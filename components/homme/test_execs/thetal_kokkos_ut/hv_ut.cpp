#include <catch2/catch.hpp>

#include <random>

#include "Types.hpp"
#include "Context.hpp"
#include "ElementsGeometry.hpp"
#include "ElementsState.hpp"
#include "ElementsDerivedState.hpp"
#include "FunctorsBuffersManager.hpp"
#include "HybridVCoord.hpp"
#include "HyperviscosityFunctorImpl.hpp"
#include "SimulationParams.hpp"
#include "SphereOperators.hpp"
#include "mpi/MpiBuffersManager.hpp"
#include "mpi/Connectivity.hpp"

#include "utilities/TestUtils.hpp"
#include "utilities/SyncUtils.hpp"
#include "utilities/ViewUtils.hpp"

using namespace Homme;

extern "C" {
void init_f90 (const int& ne,
               const Real* hyai_ptr, const Real* hybi_ptr,
               const Real* hyam_ptr, const Real* hybm_ptr,
               Real* dvv, Real* mp,
               const Real& ps0, const int& hypervis_subcycle,
               const Real& nu, const Real& nu_div, const Real& nu_top,
               const Real& nu_p, const Real& nu_s);
void init_geo_views_f90 (Real*& d_ptr,Real*& dinv_ptr,
               const Real*& phis_ptr, const Real*& gradphis_ptr,
               Real*& sphmp_ptr, Real*& rspmp_ptr,
               Real*& tVisc_ptr, Real*& sph2c_ptr,
               Real*& metdet_ptr, Real*& metinv_ptr);
void biharmonic_wk_theta_f90 (const int& np1, const Real& hv_scaling, const bool& hydrostatic,
                              const Real*& dp, const Real*& vtheta_dp,
                              const Real*& w,  const Real*& phi, const Real*& v,
                              Real*& dptens, Real*& ttens, Real*& wtens,
                              Real*& phitens, Real*& vtens);
void advance_hypervis_f90 (const int& np1, const Real& dt, const Real& eta_ave_w,
                           const Real& hv_scaling, const bool& hydrostatic,
                           Real*& v_state, Real*& w_state, Real*& vtheta_state,
                           Real*& dp_state, Real*& phinh_state);
void cleanup_f90();
} // extern "C"

// This class is basically a HyperviscosityFunctorImpl, but
// using this instead of HVF we can access protected members
// of HVF (e.g., for initialization) without exposing them in HVF.

class HVFTester : public HyperviscosityFunctorImpl {
public:
  HVFTester (const SimulationParams&     params,
             const ElementsGeometry&     geometry,
             const ElementsState&        state,
             const ElementsDerivedState& derived)
   : HyperviscosityFunctorImpl(params,geometry,state,derived)
  {
    // Zero out ref states
    Kokkos::deep_copy(m_buffers.dp_ref,0.0);
    Kokkos::deep_copy(m_buffers.theta_ref,0.0);
    Kokkos::deep_copy(m_buffers.phi_i_ref,0.0);
  }

  ~HVFTester () = default;

  void set_timestep_data (const int np1, const Real dt, const Real eta_ave_w, const bool hydrostatic)
  {
    m_data.np1 = np1;
    m_data.dt = dt/m_data.hypervis_subcycle;
    m_data.eta_ave_w = eta_ave_w;
    this->m_theta_hydrostatic_mode = hydrostatic;
  }

  void set_hv_data (const Real hv_scaling, const Real nu_ratio1, const Real nu_ratio2)
  {
    m_data.nu_ratio1 = nu_ratio1;
    m_data.nu_ratio2 = nu_ratio2;
    m_data.consthv = (hv_scaling==0.0);
  }

  using ScalarTens = ExecViewUnmanaged<Scalar*   [NP][NP][NUM_LEV]>;
#ifdef XX_NONBFB_COMING
  using ScalarTensInt = ExecViewUnmanaged<Scalar*   [NP][NP][NUM_LEV_P]>;
#else
  using ScalarTensInt = ScalarTens;
#endif
  using VectorTens = ExecViewUnmanaged<Scalar*[2][NP][NP][NUM_LEV]>;

  ScalarTens get_dptens () const { return m_buffers.dptens; }
  ScalarTens get_ttens ()  const { return m_buffers.ttens; }
  ScalarTensInt get_wtens ()  const { return m_buffers.wtens; }
  ScalarTens get_phitens ()  const { return m_buffers.phitens; }
  VectorTens get_vtens ()  const { return m_buffers.vtens; }
};

TEST_CASE("hvf", "biharmonic") {

  // Catch runs these blocks of code multiple times, namely once per each
  // session within each test case. This is problematic for Context, which
  // is a static singleton.
  // We cannot call 'create' unless we are sure the object is not already stored
  // in the context. One solution is to call 'create_if_not_there', but that's not what
  // happens in mpi_cxx_f90_interface, which is called by the geometry_interface
  // fortran module.
  // Two solutions:
  //  - cleaning up the context at the end of TEST_CASE: this would also delete
  //    the comm object in the context, so you have to re-create it.
  //    Notice, however, that the comm would *already be there* when this block
  //    of code is executed for the first time (is created in tester.cpp),
  //    so you need to check if it's there first.
  //  - change mpi_cxx_f90_interface, to create the Connectivity only if not
  //    already present.
  //
  // Among the two, the former seems cleaner, since it does not affect the
  // src folder of Homme, only the test one. So I'm going with that.
  // More precisely, I'm getting a copy of the existing Comm from the context,
  // and reset it back in it after the cleanup

  constexpr int ne = 2;

  // The random numbers generator
  std::random_device rd;
  using rngAlg = std::mt19937_64;
  // const int seed = 1; // Change to the following line after debugging
  const int seed = rd();
  rngAlg engine(seed);
  using RPDF = std::uniform_real_distribution<Real>;
  using IPDF = std::uniform_int_distribution<int>;

  // Use stuff from Context, to increase similarity with actual runs
  auto& c = Context::singleton();
  auto old_comm = c.get_ptr<Comm>();

  // Init parameters
  auto& params = c.create<SimulationParams>();
  params.nu_top            = RPDF(1e-6,1e-3)(engine);
  params.nu                = RPDF(1e-6,1e-3)(engine);
  params.nu_p              = RPDF(1e-6,1e-3)(engine);
  params.nu_s              = RPDF(1e-6,1e-3)(engine);
  params.nu_div            = RPDF(1e-6,1e-3)(engine);
  params.hypervis_scaling  = RPDF(0.1,1.0)(engine);
  params.hypervis_subcycle = IPDF(1,3)(engine);
  params.params_set = true;

  // Create and init hvcoord and ref_elem, needed to init the fortran interface
  auto& hvcoord = c.create<HybridVCoord>();
  auto& ref_FE  = c.create<ReferenceElement>();
  hvcoord.random_init(seed);

  auto hyai = Kokkos::create_mirror_view(hvcoord.hybrid_ai);
  auto hybi = Kokkos::create_mirror_view(hvcoord.hybrid_bi);
  auto hyam = Kokkos::create_mirror_view(hvcoord.hybrid_am);
  auto hybm = Kokkos::create_mirror_view(hvcoord.hybrid_bm);
  Kokkos::deep_copy(hyai,hvcoord.hybrid_ai);
  Kokkos::deep_copy(hybi,hvcoord.hybrid_bi);
  Kokkos::deep_copy(hyam,hvcoord.hybrid_am);
  Kokkos::deep_copy(hybm,hvcoord.hybrid_bm);
  const Real* hyai_ptr  = hyai.data();
  const Real* hybi_ptr  = hybi.data();
  const Real* hyam_ptr  = reinterpret_cast<Real*>(hyam.data());
  const Real* hybm_ptr  = reinterpret_cast<Real*>(hybm.data());

  std::vector<Real> dvv(NP*NP);
  std::vector<Real> mp(NP*NP);

  // This will also init the c connectivity.
  init_f90(ne,hyai_ptr,hybi_ptr,hyam_ptr,hybm_ptr,dvv.data(),mp.data(),
           hvcoord.ps0,params.hypervis_subcycle,
           params.nu, params.nu_div, params.nu_top,
           params.nu_p, params.nu_s);
  ref_FE.init_mass(mp.data());
  ref_FE.init_deriv(dvv.data());

  // Create and init elements
  const int num_elems = c.get<Connectivity>().get_num_local_elements();

  auto& geo = c.create<ElementsGeometry>();
  geo.init(num_elems,false,true);
  geo.randomize(seed);

  auto& state = c.create<ElementsState>();
  state.init(num_elems);
  const auto max_pressure = 1000.0 + hvcoord.ps0; // This ensures max_p > ps0
  state.randomize(seed,max_pressure,hvcoord.ps0);

  auto& derived = c.create<ElementsDerivedState>();
  derived.init(num_elems);
  derived.randomize(seed,RPDF(0.1,1.0)(engine));

  // Init f90
  auto d        = Kokkos::create_mirror_view(geo.m_d);
  auto dinv     = Kokkos::create_mirror_view(geo.m_dinv);
  auto phis     = Kokkos::create_mirror_view(geo.m_phis);
  auto gradphis = Kokkos::create_mirror_view(geo.m_gradphis);
  auto spmp     = Kokkos::create_mirror_view(geo.m_spheremp);
  auto rspmp    = Kokkos::create_mirror_view(geo.m_rspheremp);
  auto tVisc    = Kokkos::create_mirror_view(geo.m_tensorvisc);
  auto sph2c    = Kokkos::create_mirror_view(geo.m_vec_sph2cart);
  auto mdet     = Kokkos::create_mirror_view(geo.m_metdet);
  auto minv     = Kokkos::create_mirror_view(geo.m_metinv);
  Kokkos::deep_copy(phis,geo.m_phis);
  Kokkos::deep_copy(gradphis,geo.m_gradphis);

  Real* d_ptr        = d.data();
  Real* dinv_ptr     = dinv.data();
  Real* spmp_ptr     = spmp.data();
  Real* rspmp_ptr    = rspmp.data();
  Real* tVisc_ptr    = tVisc.data();
  Real* sph2c_ptr    = sph2c.data();
  Real* mdet_ptr     = mdet.data();
  Real* minv_ptr     = minv.data();
  const Real* phis_ptr     = phis.data();
  const Real* gradphis_ptr = gradphis.data();

  // This will also init the c connectivity.
  init_geo_views_f90(d_ptr,dinv_ptr,phis_ptr,gradphis_ptr,
                     spmp_ptr,rspmp_ptr,tVisc_ptr,
                     sph2c_ptr,mdet_ptr,minv_ptr);

  Kokkos::deep_copy(geo.m_d,d);
  Kokkos::deep_copy(geo.m_dinv,dinv);
  Kokkos::deep_copy(geo.m_spheremp,spmp);
  Kokkos::deep_copy(geo.m_rspheremp,rspmp);
  Kokkos::deep_copy(geo.m_tensorvisc,tVisc);
  Kokkos::deep_copy(geo.m_vec_sph2cart,sph2c);
  Kokkos::deep_copy(geo.m_metdet,mdet);
  Kokkos::deep_copy(geo.m_metinv,minv);

  // Get or create and init other structures needed by HVF
  auto& bmm = c.create<MpiBuffersManagerMap>();
  auto& sphop = c.create<SphereOperators>();
  FunctorsBuffersManager fbm;

  sphop.setup(geo,ref_FE);
  if (!bmm.is_connectivity_set ()) {
    bmm.set_connectivity(c.get_ptr<Connectivity>());
  }

  // Create the HVF tester
  HVFTester hvf(params,geo,state,derived);
  fbm.request_size( hvf.requested_buffer_size() );
  fbm.allocate();
  hvf.init_buffers(fbm);

  SECTION ("biharmonic_wk_theta") {
    for (Real hv_scaling : {0.0, RPDF(0.5,5.0)(engine)}) {
      for (const bool hydrostatic : {true, false}) {
        params.theta_hydrostatic_mode = hydrostatic;

        // Generate timestep settings
        const Real dt = 1.0;//rpdf(engine);
        const Real eta_ave_w = 1.0;//rpdf(engine);
        const int  np1 = 0;//ipdf(engine);
        hvf.set_timestep_data(np1,dt,eta_ave_w, hydrostatic);

        // The be needs to be inited after the hydrostatic option has been set
        hvf.init_boundary_exchanges();

        // Update hv settings
        params.hypervis_scaling = hv_scaling;
        if (params.nu != params.nu_div) {
          Real ratio = params.nu_div / params.nu;
          if (params.hypervis_scaling != 0.0) {
            params.nu_ratio1 = ratio;
            params.nu_ratio2 = 1.0;
          }else{
            params.nu_ratio1 = ratio;
            params.nu_ratio2 = 1.0;
          }
        }else{
          params.nu_ratio1 = 1.0;
          params.nu_ratio2 = 1.0;
        }

        // Set the hv scaling
        hvf.set_hv_data(hv_scaling,params.nu_ratio1,params.nu_ratio2);

        // Run kokkos version
        hvf.biharmonic_wk_theta();

        // Run fortran version
        using ScalarStateF90    = HostViewManaged<Real*[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV][NP][NP]>;
        using ScalarStateIntF90 = HostViewManaged<Real*[NUM_TIME_LEVELS][NUM_INTERFACE_LEV][NP][NP]>;
        using VectorStateF90    = HostViewManaged<Real*[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV][2][NP][NP]>;

        using ScalarViewF90 = HostViewManaged<Real*[NUM_PHYSICAL_LEV][NP][NP]>;
        using VectorViewF90 = HostViewManaged<Real*[NUM_PHYSICAL_LEV][2][NP][NP]>;

        ScalarStateF90 dp3d_f90("",num_elems);
        ScalarStateF90 vtheta_f90("",num_elems);
        ScalarStateIntF90 w_f90("",num_elems);
        ScalarStateIntF90 phi_f90("",num_elems);
        VectorStateF90 v_f90("",num_elems);

        sync_to_host(state.m_dp3d,dp3d_f90);
        sync_to_host(state.m_vtheta_dp,vtheta_f90);
        sync_to_host(state.m_w_i,w_f90);
        sync_to_host(state.m_phinh_i,phi_f90);
        sync_to_host(state.m_v,v_f90);

        const Real* dp_ptr     = reinterpret_cast<const Real*>(dp3d_f90.data());
        const Real* vtheta_ptr = reinterpret_cast<const Real*>(vtheta_f90.data());
        const Real* w_ptr      = reinterpret_cast<const Real*>(w_f90.data());
        const Real* phi_ptr    = reinterpret_cast<const Real*>(phi_f90.data());
        const Real* v_ptr      = reinterpret_cast<const Real*>(v_f90.data());

        ScalarViewF90 dptens_f90("",num_elems);
        ScalarViewF90 ttens_f90("",num_elems);
        ScalarViewF90 wtens_f90("",num_elems);
        ScalarViewF90 phitens_f90("",num_elems);
        VectorViewF90 vtens_f90("",num_elems);

        auto dptens_ptr  = dptens_f90.data();
        auto ttens_ptr   = ttens_f90.data();
        auto wtens_ptr   = wtens_f90.data();
        auto phitens_ptr = phitens_f90.data();
        auto vtens_ptr   = vtens_f90.data();
        biharmonic_wk_theta_f90(np1+1, params.hypervis_scaling, hydrostatic,
                                dp_ptr,vtheta_ptr,w_ptr,phi_ptr,v_ptr,
                                dptens_ptr,ttens_ptr,wtens_ptr, phitens_ptr,vtens_ptr);

        // Compare answers
        auto h_dptens = Kokkos::create_mirror_view(hvf.get_dptens());
        auto h_ttens = Kokkos::create_mirror_view(hvf.get_ttens());
        auto h_wtens = Kokkos::create_mirror_view(hvf.get_wtens());
        auto h_phitens = Kokkos::create_mirror_view(hvf.get_phitens());
        auto h_vtens = Kokkos::create_mirror_view(hvf.get_vtens());

        Kokkos::deep_copy(h_dptens,hvf.get_dptens());
        Kokkos::deep_copy(h_ttens,hvf.get_ttens());
        Kokkos::deep_copy(h_wtens,hvf.get_wtens());
        Kokkos::deep_copy(h_phitens,hvf.get_phitens());
        Kokkos::deep_copy(h_vtens,hvf.get_vtens());
        for (int ie=0; ie<num_elems; ++ie) {
          auto dptens_cxx = viewAsReal(Homme::subview(h_dptens,ie));
          auto ttens_cxx = viewAsReal(Homme::subview(h_ttens,ie));
          auto wtens_cxx = viewAsReal(Homme::subview(h_wtens,ie));
          auto phitens_cxx = viewAsReal(Homme::subview(h_phitens,ie));
          auto vtens_cxx = viewAsReal(Homme::subview(h_vtens,ie));
          for (int igp=0; igp<NP; ++igp) {
            for (int jgp=0; jgp<NP; ++jgp) {
              for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
                if(dptens_cxx(igp,jgp,k)!=dptens_f90(ie,k,igp,jgp)) {
                  printf("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                  printf("hv_scaling: %3.17f\n",hv_scaling);
                  printf("dptens cxx: %3.40f\n",dptens_cxx(igp,jgp,k));
                  printf("dptens f90: %3.40f\n",dptens_f90(ie,k,igp,jgp));
                }
                REQUIRE(dptens_cxx(igp,jgp,k)==dptens_f90(ie,k,igp,jgp));
                if(ttens_cxx(igp,jgp,k)!=ttens_f90(ie,k,igp,jgp)) {
                  printf("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                  printf("hv_scaling: %3.17f\n",hv_scaling);
                  printf("ttens cxx: %3.17f\n",ttens_cxx(igp,jgp,k));
                  printf("ttens f90: %3.17f\n",ttens_f90(ie,k,igp,jgp));
                }
                REQUIRE(ttens_cxx(igp,jgp,k)==ttens_f90(ie,k,igp,jgp));
                if(wtens_cxx(igp,jgp,k)!=wtens_f90(ie,k,igp,jgp)) {
                  printf("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                  printf("hv_scaling: %3.17f\n",hv_scaling);
                  printf("wtens cxx: %3.17f\n",wtens_cxx(igp,jgp,k));
                  printf("wtens f90: %3.17f\n",wtens_f90(ie,k,igp,jgp));
                }
                REQUIRE(wtens_cxx(igp,jgp,k)==wtens_f90(ie,k,igp,jgp));
                if(phitens_cxx(igp,jgp,k)!=phitens_f90(ie,k,igp,jgp)) {
                  printf("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                  printf("hv_scaling: %3.17f\n",hv_scaling);
                  printf("phitens cxx: %3.17f\n",phitens_cxx(igp,jgp,k));
                  printf("phitens f90: %3.17f\n",phitens_f90(ie,k,igp,jgp));
                }
                REQUIRE(phitens_cxx(igp,jgp,k)==phitens_f90(ie,k,igp,jgp));

                if(vtens_cxx(0,igp,jgp,k)!=vtens_f90(ie,k,0,igp,jgp)) {
                  printf("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                  printf("hv_scaling: %3.17f\n",hv_scaling);
                  printf("vtens cxx: %3.17f\n",vtens_cxx(0,igp,jgp,k));
                  printf("vtens f90: %3.17f\n",vtens_f90(ie,k,0,igp,jgp));
                }
                REQUIRE(vtens_cxx(0,igp,jgp,k)==vtens_f90(ie,k,0,igp,jgp));
                if(vtens_cxx(1,igp,jgp,k)!=vtens_f90(ie,k,1,igp,jgp)) {
                  printf("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                  printf("hv_scaling: %3.17f\n",hv_scaling);
                  printf("vtens cxx: %3.17f\n",vtens_cxx(1,igp,jgp,k));
                  printf("vtens f90: %3.17f\n",vtens_f90(ie,k,1,igp,jgp));
                }
                REQUIRE(vtens_cxx(1,igp,jgp,k)==vtens_f90(ie,k,1,igp,jgp));
              }
            }
          }
        }
      }
    }
  }

  SECTION ("hypervis") {
    for (Real hv_scaling : {0.0, RPDF(0.5,5.0)(engine)}) {
      for (const bool hydrostatic : {true, false}) {
        params.theta_hydrostatic_mode = hydrostatic;

        // Generate timestep settings
        const Real dt = RPDF(1e-5,1e-3)(engine);
        const Real eta_ave_w = 1.0;
        const int  np1 = IPDF(0,2)(engine);
        hvf.set_timestep_data(np1,dt,eta_ave_w, hydrostatic);

        // Generate random states
        state.randomize(seed);

        // The HV functor as a whole is more delicate than biharmonic_wk.
        // In particular, the EOS is used a couple of times. This means
        // that inputs *must* satisfy some minimum requirements, like
        // dp>0, vtheta>0, and d(phi)>0. This is very unlikely with random
        // inputs coming from state.randomize(seed), so we generate data
        // as "realistic" as possible, and perturb it.
        using PDF = std::uniform_real_distribution<Real>;
        ExecViewManaged<Scalar*[NP][NP][NUM_LEV_P]> perturb("",num_elems);
        genRandArray(perturb,engine,PDF(-0.05,0.05));

        EquationOfState eos;
        eos.init(hydrostatic,hvcoord);

        ElementOps elem_ops;
        elem_ops.init(hvcoord);

        auto s = state;
        auto g = geo;
        ExecViewManaged<Scalar[NUM_LEV]> buf_m("");
        ExecViewManaged<Scalar[NUM_LEV_P]> buf_i("");
        Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                             KOKKOS_LAMBDA(const TeamMember& team){
          KernelVariables kv(team);
          Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                               [&](const int idx){
            const int igp = idx / NP;
            const int jgp = idx % NP;

            auto noise = Homme::subview(perturb,kv.ie,igp,jgp);
            auto dp = Homme::subview(s.m_dp3d,kv.ie,np1,igp,jgp);
            auto theta = Homme::subview(s.m_vtheta_dp,kv.ie,np1,igp,jgp);
            auto phi = Homme::subview(s.m_phinh_i,kv.ie,np1,igp,jgp);

            // First, compute dp = dp_ref+noise
            hvcoord.compute_dp_ref(kv,s.m_ps_v(kv.ie,np1,igp,jgp),dp);
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                                 [&](const int ilev){
              dp(ilev) *= 1.0 + noise(ilev);
            });
            // Compute pressure
            elem_ops.compute_hydrostatic_p(kv,dp,buf_i,buf_m);

            // Compute vtheta_dp = theta_ref*dp
            elem_ops.compute_theta_ref(kv,buf_m,theta);
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                                 [&](const int ilev){
              theta(ilev) *= dp(ilev);
            });

            // Compute phi
            eos.compute_phi_i(kv,g.m_phis(kv.ie,igp,jgp),
                                 theta,buf_m,phi);
          });
        });

        // The be needs to be inited after the hydrostatic option has been set
        hvf.init_boundary_exchanges();

        // Copy states into f90 pointers
        HostViewManaged<Real*[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV][2][NP][NP]> v_f90("",num_elems);
        HostViewManaged<Real*[NUM_TIME_LEVELS][NUM_INTERFACE_LEV][NP][NP]>   w_f90("",num_elems);
        HostViewManaged<Real*[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV][NP][NP]>    dp_f90("",num_elems);
        HostViewManaged<Real*[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV][NP][NP]>    vtheta_f90("",num_elems);
        HostViewManaged<Real*[NUM_TIME_LEVELS][NUM_INTERFACE_LEV][NP][NP]>   phinh_f90("",num_elems);

        sync_to_host(state.m_v,v_f90);
        sync_to_host(state.m_w_i,w_f90);
        sync_to_host(state.m_dp3d,dp_f90);
        sync_to_host(state.m_vtheta_dp,vtheta_f90);
        sync_to_host(state.m_phinh_i,phinh_f90);

        Real* v_f90_ptr      = v_f90.data();
        Real* w_f90_ptr      = w_f90.data();
        Real* dp_f90_ptr     = dp_f90.data();
        Real* vtheta_f90_ptr = vtheta_f90.data();
        Real* phinh_f90_ptr  = phinh_f90.data();

        // Update hv settings
        params.hypervis_scaling = hv_scaling;
        if (params.nu != params.nu_div) {
          Real ratio = params.nu_div / params.nu;
          if (params.hypervis_scaling != 0.0) {
            params.nu_ratio1 = ratio;
            params.nu_ratio2 = 1.0;
          }else{
            params.nu_ratio1 = ratio;
            params.nu_ratio2 = 1.0;
          }
        }else{
          params.nu_ratio1 = 1.0;
          params.nu_ratio2 = 1.0;
        }

        // Set the viscosity params
        hvf.set_hv_data(hv_scaling,params.nu_ratio1,params.nu_ratio2);

        // Run the cxx functor
        hvf.run(np1,dt,eta_ave_w);

        // Run the f90 functor
        advance_hypervis_f90(np1+1,dt,eta_ave_w, hv_scaling, hydrostatic,
                             v_f90_ptr, w_f90_ptr, vtheta_f90_ptr, dp_f90_ptr, phinh_f90_ptr);


        // Compare answers
        auto v_cxx      = Kokkos::create_mirror_view(state.m_v);
        auto w_cxx      = Kokkos::create_mirror_view(state.m_w_i);
        auto vtheta_cxx = Kokkos::create_mirror_view(state.m_vtheta_dp);
        auto dp_cxx     = Kokkos::create_mirror_view(state.m_dp3d);
        auto phinh_cxx  = Kokkos::create_mirror_view(state.m_phinh_i);

        Kokkos::deep_copy(v_cxx,      state.m_v);
        Kokkos::deep_copy(w_cxx,      state.m_w_i);
        Kokkos::deep_copy(vtheta_cxx, state.m_vtheta_dp);
        Kokkos::deep_copy(dp_cxx,     state.m_dp3d);
        Kokkos::deep_copy(phinh_cxx,  state.m_phinh_i);

        for (int ie=0; ie<num_elems; ++ie) {
          for (int igp=0; igp<NP; ++igp) {
            for (int jgp=0; jgp<NP; ++jgp) {
              for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
                const int ilev = k / VECTOR_SIZE;
                const int ivec = k % VECTOR_SIZE;

                if (v_cxx(ie,np1,0,igp,jgp,ilev)[ivec]!=v_f90(ie,np1,k,0,igp,jgp)) {
                  printf ("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                  printf ("v_cxx: %3.16f\n",v_cxx(ie,np1,0,igp,jgp,ilev)[ivec]);
                  printf ("v_f90: %3.16f\n",v_f90(ie,np1,k,0,igp,jgp));
                }
                REQUIRE (v_cxx(ie,np1,0,igp,jgp,ilev)[ivec]==v_f90(ie,np1,k,0,igp,jgp));

                if (v_cxx(ie,np1,1,igp,jgp,ilev)[ivec]!=v_f90(ie,np1,k,1,igp,jgp)) {
                  printf ("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                  printf ("v_cxx: %3.16f\n",v_cxx(ie,np1,1,igp,jgp,ilev)[ivec]);
                  printf ("v_f90: %3.16f\n",v_f90(ie,np1,k,1,igp,jgp));
                }
                REQUIRE (v_cxx(ie,np1,1,igp,jgp,ilev)[ivec]==v_f90(ie,np1,k,1,igp,jgp));

                if (dp_cxx(ie,np1,igp,jgp,ilev)[ivec]!=dp_f90(ie,np1,k,igp,jgp)) {
                  printf ("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                  printf ("dp_cxx: %3.16f\n",dp_cxx(ie,np1,igp,jgp,ilev)[ivec]);
                  printf ("dp_f90: %3.16f\n",dp_f90(ie,np1,k,igp,jgp));
                }
                REQUIRE (dp_cxx(ie,np1,igp,jgp,ilev)[ivec]==dp_f90(ie,np1,k,igp,jgp));

                if (w_cxx(ie,np1,igp,jgp,ilev)[ivec]!=w_f90(ie,np1,k,igp,jgp)) {
                  printf ("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                  printf ("w_cxx: %3.16f\n",w_cxx(ie,np1,igp,jgp,ilev)[ivec]);
                  printf ("w_f90: %3.16f\n",w_f90(ie,np1,k,igp,jgp));
                }
                REQUIRE (w_cxx(ie,np1,igp,jgp,ilev)[ivec]==w_f90(ie,np1,k,igp,jgp));

                if (phinh_cxx(ie,np1,igp,jgp,ilev)[ivec]!=phinh_f90(ie,np1,k,igp,jgp)) {
                  printf ("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                  printf ("phinh_cxx: %3.16f\n",phinh_cxx(ie,np1,igp,jgp,ilev)[ivec]);
                  printf ("phinh_f90: %3.16f\n",phinh_f90(ie,np1,k,igp,jgp));
                }
                REQUIRE (phinh_cxx(ie,np1,igp,jgp,ilev)[ivec]==phinh_f90(ie,np1,k,igp,jgp));

                if (vtheta_cxx(ie,np1,igp,jgp,ilev)[ivec]!=vtheta_f90(ie,np1,k,igp,jgp)) {
                  printf ("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                  printf ("vtheta_cxx: %3.16f\n",vtheta_cxx(ie,np1,igp,jgp,ilev)[ivec]);
                  printf ("vtheta_f90: %3.16f\n",vtheta_f90(ie,np1,k,igp,jgp));
                }
                REQUIRE (vtheta_cxx(ie,np1,igp,jgp,ilev)[ivec]==vtheta_f90(ie,np1,k,igp,jgp));

              }

              // Last interface
              const int k = NUM_INTERFACE_LEV-1;
              const int ilev = k / VECTOR_SIZE;
              const int ivec = k % VECTOR_SIZE;

              if (w_cxx(ie,np1,igp,jgp,ilev)[ivec]!=w_f90(ie,np1,k,igp,jgp)) {
                printf ("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                printf ("w_cxx: %3.16f\n",w_cxx(ie,np1,igp,jgp,ilev)[ivec]);
                printf ("w_f90: %3.16f\n",w_f90(ie,np1,k,igp,jgp));
              }
              REQUIRE (w_cxx(ie,np1,igp,jgp,ilev)[ivec]==w_f90(ie,np1,k,igp,jgp));

              if (phinh_cxx(ie,np1,igp,jgp,ilev)[ivec]!=phinh_f90(ie,np1,k,igp,jgp)) {
                printf ("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                printf ("phinh_cxx: %3.16f\n",phinh_cxx(ie,np1,igp,jgp,ilev)[ivec]);
                printf ("phinh_f90: %3.16f\n",phinh_f90(ie,np1,k,igp,jgp));
              }
              REQUIRE (phinh_cxx(ie,np1,igp,jgp,ilev)[ivec]==phinh_f90(ie,np1,k,igp,jgp));
            }
          }
        }
      }
    }
  }

  c.finalize_singleton();
  auto& new_comm = c.create<Comm>();
  new_comm = *old_comm;
  cleanup_f90();
}
