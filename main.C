// Copyright (c) 2018-2019, SofTMech.org
// All rights reserved.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE


// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>
#include <stdio.h>
#include <time.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFECentroidPostProcessor.h>
#include <ibamr/IBFEMethod.h>

#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredNonConservativeHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>

#include <ibamr/app_namespaces.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petsc.h>  //This is from wchen YJ
#include <petscsys.h>  //This is from Scott YJ

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

#include <ibamr/INSVCStaggeredHierarchyIntegrator.h>
#include <ibamr/INSVCStaggeredNonConservativeHierarchyIntegrator.h>  //This is from Scott YJ

#include <ibamr/app_namespaces.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Headers for basic libMesh objects
#include <libmesh/exodusII_io.h>
#include <libmesh/mesh.h>
#include <libmesh/equation_systems.h>

// Headers for application-specific algorithm/data structure objects
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBFEMethod.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>
#include <ibamr/IBFECentroidPostProcessor.h>


// Application
#include "SetFluidPropertiesBySolid.h"
#include "ModelInitialization.h"
#include "BoundaryConditions.h"
#include "ExternalSources.h"
#include "MechanicsModel.h"
#include "sourceToIB.h"
#include "DarcyMass.h"
#include "DarcyPressure.h"
#include "DarcyVelocity.h"
#include "MapMeshToFluid.h"
#include "MapmAddedToFluid.h"
#include "ActiveContraction.h"
#include "FeedbackToNetwork.h"
#include "MeshConversion.h"
#include "ConoraryFlowInterface.h"
#include "AHAdivision.h"
#include "sourceToDarcy.h"
#include "AdvanceCoronaryCode.h"

#include "Coro.h"
#include "tools.h"
#include "arteriesPD.h"

//for linux system call
#include <stdlib.h>
#include <stdio.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>
#include <chrono>

// to calculate the volume
#include "tetVolumeCalculation.h"
#include "tetgen.h"

using namespace std;

extern "C" void impedance_driver_(int *tmstps, double *Period, double *rho, double *mu,
                                  double *r_root, double *rmin,
                                  double *y11, double *y12, double *y21, double *y22,
                                  double *ye11,double *ye12,double *ye21,double *ye22,
                                  double *Lr, double *q, double *g,
                                  double *fa1, double *fa2, double *fa3,
                                  double *fv1, double *fv2, double *fv3,
                                  double *asym, double *expo,
                                  double *lrrA, double *lrrV);


/*******************************************************************************

 *******************************************************************************/

// Function prototypes
void output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 Pointer<IBHierarchyIntegrator> time_integrator,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname);

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/
bool
run_example(int argc, char* argv[], Tube *Arteries[])
{
    // Initialize PETSc, MPI, and SAMRAI.
    //PetscInitialize(&argc, &argv, NULL, NULL);

    LibMeshInit init(argc, argv); //initialize libMesh, which will initlize Petsc

    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    // Increase maximum patch data component indices
    SAMRAIManager::setMaxNumberPatchDataEntries(2500);

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "INS.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

    // Get various standard options set in the input file.
    const bool dump_viz_data = app_initializer->dumpVizData();
    const int viz_dump_interval = app_initializer->getVizDumpInterval();
    const bool uses_visit = dump_viz_data && !app_initializer->getVisItDataWriter().isNull();
    const bool uses_exodus = dump_viz_data && !app_initializer->getExodusIIFilename().empty();
    const string exodus_filename = app_initializer->getExodusIIFilename();

    const bool dump_restart_data = app_initializer->dumpRestartData();
    const int restart_dump_interval = app_initializer->getRestartDumpInterval();
    const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();
    const string restart_read_dirname = app_initializer->getRestartReadDirectory();
    const int restart_restore_num = app_initializer->getRestartRestoreNumber();

    const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
    const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
    const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();

        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();
        // till here everything is similar as the IB/FE code



        // Configure the passive stress model.
        MechanicsModel::normalize_stress = input_db->getBool("NORMALIZE_STRESS");


        // Configure the active tension model.
        MechanicsModel::enable_active_tension = input_db->getBool("ENABLE_ACTIVE_TENSION");
        MechanicsModel::T_scale = input_db->getDouble("T_SCALE");
        MechanicsModel::beta_s = input_db->getDouble("BETA_S");


        BoundaryConditions::P_load = input_db->getDouble("P_LOAD");
        BoundaryConditions::t_load = input_db->getDouble("T_LOAD");
        BoundaryConditions::kappa = input_db->getDouble("KAPPA");
        BoundaryConditions::t_end_diastole = input_db->getDouble("T_END_DIASTOLE");
        BoundaryConditions::diastole_systole_difference = input_db->getDouble("DIASTOLE_SYSTOLE_DIFFERENCE");

        BoundaryConditions::cycle_no = input_db->getDouble("NO_OF_CYCLES");
        BoundaryConditions::t_period = input_db->getDouble("TIME_PERIOD");
        BoundaryConditions::LV_CIsoContract = input_db->getDouble("CIsoContract");
        BoundaryConditions::P_AV_open = input_db->getDouble("P_AVOPEN");
        BoundaryConditions::Pwk_AV = input_db->getDouble("P_AVOPEN");

        BoundaryConditions::pao_initial = input_db->getDouble("P_AO_OPEN");
        BoundaryConditions::pwk_initial = input_db->getDouble("P_AO_OPEN");

        BoundaryConditions::P_Ao = input_db->getDouble("P_AO_OPEN");

        BoundaryConditions::LV_CIsoRelax = input_db->getDouble("LV_CIsoRelax");
        BoundaryConditions::current_cycle=1;
        BoundaryConditions::cycle_time_saved=0.0;
        BoundaryConditions::single_cycle_length=input_db->getDouble("SINGLE_CYCLE_LENGTH");
        BoundaryConditions::total_node_num = input_db->getInteger("TOTAL_NODE_NUM");
        BoundaryConditions::total_elem_num = input_db->getInteger("TOTAL_ELEM_NUM");
        BoundaryConditions::total_endo_nodes = input_db->getInteger("TOTAL_ENDO_NODES");

        DarcyMass::phi0=input_db->getDouble("PHI_0");
        DarcyMass::biot_modulus=input_db->getDouble("BIOT_MODULUS");
        DarcyMass::kappa0=input_db->getDouble("KAPPA_0");


        DarcyPressure::phi0=input_db->getDouble("PHI_0");
        DarcyPressure::biot_modulus=input_db->getDouble("BIOT_MODULUS");
        DarcyPressure::kappa0=input_db->getDouble("KAPPA_0");
        DarcyPressure::cap_K=input_db->getDouble("CAP_K");


        sourceToDarcy::beta_a=input_db->getDouble("BETA_A");
        sourceToDarcy::beta_v=input_db->getDouble("BETA_V");
        sourceToDarcy::p_a=input_db->getDouble("P_A");
        sourceToDarcy::p_v=input_db->getDouble("P_V");

        ConoraryFlowInterface::send_source_from_coro = input_db->getBool("SEND_SOURCE_FROM_CORO");
        ConoraryFlowInterface::visualise_coro_data = input_db->getBool("VISUALISE_CORO_DATA");
        ConoraryFlowInterface::output_coro_data_to_screen = input_db->getBool("OUTPUT_CORO_DATA_TO_SCREEN");
        ConoraryFlowInterface::output_coro_data_to_file = input_db->getBool("OUTPUT_CORO_DATA_TO_FILE");
        ConoraryFlowInterface::aha_to_subdomain = input_db->getBool("AHA_TO_SUBDOMAINS");
        ConoraryFlowInterface::exp_pao_to_file = input_db->getBool("EXP_PAO_TO_FILE");


        //load and initialize the FE mesh
        Mesh mesh(init.comm(), NDIM);
        ModelInitialization::initialize_mesh(mesh, input_db);

    // Create major IBAMR solver objects.
    Pointer<IBFEMethod> ib_method_ops = new IBFEMethod("IBFEMethod",
                     app_initializer->getComponentDatabase("IBFEMethod"),
                     &mesh,
                     app_initializer->getComponentDatabase("GriddingAlgorithm")->getInteger("max_levels"),
                    /*register_for_restart*/  true, restart_read_dirname, restart_restore_num);
    Pointer<INSHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
                     "INSStaggeredHierarchyIntegrator",
                     app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));
    Pointer<IBHierarchyIntegrator> time_integrator = new IBExplicitHierarchyIntegrator(
                     "IBExplicitHierarchyIntegrator",
                     app_initializer->getComponentDatabase("IBExplicitHierarchyIntegrator"),
                     ib_method_ops,
                     navier_stokes_integrator);

    // Create major SAMRAI algorithm and data objects.
    Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
                      "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>(
                      "PatchHierarchy", grid_geometry);
    Pointer<StandardTagAndInitialize<NDIM> > error_detector =
                      new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                      time_integrator,
                      app_initializer->getComponentDatabase("StandardTagAndInitialize"));
    Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
    Pointer<LoadBalancer<NDIM> > load_balancer = new LoadBalancer<NDIM>(
                       "LoadBalancer",
                       app_initializer->getComponentDatabase("LoadBalancer"));
    Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm = new GriddingAlgorithm<NDIM>(
                        "GriddingAlgorithm",
                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                        error_detector, box_generator, load_balancer);

        // Configure the IBFE solver
        pout << "Configuring the solver...\n";
        //figure out the integration order
//        const Order full_quad_order =
//        Utility::string_to_enum<Order>(input_db->getString("FULL_QUAD_ORDER"));
//        const Order reduced_quad_order =
//        Utility::string_to_enum<Order>(input_db->getString("REDUCED_QUAD_ORDER"));

        IBFEMethod::PK1StressFcnData PK1_dev_stress_data;
        // prepare the data transfered to PK1 function
        std::vector<int> vars_fibre(3);
        vars_fibre[0] = 0;
        vars_fibre[1] = 1;
        vars_fibre[2] = 2;
        std::vector<int> vars_sheet(3);
        vars_sheet[0] = 0;
        vars_sheet[1] = 1;
        vars_sheet[2] = 2;
        std::vector<int> vars_acT(1);
        vars_acT[0] = 0;
        std::vector<SystemData> fsnAct_data(3);
        fsnAct_data[0] = SystemData("fiber direction", vars_fibre);
        fsnAct_data[1] = SystemData("sheet direction", vars_sheet);
        fsnAct_data[2] = SystemData("active tension", vars_acT);

        PK1_dev_stress_data.fcn = MechanicsModel::PK1_dev_stress_function;
        PK1_dev_stress_data.quad_type = QGAUSS;
        PK1_dev_stress_data.quad_order = THIRD;  // full-order integration
        PK1_dev_stress_data.system_data = fsnAct_data;
        ib_method_ops->registerPK1StressFunction(PK1_dev_stress_data);

    std::vector<int> vars_dis(NDIM);
    for (unsigned int d = 0; d < NDIM; ++d) vars_dis[d] = d;
    vector<SystemData> sys_data_dis(1, SystemData(IBFEMethod::COORD_MAPPING_SYSTEM_NAME, vars_dis));
    IBFEMethod::LagSurfacePressureFcnData surface_pressure_data(
              BoundaryConditions::loading_force_function, sys_data_dis, NULL);
    ib_method_ops->registerLagSurfacePressureFunction(surface_pressure_data);
    IBFEMethod::LagSurfaceForceFcnData surface_fcn_data(
              BoundaryConditions::tether_force_function,sys_data_dis,NULL);
    ib_method_ops->registerLagSurfaceForceFunction(surface_fcn_data);

    ib_method_ops->initializeFEEquationSystems();
    EquationSystems* equation_systems =
                         ib_method_ops->getFEDataManager()->getEquationSystems();
    ModelInitialization::initialize_equation_systems(equation_systems);

        //adding the system data for PK1_dil stress, the pore pressure and the added M

        std::vector<int> vars_P(1);
        vars_P[0] = 0;
        std::vector<int> vars_M(1);
        vars_M[0] = 0;

        std::vector<int> vars_fibres(3);
        vars_fibres[0] = 0;
        vars_fibres[1] = 1;
        vars_fibres[2] = 2;
        std::vector<int> vars_sheets(3);
        vars_sheets[0] = 0;
        vars_sheets[1] = 1;
        vars_sheets[2] = 2;
        std::vector<int> vars_acTs(1);
        vars_acTs[0] = 0;

        std::vector<SystemData> darcy_P_M_data(5);

        darcy_P_M_data[0] = SystemData(DarcyPressure::pSystemName, vars_P);
        darcy_P_M_data[1] = SystemData(DarcyMass::mSystemName, vars_M);

        darcy_P_M_data[2] = SystemData("fiber direction", vars_fibres);
        darcy_P_M_data[3] = SystemData("sheet direction", vars_sheets);
        darcy_P_M_data[4] = SystemData("active tension", vars_acTs);

        IBFEMethod::PK1StressFcnData PK1_dil_stress_data;
        PK1_dil_stress_data.fcn = MechanicsModel::PK1_dil_stress_function;
        PK1_dil_stress_data.quad_type = QGAUSS;
        PK1_dil_stress_data.quad_order = FIRST;  // reduced-order integration
        PK1_dil_stress_data.system_data = darcy_P_M_data;
        ib_method_ops->registerPK1StressFunction(PK1_dil_stress_data); //

        //prepare the system data for transfering the variables, but not needed
        std::vector<int> vars(NDIM);
        for (unsigned int d = 0; d < NDIM; ++d) vars[d] = d;
        vector<SystemData> sys_data(1, SystemData(IBFEMethod::COORD_MAPPING_SYSTEM_NAME, vars));
//        IBFEMethod::LagSurfacePressureFcnData surface_pressure_data(BoundaryConditions::loading_force_function, sys_data, NULL);
        ib_method_ops->registerLagSurfacePressureFunction(surface_pressure_data);
//        IBFEMethod::LagSurfaceForceFcnData surface_fcn_data(BoundaryConditions::tether_force_function,sys_data,NULL);
        ib_method_ops->registerLagSurfaceForceFunction(surface_fcn_data);

        ib_method_ops->initializeFEEquationSystems();

        //*************************************************************************************************
        //set up the link with darcy solver for interpolating to the IB solver
        std::vector<int> vars_Sstar(1);
        vars_Sstar[0] = 0;
        std::vector<SystemData> darcy_Sstar_data(1);
        darcy_Sstar_data[0] = SystemData(sourceToIB::SstarSystemName, vars_Sstar);

        IBFEMethod::LagBodySourceFcnData lag_body_source_data;

        lag_body_source_data.fcn = ExternalSources::external_source_function;
        lag_body_source_data.system_data = darcy_Sstar_data;
        ib_method_ops->registerLagBodySourceFunction(lag_body_source_data);
        //*************************************************************************************************
 //       EquationSystems* equation_systems = ib_method_ops->getFEDataManager()->getEquationSystems();

        //setup from ModelInitialization
        ModelInitialization::initialize_equation_systems(equation_systems);
        //set up the phi_mesh system in the solid region
        MapMeshToFluid::initialize_equation_systems(equation_systems);
        //set up the m_added system in the solid region
        MapmAddedToFluid::initialize_equation_systems(equation_systems);

        DarcyPressure::setup_pstar(equation_systems);
        DarcyPressure::setup_p(equation_systems);

        DarcyMass::setup_m(equation_systems);
        DarcyMass::setup_mstar(equation_systems);
        DarcyMass::setup_mextra(equation_systems);
        DarcyMass::setup_maverage(equation_systems);

        DarcyVelocity::setup_wx(equation_systems);
        DarcyVelocity::setup_wy(equation_systems);
        DarcyVelocity::setup_wz(equation_systems);
        DarcyVelocity::setup_wxstar(equation_systems);
        DarcyVelocity::setup_wystar(equation_systems);
        DarcyVelocity::setup_wzstar(equation_systems);

        sourceToIB::setup_Sstar(equation_systems);

        sourceToDarcy::setup_std(equation_systems);

        FeedbackToNetwork::setup_feedback(equation_systems);

      	ConoraryFlowInterface::setup_Q_elem(equation_systems);
        ConoraryFlowInterface::setup_Pext(equation_systems);
        ConoraryFlowInterface::setup_id_test(equation_systems);
        ConoraryFlowInterface::setup_visualise_coro_data(equation_systems);

        AHAdivision::setup_mesh_division(equation_systems);

        MeshConversion::setup_mesh_long(equation_systems);
        MeshConversion::setup_mesh_radial(equation_systems);
        MeshConversion::setup_mesh_circum(equation_systems);

        // Set up a post-processor to reconstruct the various quantities of interest.
        Pointer<IBFEPostProcessor> ib_post_processor =
        new IBFECentroidPostProcessor("IBFEPostProcessor", ib_method_ops->getFEDataManager());

        //1) Deformation gradient tensor FF = dX/ds.
        ib_post_processor->registerTensorVariable(
            "FF", MONOMIAL, CONSTANT, IBFEPostProcessor::FF_fcn);

        // 2) Deformed fiber and sheet axes.
        vector<SystemData> f_system_data(1);
        f_system_data[0] = fsnAct_data[0];
        ib_post_processor->registerVectorVariable(
           "f", MONOMIAL, CONSTANT, IBFEPostProcessor::deformed_material_axis_fcn,f_system_data);
        vector<SystemData> s_system_data(1);
        s_system_data[0] = fsnAct_data[1];
        ib_post_processor->registerVectorVariable(
               "s", MONOMIAL, CONSTANT, IBFEPostProcessor::deformed_material_axis_fcn,s_system_data);

        // 3) Fiber and sheet stretches.
        ib_post_processor->registerScalarVariable(
                   "lambda_f", MONOMIAL, CONSTANT, IBFEPostProcessor::material_axis_stretch_fcn,f_system_data);
        ib_post_processor->registerScalarVariable(
                   "lambda_s", MONOMIAL, CONSTANT, IBFEPostProcessor::material_axis_stretch_fcn,s_system_data);


        // 4) Cauchy stress sigma = (1/J) PP FF^T.
        ib_post_processor->registerTensorVariable(
                       "sigma_dev", MONOMIAL, CONSTANT,
                       IBFEPostProcessor::cauchy_stress_from_PK1_stress_fcn,
                        PK1_dev_stress_data.system_data, &PK1_dev_stress_data);
        ib_post_processor->registerTensorVariable(
                        "sigma_dil", MONOMIAL, CONSTANT,
                        IBFEPostProcessor::cauchy_stress_from_PK1_stress_fcn,
                        PK1_dil_stress_data.system_data, &PK1_dil_stress_data);


        //5) Eulerian pressure p_f.
        Pointer<hier::Variable<NDIM> > p_var = time_integrator->getPressureVariable();
        Pointer<VariableContext> p_current_ctx = time_integrator->getCurrentContext();
        HierarchyGhostCellInterpolation::InterpolationTransactionComponent p_ghostfill(
                            /*data_idx*/ -1,
                            "LINEAR_REFINE",
                            /*use_cf_bdry_interpolation*/ false,
                            "CONSERVATIVE_COARSEN",
                            "LINEAR");
        FEDataManager::InterpSpec p_interp_spec("PIECEWISE_LINEAR",
                              QGAUSS,
                              FIFTH,
                              /*use_adaptive_quadrature*/ false,
                              /*point_density*/ 2.0,
                              /*use_consistent_mass_matrix*/ true,
                              /*use_nodal_quadrature*/ false);

        ib_post_processor->registerInterpolatedScalarEulerianVariable("p_f", LAGRANGE, FIRST, p_var, p_current_ctx,
                                                                                  p_ghostfill, p_interp_spec);

        //if we donot use the level set, just to reinitliaze the regions according to the mesh position
        VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
        Pointer<VariableContext> ctx = var_db->getContext("level_set_mesh_context");
        Pointer<CellVariable<NDIM, double> > phi_mesh_var = new CellVariable<NDIM, double>(MapMeshToFluid::phi_mesh_var_name);
        const IntVector<NDIM> ib_ghosts = ib_method_ops->getMinimumGhostCellWidth();
        const int phi_mesh_idx = var_db->registerVariableAndContext(phi_mesh_var, ctx, ib_ghosts);
        pout << "declare the phi_mesh_variable\n";

        Pointer<VariableContext> ctx_m_added = var_db->getContext("m_added_context");
        Pointer<CellVariable<NDIM, double> > m_added_var = new CellVariable<NDIM, double>(MapmAddedToFluid::m_added_var_name);
        const int m_added_idx = var_db->registerVariableAndContext(m_added_var, ctx_m_added, ib_ghosts);
        pout << "declare the m_added variable\n";

		// Create Eulerian boundary condition specification objects when needed.
    const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
    std::vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
    const bool periodic_boundaries = periodic_shift.min() > 0;
    if (!periodic_boundaries)
    {
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            ostringstream bc_coefs_name_stream;
            bc_coefs_name_stream << "u_bc_coefs_" << d;
            const string bc_coefs_name = bc_coefs_name_stream.str();
            ostringstream bc_coefs_db_name_stream;
            bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
            const string bc_coefs_db_name = bc_coefs_db_name_stream.str();
            u_bc_coefs[d] = new muParserRobinBcCoefs(bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
        }
        navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
    }

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        visit_data_writer->registerPlotQuantity(phi_mesh_var->getName(), "SCALAR", phi_mesh_idx);
        visit_data_writer->registerPlotQuantity(m_added_var->getName(), "SCALAR", m_added_idx);
        if (uses_visit)
        {
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }
        UniquePtr<ExodusII_IO> exodus_io(uses_exodus? new ExodusII_IO(mesh): NULL);

        // Initialize FE data.
        ib_method_ops->initializeFEData();
        ib_post_processor->initializeFEData();

        // reading in the fibre and sheet directions
        ModelInitialization::initialize_material_axes(mesh, equation_systems,input_db);

        //initialize the Darcy system here
        DarcyPressure::initialize_p(equation_systems);

        DarcyMass::initialize_m(equation_systems);
        DarcyMass::initialize_mstar(equation_systems);
        DarcyMass::initialize_mextra(equation_systems);
        DarcyMass::initialize_maverage(equation_systems);

        DarcyVelocity::initialize_wx(equation_systems);
        DarcyVelocity::initialize_wy(equation_systems);
        DarcyVelocity::initialize_wz(equation_systems);

        DarcyVelocity::initialize_wxstar(equation_systems);
        DarcyVelocity::initialize_wystar(equation_systems);
        DarcyVelocity::initialize_wzstar(equation_systems);

        sourceToIB::initialize_Sstar(equation_systems);
        DarcyPressure::initialize_pstar(equation_systems);

        FeedbackToNetwork::initialize_feedback(equation_systems);

        sourceToDarcy::initialize_std(equation_systems);

      	ConoraryFlowInterface::initialize_Q_elem(equation_systems);
        ConoraryFlowInterface::initialize_Pext(equation_systems);
        ConoraryFlowInterface::initialize_id_test(equation_systems);
        ConoraryFlowInterface::initialize_visualise_coro_data(equation_systems);


        AHAdivision::initialize_mesh_division(equation_systems);

        MeshConversion::initialize_mesh_long(equation_systems);
        MeshConversion::initialize_mesh_radial(equation_systems);
        MeshConversion::initialize_mesh_circum(equation_systems);

        //initialize the phi system
        MapMeshToFluid::initialize_phi(equation_systems);
        //intiialize the m_added_system
        MapmAddedToFluid::initialize_m_added(equation_systems);
        //initialize the J_system to 1
        MechanicsModel::initialize_J_system(equation_systems);
        //initialize the FF_system to an identity add_matrix
        MechanicsModel::initialize_FF_system(equation_systems);


        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);

        //allocate data on each level for phi_mesh_idx, phi_mesh_var, and m_added var
        for (int ln=0; ln <= patch_hierarchy->getFinestLevelNumber(); ln++)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            level->allocatePatchData(phi_mesh_idx, 0.0);
            level->allocatePatchData(m_added_idx, 0.0);
        }
        //initialize to be zero for phi_mesh_idx
        for (int ln=0; ln<= patch_hierarchy->getFinestLevelNumber(); ln++)
        {
            Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
            for(PatchLevel<NDIM>::Iterator p(level); p; p++)
            {
                Pointer<Patch<NDIM> > patch = level->getPatch(p());
                const Box<NDIM>& patch_box = patch->getBox();
                Pointer<CellData<NDIM, double> > D_data = patch->getPatchData(phi_mesh_idx);
                Pointer<CellData<NDIM, double> > m_added_data = patch->getPatchData(m_added_idx);
                for(Box<NDIM>::Iterator it(patch_box); it; it++ )
                {
                    CellIndex<NDIM> ci(it());
                    (*D_data)(ci) = -1;
                    (*m_added_data)(ci) = -1.0e-6;
                }
            }
        }

        // Remove the AppInitializer
        app_initializer.setNull();

        //update J_system
        //MechanicsModel::update_J_system(equation_systems);
        //update the FF systems
        //MechanicsModel::update_FF_system(equation_systems);

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data)
        {
            pout << "\n\nWriting visualization files...\n\n";
            if (uses_visit)
            {
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            }
            if (uses_exodus)
            {
                pout<< "write out the exodus file\n";
                ib_post_processor->postProcessData(loop_time);
                exodus_io->write_timestep(exodus_filename, *equation_systems,
                                          iteration_num/viz_dump_interval+1, loop_time);
            }
        }



BoundaryConditions::readingPointsGeneral(mesh,
                                BoundaryConditions::LV_endo_points_list,
                                BoundaryConditions::LV_endo_points,
                                BoundaryConditions::LV_all_points,
                                BoundaryConditions::LV_NoOfEndoNode,
                                input_db->getString("ENDO_POINTS_LIST"));


BoundaryConditions::updatePointsPositionGeneral(equation_systems,
                            BoundaryConditions::LV_endo_points_list,
                            BoundaryConditions::LV_endo_points,
                            BoundaryConditions::LV_all_points,
                            BoundaryConditions::LV_NoOfEndoNode);


                            BoundaryConditions::solve_wall_volume(equation_systems);


                            int my_rank20;
                            MPI_Comm_rank(MPI_COMM_WORLD, &my_rank20);

                            // Determine the rank of the broadcast emitter process
                            int broadcast_root20 = 0;
                             MPI_Barrier(init.comm().get());
                            if(my_rank20 == broadcast_root20)
                            {

                        		BoundaryConditions::whole_volume = tetVolumeCalculationByPoints(BoundaryConditions::LV_all_points);
                            BoundaryConditions::LV_volume=BoundaryConditions::whole_volume-BoundaryConditions::wall_volume;

                        		std::ofstream pv_data_file;
                        		pv_data_file.open("pressure_volume_data.dat", std::ofstream::out);
                        		pv_data_file << "#time \t" <<"pressure \t"<<"LV volume\t"<<std::endl;

                            pv_data_file <<loop_time<<"    "<<BoundaryConditions::P_current_loading<<"    "<<BoundaryConditions::LV_volume<<std::endl;
                        		pv_data_file.close();

            pout << "checking the lv volume on proc 0 " << mesh.processor_id() << "  " <<  BoundaryConditions::LV_volume << "\n";

                        	  }

                           MPI_Bcast(&BoundaryConditions::LV_volume, 1, MPI_DOUBLE, broadcast_root20, MPI_COMM_WORLD);
MPI_Barrier(init.comm().get());


                            printf("processor %d, ini volume is %f\n", mesh.processor_id(), BoundaryConditions::LV_volume);


double systemtime = loop_time;

ActiveContraction::time_to_enter=0;


int my_rank;
MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
int broadcast_root = 0;


if(my_rank == broadcast_root)
{
ofstream outdata_coro;
ofstream outdata_coro_solo;
ofstream outdata_pressures;
ofstream outdata_external_pressures;
ofstream outdata_volume_data;
ofstream outdata_average_data;
ofstream outdata_imp_holder_period;
ofstream outdata_flow_different_levels;
ofstream outdata_pressure_different_levels;
ofstream outdata_coupled_flow_different_levels;
ofstream outdata_coupled_pressure_different_levels;
outdata_coro.open("/xlwork1/scott/full_model_nov_22/coronary_data_coupled_run_UPDATED_TREE_RUNS_FOR_PAPER_MECHANICS_60PERCENT.dat", std::ofstream::trunc); // opens the file
outdata_coro_solo.open("/xlwork1/scott/full_model_nov_22/coronary_data_solo_run_UPDATED_TREE_RUNS_FOR_PAPER_MECHANICS_60PERCENT.dat", std::ofstream::trunc); // opens the file
outdata_pressures.open("/xlwork1/scott/full_model_nov_22/pressure_data_coupled_run_UPDATED_TREE_RUNS_FOR_PAPER_MECHANICS_60PERCENT.dat", std::ofstream::trunc); // opens the file
outdata_external_pressures.open("/xlwork1/scott/full_model_nov_22/external_pressure_UPDATED_TREE_RUNS_FOR_PAPER_MECHANICS_60PERCENT.dat", std::ofstream::trunc); // opens the file
outdata_volume_data.open("/xlwork1/scott/full_model_nov_22/volume_data_UPDATED_TREE_RUNS_FOR_PAPER_MECHANICS_60PERCENT.dat", std::ofstream::trunc); // opens the file
outdata_average_data.open("/xlwork1/scott/full_model_nov_22/average_data_UPDATED_TREE_RUNS_FOR_PAPER_MECHANICS_60PERCENT.dat", std::ofstream::trunc); // opens the file
outdata_imp_holder_period.open("/xlwork1/scott/full_model_nov_22/imp_holder_period_UPDATED_TREE_RUNS_FOR_PAPER_MECHANICS_60PERCENT.dat", std::ofstream::trunc); // opens the file
outdata_flow_different_levels.open("/xlwork1/scott/full_model_nov_22/flow_different_levels_UPDATED_TREE_RUNS_FOR_PAPER_MECHANICS_60PERCENT.dat", std::ofstream::trunc);
outdata_pressure_different_levels.open("/xlwork1/scott/full_model_nov_22/pressure_different_levels_UPDATED_TREE_RUNS_FOR_PAPER_MECHANICS_60PERCENT.dat", std::ofstream::trunc);
outdata_coupled_flow_different_levels.open("/xlwork1/scott/full_model_nov_22/coupled_flow_different_levels_UPDATED_TREE_RUNS_FOR_PAPER_MECHANICS_60PERCENT.dat", std::ofstream::trunc);
outdata_coupled_pressure_different_levels.open("/xlwork1/scott/full_model_nov_22/coupled_pressure_different_levels_UPDATED_TREE_RUNS_FOR_PAPER_MECHANICS_60PERCENT.dat", std::ofstream::trunc);
outdata_coro.close();
outdata_coro_solo.close();
outdata_pressures.close();
outdata_external_pressures.close();
outdata_volume_data.close();
outdata_average_data.close();
outdata_imp_holder_period.close();
outdata_flow_different_levels.close();
outdata_pressure_different_levels.close();
outdata_coupled_flow_different_levels.close();
outdata_coupled_pressure_different_levels.close();
}



//--------------------------------------------------------------------------------------------------------

        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {


//print out some information to the log file about how long each iteration takes
          if(my_rank == broadcast_root)
          {
          time_t rawtime;
          struct tm * timeinfo;

          time( &rawtime );
          timeinfo = localtime( &rawtime );

          pout<< "information about time ... "<<asctime (timeinfo)<< endl;
          }


            iteration_num = time_integrator->getIntegratorStep()+1;


            if(iteration_num==1)
            {
            loop_time = time_integrator->getIntegratorTime();
            }


            if(my_rank == broadcast_root)
            {
              if(iteration_num==1)
              {

              Tube::initial_or_integrated=0;
              Tube::tstart1 = 0.0;
              Tube::finaltime1 = numHeartBeats*Tube::Period;      // Final end-time during a simulation.
              Tube::tend1 = Tube::tstart1 + Tube::onePeriodDt;
              Tube::current_position_in_coronary=0;

              Tube::counter_sec=0;

                  double t0 =0.0;



              while (Tube::tend1 <= Tube::finaltime1)
              {

              solver (Arteries, Tube::tstart1, Tube::tend1, Tube::onePeriodDt,t0);

              Tube::tstart1 = Tube::tend1;
              Tube::tend1 = Tube::tstart1 + Tube::onePeriodDt;
              Tube::counter_sec=Tube::counter_sec+1;


              if (Tube::counter_sec % 100 ==0)
              {
                ofstream outdata;
                outdata.open("/xlwork1/scott/full_model_nov_22/coronary_data_solo_run_UPDATED_TREE_RUNS_FOR_PAPER_MECHANICS_60PERCENT.dat", std::ofstream::app); // opens the file
                for (int i=0; i < nbrves; i ++)
               {
              outdata << Tube::tstart1 <<"\t"<< i <<"\t"<< Arteries[i]->Qnew[Arteries[i]->N]*Tube::q << endl;
               }
                outdata.close();

                int arteries_0_position=floor(0.5*Arteries[0]->N);
                int arteries_1_position=floor(0.5*Arteries[1]->N);
                int arteries_2_position=floor(0.5*Arteries[2]->N);
                int arteries_3_position=floor(0.5*Arteries[3]->N);
                int arteries_4_position=floor(0.5*Arteries[4]->N);
                int arteries_5_position=floor(0.5*Arteries[5]->N);
                int arteries_6_position=floor(0.5*Arteries[6]->N);
                int arteries_7_position=floor(0.5*Arteries[7]->N);
                int arteries_8_position=floor(0.5*Arteries[8]->N);
                int arteries_9_position=floor(0.5*Arteries[9]->N);
                int arteries_10_position=floor(0.5*Arteries[10]->N);
                int arteries_11_position=floor(0.5*Arteries[11]->N);
                int arteries_12_position=floor(0.5*Arteries[12]->N);
                int arteries_13_position=floor(0.5*Arteries[13]->N);
                int arteries_14_position=floor(0.5*Arteries[14]->N);
                int arteries_15_position=floor(0.5*Arteries[15]->N);
                int arteries_16_position=floor(0.5*Arteries[16]->N);
                int arteries_17_position=floor(0.5*Arteries[17]->N);
                int arteries_18_position=floor(0.5*Arteries[18]->N);
                int arteries_19_position=floor(0.5*Arteries[19]->N);
                int arteries_20_position=floor(0.5*Arteries[20]->N);
                int arteries_21_position=floor(0.5*Arteries[21]->N);
                int arteries_22_position=floor(0.5*Arteries[22]->N);
                int arteries_23_position=floor(0.5*Arteries[23]->N);

                ofstream outdata_flow_different_levels;
                outdata_flow_different_levels.open("/xlwork1/scott/full_model_nov_22/flow_different_levels_UPDATED_TREE_RUNS_FOR_PAPER_MECHANICS_60PERCENT.dat", std::ofstream::app);
                outdata_flow_different_levels<< Tube::tstart1 <<"\t" << Tube::current_position_in_coronary <<"\t"
                << Arteries[0]->Qnew[arteries_0_position]*Tube::q << "\t"
                << Arteries[1]->Qnew[arteries_1_position]*Tube::q <<"\t"
                << Arteries[2]->Qnew[arteries_2_position]*Tube::q <<"\t"
                << Arteries[3]->Qnew[arteries_3_position]*Tube::q <<"\t"
                << Arteries[4]->Qnew[arteries_4_position]*Tube::q <<"\t"
                << Arteries[5]->Qnew[arteries_5_position]*Tube::q <<"\t"
                << Arteries[6]->Qnew[arteries_6_position]*Tube::q <<"\t"
                << Arteries[7]->Qnew[arteries_7_position]*Tube::q <<"\t"
                << Arteries[8]->Qnew[arteries_8_position]*Tube::q <<"\t"
                << Arteries[9]->Qnew[arteries_9_position]*Tube::q <<"\t"
                << Arteries[10]->Qnew[arteries_10_position]*Tube::q << "\t"
                << Arteries[11]->Qnew[arteries_11_position]*Tube::q <<"\t"
                << Arteries[12]->Qnew[arteries_12_position]*Tube::q <<"\t"
                << Arteries[13]->Qnew[arteries_13_position]*Tube::q <<"\t"
                << Arteries[14]->Qnew[arteries_14_position]*Tube::q <<"\t"
                << Arteries[15]->Qnew[arteries_15_position]*Tube::q <<"\t"
                << Arteries[16]->Qnew[arteries_16_position]*Tube::q <<"\t"
                << Arteries[17]->Qnew[arteries_17_position]*Tube::q <<"\t"
                << Arteries[18]->Qnew[arteries_18_position]*Tube::q <<"\t"
                << Arteries[19]->Qnew[arteries_19_position]*Tube::q <<"\t"
                << Arteries[20]->Qnew[arteries_20_position]*Tube::q << "\t"
                << Arteries[21]->Qnew[arteries_21_position]*Tube::q <<"\t"
                << Arteries[22]->Qnew[arteries_22_position]*Tube::q <<"\t"
                << Arteries[23]->Qnew[arteries_23_position]*Tube::q << endl;
                outdata_flow_different_levels.close();

                double anew_arteries_0_position=Arteries[0]->Anew[arteries_0_position];
                double anew_arteries_1_position=Arteries[1]->Anew[arteries_1_position];
                double anew_arteries_2_position=Arteries[2]->Anew[arteries_2_position];
                double anew_arteries_3_position=Arteries[3]->Anew[arteries_3_position];
                double anew_arteries_4_position=Arteries[4]->Anew[arteries_4_position];
                double anew_arteries_5_position=Arteries[5]->Anew[arteries_5_position];
                double anew_arteries_6_position=Arteries[6]->Anew[arteries_6_position];
                double anew_arteries_7_position=Arteries[7]->Anew[arteries_7_position];
                double anew_arteries_8_position=Arteries[8]->Anew[arteries_8_position];
                double anew_arteries_9_position=Arteries[9]->Anew[arteries_9_position];
                double anew_arteries_10_position=Arteries[10]->Anew[arteries_10_position];
                double anew_arteries_11_position=Arteries[11]->Anew[arteries_11_position];
                double anew_arteries_12_position=Arteries[12]->Anew[arteries_12_position];
                double anew_arteries_13_position=Arteries[13]->Anew[arteries_13_position];
                double anew_arteries_14_position=Arteries[14]->Anew[arteries_14_position];
                double anew_arteries_15_position=Arteries[15]->Anew[arteries_15_position];
                double anew_arteries_16_position=Arteries[16]->Anew[arteries_16_position];
                double anew_arteries_17_position=Arteries[17]->Anew[arteries_17_position];
                double anew_arteries_18_position=Arteries[18]->Anew[arteries_18_position];
                double anew_arteries_19_position=Arteries[19]->Anew[arteries_19_position];
                double anew_arteries_20_position=Arteries[20]->Anew[arteries_20_position];
                double anew_arteries_21_position=Arteries[21]->Anew[arteries_21_position];
                double anew_arteries_22_position=Arteries[22]->Anew[arteries_22_position];
                double anew_arteries_23_position=Arteries[23]->Anew[arteries_23_position];

                double val_to_multiply=Tube::rho*Tube::g*Tube::Lr/Tube::conv;

               ofstream outdata_pressure_different_levels;
               outdata_pressure_different_levels.open("/xlwork1/scott/full_model_nov_22/pressure_different_levels_UPDATED_TREE_RUNS_FOR_PAPER_MECHANICS_60PERCENT.dat", std::ofstream::app);
               outdata_pressure_different_levels << Tube::tstart1 <<"\t"<< Tube::current_position_in_coronary <<"\t"
               << Arteries[0]->P(Tube::current_position_in_coronary,arteries_0_position,anew_arteries_0_position)*val_to_multiply << "\t"                            << Arteries[1]->P(Tube::current_position_in_coronary,arteries_1_position,anew_arteries_1_position)*val_to_multiply << "\t"
               << Arteries[2]->P(Tube::current_position_in_coronary,arteries_2_position,anew_arteries_2_position)*val_to_multiply << "\t"
               << Arteries[3]->P(Tube::current_position_in_coronary,arteries_3_position,anew_arteries_3_position)*val_to_multiply << "\t"
               << Arteries[4]->P(Tube::current_position_in_coronary,arteries_4_position,anew_arteries_4_position)*val_to_multiply << "\t"
               << Arteries[5]->P(Tube::current_position_in_coronary,arteries_5_position,anew_arteries_5_position)*val_to_multiply << "\t"
               << Arteries[6]->P(Tube::current_position_in_coronary,arteries_6_position,anew_arteries_6_position)*val_to_multiply << "\t"
               << Arteries[7]->P(Tube::current_position_in_coronary,arteries_7_position,anew_arteries_7_position)*val_to_multiply << "\t"
               << Arteries[8]->P(Tube::current_position_in_coronary,arteries_8_position,anew_arteries_8_position)*val_to_multiply <<"\t"
               << Arteries[9]->P(Tube::current_position_in_coronary,arteries_9_position,anew_arteries_9_position)*val_to_multiply << "\t"
               << Arteries[10]->P(Tube::current_position_in_coronary,arteries_10_position,anew_arteries_10_position)*val_to_multiply << "\t"                         << Arteries[11]->P(Tube::current_position_in_coronary,arteries_11_position,anew_arteries_11_position)*val_to_multiply << "\t"
               << Arteries[12]->P(Tube::current_position_in_coronary,arteries_12_position,anew_arteries_12_position)*val_to_multiply << "\t"
               << Arteries[13]->P(Tube::current_position_in_coronary,arteries_13_position,anew_arteries_13_position)*val_to_multiply << "\t"
               << Arteries[14]->P(Tube::current_position_in_coronary,arteries_14_position,anew_arteries_14_position)*val_to_multiply << "\t"
               << Arteries[15]->P(Tube::current_position_in_coronary,arteries_15_position,anew_arteries_15_position)*val_to_multiply << "\t"
               << Arteries[16]->P(Tube::current_position_in_coronary,arteries_16_position,anew_arteries_16_position)*val_to_multiply << "\t"
               << Arteries[17]->P(Tube::current_position_in_coronary,arteries_17_position,anew_arteries_17_position)*val_to_multiply << "\t"
               << Arteries[18]->P(Tube::current_position_in_coronary,arteries_18_position,anew_arteries_18_position)*val_to_multiply <<"\t"
               << Arteries[19]->P(Tube::current_position_in_coronary,arteries_19_position,anew_arteries_19_position)*val_to_multiply << "\t"                         << Arteries[20]->P(Tube::current_position_in_coronary,arteries_20_position,anew_arteries_20_position)*val_to_multiply << "\t"                         << Arteries[21]->P(Tube::current_position_in_coronary,arteries_21_position,anew_arteries_21_position)*val_to_multiply << "\t"
               << Arteries[22]->P(Tube::current_position_in_coronary,arteries_22_position,anew_arteries_22_position)*val_to_multiply << "\t"
               << Arteries[23]->P(Tube::current_position_in_coronary,arteries_23_position,anew_arteries_23_position)*val_to_multiply << endl;
               outdata_pressure_different_levels.close();


              }

Tube::current_position_in_coronary=Tube::current_position_in_coronary+1;

if (Tube::current_position_in_coronary == 65537)
{
  Tube::current_position_in_coronary=0;
}

              }

            }

            }


            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
      //      pout << "loading pressure is " << SAMRAI_MPI::maxReduction(BoundaryConditions::P_current_load) << "\n";
            // since CPU 0 may not contain any pressure boundary, the P_current_load may not output properly;
            //SAMRAI_MPI::maxReduction(T_max), works for positive pressure

            //pout << "\nloading pressure beginning: " << BoundaryConditions::loading_pressure(loop_time) << " mmHg" << endl << endl;
            pout << "lump pressure beginning: " << BoundaryConditions::P_old_loading<< " mmHg" << endl;
    	      pout << "LV cavity volume beginning: " << BoundaryConditions::LV_volume<< "  mL"<<endl;
            pout << "cardiac phase in: " << BoundaryConditions::LV_states<<endl;
            BoundaryConditions::P_old_loading = BoundaryConditions::P_current_loading;


            BoundaryConditions::LV_volume_old = BoundaryConditions::LV_volume;
            dt = time_integrator->getMaximumTimeStepSize();


            BoundaryConditions::dt = dt;
            BoundaryConditions::boundary_info = mesh.boundary_info.get();

            time_integrator->advanceHierarchy(dt);


//*******************************************************************************
                    ib_post_processor->postProcessData(loop_time);

                    System& lamda_system = equation_systems->get_system<System>("p_f interpolation system");
                    pout<< "lambda_system number: "<< lamda_system.number()<<"\n";
                    NumericVector<double>* pf_vec = lamda_system.solution.get();
                    NumericVector<double>* pf_ghost_vec = lamda_system.current_local_solution.get();
                    pf_vec->localize(*pf_ghost_vec);

                    System& sigma_dil_system = equation_systems->get_system<System>("sigma_dil reconstruction system");
                    NumericVector<double>* sigma_dil_vec = sigma_dil_system.solution.get();
                    NumericVector<double>* sigma_dil_ghost_vec = sigma_dil_system.current_local_solution.get();
                    sigma_dil_vec->localize(*sigma_dil_ghost_vec);

                    System& sigma_dev_system = equation_systems->get_system<System>("sigma_dev reconstruction system");
                    NumericVector<double>* sigma_dev_vec = sigma_dev_system.solution.get();
                    NumericVector<double>* sigma_dev_ghost_vec = sigma_dev_system.current_local_solution.get();
                    sigma_dev_vec->localize(*sigma_dev_ghost_vec);
//******************************************************************************


if(iteration_num==1)
{

if(my_rank == broadcast_root)
{

      std::ifstream ifsendo("nodes_9_subdomains.csv");

      int current_subdomain_id;
      int IDtemp=1;
      int pIndex = 0;

      while (!ifsendo.eof() & IDtemp<=BoundaryConditions::total_node_num)
      {
       ifsendo>>current_subdomain_id;
       IDtemp++;
       ConoraryFlowInterface::node_current_subdomain_lists[pIndex]=current_subdomain_id;
       pIndex = pIndex + 1;

Tube::number_of_nodes_in_subdomain[current_subdomain_id-1]=Tube::number_of_nodes_in_subdomain[current_subdomain_id-1]+1;

      }


Tube::initial_or_integrated=1;
Tube::current_position_in_coronary=0;

}

}

if(iteration_num==1)
{
ConoraryFlowInterface:: update_id_test(equation_systems);
}
//******************************************************************************



//compute the subdomain volumes in the following.
ConoraryFlowInterface::update_Q_elem(equation_systems,iteration_num);



//assign the sinks and sources appropriately.
if(my_rank == broadcast_root)
{


  for (int q_coro = 0; q_coro <Tube::number_of_subdomains; q_coro++)
  {
    Tube::source_from_coro[q_coro] = Tube::q*(Arteries[Tube::source_number[q_coro]]->Qnew[Arteries[Tube::source_number[q_coro]]->N]);
    Tube::sink_from_coro[q_coro] = Tube::q*(Arteries[Tube::sink_number[q_coro]]->Qnew[Arteries[Tube::sink_number[q_coro]]->N]);

    Tube::net_from_coro[q_coro]=Tube::source_from_coro[q_coro]+Tube::sink_from_coro[q_coro];
    Tube::net_from_coro_averaged[q_coro]=(Tube::net_from_coro[q_coro])/(ConoraryFlowInterface::total_subdomain_volumes[q_coro]);

    Tube::net_from_coro_in_dimensional_space[q_coro]=Tube::net_from_coro[q_coro]*((Tube::onePeriodDt)/(dt));
    Tube::net_from_coro_averaged_in_dimensional_space[q_coro]=(Tube::net_from_coro_in_dimensional_space[q_coro])/(ConoraryFlowInterface::total_subdomain_volumes[q_coro]);
  }

  Tube::flow_added_this_timestep=(Tube::net_from_coro[0]+Tube::net_from_coro[1]+Tube::net_from_coro[2]+Tube::net_from_coro[3]+Tube::net_from_coro[4]+Tube::net_from_coro[5]+Tube::net_from_coro[6]+Tube::net_from_coro[7]+Tube::net_from_coro[8])*(Tube::onePeriodDt);
  Tube::total_flow_added=Tube::total_flow_added+Tube::flow_added_this_timestep;

  Tube::flow_added_this_timestep_in_dimensional_space=(Tube::net_from_coro_in_dimensional_space[0]+Tube::net_from_coro_in_dimensional_space[1]+Tube::net_from_coro_in_dimensional_space[2]+Tube::net_from_coro_in_dimensional_space[3]+Tube::net_from_coro_in_dimensional_space[4]+Tube::net_from_coro_in_dimensional_space[5]+Tube::net_from_coro_in_dimensional_space[6]+Tube::net_from_coro_in_dimensional_space[7]+Tube::net_from_coro_in_dimensional_space[8])*(dt);
  Tube::total_flow_added_in_dimensional_space=Tube::total_flow_added_in_dimensional_space+Tube::flow_added_this_timestep_in_dimensional_space;

}


MPI_Bcast(&Tube::net_from_coro_averaged_in_dimensional_space[0], 1, MPI_DOUBLE, broadcast_root, MPI_COMM_WORLD);
MPI_Bcast(&Tube::net_from_coro_averaged_in_dimensional_space[1], 1, MPI_DOUBLE, broadcast_root, MPI_COMM_WORLD);
MPI_Bcast(&Tube::net_from_coro_averaged_in_dimensional_space[2], 1, MPI_DOUBLE, broadcast_root, MPI_COMM_WORLD);
MPI_Bcast(&Tube::net_from_coro_averaged_in_dimensional_space[3], 1, MPI_DOUBLE, broadcast_root, MPI_COMM_WORLD);
MPI_Bcast(&Tube::net_from_coro_averaged_in_dimensional_space[4], 1, MPI_DOUBLE, broadcast_root, MPI_COMM_WORLD);
MPI_Bcast(&Tube::net_from_coro_averaged_in_dimensional_space[5], 1, MPI_DOUBLE, broadcast_root, MPI_COMM_WORLD);
MPI_Bcast(&Tube::net_from_coro_averaged_in_dimensional_space[6], 1, MPI_DOUBLE, broadcast_root, MPI_COMM_WORLD);
MPI_Bcast(&Tube::net_from_coro_averaged_in_dimensional_space[7], 1, MPI_DOUBLE, broadcast_root, MPI_COMM_WORLD);
MPI_Bcast(&Tube::net_from_coro_averaged_in_dimensional_space[8], 1, MPI_DOUBLE, broadcast_root, MPI_COMM_WORLD);
MPI_Barrier(init.comm().get());



          sourceToDarcy::solve_std(equation_systems);

//******************************************************************************

          DarcyPressure::solve_pstar(equation_systems,loop_time,dt);
          DarcyMass::solve_mstar(equation_systems,loop_time,dt);

          MechanicsModel::update_J_system(equation_systems, dt);
          MechanicsModel::update_FF_system(equation_systems);

          DarcyPressure::solve_p(equation_systems,loop_time,dt);
          DarcyVelocity::solve_wx(equation_systems);
          DarcyVelocity::solve_wy(equation_systems);
          DarcyVelocity::solve_wz(equation_systems);

          DarcyMass::solve_m(equation_systems,loop_time,dt);
          DarcyMass::solve_maverage(equation_systems,loop_time,dt);

          sourceToIB::solve_Sstar(equation_systems,loop_time,dt);




//******************************************************************************


AdvanceCoronaryCode::advance_variables(equation_systems);

if (BoundaryConditions::LV_states == 5 && ActiveContraction::time_to_enter == 0)
{
  ActiveContraction::save_reduction_value(equation_systems, loop_time, dt);
  ActiveContraction::time_to_enter=1;
}

if (BoundaryConditions::LV_states == 1)
{
  ActiveContraction::time_to_enter=0;
}



if (MechanicsModel::enable_active_tension)
{
    ActiveContraction::update_active_tension_model_state_variables(equation_systems, loop_time, dt);
}



//calculate the LV cavity volume for pressure update

   BoundaryConditions::updatePointsPositionGeneral(equation_systems,
                                BoundaryConditions::LV_endo_points_list,
                                BoundaryConditions::LV_endo_points,
                                BoundaryConditions::LV_all_points,
                                BoundaryConditions::LV_NoOfEndoNode);


   BoundaryConditions::solve_wall_volume(equation_systems);



   if(my_rank == broadcast_root)
   {
   BoundaryConditions::whole_volume = tetVolumeCalculationByPoints(BoundaryConditions::LV_all_points);
   BoundaryConditions::LV_volume=BoundaryConditions::whole_volume-BoundaryConditions::wall_volume;
   }
   MPI_Bcast(&BoundaryConditions::whole_volume, 1, MPI_DOUBLE, broadcast_root, MPI_COMM_WORLD);
   MPI_Bcast(&BoundaryConditions::LV_volume, 1, MPI_DOUBLE, broadcast_root, MPI_COMM_WORLD);
   MPI_Barrier(init.comm().get());

    BoundaryConditions::phase_determination(loop_time); //determine what phase we are in
    BoundaryConditions::lump_para_pressure(equation_systems,loop_time);  //update the pressure


    //for checking purposes
    BoundaryConditions::print_info(loop_time);


//---------------------------------------------------------------------------------------------





if(my_rank == broadcast_root)
{

for(int qq=0; qq<9; qq++)
{
  Tube::pressure_to_coronary_current[qq]=0.0;
}

}

// find the new external pressure to be passed to the coronary code

ConoraryFlowInterface::update_Pext(equation_systems,loop_time);




if(my_rank == broadcast_root)
{



  if(BoundaryConditions::current_cycle==2 && Tube::second_period_exp_update==0)
  {
  for(int hhh=0; hhh<65537; hhh++)
  {
Tube::imp_holder_1[hhh]=Tube::imp_holder_1_period[hhh];
Tube::imp_holder_2[hhh]=Tube::imp_holder_2_period[hhh];
Tube::imp_holder_3[hhh]=Tube::imp_holder_3_period[hhh];
Tube::imp_holder_4[hhh]=Tube::imp_holder_4_period[hhh];
Tube::imp_holder_5[hhh]=Tube::imp_holder_5_period[hhh];
Tube::imp_holder_6[hhh]=Tube::imp_holder_6_period[hhh];
Tube::imp_holder_7[hhh]=Tube::imp_holder_7_period[hhh];
Tube::imp_holder_8[hhh]=Tube::imp_holder_8_period[hhh];
Tube::imp_holder_9[hhh]=Tube::imp_holder_9_period[hhh];
  }
 Tube::second_period_exp_update=1;
  }




for(int qq=0; qq<9; qq++)
{
Tube::pressure_to_coronary_current[qq]=(Tube::dyne_to_mmhg)*(Tube::pressure_to_coronary_current[qq]/(ConoraryFlowInterface::total_subdomain_volumes[qq]));
}


if(loop_time < 0.001)
{
  for(int qq=0; qq<9; qq++)
  {
  Tube::pressure_to_coronary_current[qq]=0.0;
  }
}


Tube::imp_holder_1_period[Tube::exp_position]=Tube::pressure_to_coronary_current[0];
Tube::imp_holder_2_period[Tube::exp_position]=Tube::pressure_to_coronary_current[1];
Tube::imp_holder_3_period[Tube::exp_position]=Tube::pressure_to_coronary_current[2];
Tube::imp_holder_4_period[Tube::exp_position]=Tube::pressure_to_coronary_current[3];
Tube::imp_holder_5_period[Tube::exp_position]=Tube::pressure_to_coronary_current[4];
Tube::imp_holder_6_period[Tube::exp_position]=Tube::pressure_to_coronary_current[5];
Tube::imp_holder_7_period[Tube::exp_position]=Tube::pressure_to_coronary_current[6];
Tube::imp_holder_8_period[Tube::exp_position]=Tube::pressure_to_coronary_current[7];
Tube::imp_holder_9_period[Tube::exp_position]=Tube::pressure_to_coronary_current[8];

  ofstream outdata_imp_holder_period;
  if (Tube::counter_sec % 1 ==0)
  {
  outdata_imp_holder_period.open("/xlwork1/scott/full_model_nov_22/imp_holder_period_UPDATED_TREE_RUNS_FOR_PAPER_MECHANICS_60PERCENT.dat", std::ofstream::app); // opens the file
  outdata_imp_holder_period <<
  loop_time << "\t" <<
  Tube::exp_position << "\t" <<
  Tube::imp_holder_1_period[Tube::exp_position] << "\t" <<
  Tube::imp_holder_2_period[Tube::exp_position] << "\t" <<
  Tube::imp_holder_3_period[Tube::exp_position] << "\t" <<
  Tube::imp_holder_4_period[Tube::exp_position] << "\t" <<
  Tube::imp_holder_5_period[Tube::exp_position] << "\t" <<
  Tube::imp_holder_6_period[Tube::exp_position] << "\t" <<
  Tube::imp_holder_7_period[Tube::exp_position] << "\t" <<
  Tube::imp_holder_8_period[Tube::exp_position] << "\t" <<
  Tube::imp_holder_9_period[Tube::exp_position] << endl;
  outdata_imp_holder_period.close();
  }


  if(BoundaryConditions::current_cycle>2)
{

Tube::imp_holder_1[Tube::exp_position]=Tube::imp_holder_1_period[Tube::exp_position];
Tube::imp_holder_2[Tube::exp_position]=Tube::imp_holder_2_period[Tube::exp_position];
Tube::imp_holder_3[Tube::exp_position]=Tube::imp_holder_3_period[Tube::exp_position];
Tube::imp_holder_4[Tube::exp_position]=Tube::imp_holder_4_period[Tube::exp_position];
Tube::imp_holder_5[Tube::exp_position]=Tube::imp_holder_5_period[Tube::exp_position];
Tube::imp_holder_6[Tube::exp_position]=Tube::imp_holder_6_period[Tube::exp_position];
Tube::imp_holder_7[Tube::exp_position]=Tube::imp_holder_7_period[Tube::exp_position];
Tube::imp_holder_8[Tube::exp_position]=Tube::imp_holder_8_period[Tube::exp_position];
Tube::imp_holder_9[Tube::exp_position]=Tube::imp_holder_9_period[Tube::exp_position];

}


Tube::exp_position=Tube::exp_position+1;

if(Tube::exp_position == 65537)
{
Tube::exp_position=0;
}




}





//---------------------------------------------------------------------------------------------
//now run the coronary model again


    if(my_rank == broadcast_root)
    {

    solver (Arteries, Tube::tstart1, Tube::tend1, Tube::onePeriodDt, 0.0);

    Tube::tstart1 = Tube::tend1;
    Tube::tend1 = Tube::tstart1 + Tube::onePeriodDt;
    Tube::counter_sec=Tube::counter_sec+1;

    ofstream outdata_coro;
    if (Tube::counter_sec % 100 ==0)
    {
    outdata_coro.open("/xlwork1/scott/full_model_nov_22/coronary_data_coupled_run_UPDATED_TREE_RUNS_FOR_PAPER_MECHANICS_60PERCENT.dat", std::ofstream::app); // opens the file
      for (int i=0; i < Tube::nbrves_passer; i ++)
     {
      outdata_coro << Tube::tstart1 <<"\t"<< i <<"\t"<< Arteries[i]->Qnew[Arteries[i]->N]*Tube::q << endl;
     }
      outdata_coro.close();
    }

      ofstream outdata_pressures;
      if (Tube::counter_sec % 100 ==0)
      {
      outdata_pressures.open("/xlwork1/scott/full_model_nov_22/pressure_data_coupled_run_UPDATED_TREE_RUNS_FOR_PAPER_MECHANICS_60PERCENT.dat", std::ofstream::app); // opens the file
      outdata_pressures << loop_time << "\t"<< BoundaryConditions::P_Ao <<"\t"<< BoundaryConditions::Pwk_AV <<"\t"<< BoundaryConditions::P_current_loading  << endl;
      outdata_pressures.close();
      }


      ofstream outdata_external_pressures;
      if (Tube::counter_sec % 100 ==0)
      {
      outdata_external_pressures.open("/xlwork1/scott/full_model_nov_22/external_pressure_UPDATED_TREE_RUNS_FOR_PAPER_MECHANICS_60PERCENT.dat", std::ofstream::app); // opens the file
        for (int i=0; i < Tube::number_of_subdomains; i ++)
       {
      outdata_external_pressures << Tube::tstart1 <<"\t"<< i <<"\t"<< Tube::pressure_to_coronary_current[i] <<"\t"<< Tube::net_from_coro[i] << "\t" << Tube::net_from_coro_averaged[i] << "\t" << Tube::source_from_coro[i] << "\t" << Tube::sink_from_coro[i] << endl;
       }
      outdata_external_pressures.close();
       }


      ofstream outdata_volume_data;
      if (Tube::counter_sec % 100 ==0)
      {
      outdata_volume_data.open("/xlwork1/scott/full_model_nov_22/volume_data_UPDATED_TREE_RUNS_FOR_PAPER_MECHANICS_60PERCENT.dat", std::ofstream::app); // opens the file
      outdata_volume_data << loop_time << "\t"<< BoundaryConditions::whole_volume <<"\t"<< BoundaryConditions::wall_volume <<"\t"<< BoundaryConditions::LV_volume <<"\t"<< Tube::total_flow_added << endl;
      outdata_volume_data.close();
      }

      ofstream outdata_average_data;
      if (Tube::counter_sec % 100 ==0)
      {
      outdata_average_data.open("/xlwork1/scott/full_model_nov_22/average_data_UPDATED_TREE_RUNS_FOR_PAPER_MECHANICS_60PERCENT.dat", std::ofstream::app); // opens the file
      outdata_average_data << loop_time << "\t"<< Tube::flow_added_this_timestep <<"\t"<< Tube::summed_source <<"\t"<< Tube::total_flow_added <<"\t"<< Tube::summed_source_all_time <<"\t"<< Tube::added_mass_sum << "\t" << Tube::average_added_mass_sum << "\t" << Tube::flow_added_this_timestep_in_dimensional_space << "\t" << Tube::total_flow_added_in_dimensional_space << endl;
      outdata_average_data.close();
      }


                int arteries_0_position=floor(0.5*Arteries[0]->N);
                int arteries_1_position=floor(0.5*Arteries[1]->N);
                int arteries_2_position=floor(0.5*Arteries[2]->N);
                int arteries_3_position=floor(0.5*Arteries[3]->N);
                int arteries_4_position=floor(0.5*Arteries[4]->N);
                int arteries_5_position=floor(0.5*Arteries[5]->N);
                int arteries_6_position=floor(0.5*Arteries[6]->N);
                int arteries_7_position=floor(0.5*Arteries[7]->N);
                int arteries_8_position=floor(0.5*Arteries[8]->N);
                int arteries_9_position=floor(0.5*Arteries[9]->N);
                int arteries_10_position=floor(0.5*Arteries[10]->N);
                int arteries_11_position=floor(0.5*Arteries[11]->N);
                int arteries_12_position=floor(0.5*Arteries[12]->N);
                int arteries_13_position=floor(0.5*Arteries[13]->N);
                int arteries_14_position=floor(0.5*Arteries[14]->N);
                int arteries_15_position=floor(0.5*Arteries[15]->N);
                int arteries_16_position=floor(0.5*Arteries[16]->N);
                int arteries_17_position=floor(0.5*Arteries[17]->N);
                int arteries_18_position=floor(0.5*Arteries[18]->N);
                int arteries_19_position=floor(0.5*Arteries[19]->N);
                int arteries_20_position=floor(0.5*Arteries[20]->N);
                int arteries_21_position=floor(0.5*Arteries[21]->N);
                int arteries_22_position=floor(0.5*Arteries[22]->N);
                int arteries_23_position=floor(0.5*Arteries[23]->N);

                ofstream outdata_coupled_flow_different_levels;
      outdata_coupled_flow_different_levels.open("/xlwork1/scott/full_model_nov_22/coupled_flow_different_levels_UPDATED_TREE_RUNS_FOR_PAPER_MECHANICS_60PERCENT.dat", std::ofstream::app);
      outdata_coupled_flow_different_levels<< Tube::tstart1 <<"\t" << Tube::current_position_in_coronary <<"\t"
                << Arteries[0]->Qnew[arteries_0_position]*Tube::q << "\t"
                << Arteries[1]->Qnew[arteries_1_position]*Tube::q <<"\t"
                << Arteries[2]->Qnew[arteries_2_position]*Tube::q <<"\t"
                << Arteries[3]->Qnew[arteries_3_position]*Tube::q <<"\t"
                << Arteries[4]->Qnew[arteries_4_position]*Tube::q <<"\t"
                << Arteries[5]->Qnew[arteries_5_position]*Tube::q <<"\t"
                << Arteries[6]->Qnew[arteries_6_position]*Tube::q <<"\t"
                << Arteries[7]->Qnew[arteries_7_position]*Tube::q <<"\t"
                << Arteries[8]->Qnew[arteries_8_position]*Tube::q <<"\t"
                << Arteries[9]->Qnew[arteries_9_position]*Tube::q <<"\t"
                << Arteries[10]->Qnew[arteries_10_position]*Tube::q << "\t"
                << Arteries[11]->Qnew[arteries_11_position]*Tube::q <<"\t"
                << Arteries[12]->Qnew[arteries_12_position]*Tube::q <<"\t"
                << Arteries[13]->Qnew[arteries_13_position]*Tube::q <<"\t"
                << Arteries[14]->Qnew[arteries_14_position]*Tube::q <<"\t"
                << Arteries[15]->Qnew[arteries_15_position]*Tube::q <<"\t"
                << Arteries[16]->Qnew[arteries_16_position]*Tube::q <<"\t"
                << Arteries[17]->Qnew[arteries_17_position]*Tube::q <<"\t"
                << Arteries[18]->Qnew[arteries_18_position]*Tube::q <<"\t"
                << Arteries[19]->Qnew[arteries_19_position]*Tube::q <<"\t"
                << Arteries[20]->Qnew[arteries_20_position]*Tube::q << "\t"
                << Arteries[21]->Qnew[arteries_21_position]*Tube::q <<"\t"
                << Arteries[22]->Qnew[arteries_22_position]*Tube::q <<"\t"
                << Arteries[23]->Qnew[arteries_23_position]*Tube::q << endl;
                outdata_coupled_flow_different_levels.close();

                double anew_arteries_0_position=Arteries[0]->Anew[arteries_0_position];
                double anew_arteries_1_position=Arteries[1]->Anew[arteries_1_position];
                double anew_arteries_2_position=Arteries[2]->Anew[arteries_2_position];
                double anew_arteries_3_position=Arteries[3]->Anew[arteries_3_position];
                double anew_arteries_4_position=Arteries[4]->Anew[arteries_4_position];
                double anew_arteries_5_position=Arteries[5]->Anew[arteries_5_position];
                double anew_arteries_6_position=Arteries[6]->Anew[arteries_6_position];
                double anew_arteries_7_position=Arteries[7]->Anew[arteries_7_position];
                double anew_arteries_8_position=Arteries[8]->Anew[arteries_8_position];
                double anew_arteries_9_position=Arteries[9]->Anew[arteries_9_position];
                double anew_arteries_10_position=Arteries[10]->Anew[arteries_10_position];
                double anew_arteries_11_position=Arteries[11]->Anew[arteries_11_position];
                double anew_arteries_12_position=Arteries[12]->Anew[arteries_12_position];
                double anew_arteries_13_position=Arteries[13]->Anew[arteries_13_position];
                double anew_arteries_14_position=Arteries[14]->Anew[arteries_14_position];
                double anew_arteries_15_position=Arteries[15]->Anew[arteries_15_position];
                double anew_arteries_16_position=Arteries[16]->Anew[arteries_16_position];
                double anew_arteries_17_position=Arteries[17]->Anew[arteries_17_position];
                double anew_arteries_18_position=Arteries[18]->Anew[arteries_18_position];
                double anew_arteries_19_position=Arteries[19]->Anew[arteries_19_position];
                double anew_arteries_20_position=Arteries[20]->Anew[arteries_20_position];
                double anew_arteries_21_position=Arteries[21]->Anew[arteries_21_position];
                double anew_arteries_22_position=Arteries[22]->Anew[arteries_22_position];
                double anew_arteries_23_position=Arteries[23]->Anew[arteries_23_position];

                double val_to_multiply=Tube::rho*Tube::g*Tube::Lr/Tube::conv;

                ofstream outdata_coupled_pressure_different_levels;
     outdata_coupled_pressure_different_levels.open("/xlwork1/scott/full_model_nov_22/coupled_pressure_different_levels_UPDATED_TREE_RUNS_FOR_PAPER_MECHANICS_60PERCENT.dat", std::ofstream::app);
     outdata_coupled_pressure_different_levels << Tube::tstart1 <<"\t"<< Tube::current_position_in_coronary <<"\t"
               << Arteries[0]->P(Tube::current_position_in_coronary,arteries_0_position,anew_arteries_0_position)*val_to_multiply << "\t"                            << Arteries[1]->P(Tube::current_position_in_coronary,arteries_1_position,anew_arteries_1_position)*val_to_multiply << "\t"
               << Arteries[2]->P(Tube::current_position_in_coronary,arteries_2_position,anew_arteries_2_position)*val_to_multiply << "\t"
               << Arteries[3]->P(Tube::current_position_in_coronary,arteries_3_position,anew_arteries_3_position)*val_to_multiply << "\t"
               << Arteries[4]->P(Tube::current_position_in_coronary,arteries_4_position,anew_arteries_4_position)*val_to_multiply << "\t"
               << Arteries[5]->P(Tube::current_position_in_coronary,arteries_5_position,anew_arteries_5_position)*val_to_multiply << "\t"
               << Arteries[6]->P(Tube::current_position_in_coronary,arteries_6_position,anew_arteries_6_position)*val_to_multiply << "\t"
               << Arteries[7]->P(Tube::current_position_in_coronary,arteries_7_position,anew_arteries_7_position)*val_to_multiply << "\t"
               << Arteries[8]->P(Tube::current_position_in_coronary,arteries_8_position,anew_arteries_8_position)*val_to_multiply <<"\t"
               << Arteries[9]->P(Tube::current_position_in_coronary,arteries_9_position,anew_arteries_9_position)*val_to_multiply << "\t"
               << Arteries[10]->P(Tube::current_position_in_coronary,arteries_10_position,anew_arteries_10_position)*val_to_multiply << "\t"                         << Arteries[11]->P(Tube::current_position_in_coronary,arteries_11_position,anew_arteries_11_position)*val_to_multiply << "\t"
               << Arteries[12]->P(Tube::current_position_in_coronary,arteries_12_position,anew_arteries_12_position)*val_to_multiply << "\t"
               << Arteries[13]->P(Tube::current_position_in_coronary,arteries_13_position,anew_arteries_13_position)*val_to_multiply << "\t"
               << Arteries[14]->P(Tube::current_position_in_coronary,arteries_14_position,anew_arteries_14_position)*val_to_multiply << "\t"
               << Arteries[15]->P(Tube::current_position_in_coronary,arteries_15_position,anew_arteries_15_position)*val_to_multiply << "\t"
               << Arteries[16]->P(Tube::current_position_in_coronary,arteries_16_position,anew_arteries_16_position)*val_to_multiply << "\t"
               << Arteries[17]->P(Tube::current_position_in_coronary,arteries_17_position,anew_arteries_17_position)*val_to_multiply << "\t"
               << Arteries[18]->P(Tube::current_position_in_coronary,arteries_18_position,anew_arteries_18_position)*val_to_multiply <<"\t"
               << Arteries[19]->P(Tube::current_position_in_coronary,arteries_19_position,anew_arteries_19_position)*val_to_multiply << "\t"                         << Arteries[20]->P(Tube::current_position_in_coronary,arteries_20_position,anew_arteries_20_position)*val_to_multiply << "\t"                         << Arteries[21]->P(Tube::current_position_in_coronary,arteries_21_position,anew_arteries_21_position)*val_to_multiply << "\t"
               << Arteries[22]->P(Tube::current_position_in_coronary,arteries_22_position,anew_arteries_22_position)*val_to_multiply << "\t"
               << Arteries[23]->P(Tube::current_position_in_coronary,arteries_23_position,anew_arteries_23_position)*val_to_multiply << endl;
     outdata_coupled_pressure_different_levels.close();





      Tube::current_position_in_coronary=Tube::current_position_in_coronary+1;

      if (Tube::current_position_in_coronary == 65537)
      {
        Tube::current_position_in_coronary=0;
      }




    }




//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------


            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time in the ventricle model is " << loop_time << "\n";
            pout << "Simulation time in the coronary model (dimensionless) is " << Tube::tstart1 << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";


            loop_time += dt;


            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
    //        iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data  && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                if (uses_visit)
                {
                    time_integrator->setupPlotData();
                    visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                }
                if (uses_exodus)
                {
                    ib_post_processor->postProcessData(loop_time);
                    exodus_io->write_timestep(exodus_filename, *equation_systems,
                                              iteration_num/viz_dump_interval+1, loop_time);
                }
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
                ib_method_ops->writeFEDataToRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
            if (dump_postproc_data && (iteration_num % postproc_data_dump_interval == 0 || last_step))
            {
                output_data(patch_hierarchy, time_integrator, iteration_num, loop_time, postproc_data_dump_dirname);
            }
        }

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    //PetscFinalize(); //seems there is no need to use PetscFinalize() since using LibMeshInit to manage
    return true;
} // run_example





void
output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
            Pointer<IBHierarchyIntegrator> time_integrator,
            const int iteration_num,
            const double loop_time,
            const string& data_dump_dirname)
{
    plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    plog << "simulation time is " << loop_time << endl;

    // Write Cartesian data.
    string file_name = data_dump_dirname + "/" + "hier_data.";
    char temp_buf[128];
    sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, SAMRAI_MPI::getRank());
    file_name += temp_buf;
    Pointer<HDFDatabase> hier_db = new HDFDatabase("hier_db");
    hier_db->create(file_name);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    ComponentSelector hier_data;
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(time_integrator->getVelocityVariable(),
                                                           time_integrator->getCurrentContext()));
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(time_integrator->getPressureVariable(),
                                                           time_integrator->getCurrentContext()));
    patch_hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
    hier_db->putDouble("loop_time", loop_time);
    hier_db->putInteger("iteration_num", iteration_num);
    hier_db->close();
    return;
} // output_data




















int main(int argc, char** argv)
{

  verbosity=1;

	//Initialisation for the BV codes

  LibMeshInit init(argc, argv); //initialize libMesh, which will initlize Petsc

  SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
  SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
  SAMRAIManager::startup();

  // Increase maximum patch data component indices
  SAMRAIManager::setMaxNumberPatchDataEntries(2500);


      Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "INS.log");
      Pointer<Database> input_db = app_initializer->getInputDatabase();




Tube::exp_reducer=0;
Tube::counter_sec=0;
Tube::exp_position=0;

Tube::plts=input_db->getInteger("PLTS");
Tube::tmstps=input_db->getInteger("TMSTPS");
Tube::time_step_discrepency=input_db->getInteger("TIMESTEPDISCREPENCY");
Tube::number_of_subdomains=input_db->getInteger("NUMBER_OF_SUBDOMAINS");

Tube::rho=input_db->getDouble("RHO_CNM");
Tube::conv=input_db->getDouble("CONV_CNM");
Tube::mu1=input_db->getDouble("MU1_CNM");
Tube::mu_pl=input_db->getDouble("MU_PL_CNM");
Tube::Tper=input_db->getDouble("TPER_CNM");
Tube::nu=input_db->getDouble("NU_CNM");
Tube::Fcst=input_db->getDouble("FCST");
Tube::Lr=input_db->getDouble("LR");
Tube::Lr2=input_db->getDouble("LR2");
Tube::Lr3=input_db->getDouble("LR3");
Tube::g=input_db->getDouble("GG");
Tube::Re=input_db->getDouble("RE");
Tube::q=input_db->getDouble("QQ");
Tube::Fr2=input_db->getDouble("FR2");
Tube::Period=input_db->getDouble("PERIOD");
Tube::p0=input_db->getDouble("P0");
Tube::Deltat=input_db->getDouble("DELTAT");
Tube::dyne_to_mmhg=input_db->getDouble("DYNE_TO_MMHG");

Tube::tstart1=input_db->getDouble("TSTART1");
Tube::tend1=input_db->getDouble("TEND1");
Tube::finaltime1=input_db->getDouble("FINALTIME1");
Tube::onePeriodDt=input_db->getDouble("ONEPERIODDT");


for(int qq=0; qq<65537; qq++)
{
  Tube::imp_holder_1[qq]=0.0;
  Tube::imp_holder_2[qq]=0.0;
  Tube::imp_holder_3[qq]=0.0;
  Tube::imp_holder_4[qq]=0.0;
  Tube::imp_holder_5[qq]=0.0;
  Tube::imp_holder_6[qq]=0.0;
  Tube::imp_holder_7[qq]=0.0;
  Tube::imp_holder_8[qq]=0.0;
  Tube::imp_holder_9[qq]=0.0;

  Tube::imp_holder_1_period[qq]=0.0;
  Tube::imp_holder_2_period[qq]=0.0;
  Tube::imp_holder_3_period[qq]=0.0;
  Tube::imp_holder_4_period[qq]=0.0;
  Tube::imp_holder_5_period[qq]=0.0;
  Tube::imp_holder_6_period[qq]=0.0;
  Tube::imp_holder_7_period[qq]=0.0;
  Tube::imp_holder_8_period[qq]=0.0;
  Tube::imp_holder_9_period[qq]=0.0;
}



for(int qq=0; qq<Tube::number_of_subdomains; qq++)
{
  Tube::source_from_coro[qq]=0.0;
  Tube::sink_from_coro[qq]=0.0;
  Tube::net_from_coro[qq]=0.0;
  Tube::net_from_coro_averaged[qq]=0.0;

  Tube::net_from_coro_in_dimensional_space[qq]=0.0;
  Tube::net_from_coro_averaged_in_dimensional_space[qq]=0.0;

  Tube::pressure_to_coronary_current[qq]=0.0;
  Tube::pressure_to_coronary_previous[qq]=0.0;
  Tube::pressure_to_coronary_reduction[qq]=0.0;

  Tube::source_number[qq]=0;
  Tube::sink_number[qq]=0;
  Tube::number_of_nodes_in_subdomain[qq]=0;
  Tube::number_of_elements_in_subdomain[qq]=0;

}

Tube::flow_added_this_timestep=0.0;
Tube::total_flow_added=0.0;
Tube::flow_added_this_timestep_in_dimensional_space=0.0;
Tube::total_flow_added_in_dimensional_space=0.0;
Tube::summed_source=0.0;
Tube::summed_source_all_time=0.0;
Tube::added_mass_sum=0.0;
Tube::average_added_mass_sum=0.0;

Tube::transfer_aortic_previous=80.0;
Tube::transfer_aortic=80.0;
Tube::loop_time_arteries=0.0;
Tube::pass_val=0.0;
Tube::initial_or_integrated=0;
Tube::current_or_previous=0;
Tube::nbrves_passer=0;
Tube::current_position_in_coronary=0;

Tube::total_source=0.0;
Tube::total_sink=0.0;
Tube::net_sink=0.0;

Tube::second_period_exp_update=0;
Tube::conv_compute=0;
Tube::beat_tracker=0;
Tube::time_holder=0.0;


rm = 0.002;
f1 = 0; f2 = 1; f3 = 2530000;
fa1 = 0; fa2 = 1; fa3 = 1030000;

//fa1 = 1.99925e09; fa2 = -23.767; fa3 = 1030000;

asym = 0.41; expo = 2.76; lrrA = 36.0;
lrrV = 36.0;
numHeartBeats = 8;
id = 1;
SVPA = 0;
LVPA = 0;
gridPoints=10;

fv1 = fa1;
fv2 = fa2;
fv3 = fa3;


// Parameters required to initiate class Tube (Length,topradius,botradius,LeftDaughter,RightDaughter,rmin, points,
                                               //init,K,f1,f2,f3,fa1,fa2,fa3,fv1,fv2,fv3,asym,expo,lrrA,lrrV);

char str[20]; int char_count = 0; FILE * treefile;
treefile = fopen("/xlwork1/scott/full_model_nov_22/NewSumTab/scott_tree.tsv", "r");

if (treefile) {while (fscanf(treefile, "%s", str) != EOF) { char_count++; }}; fclose(treefile);
int line_len = 9; nbrves = int(char_count / line_len);
Tube * Arteries[nbrves];

Tube::nbrves_passer=nbrves;

double * array[nbrves];
for (int i = 0; i < nbrves; i++) {array[i] = new double[line_len];}

for (int row = 0; row < nbrves; row++) {
    for (int col = 0; col < line_len; col++) {
        array[row][col] = 0.0;
    }
}

treefile = fopen("/xlwork1/scott/full_model_nov_22/NewSumTab/scott_tree.tsv", "r");

for (int i = 0; i < char_count; i++) {
    int row = float(i) / float(line_len), col = i % line_len;
    fscanf(treefile, "%s", str); array[row][col] = atof(str);
}
fclose(treefile); printf("loaded tube data from tree file\n");

for (int ii = 0; ii < nbrves; ii++) {


    for (int j = 0; j <= 6; j++)
    {
        if (int(array[ii][j]) == array[ii][j])
        {pout << int(array[ii][j]) << "\t";}
        else    {pout << array[ii][j] << "\t";}
    }




double rm_val;
if (int(array[ii][4]) != 0){rm_val =0;}else {rm_val = rm;}
double init_val;
if (int(array[ii][4]) == 0 && int(array[ii][6]) != 0){init_val =3;} else{init_val =  array[ii][7];}

//  in order, handle trifs
if (int(array[ii][4]) != 0 && int(array[ii][5]) != 0 && int(array[ii][6])!=0)
{
    Arteries[int(array[ii][0])] = new Tube(array[ii][1], array[ii][2], array[ii][3], Arteries[int(array[ii][4])],Arteries[int(array[ii][5])], Arteries[int(array[ii][6])],rm_val,  gridPoints,  init_val, f1, f2, f3, fa1, fa2, fa3, fv1, fv2, fv3, asym, expo, lrrA, lrrV, array[ii][8], LVPA);
}
else if (int(array[ii][4]) != 0 && int(array[ii][5]) == 0 && int(array[ii][6])!=0)
{
    Arteries[int(array[ii][0])] = new Tube(array[ii][1], array[ii][2], array[ii][3], Arteries[int(array[ii][4])],0, Arteries[int(array[ii][6])],rm_val,  gridPoints,  init_val, f1, f2, f3, fa1, fa2, fa3, fv1, fv2, fv3, asym, expo, lrrA, lrrV, array[ii][8], LVPA);
}
else if (int(array[ii][4]) != 0 && int(array[ii][5]) == 0 && int(array[ii][6])==0)
{
    Arteries[int(array[ii][0])] = new Tube(array[ii][1], array[ii][2], array[ii][3], Arteries[int(array[ii][4])],0, 0,rm_val,  gridPoints,  init_val, f1, f2, f3, fa1, fa2, fa3, fv1, fv2, fv3, asym, expo, lrrA, lrrV, array[ii][8], LVPA);
}
else if (int(array[ii][4]) == 0 && int(array[ii][5]) == 0 && int(array[ii][6])!=0)
{
    Arteries[int(array[ii][0])] = new Tube(array[ii][1], array[ii][2], array[ii][3], 0,0, Arteries[int(array[ii][6])],rm_val,  gridPoints,  init_val, f1, f2, f3, fa1, fa2, fa3, fv1, fv2, fv3, asym, expo, lrrA, lrrV, array[ii][8], LVPA);
}
else if (int(array[ii][4]) == 0 && int(array[ii][5]) == 0 && int(array[ii][6])==0)
{
    Arteries[int(array[ii][0])] = new Tube(array[ii][1], array[ii][2], array[ii][3], 0,0, 0,rm_val,  gridPoints,  init_val, f1, f2, f3, fa1, fa2, fa3, fv1, fv2, fv3, asym, expo, lrrA, lrrV, array[ii][8], LVPA);
}
else
{
  printf("didn't enter anything \n");
}

// Arteries[int(array[ii][0])] = new Tube(array[ii][1], array[ii][2], array[ii][3], Arteries[int(array[ii][4])], Arteries[int(array[ii][5])], Arteries[int(array[ii][6])],rm_val,  gridPoints,  init_val, f1, f2, f3, fa1, fa2, fa3, fv1, fv2, fv3, asym, expo, lrrA, lrrV, SVPA, LVPA);


Tube::source_number[0]=4;
Tube::source_number[1]=6;
Tube::source_number[2]=7;
Tube::source_number[3]=10;
Tube::source_number[4]=11;
Tube::source_number[5]=12;
Tube::source_number[6]=13;
Tube::source_number[7]=14;
Tube::source_number[8]=22;

Tube::sink_number[0]=24;
Tube::sink_number[1]=25;
Tube::sink_number[2]=26;
Tube::sink_number[3]=27;
Tube::sink_number[4]=28;
Tube::sink_number[5]=29;
Tube::sink_number[6]=30;
Tube::sink_number[7]=31;
Tube::sink_number[8]=35;





if (rm_val == 0){pout << "rm = " << Arteries[int(array[ii][0])]->rm << "\t\t";}
else{pout << "rm = " << Arteries[int(array[ii][0])]->rm << "\t";}


pout << "init = "<<Arteries[int(array[ii][0])]->init<<"\n";
}
pout<<"\nVessels are made\n\n";


	// before real calculation, run the BV codes first to initialize

  // Workspace used by bound_match
  for(int i=0; i<8; i++) fj[i] = new double[8];

  // Workspace used by bound_bif
  for(int i=0; i<18; i++) fjac[i] = new double[18];
  for(int i=0; i<12; i++) fmono[i] = new double[12];
  for(int i=0; i<24; i++) ftrif[i] = new double[24];


	printf("Period:  %f \n numHeartBeats:  %d \n tmstps:  %d \n", Tube::Period, numHeartBeats,Tube::tmstps);

	Tube::onePeriodDt = Tube::Period/Tube::tmstps;


	printf("tstart:  %f \n finaltime:  %f \n tend:  %f \n onePeriodDt:  %f tmstps: %i \n", Tube::tstart1, Tube::finaltime1, Tube::tend1, Tube::onePeriodDt, Tube::tmstps);


















// make the call to the main part of the code
    run_example(argc, argv, Arteries);


        if (verbosity)
    {fprintf(stdout,"\n");}

  // In order to termate the program correctly the vessel network and hence
  // all the vessels and their workspace are deleted.
  for (int i=0; i<nbrves; i++) delete Arteries[i];

  // Matrices and arrays are deleted
  for (int i=0; i<18; i++) delete[] fjac[i];
  for (int i=0; i<8;  i++) delete[] fj[i];

  for (int i=0; i<12; i++) delete[] fmono[i];
  for (int i=0; i<24; i++) delete[] ftrif[i];







    return 0;
}
