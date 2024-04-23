// APPLICATION INCLUDES
#include <MapMeshToFluid.h>

//libmesh includes
#include <libmesh/equation_systems.h>
#include <libmesh/explicit_system.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_tools.h>

//IBFE includes
#include <ibamr/IBFEMethod.h>



// STATIC VARIABLES
std::string MapMeshToFluid::phi_mesh_var_name = "phi_mesh";
std::string MapMeshToFluid::phi_mesh_ctx_name = "phi_mesh_context";

std::string MapMeshToFluid::phiMeshSystemName = "phi_mesh_system";
std::string MapMeshToFluid::phiMeshVariableName = "phi_mesh";

// CLASS IMPLEMENTATION
void MapMeshToFluid::initialize_equation_systems(
                                                 EquationSystems* equation_systems)
{
    
    System& phi_system = equation_systems->add_system<System>(phiMeshSystemName);
    phi_system.add_variable(phiMeshVariableName, FIRST, LAGRANGE);
    phi_system.init();
    
    return;
}


void MapMeshToFluid::initialize_phi(EquationSystems *es)
{
    System& phi_system = es->get_system<System>(phiMeshSystemName);
    unsigned int system_num = phi_system.number();
    libMesh::MeshBase& mesh = es->get_mesh();
    
    const libMesh::DofMap& dof_map = phi_system.get_dof_map();
    libMesh::NumericVector<double>& phi_soln = *phi_system.solution;
    
    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();
    
    double phi_0 = 1.0;
    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;
        const unsigned int dof_num = nd->dof_number(system_num, 0,0);
        phi_soln.set(dof_num, phi_0);
    }
    phi_soln.close();
    phi_soln.localize(*phi_system.current_local_solution);
    
    return;
    
}

void MapMeshToFluid::spread_phi_solid_fluid(
                                            SAMRAI::tbox::Pointer<IBFEMethod> ib_method_ops,
                                            SAMRAI::tbox::Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                            const int phi_data_idx,
                                            const double data_time)
{
    
    //redefine all values to be zero in patch_hierarchy fpr phi_data_idx
    for(int ln=0; ln<=patch_hierarchy->getFinestLevelNumber(); ln++)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for( PatchLevel<NDIM>::Iterator p(level); p; p++ )
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            const Box<NDIM>& patch_box = patch->getBox();
            Pointer< CellData<NDIM, double> > D_data = patch->getPatchData(phi_data_idx);
            D_data->fillAll(0.0);
        }
    }
    
    //spreading phi from solid to fluid
    pout<<"speading phi from structure to fluid\n";
    
    libMesh::PetscVector<double> * d_X_current_vec;
    libMesh::PetscVector<double> * d_X_IB_ghost_vec;
    libMesh::PetscVector<double> * d_phi_current_vec;
    libMesh::PetscVector<double> * d_phi_IB_ghost_Vec;
    
    
    FEDataManager* d_fe_data_manager = ib_method_ops->getFEDataManager();
    EquationSystems* es = d_fe_data_manager->getEquationSystems();
    
    //initializate from the system
    
    System& d_X_system = es->get_system<System>(IBFEMethod::COORDS_SYSTEM_NAME);
    d_X_current_vec = dynamic_cast<PetscVector<double>* >(d_X_system.current_local_solution.get());
    d_X_IB_ghost_vec = dynamic_cast<PetscVector<double>* >(
                                                           d_fe_data_manager->buildGhostedCoordsVector(/*localize_data*/ false));
    
    
    System& d_phi_system = es->get_system<System>(phiMeshSystemName);
    d_phi_current_vec = dynamic_cast<PetscVector<double>* >(d_phi_system.solution.get());
    d_phi_IB_ghost_Vec = dynamic_cast<PetscVector<double>* >(
                                                             d_fe_data_manager->buildGhostedSolutionVector(phiMeshSystemName, /*localized_data*/ false) );
    
    PetscVector<double>* X_vec = d_X_current_vec;
    PetscVector<double>* phi_vec = d_phi_current_vec;
    PetscVector<double>* X_ghost_vec = d_X_IB_ghost_vec;
    PetscVector<double>* phi_ghost_vec = d_phi_IB_ghost_Vec;
    X_vec->localize(*X_ghost_vec);
    phi_vec->localize(*phi_ghost_vec);
    
    //speading now
    d_fe_data_manager->spread(
                              phi_data_idx, *phi_ghost_vec, *X_ghost_vec, phiMeshSystemName, NULL, data_time);
    
    //may need some postprocessing
    //spreading, then we can decide, if D_data > 1, then D_data = 1; if D_data<0, then D_data = 0;
    
    return;
}

