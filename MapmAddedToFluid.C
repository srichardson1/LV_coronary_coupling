// APPLICATION INCLUDES
#include <MapmAddedToFluid.h>

//libmesh includes
#include <libmesh/equation_systems.h>
#include <libmesh/explicit_system.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_tools.h>

//IBFE includes
#include <ibamr/IBFEMethod.h>



// STATIC VARIABLES
std::string MapmAddedToFluid::m_added_var_name = "m_added";
std::string MapmAddedToFluid::m_added_ctx_name = "m_added_context";

std::string MapmAddedToFluid::mAddedSystemName = "m_added_system";
std::string MapmAddedToFluid::mAddedVariableName = "m_added";

// CLASS IMPLEMENTATION
void MapmAddedToFluid::initialize_equation_systems(
                                                   EquationSystems* equation_systems)
{
    
    System& m_added_system = equation_systems->add_system<System>(mAddedSystemName);
    m_added_system.add_variable(mAddedVariableName, FIRST, LAGRANGE);
    m_added_system.init();
    
    return;
}


void MapmAddedToFluid::initialize_m_added(EquationSystems *es)
{
    System& m_added_system = es->get_system<System>(mAddedSystemName);
    unsigned int system_num = m_added_system.number();
    libMesh::MeshBase& mesh = es->get_mesh();
    
    const libMesh::DofMap& dof_map = m_added_system.get_dof_map();
    libMesh::NumericVector<double>& m_added_soln = *m_added_system.solution;
    
    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();
    
    double m_add_0 = 0.5;
    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;
        const unsigned int dof_num = nd->dof_number(system_num, 0,0);
        m_added_soln.set(dof_num, m_add_0);
    }
    m_added_soln.close();
    m_added_soln.localize(*m_added_system.current_local_solution);
    
    return;
    
}

void MapmAddedToFluid::spread_m_added_solid_fluid(
                                                  SAMRAI::tbox::Pointer<IBFEMethod> ib_method_ops,
                                                  SAMRAI::tbox::Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                                  const int m_added_data_idx,
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
            Pointer< CellData<NDIM, double> > D_data = patch->getPatchData(m_added_data_idx);
            D_data->fillAll(0.0);
        }
    }
    
    //spreading m_added from solid to fluid
    pout<<"speading added m from structure to fluid\n";
    
    libMesh::PetscVector<double> * d_X_current_vec;
    libMesh::PetscVector<double> * d_X_IB_ghost_vec;
    libMesh::PetscVector<double> * d_m_added_current_vec;
    libMesh::PetscVector<double> * d_m_added_IB_ghost_Vec;
    
    
    FEDataManager* d_fe_data_manager = ib_method_ops->getFEDataManager();
    EquationSystems* es = d_fe_data_manager->getEquationSystems();
    
    //initializate from the system
    
    System& d_X_system = es->get_system<System>(IBFEMethod::COORDS_SYSTEM_NAME);
    d_X_current_vec = dynamic_cast<PetscVector<double>* >(d_X_system.current_local_solution.get());
    d_X_IB_ghost_vec = dynamic_cast<PetscVector<double>* >(
                                                           d_fe_data_manager->buildGhostedCoordsVector(/*localize_data*/ false));
    
    
    System& d_m_added_system = es->get_system<System>(mAddedSystemName);
    d_m_added_current_vec = dynamic_cast<PetscVector<double>* >(d_m_added_system.solution.get());
    d_m_added_IB_ghost_Vec = dynamic_cast<PetscVector<double>* >(
                                                                 d_fe_data_manager->buildGhostedSolutionVector(mAddedSystemName, /*localized_data*/ false) );
    
    PetscVector<double>* X_vec = d_X_current_vec;
    PetscVector<double>* m_added_vec = d_m_added_current_vec;
    PetscVector<double>* X_ghost_vec = d_X_IB_ghost_vec;
    PetscVector<double>* m_added_ghost_vec = d_m_added_IB_ghost_Vec;
    X_vec->localize(*X_ghost_vec);
    m_added_vec->localize(*m_added_ghost_vec);
    
    //speading now
    d_fe_data_manager->spread(
                              m_added_data_idx, *m_added_ghost_vec, *X_ghost_vec, mAddedSystemName, NULL, data_time);
    
    //may need some postprocessing
    //spreading, then we can decide, if D_data > 1, then D_data = 1; if D_data<0, then D_data = 0;
    
    return;
}
