#include <sourceToDarcy.h>

// LIBMESH INCLUDES
// not sure how many are needed, just pull in all of them
#include <libmesh/centroid_partitioner.h>
#include <libmesh/equation_systems.h>
#include <libmesh/explicit_system.h>
#include "libmesh/linear_implicit_system.h"
#include <libmesh/mesh.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_tools.h>

#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"

#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"

#include "libmesh/perf_log.h"
#include "libmesh/elem.h"

#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/analytic_function.h"

#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"

#include <libmesh/boundary_info.h>

#include "ibtk/FEDataInterpolation.h"
#include "ibtk/FEDataManager.h"

#include "ConoraryFlowInterface.h"
#include "DarcyPressure.h"
#include "arteriesPD.h"

std::string sourceToDarcy::sourceToDarcySystemName = "source_to_darcy_system";
std::string sourceToDarcy::sourceToDarcyvariableName = "source_to_darcy_variable";
Order sourceToDarcy::sourceToDarcyvariable_order = FIRST;

BoundaryInfo* sourceToDarcy::boundary_info;

double sourceToDarcy::p_a;
double sourceToDarcy::p_v;
double sourceToDarcy::beta_a;
double sourceToDarcy::beta_v;

double sourceToDarcy::time_std;
double sourceToDarcy::time_step_std;

void sourceToDarcy::solve_std(EquationSystems* es)
{
    LinearImplicitSystem& system_std = es->get_system<LinearImplicitSystem>(sourceToDarcySystemName);
    system_std.solve();


    Tube::summed_source = system_std.solution->sum();
    Tube::summed_source_all_time=Tube::summed_source_all_time+Tube::summed_source;


    return;
}

//*****************************************************************************
//*****************************************************************************

void sourceToDarcy::setup_std(EquationSystems* es)
{
    es->add_system<LinearImplicitSystem> (sourceToDarcySystemName);
    LinearImplicitSystem& system_std = es->get_system<LinearImplicitSystem>(sourceToDarcySystemName);
    system_std.add_variable(sourceToDarcyvariableName, sourceToDarcyvariable_order, LAGRANGE);
    return;
}

//*****************************************************************************
//*****************************************************************************

void sourceToDarcy::initialize_std(EquationSystems* es)
{
    LinearImplicitSystem& system_std = es->get_system<LinearImplicitSystem>(sourceToDarcySystemName);
    system_std.attach_assemble_function(assemble_std);
    unsigned int system_num = system_std.number();

    libMesh::MeshBase& mesh = es->get_mesh();
    libMesh::NumericVector<double>& std_soln = *system_std.solution;
    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();

    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();

    double std_0 = 0.0;
    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;
        const unsigned int dof_num = nd->dof_number(system_num, 0,0);
        std_soln.set(dof_num, std_0);
    }
    std_soln.close();
    std_soln.localize(*system_std.current_local_solution);

    return;
}



//*****************************************************************************
//*****************************************************************************

void sourceToDarcy::assemble_std(EquationSystems& es,const std::string& system_name)
{

    libmesh_assert_equal_to (system_name, sourceToDarcySystemName);
    const MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    LinearImplicitSystem & system_std = es.get_system<LinearImplicitSystem> (sourceToDarcySystemName);
    LinearImplicitSystem & system_p = es.get_system<LinearImplicitSystem> (DarcyPressure::pSystemName);

    const libMesh::DofMap& dof_map_std = system_std.get_dof_map();
    const DofMap & dof_map_p = system_p.get_dof_map();

//***************************************************************************
    FEType fe_type = dof_map_std.variable_type(0);
    UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
    QGauss qrule (dim, FIFTH);
    fe->attach_quadrature_rule (&qrule);
    UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));
    QGauss qface(dim-1, FIFTH);
//***************************************************************************

    fe_face->attach_quadrature_rule (&qface);

    const std::vector<Real> & JxW = fe->get_JxW();
    const std::vector<std::vector<Real>> & phi = fe->get_phi();
    const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

std::vector<dof_id_type> dof_indices_p;
std::vector<dof_id_type> dof_indices_std;
//****************************************************************************


    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

    for (;  el != end_el; ++el)
    {

      Elem* elem = *el;

        dof_map_p.dof_indices (elem, dof_indices_p);
        dof_map_std.dof_indices (elem, dof_indices_std);



        fe->reinit (elem);

        Ke.resize (dof_indices_std.size(),
                   dof_indices_std.size());

        Fe.resize (dof_indices_std.size());


        double sink_source_from_coro=0.0;


        if(ConoraryFlowInterface::element_current_subdomain_lists[elem->id()] == 1)
        {
          sink_source_from_coro=Tube::net_from_coro_averaged_in_dimensional_space[0];
        }
        else if(ConoraryFlowInterface::element_current_subdomain_lists[elem->id()] == 2)
        {
          sink_source_from_coro=Tube::net_from_coro_averaged_in_dimensional_space[1];
        }
        else if(ConoraryFlowInterface::element_current_subdomain_lists[elem->id()] == 3)
        {
          sink_source_from_coro=Tube::net_from_coro_averaged_in_dimensional_space[2];
        }
        else if(ConoraryFlowInterface::element_current_subdomain_lists[elem->id()] == 4)
        {
          sink_source_from_coro=Tube::net_from_coro_averaged_in_dimensional_space[3];
        }
        else if(ConoraryFlowInterface::element_current_subdomain_lists[elem->id()] == 5)
        {
          sink_source_from_coro=Tube::net_from_coro_averaged_in_dimensional_space[4];
        }
        else if(ConoraryFlowInterface::element_current_subdomain_lists[elem->id()] == 6)
        {
          sink_source_from_coro=Tube::net_from_coro_averaged_in_dimensional_space[5];
        }
        else if(ConoraryFlowInterface::element_current_subdomain_lists[elem->id()] == 7)
        {
          sink_source_from_coro=Tube::net_from_coro_averaged_in_dimensional_space[6];
        }
        else if(ConoraryFlowInterface::element_current_subdomain_lists[elem->id()] == 8)
        {
          sink_source_from_coro=Tube::net_from_coro_averaged_in_dimensional_space[7];
        }
        else if(ConoraryFlowInterface::element_current_subdomain_lists[elem->id()] == 9)
        {
          sink_source_from_coro=Tube::net_from_coro_averaged_in_dimensional_space[8];
        }
        else
        {pout << "something has gone wrong, element has not been assigned a subdomain" << "\n";}




        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

            for (std::size_t i=0; i<phi.size(); i++)
            {
                for (std::size_t j=0; j<phi.size(); j++)
                {
                    Ke(i,j) += JxW[qp]*phi[j][qp]*phi[i][qp];
                }
            }

            for (std::size_t i=0; i<phi.size(); i++)
            {
                Fe(i) += JxW[qp]*(sink_source_from_coro)*(phi[i][qp]);
            }

        }


        dof_map_std.constrain_element_matrix_and_vector (Ke, Fe, dof_indices_std);
        system_std.matrix->add_matrix (Ke, dof_indices_std);
        system_std.rhs->add_vector    (Fe, dof_indices_std);
    }
    return;

}
