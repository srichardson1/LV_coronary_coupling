#include <MeshConversion.h>

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

#include "ibtk/FEDataInterpolation.h"
#include "ibtk/FEDataManager.h"


//the following variables can be overwritten in the main function
std::string MeshConversion::meshlongSystemName = "meshlong_system";
std::string MeshConversion::meshlongvariableName = "meshlong_variable";
Order MeshConversion::meshlongvariable_order = FIRST;

std::string MeshConversion::meshradialSystemName = "meshradial_system";
std::string MeshConversion::meshradialvariableName = "meshradial_variable";
Order MeshConversion::meshradialvariable_order = FIRST;

std::string MeshConversion::meshcircumSystemName = "meshcircum_system";
std::string MeshConversion::meshcircumvariableName = "meshcircum_variable";
Order MeshConversion::meshcircumvariable_order = FIRST;

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


void MeshConversion::setup_mesh_long(EquationSystems* es)
{
    es->add_system<LinearImplicitSystem> (meshlongSystemName);
    LinearImplicitSystem& system_meshlong = es->get_system<LinearImplicitSystem>(meshlongSystemName);
    system_meshlong.add_variable(meshlongvariableName, meshlongvariable_order, LAGRANGE);
    return;
}

//*****************************************************************************

void MeshConversion::setup_mesh_radial(EquationSystems* es)
{
    es->add_system<LinearImplicitSystem> (meshradialSystemName);
    LinearImplicitSystem& system_meshradial = es->get_system<LinearImplicitSystem>(meshradialSystemName);
    system_meshradial.add_variable(meshradialvariableName, meshradialvariable_order, LAGRANGE);
    return;
}

//*****************************************************************************


void MeshConversion::setup_mesh_circum(EquationSystems* es)
{
    es->add_system<LinearImplicitSystem> (meshcircumSystemName);
    LinearImplicitSystem& system_meshcircum = es->get_system<LinearImplicitSystem>(meshcircumSystemName);
    system_meshcircum.add_variable(meshcircumvariableName, meshcircumvariable_order, LAGRANGE);
    return;
}

//*****************************************************************************


void MeshConversion::initialize_mesh_long(EquationSystems* es)
{
    LinearImplicitSystem& system_meshlong = es->get_system<LinearImplicitSystem>(meshlongSystemName);
    system_meshlong.attach_assemble_function(assemble_mesh_long);
    unsigned int system_num = system_meshlong.number();

    libMesh::MeshBase& mesh = es->get_mesh();

    libMesh::NumericVector<double>& meshlong_soln = *system_meshlong.solution;

    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();

    double meshlong_0 = 0.0;
    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;
        const unsigned int dof_num = nd->dof_number(system_num, 0,0);
        meshlong_soln.set(dof_num, meshlong_0);
    }
    meshlong_soln.close();
    meshlong_soln.localize(*system_meshlong.current_local_solution);

    return;
}

//*****************************************************************************


void MeshConversion::initialize_mesh_radial(EquationSystems* es)
{
    LinearImplicitSystem& system_meshradial = es->get_system<LinearImplicitSystem>(meshradialSystemName);
    system_meshradial.attach_assemble_function(assemble_mesh_radial);
    unsigned int system_num = system_meshradial.number();

    libMesh::MeshBase& mesh = es->get_mesh();

    libMesh::NumericVector<double>& meshradial_soln = *system_meshradial.solution;

    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();

    double meshradial_0 = 0.0;
    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;
        const unsigned int dof_num = nd->dof_number(system_num, 0,0);
        meshradial_soln.set(dof_num, meshradial_0);
    }
    meshradial_soln.close();
    meshradial_soln.localize(*system_meshradial.current_local_solution);

    return;
}

//*****************************************************************************

void MeshConversion::initialize_mesh_circum(EquationSystems* es)
{
    LinearImplicitSystem& system_meshcircum = es->get_system<LinearImplicitSystem>(meshcircumSystemName);
    system_meshcircum.attach_assemble_function(assemble_mesh_circum);
    unsigned int system_num = system_meshcircum.number();

    libMesh::MeshBase& mesh = es->get_mesh();

    libMesh::NumericVector<double>& meshcircum_soln = *system_meshcircum.solution;

    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();

    double meshcircum_0 = 0.0;
    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;
        const unsigned int dof_num = nd->dof_number(system_num, 0,0);
        meshcircum_soln.set(dof_num, meshcircum_0);
    }
    meshcircum_soln.close();
    meshcircum_soln.localize(*system_meshcircum.current_local_solution);

    return;
}

//*****************************************************************************



void MeshConversion::solve_mesh_long(EquationSystems* es)
{
    LinearImplicitSystem& system_meshlong = es->get_system<LinearImplicitSystem>(meshlongSystemName);
    system_meshlong.solve();
    return;
}


void MeshConversion::solve_mesh_radial(EquationSystems* es)
{
    LinearImplicitSystem& system_meshradial = es->get_system<LinearImplicitSystem>(meshradialSystemName);
    system_meshradial.solve();
    return;
}


void MeshConversion::solve_mesh_circum(EquationSystems* es)
{
    LinearImplicitSystem& system_meshcircum = es->get_system<LinearImplicitSystem>(meshcircumSystemName);
    system_meshcircum.solve();
    return;
}

//----------------------------------------------------------------------------

void MeshConversion::assemble_mesh_long(EquationSystems& es,
                             const std::string& system_name)
{
    libmesh_assert_equal_to (system_name, meshlongSystemName);
    const MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    LinearImplicitSystem & system_meshlong = es.get_system<LinearImplicitSystem> (meshlongSystemName);

    const DofMap & dof_map_meshlong = system_meshlong.get_dof_map();


    FEType fe_type = dof_map_meshlong.variable_type(0);
    UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
    QGauss qrule (dim, FIFTH); // Not sure what the best accuracy is here
    fe->attach_quadrature_rule (&qrule);
    UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));
    QGauss qface(dim-1, FIFTH); // Not sure what the best accuracy is here
    //----------------------------------------------------------------------------------

    fe_face->attach_quadrature_rule (&qface);

    const std::vector<Real> & JxW = fe->get_JxW();
    const std::vector<std::vector<Real>> & phi = fe->get_phi();
    const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

//    const std::vector<Point> & q_point = fe->get_xyz();

    const std::vector<libMesh::Point> & q_point = fe->get_xyz();


    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    std::vector<dof_id_type> dof_indices_meshlong;

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    for (;  el != end_el; ++el)

    {

        Elem* elem = *el;

        dof_map_meshlong.dof_indices (elem, dof_indices_meshlong);

        fe->reinit (elem);

        Ke.resize (dof_indices_meshlong.size(),
                   dof_indices_meshlong.size());

        Fe.resize (dof_indices_meshlong.size());

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

            for (std::size_t i=0; i<phi.size(); i++)
            {
                for (std::size_t j=0; j<phi.size(); j++)
                {
                  Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
                }
            }




//------------------------------------------------------------------------------------------------

            for (std::size_t i=0; i<phi.size(); i++)
            {
                Fe(i) += JxW[qp]*0.0*phi[i][qp];
            }

        }



        dof_map_meshlong.constrain_element_matrix_and_vector (Ke, Fe, dof_indices_meshlong);
        system_meshlong.matrix->add_matrix (Ke, dof_indices_meshlong);
        system_meshlong.rhs->add_vector    (Fe, dof_indices_meshlong);
    }

    return;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void MeshConversion::assemble_mesh_radial(EquationSystems& es,
                             const std::string& system_name)
{
    libmesh_assert_equal_to (system_name, meshradialSystemName);
    const MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    LinearImplicitSystem & system_meshradial = es.get_system<LinearImplicitSystem> (meshradialSystemName);

    const DofMap & dof_map_meshradial = system_meshradial.get_dof_map();


    FEType fe_type = dof_map_meshradial.variable_type(0);
    UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
    QGauss qrule (dim, FIFTH); // Not sure what the best accuracy is here
    fe->attach_quadrature_rule (&qrule);
    UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));
    QGauss qface(dim-1, FIFTH); // Not sure what the best accuracy is here
    //----------------------------------------------------------------------------------

    fe_face->attach_quadrature_rule (&qface);

    const std::vector<Real> & JxW = fe->get_JxW();
    const std::vector<std::vector<Real>> & phi = fe->get_phi();
    const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

//    const std::vector<Point> & q_point = fe->get_xyz();

    const std::vector<libMesh::Point> & q_point = fe->get_xyz();


    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    std::vector<dof_id_type> dof_indices_meshradial;

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    for (;  el != end_el; ++el)

    {

        Elem* elem = *el;

        dof_map_meshradial.dof_indices (elem, dof_indices_meshradial);

        fe->reinit (elem);

        Ke.resize (dof_indices_meshradial.size(),
                   dof_indices_meshradial.size());

        Fe.resize (dof_indices_meshradial.size());

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

            for (std::size_t i=0; i<phi.size(); i++)
            {
                for (std::size_t j=0; j<phi.size(); j++)
                {
                  Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
                }
            }




//------------------------------------------------------------------------------------------------

            for (std::size_t i=0; i<phi.size(); i++)
            {
                Fe(i) += JxW[qp]*0.0*phi[i][qp];
            }

        }



        dof_map_meshradial.constrain_element_matrix_and_vector (Ke, Fe, dof_indices_meshradial);
        system_meshradial.matrix->add_matrix (Ke, dof_indices_meshradial);
        system_meshradial.rhs->add_vector    (Fe, dof_indices_meshradial);
    }

    return;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void MeshConversion::assemble_mesh_circum(EquationSystems& es,
                             const std::string& system_name)
{
    libmesh_assert_equal_to (system_name, meshcircumSystemName);
    const MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    LinearImplicitSystem & system_meshcircum = es.get_system<LinearImplicitSystem> (meshcircumSystemName);

    const DofMap & dof_map_meshcircum = system_meshcircum.get_dof_map();


    FEType fe_type = dof_map_meshcircum.variable_type(0);
    UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
    QGauss qrule (dim, FIFTH); // Not sure what the best accuracy is here
    fe->attach_quadrature_rule (&qrule);
    UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));
    QGauss qface(dim-1, FIFTH); // Not sure what the best accuracy is here
    //----------------------------------------------------------------------------------

    fe_face->attach_quadrature_rule (&qface);

    const std::vector<Real> & JxW = fe->get_JxW();
    const std::vector<std::vector<Real>> & phi = fe->get_phi();
    const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

//    const std::vector<Point> & q_point = fe->get_xyz();

    const std::vector<libMesh::Point> & q_point = fe->get_xyz();


    libMesh::NumericVector<double>& meshcircum_soln = *system_meshcircum.solution;




    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    std::vector<dof_id_type> dof_indices_meshcircum;

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();


//    meshcircum_soln.close();
//    meshcircum_soln.localize(*meshcircum_system.current_local_solution);


    return;
}
