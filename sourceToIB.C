#include <sourceToIB.h>

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

#include "sourceToDarcy.h"
#include "MechanicsModel.h"
#include "DarcyMass.h"
#include "DarcyPressure.h"
#include "DarcyVelocity.h"
#include "ConoraryFlowInterface.h"
#include "arteriesPD.h"

//the following variables can be overwritten in the main function

std::string sourceToIB::SstarSystemName = "Sstar_system";
std::string sourceToIB::SstarvariableName = "Sstar_variable";
Order sourceToIB::Sstarvariable_order = FIRST;

double sourceToIB::time;
double sourceToIB::time_step;

//*****************************************************************************

void sourceToIB::setup_Sstar(EquationSystems* es)
{
    es->add_system<LinearImplicitSystem> (SstarSystemName);
    LinearImplicitSystem& system_Sstar = es->get_system<LinearImplicitSystem>(SstarSystemName);
    system_Sstar.add_variable(SstarvariableName, Sstarvariable_order, LAGRANGE);
    return;
}

//*****************************************************************************

void sourceToIB::initialize_Sstar(EquationSystems* es)
{
    LinearImplicitSystem& system_Sstar = es->get_system<LinearImplicitSystem>(SstarSystemName);
    system_Sstar.attach_assemble_function(assemble_Sstar);
    unsigned int system_num = system_Sstar.number();

    libMesh::MeshBase& mesh = es->get_mesh();

    libMesh::NumericVector<double>& Sstar_soln = *system_Sstar.solution;

    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();

    double Sstar_0 = 0.0;
    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;
        const unsigned int dof_num = nd->dof_number(system_num, 0,0);
        Sstar_soln.set(dof_num, Sstar_0);
    }
    Sstar_soln.close();
    Sstar_soln.localize(*system_Sstar.current_local_solution);

    return;
}

//*****************************************************************************

void sourceToIB::solve_Sstar(EquationSystems* es,const double loop_time,const double dt)
{

    time=loop_time;
    LinearImplicitSystem& system_Sstar = es->get_system<LinearImplicitSystem>(SstarSystemName);
    system_Sstar.solve();

    return;
}

//*****************************************************************************

void sourceToIB::assemble_Sstar(EquationSystems& es,
                                 const std::string& system_name)
{

    libmesh_assert_equal_to (system_name, SstarSystemName);
    const MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    LinearImplicitSystem & system_Sstar = es.get_system<LinearImplicitSystem> (SstarSystemName);
    LinearImplicitSystem & system_wx = es.get_system<LinearImplicitSystem> (DarcyVelocity::wxSystemName);
    LinearImplicitSystem & system_wy = es.get_system<LinearImplicitSystem> (DarcyVelocity::wySystemName);
    LinearImplicitSystem & system_wz = es.get_system<LinearImplicitSystem> (DarcyVelocity::wzSystemName);
    LinearImplicitSystem& system_source = es.get_system<LinearImplicitSystem>(sourceToDarcy::sourceToDarcySystemName);

    const DofMap & dof_map_Sstar = system_Sstar.get_dof_map();
    const DofMap & dof_map_wx = system_wx.get_dof_map();
    const DofMap & dof_map_wy = system_wy.get_dof_map();
    const DofMap & dof_map_wz = system_wz.get_dof_map();
    const DofMap & dof_map_source = system_source.get_dof_map();

//***************************************************************************
    FEType fe_type = dof_map_Sstar.variable_type(0);
    UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
    QGauss qrule (dim, FIFTH); // Not sure what the best accuracy is here
    fe->attach_quadrature_rule (&qrule);
    UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));
    QGauss qface(dim-1, FIFTH); // Not sure what the best accuracy is here
//***************************************************************************

    fe_face->attach_quadrature_rule (&qface);

    const std::vector<Real> & JxW = fe->get_JxW();
    const std::vector<std::vector<Real>> & phi = fe->get_phi();
    const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    std::vector<dof_id_type> dof_indices_Sstar;
    std::vector<dof_id_type> dof_indices_wx;
    std::vector<dof_id_type> dof_indices_wy;
    std::vector<dof_id_type> dof_indices_wz;
    std::vector<dof_id_type> dof_indices_source;
    std::vector<dof_id_type> dof_indices_Q;
//****************************************************************************


    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    for (;  el != end_el; ++el)
    {

        Elem* elem = *el;

        dof_map_Sstar.dof_indices (elem, dof_indices_Sstar);
        dof_map_wx.dof_indices (elem, dof_indices_wx);
        dof_map_wy.dof_indices (elem, dof_indices_wy);
        dof_map_wz.dof_indices (elem, dof_indices_wz);
        dof_map_source.dof_indices (elem, dof_indices_source);
//****************************************************************************


        fe->reinit (elem);

        Ke.resize (dof_indices_Sstar.size(),
                   dof_indices_Sstar.size());

        Fe.resize (dof_indices_Sstar.size());

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

            double wx_cur = 0.0;
            double wy_cur = 0.0;
            double wz_cur = 0.0;

            double source_cur = 0.0;
            double RHS_Sstar = 0.0;

            for (unsigned int l=0; l<phi.size(); l++)
            {
                wx_cur += dphi[l][qp](0)*system_wx.current_solution(dof_indices_wx[l]);
                wy_cur += dphi[l][qp](1)*system_wy.current_solution(dof_indices_wy[l]);
                wz_cur += dphi[l][qp](2)*system_wz.current_solution(dof_indices_wz[l]);

                source_cur += phi[l][qp]*system_source.current_solution(dof_indices_source[l]);
            }


//---------------------------------------------------------------------------------------

            for (std::size_t i=0; i<phi.size(); i++)
            {
                for (std::size_t j=0; j<phi.size(); j++)
                {
                    Ke(i,j) += JxW[qp]*phi[j][qp]*phi[i][qp];
                }
            }

        //    RHS_Sstar=source_cur-wx_cur-wy_cur-wz_cur;
RHS_Sstar=source_cur;

            for (std::size_t i=0; i<phi.size(); i++)
            {
                Fe(i) += JxW[qp]*(RHS_Sstar)*(phi[i][qp]);
            }
        }
        dof_map_Sstar.constrain_element_matrix_and_vector (Ke, Fe, dof_indices_Sstar);
        system_Sstar.matrix->add_matrix (Ke, dof_indices_Sstar);
        system_Sstar.rhs->add_vector    (Fe, dof_indices_Sstar);
    }

    return;

}



//*****************************************************************************
