#include <FeedbackToNetwork.h>

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

#include "ConoraryFlowInterface.h"
#include "MechanicsModel.h"
#include "DarcyPressure.h"

std::string FeedbackToNetwork::FeedbackSystemName = "feedback_system";
std::string FeedbackToNetwork::FeedbackvariableName = "feedback_variable";
Order FeedbackToNetwork::Feedbackvariable_order = FIRST;


namespace //private name space
{

static const double kappa_1 = 20000.0;
static const double kappa_2 = 330.0;

}


void FeedbackToNetwork::solve_feedback(EquationSystems* es)
{
    LinearImplicitSystem& system_feedback = es->get_system<LinearImplicitSystem>(FeedbackSystemName);
    system_feedback.solve();
    return;
}

//*****************************************************************************
//*****************************************************************************

void FeedbackToNetwork::setup_feedback(EquationSystems* es)
{
    es->add_system<LinearImplicitSystem> (FeedbackSystemName);
    LinearImplicitSystem& system_feedback = es->get_system<LinearImplicitSystem>(FeedbackSystemName);
    system_feedback.add_variable(FeedbackvariableName, Feedbackvariable_order, LAGRANGE);
    return;
}

//*****************************************************************************
//*****************************************************************************

void FeedbackToNetwork::initialize_feedback(EquationSystems* es)
{
    LinearImplicitSystem& system_feedback = es->get_system<LinearImplicitSystem>(FeedbackSystemName);
    system_feedback.attach_assemble_function(assemble_feedback);
    unsigned int system_num = system_feedback.number();

    libMesh::MeshBase& mesh = es->get_mesh();
    libMesh::NumericVector<double>& feedback_soln = *system_feedback.solution;
    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();

    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();

    double feedback_0 = 0.0;
    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;
        const unsigned int dof_num = nd->dof_number(system_num, 0,0);
        feedback_soln.set(dof_num, feedback_0);
    }
    feedback_soln.close();
    feedback_soln.localize(*system_feedback.current_local_solution);

    return;
}



//*****************************************************************************
//*****************************************************************************

void FeedbackToNetwork::assemble_feedback(EquationSystems& es,const std::string& system_name)
{

    libmesh_assert_equal_to (system_name, FeedbackSystemName);
    const MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    LinearImplicitSystem & system_feedback = es.get_system<LinearImplicitSystem> (FeedbackSystemName);
    System& system_Q_elem = es.get_system<System>(ConoraryFlowInterface::Q_elem_systemName);
    LinearImplicitSystem & system_p = es.get_system<LinearImplicitSystem> (DarcyPressure::pSystemName);

    LinearImplicitSystem& J_system = es.get_system<LinearImplicitSystem>(MechanicsModel::J_system_name);
    LinearImplicitSystem& F_system = es.get_system<LinearImplicitSystem>(MechanicsModel::FF_system_name);

    const DofMap & dof_map_feedback = system_feedback.get_dof_map();
    const libMesh::DofMap& dof_map_Q = system_Q_elem.get_dof_map();
    const DofMap & dof_map_p = system_p.get_dof_map();
    const DofMap & dof_map_J = J_system.get_dof_map();
    const DofMap & dof_map_F = F_system.get_dof_map();

//***************************************************************************
    FEType fe_type = dof_map_feedback.variable_type(0);
    UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
    QGauss qrule (dim, FIFTH);
    fe->attach_quadrature_rule (&qrule);
    UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));
    QGauss qface(dim-1, FIFTH);
//***************************************************************************

    const unsigned int F00_var = 0;
    const unsigned int F01_var = 1;
    const unsigned int F02_var = 2;
    const unsigned int F10_var = 3;
    const unsigned int F11_var = 4;
    const unsigned int F12_var = 5;
    const unsigned int F20_var = 6;
    const unsigned int F21_var = 7;
    const unsigned int F22_var = 8;

    TensorValue<double> CC;
    TensorValue<double> BB;
    TensorValue<double> FF;
    TensorValue<double> FF_inv_tr;
    TensorValue<double> FF_tr;
    TensorValue<double> sigma_total;

    fe_face->attach_quadrature_rule (&qface);

    const std::vector<Real> & JxW = fe->get_JxW();

    const std::vector<std::vector<Real>> & phi = fe->get_phi();
    const std::vector<std::vector<Real>> & phi_in = fe->get_phi();
    const std::vector<std::vector<Real>> & phi_out = fe->get_phi();


    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    std::vector<dof_id_type> dof_indices_feedback;
    std::vector<dof_id_type> dof_indices_Q;
    std::vector<dof_id_type> dof_indices_p;

    std::vector<dof_id_type> dof_indices_J;

    std::vector<dof_id_type> dof_indices_F00;
    std::vector<dof_id_type> dof_indices_F01;
    std::vector<dof_id_type> dof_indices_F02;

    std::vector<dof_id_type> dof_indices_F10;
    std::vector<dof_id_type> dof_indices_F11;
    std::vector<dof_id_type> dof_indices_F12;

    std::vector<dof_id_type> dof_indices_F20;
    std::vector<dof_id_type> dof_indices_F21;
    std::vector<dof_id_type> dof_indices_F22;

//****************************************************************************

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    for (;  el != end_el; ++el)
    {

        Elem* elem = *el;

        dof_map_feedback.dof_indices (elem, dof_indices_feedback);
        dof_map_Q.dof_indices (elem, dof_indices_Q, 0);
        dof_map_p.dof_indices (elem, dof_indices_p);

        dof_map_J.dof_indices (elem, dof_indices_J);

        dof_map_F.dof_indices (elem, dof_indices_F00, F00_var);
        dof_map_F.dof_indices (elem, dof_indices_F01, F01_var);
        dof_map_F.dof_indices (elem, dof_indices_F02, F02_var);

        dof_map_F.dof_indices (elem, dof_indices_F10, F10_var);
        dof_map_F.dof_indices (elem, dof_indices_F11, F11_var);
        dof_map_F.dof_indices (elem, dof_indices_F12, F12_var);

        dof_map_F.dof_indices (elem, dof_indices_F20, F20_var);
        dof_map_F.dof_indices (elem, dof_indices_F21, F21_var);
        dof_map_F.dof_indices (elem, dof_indices_F22, F22_var);



        fe->reinit (elem);

        Ke.resize (dof_indices_feedback.size(),
                   dof_indices_feedback.size());

        Fe.resize (dof_indices_feedback.size());

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

          double p_cur = 0.0;

          double J_cur = 0.0;

          double F00_cur = 0.0;
          double F01_cur = 0.0;
          double F02_cur = 0.0;

          double F10_cur = 0.0;
          double F11_cur = 0.0;
          double F12_cur = 0.0;

          double F20_cur = 0.0;
          double F21_cur = 0.0;
          double F22_cur = 0.0;

          double rhs_total = 0.0;

          for (unsigned int l=0; l<phi.size(); l++)
          {
              p_cur += phi[l][qp]*system_p.current_solution(dof_indices_p[l]);

              F00_cur += phi[l][qp]*F_system.current_solution(dof_indices_F00[l]);
              F01_cur += phi[l][qp]*F_system.current_solution(dof_indices_F01[l]);
              F02_cur += phi[l][qp]*F_system.current_solution(dof_indices_F02[l]);

              F10_cur += phi[l][qp]*F_system.current_solution(dof_indices_F10[l]);
              F11_cur += phi[l][qp]*F_system.current_solution(dof_indices_F11[l]);
              F12_cur += phi[l][qp]*F_system.current_solution(dof_indices_F12[l]);

              F20_cur += phi[l][qp]*F_system.current_solution(dof_indices_F20[l]);
              F21_cur += phi[l][qp]*F_system.current_solution(dof_indices_F21[l]);
              F22_cur += phi[l][qp]*F_system.current_solution(dof_indices_F22[l]);

              J_cur += phi[l][qp]*J_system.current_solution(dof_indices_J[l]);

          }

          FF(0,0)=F00_cur;
          FF(0,1)=F01_cur;
          FF(0,2)=F02_cur;
          FF(1,0)=F10_cur;
          FF(1,1)=F11_cur;
          FF(1,2)=F12_cur;
          FF(2,0)=F20_cur;
          FF(2,1)=F21_cur;
          FF(2,2)=F22_cur;

          FF_tr(0,0)=F00_cur;
        	FF_tr(0,1)=F10_cur;
        	FF_tr(0,2)=F20_cur;
        	FF_tr(1,0)=F01_cur;
        	FF_tr(1,1)=F11_cur;
        	FF_tr(1,2)=F21_cur;
        	FF_tr(2,0)=F02_cur;
        	FF_tr(2,1)=F12_cur;
        	FF_tr(2,2)=F22_cur;

          FF_inv_tr=tensor_inverse_transpose(FF,3);

          CC = FF_tr * FF;
          BB = FF * FF_tr;

          TensorValue<double> II(1.0, 0.0, 0.0, 0.0, 1.0,0.0, 0.0, 0.0, 1.0);
          double I1 = CC.tr();

          sigma_total= (2.0*kappa_1*(FF - FF_inv_tr) + 2.0*kappa_2*(I1*II - BB)*(FF - FF_inv_tr))*(1.0/J_cur)*FF_tr;

          double pf_cur = 0.0;
          pf_cur=(-1.0/3.0)*(sigma_total(0,0)+sigma_total(1,1)+sigma_total(2,2));

          rhs_total=pf_cur+p_cur;

            for (std::size_t i=0; i<phi_out.size(); i++)
            {
                for (std::size_t j=0; j<phi_out.size(); j++)
                {
                    Ke(i,j) += JxW[qp]*phi_out[j][qp]*phi_out[i][qp];
                }
            }

            for (std::size_t i=0; i<phi_out.size(); i++)
            {
                Fe(i) += JxW[qp]*(rhs_total)*(phi_out[i][qp]);
            }

        }

        dof_map_feedback.constrain_element_matrix_and_vector (Ke, Fe, dof_indices_feedback);
        system_feedback.matrix->add_matrix (Ke, dof_indices_feedback);
        system_feedback.rhs->add_vector    (Fe, dof_indices_feedback);

    }


    return;


}
