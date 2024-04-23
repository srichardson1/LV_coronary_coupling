#include <DarcyPressure.h>

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
#include "sourceToIB.h"
#include "ConoraryFlowInterface.h"
#include "arteriesPD.h"

double DarcyPressure::phi0;
double DarcyPressure::biot_modulus;
double DarcyPressure::kappa0;
double DarcyPressure::cap_K;

//the following variables can be overwritten in the main function
std::string DarcyPressure::pSystemName = "p_system";
std::string DarcyPressure::pvariableName = "p_variable";
Order DarcyPressure::pvariable_order = FIRST;

std::string DarcyPressure::pstarSystemName = "pstar_system";
std::string DarcyPressure::pstarvariableName = "pstar_variable";
Order DarcyPressure::pstarvariable_order = FIRST;


double DarcyPressure::time;
double DarcyPressure::time_step;

namespace //private name space
{
    static const TensorValue<double> Perm_Tensor(0.000002,0.0,0.0,0.0,0.000002,0.0,0.0,0.0,0.000002);
    static const double p_0=0.15;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


void DarcyPressure::setup_p(EquationSystems* es)
{
    es->add_system<LinearImplicitSystem> (pSystemName);
    LinearImplicitSystem& system_p = es->get_system<LinearImplicitSystem>(pSystemName);
    system_p.add_variable(pvariableName, pvariable_order, LAGRANGE);
    return;
}

//*****************************************************************************

void DarcyPressure::setup_pstar(EquationSystems* es)
{
    es->add_system<LinearImplicitSystem> (pstarSystemName);
    LinearImplicitSystem& system_pstar = es->get_system<LinearImplicitSystem>(pstarSystemName);
    system_pstar.add_variable(pstarvariableName, pstarvariable_order, LAGRANGE);
    return;
}




void DarcyPressure::initialize_p(EquationSystems* es)
{
    LinearImplicitSystem& system_p = es->get_system<LinearImplicitSystem>(pSystemName);
    system_p.attach_assemble_function(assemble_p);
    unsigned int system_num = system_p.number();

    libMesh::MeshBase& mesh = es->get_mesh();
//    const unsigned int dim = mesh.mesh_dimension();

//    const libMesh::DofMap& dof_map = system_p.get_dof_map();
    libMesh::NumericVector<double>& p_soln = *system_p.solution;

    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();

    double p_0 = 0.0;
    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;
        const unsigned int dof_num = nd->dof_number(system_num, 0,0);
        p_soln.set(dof_num, p_0);
    }
    p_soln.close();
    p_soln.localize(*system_p.current_local_solution);

    return;
}

//*****************************************************************************


void DarcyPressure::initialize_pstar(EquationSystems* es)
{
    LinearImplicitSystem& system_pstar = es->get_system<LinearImplicitSystem>(pstarSystemName);
    system_pstar.attach_assemble_function(assemble_pstar);
    unsigned int system_num = system_pstar.number();

    libMesh::MeshBase& mesh = es->get_mesh();
//    const unsigned int dim = mesh.mesh_dimension();

//    const libMesh::DofMap& dof_map = system_pstar.get_dof_map();
    libMesh::NumericVector<double>& pstar_soln = *system_pstar.solution;

    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();

    double pstar_0 = 0.0;
    for(; node_it != node_end; ++node_it)
    {
     	const libMesh::Node* nd = *node_it;
        const unsigned int dof_num = nd->dof_number(system_num, 0,0);
        pstar_soln.set(dof_num, pstar_0);
    }
    pstar_soln.close();
    pstar_soln.localize(*system_pstar.current_local_solution);

    return;
}

//*****************************************************************************


void DarcyPressure::solve_p(EquationSystems* es,const double loop_time,const double dt)
{
    time=loop_time;
    time_step=dt;
    LinearImplicitSystem& system_p = es->get_system<LinearImplicitSystem>(pSystemName);
    system_p.solve();

    return;
}

//*****************************************************************************

void DarcyPressure::solve_pstar(EquationSystems* es,const double loop_time,const double dt)
{
    time=loop_time;
    time_step=dt;
    LinearImplicitSystem& system_pstar = es->get_system<LinearImplicitSystem>(pstarSystemName);
    system_pstar.solve();

    return;
}

//----------------------------------------------------------------------------

void DarcyPressure::assemble_p(EquationSystems& es,
                             const std::string& system_name)
{
    libmesh_assert_equal_to (system_name, pSystemName);
    const MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    LinearImplicitSystem & system_p = es.get_system<LinearImplicitSystem> (pSystemName);
    LinearImplicitSystem & system_pstar = es.get_system<LinearImplicitSystem> (pstarSystemName);
    LinearImplicitSystem & system_m = es.get_system<LinearImplicitSystem> (DarcyMass::mSystemName);
    LinearImplicitSystem & system_mstar = es.get_system<LinearImplicitSystem> (DarcyMass::mstarSystemName);
    LinearImplicitSystem & system_mold = es.get_system<LinearImplicitSystem> (DarcyMass::moldSystemName);
    LinearImplicitSystem& F_system = es.get_system<LinearImplicitSystem>(MechanicsModel::FF_system_name);
    LinearImplicitSystem& dJdt_system = es.get_system<LinearImplicitSystem>(MechanicsModel::dJdt_system_name);
    LinearImplicitSystem& system_source = es.get_system<LinearImplicitSystem>(sourceToDarcy::sourceToDarcySystemName);

    LinearImplicitSystem& J_system = es.get_system<LinearImplicitSystem>(MechanicsModel::J_system_name);
    LinearImplicitSystem & system_jold = es.get_system<LinearImplicitSystem> (MechanicsModel::Jold_system_name);


    const DofMap & dof_map_m = system_m.get_dof_map();
    const DofMap & dof_map_mstar = system_mstar.get_dof_map();
    const DofMap & dof_map_jold = system_jold.get_dof_map();
    const DofMap & dof_map_mold = system_mold.get_dof_map();
    const DofMap & dof_map_pstar = system_pstar.get_dof_map();
    const DofMap & dof_map_J = J_system.get_dof_map();
    const DofMap & dof_map_F = F_system.get_dof_map();
    const DofMap & dof_map_dJdt = dJdt_system.get_dof_map();
    const DofMap & dof_map_p = system_p.get_dof_map();
    const DofMap & dof_map_source = system_source.get_dof_map();


    FEType fe_type = dof_map_p.variable_type(0);
    UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
    QGauss qrule (dim, FIFTH); // Not sure what the best accuracy is here
    fe->attach_quadrature_rule (&qrule);
    UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));
    QGauss qface(dim-1, FIFTH); // Not sure what the best accuracy is here
    //----------------------------------------------------------------------------------

    const unsigned int F00_var = 0;
    const unsigned int F01_var = 1;
    const unsigned int F02_var = 2;
    const unsigned int F10_var = 3;
    const unsigned int F11_var = 4;
    const unsigned int F12_var = 5;
    const unsigned int F20_var = 6;
    const unsigned int F21_var = 7;
    const unsigned int F22_var = 8;


    TensorValue<double> final_input; // The assembled LHS matrices included for convenience
    TensorValue<double> FF_def; //FF_def contains the deformation gradient
    TensorValue<double> FF_def_inv; // This is the inverse of the deformation gradient
    TensorValue<double> FF_def_inv_tr; // The is the transpose of the deformation gradient inverse
    double RHS_darcy; // The assembled RHS value included for convenience
    double cap_B;

    //----------------------------------------------------------------------------------

    fe_face->attach_quadrature_rule (&qface);

    const std::vector<Real> & JxW = fe->get_JxW();
    const std::vector<std::vector<Real>> & phi = fe->get_phi();
    const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    std::vector<dof_id_type> dof_indices_pstar;

    std::vector<dof_id_type> dof_indices_source;

    std::vector<dof_id_type> dof_indices_jold;
    std::vector<dof_id_type> dof_indices_mold;
    std::vector<dof_id_type> dof_indices_mstar;
    std::vector<dof_id_type> dof_indices_m;
    std::vector<dof_id_type> dof_indices_p;

    std::vector<dof_id_type> dof_indices_J;
    std::vector<dof_id_type> dof_indices_dJdt;

    std::vector<dof_id_type> dof_indices_F00;
    std::vector<dof_id_type> dof_indices_F01;
    std::vector<dof_id_type> dof_indices_F02;

    std::vector<dof_id_type> dof_indices_F10;
    std::vector<dof_id_type> dof_indices_F11;
    std::vector<dof_id_type> dof_indices_F12;

    std::vector<dof_id_type> dof_indices_F20;
    std::vector<dof_id_type> dof_indices_F21;
    std::vector<dof_id_type> dof_indices_F22;
    std::vector<dof_id_type> dof_indices_Q;



//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

    //traditional way of looping all elements to be compatible with 1.0 version
    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    for (;  el != end_el; ++el)
        //   for (const auto & elem : mesh.active_local_element_ptr_range()) //this is a new feature
    {
        Elem* elem = *el;


        dof_map_pstar.dof_indices (elem, dof_indices_pstar);
        dof_map_source.dof_indices (elem, dof_indices_source);

        dof_map_jold.dof_indices (elem, dof_indices_jold);
        dof_map_mstar.dof_indices (elem, dof_indices_mstar);
        dof_map_m.dof_indices (elem, dof_indices_m);
        dof_map_mold.dof_indices (elem, dof_indices_mold);
        dof_map_J.dof_indices (elem, dof_indices_J);
        dof_map_dJdt.dof_indices(elem, dof_indices_dJdt);
        dof_map_p.dof_indices (elem, dof_indices_p);

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

        Ke.resize (dof_indices_p.size(),
                   dof_indices_p.size());

        Fe.resize (dof_indices_p.size());


        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

            double p_cur = 0.0;
            double j_old = 0.0;
            double m_old = 0.0;
            double mstar_cur = 0.0;
            double pstar_cur = 0.0;
            double m_cur = 0.0;

            double J_cur = 0.0;
            double dJdt_cur = 0.0;
            double source2=0.0;

            double F00_cur = 0.0;
            double F01_cur = 0.0;
            double F02_cur = 0.0;

            double F10_cur = 0.0;
            double F11_cur = 0.0;
            double F12_cur = 0.0;

            double F20_cur = 0.0;
            double F21_cur = 0.0;
            double F22_cur = 0.0;

            double source_cur = 0.0;
            double source_other = 0.0;

            double penalty=0.0;
            double penalty_old=0.0;
            double penalty_change=0.0;

            double djdt_cur=0.0;
            double f_of_j=0.0;
            double f_of_j_old=0.0;
            double f_prime_of_j=0.0;
            double f_prime_of_j_alternative=0.0;
            double RHS_extra=0.0;
            double RHS_extra_coro=0.0;


            for (unsigned int l=0; l<phi.size(); l++)
            {

                F00_cur += phi[l][qp]*F_system.current_solution(dof_indices_F00[l]);
                F01_cur += phi[l][qp]*F_system.current_solution(dof_indices_F01[l]);
                F02_cur += phi[l][qp]*F_system.current_solution(dof_indices_F02[l]);

                F10_cur += phi[l][qp]*F_system.current_solution(dof_indices_F10[l]);
                F11_cur += phi[l][qp]*F_system.current_solution(dof_indices_F11[l]);
                F12_cur += phi[l][qp]*F_system.current_solution(dof_indices_F12[l]);

                F20_cur += phi[l][qp]*F_system.current_solution(dof_indices_F20[l]);
                F21_cur += phi[l][qp]*F_system.current_solution(dof_indices_F21[l]);
                F22_cur += phi[l][qp]*F_system.current_solution(dof_indices_F22[l]);

                j_old += phi[l][qp]*system_jold.current_solution(dof_indices_jold[l]);
                J_cur += phi[l][qp]*J_system.current_solution(dof_indices_J[l]);

                pstar_cur += phi[l][qp]*system_pstar.current_solution(dof_indices_pstar[l]);
                p_cur += phi[l][qp]*system_p.current_solution(dof_indices_p[l]);

                mstar_cur += phi[l][qp]*system_mstar.current_solution(dof_indices_mstar[l]);
                m_old += phi[l][qp]*system_mold.current_solution(dof_indices_mold[l]);
                m_cur += phi[l][qp]*system_m.current_solution(dof_indices_m[l]);

                source_cur += phi[l][qp]*system_source.current_solution(dof_indices_source[l]);

                dJdt_cur += phi[l][qp]*dJdt_system.current_solution(dof_indices_dJdt[l]);
            }

//---------------------------------------------------------------------------------------

            FF_def(0,0)=F00_cur;
            FF_def(0,1)=F01_cur;
            FF_def(0,2)=F02_cur;
            FF_def(1,0)=F10_cur;
            FF_def(1,1)=F11_cur;
            FF_def(1,2)=F12_cur;
            FF_def(2,0)=F20_cur;
            FF_def(2,1)=F21_cur;
            FF_def(2,2)=F22_cur;


            if (J_cur==1.0)
            {
                f_of_j=1.0;
                f_of_j_old=1.0;
                f_prime_of_j=(-2.0/3.0)*(dJdt_cur);
            }
            else
            {
                f_of_j=(2.0*(J_cur-1.0-log(J_cur)))/((J_cur-1.0)*(J_cur-1.0));
                f_of_j_old=(2.0*(j_old-1.0-log(j_old)))/((j_old-1.0)*(j_old-1.0));
                f_prime_of_j=((((2.0)*(J_cur-1.0)*(1.0-((1.0)/(J_cur))))-(4.0*(J_cur-1.0-log(J_cur))))/((J_cur-1.0)*(J_cur-1.0)*(J_cur-1.0)))*(dJdt_cur);
            }

            f_prime_of_j_alternative=(f_of_j-f_of_j_old)/(time_step);

            if (mstar_cur+phi0<0.0001)
            {
                penalty=0.0;
                penalty_old=0.0;
            }
            else
            {
                penalty=(1.0)/(m_cur+phi0);
                penalty_old=(1.0)/(m_old+phi0);
            }

            penalty_change=(penalty-penalty_old)/(time_step);

            FF_def_inv.zero();
            FF_def_inv_tr.zero();

            FF_def_inv=tensor_inverse(FF_def,NDIM);
            FF_def_inv_tr=tensor_inverse_transpose(FF_def,NDIM);

            cap_B=biot_modulus*f_of_j;

            final_input.zero();
            final_input=(cap_B)*(1.0)*(J_cur)*(FF_def_inv)*(Perm_Tensor)*(FF_def_inv_tr);

            double a1=final_input(0,0);
            double a2=final_input(0,1);
            double a3=final_input(0,2);

            double a4=final_input(1,0);
            double a5=final_input(1,1);
            double a6=final_input(1,2);

            double a7=final_input(2,0);
            double a8=final_input(2,1);
            double a9=final_input(2,2);

            //------------------------------------------------------------------------------------------------
            for (std::size_t i=0; i<phi.size(); i++)
            {
                for (std::size_t j=0; j<phi.size(); j++)
                {

                    Ke(i,j) += (JxW[qp]*(phi[i][qp]*phi[j][qp]))+((0.5*time_step)*(JxW[qp])*(((dphi[i][qp](0))*(((a1)*(dphi[j][qp](0)))+((a2)*(dphi[j][qp](1)))+((a3)*(dphi[j][qp](2)))))+((dphi[i][qp](1))*(((a4)*(dphi[j][qp](0)))+((a5)*(dphi[j][qp](1)))+((a6)*(dphi[j][qp](2)))))+((dphi[i][qp](2))*(((a7)*(dphi[j][qp](0)))+((a8)*(dphi[j][qp](1)))+((a9)*(dphi[j][qp](2)))))));

                }
            }
            //------------------------------------------------------------------------------------------------

            RHS_extra=(f_of_j*source_cur*biot_modulus)-(f_of_j*dJdt_cur*biot_modulus)+(f_prime_of_j*biot_modulus*(1.0-J_cur+mstar_cur))-(kappa0*penalty_change);

        //    RHS_extra=-(f_of_j*dJdt_cur*biot_modulus)+(f_prime_of_j*biot_modulus*(1.0-J_cur+mstar_cur))-(kappa0*penalty_change);
        //    RHS_extra_coro=(f_of_j*source_cur*biot_modulus);


            RHS_darcy=pstar_cur+(0.5*time_step*RHS_extra);

            for (std::size_t i=0; i<phi.size(); i++)
            {
                Fe(i) += JxW[qp]*RHS_darcy*phi[i][qp];
            }

        }
        //--------------------------------------------------------------------------------------------


        dof_map_p.constrain_element_matrix_and_vector (Ke, Fe, dof_indices_p);
        system_p.matrix->add_matrix (Ke, dof_indices_p);
        system_p.rhs->add_vector    (Fe, dof_indices_p);
    }



    return;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void DarcyPressure::assemble_pstar(EquationSystems& es,
                             const std::string& system_name)
{
    libmesh_assert_equal_to (system_name, pstarSystemName);
    const MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    LinearImplicitSystem & system_p = es.get_system<LinearImplicitSystem> (pSystemName);
    LinearImplicitSystem & system_pstar = es.get_system<LinearImplicitSystem> (pstarSystemName);
    LinearImplicitSystem & system_m = es.get_system<LinearImplicitSystem> (DarcyMass::mSystemName);
    LinearImplicitSystem & system_mold = es.get_system<LinearImplicitSystem> (DarcyMass::moldSystemName);
    LinearImplicitSystem & system_jold = es.get_system<LinearImplicitSystem> (MechanicsModel::Jold_system_name);
    LinearImplicitSystem& J_system = es.get_system<LinearImplicitSystem>(MechanicsModel::J_system_name);
    LinearImplicitSystem& F_system = es.get_system<LinearImplicitSystem>(MechanicsModel::FF_system_name);
    LinearImplicitSystem& dJdt_system = es.get_system<LinearImplicitSystem>(MechanicsModel::dJdt_system_name);
    LinearImplicitSystem& system_source = es.get_system<LinearImplicitSystem>(sourceToDarcy::sourceToDarcySystemName);

    const DofMap & dof_map_m = system_m.get_dof_map();
    const DofMap & dof_map_mold = system_mold.get_dof_map();
    const DofMap & dof_map_jold = system_jold.get_dof_map();
    const DofMap & dof_map_pstar = system_pstar.get_dof_map();
    const DofMap & dof_map_J = J_system.get_dof_map();
    const DofMap & dof_map_F = F_system.get_dof_map();
    const DofMap & dof_map_dJdt = dJdt_system.get_dof_map();
    const DofMap & dof_map_p = system_p.get_dof_map();
    const DofMap & dof_map_source = system_source.get_dof_map();

    FEType fe_type = dof_map_pstar.variable_type(0);
    UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
    QGauss qrule (dim, FIFTH); // Not sure what the best accuracy is here
    fe->attach_quadrature_rule (&qrule);
    UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));
    QGauss qface(dim-1, FIFTH); // Not sure what the best accuracy is here
    //----------------------------------------------------------------------------------

    const unsigned int F00_var = 0;
    const unsigned int F01_var = 1;
    const unsigned int F02_var = 2;
    const unsigned int F10_var = 3;
    const unsigned int F11_var = 4;
    const unsigned int F12_var = 5;
    const unsigned int F20_var = 6;
    const unsigned int F21_var = 7;
    const unsigned int F22_var = 8;

    TensorValue<double> final_input; // The assembled LHS matrices included for convenience
    TensorValue<double> FF_def; //FF_def contains the deformation gradient
    TensorValue<double> FF_def_inv; // This is the inverse of the deformation gradient
    TensorValue<double> FF_def_inv_tr; // The is the transpose of the deformation gradient inverse
    double RHS_darcy; // The assembled RHS value included for convenience
    double cap_B;

    //----------------------------------------------------------------------------------

    fe_face->attach_quadrature_rule (&qface);

    const std::vector<Real> & JxW = fe->get_JxW();
    const std::vector<std::vector<Real>> & phi = fe->get_phi();
    const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    std::vector<dof_id_type> dof_indices_source;

    std::vector<dof_id_type> dof_indices_pstar;

    std::vector<dof_id_type> dof_indices_jold;
    std::vector<dof_id_type> dof_indices_mold;
    std::vector<dof_id_type> dof_indices_m;
    std::vector<dof_id_type> dof_indices_p;

    std::vector<dof_id_type> dof_indices_J;
    std::vector<dof_id_type> dof_indices_dJdt;

    std::vector<dof_id_type> dof_indices_F00;
    std::vector<dof_id_type> dof_indices_F01;
    std::vector<dof_id_type> dof_indices_F02;

    std::vector<dof_id_type> dof_indices_F10;
    std::vector<dof_id_type> dof_indices_F11;
    std::vector<dof_id_type> dof_indices_F12;

    std::vector<dof_id_type> dof_indices_F20;
    std::vector<dof_id_type> dof_indices_F21;
    std::vector<dof_id_type> dof_indices_F22;

    //----------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------

    //traditional way of looping all elements to be compatible with 1.0 version
    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    for (;  el != end_el; ++el)
        //   for (const auto & elem : mesh.active_local_element_ptr_range()) //this is a new feature
    {

        Elem* elem = *el;

        dof_map_source.dof_indices (elem, dof_indices_source);
        dof_map_pstar.dof_indices (elem, dof_indices_pstar);

        dof_map_mold.dof_indices (elem, dof_indices_mold);
        dof_map_jold.dof_indices (elem, dof_indices_jold);
        dof_map_m.dof_indices (elem, dof_indices_m);
        dof_map_J.dof_indices (elem, dof_indices_J);
        dof_map_dJdt.dof_indices(elem, dof_indices_dJdt);
        dof_map_p.dof_indices (elem, dof_indices_p);

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

        Ke.resize (dof_indices_pstar.size(),
                   dof_indices_pstar.size());

        Fe.resize (dof_indices_pstar.size());


        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

            double p_cur = 0.0;
            double j_old = 0.0;
            double m_old = 0.0;
            double m_cur = 0.0;
            double pstar_cur = 0.0;

            double source_cur = 0.0;
            double source_other = 0.0;

            double J_cur = 0.0;
            double dJdt_cur = 0.0;

            double F00_cur = 0.0;
            double F01_cur = 0.0;
            double F02_cur = 0.0;

            double F10_cur = 0.0;
            double F11_cur = 0.0;
            double F12_cur = 0.0;

            double F20_cur = 0.0;
            double F21_cur = 0.0;
            double F22_cur = 0.0;

            double penalty=0.0;
            double penalty_old=0.0;
            double penalty_change=0.0;

            double djdt_cur=0.0;
            double f_of_j=0.0;
            double f_of_j_old=0.0;
            double f_prime_of_j=0.0;
            double f_prime_of_j_alternative=0.0;
            double RHS_extra=0.0;
            double source2=0.0;

            for (unsigned int l=0; l<phi.size(); l++)
            {
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
                m_cur += phi[l][qp]*system_m.current_solution(dof_indices_m[l]);
                p_cur += phi[l][qp]*system_p.current_solution(dof_indices_p[l]);

                m_old += phi[l][qp]*system_mold.current_solution(dof_indices_mold[l]);
                j_old += phi[l][qp]*system_jold.current_solution(dof_indices_jold[l]);

                source_cur += phi[l][qp]*system_source.current_solution(dof_indices_source[l]);
                pstar_cur += phi[l][qp]*system_pstar.current_solution(dof_indices_pstar[l]);
                dJdt_cur += phi[l][qp]*dJdt_system.current_solution(dof_indices_dJdt[l]);

            }

//---------------------------------------------------------------------------------------

            FF_def(0,0)=F00_cur;
            FF_def(0,1)=F01_cur;
            FF_def(0,2)=F02_cur;
            FF_def(1,0)=F10_cur;
            FF_def(1,1)=F11_cur;
            FF_def(1,2)=F12_cur;
            FF_def(2,0)=F20_cur;
            FF_def(2,1)=F21_cur;
            FF_def(2,2)=F22_cur;


            if (J_cur==1.0)
            {
             	  f_of_j=1.0;
                f_of_j_old=1.0;
                f_prime_of_j=(-2.0/3.0)*(dJdt_cur);
            }
            else
            {
              	f_of_j=(2.0*(J_cur-1.0-log(J_cur)))/((J_cur-1.0)*(J_cur-1.0));
                f_of_j_old=(2.0*(j_old-1.0-log(j_old)))/((j_old-1.0)*(j_old-1.0));
                f_prime_of_j=((((2.0)*(J_cur-1.0)*(1.0-((1.0)/(J_cur))))-(4.0*(J_cur-1.0-log(J_cur))))/((J_cur-1.0)*(J_cur-1.0)*(J_cur-1.0)))*(dJdt_cur);
            }

            f_prime_of_j_alternative=(f_of_j-f_of_j_old)/(time_step);

            if (m_cur+phi0<0.0001)
            {
             	penalty=0.0;
              penalty_old=0.0;
            }
            else
            {
                penalty=(1.0)/(m_cur+phi0);
                penalty_old=(1.0)/(m_old+phi0);
            }


            penalty_change=(penalty-penalty_old)/(time_step);


            FF_def_inv.zero();
            FF_def_inv_tr.zero();

            FF_def_inv=tensor_inverse(FF_def,NDIM);
            FF_def_inv_tr=tensor_inverse_transpose(FF_def,NDIM);


            cap_B=biot_modulus*f_of_j;


            final_input.zero();
            final_input=(cap_B)*(J_cur)*(FF_def_inv)*(Perm_Tensor)*(FF_def_inv_tr);

            double a1=final_input(0,0);
            double a2=final_input(0,1);
            double a3=final_input(0,2);

            double a4=final_input(1,0);
            double a5=final_input(1,1);
            double a6=final_input(1,2);

            double a7=final_input(2,0);
            double a8=final_input(2,1);
            double a9=final_input(2,2);



//------------------------------------------------------------------------------------------------
            for (std::size_t i=0; i<phi.size(); i++)
            {
                for (std::size_t j=0; j<phi.size(); j++)
                {

                    Ke(i,j) += (JxW[qp]*(phi[i][qp]*phi[j][qp]))+((0.5*time_step)*(JxW[qp])*(((dphi[i][qp](0))*(((a1)*(dphi[j][qp](0)))+((a2)*(dphi[j][qp](1)))+((a3)*(dphi[j][qp](2)))))+((dphi[i][qp](1))*(((a4)*(dphi[j][qp](0)))+((a5)*(dphi[j][qp](1)))+((a6)*(dphi[j][qp](2)))))+((dphi[i][qp](2))*(((a7)*(dphi[j][qp](0)))+((a8)*(dphi[j][qp](1)))+((a9)*(dphi[j][qp](2)))))));

                }
            }

//------------------------------------------------------------------------------------------------

            RHS_extra=(f_of_j*source_cur*biot_modulus)-(f_of_j*dJdt_cur*biot_modulus)+(f_prime_of_j*biot_modulus*(1.0-J_cur+m_cur))-(kappa0*penalty_change);

            RHS_darcy=p_cur+(0.5*time_step*RHS_extra);


            for (std::size_t i=0; i<phi.size(); i++)
            {
                Fe(i) += JxW[qp]*RHS_darcy*phi[i][qp];
            }


        }


        dof_map_pstar.constrain_element_matrix_and_vector (Ke, Fe, dof_indices_pstar);
        system_pstar.matrix->add_matrix (Ke, dof_indices_pstar);
        system_pstar.rhs->add_vector    (Fe, dof_indices_pstar);
    }

    return;
}
