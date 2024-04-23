#include <DarcyVelocity.h>

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

#include "MechanicsModel.h"
#include "DarcyPressure.h"

std::string DarcyVelocity::wxSystemName = "wx_system";
std::string DarcyVelocity::wxvariableName = "wx_variable";
Order DarcyVelocity::wxvariable_order = FIRST;

std::string DarcyVelocity::wySystemName = "wy_system";
std::string DarcyVelocity::wyvariableName = "wy_variable";
Order DarcyVelocity::wyvariable_order = FIRST;

std::string DarcyVelocity::wzSystemName = "wz_system";
std::string DarcyVelocity::wzvariableName = "wz_variable";
Order DarcyVelocity::wzvariable_order = FIRST;


std::string DarcyVelocity::wxstarSystemName = "wxstar_system";
std::string DarcyVelocity::wxstarvariableName = "wxstar_variable";
Order DarcyVelocity::wxstarvariable_order = FIRST;

std::string DarcyVelocity::wystarSystemName = "wystar_system";
std::string DarcyVelocity::wystarvariableName = "wystar_variable";
Order DarcyVelocity::wystarvariable_order = FIRST;

std::string DarcyVelocity::wzstarSystemName = "wzstar_system";
std::string DarcyVelocity::wzstarvariableName = "wzstar_variable";
Order DarcyVelocity::wzstarvariable_order = FIRST;


double DarcyVelocity::time;
double DarcyVelocity::time_step;

namespace //private name space
{
    static const TensorValue<double> Perm_Tensor(0.000002,0.0,0.0,0.0,0.000002,0.0,0.0,0.0,0.000002);
}


void DarcyVelocity::setup_wx(EquationSystems* es)
{
    es->add_system<LinearImplicitSystem> (wxSystemName);
    LinearImplicitSystem& system_wx = es->get_system<LinearImplicitSystem>(wxSystemName);
    system_wx.add_variable(wxvariableName, wxvariable_order, LAGRANGE);
    return;
}

//*****************************************************************************

void DarcyVelocity::setup_wy(EquationSystems* es)
{
    es->add_system<LinearImplicitSystem> (wySystemName);
    LinearImplicitSystem& system_wy = es->get_system<LinearImplicitSystem>(wySystemName);
    system_wy.add_variable(wyvariableName, wyvariable_order, LAGRANGE);
    return;
}

//*****************************************************************************

void DarcyVelocity::setup_wz(EquationSystems* es)
{
    es->add_system<LinearImplicitSystem> (wzSystemName);
    LinearImplicitSystem& system_wz = es->get_system<LinearImplicitSystem>(wzSystemName);
    system_wz.add_variable(wzvariableName, wzvariable_order, LAGRANGE);
    return;
}

//*****************************************************************************

void DarcyVelocity::setup_wxstar(EquationSystems* es)
{
    es->add_system<LinearImplicitSystem> (wxstarSystemName);
    LinearImplicitSystem& system_wxstar = es->get_system<LinearImplicitSystem>(wxstarSystemName);
    system_wxstar.add_variable(wxstarvariableName, wxstarvariable_order, LAGRANGE);
    return;
}

//*****************************************************************************

void DarcyVelocity::setup_wystar(EquationSystems* es)
{
    es->add_system<LinearImplicitSystem> (wystarSystemName);
    LinearImplicitSystem& system_wystar = es->get_system<LinearImplicitSystem>(wystarSystemName);
    system_wystar.add_variable(wystarvariableName, wystarvariable_order, LAGRANGE);
    return;
}

//*****************************************************************************

void DarcyVelocity::setup_wzstar(EquationSystems* es)
{
    es->add_system<LinearImplicitSystem> (wzstarSystemName);
    LinearImplicitSystem& system_wzstar = es->get_system<LinearImplicitSystem>(wzstarSystemName);
    system_wzstar.add_variable(wzstarvariableName, wzstarvariable_order, LAGRANGE);
    return;
}



//*****************************************************************************

void DarcyVelocity::initialize_wx(EquationSystems* es)
{
    LinearImplicitSystem& system_wx = es->get_system<LinearImplicitSystem>(wxSystemName);
    system_wx.attach_assemble_function(assemble_wx);
    unsigned int system_num = system_wx.number();

    libMesh::MeshBase& mesh = es->get_mesh();
//    const unsigned int dim = mesh.mesh_dimension();

//    const libMesh::DofMap& dof_map = system_wx.get_dof_map();
    libMesh::NumericVector<double>& wx_soln = *system_wx.solution;

    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();

    double wx_0 = 0.0;
    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;
        const unsigned int dof_num = nd->dof_number(system_num, 0,0);
        wx_soln.set(dof_num, wx_0);
    }
    wx_soln.close();
    wx_soln.localize(*system_wx.current_local_solution);

    return;
}

//*****************************************************************************

void DarcyVelocity::initialize_wy(EquationSystems* es)
{
    LinearImplicitSystem& system_wy = es->get_system<LinearImplicitSystem>(wySystemName);
    system_wy.attach_assemble_function(assemble_wy);
    unsigned int system_num = system_wy.number();

    libMesh::MeshBase& mesh = es->get_mesh();
//    const unsigned int dim = mesh.mesh_dimension();

//    const libMesh::DofMap& dof_map = system_wy.get_dof_map();
    libMesh::NumericVector<double>& wy_soln = *system_wy.solution;

    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();

    double wy_0 = 0.0;
    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;
        const unsigned int dof_num = nd->dof_number(system_num, 0,0);
        wy_soln.set(dof_num, wy_0);
    }
    wy_soln.close();
    wy_soln.localize(*system_wy.current_local_solution);

    return;
}

//*****************************************************************************

void DarcyVelocity::initialize_wz(EquationSystems* es)
{
    LinearImplicitSystem& system_wz = es->get_system<LinearImplicitSystem>(wzSystemName);
    system_wz.attach_assemble_function(assemble_wz);
    unsigned int system_num = system_wz.number();

    libMesh::MeshBase& mesh = es->get_mesh();
//    const unsigned int dim = mesh.mesh_dimension();

//    const libMesh::DofMap& dof_map = system_wz.get_dof_map();
    libMesh::NumericVector<double>& wz_soln = *system_wz.solution;

    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();

    double wz_0 = 0.0;
    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;
        const unsigned int dof_num = nd->dof_number(system_num, 0,0);
        wz_soln.set(dof_num, wz_0);
    }
    wz_soln.close();
    wz_soln.localize(*system_wz.current_local_solution);

    return;
}


//*****************************************************************************

void DarcyVelocity::initialize_wxstar(EquationSystems* es)
{
    LinearImplicitSystem& system_wxstar = es->get_system<LinearImplicitSystem>(wxstarSystemName);
    system_wxstar.attach_assemble_function(assemble_wxstar);
    unsigned int system_num = system_wxstar.number();

    libMesh::MeshBase& mesh = es->get_mesh();
//    const unsigned int dim = mesh.mesh_dimension();

//    const libMesh::DofMap& dof_map = system_wxstar.get_dof_map();
    libMesh::NumericVector<double>& wxstar_soln = *system_wxstar.solution;

    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();

    double wxstar_0 = 0.0;
    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;
        const unsigned int dof_num = nd->dof_number(system_num, 0,0);
        wxstar_soln.set(dof_num, wxstar_0);
    }
    wxstar_soln.close();
    wxstar_soln.localize(*system_wxstar.current_local_solution);

    return;
}

//*****************************************************************************

void DarcyVelocity::initialize_wystar(EquationSystems* es)
{
    LinearImplicitSystem& system_wystar = es->get_system<LinearImplicitSystem>(wystarSystemName);
    system_wystar.attach_assemble_function(assemble_wystar);
    unsigned int system_num = system_wystar.number();

    libMesh::MeshBase& mesh = es->get_mesh();
//    const unsigned int dim = mesh.mesh_dimension();

//    const libMesh::DofMap& dof_map = system_wystar.get_dof_map();
    libMesh::NumericVector<double>& wystar_soln = *system_wystar.solution;

    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();

    double wystar_0 = 0.0;
    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;
        const unsigned int dof_num = nd->dof_number(system_num, 0,0);
        wystar_soln.set(dof_num, wystar_0);
    }
    wystar_soln.close();
    wystar_soln.localize(*system_wystar.current_local_solution);

    return;
}

//*****************************************************************************

void DarcyVelocity::initialize_wzstar(EquationSystems* es)
{
    LinearImplicitSystem& system_wzstar = es->get_system<LinearImplicitSystem>(wzstarSystemName);
    system_wzstar.attach_assemble_function(assemble_wzstar);
    unsigned int system_num = system_wzstar.number();

    libMesh::MeshBase& mesh = es->get_mesh();
//    const unsigned int dim = mesh.mesh_dimension();

//    const libMesh::DofMap& dof_map = system_wzstar.get_dof_map();
    libMesh::NumericVector<double>& wzstar_soln = *system_wzstar.solution;

    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();

    double wzstar_0 = 0.0;
    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;
        const unsigned int dof_num = nd->dof_number(system_num, 0,0);
        wzstar_soln.set(dof_num, wzstar_0);
    }
    wzstar_soln.close();
    wzstar_soln.localize(*system_wzstar.current_local_solution);

    return;
}


//*****************************************************************************

void DarcyVelocity::solve_wx(EquationSystems* es)
{
    LinearImplicitSystem& system_wx = es->get_system<LinearImplicitSystem>(wxSystemName);
    system_wx.solve();
    return;
}

//*****************************************************************************

void DarcyVelocity::solve_wy(EquationSystems* es)
{
    LinearImplicitSystem& system_wy = es->get_system<LinearImplicitSystem>(wySystemName);
    system_wy.solve();
    return;
}

//*****************************************************************************

void DarcyVelocity::solve_wz(EquationSystems* es)
{
    LinearImplicitSystem& system_wz = es->get_system<LinearImplicitSystem>(wzSystemName);
    system_wz.solve();
    return;
}

//*****************************************************************************

void DarcyVelocity::solve_wxstar(EquationSystems* es)
{
    LinearImplicitSystem& system_wxstar = es->get_system<LinearImplicitSystem>(wxstarSystemName);
    system_wxstar.solve();
    return;
}

//*****************************************************************************

void DarcyVelocity::solve_wystar(EquationSystems* es)
{
    LinearImplicitSystem& system_wystar = es->get_system<LinearImplicitSystem>(wystarSystemName);
    system_wystar.solve();
    return;
}

//*****************************************************************************

void DarcyVelocity::solve_wzstar(EquationSystems* es)
{
    LinearImplicitSystem& system_wzstar = es->get_system<LinearImplicitSystem>(wzstarSystemName);
    system_wzstar.solve();
    return;
}




//*****************************************************************************

void DarcyVelocity::assemble_wx(EquationSystems& es,
                              const std::string& system_name)
{

    libmesh_assert_equal_to (system_name, wxSystemName);
    const MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    LinearImplicitSystem & system_wx = es.get_system<LinearImplicitSystem> (wxSystemName);
    LinearImplicitSystem & system_p = es.get_system<LinearImplicitSystem> (DarcyPressure::pSystemName);
    LinearImplicitSystem& J_system = es.get_system<LinearImplicitSystem>(MechanicsModel::J_system_name);
    LinearImplicitSystem& F_system = es.get_system<LinearImplicitSystem>(MechanicsModel::FF_system_name);

    const DofMap & dof_map_wx = system_wx.get_dof_map();
    const DofMap & dof_map_p = system_p.get_dof_map();
    const DofMap & dof_map_J = J_system.get_dof_map();
    const DofMap & dof_map_F = F_system.get_dof_map();

    //***************************************************************************
    FEType fe_type = dof_map_wx.variable_type(0);
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

    TensorValue<double> final_input; // The assembled RHS matrices
    TensorValue<double> FF_def; //FF_def contains the deformation gradient
    TensorValue<double> FF_def_inv; // This is the inverse of the deformation gradient
    TensorValue<double> FF_def_inv_tr; // The is the transpose of the deformation gradient inverse
 //   double f_of_J;

    fe_face->attach_quadrature_rule (&qface);

    const std::vector<Real> & JxW = fe->get_JxW();
//    const std::vector<libMesh::Point> & q_point = fe->get_xyz();
    const std::vector<std::vector<Real>> & phi = fe->get_phi();
    const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    std::vector<dof_id_type> dof_indices_p;
    std::vector<dof_id_type> dof_indices_wx;
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


        dof_map_p.dof_indices (elem, dof_indices_p);
        dof_map_wx.dof_indices (elem, dof_indices_wx);
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


        //****************************************************************************

        fe->reinit (elem);

        Ke.resize (dof_indices_wx.size(),
                   dof_indices_wx.size());

        Fe.resize (dof_indices_wx.size());

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

            //****************************************************************************

  //          double A=0.0;
            double RHS_wx=0.0;

            double px_cur = 0.0;
            double py_cur = 0.0;
            double pz_cur = 0.0;

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

                px_cur += dphi[l][qp](0)*system_p.current_solution(dof_indices_p[l]);
                py_cur += dphi[l][qp](1)*system_p.current_solution(dof_indices_p[l]);
                pz_cur += dphi[l][qp](2)*system_p.current_solution(dof_indices_p[l]);
            }


            FF_def(0,0)=F00_cur;
            FF_def(0,1)=F01_cur;
            FF_def(0,2)=F02_cur;
            FF_def(1,0)=F10_cur;
            FF_def(1,1)=F11_cur;
            FF_def(1,2)=F12_cur;
            FF_def(2,0)=F20_cur;
            FF_def(2,1)=F21_cur;
            FF_def(2,2)=F22_cur;

            FF_def_inv.zero();
            FF_def_inv_tr.zero();

            FF_def_inv=tensor_inverse(FF_def,NDIM);
            FF_def_inv_tr=tensor_inverse_transpose(FF_def,NDIM);


            for (std::size_t i=0; i<phi.size(); i++)
            {
                for (std::size_t j=0; j<phi.size(); j++)
                {
                    Ke(i,j) += JxW[qp]*phi[j][qp]*phi[i][qp];
                }
            }

            final_input=-(J_cur)*(FF_def_inv)*(Perm_Tensor)*(FF_def_inv_tr);
            RHS_wx=((final_input(0,0))*(px_cur))+((final_input(0,1))*(py_cur))+((final_input(0,2))*(pz_cur));

            for (std::size_t i=0; i<phi.size(); i++)
            {
                Fe(i) += JxW[qp]*(RHS_wx)*(phi[i][qp]);
            }
        }
        dof_map_wx.constrain_element_matrix_and_vector (Ke, Fe, dof_indices_wx);
        system_wx.matrix->add_matrix (Ke, dof_indices_wx);
        system_wx.rhs->add_vector    (Fe, dof_indices_wx);
    }


    return;
}


//*****************************************************************************

void DarcyVelocity::assemble_wy(EquationSystems& es,
                              const std::string& system_name)
{


    libmesh_assert_equal_to (system_name, wySystemName);
    const MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    LinearImplicitSystem & system_wy = es.get_system<LinearImplicitSystem> (wySystemName);
    LinearImplicitSystem & system_p = es.get_system<LinearImplicitSystem> (DarcyPressure::pSystemName);
    LinearImplicitSystem& J_system = es.get_system<LinearImplicitSystem>(MechanicsModel::J_system_name);
    LinearImplicitSystem& F_system = es.get_system<LinearImplicitSystem>(MechanicsModel::FF_system_name);

    const DofMap & dof_map_wy = system_wy.get_dof_map();
    const DofMap & dof_map_p = system_p.get_dof_map();
    const DofMap & dof_map_J = J_system.get_dof_map();
    const DofMap & dof_map_F = F_system.get_dof_map();

    //***************************************************************************
    FEType fe_type = dof_map_wy.variable_type(0);
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

    TensorValue<double> final_input; // The assembled RHS matrices
    TensorValue<double> FF_def; //FF_def contains the deformation gradient
    TensorValue<double> FF_def_inv; // This is the inverse of the deformation gradient
    TensorValue<double> FF_def_inv_tr; // The is the transpose of the deformation gradient inverse

    fe_face->attach_quadrature_rule (&qface);

    const std::vector<Real> & JxW = fe->get_JxW();
//    const std::vector<libMesh::Point> & q_point = fe->get_xyz();
    const std::vector<std::vector<Real>> & phi = fe->get_phi();
    const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    std::vector<dof_id_type> dof_indices_p;
    std::vector<dof_id_type> dof_indices_wy;
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


        dof_map_p.dof_indices (elem, dof_indices_p);
        dof_map_wy.dof_indices (elem, dof_indices_wy);
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


        //****************************************************************************

        fe->reinit (elem);

        Ke.resize (dof_indices_wy.size(),
                   dof_indices_wy.size());

        Fe.resize (dof_indices_wy.size());

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

            //****************************************************************************

  //          double A=0.0;
            double RHS_wy=0.0;

            double px_cur = 0.0;
            double py_cur = 0.0;
            double pz_cur = 0.0;

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

                px_cur += dphi[l][qp](0)*system_p.current_solution(dof_indices_p[l]);
                py_cur += dphi[l][qp](1)*system_p.current_solution(dof_indices_p[l]);
                pz_cur += dphi[l][qp](2)*system_p.current_solution(dof_indices_p[l]);
            }


            FF_def(0,0)=F00_cur;
            FF_def(0,1)=F01_cur;
            FF_def(0,2)=F02_cur;
            FF_def(1,0)=F10_cur;
            FF_def(1,1)=F11_cur;
            FF_def(1,2)=F12_cur;
            FF_def(2,0)=F20_cur;
            FF_def(2,1)=F21_cur;
            FF_def(2,2)=F22_cur;

            FF_def_inv.zero();
            FF_def_inv_tr.zero();

            FF_def_inv=tensor_inverse(FF_def,NDIM);
            FF_def_inv_tr=tensor_inverse_transpose(FF_def,NDIM);


            for (std::size_t i=0; i<phi.size(); i++)
            {
                for (std::size_t j=0; j<phi.size(); j++)
                {
                    Ke(i,j) += JxW[qp]*phi[j][qp]*phi[i][qp];
                }
            }

            final_input=-(J_cur)*(FF_def_inv)*(Perm_Tensor)*(FF_def_inv_tr);
            RHS_wy=((final_input(1,0))*(px_cur))+((final_input(1,1))*(py_cur))+((final_input(1,2))*(pz_cur));

            for (std::size_t i=0; i<phi.size(); i++)
            {

                Fe(i) += JxW[qp]*(RHS_wy)*(phi[i][qp]);

            }
        }
        dof_map_wy.constrain_element_matrix_and_vector (Ke, Fe, dof_indices_wy);
        system_wy.matrix->add_matrix (Ke, dof_indices_wy);
        system_wy.rhs->add_vector    (Fe, dof_indices_wy);
    }


    return;
}


//*****************************************************************************

void DarcyVelocity::assemble_wz(EquationSystems& es,
                              const std::string& system_name)
{

    libmesh_assert_equal_to (system_name, wzSystemName);
    const MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    LinearImplicitSystem & system_wz = es.get_system<LinearImplicitSystem> (wzSystemName);
    LinearImplicitSystem & system_p = es.get_system<LinearImplicitSystem> (DarcyPressure::pSystemName);
    LinearImplicitSystem& J_system = es.get_system<LinearImplicitSystem>(MechanicsModel::J_system_name);
    LinearImplicitSystem& F_system = es.get_system<LinearImplicitSystem>(MechanicsModel::FF_system_name);

    const DofMap & dof_map_wz = system_wz.get_dof_map();
    const DofMap & dof_map_p = system_p.get_dof_map();
    const DofMap & dof_map_J = J_system.get_dof_map();
    const DofMap & dof_map_F = F_system.get_dof_map();

    //***************************************************************************
    FEType fe_type = dof_map_wz.variable_type(0);
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

    TensorValue<double> final_input; // The assembled RHS matrices
    TensorValue<double> FF_def; //FF_def contains the deformation gradient
    TensorValue<double> FF_def_inv; // This is the inverse of the deformation gradient
    TensorValue<double> FF_def_inv_tr; // The is the transpose of the deformation gradient inverse


    fe_face->attach_quadrature_rule (&qface);

    const std::vector<Real> & JxW = fe->get_JxW();
 //   const std::vector<libMesh::Point> & q_point = fe->get_xyz();
    const std::vector<std::vector<Real>> & phi = fe->get_phi();
    const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    std::vector<dof_id_type> dof_indices_p;
    std::vector<dof_id_type> dof_indices_wz;
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

        dof_map_p.dof_indices (elem, dof_indices_p);
        dof_map_wz.dof_indices (elem, dof_indices_wz);
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


        //****************************************************************************

        fe->reinit (elem);

        Ke.resize (dof_indices_wz.size(),
                   dof_indices_wz.size());

        Fe.resize (dof_indices_wz.size());

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

            //****************************************************************************

  //          double A=0.0;
            double RHS_wz=0.0;

            double px_cur = 0.0;
            double py_cur = 0.0;
            double pz_cur = 0.0;

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

                px_cur += dphi[l][qp](0)*system_p.current_solution(dof_indices_p[l]);
                py_cur += dphi[l][qp](1)*system_p.current_solution(dof_indices_p[l]);
                pz_cur += dphi[l][qp](2)*system_p.current_solution(dof_indices_p[l]);
            }


            FF_def(0,0)=F00_cur;
            FF_def(0,1)=F01_cur;
            FF_def(0,2)=F02_cur;
            FF_def(1,0)=F10_cur;
            FF_def(1,1)=F11_cur;
            FF_def(1,2)=F12_cur;
            FF_def(2,0)=F20_cur;
            FF_def(2,1)=F21_cur;
            FF_def(2,2)=F22_cur;

            FF_def_inv.zero();
            FF_def_inv_tr.zero();

            FF_def_inv=tensor_inverse(FF_def,NDIM);
            FF_def_inv_tr=tensor_inverse_transpose(FF_def,NDIM);


            for (std::size_t i=0; i<phi.size(); i++)
            {
                for (std::size_t j=0; j<phi.size(); j++)
                {
                    Ke(i,j) += JxW[qp]*phi[j][qp]*phi[i][qp];
                }
            }

            final_input=-(J_cur)*(FF_def_inv)*(Perm_Tensor)*(FF_def_inv_tr);
            RHS_wz=((final_input(2,0))*(px_cur))+((final_input(2,1))*(py_cur))+((final_input(2,2))*(pz_cur));

            for (std::size_t i=0; i<phi.size(); i++)
            {

                Fe(i) += JxW[qp]*(RHS_wz)*(phi[i][qp]);

            }
        }
        dof_map_wz.constrain_element_matrix_and_vector (Ke, Fe, dof_indices_wz);
        system_wz.matrix->add_matrix (Ke, dof_indices_wz);
        system_wz.rhs->add_vector    (Fe, dof_indices_wz);
    }


    return;

}



//*****************************************************************************

void DarcyVelocity::assemble_wxstar(EquationSystems& es,
                              const std::string& system_name)
{

    libmesh_assert_equal_to (system_name, wxstarSystemName);
    const MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    LinearImplicitSystem & system_wxstar = es.get_system<LinearImplicitSystem> (wxstarSystemName);
    LinearImplicitSystem & system_pstar = es.get_system<LinearImplicitSystem> (DarcyPressure::pstarSystemName);
    LinearImplicitSystem& J_system = es.get_system<LinearImplicitSystem>(MechanicsModel::J_system_name);
    LinearImplicitSystem& F_system = es.get_system<LinearImplicitSystem>(MechanicsModel::FF_system_name);

    const DofMap & dof_map_wxstar = system_wxstar.get_dof_map();
    const DofMap & dof_map_pstar = system_pstar.get_dof_map();
    const DofMap & dof_map_J = J_system.get_dof_map();
    const DofMap & dof_map_F = F_system.get_dof_map();

    //***************************************************************************
    FEType fe_type = dof_map_wxstar.variable_type(0);
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

    TensorValue<double> final_input; // The assembled RHS matrices
    TensorValue<double> FF_def; //FF_def contains the deformation gradient
    TensorValue<double> FF_def_inv; // This is the inverse of the deformation gradient
    TensorValue<double> FF_def_inv_tr; // The is the transpose of the deformation gradient inverse
//    double f_of_J;

    fe_face->attach_quadrature_rule (&qface);

    const std::vector<Real> & JxW = fe->get_JxW();
//    const std::vector<libMesh::Point> & q_point = fe->get_xyz();
    const std::vector<std::vector<Real>> & phi = fe->get_phi();
    const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    std::vector<dof_id_type> dof_indices_pstar;
    std::vector<dof_id_type> dof_indices_wxstar;
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


        dof_map_pstar.dof_indices (elem, dof_indices_pstar);
        dof_map_wxstar.dof_indices (elem, dof_indices_wxstar);
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


        //****************************************************************************

        fe->reinit (elem);

        Ke.resize (dof_indices_wxstar.size(),
                   dof_indices_wxstar.size());

        Fe.resize (dof_indices_wxstar.size());

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

            //****************************************************************************

//            double A=0.0;
            double RHS_wxstar=0.0;

            double pxstar_cur = 0.0;
            double pystar_cur = 0.0;
            double pzstar_cur = 0.0;

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

                pxstar_cur += dphi[l][qp](0)*system_pstar.current_solution(dof_indices_pstar[l]);
                pystar_cur += dphi[l][qp](1)*system_pstar.current_solution(dof_indices_pstar[l]);
                pzstar_cur += dphi[l][qp](2)*system_pstar.current_solution(dof_indices_pstar[l]);
            }


            FF_def(0,0)=F00_cur;
            FF_def(0,1)=F01_cur;
            FF_def(0,2)=F02_cur;
            FF_def(1,0)=F10_cur;
            FF_def(1,1)=F11_cur;
            FF_def(1,2)=F12_cur;
            FF_def(2,0)=F20_cur;
            FF_def(2,1)=F21_cur;
            FF_def(2,2)=F22_cur;

            FF_def_inv.zero();
            FF_def_inv_tr.zero();

            FF_def_inv=tensor_inverse(FF_def,NDIM);
            FF_def_inv_tr=tensor_inverse_transpose(FF_def,NDIM);


            for (std::size_t i=0; i<phi.size(); i++)
            {
                for (std::size_t j=0; j<phi.size(); j++)
                {
                    Ke(i,j) += JxW[qp]*phi[j][qp]*phi[i][qp];
                }
            }

            final_input=-(J_cur)*(FF_def_inv)*(Perm_Tensor)*(FF_def_inv_tr);
            RHS_wxstar=((final_input(0,0))*(pxstar_cur))+((final_input(0,1))*(pystar_cur))+((final_input(0,2))*(pzstar_cur));



            for (std::size_t i=0; i<phi.size(); i++)
            {
                Fe(i) += JxW[qp]*(RHS_wxstar)*(phi[i][qp]);
            }
        }
        dof_map_wxstar.constrain_element_matrix_and_vector (Ke, Fe, dof_indices_wxstar);
        system_wxstar.matrix->add_matrix (Ke, dof_indices_wxstar);
        system_wxstar.rhs->add_vector    (Fe, dof_indices_wxstar);
    }


    return;
}


//*****************************************************************************

void DarcyVelocity::assemble_wystar(EquationSystems& es,
                              const std::string& system_name)
{


    libmesh_assert_equal_to (system_name, wystarSystemName);
    const MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    LinearImplicitSystem & system_wystar = es.get_system<LinearImplicitSystem> (wystarSystemName);
    LinearImplicitSystem & system_pstar = es.get_system<LinearImplicitSystem> (DarcyPressure::pstarSystemName);
    LinearImplicitSystem& J_system = es.get_system<LinearImplicitSystem>(MechanicsModel::J_system_name);
    LinearImplicitSystem& F_system = es.get_system<LinearImplicitSystem>(MechanicsModel::FF_system_name);

    const DofMap & dof_map_wystar = system_wystar.get_dof_map();
    const DofMap & dof_map_pstar = system_pstar.get_dof_map();
    const DofMap & dof_map_J = J_system.get_dof_map();
    const DofMap & dof_map_F = F_system.get_dof_map();

    //***************************************************************************
    FEType fe_type = dof_map_wystar.variable_type(0);
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

    TensorValue<double> final_input; // The assembled RHS matrices
    TensorValue<double> FF_def; //FF_def contains the deformation gradient
    TensorValue<double> FF_def_inv; // This is the inverse of the deformation gradient
    TensorValue<double> FF_def_inv_tr; // The is the transpose of the deformation gradient inverse

    fe_face->attach_quadrature_rule (&qface);

    const std::vector<Real> & JxW = fe->get_JxW();
//    const std::vector<libMesh::Point> & q_point = fe->get_xyz();
    const std::vector<std::vector<Real>> & phi = fe->get_phi();
    const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    std::vector<dof_id_type> dof_indices_pstar;
    std::vector<dof_id_type> dof_indices_wystar;
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


        dof_map_pstar.dof_indices (elem, dof_indices_pstar);
        dof_map_wystar.dof_indices (elem, dof_indices_wystar);
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


        //****************************************************************************

        fe->reinit (elem);

        Ke.resize (dof_indices_wystar.size(),
                   dof_indices_wystar.size());

        Fe.resize (dof_indices_wystar.size());

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

            //****************************************************************************

//            double A=0.0;
            double RHS_wystar=0.0;

            double pxstar_cur = 0.0;
            double pystar_cur = 0.0;
            double pzstar_cur = 0.0;

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

                pxstar_cur += dphi[l][qp](0)*system_pstar.current_solution(dof_indices_pstar[l]);
                pystar_cur += dphi[l][qp](1)*system_pstar.current_solution(dof_indices_pstar[l]);
                pzstar_cur += dphi[l][qp](2)*system_pstar.current_solution(dof_indices_pstar[l]);
            }


            FF_def(0,0)=F00_cur;
            FF_def(0,1)=F01_cur;
            FF_def(0,2)=F02_cur;
            FF_def(1,0)=F10_cur;
            FF_def(1,1)=F11_cur;
            FF_def(1,2)=F12_cur;
            FF_def(2,0)=F20_cur;
            FF_def(2,1)=F21_cur;
            FF_def(2,2)=F22_cur;

            FF_def_inv.zero();
            FF_def_inv_tr.zero();

            FF_def_inv=tensor_inverse(FF_def,NDIM);
            FF_def_inv_tr=tensor_inverse_transpose(FF_def,NDIM);


            for (std::size_t i=0; i<phi.size(); i++)
            {
                for (std::size_t j=0; j<phi.size(); j++)
                {
                    Ke(i,j) += JxW[qp]*phi[j][qp]*phi[i][qp];
                }
            }

            final_input=-(J_cur)*(FF_def_inv)*(Perm_Tensor)*(FF_def_inv_tr);
            RHS_wystar=((final_input(1,0))*(pxstar_cur))+((final_input(1,1))*(pystar_cur))+((final_input(1,2))*(pzstar_cur));

            for (std::size_t i=0; i<phi.size(); i++)
            {

                Fe(i) += JxW[qp]*(RHS_wystar)*(phi[i][qp]);

            }
        }
        dof_map_wystar.constrain_element_matrix_and_vector (Ke, Fe, dof_indices_wystar);
        system_wystar.matrix->add_matrix (Ke, dof_indices_wystar);
        system_wystar.rhs->add_vector    (Fe, dof_indices_wystar);
    }


    return;
}


//*****************************************************************************

void DarcyVelocity::assemble_wzstar(EquationSystems& es,
                              const std::string& system_name)
{

    libmesh_assert_equal_to (system_name, wzstarSystemName);
    const MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    LinearImplicitSystem & system_wzstar = es.get_system<LinearImplicitSystem> (wzstarSystemName);
    LinearImplicitSystem & system_pstar = es.get_system<LinearImplicitSystem> (DarcyPressure::pstarSystemName);
    LinearImplicitSystem& J_system = es.get_system<LinearImplicitSystem>(MechanicsModel::J_system_name);
    LinearImplicitSystem& F_system = es.get_system<LinearImplicitSystem>(MechanicsModel::FF_system_name);

    const DofMap & dof_map_wzstar = system_wzstar.get_dof_map();
    const DofMap & dof_map_pstar = system_pstar.get_dof_map();
    const DofMap & dof_map_J = J_system.get_dof_map();
    const DofMap & dof_map_F = F_system.get_dof_map();

    //***************************************************************************
    FEType fe_type = dof_map_wzstar.variable_type(0);
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

    TensorValue<double> final_input; // The assembled RHS matrices
    TensorValue<double> FF_def; //FF_def contains the deformation gradient
    TensorValue<double> FF_def_inv; // This is the inverse of the deformation gradient
    TensorValue<double> FF_def_inv_tr; // The is the transpose of the deformation gradient inverse


    fe_face->attach_quadrature_rule (&qface);

    const std::vector<Real> & JxW = fe->get_JxW();
//    const std::vector<libMesh::Point> & q_point = fe->get_xyz();
    const std::vector<std::vector<Real>> & phi = fe->get_phi();
    const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    std::vector<dof_id_type> dof_indices_pstar;
    std::vector<dof_id_type> dof_indices_wzstar;
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


        dof_map_pstar.dof_indices (elem, dof_indices_pstar);
        dof_map_wzstar.dof_indices (elem, dof_indices_wzstar);
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


        //****************************************************************************

        fe->reinit (elem);

        Ke.resize (dof_indices_wzstar.size(),
                   dof_indices_wzstar.size());

        Fe.resize (dof_indices_wzstar.size());

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

            //****************************************************************************

//            double A=0.0;
            double RHS_wzstar=0.0;

            double pxstar_cur = 0.0;
            double pystar_cur = 0.0;
            double pzstar_cur = 0.0;

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

                pxstar_cur += dphi[l][qp](0)*system_pstar.current_solution(dof_indices_pstar[l]);
                pystar_cur += dphi[l][qp](1)*system_pstar.current_solution(dof_indices_pstar[l]);
                pzstar_cur += dphi[l][qp](2)*system_pstar.current_solution(dof_indices_pstar[l]);
            }


            FF_def(0,0)=F00_cur;
            FF_def(0,1)=F01_cur;
            FF_def(0,2)=F02_cur;
            FF_def(1,0)=F10_cur;
            FF_def(1,1)=F11_cur;
            FF_def(1,2)=F12_cur;
            FF_def(2,0)=F20_cur;
            FF_def(2,1)=F21_cur;
            FF_def(2,2)=F22_cur;

            FF_def_inv.zero();
            FF_def_inv_tr.zero();

            FF_def_inv=tensor_inverse(FF_def,NDIM);
            FF_def_inv_tr=tensor_inverse_transpose(FF_def,NDIM);


            for (std::size_t i=0; i<phi.size(); i++)
            {
                for (std::size_t j=0; j<phi.size(); j++)
                {
                    Ke(i,j) += JxW[qp]*phi[j][qp]*phi[i][qp];
                }
            }

            final_input=-(J_cur)*(FF_def_inv)*(Perm_Tensor)*(FF_def_inv_tr);
            RHS_wzstar=((final_input(2,0))*(pxstar_cur))+((final_input(2,1))*(pystar_cur))+((final_input(2,2))*(pzstar_cur));

            for (std::size_t i=0; i<phi.size(); i++)
            {

                Fe(i) += JxW[qp]*(RHS_wzstar)*(phi[i][qp]);

            }
        }
        dof_map_wzstar.constrain_element_matrix_and_vector (Ke, Fe, dof_indices_wzstar);
        system_wzstar.matrix->add_matrix (Ke, dof_indices_wzstar);
        system_wzstar.rhs->add_vector    (Fe, dof_indices_wzstar);
    }


    return;

}
