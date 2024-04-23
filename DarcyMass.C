#include <DarcyMass.h>

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
#include "DarcyPressure.h"
#include "DarcyVelocity.h"
#include "ConoraryFlowInterface.h"
#include "arteriesPD.h"

double DarcyMass::phi0;
double DarcyMass::biot_modulus;
double DarcyMass::kappa0;

std::string DarcyMass::mSystemName = "m_system";
std::string DarcyMass::mvariableName = "m_variable";
Order DarcyMass::mvariable_order = FIRST;

std::string DarcyMass::maverageSystemName = "maverage_system";
std::string DarcyMass::maveragevariableName = "maverage_variable";
Order DarcyMass::maveragevariable_order = FIRST;

std::string DarcyMass::mstarSystemName = "mstar_system";
std::string DarcyMass::mstarvariableName = "mstar_variable";
Order DarcyMass::mstarvariable_order = FIRST;

std::string DarcyMass::moldSystemName = "mold_system";
std::string DarcyMass::moldvariableName = "mold_variable";
Order DarcyMass::moldvariable_order = FIRST;

std::string DarcyMass::mextraSystemName = "mextra_system";
std::string DarcyMass::mextravariableName = "mextra_variable";
Order DarcyMass::mextravariable_order = FIRST;


double DarcyMass::time;
double DarcyMass::time_step;

int DarcyMass::m_system_num;
int DarcyMass::mold_system_num;

void DarcyMass::setup_mstar(EquationSystems* es)
{
    es->add_system<LinearImplicitSystem> (mstarSystemName);
    LinearImplicitSystem& system_mstar = es->get_system<LinearImplicitSystem>(mstarSystemName);
    system_mstar.add_variable(mstarvariableName, mstarvariable_order, LAGRANGE);
    return;
}


void DarcyMass::setup_maverage(EquationSystems* es)
{
    es->add_system<LinearImplicitSystem> (maverageSystemName);
    LinearImplicitSystem& system_maverage = es->get_system<LinearImplicitSystem>(maverageSystemName);
    system_maverage.add_variable(maveragevariableName, maveragevariable_order, LAGRANGE);
    return;
}





void DarcyMass::setup_mextra(EquationSystems* es)
{
    es->add_system<LinearImplicitSystem> (mextraSystemName);
    LinearImplicitSystem& system_mextra = es->get_system<LinearImplicitSystem>(mextraSystemName);
    system_mextra.add_variable(mextravariableName, mextravariable_order, LAGRANGE);
    return;
}



//*****************************************************************************

void DarcyMass::setup_m(EquationSystems* es)
{
    es->add_system<LinearImplicitSystem> (mSystemName);
    LinearImplicitSystem& system_m = es->get_system<LinearImplicitSystem>(mSystemName);
    system_m.add_variable(mvariableName, mvariable_order, LAGRANGE);

    es->add_system<LinearImplicitSystem> (moldSystemName);
    LinearImplicitSystem& system_mold = es->get_system<LinearImplicitSystem>(moldSystemName);
    system_mold.add_variable(moldvariableName, moldvariable_order, LAGRANGE);

    return;
}


//*****************************************************************************


void DarcyMass::initialize_mstar(EquationSystems* es)
{
    LinearImplicitSystem& system_mstar = es->get_system<LinearImplicitSystem>(mstarSystemName);
    system_mstar.attach_assemble_function(assemble_mstar);
    unsigned int system_num = system_mstar.number();

    libMesh::MeshBase& mesh = es->get_mesh();

    libMesh::NumericVector<double>& mstar_soln = *system_mstar.solution;

    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();

    double mstar_0 = 0.0;
    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;
        const unsigned int dof_num = nd->dof_number(system_num, 0,0);
        mstar_soln.set(dof_num, mstar_0);
    }
    mstar_soln.close();
    mstar_soln.localize(*system_mstar.current_local_solution);

    return;
}




void DarcyMass::initialize_maverage(EquationSystems* es)
{
    LinearImplicitSystem& system_maverage = es->get_system<LinearImplicitSystem>(maverageSystemName);
    system_maverage.attach_assemble_function(assemble_maverage);
    unsigned int system_num = system_maverage.number();

    libMesh::MeshBase& mesh = es->get_mesh();

    libMesh::NumericVector<double>& maverage_soln = *system_maverage.solution;

    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();

    double maverage_0 = 0.0;
    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;
        const unsigned int dof_num = nd->dof_number(system_num, 0,0);
        maverage_soln.set(dof_num, maverage_0);
    }
    maverage_soln.close();
    maverage_soln.localize(*system_maverage.current_local_solution);

    return;
}












void DarcyMass::initialize_mextra(EquationSystems* es)
{
    LinearImplicitSystem& system_mextra = es->get_system<LinearImplicitSystem>(mextraSystemName);
    unsigned int system_num = system_mextra.number();

    libMesh::MeshBase& mesh = es->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    const libMesh::DofMap& dof_map = system_mextra.get_dof_map();
    libMesh::NumericVector<double>& mextra_soln = *system_mextra.solution;

    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();

    double mext_0 = 0.0;
    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;
        const unsigned int dof_num = nd->dof_number(system_num, 0,0);
        mextra_soln.set(dof_num, mext_0);
    }
    mextra_soln.close();
    mextra_soln.localize(*system_mextra.current_local_solution);

    return;
}





//*****************************************************************************

void DarcyMass::initialize_m(EquationSystems* es)
{
    LinearImplicitSystem& system_m = es->get_system<LinearImplicitSystem>(mSystemName);
    system_m.attach_assemble_function(assemble_m);

    LinearImplicitSystem & system_mold = es->get_system<LinearImplicitSystem> (moldSystemName);

    unsigned int m_system_num = system_m.number();
    unsigned int m_old_system_num = system_mold.number();

    libMesh::MeshBase& mesh = es->get_mesh();


    libMesh::NumericVector<double>& m_soln = *system_m.solution;

    libMesh::NumericVector<double>& mold_soln = *system_mold.solution;



    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();

    double m_0 = 0.0;
    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;
        const unsigned int dof_num = nd->dof_number(m_system_num, 0,0);
        m_soln.set(dof_num, m_0);

        const unsigned int mold_dof_num = nd->dof_number(m_old_system_num, 0,0);
        mold_soln.set(mold_dof_num, m_0);

    }
    m_soln.close();
    m_soln.localize(*system_m.current_local_solution);

    mold_soln.close();
    mold_soln.localize(*system_mold.current_local_solution);

    return;
}

//*****************************************************************************


void DarcyMass::solve_mstar(EquationSystems* es,const double loop_time,const double dt)
{
    time=loop_time;
    time_step=dt;
    LinearImplicitSystem& system_mstar = es->get_system<LinearImplicitSystem>(mstarSystemName);
    system_mstar.solve();

    return;
}

//*****************************************************************************


void DarcyMass::solve_maverage(EquationSystems* es,const double loop_time,const double dt)
{
    time=loop_time;
    time_step=dt;
    LinearImplicitSystem& system_maverage = es->get_system<LinearImplicitSystem>(maverageSystemName);
    system_maverage.solve();
    Tube::average_added_mass_sum = system_maverage.solution->sum();
    return;
}

//*****************************************************************************



void DarcyMass::solve_mextra(EquationSystems* es,const double loop_time,const double dt)
{
    time=loop_time;
    time_step=dt;
    LinearImplicitSystem& system_mextra = es->get_system<LinearImplicitSystem>(mextraSystemName);
    system_mextra.solve();

    return;
}

//*****************************************************************************






void DarcyMass::solve_m(EquationSystems* es,const double loop_time,const double dt)
{
    time=loop_time;
    time_step=dt;
    LinearImplicitSystem& system_m = es->get_system<LinearImplicitSystem>(mSystemName);

    LinearImplicitSystem& system_mold = es->get_system<LinearImplicitSystem>(moldSystemName);
    libMesh::NumericVector<double>& mold_solution = *(system_mold.solution);

    libMesh::NumericVector<double>& m_cur_solution = *(system_m.solution);
    m_cur_solution.close();
    m_cur_solution.localize(*system_m.current_local_solution);


    libMesh::MeshBase& mesh = es->get_mesh();
    libMesh::MeshBase::const_node_iterator node_it = mesh.local_nodes_begin();
    libMesh::MeshBase::const_node_iterator node_end = mesh.local_nodes_end();

    unsigned int m_system_num = system_m.number();
    unsigned int m_old_system_num = system_mold.number();

    for(; node_it != node_end; ++node_it)
    {
      const libMesh::Node* nd = *node_it;
      const unsigned int dof_num = nd->dof_number(m_system_num, 0, 0);
      const unsigned int mold_dof_num = nd->dof_number(m_old_system_num, 0, 0);

      double mold_t = mold_solution(mold_dof_num);
      double mcur_t = m_cur_solution(dof_num);

      mold_solution.set(mold_dof_num, mcur_t);

    }

   mold_solution.close();
   mold_solution.localize(*system_mold.current_local_solution);


   system_m.solve();

Tube::added_mass_sum = system_m.solution->sum();

    return;
}

//*****************************************************************************

void DarcyMass::assemble_mstar(EquationSystems& es,
                             const std::string& system_name)
{

    libmesh_assert_equal_to (system_name, mstarSystemName);
    const MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    LinearImplicitSystem & system_mstar = es.get_system<LinearImplicitSystem> (mstarSystemName);
    LinearImplicitSystem & system_pstar = es.get_system<LinearImplicitSystem> (DarcyPressure::pstarSystemName);
    LinearImplicitSystem& J_system = es.get_system<LinearImplicitSystem>(MechanicsModel::J_system_name);
    LinearImplicitSystem & system_m = es.get_system<LinearImplicitSystem> (mSystemName);
    LinearImplicitSystem & system_wx = es.get_system<LinearImplicitSystem> (DarcyVelocity::wxstarSystemName);
    LinearImplicitSystem & system_wy = es.get_system<LinearImplicitSystem> (DarcyVelocity::wystarSystemName);
    LinearImplicitSystem & system_wz = es.get_system<LinearImplicitSystem> (DarcyVelocity::wzstarSystemName);
    LinearImplicitSystem& system_source = es.get_system<LinearImplicitSystem>(sourceToDarcy::sourceToDarcySystemName);

    const DofMap & dof_map_mstar = system_mstar.get_dof_map();
    const DofMap & dof_map_pstar = system_pstar.get_dof_map();
    const DofMap & dof_map_J = J_system.get_dof_map();
    const DofMap & dof_map_m = system_m.get_dof_map();
const DofMap & dof_map_wx = system_wx.get_dof_map();
const DofMap & dof_map_wy = system_wy.get_dof_map();
const DofMap & dof_map_wz = system_wz.get_dof_map();
const DofMap & dof_map_source = system_source.get_dof_map();
    //***************************************************************************
    FEType fe_type = dof_map_mstar.variable_type(0);
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

    std::vector<dof_id_type> dof_indices_pstar;
    std::vector<dof_id_type> dof_indices_mstar;
    std::vector<dof_id_type> dof_indices_J;
    std::vector<dof_id_type> dof_indices_m;
    std::vector<dof_id_type> dof_indices_wx;
    std::vector<dof_id_type> dof_indices_wy;
    std::vector<dof_id_type> dof_indices_wz;
    std::vector<dof_id_type> dof_indices_source;

    //****************************************************************************

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    for (;  el != end_el; ++el)
    {

        Elem* elem = *el;

        dof_map_pstar.dof_indices (elem, dof_indices_pstar);
        dof_map_mstar.dof_indices (elem, dof_indices_mstar);
        dof_map_J.dof_indices (elem, dof_indices_J);
        dof_map_m.dof_indices (elem, dof_indices_m);
        dof_map_wx.dof_indices (elem, dof_indices_wx);
        dof_map_wy.dof_indices (elem, dof_indices_wy);
        dof_map_wz.dof_indices (elem, dof_indices_wz);
        dof_map_source.dof_indices (elem, dof_indices_source);

        //****************************************************************************

        fe->reinit (elem);

        Ke.resize (dof_indices_mstar.size(),
                   dof_indices_mstar.size());

        Fe.resize (dof_indices_mstar.size());

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
            //****************************************************************************

            double m_cur = 0.0;
            double mstar_cur = 0.0;
            double pstar_cur = 0.0;
            double J_cur = 0.0;
            double RHS_mstar = 0.0;
            double f_of_j = 0.0;
            double RHS_extra=0.0;
            double source_cur=0.0;
            double penalty=0.0;
            double wx_cur=0.0;
            double wy_cur=0.0;
            double wz_cur=0.0;

            for (unsigned int l=0; l<phi.size(); l++)
            {
                J_cur += phi[l][qp]*J_system.current_solution(dof_indices_J[l]);
                pstar_cur += phi[l][qp]*system_pstar.current_solution(dof_indices_pstar[l]);
                m_cur += phi[l][qp]*system_m.current_solution(dof_indices_m[l]);
                mstar_cur += phi[l][qp]*system_mstar.current_solution(dof_indices_mstar[l]);
                source_cur += phi[l][qp]*system_source.current_solution(dof_indices_source[l]);
                wx_cur += dphi[l][qp](0)*system_wx.current_solution(dof_indices_wx[l]);
                wy_cur += dphi[l][qp](1)*system_wy.current_solution(dof_indices_wy[l]);
                wz_cur += dphi[l][qp](2)*system_wz.current_solution(dof_indices_wz[l]);
            }

            if (J_cur==1.0)
            {
             	f_of_j=1;
            }
            else
            {
             	f_of_j=(2.0*(J_cur-1.0-log(J_cur)))/((J_cur-1.0)*(J_cur-1.0));
            }

            if (m_cur+phi0<0.0001)
            {
                penalty=0.0;
            }
            else
            {
                penalty=(kappa0)/(m_cur+phi0);
            }


            for (std::size_t i=0; i<phi.size(); i++)
            {
                for (std::size_t j=0; j<phi.size(); j++)
                {
                    Ke(i,j) += JxW[qp]*phi[j][qp]*phi[i][qp];
                }
            }

  RHS_mstar=((pstar_cur+penalty)/(biot_modulus*f_of_j))+J_cur-1.0;

//  RHS_mstar=m_cur+((0.5*time_step)*(source_cur-wx_cur-wy_cur-wz_cur));



            for (std::size_t i=0; i<phi.size(); i++)
            {
                Fe(i) += JxW[qp]*(RHS_mstar)*(phi[i][qp]);
            }
        }
        dof_map_mstar.constrain_element_matrix_and_vector (Ke, Fe, dof_indices_mstar);
        system_mstar.matrix->add_matrix (Ke, dof_indices_mstar);
        system_mstar.rhs->add_vector    (Fe, dof_indices_mstar);
    }


    return;


}

//*****************************************************************************

void DarcyMass::assemble_m(EquationSystems& es,
                             const std::string& system_name)
{

    libmesh_assert_equal_to (system_name, mSystemName);
    const MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    LinearImplicitSystem & system_mstar = es.get_system<LinearImplicitSystem> (mstarSystemName);
    LinearImplicitSystem & system_m = es.get_system<LinearImplicitSystem> (mSystemName);
    LinearImplicitSystem & system_wx = es.get_system<LinearImplicitSystem> (DarcyVelocity::wxSystemName);
    LinearImplicitSystem & system_wy = es.get_system<LinearImplicitSystem> (DarcyVelocity::wySystemName);
    LinearImplicitSystem & system_wz = es.get_system<LinearImplicitSystem> (DarcyVelocity::wzSystemName);
    LinearImplicitSystem& system_source = es.get_system<LinearImplicitSystem>(sourceToDarcy::sourceToDarcySystemName);


    const DofMap & dof_map_mstar = system_mstar.get_dof_map();
    const DofMap & dof_map_m = system_m.get_dof_map();
    const DofMap & dof_map_wx = system_wx.get_dof_map();
    const DofMap & dof_map_wy = system_wy.get_dof_map();
    const DofMap & dof_map_wz = system_wz.get_dof_map();
    const DofMap & dof_map_source = system_source.get_dof_map();

    //***************************************************************************
    FEType fe_type = dof_map_m.variable_type(0);
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

    std::vector<dof_id_type> dof_indices_m;
    std::vector<dof_id_type> dof_indices_mstar;
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

        dof_map_wx.dof_indices (elem, dof_indices_wx);
        dof_map_wy.dof_indices (elem, dof_indices_wy);
        dof_map_wz.dof_indices (elem, dof_indices_wz);
        dof_map_m.dof_indices (elem, dof_indices_m);
        dof_map_mstar.dof_indices (elem, dof_indices_mstar);
        dof_map_source.dof_indices (elem, dof_indices_source);
        //****************************************************************************

        fe->reinit (elem);

        Ke.resize (dof_indices_m.size(),
                   dof_indices_m.size());

        Fe.resize (dof_indices_m.size());

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
            //****************************************************************************

            double mstar_cur = 0.0;
            double RHS_m = 0.0;
            double wx_cur=0.0;
            double wy_cur=0.0;
            double wz_cur=0.0;
            double RHS_extra=0.0;
            double source_cur=0.0;
            double source_other=0.0;

            for (unsigned int l=0; l<phi.size(); l++)
            {
                mstar_cur += phi[l][qp]*system_mstar.current_solution(dof_indices_mstar[l]);
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

            RHS_extra=source_cur-wx_cur-wy_cur-wz_cur;
            RHS_m=mstar_cur+((0.5*time_step)*(RHS_extra));


            for (std::size_t i=0; i<phi.size(); i++)
            {
                Fe(i) += JxW[qp]*(RHS_m)*(phi[i][qp]);
            }
        }

        dof_map_m.constrain_element_matrix_and_vector (Ke, Fe, dof_indices_m);
        system_m.matrix->add_matrix (Ke, dof_indices_m);
        system_m.rhs->add_vector    (Fe, dof_indices_m);
    }


    return;


}








//*****************************************************************************

void DarcyMass::assemble_maverage(EquationSystems& es,
                             const std::string& system_name)
{

    libmesh_assert_equal_to (system_name, mSystemName);
    const MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    LinearImplicitSystem & system_m = es.get_system<LinearImplicitSystem> (mSystemName);
    LinearImplicitSystem & system_maverage = es.get_system<LinearImplicitSystem> (maverageSystemName);



    const DofMap & dof_map_m = system_m.get_dof_map();
    const DofMap & dof_map_maverage = system_maverage.get_dof_map();


    //***************************************************************************
    FEType fe_type = dof_map_maverage.variable_type(0);
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

    std::vector<dof_id_type> dof_indices_maverage;
    std::vector<dof_id_type> dof_indices_m;
    //****************************************************************************



    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    for (;  el != end_el; ++el)
    {

        Elem* elem = *el;


        dof_map_maverage.dof_indices (elem, dof_indices_maverage);
        dof_map_m.dof_indices (elem, dof_indices_m);
        //****************************************************************************

        fe->reinit (elem);

        Ke.resize (dof_indices_maverage.size(),
                   dof_indices_maverage.size());

        Fe.resize (dof_indices_maverage.size());


        double element_volume=elem->volume();


        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
            //****************************************************************************

            double m_cur = 0.0;
            double RHS_m = 0.0;

            for (unsigned int l=0; l<phi.size(); l++)
            {
                m_cur += phi[l][qp]*system_m.current_solution(dof_indices_m[l]);
            }

//---------------------------------------------------------------------------------------

            for (std::size_t i=0; i<phi.size(); i++)
            {
                for (std::size_t j=0; j<phi.size(); j++)
                {
                    Ke(i,j) += JxW[qp]*phi[j][qp]*phi[i][qp];
                }
            }


            RHS_m=(element_volume)*(m_cur);

            for (std::size_t i=0; i<phi.size(); i++)
            {
                Fe(i) += JxW[qp]*(RHS_m)*(phi[i][qp]);
            }
        }

        dof_map_maverage.constrain_element_matrix_and_vector (Ke, Fe, dof_indices_maverage);
        system_maverage.matrix->add_matrix (Ke, dof_indices_maverage);
        system_maverage.rhs->add_vector    (Fe, dof_indices_maverage);
    }


    return;


}


















//*****************************************************************************

void DarcyMass::assemble_mextra(EquationSystems& es,
                             const std::string& system_name)
{

    libmesh_assert_equal_to (system_name, mextraSystemName);
    const MeshBase & mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    LinearImplicitSystem & system_mextra = es.get_system<LinearImplicitSystem> (mextraSystemName);
    LinearImplicitSystem& system_source = es.get_system<LinearImplicitSystem>(sourceToDarcy::sourceToDarcySystemName);

    const DofMap & dof_map_mextra = system_mextra.get_dof_map();
const DofMap & dof_map_source = system_source.get_dof_map();
    //***************************************************************************
    FEType fe_type = dof_map_mextra.variable_type(0);
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

    std::vector<dof_id_type> dof_indices_mextra;
    std::vector<dof_id_type> dof_indices_source;

    //****************************************************************************

    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    for (;  el != end_el; ++el)
    {

        Elem* elem = *el;

        dof_map_mextra.dof_indices (elem, dof_indices_mextra);
        dof_map_source.dof_indices (elem, dof_indices_source);

        //****************************************************************************

        fe->reinit (elem);

        Ke.resize (dof_indices_mextra.size(),
                   dof_indices_mextra.size());

        Fe.resize (dof_indices_mextra.size());

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
            //****************************************************************************

            double mextra_cur = 0.0;
            double RHS_mextra=0.0;
            double source_cur=0.0;


            for (unsigned int l=0; l<phi.size(); l++)
            {
                source_cur += phi[l][qp]*system_source.current_solution(dof_indices_source[l]);
                mextra_cur += phi[l][qp]*system_mextra.current_solution(dof_indices_mextra[l]);
            }




            for (std::size_t i=0; i<phi.size(); i++)
            {
                for (std::size_t j=0; j<phi.size(); j++)
                {
                    Ke(i,j) += JxW[qp]*phi[j][qp]*phi[i][qp];
                }
            }

            RHS_mextra=mextra_cur+((time_step)*(source_cur));

            for (std::size_t i=0; i<phi.size(); i++)
            {
                Fe(i) += JxW[qp]*(RHS_mextra)*(phi[i][qp]);
            }
        }
        dof_map_mextra.constrain_element_matrix_and_vector (Ke, Fe, dof_indices_mextra);
        system_mextra.matrix->add_matrix (Ke, dof_indices_mextra);
        system_mextra.rhs->add_vector    (Fe, dof_indices_mextra);
    }


    return;


}
