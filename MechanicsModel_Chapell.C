// Copyright (c) 2011-2013, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
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
// POSSIBILITY OF SUCH DAMAGE.

// IBAMR INCLUDES
#include <ibtk/libmesh_utilities.h>
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"


// APPLICATION INCLUDES
#include <MechanicsModel.h>
#include <ActiveContraction.h>

// STATIC VARIABLES
int MechanicsModel::f0_system_num, MechanicsModel::s0_system_num;
bool MechanicsModel::enable_active_tension, MechanicsModel::normalize_stress;
double MechanicsModel::T_scale, MechanicsModel::beta_s;
double MechanicsModel::I1_dev_max;
double MechanicsModel::I1_dev_min;
double MechanicsModel::J_dev_max;
double MechanicsModel::J_dev_min;
double MechanicsModel::I1_dil_max;
double MechanicsModel::I1_dil_min;
double MechanicsModel::J_dil_max;
double MechanicsModel::J_dil_min;





//this is for J and FF systems
int MechanicsModel::J_system_num, MechanicsModel::FF_system_num;
std::string MechanicsModel::J_system_name = "J_detF";
std::string MechanicsModel::FF_system_name = "FF_darcy";

std::string MechanicsModel::dJdt_system_name = "dJdt";
int MechanicsModel::dJdt_system_num;

std::string MechanicsModel::Jold_system_name = "Jold";
int MechanicsModel::Jold_system_num;


namespace // private namespace
{
    // Material parameters (from Chapelle et al. 2010).
    static const double b  = 1.0;  // nondimensiona
    static const double rho0 = 1.0; // g/cm^2
    static const double kappa_1 = 20000.0;
    static const double kappa_2 = 330.0;
    static const double cap_K= 2200000.0;
    static const double M=2180000.0;
    static const double kappa_0 = 0.1;
    static const double phi_0 = 0.15;


   static const double a_coefficient=1.0; //for parameter sweeping


   static const double a_HO   = a_coefficient*2244.87;  // dyne/cm^2
   static const double b_HO   = 1.6215;  // nondimensional
   static const double af_HO  = a_coefficient*2.4267e4;  // dyne/cm^2
   static const double bf_HO  = 1.8268;  // nondimensional
   static const double as_HO  = a_coefficient*5562.38;  // dyne/cm^2
   static const double bs_HO  = 0.7746;  // nondimensional
   static const double afs_HO = a_coefficient*3905.16;  // dyne/cm^2
   static const double bfs_HO = 1.695;  // nondimensional

//   static const double a_HO   = a_coefficient*22448.7;  // dyne/cm^2
//   static const double b_HO   = 1.6215;  // nondimensional
//   static const double af_HO  = a_coefficient*24267.0;  // dyne/cm^2
//   static const double bf_HO  = 1.8268;  // nondimensional
//   static const double as_HO  = a_coefficient*5562.38;  // dyne/cm^2
//   static const double bs_HO  = 0.7746;  // nondimensional
//   static const double afs_HO = a_coefficient*3905.16;  // dyne/cm^2
//   static const double bfs_HO = 1.695;  // nondimensional
}

// CLASS IMPLEMENTATION

void
MechanicsModel::PK1_dev_stress_function(
                                        TensorValue<double>& PP,
                                        const TensorValue<double>& FF,
                                        const libMesh::Point& /* X */,  // current   location
                                        const libMesh::Point& /* s */,  // reference location
                                        Elem* const elem,
                                        const vector<const vector<double>*>& var_data,
                                        const vector<const vector<libMesh::VectorValue<double> >*>& /*system_grad_var_data*/,
                                        double /* data_time */,
                                        void* /* ctx */)
{
    const TensorValue<double> CC = FF.transpose() * FF;
    const TensorValue<double> BB = FF * FF.transpose();
    const TensorValue<double> FF_inv_trans = tensor_inverse_transpose(FF,NDIM);

    const TensorValue<double> II(1.0, 0.0, 0.0, 0.0, 1.0,0.0, 0.0, 0.0, 1.0);


    const double I1 = CC.tr();
    const double J = FF.det();

    I1_dev_max = max(I1_dev_max, I1);
    I1_dev_min = min(I1_dev_min, I1);

    J_dev_max = max(J_dev_max, J);
    J_dev_min = min(J_dev_min, J);

    const std::vector<double>& f0_var = *var_data[0];
    const std::vector<double>& s0_var = *var_data[1];
    const VectorValue<double> f0(f0_var[0], f0_var[1], f0_var[2]);
    const VectorValue<double> s0(s0_var[0], s0_var[1], s0_var[2]);

     PP = a_HO*exp(b_HO*(I1-3.0))*FF;

    //fibre contribution
    const double I4f = f0*(CC*f0);
    if (I4f > 1.0)
    {
        PP += 2.0*af_HO*(I4f-1.0)*exp(bf_HO*(I4f-1.0)*(I4f-1.0))*FF*outer_product(f0,f0);
    }


    // sheet contribution
    const double I4s = s0*(CC*s0);
    if (I4s > 1.0)
    {
        PP += 2.0*as_HO*(I4s-1.0)*exp(bs_HO*(I4s-1.0)*(I4s-1.0))*FF*outer_product(s0,s0);
    }


    // Fiber-sheet contribution.
    const double I8fs = f0*(CC*s0);
    PP += afs_HO*I8fs*exp(bfs_HO*I8fs*I8fs)*FF*(outer_product(f0,s0)+outer_product(s0,f0));



    // Active contribution
    if (enable_active_tension)
    {
      const std::vector<double>& acT_var = *var_data[2];
      const double acT = acT_var[0];

// here acT is redefined to 0 less than 0.4, 60*((t-0.4)/(0.2)) between 0.4 and 0.6 and 60.0 afterwards

      const double T = T_scale*acT*1.0e4;  // NOTE: T is stored in kPa; must be converted to dyne/cm^2

      if (T > 0.0)
      {
          PP += J*T*FF*outer_product(f0,f0);
      }
    }



    return;
}// PK1_dev_stress_function

void
MechanicsModel::PK1_dil_stress_function(
                                        TensorValue<double>& PP,
                                        const TensorValue<double>& FF,
                                        const libMesh::Point& /* X */,  // current   location
                                        const libMesh::Point& /* s */,  // reference location
                                        Elem* const /*elem*/,
                                        const vector<const vector<double>*>& system_data,
                                        const vector<const vector<libMesh::VectorValue<double> >*>& /*system_grad_var_data*/,
                                        double /* data_time */,
                                        void* /* ctx */)
{
    const TensorValue<double> CC = FF.transpose() * FF;
    const TensorValue<double> FF_inv_trans = tensor_inverse_transpose(FF,NDIM);
    const TensorValue<double> BB = FF * FF.transpose();

    const double I1 = CC.tr();
    const double J = FF.det();

    I1_dil_max = max(I1_dil_max, I1);
    I1_dil_min = min(I1_dil_min, I1);

    J_dil_max = max(J_dil_max, J);
    J_dil_min = min(J_dil_min, J);

    const std::vector<double>& f0_var = *system_data[2];
    const std::vector<double>& s0_var = *system_data[3];
    const VectorValue<double> f0(f0_var[0], f0_var[1], f0_var[2]);
    const VectorValue<double> s0(s0_var[0], s0_var[1], s0_var[2]);


    PP.zero();

   //this term is to ensure when F=I, PP = 0
    PP = -a_HO*exp(b_HO*(I1-3.0))*FF_inv_trans;

    //access the M system from Darcy Solver
    const std::vector<double> & M_vec = *system_data[1];
    const double M_darcy = M_vec[0];

    // decide whether to include the penalty term on J - m/rho0
PP += (  (2.0*cap_K*J*log(J-M_darcy))/(J-M_darcy)  )*FF_inv_trans;


    return;
}// PK0_dil_stress_function













void MechanicsModel::assemble_J_system( EquationSystems& es,
                                       const std::string& libmesh_dbg_var(system_name))
{
    std::cout<< "assemble J_system\n";

    LinearImplicitSystem& J_system = es.get_system<LinearImplicitSystem>(J_system_name);
  //  unsigned int J_system_num = J_system.number();

    libMesh::MeshBase& mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();


    //access the COORD_MAPPING_SYSTEM
    System& dX_system = es.get_system<System>(IBFEMethod::COORD_MAPPING_SYSTEM_NAME);
    dX_system.solution->localize(*dX_system.current_local_solution);
    NumericVector<double>& dX_data = *(dX_system.current_local_solution);
    dX_data.close();
    //for displacement vector
    std::vector<int> vars(NDIM);
    for (unsigned int d=0; d < NDIM; ++d) vars[d] = d;
    const DofMap & dX_dof_map = dX_system.get_dof_map();

    const libMesh::DofMap& dof_map = J_system.get_dof_map();
    FEType fe_type = dof_map.variable_type(0);
    UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
    QGauss qrule (dim, FIFTH);
    fe->attach_quadrature_rule(&qrule);

    UniquePtr<FEBase> fe_face (FEBase::build(dim,fe_type));
    QGauss qface(dim-1, FIFTH);
    fe_face->attach_quadrature_rule(&qface);

//    libMesh::NumericVector<double>& J_soln = *(J_system.solution);

    const std::vector<Real>& JxW = fe->get_JxW();
//    const std::vector<libMesh::Point>& q_point = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();


    TensorValue<double> II(1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0);
    TensorValue<double> FF;
    double FFdet = 0.0;
    std::vector<std::vector<unsigned int> > dof_indices_dX(NDIM);
    std::vector<dof_id_type> dof_indices_dx;
    std::vector<dof_id_type> dof_indices_dy;
    std::vector<dof_id_type> dof_indices_dz;




    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    std::vector<dof_id_type> dof_indices;

    MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el =  mesh.active_local_elements_end();

    for (; el != end_el; ++el)
    {
        Elem* elem = *el;
        dof_map.dof_indices (elem, dof_indices);

        dX_dof_map.dof_indices(elem, dof_indices_dx, 0); // for dx
        dX_dof_map.dof_indices(elem, dof_indices_dy, 1); // for dy
        dX_dof_map.dof_indices(elem, dof_indices_dz, 2); // for dz



        fe->reinit(elem);

        Ke.resize (dof_indices.size(),
                   dof_indices.size());
        Fe.resize (dof_indices.size());

        for(unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
            for (unsigned int i = 0; i<phi.size(); i++)
                for(unsigned int j =0; j<phi.size(); j++)
                {
                    Ke(i,j) += JxW[qp]*phi[i][qp]*phi[j][qp];
                }


            Gradient grad_dx, grad_dy, grad_dz;
            grad_dx.zero();
            grad_dy.zero();
            grad_dz.zero();

            for (unsigned int i=0; i<dof_indices_dx.size(); i++)
            {
                double dx = dX_data(dof_indices_dx[i]);
                double dy = dX_data(dof_indices_dy[i]);
                double dz = dX_data(dof_indices_dz[i]);

                grad_dx.add_scaled(dphi[i][qp], dx);
                grad_dy.add_scaled(dphi[i][qp], dy);
                grad_dz.add_scaled(dphi[i][qp], dz);
            }
            //calculating FF
            FF.zero();
            FF(0,0) = grad_dx(0);
            FF(0,1) = grad_dx(1);
            FF(0,2) = grad_dx(2);

            FF(1,0) = grad_dy(0);
            FF(1,1) = grad_dy(1);
            FF(1,2) = grad_dy(2);

            FF(2,0) = grad_dz(0);
            FF(2,1) = grad_dz(1);
            FF(2,2) = grad_dz(2);
            FF = FF + II;
            FFdet = FF.det();

            for (std::size_t i = 0; i<phi.size(); i++)
            {
                Fe(i) += JxW[qp]*phi[i][qp]*FFdet;
            }

        }




        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
        J_system.matrix->add_matrix(Ke,dof_indices);
        J_system.rhs->add_vector   (Fe, dof_indices);


    } //end for el

    return;
}


void MechanicsModel::initialize_J_system(EquationSystems* es)
{
    LinearImplicitSystem& J_system = es->get_system<LinearImplicitSystem>(MechanicsModel::J_system_name);
    std::cout << MechanicsModel::J_system_name<<std::endl;
    J_system.attach_assemble_function(MechanicsModel::assemble_J_system);


    LinearImplicitSystem& dJdt_system = es->get_system<LinearImplicitSystem>(MechanicsModel::dJdt_system_name);
    LinearImplicitSystem& Jold_system = es->get_system<LinearImplicitSystem>(MechanicsModel::Jold_system_name);
    std::cout << MechanicsModel::dJdt_system_name<<std::endl;
    std::cout << MechanicsModel::Jold_system_name<<std::endl;


    //for initialization
    unsigned int system_num = J_system.number();

    libMesh::MeshBase& mesh = es->get_mesh();
//    const unsigned int dim = mesh.mesh_dimension();

//    const libMesh::DofMap& dof_map = J_system.get_dof_map();
    libMesh::NumericVector<double>& soln = *(J_system.solution);

    libMesh::NumericVector<double>& dJdt_soln = *(dJdt_system.solution);
    libMesh::NumericVector<double>& Jold_soln = *(Jold_system.solution);

    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();


    double JTemp = 1.0;
    double dJdtTemp = 0.0;
    for (; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;
        const unsigned int dof_num = nd->dof_number(system_num, 0, 0);
        soln.set(dof_num, JTemp);

        const unsigned int dJdt_dof_num = nd->dof_number(dJdt_system_num, 0, 0);
        dJdt_soln.set(dJdt_dof_num, dJdtTemp);

        const unsigned int Jold_dof_num = nd->dof_number(Jold_system_num, 0, 0);
        Jold_soln.set(Jold_dof_num, JTemp);

    }
    soln.close();
    soln.localize(*J_system.current_local_solution);

    dJdt_soln.close();
    dJdt_soln.localize(*dJdt_system.current_local_solution);

    Jold_soln.close();
    Jold_soln.localize(*Jold_system.current_local_solution);

    return;

}


void MechanicsModel::update_J_system(EquationSystems* es, double dt)
{
    std::cout<< "updating J_system\n";
    LinearImplicitSystem& J_system = es->get_system<LinearImplicitSystem>(J_system_name);

    LinearImplicitSystem& Jold_system = es->get_system<LinearImplicitSystem>(Jold_system_name);
    libMesh::NumericVector<double>& J_old_solution = *(Jold_system.solution);

    libMesh::NumericVector<double>& J_cur_solution = *(J_system.solution);
    J_cur_solution.close();
    J_cur_solution.localize(*J_system.current_local_solution);

    // now calculate the dJ/dt
    LinearImplicitSystem& dJdt_system = es->get_system<LinearImplicitSystem>(dJdt_system_name);
    libMesh::NumericVector<double>& dJdt_soln = *(dJdt_system.solution);

    libMesh::MeshBase& mesh = es->get_mesh();
    libMesh::MeshBase::const_node_iterator node_it = mesh.local_nodes_begin();
    libMesh::MeshBase::const_node_iterator node_end = mesh.local_nodes_end();

    for(; node_it != node_end; ++node_it)
    {
      const libMesh::Node* nd = *node_it;
      const unsigned int dof_num = nd->dof_number(J_system_num, 0, 0);
      const unsigned int Jold_dof_num = nd->dof_number(Jold_system_num, 0, 0);
      const unsigned int dJdt_dof_num = nd->dof_number(dJdt_system_num, 0, 0);

      double Jold_t = J_old_solution(Jold_dof_num);
      double Jcur_t = J_cur_solution(dof_num);

      if (abs(dt) < 1.0e-9) //virtual zero
          dJdt_soln.set(dJdt_dof_num, (Jcur_t - Jold_t) );
      else
          dJdt_soln.set(dJdt_dof_num, (Jcur_t - Jold_t)/dt );

     //update Jold solution
     J_old_solution.set(Jold_dof_num, Jcur_t);

    }

   dJdt_soln.close();
   dJdt_soln.localize(*dJdt_system.current_local_solution);

   J_old_solution.close();
   J_old_solution.localize(*Jold_system.current_local_solution);

   J_system.solve();

    return;
}


void MechanicsModel::initialize_FF_system(EquationSystems* es)
{
    LinearImplicitSystem& FF_system = es->get_system<LinearImplicitSystem>(MechanicsModel::FF_system_name);
    std::cout << "intialize FF to identity matrix "<< MechanicsModel::FF_system_name<<std::endl;
    FF_system.attach_assemble_function(MechanicsModel::assemble_FF_system);

    unsigned int system_num = FF_system.number();
    libMesh::MeshBase& mesh = es->get_mesh();
//    const unsigned int dim = mesh.mesh_dimension();

//    const libMesh::DofMap& dof_map = FF_system.get_dof_map();
    libMesh::NumericVector<double>& soln = *(FF_system.solution);

    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();

    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;
        const unsigned int dof_num_00 = nd->dof_number(system_num, 0, 0);
        const unsigned int dof_num_01 = nd->dof_number(system_num, 1, 0);
        const unsigned int dof_num_02 = nd->dof_number(system_num, 2, 0);

        const unsigned int dof_num_10 = nd->dof_number(system_num, 3, 0);
        const unsigned int dof_num_11 = nd->dof_number(system_num, 4, 0);
        const unsigned int dof_num_12 = nd->dof_number(system_num, 5, 0);

        const unsigned int dof_num_20 = nd->dof_number(system_num, 6, 0);
        const unsigned int dof_num_21 = nd->dof_number(system_num, 7, 0);
        const unsigned int dof_num_22 = nd->dof_number(system_num, 8, 0);

        soln.set(dof_num_00, 1.0);
        soln.set(dof_num_01, 0.0);
        soln.set(dof_num_02, 0.0);

        soln.set(dof_num_10, 0.0);
        soln.set(dof_num_11, 1.0);
        soln.set(dof_num_12, 0.0);

        soln.set(dof_num_20, 0.0);
        soln.set(dof_num_21, 0.0);
        soln.set(dof_num_22, 1.0);
    }
    soln.close();
    soln.localize(*FF_system.current_local_solution);
    return;
}


void MechanicsModel::assemble_FF_system(EquationSystems& es,
                                        const std::string& system_name)
{
    std::cout << "assemble FF system\n";
    LinearImplicitSystem& FF_system = es.get_system<LinearImplicitSystem>(MechanicsModel::FF_system_name);
//    unsigned int FF_system_num = FF_system.number();

    libMesh::MeshBase& mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    //access the COORD_MAPPING_SYSTEM
    System& dX_system = es.get_system<System>(IBFEMethod::COORD_MAPPING_SYSTEM_NAME);
    dX_system.solution->localize(*dX_system.current_local_solution);
    NumericVector<double>& dX_data = *(dX_system.current_local_solution);
    dX_data.close();
    //for displacement vector
    std::vector<int> vars_dX(NDIM);
    for (unsigned int d=0; d < NDIM; ++d) vars_dX[d] = d;
    const DofMap & dX_dof_map = dX_system.get_dof_map();

    const libMesh::DofMap& dof_map = FF_system.get_dof_map();
    FEType fe_type = dof_map.variable_type(0);
    UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
    QGauss qrule (dim, FIFTH);
    fe->attach_quadrature_rule(&qrule);

    UniquePtr<FEBase> fe_face (FEBase::build(dim,fe_type));
    QGauss qface(dim-1, FIFTH);
    fe_face->attach_quadrature_rule(&qface);

//    libMesh::NumericVector<double>& FF_soln = *(FF_system.solution);

    const std::vector<Real>& JxW = fe->get_JxW();
//    const std::vector<libMesh::Point>& q_point = fe->get_xyz();
    const std::vector<std::vector<Real> >& phi = fe->get_phi();
    const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

    TensorValue<double> II(1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0);
    TensorValue<double> FF;
//    double FFdet = 0.0;
    std::vector<std::vector<unsigned int> > dof_indices_dX(NDIM);
    std::vector<dof_id_type> dof_indices_dx;
    std::vector<dof_id_type> dof_indices_dy;
    std::vector<dof_id_type> dof_indices_dz;

    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;

    DenseSubMatrix<Number>
    K00(Ke), K01(Ke), K02(Ke),
    K10(Ke), K11(Ke), K12(Ke),
    K20(Ke), K21(Ke), K22(Ke);

    DenseSubVector<Number>
    F00(Fe), F01(Fe), F02(Fe),
    F10(Fe), F11(Fe), F12(Fe),
    F20(Fe), F21(Fe), F22(Fe);

    std::vector<dof_id_type> dof_indices;
    std::vector<dof_id_type> dof_indices_F00;
    std::vector<dof_id_type> dof_indices_F01;
    std::vector<dof_id_type> dof_indices_F02;

    std::vector<dof_id_type> dof_indices_F10;
    std::vector<dof_id_type> dof_indices_F11;
    std::vector<dof_id_type> dof_indices_F12;

    std::vector<dof_id_type> dof_indices_F20;
    std::vector<dof_id_type> dof_indices_F21;
    std::vector<dof_id_type> dof_indices_F22;

    MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el =  mesh.active_local_elements_end();

    for(; el != end_el; ++el)
    {
        Elem* elem = *el;
        dof_map.dof_indices(elem, dof_indices);

        dX_dof_map.dof_indices(elem, dof_indices_dx, 0); // for dx
        dX_dof_map.dof_indices(elem, dof_indices_dy, 1); // for dy
        dX_dof_map.dof_indices(elem, dof_indices_dz, 2); // for dz

        dof_map.dof_indices(elem, dof_indices_F00, 0);
        dof_map.dof_indices(elem, dof_indices_F01, 1);
        dof_map.dof_indices(elem, dof_indices_F02, 2);

        dof_map.dof_indices(elem, dof_indices_F10, 3);
        dof_map.dof_indices(elem, dof_indices_F11, 4);
        dof_map.dof_indices(elem, dof_indices_F12, 5);

        dof_map.dof_indices(elem, dof_indices_F20, 6);
        dof_map.dof_indices(elem, dof_indices_F21, 7);
        dof_map.dof_indices(elem, dof_indices_F22, 8);

        const unsigned int n_dofs = dof_indices.size();
        const unsigned int n_00_dofs = dof_indices_F00.size();

        fe->reinit(elem);

        Ke.resize(n_dofs, n_dofs);
        Fe.resize(n_dofs);

        K00.reposition(0*n_00_dofs, 0*n_00_dofs, n_00_dofs, n_00_dofs);
        K01.reposition(1*n_00_dofs, 1*n_00_dofs, n_00_dofs, n_00_dofs);
        K02.reposition(2*n_00_dofs, 2*n_00_dofs, n_00_dofs, n_00_dofs);

        K10.reposition(3*n_00_dofs, 3*n_00_dofs, n_00_dofs, n_00_dofs);
        K11.reposition(4*n_00_dofs, 4*n_00_dofs, n_00_dofs, n_00_dofs);
        K12.reposition(5*n_00_dofs, 5*n_00_dofs, n_00_dofs, n_00_dofs);

        K20.reposition(6*n_00_dofs, 6*n_00_dofs, n_00_dofs, n_00_dofs);
        K21.reposition(7*n_00_dofs, 7*n_00_dofs, n_00_dofs, n_00_dofs);
        K22.reposition(8*n_00_dofs, 8*n_00_dofs, n_00_dofs, n_00_dofs);

        F00.reposition(0*n_00_dofs, n_00_dofs);
        F01.reposition(1*n_00_dofs, n_00_dofs);
        F02.reposition(2*n_00_dofs, n_00_dofs);

        F10.reposition(3*n_00_dofs, n_00_dofs);
        F11.reposition(4*n_00_dofs, n_00_dofs);
        F12.reposition(5*n_00_dofs, n_00_dofs);

        F20.reposition(6*n_00_dofs, n_00_dofs);
        F21.reposition(7*n_00_dofs, n_00_dofs);
        F22.reposition(8*n_00_dofs, n_00_dofs);

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {

            for (unsigned int i = 0; i<phi.size(); i++)
                for (unsigned int j = 0; j<phi.size(); j++)
                {
                    K00(i,j) += JxW[qp]*phi[i][qp]*phi[j][qp];
                    K01(i,j) += JxW[qp]*phi[i][qp]*phi[j][qp];
                    K02(i,j) += JxW[qp]*phi[i][qp]*phi[j][qp];

                    K10(i,j) += JxW[qp]*phi[i][qp]*phi[j][qp];
                    K11(i,j) += JxW[qp]*phi[i][qp]*phi[j][qp];
                    K12(i,j) += JxW[qp]*phi[i][qp]*phi[j][qp];

                    K20(i,j) += JxW[qp]*phi[i][qp]*phi[j][qp];
                    K21(i,j) += JxW[qp]*phi[i][qp]*phi[j][qp];
                    K22(i,j) += JxW[qp]*phi[i][qp]*phi[j][qp];
                }

            //we now calculate the FF
            Gradient grad_dx, grad_dy, grad_dz;
            grad_dx.zero();
            grad_dy.zero();
            grad_dz.zero();

            for (unsigned int i=0; i<dof_indices_dx.size(); i++)
            {
                double dx = dX_data(dof_indices_dx[i]);
                double dy = dX_data(dof_indices_dy[i]);
                double dz = dX_data(dof_indices_dz[i]);

                grad_dx.add_scaled(dphi[i][qp], dx);
                grad_dy.add_scaled(dphi[i][qp], dy);
                grad_dz.add_scaled(dphi[i][qp], dz);
            }
            //calculating FF
            FF.zero();
            FF(0,0) = grad_dx(0);
            FF(0,1) = grad_dx(1);
            FF(0,2) = grad_dx(2);

            FF(1,0) = grad_dy(0);
            FF(1,1) = grad_dy(1);
            FF(1,2) = grad_dy(2);

            FF(2,0) = grad_dz(0);
            FF(2,1) = grad_dz(1);
            FF(2,2) = grad_dz(2);
            FF = FF + II;

            for (std::size_t i = 0; i < phi.size(); i++)
            {
                F00(i) += JxW[qp]*phi[i][qp]*FF(0,0);
                F01(i) += JxW[qp]*phi[i][qp]*FF(0,1);
                F02(i) += JxW[qp]*phi[i][qp]*FF(0,2);

                F10(i) += JxW[qp]*phi[i][qp]*FF(1,0);
                F11(i) += JxW[qp]*phi[i][qp]*FF(1,1);
                F12(i) += JxW[qp]*phi[i][qp]*FF(1,2);

                F20(i) += JxW[qp]*phi[i][qp]*FF(2,0);
                F21(i) += JxW[qp]*phi[i][qp]*FF(2,1);
                F22(i) += JxW[qp]*phi[i][qp]*FF(2,2);

            }

        } //end for qp

        dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
        FF_system.matrix->add_matrix(Ke, dof_indices);
        FF_system.rhs->add_vector(Fe, dof_indices);


    } //end for el

    return;
} //end for the function


void MechanicsModel::update_FF_system(EquationSystems* es)
{
    std::cout<<"updating FF_system\n";
    LinearImplicitSystem& FF_system = es->get_system<LinearImplicitSystem>(MechanicsModel::FF_system_name);
    FF_system.solve();
    return;
}
