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

#include "DarcyMass.h"
#include "AHAdivision.h"
#include "MechanicsModel.h"
#include "arteriesPD.h"

#include <stdio.h>

std::vector<int>  AHAdivision::node_aha_list(16403);
std::vector<int>  AHAdivision::element_aha_list(82290);

double AHAdivision::total_aha_volumes[17];

std::string AHAdivision::mesh_division_systemName = "mesh_division_system";
std::string AHAdivision::mesh_division_variableName = "mesh_division_variable";
Order AHAdivision::mesh_division_variable_order = FIRST;



void AHAdivision::setup_mesh_division(EquationSystems* es)
{
    es->add_system<System> (AHAdivision::mesh_division_systemName);
    System& system_mesh_division = es->get_system<System>(AHAdivision::mesh_division_systemName);
    system_mesh_division.add_variable(mesh_division_variableName, mesh_division_variable_order, LAGRANGE);
    return;
}




void AHAdivision::initialize_mesh_division(EquationSystems* es)
{
    System& system_mesh_division = es->get_system<System>(AHAdivision::mesh_division_systemName);
    unsigned int system_num = system_mesh_division.number();

    libMesh::MeshBase& mesh = es->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    const libMesh::DofMap& dof_map = system_mesh_division.get_dof_map();
    libMesh::NumericVector<double>& mesh_division_soln = *system_mesh_division.solution;

    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();

    double test_0 = 0.0;
    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;
        const unsigned int dof_num = nd->dof_number(system_num, 0,0);
        mesh_division_soln.set(dof_num, test_0);
    }
    mesh_division_soln.close();
    mesh_division_soln.localize(*system_mesh_division.current_local_solution);

    return;
}



void AHAdivision:: update_mesh_division(EquationSystems* es)
{

    System& system_mesh_division = es->get_system<System>(AHAdivision::mesh_division_systemName);
    unsigned int system_num = system_mesh_division.number();

    libMesh::MeshBase& mesh = es->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    libMesh::NumericVector<double>& mesh_division_soln = *(system_mesh_division.solution);
    mesh_division_soln.close();
    mesh_division_soln.localize(*system_mesh_division.current_local_solution);

    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();


    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;

        int current_node=nd->id();


        double x_coord=mesh.point(current_node)(0);
        double y_coord=mesh.point(current_node)(1);
        double z_coord=mesh.point(current_node)(2);



        double angle=0.0;
        double angle_shift=0.0;
        double a1=0.0;
        double a2=0.0;
        double b1=0.0;
        double b2=0.0;
        double c1=0.0;
        double c2=0.0;
        double d1=0.0;
        double d2=0.0;
        double e1=0.0;
        double e2=0.0;
        double f1=0.0;
        double f2=0.0;

const unsigned int dof_num = nd->dof_number(system_num, 0, 0);

         if (z_coord > 4.4544)
         {

angle=(atan2(y_coord,x_coord)*(180/3.14159))+180.0;

//specify at what angle region 1 should start (moving from left to right)
//note that regions are defined in a counter clockwise direction
angle_shift=30.0;

a1=0.0;
a2=60.0;

b1=60.0;
b2=120.0;

c1=120.0;
c2=180.0;

d1=180.0;
d2=240.0;

e1=240.0;
e2=300.0;

f1=300.0;
f2=360.0;


if (angle >= a1 && angle < a2)
{
  node_aha_list[current_node]=1;
  mesh_division_soln.set(dof_num, 1.0);
}
else if (angle >= b1 && angle < b2)
{
  node_aha_list[current_node]=2;
  mesh_division_soln.set(dof_num, 2.0);
}
else if (angle >= c1 && angle < c2)
{
  node_aha_list[current_node]=3;
  mesh_division_soln.set(dof_num, 3.0);
}
else if (angle >= d1 && angle < d2)
{
  node_aha_list[current_node]=4;
  mesh_division_soln.set(dof_num, 4.0);
}
else if (angle >= e1 && angle < e2)
{
  node_aha_list[current_node]=5;
  mesh_division_soln.set(dof_num, 5.0);
}
else if (angle >= f1 && angle<=f2)
{
  node_aha_list[current_node]=6;
  mesh_division_soln.set(dof_num, 6.0);
}

         }
         else if(z_coord <= 4.4544 && z_coord > 2.8274)
         {


           angle=(atan2(y_coord,x_coord)*(180/3.14159))+180.0;


   if (angle < 60.0)
   {
     node_aha_list[current_node]=7;
     mesh_division_soln.set(dof_num, 7.0);
   }
   else if (angle >= 60.0 && angle < 120.0)
   {
     node_aha_list[current_node]=8;
     mesh_division_soln.set(dof_num, 8.0);
   }
   else if (angle >= 120.0 && angle < 180.0)
   {
     node_aha_list[current_node]=9;
     mesh_division_soln.set(dof_num, 9.0);
   }
   else if (angle >= 180.0 && angle < 240.0)
   {
     node_aha_list[current_node]=10;
     mesh_division_soln.set(dof_num, 10.0);
   }
   else if (angle >= 240.0 && angle < 300.0)
   {
     node_aha_list[current_node]=11;
     mesh_division_soln.set(dof_num, 11.0);
   }
   else if (angle >= 300.0)
   {
     node_aha_list[current_node]=12;
     mesh_division_soln.set(dof_num, 12.0);
   }


         }
        else if(z_coord <= 2.8274 && z_coord > 1.2004)
        {


          angle=(atan2(y_coord,x_coord)*(180/3.14159))+180.0;


  if (angle < 90.0)
  {
    node_aha_list[current_node]=13;
    mesh_division_soln.set(dof_num, 13.0);
  }
  else if (angle >= 90.0 && angle < 180.0)
  {
    node_aha_list[current_node]=14;
    mesh_division_soln.set(dof_num, 14.0);
  }
  else if (angle >= 180.0 && angle < 270.0)
  {
    node_aha_list[current_node]=15;
    mesh_division_soln.set(dof_num, 15.0);
  }
  else if (angle >= 270.0)
  {
    node_aha_list[current_node]=16;
    mesh_division_soln.set(dof_num, 16.0);
  }



        }
        else if(z_coord<=1.2004)
        {
          node_aha_list[current_node]=17;
          mesh_division_soln.set(dof_num, 17.0);
        }



    }




    mesh_division_soln.close();
    mesh_division_soln.localize(*system_mesh_division.current_local_solution);


    return;
}







void AHAdivision:: find_aha_volumes(EquationSystems* es)
{

  libMesh::MeshBase& mesh = es->get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  //----------------------------------------------------------------------------------
  LinearImplicitSystem& J_system = es->get_system<LinearImplicitSystem>(MechanicsModel::J_system_name);

  J_system.solution->localize(*J_system.current_local_solution);
  NumericVector<double>& J_data = *(J_system.current_local_solution);
  J_data.close();

  const DofMap & dof_map_J = J_system.get_dof_map();

  //----------------------------------------------------------------------------------
       FEType fe_type = dof_map_J.variable_type(0);
       UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
       QGauss qrule (dim, CONSTANT); // Not sure what the best accuracy is here
       fe->attach_quadrature_rule (&qrule);
       UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));
       QGauss qface(dim-1, CONSTANT); // Not sure what the best accuracy is here
  //----------------------------------------------------------------------------------
       fe_face->attach_quadrature_rule (&qface);

       const std::vector<std::vector<Real>> & phi = fe->get_phi();
  //----------------------------------------------------------------------------------

  vector<double> J_vec;
  J_data.localize(J_vec);

  std::vector<dof_id_type> dof_indices_J;



  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();


  for(; el != end_el; ++el)
  {

    Elem* elem = *el;

    int current_element=elem->id();

    fe->reinit (elem);
    dof_map_J.dof_indices (elem, dof_indices_J);
    const unsigned int n_var_dofs = dof_indices_J.size();

    double J_cur = 0.0;
       for (unsigned int j=0; j<n_var_dofs; j++)
        {
         J_cur += phi[j][0]*J_vec[dof_indices_J[j]];
        }

//-------------------------------------------------------------------------------------
int first_node=elem->node_id(0);
int second_node=elem->node_id(1);
int third_node=elem->node_id(2);
int fourth_node=elem->node_id(3);

double x_first_node=mesh.point(first_node)(0);
double y_first_node=mesh.point(first_node)(1);
double z_first_node=mesh.point(first_node)(2);

double x_second_node=mesh.point(second_node)(0);
double y_second_node=mesh.point(second_node)(1);
double z_second_node=mesh.point(second_node)(2);

double x_third_node=mesh.point(third_node)(0);
double y_third_node=mesh.point(third_node)(1);
double z_third_node=mesh.point(third_node)(2);

double x_fourth_node=mesh.point(fourth_node)(0);
double y_fourth_node=mesh.point(fourth_node)(1);
double z_fourth_node=mesh.point(fourth_node)(2);

int first_aha=AHAdivision::node_aha_list[first_node];
int second_aha=AHAdivision::node_aha_list[second_node];
int third_aha=AHAdivision::node_aha_list[third_node];
int fourth_aha=AHAdivision::node_aha_list[fourth_node];

      double x_centroid=elem->centroid()(0);
      double y_centroid=elem->centroid()(1);
      double z_centroid=elem->centroid()(2);

vector<int> subdomain_holder(18, 0);

subdomain_holder[first_aha]=subdomain_holder[first_aha]+1;
subdomain_holder[second_aha]=subdomain_holder[second_aha]+1;
subdomain_holder[third_aha]=subdomain_holder[third_aha]+1;
subdomain_holder[fourth_aha]=subdomain_holder[fourth_aha]+1;


int maxValue = *max_element(subdomain_holder.begin(),subdomain_holder.end());

//if the maxiumum value is 3 (or 4) then the element will be assigned the subdomain id associated
//with these 3 (or 4) nodes.
if(maxValue>2)
{

int MaxElementIndex = std::max_element(subdomain_holder.begin(),subdomain_holder.end()) - subdomain_holder.begin();
AHAdivision::element_aha_list[current_element]=MaxElementIndex;

}
else if(maxValue==2)
{

int MaxElementIndex = std::max_element(subdomain_holder.begin(),subdomain_holder.end()) - subdomain_holder.begin();

//The index where the current maximum is occuring is set to be 0, this will then allow us to find
//out whether this element is a ..,2,2... or a ...,2,1,1,.... and we can then proceed accordingly
//subdomain_holder[MaxElementIndex]=0;

int maxValueAgain = *max_element(subdomain_holder.begin(),subdomain_holder.end());


if(maxValueAgain==1)
{
//this is the case where the element has 2 nodes in 1 subdomains and each of the
//other 2 nodes are in different subdomains, Therefore assign this element the subdomain
//from the '2' nodes
AHAdivision::element_aha_list[current_element]=MaxElementIndex;
}
else if(maxValueAgain==2)
{
//this is the case where the element has 2 nodes in 2 different subdomains (...,2,2....),
//for this we compute the distance between the element centroid and the 4 respective nodes.
//Whichever of these is closest to the element centroid is the subdomain that will be used

  std::vector<double> distance(4, 0.0);

  distance[0]=sqrt(((x_centroid-x_first_node)*(x_centroid-x_first_node))+((y_centroid-y_first_node)*(y_centroid-y_first_node))+((z_centroid-z_first_node)*(z_centroid-z_first_node)));
  distance[1]=sqrt(((x_centroid-x_second_node)*(x_centroid-x_second_node))+((y_centroid-y_second_node)*(y_centroid-y_second_node))+((z_centroid-z_second_node)*(z_centroid-z_second_node)));
  distance[2]=sqrt(((x_centroid-x_third_node)*(x_centroid-x_third_node))+((y_centroid-y_third_node)*(y_centroid-y_third_node))+((z_centroid-z_third_node)*(z_centroid-z_third_node)));
  distance[3]=sqrt(((x_centroid-x_fourth_node)*(x_centroid-x_fourth_node))+((y_centroid-y_fourth_node)*(y_centroid-y_fourth_node))+((z_centroid-z_fourth_node)*(z_centroid-z_fourth_node)));

  int MinElementIndex = std::min_element(distance.begin(),distance.end()) - distance.begin();

  if(MinElementIndex==0)
  {
  AHAdivision::element_aha_list[current_element]=first_aha;
  }
  else if(MinElementIndex==1)
  {
  AHAdivision::element_aha_list[current_element]=second_aha;
  }
  else if(MinElementIndex==2)
  {
  AHAdivision::element_aha_list[current_element]=third_aha;
  }
  else if(MinElementIndex==3)
  {
  AHAdivision::element_aha_list[current_element]=fourth_aha;
  }
  else
  {
    cout << "something has gone wrong" << endl;
  }

}
else
{
  cout << "something has gone wrong" << endl;
}

}
else if(maxValue==1)
{
//this is the case where all 4 nodes are in different subdomains. In this case compute the distance
//between the element centroid and each of the 4 nodes and assign the subdomain according to which node is
//the closest.
std::vector<double> distances(4, 0.0);

distances[0]=sqrt(((x_centroid-x_first_node)*(x_centroid-x_first_node))+((y_centroid-y_first_node)*(y_centroid-y_first_node))+((z_centroid-z_first_node)*(z_centroid-z_first_node)));
distances[1]=sqrt(((x_centroid-x_second_node)*(x_centroid-x_second_node))+((y_centroid-y_second_node)*(y_centroid-y_second_node))+((z_centroid-z_second_node)*(z_centroid-z_second_node)));
distances[2]=sqrt(((x_centroid-x_third_node)*(x_centroid-x_third_node))+((y_centroid-y_third_node)*(y_centroid-y_third_node))+((z_centroid-z_third_node)*(z_centroid-z_third_node)));
distances[3]=sqrt(((x_centroid-x_fourth_node)*(x_centroid-x_fourth_node))+((y_centroid-y_fourth_node)*(y_centroid-y_fourth_node))+((z_centroid-z_fourth_node)*(z_centroid-z_fourth_node)));

int MinElementIndex = std::min_element(distances.begin(),distances.end()) - distances.begin();

if(MinElementIndex==0)
{
AHAdivision::element_aha_list[current_element]=first_aha;
}
else if(MinElementIndex==1)
{
AHAdivision::element_aha_list[current_element]=second_aha;
}
else if(MinElementIndex==2)
{
AHAdivision::element_aha_list[current_element]=third_aha;
}
else if(MinElementIndex==3)
{
AHAdivision::element_aha_list[current_element]=fourth_aha;
}
else
{
  cout << "something has gone wrong" << endl;
}

}

int element_subdomain=AHAdivision::element_aha_list[current_element];
double element_volume=(J_cur)*(elem->volume());

AHAdivision::total_aha_volumes[element_subdomain-1]=AHAdivision::total_aha_volumes[element_subdomain-1]+element_volume;



  }


  for(int qq=0; qq<17; qq++)
   {
       pout << AHAdivision::total_aha_volumes[qq] << "\n";
   }




  return;





}
