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
#include "ConoraryFlowInterface.h"
#include "MechanicsModel.h"
#include "arteriesPD.h"
#include "BoundaryConditions.h"

#include <stdio.h>

double ConoraryFlowInterface::total_subdomain_volumes[100];
double ConoraryFlowInterface::flow_change[100];
std::vector<int>  ConoraryFlowInterface::node_current_subdomain_lists(16403);
std::vector<int>  ConoraryFlowInterface::node_current_aha_lists(16403);
std::vector<int>  ConoraryFlowInterface::element_current_subdomain_lists(82290);

double ConoraryFlowInterface::total_flow_in;
double ConoraryFlowInterface::total_flow_out;
double ConoraryFlowInterface::total_flow_difference;
double ConoraryFlowInterface::total_flow_difference_averaged_over_time;

double ConoraryFlowInterface::total_flow_change;


double ConoraryFlowInterface::total_flow_difference_averaged_over_time_sub1;
double ConoraryFlowInterface::total_flow_difference_averaged_over_time_sub2;
double ConoraryFlowInterface::total_flow_difference_averaged_over_time_sub3;
double ConoraryFlowInterface::total_flow_difference_averaged_over_time_sub4;
double ConoraryFlowInterface::total_flow_difference_averaged_over_time_sub5;
double ConoraryFlowInterface::total_flow_difference_averaged_over_time_sub6;
double ConoraryFlowInterface::total_flow_difference_averaged_over_time_sub7;
double ConoraryFlowInterface::total_flow_difference_averaged_over_time_sub8;
double ConoraryFlowInterface::total_flow_difference_averaged_over_time_sub9;

std::string ConoraryFlowInterface::Pext_systemName = "Cono_Pext_system";
std::string ConoraryFlowInterface::Pext_variableName = "Cono_Pext_variable";
Order ConoraryFlowInterface::Pext_variable_order = FIRST;

std::string ConoraryFlowInterface::subdomain_test_systemName = "subdomain_test_system";
std::string ConoraryFlowInterface::subdomain_test_variableName = "subdomain_test_variable";
Order ConoraryFlowInterface::subdomain_test_variable_order = FIRST;

std::string ConoraryFlowInterface::visualise_coro_data_systemName = "visualise_coro_data_system";
std::string ConoraryFlowInterface::visualise_coro_data_variableName = "visualise_coro_data_variable";
Order ConoraryFlowInterface::visualise_coro_data_variable_order = FIRST;

std::string ConoraryFlowInterface::visualise_source_data_systemName = "visualise_source_data_system";
std::string ConoraryFlowInterface::visualise_source_data_variableName = "visualise_source_data_variable";
Order ConoraryFlowInterface::visualise_source_data_variable_order = FIRST;

std::string ConoraryFlowInterface::visualise_sink_data_systemName = "visualise_sink_data_system";
std::string ConoraryFlowInterface::visualise_sink_data_variableName = "visualise_sink_data_variable";
Order ConoraryFlowInterface::visualise_sink_data_variable_order = FIRST;


std::string ConoraryFlowInterface::Q_elem_systemName = "Cono_Q_elem_system";
std::string ConoraryFlowInterface::Q_elem_variableName = "Cono_Q_elem_variable";
Order ConoraryFlowInterface::Q_elem_variable_order = CONSTANT;


std::vector<double> ConoraryFlowInterface::totalVolume(63,1.0);

bool ConoraryFlowInterface::send_source_from_coro;
bool ConoraryFlowInterface::visualise_coro_data;
bool ConoraryFlowInterface::output_coro_data_to_screen;
bool ConoraryFlowInterface::output_coro_data_to_file;
bool ConoraryFlowInterface::aha_to_subdomain;
bool ConoraryFlowInterface::exp_pao_to_file;


namespace //private name space
{
static const double kappa_1 = 20000.0;
static const double kappa_2 = 330.0;
}


void ConoraryFlowInterface::setup_Pext(EquationSystems* es)
{
    es->add_system<System> (ConoraryFlowInterface::Pext_systemName);
    System& system_Pext = es->get_system<System>(ConoraryFlowInterface::Pext_systemName);
    system_Pext.add_variable(Pext_variableName, Pext_variable_order, LAGRANGE);
    return;
}



void ConoraryFlowInterface::setup_id_test(EquationSystems* es)
{
    es->add_system<System> (ConoraryFlowInterface::subdomain_test_systemName);
    System& system_test = es->get_system<System>(ConoraryFlowInterface::subdomain_test_systemName);
    system_test.add_variable(subdomain_test_variableName, subdomain_test_variable_order, LAGRANGE);
    return;
}



void ConoraryFlowInterface::setup_visualise_coro_data(EquationSystems* es)
{
    es->add_system<System> (ConoraryFlowInterface::visualise_coro_data_systemName);
    System& system_visualise = es->get_system<System>(ConoraryFlowInterface::visualise_coro_data_systemName);
    system_visualise.add_variable(visualise_coro_data_variableName, visualise_coro_data_variable_order, LAGRANGE);

    es->add_system<System> (ConoraryFlowInterface::visualise_source_data_systemName);
    System& system_source = es->get_system<System>(ConoraryFlowInterface::visualise_source_data_systemName);
    system_source.add_variable(visualise_source_data_variableName, visualise_source_data_variable_order, LAGRANGE);

    es->add_system<System> (ConoraryFlowInterface::visualise_sink_data_systemName);
    System& system_sink = es->get_system<System>(ConoraryFlowInterface::visualise_sink_data_systemName);
    system_sink.add_variable(visualise_sink_data_variableName, visualise_sink_data_variable_order, LAGRANGE);

    return;
}



void ConoraryFlowInterface::setup_Q_elem(EquationSystems* es)
{
  es->add_system<System> (ConoraryFlowInterface::Q_elem_systemName);
  System& system_Q_elem = es->get_system<System>(ConoraryFlowInterface::Q_elem_systemName);
  system_Q_elem.add_variable(Q_elem_variableName, Q_elem_variable_order, MONOMIAL);
}




void ConoraryFlowInterface::initialize_Pext(EquationSystems* es)
{
    System& system_Pext = es->get_system<System>(ConoraryFlowInterface::Pext_systemName);
    unsigned int system_num = system_Pext.number();

    libMesh::MeshBase& mesh = es->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    const libMesh::DofMap& dof_map = system_Pext.get_dof_map();
    libMesh::NumericVector<double>& Pext_soln = *system_Pext.solution;

    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();

    double p_0 = 0.0;
    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;
        const unsigned int dof_num = nd->dof_number(system_num, 0,0);
        Pext_soln.set(dof_num, p_0);
    }
    Pext_soln.close();
    Pext_soln.localize(*system_Pext.current_local_solution);

    return;
}




void ConoraryFlowInterface::initialize_id_test(EquationSystems* es)
{
    System& system_test = es->get_system<System>(ConoraryFlowInterface::subdomain_test_systemName);
    unsigned int system_num = system_test.number();

    libMesh::MeshBase& mesh = es->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    const libMesh::DofMap& dof_map = system_test.get_dof_map();
    libMesh::NumericVector<double>& test_soln = *system_test.solution;

    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();

    double test_0 = 0.0;
    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;
        const unsigned int dof_num = nd->dof_number(system_num, 0,0);
        test_soln.set(dof_num, test_0);
    }
    test_soln.close();
    test_soln.localize(*system_test.current_local_solution);

    return;
}



void ConoraryFlowInterface::initialize_visualise_coro_data(EquationSystems* es)
{
    System& system_test = es->get_system<System>(ConoraryFlowInterface::visualise_coro_data_systemName);
    unsigned int system_num = system_test.number();

    System& system_source = es->get_system<System>(ConoraryFlowInterface::visualise_source_data_systemName);
    unsigned int system_num_1 = system_source.number();

    System& system_sink = es->get_system<System>(ConoraryFlowInterface::visualise_sink_data_systemName);
    unsigned int system_num_2 = system_sink.number();

    libMesh::MeshBase& mesh = es->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    const libMesh::DofMap& dof_map = system_test.get_dof_map();
    libMesh::NumericVector<double>& test_soln = *system_test.solution;

    const libMesh::DofMap& dof_map_1 = system_source.get_dof_map();
    libMesh::NumericVector<double>& source_soln = *system_source.solution;

    const libMesh::DofMap& dof_map_2 = system_sink.get_dof_map();
    libMesh::NumericVector<double>& sink_soln = *system_sink.solution;

    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();

    double test_0 = 0.0;
    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;

        const unsigned int dof_num = nd->dof_number(system_num, 0,0);
        const unsigned int dof_num_1 = nd->dof_number(system_num_1, 0,0);
        const unsigned int dof_num_2 = nd->dof_number(system_num_2, 0,0);

        test_soln.set(dof_num, test_0);
        source_soln.set(dof_num_1, test_0);
        sink_soln.set(dof_num_2, test_0);
    }
    test_soln.close();
    source_soln.close();
    sink_soln.close();

    test_soln.localize(*system_test.current_local_solution);
    source_soln.localize(*system_source.current_local_solution);
    sink_soln.localize(*system_sink.current_local_solution);

    return;
}




void ConoraryFlowInterface::initialize_Q_elem(EquationSystems* es)
{

    System& system_Q_elem = es->get_system<System>(ConoraryFlowInterface::Q_elem_systemName);
    std::cout << ConoraryFlowInterface::Q_elem_systemName<<std::endl;
    unsigned int system_num = system_Q_elem.number();

    libMesh::MeshBase& mesh = es->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    const libMesh::DofMap& dof_map = system_Q_elem.get_dof_map();

    libMesh::UniquePtr<QBase> qrule = QBase::build(QGAUSS, NDIM, CONSTANT);
    libMesh::UniquePtr<FEBase> Q_fe(FEBase::build(dim, dof_map.variable_type(0)));
    Q_fe->attach_quadrature_rule(qrule.get());

    libMesh::NumericVector<double>& Q_elem_soln = *system_Q_elem.solution;

    libMesh::MeshBase::const_element_iterator elem_it = mesh.active_local_elements_begin();
    const libMesh::MeshBase::const_element_iterator elem_end = mesh.active_local_elements_end();

    std::vector<unsigned int> Q_dof_indices;

    double Q0 = 0.0;
    for(; elem_it != elem_end; ++elem_it)
    {
        const libMesh::Elem* elem = *elem_it;
        Q_fe->reinit(elem);

        dof_map.dof_indices(elem, Q_dof_indices, 0);

        // initialzed to be some value
        for (unsigned int i=0; i < Q_dof_indices.size(); ++i)
        {
           Q_elem_soln.set(Q_dof_indices[i], Q0);
        }

    }
    Q_elem_soln.close();
    Q_elem_soln.localize(*system_Q_elem.current_local_solution);

    return;

}








void ConoraryFlowInterface:: update_id_test(EquationSystems* es)
{


    System& system_test = es->get_system<System>(ConoraryFlowInterface::subdomain_test_systemName);
    unsigned int system_num = system_test.number();


    libMesh::MeshBase& mesh = es->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();


    libMesh::NumericVector<double>& test_soln = *(system_test.solution);


    test_soln.close();
    test_soln.localize(*system_test.current_local_solution);


    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    int broadcast_root = 0;


if(my_rank == broadcast_root)
{

    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();



    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;
        int current_node=nd->id();
        const unsigned int dof_num = nd->dof_number(system_num, 0, 0);
        test_soln.set(dof_num, ConoraryFlowInterface::node_current_subdomain_lists[current_node]);

        double x_node=mesh.point(current_node)(0);
        double y_node=mesh.point(current_node)(1);
        double z_node=mesh.point(current_node)(2);

        ofstream outdata_xyz_subdomain;
        outdata_xyz_subdomain.open("/xlwork1/scott/full_model_nov_22/xyz_subdomain.dat", std::ofstream::app); // opens the file
        outdata_xyz_subdomain <<
        ConoraryFlowInterface::node_current_subdomain_lists[current_node] << "\t" <<
        x_node << "\t" <<
        y_node << "\t" <<
        z_node << endl;
        outdata_xyz_subdomain.close();
        
    }

}


    test_soln.close();
    test_soln.localize(*system_test.current_local_solution);


    return;
}





void ConoraryFlowInterface:: update_visualise_coro_data(EquationSystems* es)
{


    System& system_test = es->get_system<System>(ConoraryFlowInterface::visualise_coro_data_systemName);
    unsigned int system_num = system_test.number();

    System& system_source = es->get_system<System>(ConoraryFlowInterface::visualise_source_data_systemName);
    unsigned int system_num_1 = system_source.number();

    System& system_sink = es->get_system<System>(ConoraryFlowInterface::visualise_sink_data_systemName);
    unsigned int system_num_2 = system_sink.number();


    libMesh::MeshBase& mesh = es->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    const libMesh::DofMap& dof_map = system_test.get_dof_map();
    libMesh::NumericVector<double>& test_soln = *(system_test.solution);

    const libMesh::DofMap& dof_map_1 = system_source.get_dof_map();
    libMesh::NumericVector<double>& source_soln = *(system_source.solution);

    const libMesh::DofMap& dof_map_2 = system_sink.get_dof_map();
    libMesh::NumericVector<double>& sink_soln = *(system_sink.solution);


    libMesh::MeshBase::const_node_iterator node_it = mesh.active_nodes_begin();
    const libMesh::MeshBase::const_node_iterator node_end = mesh.active_nodes_end();


    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;

        int current_node=nd->id();
        int subdomain=ConoraryFlowInterface::node_current_subdomain_lists[current_node];

        double subdomain_sink_source=Tube::net_from_coro_averaged[subdomain-1];
        double subdomain_source=Tube::source_from_coro[subdomain-1];
        double subdomain_sink=Tube::sink_from_coro[subdomain-1];

        const unsigned int dof_num = nd->dof_number(system_num, 0, 0);
        test_soln.set(dof_num, subdomain_sink_source);

        const unsigned int dof_num_1 = nd->dof_number(system_num_1, 0, 0);
        source_soln.set(dof_num_1, subdomain_source);

        const unsigned int dof_num_2 = nd->dof_number(system_num_2, 0, 0);
        sink_soln.set(dof_num_2, subdomain_sink);


    }
    test_soln.close();
    test_soln.localize(*system_test.current_local_solution);

    source_soln.close();
    source_soln.localize(*system_source.current_local_solution);

    sink_soln.close();
    sink_soln.localize(*system_sink.current_local_solution);

    return;
}




//the following takes the sink/source data which has been generated and saved into .txt files
//by the coronary network code and distributes it into the different subdomains around the
//ventricle

void ConoraryFlowInterface:: update_Q_elem(EquationSystems* es,const int iter_num)
{

    System& system_Q_elem = es->get_system<System>(Q_elem_systemName);
    unsigned int system_num = system_Q_elem.number();

    libMesh::MeshBase& mesh = es->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    const libMesh::DofMap& dof_map = system_Q_elem.get_dof_map();
    libMesh::NumericVector<double>& Q_soln_elem = *system_Q_elem.solution;

     std::vector<dof_id_type> dof_indices;

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

     const std::vector<Real> & JxW = fe->get_JxW();
     const std::vector<std::vector<Real>> & phi = fe->get_phi();
//----------------------------------------------------------------------------------

     vector<double> J_vec;
     J_data.localize(J_vec);

     std::vector<dof_id_type> dof_indices_J;



     for (int qq = 0; qq < 9; qq++)
     {
      ConoraryFlowInterface::total_subdomain_volumes[qq]=0.0;
     }


if(iter_num==1)
{
//we only need to assign a subdomain number to each element at the first time step,
//otherwise it is just the volume of the subdomains that needs to be computed


int my_rank;
MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
int broadcast_root = 0;

if(my_rank == broadcast_root)
{


MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

for(; el != end_el; ++el)

{

  Elem* elem = *el;
  dof_map.dof_indices(elem, dof_indices, 0);

  //get a reference to the element we are currently on
  int current_element=elem->id();

  //-------------------------------------------------------------------------
  //find out which 4 nodes are positioned on the element, subsequently
  //find their (x,y,z) co-ordinates which will be needed for later
  //calculations and also use ConoraryFlowInterface::node_current_subdomain_lists to find
  //the subdomain which has been assigned to each of these nodes.
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


  int first_subdomain=ConoraryFlowInterface::node_current_subdomain_lists[first_node];
  int second_subdomain=ConoraryFlowInterface::node_current_subdomain_lists[second_node];
  int third_subdomain=ConoraryFlowInterface::node_current_subdomain_lists[third_node];
  int fourth_subdomain=ConoraryFlowInterface::node_current_subdomain_lists[fourth_node];
  //-------------------------------------------------------------------------

  //find the co-ordinates of the element centroid for the current element
  double x_centroid=elem->centroid()(0);
  double y_centroid=elem->centroid()(1);
  double z_centroid=elem->centroid()(2);

  //set up a zero vector which will add 1 to the subdomain index, note this for ease this is set to have length 10
  //but the subdomains only run from 1 to 9 so that the 0th entry of the vector can always be
  //missed out. As an example of how this works (0,2,1,0,1,0,0,0,0,0) would be
  //an element that has 2 nodes in subdomain 1, 1 node in subdomain 2 and 1 node
  //in subdomain 4.

  vector<int> subdomain_holder(10, 0);

  //set up this vector as described above
  subdomain_holder[first_subdomain]=subdomain_holder[first_subdomain]+1;
  subdomain_holder[second_subdomain]=subdomain_holder[second_subdomain]+1;
  subdomain_holder[third_subdomain]=subdomain_holder[third_subdomain]+1;
  subdomain_holder[fourth_subdomain]=subdomain_holder[fourth_subdomain]+1;

  //find the maxiumum value from subdomain holder. If the maximum value is 3 or 4
  //then this means that at least 3 nodes lie in the same subdomain, if this is the case
  //then this element will be assigned that subdomain_id. If the maxiumum value is 2 then the
  //options are either (..,2,1,1..) or (...,2,2,....). In the case of a (2,1,1) this element will be assigned
  //the subdomain id which 2 nodes have. In the case of a (2,2) the distance is computed between each
  //of the 4 nodes and the element centroid, the node to which the centroid is nearest is the subdomain
  //id that will be used. Finally if the maximum value is 1 that means that all 4 nodes lie in different subdomains.
  //in this case, the distance is computed between the node and the element centroid and whichever node is nearest
  //is the one that is used for the subdomain id.

  int maxValue = *max_element(subdomain_holder.begin(),subdomain_holder.end());

  //if the maxiumum value is 3 (or 4) then the element will be assigned the subdomain id associated
  //with these 3 (or 4) nodes.
  if(maxValue>2)
  {

  int MaxElementIndex = std::max_element(subdomain_holder.begin(),subdomain_holder.end()) - subdomain_holder.begin();
  ConoraryFlowInterface::element_current_subdomain_lists[current_element]=MaxElementIndex;

  }

  else if(maxValue==2)
  {

  int MaxElementIndex = std::max_element(subdomain_holder.begin(),subdomain_holder.end()) - subdomain_holder.begin();

  //The index where the current maximum is occuring is set to be 0, this will then allow us to find
  //out whether this element is a ..,2,2... or a ...,2,1,1,.... and we can then proceed accordingly

  subdomain_holder[MaxElementIndex]=0;
  int maxValueAgain = *max_element(subdomain_holder.begin(),subdomain_holder.end());



  if(maxValueAgain==1)
  {
  //this is the case where the element has 2 nodes in 1 subdomains and each of the
  //other 2 nodes are in different subdomains, Therefore assign this element the subdomain
  //from the '2' nodes
  ConoraryFlowInterface::element_current_subdomain_lists[current_element]=MaxElementIndex;
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
    ConoraryFlowInterface::element_current_subdomain_lists[current_element]=first_subdomain;
    }
    else if(MinElementIndex==1)
    {
    ConoraryFlowInterface::element_current_subdomain_lists[current_element]=second_subdomain;
    }
    else if(MinElementIndex==2)
    {
    ConoraryFlowInterface::element_current_subdomain_lists[current_element]=third_subdomain;
    }
    else if(MinElementIndex==3)
    {
    ConoraryFlowInterface::element_current_subdomain_lists[current_element]=fourth_subdomain;
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
  ConoraryFlowInterface::element_current_subdomain_lists[current_element]=first_subdomain;
  }
  else if(MinElementIndex==1)
  {
  ConoraryFlowInterface::element_current_subdomain_lists[current_element]=second_subdomain;
  }
  else if(MinElementIndex==2)
  {
  ConoraryFlowInterface::element_current_subdomain_lists[current_element]=third_subdomain;
  }
  else if(MinElementIndex==3)
  {
  ConoraryFlowInterface::element_current_subdomain_lists[current_element]=fourth_subdomain;
  }
  else
  {
    cout << "something has gone wrong" << endl;
  }

  }

  int element_subdomain=ConoraryFlowInterface::element_current_subdomain_lists[current_element];
  double element_volume=elem->volume();

  ConoraryFlowInterface::total_subdomain_volumes[element_subdomain-1]=ConoraryFlowInterface::total_subdomain_volumes[element_subdomain-1]+element_volume;
  Tube::number_of_elements_in_subdomain[element_subdomain-1]=Tube::number_of_nodes_in_subdomain[element_subdomain-1]+1;


}


}

for (int qq = 0; qq < 9; qq++)
{
MPI_Bcast(&ConoraryFlowInterface::total_subdomain_volumes[qq], 1, MPI_DOUBLE, broadcast_root, MPI_COMM_WORLD);
}

for (int qq = 0; qq < BoundaryConditions::total_elem_num; qq++)
{
MPI_Bcast(&ConoraryFlowInterface::element_current_subdomain_lists[qq], 1, MPI_INT, broadcast_root, MPI_COMM_WORLD);
}



}
else
{

  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  int broadcast_root = 0;

  if(my_rank == broadcast_root)
  {
  for (;  el != end_el; ++el)
  {
    Elem* elem = *el;
    dof_map.dof_indices(elem, dof_indices, 0);


    fe->reinit (elem);
    dof_map_J.dof_indices (elem, dof_indices_J);
    const unsigned int n_var_dofs = dof_indices_J.size();

    double J_cur = 0.0;
       for (unsigned int j=0; j<n_var_dofs; j++)
        {
         J_cur += phi[j][0]*J_vec[dof_indices_J[j]];
        }


    double element_volume=(J_cur)*(elem->volume());
    ConoraryFlowInterface::total_subdomain_volumes[ConoraryFlowInterface::element_current_subdomain_lists[elem->id()]-1]=ConoraryFlowInterface::total_subdomain_volumes[ConoraryFlowInterface::element_current_subdomain_lists[elem->id()]-1]+element_volume;
  }
  }

  for (int qq = 0; qq < 9; qq++)
  {
  MPI_Bcast(&ConoraryFlowInterface::total_subdomain_volumes[qq], 1, MPI_DOUBLE, broadcast_root, MPI_COMM_WORLD);
  }

}

if (ConoraryFlowInterface::output_coro_data_to_screen)
  {
    pout << "sub 1 volume   " << ConoraryFlowInterface::total_subdomain_volumes[0] << endl;
    pout << "sub 2 volume   " << ConoraryFlowInterface::total_subdomain_volumes[1] << endl;
    pout << "sub 3 volume   " << ConoraryFlowInterface::total_subdomain_volumes[2] << endl;
    pout << "sub 4 volume   " << ConoraryFlowInterface::total_subdomain_volumes[3] << endl;
    pout << "sub 5 volume   " << ConoraryFlowInterface::total_subdomain_volumes[4] << endl;
    pout << "sub 6 volume   " << ConoraryFlowInterface::total_subdomain_volumes[5] << endl;
    pout << "sub 7 volume   " << ConoraryFlowInterface::total_subdomain_volumes[6] << endl;
    pout << "sub 8 volume   " << ConoraryFlowInterface::total_subdomain_volumes[7] << endl;
    pout << "sub 9 volume   " << ConoraryFlowInterface::total_subdomain_volumes[8] << endl;
 }

    return;
}



void ConoraryFlowInterface:: update_Pext(EquationSystems* es,const double time)
{

    const MeshBase & mesh = es->get_mesh();
    const int dim = mesh.mesh_dimension();

//This calls the pressure from the ib code, this is a first order variable and
//is saved on the nodes, needs to be shifted from the nodes to the element centroid

    System & system_lambda = es->get_system<System> ("p_f interpolation system");
    system_lambda.solution->localize(*system_lambda.current_local_solution);
    NumericVector<double>& lambda_data = *(system_lambda.current_local_solution);
    lambda_data.close();
    const DofMap & dof_map_lambda = system_lambda.get_dof_map();


//----------------------------------------------------------------------------------
    LinearImplicitSystem& J_system = es->get_system<LinearImplicitSystem>(MechanicsModel::J_system_name);

    J_system.solution->localize(*J_system.current_local_solution);
    NumericVector<double>& J_data = *(J_system.current_local_solution);
    J_data.close();

    const DofMap & dof_map_J = J_system.get_dof_map();

//***************************************************************************
FEType fe_type = dof_map_lambda.variable_type(0);
UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
QGauss qrule (dim, CONSTANT);
fe->attach_quadrature_rule (&qrule);
UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));
QGauss qface(dim-1, CONSTANT);
//***************************************************************************

fe_face->attach_quadrature_rule (&qface);
const std::vector<std::vector<Real>> & phi = fe->get_phi();


//This calls the sigma dilatory system which is computed in the mechanics model class,
//note that this is reconstructed as a zero order variable and will be saved on the element centroid
    System & system_dilatory = es->get_system<System> ("sigma_dil reconstruction system");
    system_dilatory.solution->localize(*system_dilatory.current_local_solution);
    NumericVector<double>& dilatory_data = *(system_dilatory.current_local_solution);
    dilatory_data.close();

//This calls the sigma deviatoric system which is computed in the mechanics model class,
//note that this is reconstructed as a zero order variable and will be saved on the element centroid
    System & system_deviatoric = es->get_system<System> ("sigma_dev reconstruction system");
    system_deviatoric.solution->localize(*system_deviatoric.current_local_solution);
    NumericVector<double>& deviatoric_data = *(system_deviatoric.current_local_solution);
    deviatoric_data.close();

    vector<double> lambda_vec;
    lambda_data.localize(lambda_vec);

    vector<double> dilatory_vec;
    dilatory_data.localize(dilatory_vec);

    vector<double> deviatoric_vec;
    deviatoric_data.localize(deviatoric_vec);

    int dilatory_system_num = system_dilatory.number();
    int deviatoric_system_num = system_deviatoric.number();

    std::vector<dof_id_type> dof_indices_lambda;


    vector<double> J_vec;
    J_data.localize(J_vec);
    std::vector<dof_id_type> dof_indices_J;


    int my_rank3;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank3);
    int broadcast_root3 = 0;

    if(my_rank3 == broadcast_root3)
    {

//----------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------

MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
const MeshBase::const_element_iterator end_el = mesh.active_elements_end();
  for(; el != end_el; ++el)
{
        Elem* elem = *el;

        fe->reinit (elem);

dof_map_lambda.dof_indices (elem, dof_indices_lambda);
const unsigned int n_var_dofs = dof_indices_lambda.size();



dof_map_J.dof_indices (elem, dof_indices_J);
const unsigned int n_var_dofs_J = dof_indices_J.size();

double J_cur = 0.0;
   for (unsigned int j=0; j<n_var_dofs_J; j++)
    {
     J_cur += phi[j][0]*J_vec[dof_indices_J[j]];
    }


double element_volume=(J_cur)*(elem->volume());

double lambda_centroid=0.0;
for (unsigned int j=0; j<n_var_dofs; j++)
{
lambda_centroid += phi[j][0]*lambda_vec[dof_indices_lambda[j]];
}

const int dof_dilatory_index_00 = elem->dof_number(dilatory_system_num,0,0);
const int dof_dilatory_index_11 = elem->dof_number(dilatory_system_num,4,0);
const int dof_dilatory_index_22 = elem->dof_number(dilatory_system_num,8,0);

const int dof_deviatoric_index_00 = elem->dof_number(deviatoric_system_num,0,0);
const int dof_deviatoric_index_11 = elem->dof_number(deviatoric_system_num,4,0);
const int dof_deviatoric_index_22 = elem->dof_number(deviatoric_system_num,8,0);


double dilatory_00 = dilatory_vec[dof_dilatory_index_00];
double dilatory_11 = dilatory_vec[dof_dilatory_index_11];
double dilatory_22 = dilatory_vec[dof_dilatory_index_22];

double deviatoric_00 = deviatoric_vec[dof_deviatoric_index_00];
double deviatoric_11 = deviatoric_vec[dof_deviatoric_index_11];
double deviatoric_22 = deviatoric_vec[dof_deviatoric_index_22];

double rhs_pext=0.0;

rhs_pext=((1.0/3.0)*(dilatory_00+dilatory_11+dilatory_22))+((1.0/3.0)*(deviatoric_00+deviatoric_11+deviatoric_22))-lambda_centroid;
rhs_pext=rhs_pext*element_volume;


if(ConoraryFlowInterface::element_current_subdomain_lists[elem->id()] == 1)
{
Tube::pressure_to_coronary_current[0]=Tube::pressure_to_coronary_current[0]+rhs_pext;
}
else if(ConoraryFlowInterface::element_current_subdomain_lists[elem->id()] == 2)
{
Tube::pressure_to_coronary_current[1]=Tube::pressure_to_coronary_current[1]+rhs_pext;
}
else if(ConoraryFlowInterface::element_current_subdomain_lists[elem->id()] == 3)
{
Tube::pressure_to_coronary_current[2]=Tube::pressure_to_coronary_current[2]+rhs_pext;
}
else if(ConoraryFlowInterface::element_current_subdomain_lists[elem->id()] == 4)
{
Tube::pressure_to_coronary_current[3]=Tube::pressure_to_coronary_current[3]+rhs_pext;
}
else if(ConoraryFlowInterface::element_current_subdomain_lists[elem->id()] == 5)
{
Tube::pressure_to_coronary_current[4]=Tube::pressure_to_coronary_current[4]+rhs_pext;
}
else if(ConoraryFlowInterface::element_current_subdomain_lists[elem->id()] == 6)
{
Tube::pressure_to_coronary_current[5]=Tube::pressure_to_coronary_current[5]+rhs_pext;
}
else if(ConoraryFlowInterface::element_current_subdomain_lists[elem->id()] == 7)
{
Tube::pressure_to_coronary_current[6]=Tube::pressure_to_coronary_current[6]+rhs_pext;
}
else if(ConoraryFlowInterface::element_current_subdomain_lists[elem->id()] == 8)
{
Tube::pressure_to_coronary_current[7]=Tube::pressure_to_coronary_current[7]+rhs_pext;
}
else if(ConoraryFlowInterface::element_current_subdomain_lists[elem->id()] == 9)
{
Tube::pressure_to_coronary_current[8]=Tube::pressure_to_coronary_current[8]+rhs_pext;
}
else
{
  pout << "something has gone wrong, a node has not been assigned a subdomain in Coronary flow Interface" << "\n";
}


    }

  }



    return;
}



void ConoraryFlowInterface::aha_to_subdomains(EquationSystems* es)
{

    libMesh::MeshBase& mesh = es->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    libMesh::MeshBase::const_node_iterator node_it = mesh.local_nodes_begin();
    libMesh::MeshBase::const_node_iterator node_end = mesh.local_nodes_end();


    for(; node_it != node_end; ++node_it)
    {
        const libMesh::Node* nd = *node_it;
        int current_node=nd->id();
        int current_aha_node=ConoraryFlowInterface::node_current_aha_lists[current_node];

        if(current_aha_node == 1)
        {
         ConoraryFlowInterface::node_current_subdomain_lists[current_node]=1;
        }
        else if(current_aha_node == 2)
        {
         ConoraryFlowInterface::node_current_subdomain_lists[current_node]=1;
        }
        else if(current_aha_node == 3)
        {
         ConoraryFlowInterface::node_current_subdomain_lists[current_node]=9;
        }
        else if(current_aha_node == 4)
        {
         ConoraryFlowInterface::node_current_subdomain_lists[current_node]=2;
        }
        else if(current_aha_node == 5)
        {
         ConoraryFlowInterface::node_current_subdomain_lists[current_node]=3;
        }
        else if(current_aha_node == 6)
        {
         ConoraryFlowInterface::node_current_subdomain_lists[current_node]=3;
        }
        else if(current_aha_node == 7)
        {
         ConoraryFlowInterface::node_current_subdomain_lists[current_node]=4;
        }
        else if(current_aha_node == 8)
        {
         ConoraryFlowInterface::node_current_subdomain_lists[current_node]=4;
        }
        else if(current_aha_node == 9)
        {
         ConoraryFlowInterface::node_current_subdomain_lists[current_node]=9;
        }
        else if(current_aha_node == 10)
        {
         ConoraryFlowInterface::node_current_subdomain_lists[current_node]=2;
        }
        else if(current_aha_node == 11)
        {
         ConoraryFlowInterface::node_current_subdomain_lists[current_node]=6;
        }
        else if(current_aha_node == 12)
        {
         ConoraryFlowInterface::node_current_subdomain_lists[current_node]=5;
        }
        else if(current_aha_node == 13)
        {
         ConoraryFlowInterface::node_current_subdomain_lists[current_node]=7;
        }
        else if(current_aha_node == 14)
        {
         ConoraryFlowInterface::node_current_subdomain_lists[current_node]=7;
        }
        else if(current_aha_node == 15)
        {
         ConoraryFlowInterface::node_current_subdomain_lists[current_node]=2;
        }
        else if(current_aha_node == 16)
        {
         ConoraryFlowInterface::node_current_subdomain_lists[current_node]=8;
        }
        else if(current_aha_node == 17)
        {
         ConoraryFlowInterface::node_current_subdomain_lists[current_node]=7;
        }


    }










    return;
}







//the following takes the sink/source data which has been generated and saved into .txt files
//by the coronary network code and distributes it into the different subdomains around the
//ventricle

void ConoraryFlowInterface:: update_Q_elem_hex(EquationSystems* es,const int iter_num)
{
    std::cout<< "updating Q_elem_system\n";
    System& system_Q_elem = es->get_system<System>(Q_elem_systemName);
    unsigned int system_num = system_Q_elem.number();

    libMesh::MeshBase& mesh = es->get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    const libMesh::DofMap& dof_map = system_Q_elem.get_dof_map();
    libMesh::NumericVector<double>& Q_soln_elem = *system_Q_elem.solution;

     std::vector<dof_id_type> dof_indices;

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

     const std::vector<Real> & JxW = fe->get_JxW();
     const std::vector<std::vector<Real>> & phi = fe->get_phi();
//----------------------------------------------------------------------------------

     vector<double> J_vec;
     J_data.localize(J_vec);

     std::vector<dof_id_type> dof_indices_J;


if(iter_num==1)
{
//we only need to assign a subdomain number to each element at the first time step,
//otherwise it is just the volume of the subdomains that needs to be computed


int my_rank;
MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
int broadcast_root = 0;

if(my_rank == broadcast_root)
{


MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

for(; el != end_el; ++el)

{

  Elem* elem = *el;
  dof_map.dof_indices(elem, dof_indices, 0);

  //get a reference to the element we are currently on
  int current_element=elem->id();

  //-------------------------------------------------------------------------
  //find out which 8 nodes are positioned on the element, subsequently
  //find their (x,y,z) co-ordinates which will be needed for later
  //calculations and also use ConoraryFlowInterface::node_current_subdomain_lists to find
  //the subdomain which has been assigned to each of these nodes.
  int first_node=elem->node_id(0);
  int second_node=elem->node_id(1);
  int third_node=elem->node_id(2);
  int fourth_node=elem->node_id(3);
  int fifth_node=elem->node_id(4);
  int sixth_node=elem->node_id(5);
  int seventh_node=elem->node_id(6);
  int eighth_node=elem->node_id(7);

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

  double x_fifth_node=mesh.point(fifth_node)(0);
  double y_fifth_node=mesh.point(fifth_node)(1);
  double z_fifth_node=mesh.point(fifth_node)(2);

  double x_sixth_node=mesh.point(sixth_node)(0);
  double y_sixth_node=mesh.point(sixth_node)(1);
  double z_sixth_node=mesh.point(sixth_node)(2);

  double x_seventh_node=mesh.point(seventh_node)(0);
  double y_seventh_node=mesh.point(seventh_node)(1);
  double z_seventh_node=mesh.point(seventh_node)(2);

  double x_eighth_node=mesh.point(eighth_node)(0);
  double y_eighth_node=mesh.point(eighth_node)(1);
  double z_eighth_node=mesh.point(eighth_node)(2);


  int first_subdomain=ConoraryFlowInterface::node_current_subdomain_lists[first_node];
  int second_subdomain=ConoraryFlowInterface::node_current_subdomain_lists[second_node];
  int third_subdomain=ConoraryFlowInterface::node_current_subdomain_lists[third_node];
  int fourth_subdomain=ConoraryFlowInterface::node_current_subdomain_lists[fourth_node];
  int fifth_subdomain=ConoraryFlowInterface::node_current_subdomain_lists[fifth_node];
  int sixth_subdomain=ConoraryFlowInterface::node_current_subdomain_lists[sixth_node];
  int seventh_subdomain=ConoraryFlowInterface::node_current_subdomain_lists[seventh_node];
  int eighth_subdomain=ConoraryFlowInterface::node_current_subdomain_lists[eighth_node];
  //-------------------------------------------------------------------------

  //find the co-ordinates of the element centroid for the current element
  double x_centroid=elem->centroid()(0);
  double y_centroid=elem->centroid()(1);
  double z_centroid=elem->centroid()(2);

  //set up a zero vector which will add 1 to the subdomain index, note this for ease this is set to have length 10
  //but the subdomains only run from 1 to 9 so that the 0th entry of the vector can always be
  //missed out. As an example of how this works (0,2,1,0,1,0,0,0,0,0) would be
  //an element that has 2 nodes in subdomain 1, 1 node in subdomain 2 and 1 node
  //in subdomain 4.

  vector<int> subdomain_holder(10, 0);

  //set up this vector as described above
  subdomain_holder[first_subdomain]=subdomain_holder[first_subdomain]+1;
  subdomain_holder[second_subdomain]=subdomain_holder[second_subdomain]+1;
  subdomain_holder[third_subdomain]=subdomain_holder[third_subdomain]+1;
  subdomain_holder[fourth_subdomain]=subdomain_holder[fourth_subdomain]+1;
  subdomain_holder[fifth_subdomain]=subdomain_holder[fifth_subdomain]+1;
  subdomain_holder[sixth_subdomain]=subdomain_holder[sixth_subdomain]+1;
  subdomain_holder[seventh_subdomain]=subdomain_holder[seventh_subdomain]+1;
  subdomain_holder[eighth_subdomain]=subdomain_holder[eighth_subdomain]+1;

  //find the maxiumum value from subdomain holder. If the maximum value is 3 or 4
  //then this means that at least 3 nodes lie in the same subdomain, if this is the case
  //then this element will be assigned that subdomain_id. If the maxiumum value is 2 then the
  //options are either (..,2,1,1..) or (...,2,2,....). In the case of a (2,1,1) this element will be assigned
  //the subdomain id which 2 nodes have. In the case of a (2,2) the distance is computed between each
  //of the 4 nodes and the element centroid, the node to which the centroid is nearest is the subdomain
  //id that will be used. Finally if the maximum value is 1 that means that all 4 nodes lie in different subdomains.
  //in this case, the distance is computed between the node and the element centroid and whichever node is nearest
  //is the one that is used for the subdomain id.

  int maxValue = *max_element(subdomain_holder.begin(),subdomain_holder.end());


  if(maxValue>4)
  {

  int MaxElementIndex = std::max_element(subdomain_holder.begin(),subdomain_holder.end()) - subdomain_holder.begin();
  ConoraryFlowInterface::element_current_subdomain_lists[current_element]=MaxElementIndex;

  }
  else if(maxValue==4)
  {

//this tells you the position of the current maximum value of 4.
int CurrentMaxElementIndex = std::max_element(subdomain_holder.begin(),subdomain_holder.end()) - subdomain_holder.begin();

subdomain_holder[CurrentMaxElementIndex]=0;

int maxValueAgain = *max_element(subdomain_holder.begin(),subdomain_holder.end());
int SecondMaxElementIndex = std::max_element(subdomain_holder.begin(),subdomain_holder.end()) - subdomain_holder.begin();

if(maxValueAgain == 4)
{
std::vector<double> distances(2, 0.0);


if(CurrentMaxElementIndex==1)
{
  distances[0]=sqrt(((x_centroid-x_first_node)*(x_centroid-x_first_node))+((y_centroid-y_first_node)*(y_centroid-y_first_node))+((z_centroid-z_first_node)*(z_centroid-z_first_node)));
}
else if(CurrentMaxElementIndex==2)
{
  distances[0]=sqrt(((x_centroid-x_second_node)*(x_centroid-x_second_node))+((y_centroid-y_second_node)*(y_centroid-y_second_node))+((z_centroid-z_second_node)*(z_centroid-z_second_node)));
}
else if(CurrentMaxElementIndex==3)
{
  distances[0]=sqrt(((x_centroid-x_third_node)*(x_centroid-x_third_node))+((y_centroid-y_third_node)*(y_centroid-y_third_node))+((z_centroid-z_third_node)*(z_centroid-z_third_node)));
}
else if(CurrentMaxElementIndex==4)
{
  distances[0]=sqrt(((x_centroid-x_fourth_node)*(x_centroid-x_fourth_node))+((y_centroid-y_fourth_node)*(y_centroid-y_fourth_node))+((z_centroid-z_fourth_node)*(z_centroid-z_fourth_node)));
}
else if(CurrentMaxElementIndex==5)
{
  distances[0]=sqrt(((x_centroid-x_fifth_node)*(x_centroid-x_fifth_node))+((y_centroid-y_fifth_node)*(y_centroid-y_fifth_node))+((z_centroid-z_fifth_node)*(z_centroid-z_fifth_node)));
}
else if(CurrentMaxElementIndex==6)
{
  distances[0]=sqrt(((x_centroid-x_sixth_node)*(x_centroid-x_sixth_node))+((y_centroid-y_sixth_node)*(y_centroid-y_sixth_node))+((z_centroid-z_sixth_node)*(z_centroid-z_sixth_node)));
}
else if(CurrentMaxElementIndex==7)
{
  distances[0]=sqrt(((x_centroid-x_seventh_node)*(x_centroid-x_seventh_node))+((y_centroid-y_seventh_node)*(y_centroid-y_seventh_node))+((z_centroid-z_seventh_node)*(z_centroid-z_seventh_node)));
}
else if(CurrentMaxElementIndex==8)
{
  distances[0]=sqrt(((x_centroid-x_eighth_node)*(x_centroid-x_eighth_node))+((y_centroid-y_eighth_node)*(y_centroid-y_eighth_node))+((z_centroid-z_eighth_node)*(z_centroid-z_eighth_node)));
}




if(SecondMaxElementIndex==1)
{
  distances[1]=sqrt(((x_centroid-x_first_node)*(x_centroid-x_first_node))+((y_centroid-y_first_node)*(y_centroid-y_first_node))+((z_centroid-z_first_node)*(z_centroid-z_first_node)));
}
else if(SecondMaxElementIndex==2)
{
  distances[1]=sqrt(((x_centroid-x_second_node)*(x_centroid-x_second_node))+((y_centroid-y_second_node)*(y_centroid-y_second_node))+((z_centroid-z_second_node)*(z_centroid-z_second_node)));
}
else if(SecondMaxElementIndex==3)
{
  distances[1]=sqrt(((x_centroid-x_third_node)*(x_centroid-x_third_node))+((y_centroid-y_third_node)*(y_centroid-y_third_node))+((z_centroid-z_third_node)*(z_centroid-z_third_node)));
}
else if(SecondMaxElementIndex==4)
{
  distances[1]=sqrt(((x_centroid-x_fourth_node)*(x_centroid-x_fourth_node))+((y_centroid-y_fourth_node)*(y_centroid-y_fourth_node))+((z_centroid-z_fourth_node)*(z_centroid-z_fourth_node)));
}
else if(SecondMaxElementIndex==5)
{
  distances[1]=sqrt(((x_centroid-x_fifth_node)*(x_centroid-x_fifth_node))+((y_centroid-y_fifth_node)*(y_centroid-y_fifth_node))+((z_centroid-z_fifth_node)*(z_centroid-z_fifth_node)));
}
else if(SecondMaxElementIndex==6)
{
  distances[1]=sqrt(((x_centroid-x_sixth_node)*(x_centroid-x_sixth_node))+((y_centroid-y_sixth_node)*(y_centroid-y_sixth_node))+((z_centroid-z_sixth_node)*(z_centroid-z_sixth_node)));
}
else if(SecondMaxElementIndex==7)
{
  distances[1]=sqrt(((x_centroid-x_seventh_node)*(x_centroid-x_seventh_node))+((y_centroid-y_seventh_node)*(y_centroid-y_seventh_node))+((z_centroid-z_seventh_node)*(z_centroid-z_seventh_node)));
}
else if(SecondMaxElementIndex==8)
{
  distances[1]=sqrt(((x_centroid-x_eighth_node)*(x_centroid-x_eighth_node))+((y_centroid-y_eighth_node)*(y_centroid-y_eighth_node))+((z_centroid-z_eighth_node)*(z_centroid-z_eighth_node)));
}


int MinElementIndexDistances = std::min_element(distances.begin(),distances.end()) - distances.begin();
if(MinElementIndexDistances == 0)
{
ConoraryFlowInterface::element_current_subdomain_lists[current_element]=CurrentMaxElementIndex;
}
else if(MinElementIndexDistances == 1)
{
ConoraryFlowInterface::element_current_subdomain_lists[current_element]=SecondMaxElementIndex;
}

}
else if (maxValueAgain < 4)
{
  ConoraryFlowInterface::element_current_subdomain_lists[current_element]=CurrentMaxElementIndex;
}


  }



  else if(maxValue==3)
  {

    int CurrentMaxElementIndex = std::max_element(subdomain_holder.begin(),subdomain_holder.end()) - subdomain_holder.begin();

    subdomain_holder[CurrentMaxElementIndex]=0;

    int maxValueAgain = *max_element(subdomain_holder.begin(),subdomain_holder.end());
    int SecondMaxElementIndex = std::max_element(subdomain_holder.begin(),subdomain_holder.end()) - subdomain_holder.begin();

    if(maxValueAgain == 3)
    {

      std::vector<double> distances(2, 0.0);

      if(CurrentMaxElementIndex==1)
      {
        distances[0]=sqrt(((x_centroid-x_first_node)*(x_centroid-x_first_node))+((y_centroid-y_first_node)*(y_centroid-y_first_node))+((z_centroid-z_first_node)*(z_centroid-z_first_node)));
      }
      else if(CurrentMaxElementIndex==2)
      {
        distances[0]=sqrt(((x_centroid-x_second_node)*(x_centroid-x_second_node))+((y_centroid-y_second_node)*(y_centroid-y_second_node))+((z_centroid-z_second_node)*(z_centroid-z_second_node)));
      }
      else if(CurrentMaxElementIndex==3)
      {
        distances[0]=sqrt(((x_centroid-x_third_node)*(x_centroid-x_third_node))+((y_centroid-y_third_node)*(y_centroid-y_third_node))+((z_centroid-z_third_node)*(z_centroid-z_third_node)));
      }
      else if(CurrentMaxElementIndex==4)
      {
        distances[0]=sqrt(((x_centroid-x_fourth_node)*(x_centroid-x_fourth_node))+((y_centroid-y_fourth_node)*(y_centroid-y_fourth_node))+((z_centroid-z_fourth_node)*(z_centroid-z_fourth_node)));
      }
      else if(CurrentMaxElementIndex==5)
      {
        distances[0]=sqrt(((x_centroid-x_fifth_node)*(x_centroid-x_fifth_node))+((y_centroid-y_fifth_node)*(y_centroid-y_fifth_node))+((z_centroid-z_fifth_node)*(z_centroid-z_fifth_node)));
      }
      else if(CurrentMaxElementIndex==6)
      {
        distances[0]=sqrt(((x_centroid-x_sixth_node)*(x_centroid-x_sixth_node))+((y_centroid-y_sixth_node)*(y_centroid-y_sixth_node))+((z_centroid-z_sixth_node)*(z_centroid-z_sixth_node)));
      }
      else if(CurrentMaxElementIndex==7)
      {
        distances[0]=sqrt(((x_centroid-x_seventh_node)*(x_centroid-x_seventh_node))+((y_centroid-y_seventh_node)*(y_centroid-y_seventh_node))+((z_centroid-z_seventh_node)*(z_centroid-z_seventh_node)));
      }
      else if(CurrentMaxElementIndex==8)
      {
        distances[0]=sqrt(((x_centroid-x_eighth_node)*(x_centroid-x_eighth_node))+((y_centroid-y_eighth_node)*(y_centroid-y_eighth_node))+((z_centroid-z_eighth_node)*(z_centroid-z_eighth_node)));
      }




      if(SecondMaxElementIndex==1)
      {
        distances[1]=sqrt(((x_centroid-x_first_node)*(x_centroid-x_first_node))+((y_centroid-y_first_node)*(y_centroid-y_first_node))+((z_centroid-z_first_node)*(z_centroid-z_first_node)));
      }
      else if(SecondMaxElementIndex==2)
      {
        distances[1]=sqrt(((x_centroid-x_second_node)*(x_centroid-x_second_node))+((y_centroid-y_second_node)*(y_centroid-y_second_node))+((z_centroid-z_second_node)*(z_centroid-z_second_node)));
      }
      else if(SecondMaxElementIndex==3)
      {
        distances[1]=sqrt(((x_centroid-x_third_node)*(x_centroid-x_third_node))+((y_centroid-y_third_node)*(y_centroid-y_third_node))+((z_centroid-z_third_node)*(z_centroid-z_third_node)));
      }
      else if(SecondMaxElementIndex==4)
      {
        distances[1]=sqrt(((x_centroid-x_fourth_node)*(x_centroid-x_fourth_node))+((y_centroid-y_fourth_node)*(y_centroid-y_fourth_node))+((z_centroid-z_fourth_node)*(z_centroid-z_fourth_node)));
      }
      else if(SecondMaxElementIndex==5)
      {
        distances[1]=sqrt(((x_centroid-x_fifth_node)*(x_centroid-x_fifth_node))+((y_centroid-y_fifth_node)*(y_centroid-y_fifth_node))+((z_centroid-z_fifth_node)*(z_centroid-z_fifth_node)));
      }
      else if(SecondMaxElementIndex==6)
      {
        distances[1]=sqrt(((x_centroid-x_sixth_node)*(x_centroid-x_sixth_node))+((y_centroid-y_sixth_node)*(y_centroid-y_sixth_node))+((z_centroid-z_sixth_node)*(z_centroid-z_sixth_node)));
      }
      else if(SecondMaxElementIndex==7)
      {
        distances[1]=sqrt(((x_centroid-x_seventh_node)*(x_centroid-x_seventh_node))+((y_centroid-y_seventh_node)*(y_centroid-y_seventh_node))+((z_centroid-z_seventh_node)*(z_centroid-z_seventh_node)));
      }
      else if(SecondMaxElementIndex==8)
      {
        distances[1]=sqrt(((x_centroid-x_eighth_node)*(x_centroid-x_eighth_node))+((y_centroid-y_eighth_node)*(y_centroid-y_eighth_node))+((z_centroid-z_eighth_node)*(z_centroid-z_eighth_node)));
      }


      int MinElementIndexDistances = std::min_element(distances.begin(),distances.end()) - distances.begin();
      if(MinElementIndexDistances == 0)
      {
      ConoraryFlowInterface::element_current_subdomain_lists[current_element]=CurrentMaxElementIndex;
      }
      else if(MinElementIndexDistances == 1)
      {
      ConoraryFlowInterface::element_current_subdomain_lists[current_element]=SecondMaxElementIndex;
      }




    }
    else if (maxValueAgain < 3)
    {
      ConoraryFlowInterface::element_current_subdomain_lists[current_element]=CurrentMaxElementIndex;
    }




  }

  else if(maxValue==2)
  {


    int CurrentMaxElementIndex = std::max_element(subdomain_holder.begin(),subdomain_holder.end()) - subdomain_holder.begin();

    subdomain_holder[CurrentMaxElementIndex]=0;

    int maxValueAgain = *max_element(subdomain_holder.begin(),subdomain_holder.end());
    int SecondMaxElementIndex = std::max_element(subdomain_holder.begin(),subdomain_holder.end()) - subdomain_holder.begin();

    if(maxValueAgain == 2)
    {

    std::vector<double> distances_four(8, 0.0);
    distances_four[0]=sqrt(((x_centroid-x_first_node)*(x_centroid-x_first_node))+((y_centroid-y_first_node)*(y_centroid-y_first_node))+((z_centroid-z_first_node)*(z_centroid-z_first_node)));
    distances_four[1]=sqrt(((x_centroid-x_second_node)*(x_centroid-x_second_node))+((y_centroid-y_second_node)*(y_centroid-y_second_node))+((z_centroid-z_second_node)*(z_centroid-z_second_node)));
    distances_four[2]=sqrt(((x_centroid-x_third_node)*(x_centroid-x_third_node))+((y_centroid-y_third_node)*(y_centroid-y_third_node))+((z_centroid-z_third_node)*(z_centroid-z_third_node)));
    distances_four[3]=sqrt(((x_centroid-x_fourth_node)*(x_centroid-x_fourth_node))+((y_centroid-y_fourth_node)*(y_centroid-y_fourth_node))+((z_centroid-z_fourth_node)*(z_centroid-z_fourth_node)));
    distances_four[4]=sqrt(((x_centroid-x_fifth_node)*(x_centroid-x_fifth_node))+((y_centroid-y_fifth_node)*(y_centroid-y_fifth_node))+((z_centroid-z_fifth_node)*(z_centroid-z_fifth_node)));
    distances_four[5]=sqrt(((x_centroid-x_sixth_node)*(x_centroid-x_sixth_node))+((y_centroid-y_sixth_node)*(y_centroid-y_sixth_node))+((z_centroid-z_sixth_node)*(z_centroid-z_sixth_node)));
    distances_four[6]=sqrt(((x_centroid-x_seventh_node)*(x_centroid-x_seventh_node))+((y_centroid-y_seventh_node)*(y_centroid-y_seventh_node))+((z_centroid-z_seventh_node)*(z_centroid-z_seventh_node)));
    distances_four[7]=sqrt(((x_centroid-x_eighth_node)*(x_centroid-x_eighth_node))+((y_centroid-y_eighth_node)*(y_centroid-y_eighth_node))+((z_centroid-z_eighth_node)*(z_centroid-z_eighth_node)));

    int MinElementIndex = std::min_element(distances_four.begin(),distances_four.end()) - distances_four.begin();

    if(MinElementIndex==0)
    {
    ConoraryFlowInterface::element_current_subdomain_lists[current_element]=first_subdomain;
    }
    else if(MinElementIndex==1)
    {
    ConoraryFlowInterface::element_current_subdomain_lists[current_element]=second_subdomain;
    }
    else if(MinElementIndex==2)
    {
    ConoraryFlowInterface::element_current_subdomain_lists[current_element]=third_subdomain;
    }
    else if(MinElementIndex==3)
    {
    ConoraryFlowInterface::element_current_subdomain_lists[current_element]=fourth_subdomain;
    }
    else if(MinElementIndex==4)
    {
    ConoraryFlowInterface::element_current_subdomain_lists[current_element]=fifth_subdomain;
    }
    else if(MinElementIndex==5)
    {
    ConoraryFlowInterface::element_current_subdomain_lists[current_element]=sixth_subdomain;
    }
    else if(MinElementIndex==6)
    {
    ConoraryFlowInterface::element_current_subdomain_lists[current_element]=seventh_subdomain;
    }
    else if(MinElementIndex==7)
    {
    ConoraryFlowInterface::element_current_subdomain_lists[current_element]=eighth_subdomain;
    }
    else
    {
      cout << "something has gone wrong" << endl;
    }


    }
    else if (maxValueAgain < 2)
    {
      ConoraryFlowInterface::element_current_subdomain_lists[current_element]=CurrentMaxElementIndex;
    }



  }

  else if(maxValue==1)
  {

      std::vector<double> distances_four(8, 0.0);
    distances_four[0]=sqrt(((x_centroid-x_first_node)*(x_centroid-x_first_node))+((y_centroid-y_first_node)*(y_centroid-y_first_node))+((z_centroid-z_first_node)*(z_centroid-z_first_node)));
    distances_four[1]=sqrt(((x_centroid-x_second_node)*(x_centroid-x_second_node))+((y_centroid-y_second_node)*(y_centroid-y_second_node))+((z_centroid-z_second_node)*(z_centroid-z_second_node)));
    distances_four[2]=sqrt(((x_centroid-x_third_node)*(x_centroid-x_third_node))+((y_centroid-y_third_node)*(y_centroid-y_third_node))+((z_centroid-z_third_node)*(z_centroid-z_third_node)));
    distances_four[3]=sqrt(((x_centroid-x_fourth_node)*(x_centroid-x_fourth_node))+((y_centroid-y_fourth_node)*(y_centroid-y_fourth_node))+((z_centroid-z_fourth_node)*(z_centroid-z_fourth_node)));
    distances_four[4]=sqrt(((x_centroid-x_fifth_node)*(x_centroid-x_fifth_node))+((y_centroid-y_fifth_node)*(y_centroid-y_fifth_node))+((z_centroid-z_fifth_node)*(z_centroid-z_fifth_node)));
    distances_four[5]=sqrt(((x_centroid-x_sixth_node)*(x_centroid-x_sixth_node))+((y_centroid-y_sixth_node)*(y_centroid-y_sixth_node))+((z_centroid-z_sixth_node)*(z_centroid-z_sixth_node)));
    distances_four[6]=sqrt(((x_centroid-x_seventh_node)*(x_centroid-x_seventh_node))+((y_centroid-y_seventh_node)*(y_centroid-y_seventh_node))+((z_centroid-z_seventh_node)*(z_centroid-z_seventh_node)));
    distances_four[7]=sqrt(((x_centroid-x_eighth_node)*(x_centroid-x_eighth_node))+((y_centroid-y_eighth_node)*(y_centroid-y_eighth_node))+((z_centroid-z_eighth_node)*(z_centroid-z_eighth_node)));


    int MinElementIndex = std::min_element(distances_four.begin(),distances_four.end()) - distances_four.begin();

    if(MinElementIndex==0)
    {
    ConoraryFlowInterface::element_current_subdomain_lists[current_element]=first_subdomain;
    }
    else if(MinElementIndex==1)
    {
    ConoraryFlowInterface::element_current_subdomain_lists[current_element]=second_subdomain;
    }
    else if(MinElementIndex==2)
    {
    ConoraryFlowInterface::element_current_subdomain_lists[current_element]=third_subdomain;
    }
    else if(MinElementIndex==3)
    {
    ConoraryFlowInterface::element_current_subdomain_lists[current_element]=fourth_subdomain;
    }
    else if(MinElementIndex==4)
    {
    ConoraryFlowInterface::element_current_subdomain_lists[current_element]=fifth_subdomain;
    }
    else if(MinElementIndex==5)
    {
    ConoraryFlowInterface::element_current_subdomain_lists[current_element]=sixth_subdomain;
    }
    else if(MinElementIndex==6)
    {
    ConoraryFlowInterface::element_current_subdomain_lists[current_element]=seventh_subdomain;
    }
    else if(MinElementIndex==7)
    {
    ConoraryFlowInterface::element_current_subdomain_lists[current_element]=eighth_subdomain;
    }
    else
    {
      cout << "something has gone wrong" << endl;
    }
  }

  int element_subdomain=ConoraryFlowInterface::element_current_subdomain_lists[current_element];
  double element_volume=elem->volume();

  ConoraryFlowInterface::total_subdomain_volumes[element_subdomain-1]=ConoraryFlowInterface::total_subdomain_volumes[element_subdomain-1]+element_volume;



}


}

for (int qq = 0; qq < 9; qq++)
{
MPI_Bcast(&ConoraryFlowInterface::total_subdomain_volumes[qq], 1, MPI_DOUBLE, broadcast_root, MPI_COMM_WORLD);
}

for (int qq = 0; qq < BoundaryConditions::total_elem_num; qq++)
{
MPI_Bcast(&ConoraryFlowInterface::element_current_subdomain_lists[qq], 1, MPI_INT, broadcast_root, MPI_COMM_WORLD);
}



}
else
{

  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  int broadcast_root = 0;

  if(my_rank == broadcast_root)
  {
  for (;  el != end_el; ++el)
  {
    Elem* elem = *el;
    dof_map.dof_indices(elem, dof_indices, 0);


    fe->reinit (elem);
    dof_map_J.dof_indices (elem, dof_indices_J);
    const unsigned int n_var_dofs = dof_indices_J.size();

    double J_cur = 0.0;
       for (unsigned int j=0; j<n_var_dofs; j++)
        {
         J_cur += phi[j][0]*J_vec[dof_indices_J[j]];
        }


    double element_volume=(J_cur)*(elem->volume());
    ConoraryFlowInterface::total_subdomain_volumes[ConoraryFlowInterface::element_current_subdomain_lists[elem->id()]-1]=ConoraryFlowInterface::total_subdomain_volumes[ConoraryFlowInterface::element_current_subdomain_lists[elem->id()]-1]+element_volume;
  }
  }

  for (int qq = 0; qq < 9; qq++)
  {
  MPI_Bcast(&ConoraryFlowInterface::total_subdomain_volumes[qq], 1, MPI_DOUBLE, broadcast_root, MPI_COMM_WORLD);
  }

}

if (ConoraryFlowInterface::output_coro_data_to_screen)
  {
    pout << "sub 1 volume   " << ConoraryFlowInterface::total_subdomain_volumes[0] << endl;
    pout << "sub 2 volume   " << ConoraryFlowInterface::total_subdomain_volumes[1] << endl;
    pout << "sub 3 volume   " << ConoraryFlowInterface::total_subdomain_volumes[2] << endl;
    pout << "sub 4 volume   " << ConoraryFlowInterface::total_subdomain_volumes[3] << endl;
    pout << "sub 5 volume   " << ConoraryFlowInterface::total_subdomain_volumes[4] << endl;
    pout << "sub 6 volume   " << ConoraryFlowInterface::total_subdomain_volumes[5] << endl;
    pout << "sub 7 volume   " << ConoraryFlowInterface::total_subdomain_volumes[6] << endl;
    pout << "sub 8 volume   " << ConoraryFlowInterface::total_subdomain_volumes[7] << endl;
    pout << "sub 9 volume   " << ConoraryFlowInterface::total_subdomain_volumes[8] << endl;
 }

    return;
}











void ConoraryFlowInterface::return_Pext_to_zero(EquationSystems* es,const double time,const double dt)
{


  int my_rank3;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank3);
  int broadcast_root3 = 0;

  if(my_rank3 == broadcast_root3)
  {

int steps_left;

  if (Tube::exp_reducer==0)
  {
    Tube::exp_reducer=1;
    steps_left=(1.1-time)/(dt);

    Tube::pressure_to_coronary_reduction[0]=(fabs(Tube::pressure_to_coronary_current[0]))/(steps_left);
    Tube::pressure_to_coronary_reduction[1]=(fabs(Tube::pressure_to_coronary_current[1]))/(steps_left);
    Tube::pressure_to_coronary_reduction[2]=(fabs(Tube::pressure_to_coronary_current[2]))/(steps_left);
    Tube::pressure_to_coronary_reduction[3]=(fabs(Tube::pressure_to_coronary_current[3]))/(steps_left);
    Tube::pressure_to_coronary_reduction[4]=(fabs(Tube::pressure_to_coronary_current[4]))/(steps_left);
    Tube::pressure_to_coronary_reduction[5]=(fabs(Tube::pressure_to_coronary_current[5]))/(steps_left);
    Tube::pressure_to_coronary_reduction[6]=(fabs(Tube::pressure_to_coronary_current[6]))/(steps_left);
    Tube::pressure_to_coronary_reduction[7]=(fabs(Tube::pressure_to_coronary_current[7]))/(steps_left);
    Tube::pressure_to_coronary_reduction[8]=(fabs(Tube::pressure_to_coronary_current[8]))/(steps_left);


  }






}


}
