// Copyright (c) 2019-2020, Hao GAO
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

// APPLICATION INCLUDES
#include <BoundaryConditions.h>
#include <arteriesPD.h>
#include "MechanicsModel.h"
#include "ActiveContraction.h"

// LIBMESH INCLUDES
#include "libmesh/linear_implicit_system.h"
#include <libmesh/boundary_info.h>
#include <libmesh/point.h>
#include <libmesh/equation_systems.h>


#include <math.h>

// STATIC VARIABLES
// t_load and t_end_diastole is the same period
double BoundaryConditions::P_load, BoundaryConditions::t_load;
double BoundaryConditions::kappa;
BoundaryInfo* BoundaryConditions::boundary_info;
double BoundaryConditions::cycle_time_saved;
int BoundaryConditions::current_cycle;
double BoundaryConditions::P_current_loading;
double BoundaryConditions::P_current_loading_RV;
double BoundaryConditions::t_end_diastole;
double BoundaryConditions::t_period;
double BoundaryConditions::t_isovolumetric_relaxation_end = 0.0; //determined in phase determination
double BoundaryConditions::t_end_systole = 0.0; // determined in phase determination
double BoundaryConditions::cycle_no;
double BoundaryConditions::single_cycle_length;
double BoundaryConditions::diastole_systole_difference;
int BoundaryConditions::total_node_num;
int BoundaryConditions::total_elem_num;
int BoundaryConditions::total_endo_nodes;

// reading the endo surface nodes for volume calculation
std::vector<int>  BoundaryConditions::LV_endo_points_list;
int BoundaryConditions::LV_NoOfEndoNode;
std::vector< std::vector<double>  > BoundaryConditions::LV_endo_points;
std::vector< std::vector<double>  > BoundaryConditions::LV_all_points;

double BoundaryConditions::p_loading_bcast_from_root; //seems redundant
double BoundaryConditions::LV_volume; // calculated by CPU0 and then broadcast to others
double BoundaryConditions::whole_volume; // calculated by CPU0 and then broadcast to others
double BoundaryConditions::wall_volume; // calculated by CPU0 and then broadcast to others
double BoundaryConditions::LV_volume_old;
double BoundaryConditions::P_old_loading;
double BoundaryConditions::dP;

double BoundaryConditions::LV_CIsoContract;
double BoundaryConditions::LV_CIsoRelax;

double BoundaryConditions::dt; // will be updated at each step in the main loop

int BoundaryConditions::LV_states;
/* * LV_states:
 * 1: diastolic filling
 * 2: isovolumetric contraction
 * 3: systolic contraction
 * 4: isovolumetric relaxation
 * 5: relaxation to approach stress-free state
 * */
bool BoundaryConditions::diastolic_filling = false;
bool BoundaryConditions::isovolumetric_contraction = false;
bool BoundaryConditions::systolic_ejection = false;
bool BoundaryConditions::isovolumetric_relaxation = false;
bool BoundaryConditions::early_diastolic_relaxation = false;

double BoundaryConditions::pao_initial;
double BoundaryConditions::pwk_initial;
double BoundaryConditions::P_Ao;
double BoundaryConditions::P_AV_open;
double BoundaryConditions::P_AV_close;
double BoundaryConditions::LV_volume_ED;
double BoundaryConditions::LV_volume_ES;
double BoundaryConditions::Pwk_AV;
double BoundaryConditions::Pwk_AV_debug;
double BoundaryConditions::Pwk_AV_debug_old;
double BoundaryConditions::Pwk_LV_debug;
double BoundaryConditions::QLV_debug;

double BoundaryConditions::Ta_relax_factor = 1.0;
double BoundaryConditions::p_aortic_hold = 0.0;

double BoundaryConditions::remaining_steps;
double BoundaryConditions::p_ao_step_difference;
double BoundaryConditions::p_wk_step_difference;
double BoundaryConditions::loading_p_step_difference;

int BoundaryConditions::remaining_steps_rounded;
int BoundaryConditions::val_test;

namespace // private namespace
{

  //for aortic and loading pressure
  static const double Cwk_LV  = 1.75; //mL mmHg^-1
  static const double Rpwk_LV = 0.79*2.0; // mmHg mL^{-1} s

  //for loading pressure only
  static const double Rcwk_LV = 0.033; //mmHg mL{-1} s
}



// CLASS IMPLEMENTATION
/* to do list
double time_in_local_period = fmod(d_time, BoundaryConditions::t_period);
*/


double
BoundaryConditions::loading_pressure(
    double time)
{
    double P = 0.0;
    std::cout << "should never reach this here\n";
    return P;
}// loading_pressure

void
BoundaryConditions::loading_force_function(
    double& P,
    const libMesh::VectorValue<double>& n,
  	const libMesh::VectorValue<double>& N,
    const TensorValue<double>& /*FF*/,
    const libMesh::Point& /*x*/,
    const libMesh::Point& /*X*/,
    Elem* const elem,
    unsigned short int side,
    const vector<const vector<double>*>& /*system_data*/,
    const vector<const vector<libMesh::VectorValue<double> >*>& /*system_grad_var_data*/,
    double time,
    void* /*ctx*/)
{
    // 4096 ===>  epicardium
    // 4097 ===> endocardium
    // 4098 ===>        base
    // const vector<short int>& bdry_ids = boundary_info->boundary_ids(elem, side);
    vector<short int> bdry_ids;
    boundary_info->boundary_ids(elem, side, bdry_ids);
    if (find(bdry_ids.begin(),bdry_ids.end(),4097) != bdry_ids.end())
    {
        //P = loading_pressure(time)*1333.2239;
        P = P_current_loading*1333.2239;
    }
    else if (find(bdry_ids.begin(),bdry_ids.end(),5090) != bdry_ids.end() ||
             find(bdry_ids.begin(),bdry_ids.end(),5091) != bdry_ids.end())
    {  //for right ventricle
        P = P_current_loading*1333.2239/3.0; // only one third of the left ventricular pressure,
                                             // not used for this case
	  }
    else
    {
        P = 0.0;
    }
    return;
}// loading_force_function

void
BoundaryConditions::tether_force_function(
    VectorValue<double>& F,
    const libMesh::VectorValue<double>& n,
  	const libMesh::VectorValue<double>& N,
    const TensorValue<double>& /*FF*/,
    const libMesh::Point& x,
    const libMesh::Point& X,
    Elem* const elem,
    const unsigned short int side,
    const vector<const vector<double>*>& /*system_data*/,
    const vector<const vector<libMesh::VectorValue<double> >*>& /*system_grad_var_data*/,
    double /*time*/,
    void* /*ctx*/)
{
    F.zero();
    // 4096 ===>  epicardium
    // 4097 ===> endocardium
    // 4098 ===>        base
    // const vector<short int>& bdry_ids = boundary_info->boundary_ids(elem, side);
    vector<short int> bdry_ids;
    boundary_info->boundary_ids(elem, side, bdry_ids);
    if (find(bdry_ids.begin(),bdry_ids.end(),4098) != bdry_ids.end())
    {
        // approximately constrain motion to be in the radial direction at base
        VectorValue<double> r = x - libMesh::Point(0.0,0.0,x(2));  // radial displacement in current   configuration
        VectorValue<double> R = X - libMesh::Point(0.0,0.0,X(2));  // radial displacement in reference configuration
        VectorValue<double> r_perp = (R*r)*R/R.norm_sq() - r;
        F = kappa*r_perp;
        F(2) += kappa*(X(2)-x(2));
    }
    return;
}// tether_force_function



/*
    this is going to read endocardial surface points for volume calculation
    currently not used, instead the general one
*/
void
BoundaryConditions::readingPoints(MeshBase& mesh)
{
	   std::ifstream ifsendo("endoList.point");
	   int NoOfEndoNode, nodeEndoID;

	   ifsendo>>NoOfEndoNode;
	   LV_NoOfEndoNode = NoOfEndoNode;

	   LV_endo_points_list.resize(NoOfEndoNode);
	   LV_endo_points.resize(LV_NoOfEndoNode);
	   for (unsigned int i = 0; i< LV_NoOfEndoNode; i++)
	   {
		      LV_endo_points[i].resize(3);
	   }

	   unsigned int IDtemp=1; //reused from the initial defintion
	   unsigned int pIndex = 0;

	   //initialize end_points
	   while (!ifsendo.eof()& IDtemp<=NoOfEndoNode)
     {
			ifsendo>>nodeEndoID;
			IDtemp++;
			nodeEndoID = nodeEndoID - 1 ; //start from 0

		    LV_endo_points_list[pIndex]=nodeEndoID;

			pIndex = pIndex + 1;

		}








		printf("processor %d read %d points\n", mesh.processor_id(), pIndex);

	   return ;
}


/*
    this is going to read endocardial surface points for volume calculation
*/
void
BoundaryConditions::readingPointsGeneral(MeshBase& mesh,
                                    std::vector<int>& points_list,
                                    std::vector< std::vector<double> >& points,
                                    std::vector< std::vector<double> >& points_sec,
                                    int& NoOfPoints,
                                    std::string file_name)
{


fstream file3;
//file3.open(file_name);
file3.open("/xlwork1/scott/new_mesh_by_scott/endo_nodes.txt");
string line3;

while (getline( file3, line3,'\n'))
{
  istringstream templine(line3);
  string data;
  while (getline( templine, data,','))
  {
    points_list.push_back(atof(data.c_str()));
  }
}
file3.close();

	   points.resize(BoundaryConditions::total_endo_nodes);
	   for (unsigned int i = 0; i< BoundaryConditions::total_endo_nodes; i++)
	   {
		      points[i].resize(3);
	   }

     points_sec.resize(BoundaryConditions::total_node_num);
     for (unsigned int i = 0; i< BoundaryConditions::total_node_num; i++)
     {
          points_sec[i].resize(3);
     }


//		printf("processor %d read %d points\n", mesh.processor_id(), pIndex);
	  return ;
}


void
BoundaryConditions::updatePointsPosition(EquationSystems* equation_systems)
{
	     const MeshBase& mesh = equation_systems->get_mesh();
	     const unsigned int dim = mesh.mesh_dimension();
		   System& X_system = equation_systems->get_system<System>(IBFEMethod::COORDS_SYSTEM_NAME);

		   //libMesh::Node    node_ref;

		   for (unsigned int i = 0; i < LV_NoOfEndoNode; i++)
		   {
			   int nodeEndoID = LV_endo_points_list[i]; //start from 0

			   const libMesh::Node& node_ref = mesh.node_ref(nodeEndoID);

			   LV_endo_points[i][0] = node_ref(0);
			   LV_endo_points[i][1] = node_ref(1);
			   LV_endo_points[i][2] = node_ref(2);



		   }

		   return ;
		   //double volume = tetVolumeCalculationByPoints(endo_points, BoundaryConditions::LV_NoOfEndoNode);
}


void BoundaryConditions::updatePointsPositionGeneral(EquationSystems* equation_systems,
                                                    std::vector<int> & points_list,
                                     std::vector< std::vector<double> >& points_coor,
                                     std::vector< std::vector<double> >& points_coor_sec,
                                                    int NoOfPoints)
{
	     const MeshBase& mesh = equation_systems->get_mesh();
	     const unsigned int dim = mesh.mesh_dimension();

		   System& X_system = equation_systems->get_system<System>(IBFEMethod::COORDS_SYSTEM_NAME);
		   const DofMap& X_dof_map = X_system.get_dof_map();
		   std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);

		   X_system.solution->localize(*X_system.current_local_solution);
       NumericVector<double>& X_data = *(X_system.current_local_solution);
       X_data.close();

       //copy data to one processor, default is 0 processor
       std::vector<double> X_data_vec;
	     X_system.solution->localize_to_one(X_data_vec);
       const unsigned int X_sys_num = X_system.number();

		    //printf("working on the updating, processor: %d\n", X_system.processor_id());


       MPI_Barrier(X_system.comm().get());
       if (0 == X_system.processor_id())
       {
 			   //printf("updating node position for vol cal processor 0\n");
 			   for (unsigned int i = 0; i < BoundaryConditions::total_node_num; i++)
 			   {
 				   unsigned int nodeID = i; //start from 0

 		        const libMesh::Node& node_ref = mesh.node_ref(nodeID); //const is needed otherwise it just abort
             //printf("obtain node ref : %d\n", nodeEndoID);

 				   const int x_dof_index = node_ref.dof_number(X_sys_num, 0, 0);
 				   const int y_dof_index = node_ref.dof_number(X_sys_num, 1, 0);
 				   const int z_dof_index = node_ref.dof_number(X_sys_num, 2, 0);

 				   points_coor_sec[i][0] = X_data_vec[x_dof_index];
 				   points_coor_sec[i][1] = X_data_vec[y_dof_index];
 				   points_coor_sec[i][2] = X_data_vec[z_dof_index];


 			   }
 		   }

		   return;
}


/*
updating pressure according to the LV contracting state, the updated pressure
will be save in P_current_loading, and also return to P.
*/
double BoundaryConditions::lump_para_pressure(EquationSystems* equation_systems,double time)
{
	double P = 0.0;

	if (LV_states == 1 && (time-cycle_time_saved) <= t_load) //diastolic filling
	{
		P = P_load*(time-cycle_time_saved)/(t_load);

//------------------------------------------------------------------------------
        //three-element windkessel model
        double Pwk_next = 0.0;
     		double Pwk_bar = 0.0;
        Pwk_bar=Pwk_AV*(1.0-((dt)/(Cwk_LV*Rpwk_LV)));
        Pwk_next=Pwk_AV-(((dt)/(Cwk_LV))*((Pwk_bar+Pwk_AV)/(2.0*Rpwk_LV)));
        //updating P_Ao
     		P_Ao = Pwk_next;
     		//updating Pwk
     		Pwk_AV = Pwk_next;
//------------------------------------------------------------------------------


	}
	else if ((time-cycle_time_saved) > t_load && LV_states == 2) //isovolumetric contraction
	{
		if (LV_volume < LV_volume_old )
		 {
			 P = P_current_loading + ((LV_volume)-(LV_volume_old))/(LV_CIsoContract);
	   }
		 else
		 {
			 P = P_current_loading; //if the LV is expanding, then pressure is not change, waiting for the contraction
		 }

//------------------------------------------------------------------------------
    //three-element windkessel model
    double Pwk_next = 0.0;
 		double Pwk_bar = 0.0;
    Pwk_bar=Pwk_AV*(1.0-((dt)/(Cwk_LV*Rpwk_LV)));
    Pwk_next=Pwk_AV-(((dt)/(Cwk_LV))*((Pwk_bar+Pwk_AV)/(2.0*Rpwk_LV)));
    //updating P_Ao
 		P_Ao = Pwk_next;
 		//updating Pwk
 		Pwk_AV = Pwk_next;
//------------------------------------------------------------------------------

	}
	else if ((time-cycle_time_saved) > t_load && LV_states == 3)
	{
		//calcualte volume flow rate
		double dv = LV_volume - LV_volume_old;
		double Qlv = -dv/dt; //aortic flow rate is the oppsite side of the LV volume change rate

	  double Pwk_next = 0.0;
		double Pwk_bar = 0.0;

		//three-element windkessel model
    Pwk_bar=Pwk_AV+(((dt)/(Cwk_LV))*(Qlv-((Pwk_AV)/(Rpwk_LV))));
    Pwk_next=Pwk_AV+(((dt)/(Cwk_LV))*((Qlv)-((Pwk_bar+Pwk_AV)/(2.0*Rpwk_LV))));

    //updating P_Ao
		P_Ao = Rcwk_LV*Qlv + Pwk_next;
    P=P_Ao;

		//updating Pwk
		Pwk_AV = Pwk_next;

	}
	else if ( LV_states == 4 && (time-cycle_time_saved) >= t_end_systole)
	{
    if (LV_volume > LV_volume_old )
		  {
        P = P_current_loading + ((LV_volume) - (LV_volume_old))/LV_CIsoRelax;
      }
    else
 		 {
 			 P = P_current_loading; //if the LV is expanding, then pressure is not change,
 			                        //waiting for the contraction
 		 }

//------------------------------------------------------------------------------
      //three-element windkessel model
      double Pwk_next = 0.0;
      double Pwk_bar = 0.0;
      Pwk_bar=Pwk_AV*(1.0-((dt)/(Cwk_LV*Rpwk_LV)));
      Pwk_next=Pwk_AV-(((dt)/(Cwk_LV))*((Pwk_bar+Pwk_AV)/(2.0*Rpwk_LV)));
      //updating P_Ao
      P_Ao = Pwk_next;
      //updating Pwk
      Pwk_AV = Pwk_next;
//------------------------------------------------------------------------------

	}
  else if ((time-cycle_time_saved) > t_load && LV_states == 5 && P_current_loading >= 0)
  {
// reduce pressure P to be zero
//------------------------------------------------------------------------------
// designed to gradually reduced p_ao and p_wk to their original starting value of 80.0
//during LV_states=5. Both should always start out greater than 80.0 at the end of lvstates=4
if (val_test==0)
{
remaining_steps=(t_period-(time-cycle_time_saved))/(dt);
remaining_steps_rounded=round(remaining_steps);

p_ao_step_difference=(P_Ao-pao_initial)/remaining_steps_rounded;
p_wk_step_difference=(Pwk_AV-pwk_initial)/remaining_steps_rounded;

pout << "checking that this part of the code was only entered once at the beginning of states 5  " << "\n";

remaining_steps=(0.1)/(dt);
remaining_steps_rounded=round(remaining_steps);

loading_p_step_difference=P_current_loading/remaining_steps_rounded;


val_test=1;

}


if (P_current_loading > 0.0001)
{
P =  P_current_loading-loading_p_step_difference; //decrease P_current with a coefficient
}
else
{
P = 0.0;
}

P_Ao=P_Ao-p_ao_step_difference;
Pwk_AV=Pwk_AV-p_wk_step_difference;


//------------------------------------------------------------------------------



    //reset back to the start of the period, by this time P_current_loading will be zero so there is no need to test this
    //as a requirement to move on to the next cycle
      if (time-cycle_time_saved>=BoundaryConditions::single_cycle_length)
      {
      cycle_time_saved=current_cycle*BoundaryConditions::t_period;
      current_cycle=current_cycle+1;

    //set all the following back to start of cycle state, set diastolic_filling to be true as that is the
    //phase to be entered.
      diastolic_filling = true;
      isovolumetric_contraction = false;
      systolic_ejection = false;
      isovolumetric_relaxation = false;
      early_diastolic_relaxation = false;

      t_end_systole = 0.0;
      p_aortic_hold=0.0;

      Ta_relax_factor = 1.0;
      LV_states = 1;

      LV_volume_ED=0.0;
      LV_volume_ES=0.0;

      t_isovolumetric_relaxation_end=0.0;
      remaining_steps=0.0;
      p_ao_step_difference=0.0;
      p_wk_step_difference=0.0;
      val_test=0;
      remaining_steps_rounded=0;

      ActiveContraction::reinitialize_variables(equation_systems);

      P_Ao=pao_initial;
      Pwk_AV=pwk_initial;

      }

  }


	dP = P - P_current_loading; // for printing out
	P_current_loading = P; // saving for loading force function


//update the value into the class Tube so that it can be more easily transferred to arteriesPD.C
Tube::transfer_aortic_previous=Tube::transfer_aortic;
Tube::transfer_aortic=P_Ao;



	return P;
}


 void BoundaryConditions::phase_determination(
        double time)
{
		if ((time-cycle_time_saved) <= t_load)
		{
			LV_states = 1; // in diastolic filling
			LV_volume_ED = LV_volume; // to save the end-diastolic volume
      diastolic_filling = true;
      Ta_relax_factor = 1.0; // no change
      p_aortic_hold=Tube::transfer_aortic;
		}
		else if ((time-cycle_time_saved) > t_load && P_current_loading <= p_aortic_hold &&
             diastolic_filling == true && systolic_ejection == false)
		{
			LV_states = 2; //isotropic contraction
      isovolumetric_contraction = true;
      Ta_relax_factor = 1.0; // no change
      p_aortic_hold=Tube::transfer_aortic;
		}
		else if ( ((time-cycle_time_saved) > t_load) && (P_current_loading > p_aortic_hold) &&
              (LV_volume_old>LV_volume) && isovolumetric_contraction == true &&
              isovolumetric_relaxation == false)
		{
			LV_states = 3; // in systolic contraction
			LV_volume_ES = LV_volume;
      systolic_ejection = true;
      t_end_systole = time-cycle_time_saved;
      Ta_relax_factor = 1.0; // no change
		}
		else if ( ((time-cycle_time_saved) > t_end_systole) && early_diastolic_relaxation==false &&
               systolic_ejection == true && (time-cycle_time_saved)-t_end_systole<0.1) // will allow 0.1s for isovolumetric relaxation
		{
			LV_states = 4; //in isovolumetric relaxation
      isovolumetric_relaxation = true;
      t_isovolumetric_relaxation_end = time-cycle_time_saved;
      Ta_relax_factor = 1 - ((time-cycle_time_saved)-t_end_systole)/0.1 + 0.01; // how to update here, going down to 0.01
		}
		else if ((time-cycle_time_saved) > t_isovolumetric_relaxation_end && diastolic_filling &&
             isovolumetric_contraction && systolic_ejection && isovolumetric_relaxation)
		{
			LV_states = 5;
      early_diastolic_relaxation = true;
      Ta_relax_factor = 0.0;
		}


		return;
}


void BoundaryConditions::print_info(double time)
{
    pout <<                                                       endl;
    pout << "print out by BoundaryConditions Class"             <<endl;
    pout << "Simulation time is " << time                  <<endl;
    pout << "Time step is       " << dt                         <<endl;
    pout << "Current cycle is       " << current_cycle          <<endl;
    pout << "Checking the cycle time saved variable " << cycle_time_saved       <<endl;
    pout << "Loading pressure : " << P_current_loading<< " mmHg" << endl;
    pout << "delta pressure:    " << dP <<" mmHg"<<endl;
    pout << "Pwk pressure :" << Pwk_AV <<" mmHg"<<endl;
    pout << "P_Ao pressure :" << P_Ao <<" mmHg"<<endl;
    pout << "checking p ao step difference " << p_ao_step_difference                  <<endl;
    pout << "checking p wk step difference " << p_wk_step_difference                  <<endl;
    pout << "Current whole volume   :" << whole_volume << " mL"         << endl;
    pout << "Current wall volume   :" << wall_volume << " mL"         << endl;
    pout << "Current LV (cavity) volume   :" << LV_volume << " mL"         << endl;
    pout << "Previous LV volume  :" << LV_volume_old<<" mL"       << endl;
    pout << "LV volume change    :" << LV_volume - LV_volume_old <<" mL"<<endl;
    pout << "LV flow rate        :" << (-LV_volume - LV_volume_old)/dt<<" mL s^-1"<<endl;
    pout << "End_diastolic volume is: "<< LV_volume_ED <<endl;
    pout << "t_load is : "<<t_load<<"\t t_end_diastole is: "<<t_end_diastole<<endl;
    pout << "P_load is : "<<P_load<<"\t P_AV_open is: " <<P_AV_open<<endl;
    pout << "Cwk_LV is  : "<<Cwk_LV<<endl;
    pout << "Rcwk_LV is : "<<Rcwk_LV<<endl;
    pout << "Rpwk_LV is : "<<Rpwk_LV<<endl;
    pout << "Checking p aortic hold :" << p_aortic_hold <<endl;
    pout << "LV state is : "         << LV_states <<endl;
    pout << "diastolic_filling: "<<diastolic_filling<<";   isovolumetric_contraction: "<<isovolumetric_contraction <<endl;
    pout << "systolic_ejection: "<<systolic_ejection<<";   isovolumetric_relaxation: "<<isovolumetric_relaxation <<endl;
    pout << "early_diastolic_relaxation: "<<early_diastolic_relaxation<<endl;
    pout << "Ta_relax_factor: " << Ta_relax_factor << endl;
    if (isovolumetric_relaxation == 1)
    { pout << "t_end_systole: "<< t_end_systole <<endl;}
    //pout << "LV_CIsoContract is : "<<LV_CIsoContract<<endl;
    //pout << "LV_CSystole is: "<<LV_CSystole <<"\t LV_RSystole is: "<<LV_RSystole<<endl;
    //pout << "LV_CIsoRelax is: "<<LV_CIsoRelax<<endl;
    pout <<                                                       endl;
    pout <<                                                       endl;
}





void BoundaryConditions::solve_wall_volume(EquationSystems* equation_systems)
{

const MeshBase& mesh = equation_systems->get_mesh();
const unsigned int dim = mesh.mesh_dimension();
//----------------------------------------------------------------------------------

LinearImplicitSystem& J_system = equation_systems->get_system<LinearImplicitSystem>(MechanicsModel::J_system_name);
J_system.solution->localize(*J_system.current_local_solution);
NumericVector<double>& J_data = *(J_system.current_local_solution);
J_data.close();
const DofMap & dof_map_J = J_system.get_dof_map();

//----------------------------------------------------------------------------------
FEType fe_type = dof_map_J.variable_type(0);
UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
QGauss qrule (dim, CONSTANT);
fe->attach_quadrature_rule (&qrule);
UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));
QGauss qface(dim-1, CONSTANT);
//----------------------------------------------------------------------------------
fe_face->attach_quadrature_rule (&qface);
const std::vector<std::vector<Real>> & phi = fe->get_phi();

    double wall_volume_inter=0.0;

    vector<double> J_vec;
    J_data.localize(J_vec);

   std::vector<dof_id_type> dof_indices_J;

   int my_rank3;
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank3);
   int broadcast_root3 = 0;

   if(my_rank3 == broadcast_root3)
   {

    //traditional way of looping all elements to be compatible with 1.0 version
    MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_elements_end();
    for (;  el != end_el; ++el)
        //   for (const auto & elem : mesh.active_local_element_ptr_range()) //this is a new feature
    {
        Elem* elem = *el;

    fe->reinit (elem);
    dof_map_J.dof_indices (elem, dof_indices_J);
    const unsigned int n_var_dofs = dof_indices_J.size();

    double J_cur = 0.0;
   for (unsigned int j=0; j<n_var_dofs; j++)
    {
     J_cur += phi[j][0]*J_vec[dof_indices_J[j]];
    }


        double deformed_elem_volume=0.0;

        deformed_elem_volume=(J_cur)*(elem->volume());
        wall_volume_inter=wall_volume_inter+deformed_elem_volume;
    }


wall_volume=wall_volume_inter;

}


    return;
}





void P_Ao_compute(double pwk_prev,double q_prev)
{




		return;

}
