  /***************************************************************************/
/*                                                                         */
/*  Program: arteries.h                                                    */
/*  Version: 2.0                                                           */
/*  By: Mette Olufsen, Math-Tech                                           */
/*  Date: 14. Jan. 1997                                                    */
/*                                                                         */
/*  This header file includes the definition of the global constants, and  */
/*  the definition of the class representing one blood vessel.             */
/*                                                                         */
/***************************************************************************/

// $Id: arteries.h,v 1.9 2005/07/07 22:19:51 heine Exp $
// Last updated on October 23, 2014 by M. Umar Qureshi

#ifndef _ARTERIESPD_H
#define _ARTERIESPD_H

#include <cstdio>
#include <cmath>

//#include <ibamr/app_namespaces.h>

#include <string>


// Global parameters imported from main.h

extern double   dyne_to_mmhg,time_holder,total_source,total_sink,net_sink,total_subdomain_volumes[100],pass_val, loop_time_arteries, transfer_aortic,transfer_aortic_previous,flow_added_this_timestep,total_flow_added,flow_added_this_timestep_in_dimensional_space,total_flow_added_in_dimensional_space, conv, rho, mu_pa, mu_pl, nu, Lr, Lr2, Lr3, g, q, Fr2,
                Re, p0, pmean, tmst, Period, Fcst, CO, COm,added_mass_sum,average_added_mass_sum,summed_source,summed_source_all_time, source_from_coro[100], sink_from_coro[100], net_from_coro[100],net_from_coro_in_dimensional_space[100],
				Deltat,pressure_to_coronary_current[100],pressure_to_coronary_previous[100],pressure_to_coronary_reduction[100], net_from_coro_averaged[100], net_from_coro_averaged_in_dimensional_space[100],
         lambda_investigate[100],dilatory_investigate[100],deviatoric_investigate[100],imp_holder_1[65537],imp_holder_2[65537],imp_holder_3[65537],imp_holder_4[65537],imp_holder_5[65537],imp_holder_6[65537],imp_holder_7[65537],imp_holder_8[65537],imp_holder_9[65537],
         imp_holder_1_period[65537],imp_holder_2_period[65537],imp_holder_3_period[65537],imp_holder_4_period[65537],imp_holder_5_period[65537],imp_holder_6_period[65537],imp_holder_7_period[65537],imp_holder_8_period[65537],imp_holder_9_period[65537];

extern int time_step_discrepency,conv_compute,number_of_subdomains,number_of_elements_in_subdomain[100],number_of_nodes_in_subdomain[100],source_number[100],sink_number[100],second_period_exp_update,exp_position,time_step_discrepency_user,counter_sec,exp_reducer,initial_or_integrated,current_position_in_coronary,current_or_previous,beat_tracker,nbrves_passer;

extern double  *fjac[18], *fj[8], *fmono[12], *ftrif[24];


//const char*  COqwqw_filename    = "corPress.dat";

// The class structure.
class Tube
{
public:

static int number_of_subdomains;
static int time_step_discrepency;
static int time_step_discrepency_user;
static int counter_sec;
static int exp_reducer;
static int exp_position;
static int second_period_exp_update;

static double pressure_to_coronary_current[100];
static double pressure_to_coronary_previous[100];
static double pressure_to_coronary_reduction[100];
static double source_from_coro[100];
static double sink_from_coro[100];
static double net_from_coro[100];
static double net_from_coro_averaged[100];
static double net_from_coro_in_dimensional_space[100];
static double net_from_coro_averaged_in_dimensional_space[100];

static double imp_holder_1[65537];
static double imp_holder_2[65537];
static double imp_holder_3[65537];
static double imp_holder_4[65537];
static double imp_holder_5[65537];
static double imp_holder_6[65537];
static double imp_holder_7[65537];
static double imp_holder_8[65537];
static double imp_holder_9[65537];

static double imp_holder_1_period[65537];
static double imp_holder_2_period[65537];
static double imp_holder_3_period[65537];
static double imp_holder_4_period[65537];
static double imp_holder_5_period[65537];
static double imp_holder_6_period[65537];
static double imp_holder_7_period[65537];
static double imp_holder_8_period[65537];
static double imp_holder_9_period[65537];

static double transfer_aortic;
static double transfer_aortic_previous;
static double total_flow_added;
static double flow_added_this_timestep;
static double total_flow_added_in_dimensional_space;
static double flow_added_this_timestep_in_dimensional_space;
static double loop_time_arteries;
static double pass_val;

static double dyne_to_mmhg;
static double rho;
static double conv;
static double mu1;
static double mu_pl;
static double Tper;
static double nu;
static double Fcst;
static double Lr;
static double Lr2;
static double Lr3;
static double g;
static double Re;
static double q;
static double Fr2;
static double Period;
static double p0;
static double Deltat;

static double tstart1;
static double tend1;
static double finaltime1;
static double onePeriodDt;
static double summed_source;
static double added_mass_sum;
static double average_added_mass_sum;
static double summed_source_all_time;

static double tstart_pass;
static double tend_pass;
static double time_holder;

static double total_source;
static double total_sink;
static double net_sink;

static int conv_compute;
static int beat_tracker;
static int initial_or_integrated;
static int current_or_previous;
static int plts;
static int tmstps;
static int nbrves_passer;
static int current_position_in_coronary;

static int source_number[100];
static int sink_number[100];
static int number_of_nodes_in_subdomain[100];
static int number_of_elements_in_subdomain[100];

  double L;                    // The length of the vessel
  double rtop, rbot;           // The top and bottom radii of the vessel
  Tube *LD, *MD, *RD;               // The left and right daughter-vessels. If
   //jam                            // these are set the Peripheral resistance
                               // should be set to zero.
  double rm;
  int pts;                  // The number of grid points per cm
  int init;
  double K_loss=0.0;
  double ff1, ff2, ff3, fa1, fa2, fa3, fv1, fv2, fv3;
  double asym, expo;
  double lrrA, lrrV;
  int SVPA, LVPA;                  // external pressure amplitude, in as mmHg
  int    counter=0;

  int N;                       // The number of grid points along the vessel
  double h;                    // The interval length of delta x
  double RLrb;                 // The peripheral resistance of the vessel

  double Ah05, Qh05;

  double *Qnew, *Qold, *Qh,    // The arrays needed to store data during
         *Anew, *Aold, *Ah,    // the numerical solution of the system.
         *R1, *R2, *R1h, *R2h,
         *S1, *S2, *S1h, *S2h,
         *Qprv, *Aprv,
         *pL, *y11, *y12, *y21, *y22, *yE11, *yE12, *yE21, *yE22, *QL, *Pout, *Z,
//       *Q0,       // uncomment for flow driven / comment for pressure driven
         *Ps,       // comment for flow input / uncomment for pressure driven
         *IPP,       // small vessels external pressure
	 *r0, *r0h,
	 *dr0dx, *dr0dxh,
	 *A0, *A0h, *wom,
	 *fr, *frh,
	 *dfrdr0, *dfrdr0h,
	 *p1, *p1h,
	 *dp1dr0, *dp1dr0h;

   double *imp;

     Tube (double Length,
           double topradius, double botradius,
           Tube *LeftDaughter, Tube *MidDaughter, Tube *RightDaughter, //jam
           double rmin, int points, int init,
           double f1, double f2, double f3,
           double fa1, double fa2, double fa3,
           double fv1, double fv2, double fv3,
           double asym, double expo,
           double lrrA, double lrrV, int SVPA, int LVPA);
                                                            // Constructor.
     ~Tube ();                                              // Destructor.


     // Prints P(x,t), A(x,t), F(x,t), or Q(x,t) for all x's and t's
     // along the tube. The argument offset makes sure that the vessel
     // is located with the right offset from the inlet.
       void printPxt (FILE *fd, double t, int offset, int sim_time, int vessel_id);

     // Defines P(x,A(x,t)).
  double P (int qLnb, int i, double A);

     // Defines dPdA(x,A(x,t)).
     double dPdA (int i, double A);

     // Defines dPdx1(x,A(x,t)).
     double dPdx1 (int i, double A);

     // Defines B(x,A(x,t)).
     double B (int i, double A);

     // Defines Bh(x,A(x,t)).
     double Bh (int i, double A);

     // Defines dBdx1(x,A(x,t)).
     double dBdx1 (int i, double A);

     // Defines dBdx1h(x,A(x,t)).
     double dBdx1h (int i, double A);

     // Defines dBdAh (x,A(x,t)).
     double dBdAh (int i, double A);

     // Defines d2BdAdxh (x, A(x,t));
     double d2BdAdxh (int i, double A);


     // Tests that the CFL-condition is valid throughout the tube.
     double CFL ();

     // Finds the flux acc. to sys. eq.
     double Rvec (int k, int i, int j, double Q, double A);

     // Finds the rhs. of system eq.
     double Svec (int k, int i, int j, double Q, double A);

    // initialises external pressure
     double IPP_init (double t, double k, double Period);

    // initialises fluid pressure // I think, JAM 3-may-18
     double Ps_init (double t, double k, double Period);


     // Steps through interior points.
  void step (double k, int i, double t, double t0);


     // Updates left bndry. This should only be done for the inlet tube.
  void bound_left (int qLnb, double t, double k, double Period);

     // Updates the matching boundary.
  void bound_match (int qLnb, double t, double k, double theta, double gamma, double t0);

     // Updates right bndry. This should only be done for terminal vessels.
     double c  (int i, double A); // The wave speed through aorta.
     // double Hp (int i, double Q, double A);
     // void poschar (double theta, double &qR, double &aR, double &cR, double &HpR);
     // void bound_right (int qLnb, double k, double theta, double t);

  void call_junc(int qLnb, double theta, double gamma, Tube *Arteries[], int parent);


     // Updates outflow.
  void bound_right (int qLnb, double t, double k, double Period);

     // In order to ensure a more efficient execution of the program the following
     // functions is made as in-line functions.

   // A function returning the Friction of the system. The definition of this
   // function is given according to the derivation in the mathematical model.
   // The constant cst, determines the amount of damping in the system.
   inline double F (double Q, double A)
   {
     double tmp1 = Fcst*M_PI*Q;
     double tmp2 = A*Re;
     double tmp3 = -tmp1/tmp2;
     return(tmp3);
     // return(-Fcst*M_PI*Q/(A*Re));
   }

   inline double dFdQ (double A)
   {
     return(-Fcst*M_PI/(A*Re));
   }

   inline double dFdA (double Q, double A)
   {
     return(Fcst*M_PI*Q/((A*A)*Re));
   }

   private:
     // The private function Q0 may only be accessed from the left boundary
     // function. It ensures a certain and given CO (defined in main.h).
     double Q0_init (double t, double k, double Period);
     //double Ps_init (double t, double k, double Period);
     double PL_init (double t, double k, double Period);
   };

   void solver (Tube *Arteries[], double tstart, double tend, double k, double t0);
   void looper (int qLnb,Tube *Arteries[]);
   #endif
