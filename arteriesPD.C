// PRESSURE DRIVEN
// analytic jacobian inverse
// mono, bi, and trifurcations are permissable
/*************************************************************************/
/*  Program: arteries.C                                                  */
/*  Version: 2.?                                                         */
/*  By: Mette Olufsen, Math-Tech, M. Umar Qureshi, UofG                  */
/*  Date: 14. Jan. 1997                                                  */
/*  This version: 28. May 2020. By J. A. Mackenzie, UofG                 */
/*                                                                       */
/*  This module can predict the flow and pressure in a tree of elastic   */
/*  vessels as described in IMFUFATEKST NR 297, and D2.1-4. The          */
/*  dependencies of the vessels in the tree must be specified in the     */
/*  main module according to the tree in question (for further details   */
/*  see documentation in problem.default.cxx).                           */
/*  This module includes all the functions needed to solve the system    */
/*  of equations. That is the description of all functions in the class  */
/*  containing the vessel (for further details see arteries.h), and in   */
/*  particular the functions needed to solve the system of equations     */
/*  numerically.                                                         */
/*                                                                       */
/*  The module is dependent on the utilities in tools.C, and             */
/*  their corresponding h-files, and also arteries.h that includes the   */
/*  declaration of the vessel-object.                                    */
/*                                                                       */
/*  NEW TUBE LAW WITH EXTERNAL PRESSURE                                  */
/*  J. A. M. 20/05/2020                                                  */
/*                                                                       */
/*  All places where we computate area and pressure should have changed, */
/*  as the tube law has changed. Bifurcation and boundary matching       */
/*  conditions should also change.                                       */
/*  Transmural and external admittance matrices are computed in seperate */
/*  f90 files, so called seperately.                                     */
/*************************************************************************/
//
//  85 character lines for comfortable laptop reading. JAM
//
//  If tagged FlP, the quantity either is or uses fluid pressure,
//  if tagged TmP, it's transmural pressure.
//  If in doubt, it's probably FlP.
//  imp index in bound_match corrected by JAM july 2020.

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <assert.h>

#include "tools.h"
#include "arteriesPD.h"
#include "nr3.h"
#include "ludcmp.h"
#include "qrdcmp.h"
#include "roots_multidim.h"
#include "junction.h"

using namespace std;

extern int nbrves;
extern int tmstps;
//extern char* CO_filename;        // is the in data
//extern char* large_filename;     // is the in data // JAM. 05/2020
//extern char* small_filename;     // is the in data // JAM. 05/2020



double Tube::source_from_coro[100];
double Tube::sink_from_coro[100];
double Tube::net_from_coro[100];
double Tube::net_from_coro_averaged[100];
double Tube::net_from_coro_in_dimensional_space[100];
double Tube::net_from_coro_averaged_in_dimensional_space[100];
double Tube::pressure_to_coronary_current[100];
double Tube::pressure_to_coronary_previous[100];
double Tube::pressure_to_coronary_reduction[100];

int Tube::number_of_elements_in_subdomain[100];
int Tube::number_of_nodes_in_subdomain[100];
int Tube::source_number[100];
int Tube::sink_number[100];

double Tube::imp_holder_1[65537];
double Tube::imp_holder_2[65537];
double Tube::imp_holder_3[65537];
double Tube::imp_holder_4[65537];
double Tube::imp_holder_5[65537];
double Tube::imp_holder_6[65537];
double Tube::imp_holder_7[65537];
double Tube::imp_holder_8[65537];
double Tube::imp_holder_9[65537];

double Tube::imp_holder_1_period[65537];
double Tube::imp_holder_2_period[65537];
double Tube::imp_holder_3_period[65537];
double Tube::imp_holder_4_period[65537];
double Tube::imp_holder_5_period[65537];
double Tube::imp_holder_6_period[65537];
double Tube::imp_holder_7_period[65537];
double Tube::imp_holder_8_period[65537];
double Tube::imp_holder_9_period[65537];


int Tube::number_of_subdomains;
int Tube::exp_reducer;
int Tube::counter_sec;
int Tube::time_step_discrepency;
int Tube::time_step_discrepency_user;
int Tube::exp_position;
int Tube::second_period_exp_update;


double Tube::average_added_mass_sum;
double Tube::added_mass_sum;
double Tube::summed_source;
double Tube::summed_source_all_time;
double Tube::total_flow_added;
double Tube::flow_added_this_timestep;
double Tube::total_flow_added_in_dimensional_space;
double Tube::flow_added_this_timestep_in_dimensional_space;
double Tube::dyne_to_mmhg;
double Tube::loop_time_arteries;
double Tube::transfer_aortic;
double Tube::transfer_aortic_previous;
double Tube::mu_pl;
double Tube::rho;
double Tube::conv;
double Tube::mu1;
double Tube::Tper;
double Tube::nu;
double Tube::Fcst;
double Tube::Lr;
double Tube::Lr2;
double Tube::Lr3;
double Tube::g;
double Tube::Re;
double Tube::q;
double Tube::Fr2;
double Tube::Period;
double Tube::p0;
double Tube::Deltat;
double Tube::pass_val;

double Tube::tstart1;
double Tube::tend1;
double Tube::finaltime1;
double Tube::onePeriodDt;

double Tube::tstart_pass;
double Tube::tend_pass;
double Tube::time_holder;
double Tube::total_source;
double Tube::total_sink;
double Tube::net_sink;

int Tube::conv_compute;
int Tube::plts;
int Tube::tmstps;
int Tube::initial_or_integrated;
int Tube::current_or_previous;
int Tube::beat_tracker;
int Tube::nbrves_passer;
int Tube::current_position_in_coronary;

// A f90 subroutine that determines the impedance at the root of a
// structured tree. Used as BC for the tree consisting of the large vessels
extern "C"
void impedance_driver_(int * tmstps, double * Period, double * rho,
    double * mu,
    double * r_root, double * rmin,
    double * y11, double * y12,
    double * y21, double * y22,
    double * Lr, double * q, double * g,
    double * fa1, double * fa2, double * fa3,
    double * fv1, double * fv2, double * fv3,
    double * asym, double * expo,
    double * lrrA, double * lrrV);



//  Methods of class Tube, see arteries.h for description of this.
//
//  The constructor. When an object is made this function will
//  initialize all the attributes of the specific tube. The parameters
//  for the length of the specific vessel, the top and bottom radii,
//  and if applicable pointers to the daughter arteries will be
//  initialized according to the actual parameters passed in the call.
//  If the tube is terminal then the peripheral resistance must be set,
//  and the daughter vessels should be NIL. Otherwise the pointers to
//  the daughter vessels must be given. Further all the work arrays are
//  declared and initialized, and the initial condition for the system
//  equations is applied.

Tube :: Tube (double Length,
              double topradius, double botradius,
              Tube *LeftDaughter, Tube *MidDaughter, Tube *RightDaughter,
              double rmin, int points, int init,
              double f1, double f2, double f3,
              double fa1, double fa2, double fa3,
              double fv1, double fv2, double fv3,
              double asym, double expo,
              double lrrA, double lrrV, int SVPA, int LVPA):
    L(Length), rtop(topradius), rbot(botradius),
    LD(LeftDaughter), MD(MidDaughter), RD(RightDaughter),
    rm(rmin), pts(points), init(init),
    ff1(f1), ff2(f2), ff3(f3), SVPA(SVPA), LVPA(LVPA)
    {

      // Initialization of the basic parameters
      N = int(pts * L);
      h = 1.0 / pts / Lr;


      // Declaration and Initialization of the needed intermediate arrays.
      Qnew = new double[N + 1];
      Anew = new double[N + 1];
      Qold = new double[N + 1];
      Aold = new double[N + 1];
      Qprv = new double[N + 1];
      Aprv = new double[N + 1];
      R1 = new double[N + 1];
      R2 = new double[N + 1];
      S1 = new double[N + 1];
      S2 = new double[N + 1];
      r0 = new double[N + 1];
      r0h = new double[N + 2];
      dr0dx = new double[N + 1];
      dr0dxh = new double[N + 2];
      wom = new double[N + 1];
      A0 = new double[N + 1];
      A0h = new double[N + 2];
      fr = new double[N + 1];
      frh = new double[N + 2];
      dfrdr0 = new double[N + 1];
      dfrdr0h = new double[N + 2];
      p1 = new double[N + 1];
      p1h = new double[N + 2];
      dp1dr0 = new double[N + 1];
      dp1dr0h = new double[N + 2];
      Ah = new double[N];
      Qh = new double[N];
      R1h = new double[N];
      R2h = new double[N];
      S1h = new double[N];
      S2h = new double[N];
      pL = new double[tmstps];
      y11 = new double[tmstps];
      y12 = new double[tmstps];
      y21 = new double[tmstps];
      y22 = new double[tmstps];
      double rgLr = 4.0 / 3.0 / rho / g / Lr;
      double rgLr2 = 4.0 / 3.0 / rho / g / Lr2;

      // Vessel geometry is tabulated and initial conditions are applied
      for (int i = 0; i <= N; i++) {
          r0[i] = rtop * exp(i * log(rbot / rtop) / N) / Lr;
          r0h[i] = rtop * exp((i - 0.5) * log(rbot / rtop) / N) / Lr;
          dr0dx[i] = (log(rbot / rtop) / h / N) * r0[i];
          dr0dxh[i] = (log(rbot / rtop) / h / N) * r0h[i];
          wom[i] = r0[i] * sqrt(2.0 * M_PI / Period / nu);
          A0[i] = M_PI * sq(r0[i]);
          A0h[i] = M_PI * sq(r0h[i]);
          fr[i] = (ff1 * exp(ff2 * r0[i]) + ff3) * rgLr;
          frh[i] = (ff1 * exp(ff2 * r0h[i]) + ff3) * rgLr;
          dfrdr0[i] = ff1 * ff2 * exp(ff2 * r0[i]) * rgLr2;
          dfrdr0h[i] = ff1 * ff2 * exp(ff2 * r0h[i]) * rgLr2;
          p1[i] = fr[i] / M_PI;
          p1h[i] = frh[i] / M_PI;
          dp1dr0[i] = dfrdr0[i] / M_PI;
          dp1dr0h[i] = dfrdr0h[i] / M_PI;
          Qnew[i] = 1.0;
          Anew[i] = A0[i];
      }

      r0h[N + 1] = rtop * exp((N + 0.5) * log(rbot / rtop) / N) / Lr;
      dr0dxh[N + 1] = log(rbot / rtop) / h / N * r0h[N + 1];
      A0h[N + 1] = M_PI * sq(r0h[N + 1]);
      frh[N + 1] = (ff1 * exp(ff2 * r0h[N + 1]) + ff3) * rgLr;
      dfrdr0h[N + 1] = ff1 * ff2 * exp(ff2 * r0h[N + 1]) * rgLr2;
      p1h[N + 1] = frh[N + 1] / M_PI;
      dp1dr0h[N + 1] = dfrdr0h[N + 1] / M_PI;
      Ah05 = A0[N];
      Qh05 = 1.0;

//  =======================================================================





    //  ====================  PRESSURE AS INLET CONDITION  ====================

    if (init == 1) {
        Ps = new double[tmstps + 1];
        FILE * fi = fopen("/xlwork1/scott/full_model_nov_22/Data/aortic_128_mesh.dat", "r"); // Reading in pressure
        for (int i = 0; i <= tmstps; i++) {
            fscanf(fi, "%lf", & Ps[i]);
            Ps[i] = (Ps[i]) * conv / rho / g / Lr; // non-dimensionalise if required
        }
        fclose(fi);
    }


    //  =======================================================================



    //  =================  PRESSURE AT OUTLET OF LARGE VEINS  =================
    if (init == 2) {
        Pout = new double[tmstps + 1];
        for (int i = 0; i <= tmstps; i++) {
            Pout[i] = 5.0 * conv / rho / g / Lr;
        }
    }
    //  ====================================================================




    //  sinusoidal pressure profile applied to large vessels
    IPP = new double[tmstps + 1];
    for (int i = 0; i <= tmstps; i++) {
//        IPP[i] = LVPA * sin(i * 2 * M_PI / tmstps) * conv / rho / g / Lr; // sinusoidal
    }

    // In case of an end-tube evaluate the impedances for the boundary
    // condition. This is done by calling the f90 routine root_imp which
    // calculates the impedance at the root of a structured tree. The
    // underscores is sensitive to the compiler but can be seen at the bottom
    // of the file root_imp.o.

    if (init == 3)
    {


                                imp = new double[tmstps+1];
                                RD -> imp = new double[tmstps + 1];

                                if (SVPA == 1)
                                {
                                    FILE *fi1 = fopen ("/xlwork1/scott/full_model_nov_22/startup_files/pressure_128_mesh_sub1.dat", "r"); // Reading in pressure
                                    for (int j=0; j<=tmstps; j++)
                                    {
                                    fscanf(fi1,"%lf",&imp[j]);

                                    //imp[j] = (imp[j]/(rho*g*Lr/conv));
                                    Tube::imp_holder_1[j]=imp[j];

                                    }
                                    RD -> imp = imp;
                                    fclose(fi1);
                                }
                                if (SVPA == 2)
                                {
                                    FILE *fi1 = fopen ("/xlwork1/scott/full_model_nov_22/startup_files/pressure_128_mesh_sub2.dat", "r"); // Reading in pressure
                                    for (int j=0; j<=tmstps; j++)
                                    {
                                    fscanf(fi1,"%lf",&imp[j]);

                                    //imp[j] = (imp[j]/(rho*g*Lr/conv));
                                    Tube::imp_holder_2[j]=imp[j];
                                    }
                                    RD -> imp = imp;
                                    fclose(fi1);
                                }
                                if (SVPA == 3)
                                {
                                    FILE *fi1 = fopen ("/xlwork1/scott/full_model_nov_22/startup_files/pressure_128_mesh_sub3.dat", "r"); // Reading in pressure
                                    for (int j=0; j<=tmstps; j++)
                                    {
                                    fscanf(fi1,"%lf",&imp[j]);

                                    //imp[j] = (imp[j]/(rho*g*Lr/conv));
                                    Tube::imp_holder_3[j]=imp[j];
                                    }
                                    RD -> imp = imp;
                                    fclose(fi1);
                                }
                                if (SVPA == 4)
                                {
                                    FILE *fi1 = fopen ("/xlwork1/scott/full_model_nov_22/startup_files/pressure_128_mesh_sub4.dat", "r"); // Reading in pressure
                                    for (int j=0; j<=tmstps; j++)
                                    {
                                    fscanf(fi1,"%lf",&imp[j]);

                                    //imp[j] = (imp[j]/(rho*g*Lr/conv));
                                    Tube::imp_holder_4[j]=imp[j];
                                    }
                                    RD -> imp = imp;
                                    fclose(fi1);
                                }
                                if (SVPA == 5)
                                {
                                    FILE *fi1 = fopen ("/xlwork1/scott/full_model_nov_22/startup_files/pressure_128_mesh_sub5.dat", "r"); // Reading in pressure
                                    for (int j=0; j<=tmstps; j++)
                                    {
                                    fscanf(fi1,"%lf",&imp[j]);

                                    //imp[j] = (imp[j]/(rho*g*Lr/conv));
                                    Tube::imp_holder_5[j]=imp[j];
                                    }
                                    RD -> imp = imp;
                                    fclose(fi1);
                                }
                                if (SVPA == 6)
                                {
                                    FILE *fi1 = fopen ("/xlwork1/scott/full_model_nov_22/startup_files/pressure_128_mesh_sub6.dat", "r"); // Reading in pressure
                                    for (int j=0; j<=tmstps; j++)
                                    {
                                    fscanf(fi1,"%lf",&imp[j]);

                                    //imp[j] = (imp[j]/(rho*g*Lr/conv));
                                    Tube::imp_holder_6[j]=imp[j];
                                    }
                                    RD -> imp = imp;
                                    fclose(fi1);
                                }
                                if (SVPA == 7)
                                {
                                    FILE *fi1 = fopen ("/xlwork1/scott/full_model_nov_22/startup_files/pressure_128_mesh_sub7.dat", "r"); // Reading in pressure
                                    for (int j=0; j<=tmstps; j++)
                                    {
                                    fscanf(fi1,"%lf",&imp[j]);

                                    //imp[j] = (imp[j]/(rho*g*Lr/conv));
                                    Tube::imp_holder_7[j]=imp[j];
                                    }
                                    RD -> imp = imp;
                                    fclose(fi1);
                                }
                                if (SVPA == 8)
                                {
                                    FILE *fi1 = fopen ("/xlwork1/scott/full_model_nov_22/startup_files/pressure_128_mesh_sub8.dat", "r"); // Reading in pressure
                                    for (int j=0; j<=tmstps; j++)
                                    {
                                    fscanf(fi1,"%lf",&imp[j]);

                                    //imp[j] = (imp[j]/(rho*g*Lr/conv));
                                    Tube::imp_holder_8[j]=imp[j];
                                    }
                                    RD -> imp = imp;
                                    fclose(fi1);
                                }
                                if (SVPA == 9)
                                {
                                    FILE *fi1 = fopen ("/xlwork1/scott/full_model_nov_22/startup_files/pressure_128_mesh_sub9.dat", "r"); // Reading in pressure
                                    for (int j=0; j<=tmstps; j++)
                                    {
                                    fscanf(fi1,"%lf",&imp[j]);

                                    //imp[j] = (imp[j]/(rho*g*Lr/conv));
                                    Tube::imp_holder_9[j]=imp[j];
                                    }
                                    RD -> imp = imp;
                                    fclose(fi1);
                                }


                                if (SVPA == 0)
                                {
                                  for (int j =0; j<tmstps+1;j++)
                                  {
                                  imp[j] = 0.0;
                                  }
                                  RD -> imp = imp;
                                }










                                    impedance_driver_( & tmstps, & Period, & rho, & mu_pl, & rbot, & rmin, y11, y12, y21, y22, & Lr, & q, & g, & fa1, & fa2, & fa3, & fv1, & fv2, & fv3, & asym, & expo, & lrrA, & lrrV);



                                    for (int j=0; j<tmstps; j++)
                                    {
                                        pL[j] = p1[N];
                                    }







    }
}

// The destructor. When the tube-objects terminates, all arrays are
// deleted, in order to free the memory occupied by the object.

Tube :: ~Tube ()
{
  delete[] Anew;
  delete[] Qnew;
  delete[] Aold;
  delete[] Qold;
  delete[] Aprv;
  delete[] Qprv;
  delete[] Ah;
  delete[] Qh;
  delete[] y11;
  delete[] y12;
  delete[] y21;
  delete[] y22;
  delete[] pL;
  delete[] R1h;
  delete[] R2h;
  delete[] S1h;
  delete[] S2h;
  delete[] R1;
  delete[] R2;
  delete[] S1;
  delete[] S2;
  delete[] r0;
  delete[] r0h;
  delete[] dr0dx;
  delete[] dr0dxh;
  delete[] A0;
  delete[] A0h;
  delete[] fr;
  delete[] frh;
  delete[] dfrdr0;
  delete[] dfrdr0h;
  delete[] p1;
  delete[] p1h;
  delete[] dp1dr0;
  delete[] dp1dr0h;
}


//  The following functions prints P, Q, A, and F as functions of (x, t).
//  This is done in terms of the re-dimensionalized variables. In this case
//  the functions is plotted for a fixed time, but for all x along the
//  vessel in question. Since the doesn't have to be the first vessel in
//  the tree, it would have some offset from the heart. Which determines
//  the position for x. Therefore there are two arguments passed to this
//  function the time and the offset. Umar modified the print routines to
//  write important parameters in one file. Convenient files for WIA.

//void Tube :: printPxt (FILE *fd, double t, int offset)
void Tube::printPxt (FILE *fd, double t, int offset, int sim_time, int vessel_id){
	if (offset == 0) fprintf (fd, "\n");
	for (int i=0; i<N; i++){
		//fprintf (fd, "%d, %13.10f, %13.10f, %15.10f, %15.10f\n", vessel_id, t*Lr3/q, (i+offset)*h*Lr, (P(qLnb,i,Anew[i])+p0)*rho*g*Lr/conv,Qnew[i]*q);
		//  Time, Vessel Length, Fluid Pressure [mmHg],
		//  Flow rate [cm^3/sec], Cross-sectional Area [cm^2],
	}
}

double Tube::P(int qLnb, int i, double A) {
    double pold = fr[i] * (sqrt(A / A0[i]) - 1.0) + IPP[qLnb];
    // FLP iff A(TMP) Nov. 19.

    return pold;
}

double Tube::dPdA(int i, double A) {
    double pold = fr[i] / (sqrt(A0[i] * A) * 2.0); // Changed by JAM. FEB 2019.

    return pold;
}


double Tube::dPdx1(int i, double A)
// DERIVATIVE OF P WRT x IF WE USE A PROPERLY DEFINED AREA
{
    double pold = (dfrdr0[i] * (sqrt(A / A0[i]) - 1.0) - fr[i] * sqrt(M_PI * A) / A0[i]) *
        dr0dx[i];
    // Changed by JAM. FEB 2019.

    return pold;
}


double Tube::B(int i, double A) {
    //double pold = fr[i]*(sqrt(A0[i]*A)-A0[i])/Fr2;
    double pold = fr[i] * (sqrt(cu(A)) - sqrt(cu(A0[i]))) / (3.0 * sqrt(A0[i]) * Fr2);
    // correct assuming A(TMP)

    return pold;
}


double Tube::Bh(int i, double A) {
    int ip1 = i + 1;
    // double pold =  frh[ip1]*(sqrt(A0h[ip1]*A)-A0h[ip1])/Fr2;
    double pold = frh[ip1] * (sqrt(cu(A)) - sqrt(cu(A0h[ip1]))) /
        (3.0 * sqrt(A0h[ip1]) * Fr2);
    // correct assuming A(TMP)
    return pold;
}

double Tube::dBdx1(int i, double A) {
    double dfr = dfrdr0[i];
    double pold = ((dfr / sqrt(A0[i]) - sqrt(M_PI) * fr[i] / A0[i]) * (sqrt(cu(A)) -
            sqrt(cu(A0[i]))) / 3.0 +
        fr[i] *
        (M_PI * r0[i] * (sqrt(cu(A / A0[i])) - 1.0) - A *
            (sqrt(A / A0[i]) - 1.0) * dfr / fr[i])
    ) * dr0dx[i] / Fr2;

    // NEW TUBE LAW. JAM. JAN 2019. // LINE BREAKS. JAM. APR 2020.

    return pold;
}

double Tube::dBdx1h(int i, double A) {
    int ip1 = i + 1;

    double dfr = dfrdr0h[ip1];

    double pold = (
        (dfr / sqrt(A0h[ip1]) - sqrt(M_PI) * frh[ip1] / A0h[ip1]) *
        (sqrt(cu(A)) - sqrt(cu(A0h[ip1]))) / 3.0 +
        frh[ip1] *
        (M_PI * r0h[ip1] * (sqrt(cu(A / A0h[ip1])) - 1.0) -
            A * (sqrt(A / A0h[ip1]) - 1.0) * dfr / frh[ip1])
    ) * dr0dxh[ip1] / Fr2;

    // NEW TUBE LAW. JAM. JAN 2019.

    return pold;
}



double Tube::dBdAh(int i, double A) {
    int ip1 = i + 1;
    //  double pold = 0.5*frh[ip1]*sqrt(A0h[ip1]/A)/Fr2;
    double pold = frh[ip1] * sqrt(A / A0h[ip1]) / (2.0 * Fr2);
    // NEW TUBE LAW BY JAM. FEB 2019.

    return pold;
}

double Tube::d2BdAdxh(int i, double A) {
    int ip1 = i + 1;
    double dfr = dfrdr0h[ip1];
    //  double pold = (-dfr+1.0/sqrt(A)*(sqrt(M_PI)*frh[ip1]+
    //           sqrt(A0h[ip1])*dfr))*dr0dxh[ip1]/Fr2;
    double pold = (dfr +
        sqrt(A * M_PI) * frh[ip1] / A0h[ip1] -
        sqrt(A / A0h[ip1]) * dfr) * dr0dxh[ip1] / Fr2;

    // NEW TUBE LAW. JAM. FEB 2019. // LINE BREAKS. JAM. APR 2020.
    // CHECK AGAIN
    return pold;
}

//  When determining or checking the step-size (k) the CFL-condition is
//  applied. This is determined according to the result reached from the
//  analysis made using the method of characteristics
//  (See IMFUFATEKST no 297). In this function the minimal step-size
//  fulfilling this condition for this tube is returned.

double Tube::CFL() // The CFL-condition
{
    double minimum = 64000000.0;
    for (int i = 0; i <= N; i++) {
        double c_tmp = c(i, Anew[i]);
        double Vnew = Qnew[i] / Anew[i];
        double temp = min(h / fabs(Vnew - c_tmp),
            h / fabs(Vnew + c_tmp));
        if (temp < minimum) minimum = temp;
    }
    return (minimum);
}

// When taking a Lax-Wendroff step, the flux of the system must be
// determined. This is evaluated at i + j/2, and the prediction is given
// as described in IMFUFATEKST no 297 and D2.1-4. The integer k determines
// whether we deal with the first or the second component of the vector.

double Tube::Rvec(int k, int i, int j, double Q, double A) {
    if (k == 1) return (Q);
    else
    if (k == 2) return (sq(Q) / A + ((j == 0) ? B(i, A) : Bh(i, A)));
    else error("arteries.cxx", "Call of non-existing vector-component of R");
    return (0);
}

//  Similarly the right hand side of the system of equations must be
//  determined at i + j/2. Also in this case the function is given as
//  stated in the mathematical model, and also in this case k states the
//  needed component of the vector.

double Tube::Svec(int k, int i, int j, double Q, double A) {
    if (k == 1) return (0.0);
    else
    if (k == 2) return (F(Q, A) + ((j == 0) ? dBdx1(i, A) : dBdx1h(i, A)));
    else error("arteries.cxx", "Call of non-existing vector-component of S");
    return (0);
}

//  The solutions of Anew and Qnew are found for all interior points of the
//  vessel at (t+k), where k is the length of the current time-step. This
//  function saves the results in the arrays Anew and Qnew, and the
//  function  is made according to Lax-Wendroff's method as described in
//  IMFUFATEKST  no 297 and D2.1-4.

void Tube::step(double k, int ii, double t, double t0) {
    double theta = k / h; // Theta is determined.
    double gamma = 0.5 * k; // Gamma is determined.

    if (t0 != 0 && t == t0) {
        if (ii == 0) {printf("restarting\n");}

        // load in all quants usually found here from data files // only in restartables
        // char fAname[20]; sprintf(fAname, "./%s/sv_1_%i.2d", file_name, ii); FILE * fA = fopen(fAname, "r");
        // for (int i = 0; i <= N; i++) {fscanf(fA, "%lf,%lf,%lf,%lf,%lf,%lf",
        //     & Aold[i], & Qold[i], & R1[i], & R2[i], & S1[i], & S2[i]);}
        // fclose(fA);

        // sprintf(fAname, "./%s/sv_2_%i.2d", file_name, ii); fA = fopen(fAname, "r");
        // for (int i = 0; i < N; i++) {
        //     fscanf(fA, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
        //         &Anew[i], & Qnew[i], & Ah[i], & Qh[i], &R1h[i], & R2h[i], & S1h[i], & S2h[i]);
        // }
        // fclose(fA);
    } else { //if not restarting
        if (t == t0 + k && ii == 0) {printf("Second time step\n");}

        for (int i = 0; i <= N; i++) {
            // Remember the values at this time level.
            Qold[i] = Qnew[i];
            Aold[i] = Anew[i];

            //  Anew and Qnew are predicted at the new time level (t+k).
            R1[i] = Rvec(1, i, 0, Qold[i], Aold[i]);
            R2[i] = Rvec(2, i, 0, Qold[i], Aold[i]);
            S1[i] = Svec(1, i, 0, Qold[i], Aold[i]);
            S2[i] = Svec(2, i, 0, Qold[i], Aold[i]);
        }

        for (int i = 0; i < N; i++) {
            Ah[i] = 0.5 * (Aold[i + 1] + Aold[i]) - 0.5 * theta * (R1[i + 1] - R1[i]) + 0.5 * gamma * (S1[i + 1] + S1[i]);
            Qh[i] = 0.5 * (Qold[i + 1] + Qold[i]) - 0.5 * theta * (R2[i + 1] - R2[i]) + 0.5 * gamma * (S2[i + 1] + S2[i]);
            R1h[i] = Rvec(1, i, 1, Qh[i], Ah[i]);
            R2h[i] = Rvec(2, i, 1, Qh[i], Ah[i]);
            S1h[i] = Svec(1, i, 1, Qh[i], Ah[i]);
            S2h[i] = Svec(2, i, 1, Qh[i], Ah[i]);
        }

        for (int i = 1; i < N; i++) {
            Anew[i] = Aold[i] - theta * (R1h[i] - R1h[i - 1]) + gamma * (S1h[i] + S1h[i - 1]);
            Qnew[i] = Qold[i] - theta * (R2h[i] - R2h[i - 1]) + gamma * (S2h[i] + S2h[i - 1]);
        }
    }
}


//  ==========================  INPUT CONDITION  ==========================
//  The left boundary (x=0) uses this function to model an inflow into
//  the system. The actual parameter given to the function is the model
//  time. As stated in the mathematical model the constants of the function
//  are chosen in order to ensure a certain CO (specified in main.h). Hence
//  we have the specified value of b. Further the period (dimension-less)
//  is assumed to be Period.

//  ---------------------------  FOR PRESSURE  ----------------------------
double Tube::Ps_init(double t, double k, double Period) // FLUID PRESSURE
{
  if (t < Period) return (Ps[int(t / k)]); else // FLUID PRESSURE
  if (t >= Period) return (Ps_init((t - Period), k, Period)); // FLUID PRESSURE
  else return (0);
}
//  -----------------------------------------------------------------------

double Tube :: IPP_init (double t, double k, double Period)
{
  if (t <= Period) return (IPP[int(t/k)]); else
  if (t >  Period) return (IPP_init((t-Period),k,Period));
  else return (0);
}

// see end of arteries.h

//  Update of the left boundary at time t. This function uses Q0 to
//  determine the flow rate at the next time-step. From this the value of A
//  is predicted using Lax-Wendroff's numerical scheme. This function is
//  only relevant when the tube is an inlet vessel.


//  ---------------------------  FOR PRESSURE  ----------------------------
//  uncomment for pressure input / comment for flow

void Tube::bound_left(int qLnb, double t, double k, double Period) {


    double Pnew = 0.0; // FLUID PRESSURE
    double P05 = 0.0; // FLUID PRESSURE

    if (Tube::initial_or_integrated == 0)
    {
    Pnew   = Ps_init(t, k, Period); // FLUID PRESSURE
    P05    = Ps_init(t - k, k, Period); // FLUID PRESSURE
    }
    else if (Tube::initial_or_integrated == 1)
    {
    Pnew   = (Tube::transfer_aortic) * conv / rho / g / Lr; // FLUID PRESSURE
    P05    = (Tube::transfer_aortic_previous) * conv / rho / g / Lr; // FLUID PRESSURE
    }




    if (int(t / k) < 0) {
        printf("t/k negative in bound_right\n");
    }

    Anew[0] = A0[0] * (sq(((Pnew - IPP[qLnb]) / fr[0]) + 1.0));

    double A05;
    if (qLnb == 0) {
        A05 = A0[0] * (sq(((P05 - IPP[qLnb]) / fr[0]) + 1.0));
    } else if (qLnb != 0) {
        A05 = A0[0] * (sq(((P05 - 0.5 * (IPP[qLnb] + IPP[qLnb - 1])) / fr[0]) + 1.0));
    } //  if statements ensure that qLnb in the range [0, tmstps]

    double Ahm05 = 2.0 * A05 - Ah[0];
    double Qhm05 = R1h[0] + h * (Anew[0] - Aold[0]) / k;
    double R2hm05 = Rvec(2, -1, 1, Qhm05, Ahm05);
    double S2hm05 = Svec(2, -1, 1, Qhm05, Ahm05);
    Qnew[0] = Qold[0] - k * (R2h[0] - R2hm05) / h + k * (S2h[0] + S2hm05) / 2.0;
}
//  -----------------------------------------------------------------------
//  =======================================================================



//  The value at the right boundary at time t is predicted. NB: This should
//  only be used with terminal vessels, i.e. for vessels that don't
//  bifurcate into further branches. In that situation the bifurcation
//  boundary function should be called instead. Again the procedure
//  specified is given according to the mathematical theory presented in
//  IMFUFATEKST no 297 and D2.1-4.

double Tube::c(int i, double A) // The wave speed through aorta.
{
    //double cnst =  0.5*fr[i]*sqrt(A0[i]/A)/Fr2;
    double cnst = 0.5 * (fr[i] / sqrt(A0[i] * A)) / Fr2;
    return sqrt(cnst);
}



// The matching boundary. // All P quantites should FLP

// The matching boundary. // All P quantites should FLP

// The matching boundary. // All P quantites should FLP
void Tube::bound_match(int qLnb, double t, double k, double theta, double gamma, double t0) {

        double PA, PV, PAh, PVh; // FLP. Defined ln. 77.

        // int j = 1, ok = false, ntrial = 50000, qLnb_1 = qLnb - 1;
        int j = 1, ok = false, ntrial = 50000, qLnb_1 = qLnb + 1;

        // Make sure that qLnb_1 runs in the interval [0:tmstps-1].
        if (qLnb_1 == (int) tmstps) {qLnb_1 = 0;}

        //  In order to make a relation between P(x_L, t+dt) and Q(x_L, t+dt), and P(x_L, t+dt/2) and Q(x_L,t+dt/2) we need to extract the term involving y[0] (see mathematical derivation). This term corresponds to the peripheral. The remaining terms in the convolution present at the boundary, see mathematical derivation.

        double paterms = 0.0, pvterms = 0.0, paterms_half = 0.0, pvterms_half = 0.0;
        int M = RD -> N;
        double RDtheta = k / RD -> h;




        if (t > Period) {
            //convolution integral
            for (int m = 1; m < tmstps; m++) {
                // Corrected time-stepping from NAH.
                int pindex = (qLnb_1 + tmstps - m) % tmstps;
                int pindex1 = pindex - 1;
                if (pindex == 0) {
                    pindex1 = tmstps - 1;
                }

                double small_vessel_curr;
                double small_vessel_prev;

                if (SVPA == 0)
                {
                small_vessel_curr = 0.0;
                small_vessel_prev = 0.0;
                }


                if (SVPA == 1)
                {
                small_vessel_curr = Tube::imp_holder_1[pindex]/(rho*g*Lr/conv);
                small_vessel_prev = Tube::imp_holder_1[pindex1]/(rho*g*Lr/conv);
                }

                if (SVPA == 2)
                {
                small_vessel_curr = Tube::imp_holder_2[pindex]/(rho*g*Lr/conv);
                small_vessel_prev = Tube::imp_holder_2[pindex1]/(rho*g*Lr/conv);
                }

                if (SVPA == 3)
                {
                small_vessel_curr = Tube::imp_holder_3[pindex]/(rho*g*Lr/conv);
                small_vessel_prev = Tube::imp_holder_3[pindex1]/(rho*g*Lr/conv);
                }

                if (SVPA == 4)
                {
                small_vessel_curr = Tube::imp_holder_4[pindex]/(rho*g*Lr/conv);
                small_vessel_prev = Tube::imp_holder_4[pindex1]/(rho*g*Lr/conv);
                }

                if (SVPA == 5)
                {
                small_vessel_curr = Tube::imp_holder_5[pindex]/(rho*g*Lr/conv);
                small_vessel_prev = Tube::imp_holder_5[pindex1]/(rho*g*Lr/conv);
                }

                if (SVPA == 6)
                {
                small_vessel_curr = Tube::imp_holder_6[pindex]/(rho*g*Lr/conv);
                small_vessel_prev = Tube::imp_holder_6[pindex1]/(rho*g*Lr/conv);
                }

                if (SVPA == 7)
                {
                small_vessel_curr = Tube::imp_holder_7[pindex]/(rho*g*Lr/conv);
                small_vessel_prev = Tube::imp_holder_7[pindex1]/(rho*g*Lr/conv);
                }

                if (SVPA == 8)
                {
                small_vessel_curr = Tube::imp_holder_8[pindex]/(rho*g*Lr/conv);
                small_vessel_prev = Tube::imp_holder_8[pindex1]/(rho*g*Lr/conv);
                }

                if (SVPA == 9)
                {
                small_vessel_curr = Tube::imp_holder_9[pindex]/(rho*g*Lr/conv);
                small_vessel_prev = Tube::imp_holder_9[pindex1]/(rho*g*Lr/conv);
                }

//                double small_vessel_curr = imp[pindex];
//                double small_vessel_prev = imp[pindex1];

                paterms = paterms +
                    (pL[pindex] - small_vessel_curr) * y11[m] +
                    ((RD -> pL[pindex]) - small_vessel_curr) * y12[m];

                pvterms = pvterms +
                    (pL[pindex] - small_vessel_curr) * y21[m] +
                    ((RD -> pL[pindex]) - small_vessel_curr) * y22[m]; //YJ

                paterms_half = paterms_half +
                    0.5 * (pL[pindex] + pL[pindex1]) * y11[m] +
                    0.5 * ((RD -> pL[pindex]) + (RD -> pL[pindex1])) * y12[m] -
                    0.5 * (small_vessel_curr + small_vessel_prev) * y11[m] -
                    0.5 * (small_vessel_curr + small_vessel_prev) * y12[m];

                pvterms_half = pvterms_half +
                    0.5 * (pL[pindex] + pL[pindex1]) * y21[m] +
                    0.5 * ((RD -> pL[pindex]) + (RD -> pL[pindex1])) * y22[m] -
                    0.5 * (small_vessel_curr + small_vessel_prev) * y21[m]-
                    0.5 * (small_vessel_curr + small_vessel_prev) * y22[m];
            }
            paterms = k * paterms;
            pvterms = k * pvterms;
            paterms_half = k * paterms_half;
            pvterms_half = k * pvterms_half;
        }

        double k1 = Aold[N] + theta * R1h[N - 1];
        double k2 = Qold[N] + theta * R2h[N - 1] + gamma * S2h[N - 1];
        double k3 = paterms;
        double k4 = k * y11[0];
        double k5 = k * y12[0];
        double k6 = RD -> Aold[M] + RDtheta * (RD -> R1h[M - 1]);
        double k7 = RD -> Qold[M] + RDtheta * (RD -> R2h[M - 1]) + gamma * (RD -> S2h[M - 1]);
        double k8 = pvterms;
        double k9 = k * y21[0];
        double k10 = k * y22[0];
        double k11 = paterms_half - 0.5 * Qh[N - 1];
        double k12 = Ah[N - 1];
        double k13 = RD -> Ah[M - 1];
        double k14 = pvterms_half - 0.5 * (RD -> Qh[M - 1]);
        double k15 = k4; // JAM   //715;//
        double k16 = k5; // JAM   //716;//
        double k17 = k9; // JAM   //717;//
        double k18 = k10; // JAM   //718;//
        //  Unknowns declared, and initial guesses applied

        double xb[8];

        //  xb is a vector used to find the residuals for the boundary match
        if (t <= Period) {
            xb[0] = Ah[N - 1];
            xb[1] = Qh[N - 1];
            xb[2] = Aold[N];
            xb[3] = Qold[N];
            xb[4] = RD -> Ah[M - 1];
            xb[5] = RD -> Qh[M - 1];
            xb[6] = RD -> Aold[M];
            xb[7] = RD -> Qold[M];
        }

        if (t > Period) {
            xb[0] = Aold[N];
            xb[1] = Qold[N];
            xb[2] = Ah05; //Aold[N];
            xb[3] = Qh05; //Qold[N];
            xb[4] = RD -> Aold[M];
            xb[5] = RD -> Qold[M];
            xb[6] = RD -> Ah05; //RD->Aold[M];
            xb[7] = RD -> Qh05; //RD->Qold[M];
        }

        while (j <= ntrial && ok == false) {
            double fvec[8];
            PA = P(qLnb, N, xb[0]); // FLP iff P(i, A[i]) FLP.
            PV = RD -> P(qLnb, M, xb[4]); // FLP iff P(i, A[i]) FLP.
            PAh = P(qLnb, N, 0.5 * (k12 + xb[2])); // FLP iff P(i, A[i]) FLP.
            PVh = RD -> P(qLnb, M, 0.5 * (k13 + xb[6])); // FLP iff P(i, A[i]) FLP.

            double IMP=0.0;


            if (SVPA == 0)
            {
            IMP = 0.0;
            }


            if (SVPA == 1)
            {
            IMP = Tube::imp_holder_1[qLnb]/(rho*g*Lr/conv);
            }

            if (SVPA == 2)
            {
            IMP = Tube::imp_holder_2[qLnb]/(rho*g*Lr/conv);
            }

            if (SVPA == 3)
            {
            IMP = Tube::imp_holder_3[qLnb]/(rho*g*Lr/conv);
            }

            if (SVPA == 4)
            {
            IMP = Tube::imp_holder_4[qLnb]/(rho*g*Lr/conv);
            }

            if (SVPA == 5)
            {
            IMP = Tube::imp_holder_5[qLnb]/(rho*g*Lr/conv);
            }

            if (SVPA == 6)
            {
            IMP = Tube::imp_holder_6[qLnb]/(rho*g*Lr/conv);
            }

            if (SVPA == 7)
            {
            IMP = Tube::imp_holder_7[qLnb]/(rho*g*Lr/conv);
            }

            if (SVPA == 8)
            {
            IMP = Tube::imp_holder_8[qLnb]/(rho*g*Lr/conv);
            }

            if (SVPA == 9)
            {
            IMP = Tube::imp_holder_9[qLnb]/(rho*g*Lr/conv);
            }




//            IMP = imp[qLnb];

            //  Vector of residuals
            fvec[0] = k1 - xb[0] - theta * xb[3];

            fvec[1] = k2 - xb[1] - theta * (sq(xb[3]) / xb[2] + Bh(N, xb[2])) +
                gamma * (F(xb[3], xb[2]) + dBdx1h(N, xb[2]));

            fvec[2] = k3 - xb[1] + k4 * PA + k5 * PV - k15 * IMP - k16 * IMP;

            fvec[3] = k6 - xb[4] - RDtheta * xb[7];

            fvec[4] = k7 - xb[5] - RDtheta * (sq(xb[7]) / xb[6] + RD -> Bh(M, xb[6])) +
                gamma * (RD -> F(xb[7], xb[6]) + RD -> dBdx1h(M, xb[6]));

            fvec[5] = k8 - xb[5] + k9 * PA + k10 * PV - k17 * IMP - k18 * IMP;

            fvec[6] = k11 - xb[3] / 2.0 + k4 * PAh + k5 * PVh - k15 * IMP - k16 * IMP;

            fvec[7] = k14 - xb[7] / 2.0 + k9 * PAh + k10 * PVh - k17 * IMP - k18 * IMP;

            //  The Jacobian
            for (int row = 0; row < 8; row++) {
                for (int col = 0; col < 8; col++) {
                    fj[row][col] = 0.0;
                }
            }

            fj[0][0] = -1.0;
            fj[0][3] = -theta;

            fj[1][1] = -1.0;
            fj[1][2] = theta * (sq(xb[3] / xb[2]) - dBdAh(N, xb[2])) +
                gamma * (dFdA(xb[3], xb[2]) + d2BdAdxh(N, xb[2]));
            fj[1][3] = -2.0 * theta * xb[3] / xb[2] + gamma * dFdQ(xb[2]);

            fj[2][0] = k4 * dPdA(N, xb[0]);
            fj[2][1] = -1.0;
            fj[2][4] = k5 * RD -> dPdA(M, xb[4]);

            fj[3][4] = -1.0;
            fj[3][7] = -RDtheta;

            fj[4][5] = -1.0;
            fj[4][6] = RDtheta * (sq(xb[7] / xb[6]) - RD -> dBdAh(M, xb[6])) +
                gamma * (dFdA(xb[7], xb[6]) + RD -> d2BdAdxh(M, xb[6]));
            fj[4][7] = -2.0 * RDtheta * xb[7] / xb[6] + gamma * dFdQ(xb[6]);

            fj[5][0] = k9 * dPdA(N, xb[0]);
            fj[5][4] = k10 * RD -> dPdA(M, xb[4]);
            fj[5][5] = -1.0;

            fj[6][2] = k4 * dPdA(N, (k12 + xb[2]) / 2.0);
            fj[6][3] = -0.5;
            fj[6][6] = k5 * RD -> dPdA(M, (k13 + xb[6]) / 2.0);

            fj[7][2] = k9 * dPdA(N, (k12 + xb[2]) / 2.0);
            fj[7][6] = k10 * RD -> dPdA(M, (k13 + xb[6]) / 2.0);
            fj[7][7] = -0.5;

            // apply the zero finding function
            double tol = 1.0e-8;
            int ch = zero(xb, 8, tol, tol, fvec, fj);
            if (ch == 1) ok = true;
            j = j + 1;
        }

        Anew[N] = xb[0];
        Qnew[N] = xb[1];
        Ah05 = xb[2];
        Qh05 = xb[3];
        RD -> Anew[M] = xb[4];
        RD -> Qnew[M] = xb[5];
        RD -> Ah05 = xb[6];
        RD -> Qh05 = xb[7];
        pL[qLnb_1] = P(qLnb_1, N, Anew[N]); // FLP as P(N, Anew[N]) is FLP
        RD -> pL[qLnb_1] = RD -> P(qLnb, M, xb[4]); // FLP as P(N, Anew[N]) is FLP


        if (j >= ntrial) {

            error("arteries.C", "Root not found in the matching");
            cout << "t = "<< t * Lr3 / q << "t/(k*tmstps) = " << t / (k * tmstps) << "\n";
            printf("AM=%10.15f, QM=%10.15f \n", xb[0], xb[1]);
            printf("AL=%10.15f, QL=%10.15f \n", xb[4], xb[5]);
            printf("rbot=%10.15f \n", rbot);
  //          printf("Arteries[%i] \n", ID);
            exit(1);
        }

}


//  The value at the bifurcation point at time t is predicted. NB: This
//  should only be done for tubes that do bifurcate into further branches.
//  If this is not the case we have a terminal vessel and bound_right
//  should be called instead. The procedure operates according to the
//  specifications in the mathematical model as a link between this tube
//  and its daughters. Therefore there will be three tubes involved in
//  this function. One problem is however, that the rather complicated
//  system of equations does not converge for all choices of parameters
//  (the peripheral resistance, the top radius, and the bottom radius).



void Tube::call_junc(int qLnb, double theta, double gamma, Tube * Arteries[], int parent) {

    // cout << "junc called\n";

    Tube * D1 = Arteries[parent] -> LD;
    Tube * D2 = Arteries[parent] -> MD;
    Tube * D3 = Arteries[parent] -> RD;

    // cout << "daughters pointed to\n";

    int j = 1;
    int ok = false;
    const int ntrial = 4000;
    //    double tol = 1.0e-8;
    double tol = 1.0e-6;
    //  check how many daughter vessels are given for each junction to match, and call the appropriate function from junction.c

    if (D1 != 0 && D2 == 0 && D3 == 0) {
        // cout << "mono in\n";
        monofurcation(qLnb, theta, gamma, Arteries, parent, ntrial, tol, j, ok);
        // cout << "mono out\n";
    } else if (D1 != 0 && D2 == 0 && D3 != 0) {
        // cout << "bif in\n";
        bifurcation(qLnb, theta, gamma, Arteries, parent, ntrial, tol, j, ok);
        // cout << "bif out\n";
    } else if (D1 != 0 && D2 != 0 && D3 != 0) {
        // cout << "daughters pointed to\n";
        trifurcation(qLnb, theta, gamma, Arteries, parent, ntrial, tol, j, ok);
    }
}




//  The right boundary (x=L) uses this function to model an inflow into the
//  system. The actual parameter given to the function is the model time.
//  As stated in the mathematical model the constants of the function are
//  chosen in order to ensure a certain CO (specified in main.h). Hence we
//  have the specified value of b. Further the period (dimension-less) is
//  assumed to be Period.

double Tube::PL_init(double t, double k, double Period) {
    if (t <= Period) return (Pout[int(t / k)]);
    else
    if (t > Period) return (PL_init((t - Period), k, Period));
    else return (0);
}

//  Update of the right boundary at time t. This function uses PL to
//  determine A at the next time-step. From this the value of Q is
//  predicted using Lax-Wendroff's numerical scheme. This function is only
//  determine A at the next time-step. From this the value of Q is
//  predicted using Lax-Wendroff's numerical scheme. This function is only
//  determine A at the next time-step. From this the value of Q is
//  predicted using Lax-Wendroff's numerical scheme. This function is only
//  relevant when the tube is an outlet vessel.

void Tube::bound_right(int qLnb, double t, double k, double Period) {
    double Pnew = PL_init(t, k, Period);
    double P05 = PL_init(t - k, k, Period);

    if (int(t / k) < 0)
        printf("t/k negative in bound_right\n");

    Anew[0] = A0[0] * (sq(((Pnew - IPP[qLnb]) / fr[0]) + 1.0));

    //    double A05 = A0[0] * (sq(((P05 - IPP[qLnb]) / fr[0]) + 1.0));
    double A05;
    if (qLnb == 0) {
        A05 = A0[0] * (sq(((P05 - 0.5 * (IPP[qLnb] + IPP[qLnb])) / fr[0]) + 1.0));
    } else if (qLnb != 0) {
        A05 = A0[0] * (sq(((P05 - 0.5 * (IPP[qLnb] + IPP[qLnb - 1])) / fr[0]) + 1.0));
    }

    double Ahm05 = 2.0 * A05 - Ah[0];
    double Qhm05 = R1h[0] + h * (Anew[0] - Aold[0]) / k;
    double R2hm05 = Rvec(2, -1, 1, Qhm05, Ahm05);
    double S2hm05 = Svec(2, -1, 1, Qhm05, Ahm05);

    Qnew[0] = Qold[0] - k * (R2h[0] - R2hm05) / h + k * (S2h[0] + S2hm05) / 2.0;

}



void looper(int qLnb, Tube * Arteries[]) {
printf("in looper\n");
    double t0=0.0,tstart,tend;

    const char * nameVessel; FILE * fVessel; char v[20];

    double k = Tube::onePeriodDt,Period = Tube::Period, Deltat = Tube::Deltat;
	int tmstps = Tube::tmstps, plts = Tube::plts;

    int ves = 10, pos = Arteries[ves] -> N; // from MJC. Loop until converged in the given vessel
    int period_counter = 1 + t0 / Period; // Count the number of periods you have solved for
    double norm_sol = 1e+6, threshold = 0.01, sse, sol_p1[plts], sol_p2[plts], pct, pctmax;

    int Nmax = 0;
    for (int i = 0; i < nbrves; i++) {Nmax = int(std::fmax(Nmax, Arteries[i] -> N));}
    tstart = t0; tend = tstart + Deltat;
    int sol_ID = 0;

    while (tend <= period_counter * Period) {
        solver(Arteries, tstart, tend, k,t0);

        sol_p1[sol_ID] = Arteries[ves] -> P(qLnb,0, Arteries[ves] -> Anew[pos]);
        tstart = tend;
        tend = tend + Deltat; // The current ending time is increased by Deltat.
        sol_ID++;
    }
printf("end first solve\n");
std::cout << "end first solve\n";



sprintf(v, "V%i.2d", nbrves); nameVessel = v; fVessel = fopen(nameVessel, "w");

while (norm_sol >= threshold) {
        sol_ID = 0; sse = 0; pctmax = 0; period_counter++;

        while (tend <= period_counter * Period) {
            solver(Arteries, tstart, tend, k,t0);

        for (int j = 0; j < nbrves; j++) {Arteries[j] -> printPxt(fVessel, tend, 0, fmod(tstart / k, tmstps), j);}



            sol_p2[sol_ID] = Arteries[ves] -> P(qLnb,0, Arteries[ves] -> Anew[pos]); // for printing
            sse = sse + sq(sol_p1[sol_ID] - sol_p2[sol_ID]);

            pct = 100 * (sol_p1[sol_ID] - sol_p2[sol_ID]) / sol_p2[sol_ID];

            if (abs(pct) > abs(pctmax)) {pctmax = pct;}

            tstart = tend;
            tend = tend + Deltat; // move to the next time interval of length Deltat.
            sol_ID++;
        }
        norm_sol = abs(pctmax);
        memcpy(sol_p1, sol_p2, sizeof(sol_p2));

   printf("NORM_SOL: %f\n",sse);
            printf("pct max: %f\n\n", pctmax);
      std::cout<<"NORM_SOL: " <<sse << "\npct max: "<< pctmax << "\n\n";
    }

    std::cout<<"num_cylces: "<< int(period_counter - (t0 / Period)) <<"\n";



//    sprintf(v, "V%i.2d", nbrves); nameVessel = v; fVessel = fopen(nameVessel, "w");

    while (tend <= (period_counter + 1) * Period) {
        if (tstart == (period_counter) * Period) {
            std::cout<< "saving\n";
        }
        solver(Arteries, tstart, tend, k,t0);
        for (int j = 0; j < nbrves; j++) {Arteries[j] -> printPxt(fVessel, tend, 0, fmod(tstart / k, tmstps), j);}
        tstart = tend;
        tend = tend + Deltat; // move to the next time interval of length Deltat.
    }

fclose(fVessel);
}



//  Solves the non-linear PDE's (momentum and continuity eqn's.) from
//  t = tstart to t= tend. This function checks the maximal possible size
//  of the next time-step, reduces it to make sure that we don't walk to
//  far, and takes the next step. This is done by executing the step
//  routine, then updating the left boundary and finally updating
//  bifurcation points and the right boundaries. This is carried out as
//  long as the time hasn't passed the desired ending time (tend) which is
//  passed to the function as a parameter.

 void solver (Tube *Arteries[], double tstart, double tend, double k, double t0)
 {
   // The following definitions only used with variable time-stepping
   double t = tstart;
   int qLnb = (int) fmod(t / k, Tube::tmstps);


   // As long as we haven't passed the desired ending time do:
   while (t < tend) {
       // Check that the step we take is valid. If this error occurs when forcing a constant step-size the program must be terminated.
       if (t + k > tend) {
           double kold = k;
           k = tend - t;
           printf("ERROR (arteries.C): Step-size changed:\n");
           printf("t+k=%10.15f, tend=%10.15f k=%10.15f kold=%10.15f\n",
               t + kold, tend, k, kold);
       }

       // Check that the CFL-condition applies.
       for (int i = 0; i < nbrves; i++) {
           if (k > Arteries[i] -> CFL()) {
               cout << "Arteries[" << i << "]" << endl;
               cout << "time = " << (t * Tube::Lr3 / Tube::q) << "s" << endl;
               error("arteries.C", "Time-step size too large\n CFL-condition violated\n Now exiting.\n");
               //                exit(2); // suppress exit
           }
       }

       // solve for interior points, by calling step.
       // for (int i = 0; i < nbrves; i++) {Arteries[i] -> step(k, t0);}
       for (int i = 0; i < nbrves; i++) {
           Arteries[i] -> step(k, i, t, t0);
       }

       // Update left and right boundaries, and the bifurcation points.
       for (int i = 0; i < nbrves; i++) {
           if (Arteries[i] -> init == 1) {
               Arteries[i] -> bound_left(qLnb, t + k, k, Tube::Period);
           };
           if (Arteries[i] -> init == 2) {
               Arteries[i] -> bound_right(qLnb, t + k, k, Tube::Period);
           };
           if (Arteries[i] -> init == 3)
           {
               double theta = k / Arteries[i] -> h, gamma = k / 2.0;
               Arteries[i] -> bound_match(qLnb, t, k, theta, gamma, t0);
           }

           if (Arteries[i] -> rm == 0) {
               double theta = k / Arteries[i] -> h, gamma = k / 2.0;
               Arteries[i] -> call_junc(qLnb, theta, gamma, Arteries, i);
           }
       }
       t = t + k; // Uate the time and position within one period.
       qLnb = (qLnb + 1) % Tube::tmstps;
   }

}
