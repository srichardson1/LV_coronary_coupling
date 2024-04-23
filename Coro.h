/***************************************************************************/
/*                                                                         */
/*  Program: sor06.h                                                       */
/*  Version: 2.0                                                           */
/*  By: Mette Olufsen, Math-Tech                                           */
/*  Date: 14. Jan. 1997                                                    */
/*                                                                         */
/*  This header file defines the global parameters                         */
/*                                                                         */
/***************************************************************************/

// $Id: sor06.h,v 1.7 2005/07/07 22:19:51 heine Exp $
// Last updated on October 23, 2014 by M. Umar Qureshi
// coronary version

#ifndef _Coro_H
#define _Coro_H

#include <cmath>
#include "tools.h"

int    nbrves, N_aorta;               // Number of vessels in the tree.

//int tmstps = 131072;          // The number of timesteps per period.
//int tmstps = 32768;
//int	   plts   = 1024;                 // Number of plots per period.

	   // For pressure driven use :   const char*  CO_filename = "p0_in4096.dat";
	   //     For flow driven use :   const char*  CO_filename = "gvPA_4096.dat";


  const char*  CO_filename    = "corPress.dat";   // Input flow file at the heart. // Comment for pressure driven / Uncomment for flow driven
//const char*  PL_filename    = "zeros_4096.dat";        // Outflow (Static pressure) file.
  const char*  small_filename = "externalPressure32768.dat"; //extSin.dat"; // "Pe_16384.dat";
  const char*  large_filename = "externalPressure32768.dat"; //extSin.dat"; // "Pe_16384.dat";
  // x = linspace(0, 0.7, 32769); extSin = 0.5*(sin(x*2*pi/0.7+0.646)+1);

//double conv   = 1333.220,             // Conversion from mmHg to SI-units.
//       rho    = 1.055,                // Density of blood [g/cm^3].
//       mu1     = 0.049,                // Viscosity of blood [g/cm/s].
//       mu_pl  = 0.049,                   // Viscosity of blood [g/cm/s].
//       nu     = mu1 / rho,
//       Tper   = 0.7,                 // The period of one heart beat [s].

   //  rm     = 0.005,             // The minimum radius for the structured tree [cm].

   // Stiffness parameters for the large vessels.

   //  f1     = 0, //1.99925e07,
   //  f2     = 1,// -25.5267,
   //  f3     = 260000,//8.65e05    //which is 0.3*k3 (k3 = 8.65e05 in Olufsen's et al (1999))

   // Stiffness parameters for the small arteries.
   //  fa1    =  0.0, // 1.99925e07,
   //  fa2    =  1.0, // -25.5267,
   //  fa3    =  50000,// 465251,

   // Stiffness parameters for the small veins.
   //  fv1    =  0.0,  // 1.99925e07,
   //  fv2    =  1.0,  // -25.5267,
   //  fv3    =  50000,//165000,//50000.0, // 465251,

   // To simulate the PAH and CTEPH and HLD in Qureshi et al (2014)
   // Used f3 = 0.48*k3 for HLD plots or 340000 if starting from xi=2.4 (and for PVH)
   // For PAH use 0.45*k3 or 305000 for up to 50% increase in fa3,
   // fv3 and 336000 for up to 75% increase in fa3 and fv3

   //  asym   =  0.41,                // The asymmetry ratio of the structured tree.
   //  expo   =	 2.76,              // Exponent in radius relation. (Decrease up to 2.3 to simulate PH associated with HLD)

   //  lrrA   =  15.75,               // Length to radius ratio in small arteries.
   //  lrrV   =  14.54,               // Length to radius ratio in small veins.

   //  Fcst   = 10.0,                 // Determines the damping coeff.
//  double
//       Fcst   = 50, //17.7778,        // Determines the damping coeff.
                                        // for the friction.
//       Lr     = 0.25,                  // characteristic radius of the
                                        // vessels in the tree [cm].
//       Lr222    = sq(0.25),               // The squared radius [cm2].
//       Lr33    = cu(0.25),               // The radius to the third power [cm^3].
//       ggg      = 981.0,                // The gravitational force [cm/s^2].
//       qqq      = 10.0*Lr222,             // The characteristic flow [cm^3/s].
//       Fr212    = sq(qqq)/ggg/pow(0.25,5),    // The squared Froudes number.
//       Re     = q*1.055/0.049/0.25,          // Reynolds number.
//       Perioddd = 0.7*qqq/Lr33;           // The dimension-less period.
    // tau    = 0.08*q/Lr3,           // End of pulse, dimension-less.
    //   kkk      = Period/tmstps,        // Length of a timestep.
//       Deltat = Periodd/plts;          // Interval between each point plottet.
//       p0     = 0.0;//-10/rho/g/Lr*conv;      // Ensures a certain diastolic pressure. (Unstressed pressure)




//double tstart1=0,
//       tend1=0,
//	   finaltime1=0;

//double onePeriodDt=0.0001;


	/*
 double rm = 0.005; double f1 = 0; double f2 = 1; double f3 = 260000;
 double fa1 = 0; double fa2 = 1; double fa3 = 50000; double asym = 0.41;
 double expo = 2.76; double lrrA = 15.76; double lrrV = 14.54;
 bool verbosity = 1; double numHeartBeats = 6; int id = 1; double SVPA = 5;
 // SVPA: external pressure amplitude. Has dimensions. see sor06.h
 //double *Q0; // appears to be unused
*/


 // SVPA: external pressure amplitude. Has dimensions. see sor06.h
 //double *Q0; // appears to be unused var
 double rm; double f1; double f2; double f3; double fa1; double fa2;
 double fa3; double asym; double expo; double lrrA; double lrrV;
 int numHeartBeats; int id; double SVPA; double LVPA; int gridPoints;

static bool verbosity;


	double fv1; double fv2; double fv3;


   // variables for interaction
   double Snew[10];
   int iforsnew = 1;
   double Snewprint[10000];
   int countprint = 0;


  int total_subdomain1=63;

  // variables for interaction
  double qFromCoro[20] = {0};
  std::vector<double> Pext_coro1(total_subdomain1+2,10.0);


  double *fjac[18],*fmono[12], *ftrif[24],   // Work space used by bound_bif.
         *fj[8];      // Work space used by bound_match.

#endif
