//*******************************************************************/
/* This file contains the mono-, bi-, and trifurcation matching     */
/* routines for the 1D model code.                                  */
/* All the Jacobians were pre-inverted to skip the LU decomposition.*/
/* Developed by M J Colebank (NCSU) and J. A. Mackenzie (UoG)       */
//*******************************************************************/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <algorithm>
#include "tools.h"
#include "arteriesPD.h"
#include "junction.h"

void monofurcation(int qLnb, double theta, double gamma, Tube * Arteries[], int parent, int ntrial, double tol, int j, int ok) {
    //******************************************************************//
    //==================== MONOFURCATION MATCHING ======================//
    //******************************************************************//
    Tube * PV = Arteries[parent];
    Tube * LD = Arteries[parent] -> LD;
    int N = PV -> N;

    //    int j = 1;
    //    int ok = false;
    //    const int ntrial = 4000;
    //    double tol = 1.0e-10;

    double k1[2];
    double k2[2];
    double k3[2];
    double k4[2];
    // These are the flows at the half time step.
    k1[0] = PV -> Qold[N] + theta * (PV -> R2h[N - 1]) + gamma * (PV -> S2h[N - 1]);

    k1[1] = LD -> Qold[0] - theta * (LD -> R2h[0]) + gamma * (LD -> S2h[0]);

    k2[0] = PV -> Aold[N] + theta * (PV -> R1h[N - 1]);
    k2[1] = LD -> Aold[0] - theta * (LD -> R1h[0]);

    k3[0] = PV -> Qh[N - 1] / 2.0;
    k3[1] = LD -> Qh[0] / 2.0;

    k4[0] = PV -> Ah[N - 1] / 2.0;
    k4[1] = LD -> Ah[0] / 2.0;

    double xb[12];

    // The approximative initial guesses are applied.
    xb[0] = PV -> Qh[N - 1]; //Initial guess for Q1_xb n+1
    xb[1] = (PV -> Qold[N - 1] + PV -> Qold[N]) / 2.0; //Initial guess for Q1_xb^n+0.5
    xb[2] = PV -> Qold[N]; //Initial guess for Q1_xb+0.5 n+0.5
    xb[3] = LD -> Qh[0]; //Initial guess for Q2_xb n+1
    xb[4] = (LD -> Qold[0] + LD -> Qold[1]) / 2.0; //Initial guess for Q2_xb n+0.5
    xb[5] = LD -> Qold[0]; //Initial guess for Q2_xb+0.5 n+0.5
    xb[6] = PV -> Ah[N - 1]; //Initial guess for A1_xb n+1
    xb[7] = (PV -> Aold[N - 1] + PV -> Aold[N]) / 2.0; //Initial guess for A1_xb^n+0.5
    xb[8] = PV -> Aold[N]; //Initial guess for A1_xb+0.5 n+0.5
    xb[9] = LD -> Ah[0]; //Initial guess for A2_xb n+1
    xb[10] = (LD -> Aold[0] + LD -> Aold[1]) / 2.0; //Initial guess for A2_xb n+0.5
    xb[11] = LD -> Aold[0]; //Initial guess for A2_xb+0.5 n+0.5

    double k7nh = 0.0; // Bernoulli loss
    double k7n = 0.0; //

    while (j <= ntrial && ok == false) // Find the zero
    {
        double fvec[12];

        // The residuals

        // Characteristic Q residual at n+1
        fvec[0] = k1[0] - xb[0] -
            theta * (sq(xb[2]) / xb[8] + PV -> Bh(N, xb[8])) +
            gamma * (PV -> F(xb[2], xb[8]) + PV -> dBdx1h(N, xb[8]));

        fvec[1] = k1[1] - xb[3] +
            theta * (sq(xb[5]) / xb[11] + LD -> Bh(-1, xb[11])) +
            gamma * (PV -> F(xb[5], xb[11]) + LD -> dBdx1h(-1, xb[11]));

        // Characteristic A residual at n+1
        fvec[2] = -theta * xb[2] - xb[6] + k2[0];
        fvec[3] = theta * xb[5] - xb[9] + k2[1];

        // Flow residuals at n+1/2 (ghost points)
        fvec[4] = -xb[1] + xb[2] / 2.0 + k3[0];
        fvec[5] = -xb[4] + xb[5] / 2.0 + k3[1];

        // Area residuals at n+1/2 (ghost points)
        fvec[6] = -xb[7] + xb[8] / 2.0 + k4[0];
        fvec[7] = -xb[10] + xb[11] / 2.0 + k4[1];

        // Flow conservation residuals (n+1/2 and n+1)
        fvec[8] = -xb[1] + xb[4];
        fvec[9] = -xb[0] + xb[3];

        // Use these terms if you want to have Benoulli loss: u^2=(q/A)^2
        double sq211 = sq(xb[1] / xb[7]);
        double sq110 = sq(xb[0] / xb[6]);

        fvec[10] = -(PV -> P(qLnb, N, xb[7])) + LD -> P(qLnb, 0, xb[10]) + abs(k7nh) * sq211;

        fvec[11] = -(PV -> P(qLnb, N, xb[6])) + LD -> P(qLnb, 0, xb[9]) + abs(k7n) * sq110;

        double chi[8];

        // Here are the residuals for the characteristic matching for flow
        chi[0] = -2.0 * theta * xb[2] / xb[8] + gamma * PV -> dFdQ(xb[8]);
        chi[2] = 2.0 * theta * xb[5] / xb[11] + gamma * PV -> dFdQ(xb[11]);

        // Here are the residuals for the area characteristic matching
        chi[1] = theta * (sq(xb[2] / xb[8]) - PV -> dBdAh(N, xb[8])) +
            gamma * (PV -> dFdA(xb[2], xb[8]) + PV -> d2BdAdxh(N, xb[8]));

        chi[3] = theta * (-sq(xb[5] / xb[11]) + LD -> dBdAh(-1, xb[11])) +
            gamma * (PV -> dFdA(xb[5], xb[11]) + LD -> d2BdAdxh(-1, xb[11]));

        // Here is pressure conservation (n+1/2)
        chi[4] = -(PV -> dPdA(N, xb[7])) + sq(xb[1]) / cu(xb[7]) * (-2.0 * abs(k7n)); //Loss term
        chi[5] = LD -> dPdA(0, xb[10]);

        // Here is pressure conservation (n+1)
        chi[6] = -(PV -> dPdA(N, xb[6])) + sq(xb[0]) / cu(xb[6]) * (-2.0 * abs(k7nh)); //Loss term
        chi[7] = LD -> dPdA(0, xb[9]);

        // check the matrix is non-singular; if singular, exit
        double detFmono= (chi[1]*chi[5]*chi[7]) - (chi[1]*chi[5]*chi[6]) - (chi[3]*chi[4]*chi[6]) + (chi[3]*chi[4]*chi[7]);

        if (sqrt(sq(detFmono)) < tol) {
            error("arteries.c", "Singular Monofurcation Jacobian");
            exit(100);
        }

        for (int row = 0; row < 12; row++)
            for (int col = 0; col < 12; col++)
                fmono[row][col] = 0.0;
        // The Jacobian.
//        double chi_3_4 = chi[3]*chi[4];
//        double chi_1_5 = chi[1]*chi[5];
//        double chi_2_6 = chi[2]*chi[6];

        fmono[0][0] = -(chi[3]*chi[4])/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[0][1] = -(chi[1]*chi[5])/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[0][2] = -(chi[6]*(chi[0]*chi[3]*chi[4] + chi[1]*chi[2]*chi[5]))/(chi[1]*chi[5]*chi[6]*theta - chi[1]*chi[5]*chi[7]*theta + chi[3]*chi[4]*chi[6]*theta - chi[3]*chi[4]*chi[7]*theta);
        fmono[0][3] = -(chi[7]*(chi[0]*chi[3]*chi[4] + chi[1]*chi[2]*chi[5]))/(chi[1]*chi[5]*chi[6]*theta - chi[1]*chi[5]*chi[7]*theta + chi[3]*chi[4]*chi[6]*theta - chi[3]*chi[4]*chi[7]*theta);
        fmono[0][4] = -(2*(chi[0]*chi[3]*chi[4]*chi[7] + chi[1]*chi[2]*chi[5]*chi[6]))/(chi[1]*chi[5]*chi[6] - chi[1]*chi[5]*chi[7] + chi[3]*chi[4]*chi[6] - chi[3]*chi[4]*chi[7]);
        fmono[0][5] = (2*(chi[0]*chi[3]*chi[4]*chi[7] + chi[1]*chi[2]*chi[5]*chi[6]))/(chi[1]*chi[5]*chi[6] - chi[1]*chi[5]*chi[7] + chi[3]*chi[4]*chi[6] - chi[3]*chi[4]*chi[7]);
        fmono[0][6] = (2*chi[1]*chi[3]*chi[4])/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[0][7] = (2*chi[1]*chi[3]*chi[5])/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[0][8] = (2*(chi[0]*chi[3]*chi[4]*chi[7] + chi[1]*chi[2]*chi[5]*chi[6]))/(chi[1]*chi[5]*chi[6] - chi[1]*chi[5]*chi[7] + chi[3]*chi[4]*chi[6] - chi[3]*chi[4]*chi[7]);
        fmono[0][9] = -(chi[1]*chi[5])/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[0][10] = (2*chi[1]*chi[3])/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[0][11] = -(chi[0]*chi[3]*chi[4] + chi[1]*chi[2]*chi[5])/(chi[1]*chi[5]*chi[6]*theta - chi[1]*chi[5]*chi[7]*theta + chi[3]*chi[4]*chi[6]*theta - chi[3]*chi[4]*chi[7]*theta);


        fmono[1][0] = 0;
        fmono[1][1] = 0;
        fmono[1][2] = -chi[6]/(2*(chi[6]*theta - chi[7]*theta));
        fmono[1][3] = -chi[7]/(2*(chi[6]*theta - chi[7]*theta));
        fmono[1][4] = -chi[6]/(chi[6] - chi[7]);
        fmono[1][5] = chi[7]/(chi[6] - chi[7]);
        fmono[1][6] = 0;
        fmono[1][7] = 0;
        fmono[1][8] = chi[7]/(chi[6] - chi[7]);
        fmono[1][9] = 0;
        fmono[1][10] = 0;
        fmono[1][11] = -1/(2*(chi[6]*theta - chi[7]*theta));


        fmono[2][0] = 0;
        fmono[2][1] = 0;
        fmono[2][2] = -chi[6]/(chi[6]*theta - chi[7]*theta);
        fmono[2][3] = -chi[7]/(chi[6]*theta - chi[7]*theta);
        fmono[2][4] = -(2*chi[7])/(chi[6] - chi[7]);
        fmono[2][5] = (2*chi[7])/(chi[6] - chi[7]);
        fmono[2][6] = 0;
        fmono[2][7] = 0;
        fmono[2][8] = (2*chi[7])/(chi[6] - chi[7]);
        fmono[2][9] = 0;
        fmono[2][10] = 0;
        fmono[2][11] = -1/(chi[6]*theta - chi[7]*theta);


        fmono[3][0] = -(chi[3]*chi[4])/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[3][1] = -(chi[1]*chi[5])/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[3][2] = -(chi[6]*(chi[0]*chi[3]*chi[4] + chi[1]*chi[2]*chi[5]))/(chi[1]*chi[5]*chi[6]*theta - chi[1]*chi[5]*chi[7]*theta + chi[3]*chi[4]*chi[6]*theta - chi[3]*chi[4]*chi[7]*theta);
        fmono[3][3] = -(chi[7]*(chi[0]*chi[3]*chi[4] + chi[1]*chi[2]*chi[5]))/(chi[1]*chi[5]*chi[6]*theta - chi[1]*chi[5]*chi[7]*theta + chi[3]*chi[4]*chi[6]*theta - chi[3]*chi[4]*chi[7]*theta);
        fmono[3][4] = -(2*(chi[0]*chi[3]*chi[4]*chi[7] + chi[1]*chi[2]*chi[5]*chi[6]))/(chi[1]*chi[5]*chi[6] - chi[1]*chi[5]*chi[7] + chi[3]*chi[4]*chi[6] - chi[3]*chi[4]*chi[7]);
        fmono[3][5] = (2*(chi[0]*chi[3]*chi[4]*chi[7] + chi[1]*chi[2]*chi[5]*chi[6]))/(chi[1]*chi[5]*chi[6] - chi[1]*chi[5]*chi[7] + chi[3]*chi[4]*chi[6] - chi[3]*chi[4]*chi[7]);
        fmono[3][6] = (2*chi[1]*chi[3]*chi[4])/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[3][7] = (2*chi[1]*chi[3]*chi[5])/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[3][8] = (2*(chi[0]*chi[3]*chi[4]*chi[7] + chi[1]*chi[2]*chi[5]*chi[6]))/(chi[1]*chi[5]*chi[6] - chi[1]*chi[5]*chi[7] + chi[3]*chi[4]*chi[6] - chi[3]*chi[4]*chi[7]);
        fmono[3][9] = (chi[3]*chi[4])/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[3][10] = (2*chi[1]*chi[3])/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[3][11] = -(chi[0]*chi[3]*chi[4] + chi[1]*chi[2]*chi[5])/(chi[1]*chi[5]*chi[6]*theta - chi[1]*chi[5]*chi[7]*theta + chi[3]*chi[4]*chi[6]*theta - chi[3]*chi[4]*chi[7]*theta);


        fmono[4][0] = 0;
        fmono[4][1] = 0;
        fmono[4][2] = -chi[6]/(2*(chi[6]*theta - chi[7]*theta));
        fmono[4][3] = -chi[7]/(2*(chi[6]*theta - chi[7]*theta));
        fmono[4][4] = -chi[6]/(chi[6] - chi[7]);
        fmono[4][5] = chi[7]/(chi[6] - chi[7]);
        fmono[4][6] = 0;
        fmono[4][7] = 0;
        fmono[4][8] = chi[6]/(chi[6] - chi[7]);
        fmono[4][9] = 0;
        fmono[4][10] = 0;
        fmono[4][11] = -1/(2*(chi[6]*theta - chi[7]*theta));


        fmono[5][0] = 0;
        fmono[5][1] = 0;
        fmono[5][2] = -chi[6]/(chi[6]*theta - chi[7]*theta);
        fmono[5][3] = -chi[7]/(chi[6]*theta - chi[7]*theta);
        fmono[5][4] = -(2*chi[6])/(chi[6] - chi[7]);
        fmono[5][5] = (2*chi[6])/(chi[6] - chi[7]);
        fmono[5][6] = 0;
        fmono[5][7] = 0;
        fmono[5][8] = (2*chi[6])/(chi[6] - chi[7]);
        fmono[5][9] = 0;
        fmono[5][10] = 0;
        fmono[5][11] = -1/(chi[6]*theta - chi[7]*theta);


        fmono[6][0] = 0;
        fmono[6][1] = 0;
        fmono[6][2] = chi[7]/(chi[6] - chi[7]);
        fmono[6][3] = chi[7]/(chi[6] - chi[7]);
        fmono[6][4] = (2*chi[7]*theta)/(chi[6] - chi[7]);
        fmono[6][5] = -(2*chi[7]*theta)/(chi[6] - chi[7]);
        fmono[6][6] = 0;
        fmono[6][7] = 0;
        fmono[6][8] = -(2*chi[7]*theta)/(chi[6] - chi[7]);
        fmono[6][9] = 0;
        fmono[6][10] = 0;
        fmono[6][11] = 1/(chi[6] - chi[7]);


        fmono[7][0] = chi[5]/(2*(chi[1]*chi[5] + chi[3]*chi[4]));
        fmono[7][1] = -chi[5]/(2*(chi[1]*chi[5] + chi[3]*chi[4]));
        fmono[7][2] = (chi[5]*chi[6]*(chi[0] - chi[2]))/(2*(chi[1]*chi[5]*chi[6]*theta - chi[1]*chi[5]*chi[7]*theta + chi[3]*chi[4]*chi[6]*theta - chi[3]*chi[4]*chi[7]*theta));
        fmono[7][3] = (chi[5]*(chi[0]*chi[7] - chi[2]*chi[7]))/(2*(chi[1]*chi[5]*chi[6]*theta - chi[1]*chi[5]*chi[7]*theta + chi[3]*chi[4]*chi[6]*theta - chi[3]*chi[4]*chi[7]*theta));
        fmono[7][4] = (chi[5]*(chi[0]*chi[7] - chi[2]*chi[6]))/(chi[1]*chi[5]*chi[6] - chi[1]*chi[5]*chi[7] + chi[3]*chi[4]*chi[6] - chi[3]*chi[4]*chi[7]);
        fmono[7][5] = -(chi[5]*(chi[0]*chi[7] - chi[2]*chi[6]))/(chi[1]*chi[5]*chi[6] - chi[1]*chi[5]*chi[7] + chi[3]*chi[4]*chi[6] - chi[3]*chi[4]*chi[7]);
        fmono[7][6] = -(chi[1]*chi[5])/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[7][7] = (chi[3]*chi[5])/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[7][8] = -(chi[5]*(chi[0]*chi[7] - chi[2]*chi[6]))/(chi[1]*chi[5]*chi[6] - chi[1]*chi[5]*chi[7] + chi[3]*chi[4]*chi[6] - chi[3]*chi[4]*chi[7]);
        fmono[7][9] = -chi[5]/(2*(chi[1]*chi[5] + chi[3]*chi[4]));
        fmono[7][10] = chi[3]/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[7][11] = (chi[5]*(chi[0] - chi[2]))/(2*(chi[1]*chi[5]*chi[6]*theta - chi[1]*chi[5]*chi[7]*theta + chi[3]*chi[4]*chi[6]*theta - chi[3]*chi[4]*chi[7]*theta));


        fmono[8][0] = chi[5]/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[8][1] = -chi[5]/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[8][2] = (chi[5]*chi[6]*(chi[0] - chi[2]))/(chi[1]*chi[5]*chi[6]*theta - chi[1]*chi[5]*chi[7]*theta + chi[3]*chi[4]*chi[6]*theta - chi[3]*chi[4]*chi[7]*theta);
        fmono[8][3] = (chi[5]*(chi[0]*chi[7] - chi[2]*chi[7]))/(chi[1]*chi[5]*chi[6]*theta - chi[1]*chi[5]*chi[7]*theta + chi[3]*chi[4]*chi[6]*theta - chi[3]*chi[4]*chi[7]*theta);
        fmono[8][4] = (2*chi[5]*(chi[0]*chi[7] - chi[2]*chi[6]))/(chi[1]*chi[5]*chi[6] - chi[1]*chi[5]*chi[7] + chi[3]*chi[4]*chi[6] - chi[3]*chi[4]*chi[7]);
        fmono[8][5] = -(2*chi[5]*(chi[0]*chi[7] - chi[2]*chi[6]))/(chi[1]*chi[5]*chi[6] - chi[1]*chi[5]*chi[7] + chi[3]*chi[4]*chi[6] - chi[3]*chi[4]*chi[7]);
        fmono[8][6] = (2*chi[3]*chi[4])/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[8][7] = (2*chi[3]*chi[5])/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[8][8] = -(2*chi[5]*(chi[0]*chi[7] - chi[2]*chi[6]))/(chi[1]*chi[5]*chi[6] - chi[1]*chi[5]*chi[7] + chi[3]*chi[4]*chi[6] - chi[3]*chi[4]*chi[7]);
        fmono[8][9] = -chi[5]/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[8][10] = (2*chi[3])/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[8][11] = (chi[5]*(chi[0] - chi[2]))/(chi[1]*chi[5]*chi[6]*theta - chi[1]*chi[5]*chi[7]*theta + chi[3]*chi[4]*chi[6]*theta - chi[3]*chi[4]*chi[7]*theta);


        fmono[9][0] = 0;
        fmono[9][1] = 0;
        fmono[9][2] = -chi[6]/(chi[6] - chi[7]);
        fmono[9][3] = -chi[6]/(chi[6] - chi[7]);
        fmono[9][4] = -(2*chi[6]*theta)/(chi[6] - chi[7]);
        fmono[9][5] = (2*chi[6]*theta)/(chi[6] - chi[7]);
        fmono[9][6] = 0;
        fmono[9][7] = 0;
        fmono[9][8] = (2*chi[6]*theta)/(chi[6] - chi[7]);
        fmono[9][9] = 0;
        fmono[9][10] = 0;
        fmono[9][11] = -1/(chi[6] - chi[7]);


        fmono[10][0] = -chi[4]/(2*(chi[1]*chi[5] + chi[3]*chi[4]));
        fmono[10][1] = chi[4]/(2*(chi[1]*chi[5] + chi[3]*chi[4]));
        fmono[10][2] = -(chi[4]*chi[6]*(chi[0] - chi[2]))/(2*(chi[1]*chi[5]*chi[6]*theta - chi[1]*chi[5]*chi[7]*theta + chi[3]*chi[4]*chi[6]*theta - chi[3]*chi[4]*chi[7]*theta));
        fmono[10][3] = -(chi[0]*chi[4]*chi[7] - chi[2]*chi[4]*chi[7])/(2*(chi[1]*chi[5]*chi[6]*theta - chi[1]*chi[5]*chi[7]*theta + chi[3]*chi[4]*chi[6]*theta - chi[3]*chi[4]*chi[7]*theta));
        fmono[10][4] = -(chi[0]*chi[4]*chi[7] - chi[2]*chi[4]*chi[6])/(chi[1]*chi[5]*chi[6] - chi[1]*chi[5]*chi[7] + chi[3]*chi[4]*chi[6] - chi[3]*chi[4]*chi[7]);
        fmono[10][5] = (chi[0]*chi[4]*chi[7] - chi[2]*chi[4]*chi[6])/(chi[1]*chi[5]*chi[6] - chi[1]*chi[5]*chi[7] + chi[3]*chi[4]*chi[6] - chi[3]*chi[4]*chi[7]);
        fmono[10][6] = (chi[1]*chi[4])/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[10][7] = -(chi[3]*chi[4])/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[10][8] = (chi[0]*chi[4]*chi[7] - chi[2]*chi[4]*chi[6])/(chi[1]*chi[5]*chi[6] - chi[1]*chi[5]*chi[7] + chi[3]*chi[4]*chi[6] - chi[3]*chi[4]*chi[7]);
        fmono[10][9] = chi[4]/(2*(chi[1]*chi[5] + chi[3]*chi[4]));
        fmono[10][10] = chi[1]/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[10][11] = -(chi[4]*(chi[0] - chi[2]))/(2*(chi[1]*chi[5]*chi[6]*theta - chi[1]*chi[5]*chi[7]*theta + chi[3]*chi[4]*chi[6]*theta - chi[3]*chi[4]*chi[7]*theta));


        fmono[11][0] = -chi[4]/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[11][1] = chi[4]/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[11][2] = -(chi[4]*chi[6]*(chi[0] - chi[2]))/(chi[1]*chi[5]*chi[6]*theta - chi[1]*chi[5]*chi[7]*theta + chi[3]*chi[4]*chi[6]*theta - chi[3]*chi[4]*chi[7]*theta);
        fmono[11][3] = -(chi[0]*chi[4]*chi[7] - chi[2]*chi[4]*chi[7])/(chi[1]*chi[5]*chi[6]*theta - chi[1]*chi[5]*chi[7]*theta + chi[3]*chi[4]*chi[6]*theta - chi[3]*chi[4]*chi[7]*theta);
        fmono[11][4] = -(2*(chi[0]*chi[4]*chi[7] - chi[2]*chi[4]*chi[6]))/(chi[1]*chi[5]*chi[6] - chi[1]*chi[5]*chi[7] + chi[3]*chi[4]*chi[6] - chi[3]*chi[4]*chi[7]);
        fmono[11][5] = (2*(chi[0]*chi[4]*chi[7] - chi[2]*chi[4]*chi[6]))/(chi[1]*chi[5]*chi[6] - chi[1]*chi[5]*chi[7] + chi[3]*chi[4]*chi[6] - chi[3]*chi[4]*chi[7]);
        fmono[11][6] = (2*chi[1]*chi[4])/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[11][7] = (2*chi[1]*chi[5])/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[11][8] = (2*(chi[0]*chi[4]*chi[7] - chi[2]*chi[4]*chi[6]))/(chi[1]*chi[5]*chi[6] - chi[1]*chi[5]*chi[7] + chi[3]*chi[4]*chi[6] - chi[3]*chi[4]*chi[7]);
        fmono[11][9] = chi[4]/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[11][10] = (2*chi[1])/(chi[1]*chi[5] + chi[3]*chi[4]);
        fmono[11][11] = -(chi[4]*(chi[0] - chi[2]))/(chi[1]*chi[5]*chi[6]*theta - chi[1]*chi[5]*chi[7]*theta + chi[3]*chi[4]*chi[6]*theta - chi[3]*chi[4]*chi[7]*theta);


        // Check whether solution is close enough. If not run the loop again.
//        int ch = zero(xb, 12, tol, tol, fvec, fmono);
        int ch = zeroNoDcmp(xb, 12, tol, tol, fvec, fmono); //
        if (ch == 1) ok = true;
        j = j + 1;
    }

    // Solutions is applied, and right boundary is updated.
    PV -> Anew[N] = xb[6];
    PV -> Qnew[N] = xb[0];
    LD -> Anew[0] = xb[9];
    LD -> Qnew[0] = xb[3];

    if (j >= ntrial) {printf("Error (junction.C)\nRoot not found in the monofurcation\nArteries[%i]\nqLnb = %i", parent, qLnb); exit(200);}
}

void bifurcation(int qLnb, double theta, double gamma, Tube * Arteries[], int parent, int ntrial, double tol, int j, int ok) {
    //******************************************************************//
    //===================== BIFURCATION MATCHING =======================//
    //******************************************************************//
    Tube * PV = Arteries[parent];
    Tube * LD = Arteries[parent] -> LD;
    Tube * RD = Arteries[parent] -> RD;
    int N = PV -> N;

    double PN; // FLP. Defined otherwhere.
    double g1 = PV -> Qold[N] + theta * PV -> R2h[N - 1] + gamma * PV -> S2h[N - 1];
    double g2 = LD -> Qold[0] - theta * (LD -> R2h[0]) + gamma * (LD -> S2h[0]);
    double g2a = RD -> Qold[0] - theta * (RD -> R2h[0]) + gamma * (RD -> S2h[0]);

    double k1 = PV -> Aold[N] + theta * (PV -> R1h[N - 1]);
    double k2 = LD -> Aold[0] - theta * (LD -> R1h[0]);
    double k2a = RD -> Aold[0] - theta * (RD -> R1h[0]);

    double k3 = PV -> Qh[N - 1] / 2.0;
    double k4 = LD -> Qh[0] / 2.0;
    double k4a = RD -> Qh[0] / 2.0;

    double k5 = PV -> Ah[N - 1] / 2.0;
    double k6 = LD -> Ah[0] / 2.0;
    double k6a = RD -> Ah[0] / 2.0;

    double xb[18];

    // The approximative initial guesses are applied.
    xb[0] = PV -> Qh[N - 1]; //Initial guess for Q1_xb n+1
    xb[1] = (PV -> Qold[N - 1] + PV -> Qold[N]) / 2.0; //Initial guess for Q1_xb^n+0.5
    xb[2] = PV -> Qold[N]; //Initial guess for Q1_xb+0.5 n+0.5
    xb[3] = LD -> Qh[0]; //Initial guess for Q2_xb n+1
    xb[4] = (LD -> Qold[0] + LD -> Qold[1]) / 2.0; //Initial guess for Q2_xb n+0.5
    xb[5] = LD -> Qold[0]; //Initial guess for Q2_xb+0.5 n+0.5
    xb[6] = RD -> Qh[0]; //Initial guess for Q3_xb n+1
    xb[7] = (RD -> Qold[0] + RD -> Qold[1]) / 2.0; //Initial guess for Q3_xb n+0.5
    xb[8] = RD -> Qold[0]; //Initial guess for Q3_xb+0.5 n+0.5
    xb[9] = PV -> Ah[N - 1]; //Initial guess for A1_xb n+1
    xb[10] = (PV -> Aold[N - 1] + PV -> Aold[N]) / 2.0; //Initial guess for A1_xb^n+0.5
    xb[11] = PV -> Aold[N]; //Initial guess for A1_xb+0.5 n+0.5
    xb[12] = LD -> Ah[0]; //Initial guess for A2_xb n+1
    xb[13] = (LD -> Aold[0] + LD -> Aold[1]) / 2.0; //Initial guess for A2_xb n+0.5
    xb[14] = LD -> Aold[0]; //Initial guess for A2_xb+0.5 n+0.5
    xb[15] = RD -> Ah[0]; //Initial guess for A3_xb n+1
    xb[16] = (RD -> Aold[0] + RD -> Aold[1]) / 2.0; //Initial guess for A3_xb n+0.5
    xb[17] = RD -> Aold[0]; //Initial guess for A3_xb+0.5 n+0.5

    // Bernoulli loss terms
    double k7nh = LD -> K_loss / 2.0;
    double k7n = LD -> K_loss / 2.0;
    double k7anh = RD -> K_loss / 2.0;
    double k7an = RD -> K_loss / 2.0;

    // The residuals (fvec), and the Jacobian is determined, and if possible
    // the system of equations is solved.
    while (j <= ntrial && ok == false) // Find the zero
    {
        double fvec[18];
        // The residuals.
        fvec[0] = g1 - xb[0] - theta * (sq(xb[2]) / xb[11] + PV -> Bh(N, xb[11])) + gamma * (PV -> F(xb[2], xb[11]) + PV -> dBdx1h(N, xb[11]));
        fvec[1] = g2 - xb[3] + theta * (sq(xb[5]) / xb[14] + LD -> Bh(-1, xb[14])) + gamma * (PV -> F(xb[5], xb[14]) + LD -> dBdx1h(-1, xb[14]));
        fvec[2] = g2a - xb[6] + theta * (sq(xb[8]) / xb[17] + RD -> Bh(-1, xb[17])) + gamma * (PV -> F(xb[8], xb[17]) + RD -> dBdx1h(-1, xb[17]));
        fvec[3] = -theta * xb[2] - xb[9] + k1;
        fvec[4] = theta * xb[5] - xb[12] + k2;
        fvec[5] = theta * xb[8] - xb[15] + k2a;
        fvec[6] = -xb[1] + xb[2] / 2.0 + k3;
        fvec[7] = -xb[4] + xb[5] / 2.0 + k4;
        fvec[8] = -xb[7] + xb[8] / 2.0 + k4a;
        fvec[9] = -xb[10] + xb[11] / 2.0 + k5;
        fvec[10] = -xb[13] + xb[14] / 2.0 + k6;
        fvec[11] = -xb[16] + xb[17] / 2.0 + k6a;
        fvec[12] = -xb[1] + xb[4] + xb[7];
        fvec[13] = -xb[0] + xb[3] + xb[6];

        PN = PV -> P(qLnb, N, xb[10]); // FLP
        double sq211 = sq(xb[1] / xb[10]);
        if (xb[1] > 0) {
            fvec[14] = -PN + LD -> P(qLnb, 0, xb[13]) + k7nh * sq211;
            fvec[15] = -PN + RD -> P(qLnb, 0, xb[16]) + k7anh * sq211;
        } else {
            fvec[14] = -PN + LD -> P(qLnb, 0, xb[13]) - k7nh * sq211;
            fvec[15] = -PN + RD -> P(qLnb, 0, xb[16]) - k7anh * sq211;
        };

        PN = PV -> P(qLnb, N, xb[9]); // FLP
        double sq110 = sq(xb[0] / xb[9]);
        if (xb[0] > 0) {
            fvec[16] = -PN + LD -> P(qLnb, 0, xb[12]) + k7n * sq110;
            fvec[17] = -PN + RD -> P(qLnb, 0, xb[15]) + k7an * sq110;
        } else {
            fvec[16] = -PN + LD -> P(qLnb, 0, xb[12]) - k7n * sq110;
            fvec[17] = -PN + RD -> P(qLnb, 0, xb[15]) - k7an * sq110;
        };

        double chi[12];
        chi[0] = -2 * theta * xb[2] / xb[11] + gamma * (PV -> dFdQ(xb[11]));
        chi[1] = theta * (sq(xb[2] / xb[11]) - (PV -> dBdAh(N, xb[11]))) +
            gamma * ((PV -> dFdA(xb[2], xb[11])) + (PV -> d2BdAdxh(N, xb[11])));
        chi[2] = 2 * theta * xb[5] / xb[14] + gamma * PV -> dFdQ(xb[14]);
        chi[3] = theta * (-sq(xb[5] / xb[14]) + LD -> dBdAh(-1, xb[14])) +
            gamma * (PV -> dFdA(xb[5], xb[14]) + LD -> d2BdAdxh(-1, xb[14]));
        chi[4] = 2 * theta * xb[8] / xb[17] + gamma * PV -> dFdQ(xb[17]);
        chi[5] = theta * (-sq(xb[8] / xb[17]) + RD -> dBdAh(-1, xb[17])) +
            gamma * (PV -> dFdA(xb[8], xb[17]) + RD -> d2BdAdxh(-1, xb[17]));
        if (xb[1] > 0) {
            chi[6] = -(PV -> dPdA(N, xb[10])) + sq(xb[1]) / cu(xb[10]) * (-2 * k7nh);
        } else {
            chi[6] = -(PV -> dPdA(N, xb[10])) + sq(xb[1]) / cu(xb[10]) * (2 * k7nh);
        };
        chi[7] = LD -> dPdA(0, xb[13]);
        chi[8] = RD -> dPdA(0, xb[16]);
        if (xb[0] > 0) {
            chi[9] = -(PV -> dPdA(N, xb[9])) + sq(xb[0]) / cu(xb[9]) * (-2 * k7an);
        } else {
            chi[9] = -(PV -> dPdA(N, xb[9])) + sq(xb[0]) / cu(xb[9]) * (2 * k7n);
        };
        chi[10] = LD -> dPdA(0, xb[12]);
        chi[11] = RD -> dPdA(0, xb[15]);

        //for (int i=0; i<12; i++){
        //  fprintf(stdout, "chi(%i) = %10.15f \n",i+1, chi[i]);
        //};
        // fprintf(stdout, "theta = %10.15f \n",theta);
        // exit(1);

        // The Inverse Jacobian.
        // many elements are repeated, so only compute the first instance and duplicate the rest as needed per the list starting from ln. 1726

        double detFjac = ((chi[1] * chi[7] * chi[8] + chi[3] * chi[6] * chi[8] + chi[5] * chi[6] * chi[7]) * (chi[9] * chi[10] + chi[9] * chi[11] - chi[10] * chi[11]));
        if (sqrt(sq(detFjac)) < 1.0e-10) {
            error("arteries.c", "Singular Bifurcation Jacobian");
            exit(101);
        }

        double chi368 = chi[3] * chi[6] * chi[8];
        double chi178 = chi[1] * chi[7] * chi[8];
        double chi567 = chi[5] * chi[6] * chi[7];
        double chi910 = chi[9] * chi[10];
        double chi911 = chi[9] * chi[11];
        double chi12 = chi[1] * chi[2];
        double chi789 = chi[7] * chi[8] * chi[9];
        double chi1011 = chi[10] * chi[11];
        double chi34 = chi[3] * chi[4];

        for (int row = 0; row < 18; row++)
            for (int col = 0; col < 18; col++)
                fjac[row][col] = 0.0;

        fjac[0][0] = -(chi368 + chi567) / (chi178 + chi368 + chi567);
        fjac[0][1] = -(chi178) / (chi178 + chi368 + chi567);
        fjac[0][3] = -(chi[0] * chi368 * chi910 + chi[0] * chi368 * chi911 + chi[0] * chi567 * chi910 + chi[0] * chi567 * chi911 + chi12 * chi[7] * chi[8] * chi911 + chi[1] * chi[4] * chi[7] * chi[8] * chi910) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);
        fjac[0][4] = -(chi[0] * chi368 * chi1011 - chi12 * chi[7] * chi[8] * chi910 + chi[0] * chi567 * chi1011 + chi12 * chi[7] * chi[8] * chi1011 + chi[1] * chi[4] * chi[7] * chi[8] * chi910) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);
        fjac[0][5] = -(chi[0] * chi368 * chi1011 + chi12 * chi[7] * chi[8] * chi911 + chi[0] * chi567 * chi1011 - chi[1] * chi[4] * chi[7] * chi[8] * chi911 + chi[1] * chi[4] * chi[7] * chi[8] * chi1011) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);
        fjac[0][6] = -(2 * (chi[0] * chi368 * chi1011 + chi12 * chi[7] * chi[8] * chi911 + chi[0] * chi567 * chi1011 + chi[1] * chi[4] * chi[7] * chi[8] * chi910)) / (chi178 * chi910 + chi178 * chi911 + chi368 * chi910 - chi178 * chi1011 + chi368 * chi911 + chi567 * chi910 - chi368 * chi1011 + chi567 * chi911 - chi567 * chi1011);
        fjac[0][9] = (2 * (chi[1] * chi368 + chi[1] * chi567)) / (chi178 + chi368 + chi567);
        fjac[0][10] = (2 * chi[1] * chi[3] * chi[7] * chi[8]) / (chi178 + chi368 + chi567);
        fjac[0][11] = (2 * chi[1] * chi[5] * chi[7] * chi[8]) / (chi178 + chi368 + chi567);
        fjac[0][14] = (2 * chi[1] * chi[3] * chi[8]) / (chi178 + chi368 + chi567);
        fjac[0][15] = (2 * chi[1] * chi[5] * chi[7]) / (chi178 + chi368 + chi567);
        fjac[0][16] = -(chi[0] * chi368 * chi[11] - chi12 * chi789 + chi[0] * chi567 * chi[11] + chi12 * chi[7] * chi[8] * chi[11] + chi[1] * chi[4] * chi789) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);
        fjac[0][17] = -(chi[0] * chi368 * chi[10] + chi12 * chi789 + chi[0] * chi567 * chi[10] - chi[1] * chi[4] * chi789 + chi[1] * chi[4] * chi[7] * chi[8] * chi[10]) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);

        fjac[1][3] = -(chi910 + chi911) / (2 * (chi910 * theta + chi911 * theta - chi1011 * theta));
        fjac[1][4] = -(chi1011) / (2 * (chi910 * theta + chi911 * theta - chi1011 * theta));
        fjac[1][6] = -(chi910 + chi911) / (chi910 + chi911 - chi1011);
        fjac[1][16] = -chi[11] / (2 * (chi910 * theta + chi911 * theta - chi1011 * theta));
        fjac[1][17] = -chi[10] / (2 * (chi910 * theta + chi911 * theta - chi1011 * theta));

        fjac[2][3] = -(chi910 + chi911) / (chi910 * theta + chi911 * theta - chi1011 * theta);
        fjac[2][4] = -(chi1011) / (chi910 * theta + chi911 * theta - chi1011 * theta);
        fjac[2][6] = -(2 * chi1011) / (chi910 + chi911 - chi1011);
        fjac[2][16] = -chi[11] / (chi910 * theta + chi911 * theta - chi1011 * theta);
        fjac[2][17] = -chi[10] / (chi910 * theta + chi911 * theta - chi1011 * theta);

        fjac[3][0] = -(chi368) / (chi178 + chi368 + chi567);
        fjac[3][1] = -(chi178 + chi567) / (chi178 + chi368 + chi567);
        fjac[3][3] = -(chi[0] * chi368 * chi910 + chi[0] * chi368 * chi911 + chi12 * chi[7] * chi[8] * chi911 + chi[2] * chi567 * chi911 - chi34 * chi[6] * chi[8] * chi910) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);
        fjac[3][4] = (chi12 * chi[7] * chi[8] * chi910 - chi[0] * chi368 * chi1011 - chi12 * chi[7] * chi[8] * chi1011 + chi[2] * chi567 * chi910 + chi34 * chi[6] * chi[8] * chi910 - chi[2] * chi567 * chi1011) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);
        fjac[3][5] = -(chi[0] * chi368 * chi1011 + chi12 * chi[7] * chi[8] * chi911 + chi[2] * chi567 * chi911 + chi34 * chi[6] * chi[8] * chi911 - chi34 * chi[6] * chi[8] * chi1011) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);
        fjac[3][6] = -(2 * (chi[0] * chi368 * chi1011 + chi12 * chi[7] * chi[8] * chi911 + chi[2] * chi567 * chi911 - chi34 * chi[6] * chi[8] * chi910)) / (chi178 * chi910 + chi178 * chi911 + chi368 * chi910 - chi178 * chi1011 + chi368 * chi911 + chi567 * chi910 - chi368 * chi1011 + chi567 * chi911 - chi567 * chi1011);
        fjac[3][9] = (2 * chi[1] * chi368) / (chi178 + chi368 + chi567);
        fjac[3][10] = (2 * (chi[1] * chi[3] * chi[7] * chi[8] + chi[3] * chi567)) / (chi178 + chi368 + chi567);
        fjac[3][11] = -(2 * chi[3] * chi[5] * chi[6] * chi[8]) / (chi178 + chi368 + chi567);
        fjac[3][14] = (2 * (chi[1] * chi[3] * chi[8] + chi[3] * chi[5] * chi[6])) / (chi178 + chi368 + chi567);
        fjac[3][16] = (chi12 * chi789 - chi[0] * chi368 * chi[11] - chi12 * chi[7] * chi[8] * chi[11] + chi[2] * chi567 * chi[9] + chi34 * chi[6] * chi[8] * chi[9] - chi[2] * chi567 * chi[11]) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);
        fjac[3][17] = -(chi[0] * chi368 * chi[10] + chi12 * chi789 + chi[2] * chi567 * chi[9] + chi34 * chi[6] * chi[8] * chi[9] - chi34 * chi[6] * chi[8] * chi[10]) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);

        fjac[4][3] = -(chi911) / (2 * (chi910 * theta + chi911 * theta - chi1011 * theta));
        fjac[4][4] = (chi910 - chi1011) / (2 * (chi910 * theta + chi911 * theta - chi1011 * theta));
        fjac[4][7] = -(chi910 - chi1011) / (chi910 + chi911 - chi1011);
        fjac[4][16] = (chi[9] - chi[11]) / (2 * (chi910 * theta + chi911 * theta - chi1011 * theta));

        fjac[5][3] = -(chi911) / (chi910 * theta + chi911 * theta - chi1011 * theta);
        fjac[5][4] = (chi910 - chi1011) / (chi910 * theta + chi911 * theta - chi1011 * theta);
        fjac[5][6] = -(2 * chi911) / (chi910 + chi911 - chi1011);
        fjac[5][16] = (chi[9] - chi[11]) / (chi910 * theta + chi911 * theta - chi1011 * theta);

        fjac[6][0] = -(chi567) / (chi178 + chi368 + chi567);
        fjac[6][2] = -(chi178 + chi368) / (chi178 + chi368 + chi567);
        fjac[6][3] = -(chi[0] * chi567 * chi910 + chi[0] * chi567 * chi911 + chi[1] * chi[4] * chi[7] * chi[8] * chi910 - chi[2] * chi567 * chi911 + chi34 * chi[6] * chi[8] * chi910) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);
        fjac[6][4] = -(chi[0] * chi567 * chi1011 + chi[1] * chi[4] * chi[7] * chi[8] * chi910 + chi[2] * chi567 * chi910 + chi34 * chi[6] * chi[8] * chi910 - chi[2] * chi567 * chi1011) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);
        fjac[6][5] = -(chi[0] * chi567 * chi1011 - chi[1] * chi[4] * chi[7] * chi[8] * chi911 - chi[2] * chi567 * chi911 + chi[1] * chi[4] * chi[7] * chi[8] * chi1011 - chi34 * chi[6] * chi[8] * chi911 + chi34 * chi[6] * chi[8] * chi1011) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);
        fjac[6][6] = -(2 * (chi[0] * chi567 * chi1011 + chi[1] * chi[4] * chi[7] * chi[8] * chi910 - chi[2] * chi567 * chi911 + chi34 * chi[6] * chi[8] * chi910)) / (chi178 * chi910 + chi178 * chi911 + chi368 * chi910 - chi178 * chi1011 + chi368 * chi911 + chi567 * chi910 - chi368 * chi1011 + chi567 * chi911 - chi567 * chi1011);
        fjac[6][9] = (2 * chi[1] * chi567) / (chi178 + chi368 + chi567);
        fjac[6][10] = -(2 * chi[3] * chi567) / (chi178 + chi368 + chi567);
        fjac[6][11] = (2 * chi[5] * (chi178 + chi368)) / (chi178 + chi368 + chi567);
        fjac[6][14] = -(2 * chi[3] * chi[5] * chi[6]) / (chi178 + chi368 + chi567);
        fjac[6][15] = (2 * chi[5] * (chi[1] * chi[7] + chi[3] * chi[6])) / (chi178 + chi368 + chi567);
        fjac[6][16] = -(chi[0] * chi567 * chi[11] + chi[1] * chi[4] * chi789 + chi[2] * chi567 * chi[9] + chi34 * chi[6] * chi[8] * chi[9] - chi[2] * chi567 * chi[11]) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);
        fjac[6][17] = -(chi[0] * chi567 * chi[10] - chi[1] * chi[4] * chi789 - chi[2] * chi567 * chi[9] + chi[1] * chi[4] * chi[7] * chi[8] * chi[10] - chi34 * chi[6] * chi[8] * chi[9] + chi34 * chi[6] * chi[8] * chi[10]) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);

        fjac[7][3] = -(chi910) / (2 * (chi910 * theta + chi911 * theta - chi1011 * theta));
        fjac[7][5] = (chi911 - chi1011) / (2 * (chi910 * theta + chi911 * theta - chi1011 * theta));
        fjac[7][8] = -(chi911 - chi1011) / (chi910 + chi911 - chi1011);
        fjac[7][16] = -chi[9] / (2 * (chi910 * theta + chi911 * theta - chi1011 * theta));
        fjac[7][17] = (chi[9] - chi[10]) / (2 * (chi910 * theta + chi911 * theta - chi1011 * theta));

        fjac[8][3] = -(chi910) / (chi910 * theta + chi911 * theta - chi1011 * theta);
        fjac[8][5] = (chi911 - chi1011) / (chi910 * theta + chi911 * theta - chi1011 * theta);
        fjac[8][6] = -(2 * chi910) / (chi910 + chi911 - chi1011);
        fjac[8][16] = -chi[9] / (chi910 * theta + chi911 * theta - chi1011 * theta);
        fjac[8][17] = (chi[9] - chi[10]) / (chi910 * theta + chi911 * theta - chi1011 * theta);

        fjac[9][3] = (chi1011) / (chi910 + chi911 - chi1011);
        fjac[9][6] = (2 * chi1011 * theta) / (chi910 + chi911 - chi1011);
        fjac[9][16] = chi[11] / (chi910 + chi911 - chi1011);
        fjac[9][17] = chi[10] / (chi910 + chi911 - chi1011);

        fjac[10][0] = (chi[7] * chi[8]) / (2 * (chi178 + chi368 + chi567));
        fjac[10][3] = (chi[0] * chi[7] * chi[8] * chi910 + chi[0] * chi[7] * chi[8] * chi911 - chi[2] * chi[7] * chi[8] * chi911 - chi[4] * chi[7] * chi[8] * chi910) / (2 * (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta));
        fjac[10][4] = (chi[0] * chi[7] * chi[8] * chi1011 + chi[2] * chi[7] * chi[8] * chi910 - chi[2] * chi[7] * chi[8] * chi1011 - chi[4] * chi[7] * chi[8] * chi910) / (2 * (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta));
        fjac[10][5] = (chi[0] * chi[7] * chi[8] * chi1011 - chi[2] * chi[7] * chi[8] * chi911 + chi[4] * chi[7] * chi[8] * chi911 - chi[4] * chi[7] * chi[8] * chi1011) / (2 * (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta));
        fjac[10][6] = -(chi[2] * chi[7] * chi[8] * chi911 - chi[0] * chi[7] * chi[8] * chi1011 + chi[4] * chi[7] * chi[8] * chi910) / (chi178 * chi910 + chi178 * chi911 + chi368 * chi910 - chi178 * chi1011 + chi368 * chi911 + chi567 * chi910 - chi368 * chi1011 + chi567 * chi911 - chi567 * chi1011);
        fjac[10][10] = (chi[3] * chi[7] * chi[8]) / (chi178 + chi368 + chi567);
        fjac[10][11] = (chi[5] * chi[7] * chi[8]) / (chi178 + chi368 + chi567);
        fjac[10][14] = (chi[3] * chi[8]) / (chi178 + chi368 + chi567);
        fjac[10][15] = (chi[5] * chi[7]) / (chi178 + chi368 + chi567);
        fjac[10][16] = (chi[0] * chi[7] * chi[8] * chi[11] + chi[2] * chi789 - chi[2] * chi[7] * chi[8] * chi[11] - chi[4] * chi789) / (2 * (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta));
        fjac[10][17] = (chi[0] * chi[7] * chi[8] * chi[10] - chi[2] * chi789 + chi[4] * chi789 - chi[4] * chi[7] * chi[8] * chi[10]) / (2 * (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta));

        fjac[11][0] = (chi[7] * chi[8]) / (chi178 + chi368 + chi567);
        fjac[11][3] = (chi[0] * chi[7] * chi[8] * chi910 + chi[0] * chi[7] * chi[8] * chi911 - chi[2] * chi[7] * chi[8] * chi911 - chi[4] * chi[7] * chi[8] * chi910) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);
        fjac[11][4] = (chi[0] * chi[7] * chi[8] * chi1011 + chi[2] * chi[7] * chi[8] * chi910 - chi[2] * chi[7] * chi[8] * chi1011 - chi[4] * chi[7] * chi[8] * chi910) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);
        fjac[11][5] = (chi[0] * chi[7] * chi[8] * chi1011 - chi[2] * chi[7] * chi[8] * chi911 + chi[4] * chi[7] * chi[8] * chi911 - chi[4] * chi[7] * chi[8] * chi1011) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);
        fjac[11][6] = -(2 * (chi[2] * chi[7] * chi[8] * chi911 - chi[0] * chi[7] * chi[8] * chi1011 + chi[4] * chi[7] * chi[8] * chi910)) / (chi178 * chi910 + chi178 * chi911 + chi368 * chi910 - chi178 * chi1011 + chi368 * chi911 + chi567 * chi910 - chi368 * chi1011 + chi567 * chi911 - chi567 * chi1011);
        fjac[11][9] = (2 * (chi368 + chi567)) / (chi178 + chi368 + chi567);
        fjac[11][10] = (2 * chi[3] * chi[7] * chi[8]) / (chi178 + chi368 + chi567);
        fjac[11][11] = (2 * chi[5] * chi[7] * chi[8]) / (chi178 + chi368 + chi567);
        fjac[11][14] = (2 * chi[3] * chi[8]) / (chi178 + chi368 + chi567);
        fjac[11][15] = (2 * chi[5] * chi[7]) / (chi178 + chi368 + chi567);
        fjac[11][16] = (chi[0] * chi[7] * chi[8] * chi[11] + chi[2] * chi789 - chi[2] * chi[7] * chi[8] * chi[11] - chi[4] * chi789) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);
        fjac[11][17] = (chi[0] * chi[7] * chi[8] * chi[10] - chi[2] * chi789 + chi[4] * chi789 - chi[4] * chi[7] * chi[8] * chi[10]) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);

        fjac[12][3] = -(chi911) / (chi910 + chi911 - chi1011);
        fjac[12][6] = -(2 * chi911 * theta) / (chi910 + chi911 - chi1011);
        fjac[12][16] = (chi[9] - chi[11]) / (chi910 + chi911 - chi1011);

        fjac[13][0] = -(chi[6] * chi[8]) / (2 * (chi178 + chi368 + chi567));
        fjac[13][3] = -(chi[8] * (chi[0] * chi[6] * chi910 + chi[0] * chi[6] * chi911 - chi[2] * chi[6] * chi911 - chi[4] * chi[6] * chi910)) / (2 * (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta));
        fjac[13][4] = -(chi[8] * chi[10] * (chi[0] * chi[6] * chi[11] + chi[2] * chi[6] * chi[9] - chi[2] * chi[6] * chi[11] - chi[4] * chi[6] * chi[9])) / (2 * (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta));
        fjac[13][5] = -(chi[8] * (chi[0] * chi[6] * chi1011 - chi[2] * chi[6] * chi911 + chi[4] * chi[6] * chi911 - chi[4] * chi[6] * chi1011)) / (2 * (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta));
        fjac[13][6] = (chi[8] * (chi[2] * chi[6] * chi911 - chi[0] * chi[6] * chi1011 + chi[4] * chi[6] * chi910)) / (chi178 * chi910 + chi178 * chi911 + chi368 * chi910 - chi178 * chi1011 + chi368 * chi911 + chi567 * chi910 - chi368 * chi1011 + chi567 * chi911 - chi567 * chi1011);
        fjac[13][9] = (chi[1] * chi[6] * chi[8]) / (chi178 + chi368 + chi567);
        fjac[13][11] = -(chi[5] * chi[6] * chi[8]) / (chi178 + chi368 + chi567);
        fjac[13][14] = (chi[1] * chi[8] + chi[5] * chi[6]) / (chi178 + chi368 + chi567);
        fjac[13][15] = -(chi[5] * chi[6]) / (chi178 + chi368 + chi567);
        fjac[13][16] = -(chi[8] * (chi[0] * chi[6] * chi[11] + chi[2] * chi[6] * chi[9] - chi[2] * chi[6] * chi[11] - chi[4] * chi[6] * chi[9])) / (2 * (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta));
        fjac[13][17] = -(chi[8] * (chi[0] * chi[6] * chi[10] - chi[2] * chi[6] * chi[9] + chi[4] * chi[6] * chi[9] - chi[4] * chi[6] * chi[10])) / (2 * (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta));

        fjac[14][0] = -(chi[6] * chi[8]) / (chi178 + chi368 + chi567);
        fjac[14][3] = -(chi[8] * (chi[0] * chi[6] * chi910 + chi[0] * chi[6] * chi911 - chi[2] * chi[6] * chi911 - chi[4] * chi[6] * chi910)) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);
        fjac[14][4] = -(chi[8] * chi[10] * (chi[0] * chi[6] * chi[11] + chi[2] * chi[6] * chi[9] - chi[2] * chi[6] * chi[11] - chi[4] * chi[6] * chi[9])) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);
        fjac[14][5] = -(chi[8] * (chi[0] * chi[6] * chi1011 - chi[2] * chi[6] * chi911 + chi[4] * chi[6] * chi911 - chi[4] * chi[6] * chi1011)) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);
        fjac[14][6] = (2 * chi[8] * (chi[2] * chi[6] * chi911 - chi[0] * chi[6] * chi1011 + chi[4] * chi[6] * chi910)) / (chi178 * chi910 + chi178 * chi911 + chi368 * chi910 - chi178 * chi1011 + chi368 * chi911 + chi567 * chi910 - chi368 * chi1011 + chi567 * chi911 - chi567 * chi1011);
        fjac[14][9] = (2 * chi[1] * chi[6] * chi[8]) / (chi178 + chi368 + chi567);
        fjac[14][10] = (2 * (chi178 + chi567)) / (chi178 + chi368 + chi567);
        fjac[14][11] = -(2 * chi[5] * chi[6] * chi[8]) / (chi178 + chi368 + chi567);
        fjac[14][14] = (2 * (chi[1] * chi[8] + chi[5] * chi[6])) / (chi178 + chi368 + chi567);
        fjac[14][15] = -(2 * chi[5] * chi[6]) / (chi178 + chi368 + chi567);
        fjac[14][16] = -(chi[8] * (chi[0] * chi[6] * chi[11] + chi[2] * chi[6] * chi[9] - chi[2] * chi[6] * chi[11] - chi[4] * chi[6] * chi[9])) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);
        fjac[14][17] = -(chi[8] * (chi[0] * chi[6] * chi[10] - chi[2] * chi[6] * chi[9] + chi[4] * chi[6] * chi[9] - chi[4] * chi[6] * chi[10])) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);

        fjac[15][3] = -(chi910) / (chi910 + chi911 - chi1011);
        fjac[15][6] = -(2 * chi910 * theta) / (chi910 + chi911 - chi1011);
        fjac[15][16] = -chi[9] / (chi910 + chi911 - chi1011);
        fjac[15][17] = (chi[9] - chi[10]) / (chi910 + chi911 - chi1011);

        fjac[16][0] = -(chi[6] * chi[7]) / (2 * (chi178 + chi368 + chi567));
        fjac[16][3] = -(chi[0] * chi[6] * chi[7] * chi910 + chi[0] * chi[6] * chi[7] * chi911 - chi[2] * chi[6] * chi[7] * chi911 - chi[4] * chi[6] * chi[7] * chi910) / (2 * (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta));
        fjac[16][4] = -(chi[7] * chi[10] * (chi[0] * chi[6] * chi[11] + chi[2] * chi[6] * chi[9] - chi[2] * chi[6] * chi[11] - chi[4] * chi[6] * chi[9])) / (2 * (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta));
        fjac[16][5] = -(chi[0] * chi[6] * chi[7] * chi1011 - chi[2] * chi[6] * chi[7] * chi911 + chi[4] * chi[6] * chi[7] * chi911 - chi[4] * chi[6] * chi[7] * chi1011) / (2 * (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta));
        fjac[16][6] = (chi[2] * chi[6] * chi[7] * chi911 - chi[0] * chi[6] * chi[7] * chi1011 + chi[4] * chi[6] * chi[7] * chi910) / (chi178 * chi910 + chi178 * chi911 + chi368 * chi910 - chi178 * chi1011 + chi368 * chi911 + chi567 * chi910 - chi368 * chi1011 + chi567 * chi911 - chi567 * chi1011);

        fjac[16][9] = (chi[1] * chi[6] * chi[7]) / (chi178 + chi368 + chi567);
        fjac[16][10] = -(chi[3] * chi[6] * chi[7]) / (chi178 + chi368 + chi567);
        fjac[16][14] = -(chi[3] * chi[6]) / (chi178 + chi368 + chi567);
        fjac[16][15] = (chi[1] * chi[7] + chi[3] * chi[6]) / (chi178 + chi368 + chi567);
        fjac[16][16] = -(chi[7] * (chi[0] * chi[6] * chi[11] + chi[2] * chi[6] * chi[9] - chi[2] * chi[6] * chi[11] - chi[4] * chi[6] * chi[9])) / (2 * (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta));
        fjac[16][17] = -(chi[0] * chi[6] * chi[7] * chi[10] - chi[2] * chi[6] * chi[7] * chi[9] + chi[4] * chi[6] * chi[7] * chi[9] - chi[4] * chi[6] * chi[7] * chi[10]) / (2 * (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta));

        fjac[17][0] = -(chi[6] * chi[7]) / (chi178 + chi368 + chi567);
        fjac[17][3] = -(chi[0] * chi[6] * chi[7] * chi910 + chi[0] * chi[6] * chi[7] * chi911 - chi[2] * chi[6] * chi[7] * chi911 - chi[4] * chi[6] * chi[7] * chi910) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);
        fjac[17][4] = -(chi[7] * chi[10] * (chi[0] * chi[6] * chi[11] + chi[2] * chi[6] * chi[9] - chi[2] * chi[6] * chi[11] - chi[4] * chi[6] * chi[9])) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);
        fjac[17][5] = -(chi[0] * chi[6] * chi[7] * chi1011 - chi[2] * chi[6] * chi[7] * chi911 + chi[4] * chi[6] * chi[7] * chi911 - chi[4] * chi[6] * chi[7] * chi1011) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);
        fjac[17][6] = (2 * (chi[2] * chi[6] * chi[7] * chi911 - chi[0] * chi[6] * chi[7] * chi1011 + chi[4] * chi[6] * chi[7] * chi910)) / (chi178 * chi910 + chi178 * chi911 + chi368 * chi910 - chi178 * chi1011 + chi368 * chi911 + chi567 * chi910 - chi368 * chi1011 + chi567 * chi911 - chi567 * chi1011);
        fjac[17][9] = (2 * chi[1] * chi[6] * chi[7]) / (chi178 + chi368 + chi567);
        fjac[17][10] = -(2 * chi[3] * chi[6] * chi[7]) / (chi178 + chi368 + chi567);
        fjac[17][11] = (2 * (chi178 + chi368)) / (chi178 + chi368 + chi567);
        fjac[17][14] = -(2 * chi[3] * chi[6]) / (chi178 + chi368 + chi567);
        fjac[17][15] = (2 * (chi[1] * chi[7] + chi[3] * chi[6])) / (chi178 + chi368 + chi567);
        fjac[17][16] = -(chi[7] * (chi[0] * chi[6] * chi[11] + chi[2] * chi[6] * chi[9] - chi[2] * chi[6] * chi[11] - chi[4] * chi[6] * chi[9])) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);
        fjac[17][17] = -(chi[0] * chi[6] * chi[7] * chi[10] - chi[2] * chi[6] * chi[7] * chi[9] + chi[4] * chi[6] * chi[7] * chi[9] - chi[4] * chi[6] * chi[7] * chi[10]) / (chi178 * chi910 * theta + chi178 * chi911 * theta + chi368 * chi910 * theta - chi178 * chi1011 * theta + chi368 * chi911 * theta + chi567 * chi910 * theta - chi368 * chi1011 * theta + chi567 * chi911 * theta - chi567 * chi1011 * theta);

        fjac[0][2] = fjac[0][1];
        fjac[0][13] = fjac[0][1];
        fjac[1][7] = fjac[9][3];
        fjac[1][8] = fjac[9][3];
        fjac[1][12] = fjac[9][3];
        fjac[1][5] = fjac[1][4];
        fjac[2][5] = fjac[2][4];
        fjac[3][15] = fjac[6][14];
        fjac[4][5] = fjac[4][3];
        fjac[4][6] = fjac[12][3];
        fjac[4][17] = fjac[7][16];
        fjac[5][5] = fjac[5][3];
        fjac[5][17] = fjac[8][16];
        fjac[7][4] = fjac[7][3];
        fjac[7][6] = fjac[15][3];
        fjac[8][4] = fjac[8][3];
        fjac[9][4] = fjac[9][3];
        fjac[9][5] = fjac[9][3];
        fjac[10][9] = fjac[0][1];
        fjac[12][4] = fjac[12][3];
        fjac[12][5] = fjac[12][3];
        fjac[12][17] = fjac[15][16];
        fjac[13][10] = fjac[3][0];
        fjac[15][4] = fjac[15][3];
        fjac[15][5] = fjac[15][3];
        fjac[16][11] = fjac[6][0];

        fjac[0][7] = -fjac[0][6];
        fjac[0][8] = -fjac[0][6];
        fjac[0][12] = -fjac[0][6];
        fjac[2][7] = -fjac[2][6];
        fjac[2][8] = -fjac[2][6];
        fjac[2][12] = -fjac[2][6];
        fjac[3][2] = -fjac[3][0];
        fjac[3][13] = -fjac[3][0];
        fjac[3][7] = -fjac[3][6];
        fjac[3][8] = -fjac[3][6];
        fjac[3][12] = -fjac[3][6];
        fjac[4][8] = -fjac[12][3];
        fjac[4][12] = -fjac[12][3];
        fjac[5][7] = -fjac[5][6];
        fjac[5][8] = -fjac[5][6];
        fjac[5][12] = -fjac[5][6];
        fjac[6][1] = -fjac[6][0];
        fjac[6][13] = -fjac[6][0];
        fjac[6][7] = -fjac[6][6];
        fjac[6][8] = -fjac[6][6];
        fjac[6][12] = -fjac[6][6];
        fjac[7][7] = -fjac[15][3];
        fjac[7][12] = -fjac[15][3];
        fjac[8][7] = -fjac[8][6];
        fjac[8][8] = -fjac[8][6];
        fjac[8][12] = -fjac[8][6];
        fjac[9][7] = -fjac[9][6];
        fjac[9][8] = -fjac[9][6];
        fjac[9][12] = -fjac[9][6];
        fjac[10][1] = -fjac[10][0];
        fjac[10][2] = -fjac[10][0];
        fjac[10][13] = -fjac[10][0];
        fjac[10][7] = -fjac[10][6];
        fjac[10][8] = -fjac[10][6];
        fjac[10][12] = -fjac[10][6];
        fjac[11][1] = -fjac[11][0];
        fjac[11][2] = -fjac[11][0];
        fjac[11][13] = -fjac[11][0];
        fjac[11][7] = -fjac[11][6];
        fjac[11][8] = -fjac[11][6];
        fjac[11][12] = -fjac[11][6];
        fjac[12][7] = -fjac[12][6];
        fjac[12][8] = -fjac[12][6];
        fjac[12][12] = -fjac[12][6];
        fjac[13][1] = -fjac[13][0];
        fjac[13][2] = -fjac[13][0];
        fjac[13][13] = -fjac[13][0];
        fjac[13][7] = -fjac[13][6];
        fjac[13][8] = -fjac[13][6];
        fjac[13][12] = -fjac[13][6];
        fjac[14][1] = -fjac[14][0];
        fjac[14][2] = -fjac[14][0];
        fjac[14][13] = -fjac[14][0];
        fjac[14][7] = -fjac[14][6];
        fjac[14][8] = -fjac[14][6];
        fjac[14][12] = -fjac[14][6];
        fjac[15][7] = -fjac[15][6];
        fjac[15][8] = -fjac[15][6];
        fjac[15][12] = -fjac[15][6];
        fjac[16][1] = -fjac[16][0];
        fjac[16][2] = -fjac[16][0];
        fjac[16][13] = -fjac[16][0];
        fjac[16][7] = -fjac[16][6];
        fjac[16][8] = -fjac[16][6];
        fjac[16][12] = -fjac[16][6];
        fjac[17][1] = -fjac[17][0];
        fjac[17][2] = -fjac[17][0];
        fjac[17][13] = -fjac[17][0];
        fjac[17][7] = -fjac[17][6];
        fjac[17][8] = -fjac[17][6];
        fjac[17][12] = -fjac[17][6];

        // Check whether solution is close enough. If not run the loop again.
        //  int ch = zero(xb, 18, 1.0e-10, 1.0e-10, fvec, fjac);
        int ch = zeroNoDcmp(xb, 18, tol, tol, fvec, fjac); //
        if (ch == 1) ok = true;
        j = j + 1;
    }

    // Solutions is applied, and right boundary is updated.
    PV -> Anew[N] = xb[9];
    PV -> Qnew[N] = xb[0];
    LD -> Anew[0] = xb[12];
    LD -> Qnew[0] = xb[3];
    RD -> Anew[0] = xb[15];
    RD -> Qnew[0] = xb[6];

//    while (j <= ntrial && ok == false) // Find the zero
//    {
//        double fvec[18];
//        // The residuals.
//
//        fvec[0] = g1 - xb[0] -
//            theta * (sq(xb[2]) / xb[11] + PV ->Bh(N, xb[11])) +
//            gamma * (PV ->F(xb[2], xb[11]) + PV ->dBdx1h(N, xb[11]));
//
//        fvec[1] = g2 - xb[3] +
//            theta * (sq(xb[5]) / xb[14] + LD -> Bh(-1, xb[14])) +
//            gamma * (PV ->F(xb[5], xb[14]) + LD -> dBdx1h(-1, xb[14]));
//
//        fvec[2] = g2a - xb[6] +
//            theta * (sq(xb[8]) / xb[17] + RD -> Bh(-1, xb[17])) +
//            gamma * (PV ->F(xb[8], xb[17]) + RD -> dBdx1h(-1, xb[17]));
//
//        fvec[3] = -theta * xb[2] - xb[9] + k1;
//        fvec[4] = theta * xb[5] - xb[12] + k2;
//        fvec[5] = theta * xb[8] - xb[15] + k2a;
//        fvec[6] = -xb[1] + xb[2] / 2.0 + k3;
//        fvec[7] = -xb[4] + xb[5] / 2.0 + k4;
//        fvec[8] = -xb[7] + xb[8] / 2.0 + k4a;
//        fvec[9] = -xb[10] + xb[11] / 2.0 + k5;
//        fvec[10] = -xb[13] + xb[14] / 2.0 + k6;
//        fvec[11] = -xb[16] + xb[17] / 2.0 + k6a;
//        fvec[12] = -xb[1] + xb[4] + xb[7];
//        fvec[13] = -xb[0] + xb[3] + xb[6];
//
//        PN = PV ->P(qLnb,N, xb[10]); // FLP
//        double sq211 = sq(xb[1] / xb[10]);
//
//        if (xb[1] > 0) {
//            fvec[14] = -PN + LD -> P(qLnb,0, xb[13]) + k7nh * sq211;
//            fvec[15] = -PN + RD -> P(qLnb,0, xb[16]) + k7anh * sq211;
//        } else {
//            fvec[14] = -PN + LD -> P(qLnb,0, xb[13]) - k7nh * sq211;
//            fvec[15] = -PN + RD -> P(qLnb,0, xb[16]) - k7anh * sq211;
//        };
//
//        PN = PV ->P(qLnb,N, xb[9]); // FLP
//        double sq110 = sq(xb[0] / xb[9]);
//        if (xb[0] > 0) {
//            fvec[16] = -PN + LD -> P(qLnb,0, xb[12]) + k7n * sq110;
//            fvec[17] = -PN + RD -> P(qLnb,0, xb[15]) + k7an * sq110;
//        } else {
//            fvec[16] = -PN + LD -> P(qLnb,0, xb[12]) - k7n * sq110;
//            fvec[17] = -PN + RD -> P(qLnb,0, xb[15]) - k7an * sq110;
//        };
//
//        // define chi. JAM. March 2020.
//        double chi[13]; // 13 to make indexing the same as Umar's thesis.
//
//        chi[1] = -2 * theta * xb[2] / xb[11] + gamma * PV ->dFdQ(xb[11]);
//
//        chi[2] = theta * (sq(xb[2] / xb[11]) - PV ->dBdAh(N, xb[11])) +
//            gamma * (PV ->dFdA(xb[2], xb[11]) + PV ->d2BdAdxh(N, xb[11]));
//        chi[3] = 2 * theta * xb[5] / xb[14] + gamma * PV ->dFdQ(xb[14]);
//
//        chi[4] = theta * (-sq(xb[5] / xb[14]) + LD -> dBdAh(-1, xb[14])) +
//            gamma * (PV ->dFdA(xb[5], xb[14]) +
//                LD -> d2BdAdxh(-1, xb[14]));
//
//        chi[5] = 2 * theta * xb[8] / xb[17] + gamma * (PV ->dFdQ(xb[17]));
//
//        chi[6] = theta * (-sq(xb[8] / xb[17]) + RD -> dBdAh(-1, xb[17])) +
//            gamma * (PV ->dFdA(xb[8], xb[17]) +
//                RD -> d2BdAdxh(-1, xb[17]));
//
//        if (xb[1] > 0) {
//            chi[7] = -(PV ->dPdA(N, xb[10])) + sq(xb[1]) / cu(xb[10]) * (-2 * k7nh);
//        } else {
//            chi[7] = -(PV ->dPdA(N, xb[10])) + sq(xb[1]) / cu(xb[10]) * (2 * k7nh);
//        };
//
//        chi[8] = LD -> dPdA(0, xb[13]);
//        chi[9] = RD -> dPdA(0, xb[16]);
//
//        if (xb[0] > 0) {
//            chi[10] = -(PV ->dPdA(N, xb[9])) + sq(xb[0]) / cu(xb[9]) * (-2 * k7an);
//        } else {
//            chi[10] = -(PV ->dPdA(N, xb[9])) + sq(xb[0]) / cu(xb[9]) * (2 * k7n);
//        };
//
//        chi[11] = LD -> dPdA(0, xb[12]);
//        chi[12] = RD -> dPdA(0, xb[15]);
//
//        // The Jacobian.
//
//        for (int row = 0; row < 18; row++)
//            for (int col = 0; col < 18; col++)
//                fjac[row][col] = 0.0;
//
//        fjac[0][0] = -1.0;
//        fjac[0][2] = chi[1];
//        fjac[0][11] = chi[2];
//
//        fjac[1][3] = -1.0;
//        fjac[1][5] = chi[3];
//        fjac[1][14] = chi[4];
//
//        fjac[2][6] = -1.0;
//        fjac[2][8] = chi[5];
//        fjac[2][17] = chi[6];
//
//        fjac[3][2] = -theta;
//        fjac[3][9] = -1.0;
//        fjac[4][5] = theta;
//        fjac[4][12] = -1.0;
//        fjac[5][8] = theta;
//        fjac[5][15] = -1.0;
//
//        fjac[6][1] = -1.0; fjac[6][2] = 0.5;
//        fjac[7][4] = -1.0; fjac[7][5] = 0.5;
//        fjac[8][7] = -1.0; fjac[8][8] = 0.5;
//        fjac[9][10] = -1.0; fjac[9][11] = 0.5;
//        fjac[10][13] = -1.0; fjac[10][14] = 0.5; fjac[11][16] = -1.0; fjac[11][17] = 0.5;
//        fjac[12][1] = -1.0; fjac[12][4] = 1.0; fjac[12][7] = 1.0; fjac[13][0] = -1.0; fjac[13][3] = 1.0; fjac[13][6] = 1.0;
//
//        fjac[14][10] = chi[7]; fjac[14][13] = chi[8];
//
//        fjac[15][10] = chi[7]; fjac[15][16] = chi[9];
//
//        fjac[16][9] = chi[10]; fjac[16][12] = chi[11];
//
//        fjac[17][9] = chi[10]; fjac[17][15] = chi[12];
//
//        if (xb[1] > 0) {
//            fjac[14][1] = xb[1] / sq(xb[10]) * (2 * k7nh);
//            fjac[15][1] = xb[1] / sq(xb[10]) * (2 * k7anh);
//        } else {
//            fjac[14][1] = xb[1] / sq(xb[10]) * (-2 * k7nh);
//            fjac[15][1] = xb[1] / sq(xb[10]) * (-2 * k7anh);
//        };
//
//        if (xb[0] > 0) {
//            fjac[16][0] = xb[0] / sq(xb[9]) * (2 * k7n);
//            fjac[17][0] = xb[0] / sq(xb[9]) * (2 * k7an);
//        } else {
//            fjac[16][0] = xb[0] / sq(xb[9]) * (-2 * k7n);
//            fjac[17][0] = xb[0] / sq(xb[9]) * (-2 * k7an);
//        };
//        // Check whether solution is close enough. If not run the loop again.
//        int ch = zero(xb, 18, tol, tol, fvec, fjac);
//        if (ch == 1) ok = true;
//
//        j = j + 1;
//    }
//
//    // Solutions is applied, and right boundary is updated.
//    PV ->Anew[N] = xb[9];
//    PV ->Qnew[N] = xb[0];
//    LD -> Anew[0] = xb[12];
//    LD -> Qnew[0] = xb[3];
//    RD -> Anew[0] = xb[15];
//    RD -> Qnew[0] = xb[6];

    if (j >= ntrial) error("arteries.C", "Root not found in the bifurcation");
    if (j >= ntrial) {printf("Error (junction.C)\n Root not found in the bifurcation\n Arteries[%i]\nqLnb = %i", parent, qLnb); exit(202);}
}

void trifurcation(int qLnb, double theta, double gamma, Tube * Arteries[], int parent, int ntrial, double tol, int j, int ok) {
    //*****************************************************************//
    //===================== TRIFURCATION MATCHING =====================//
    //*****************************************************************//
    Tube * PV = Arteries[parent];
    Tube * LD = Arteries[parent] -> LD;
    Tube * MD = Arteries[parent] -> MD;
    Tube * RD = Arteries[parent] -> RD;
    int N = PV -> N;

    double PN; // FLP. Defined otherwhere.

    //    int j = 1;
    //    int ok = false;
    double k1[4];
    double k2[4];
    double k3[4];
    double k4[4];

    // These are the flows at the half time step.
    k1[0] = PV -> Qold[N] + theta * (PV -> R2h[N - 1]) + gamma * (PV -> S2h[N - 1]);

    k1[1] = LD -> Qold[0] - theta * (LD -> R2h[0]) + gamma * (LD -> S2h[0]);

    k1[2] = MD -> Qold[0] - theta * (MD -> R2h[0]) + gamma * (MD -> S2h[0]);

    k1[3] = RD -> Qold[0] - theta * (RD -> R2h[0]) + gamma * (RD -> S2h[0]);

    // These are areas at the half time step
    k2[0] = PV -> Aold[N] + theta * (PV -> R1h[N - 1]);
    k2[1] = LD -> Aold[0] - theta * (LD -> R1h[0]);
    k2[2] = MD -> Aold[0] - theta * (MD -> R1h[0]);
    k2[3] = RD -> Aold[0] - theta * (RD -> R1h[0]);

    // These are flows at the half time stpe
    k3[0] = PV -> Qh[N - 1] / 2.0;
    k3[1] = LD -> Qh[0] / 2.0;
    k3[2] = MD -> Qh[0] / 2.0;
    k3[3] = RD -> Qh[0] / 2.0;

    k4[0] = PV -> Ah[N - 1] / 2.0;
    k4[1] = LD -> Ah[0] / 2.0;
    k4[2] = MD -> Ah[0] / 2.0;
    k4[3] = RD -> Ah[0] / 2.0;

    double xb[24];

    // The approximative initial guesses are applied.
    // Initial flows
    xb[0] = PV -> Qh[N - 1]; //Initial guess for Qp_xb n+1
    xb[1] = (PV -> Qold[N - 1] + PV -> Qold[N]) / 2.0; //Initial guess for Qp_xb^n+0.5
    xb[2] = PV -> Qold[N]; //Initial guess for Qp_xb+0.5 n+0.5

    xb[3] = LD -> Qh[0]; //Initial guess for Qd1_xb n+1
    xb[4] = (LD -> Qold[0] + LD -> Qold[1]) / 2.0; //Initial guess for Qd1_xb n+0.5
    xb[5] = LD -> Qold[0]; //Initial guess for Qd1_xb+0.5 n+0.5

    xb[6] = MD -> Qh[0]; //Initial guess for Qd2_xb n+1
    xb[7] = (MD -> Qold[0] + MD -> Qold[1]) / 2.0; //Initial guess for Qd2_xb n+0.5
    xb[8] = MD -> Qold[0]; //Initial guess for Qd2_xb+0.5 n+0.5

    xb[9] = RD -> Qh[0]; //Initial guess for Qd3_xb n+1
    xb[10] = (RD -> Qold[0] + RD -> Qold[1]) / 2.0; //Initial guess for Qd3_xb n+0.5
    xb[11] = RD -> Qold[0]; //Initial guess for Qd3_xb+0.5 n+0.5

    xb[12] = PV -> Ah[N - 1]; //Initial guess for Ap_xb n+1
    xb[13] = (PV -> Aold[N - 1] + PV -> Aold[N]) / 2.0; //Initial guess for Ap_xb^n+0.5
    xb[14] = PV -> Aold[N]; //Initial guess for Ap_xb+0.5 n+0.5

    xb[15] = LD -> Ah[0]; //Initial guess for Ad1_xb n+1
    xb[16] = (LD -> Aold[0] + LD -> Aold[1]) / 2.0; //Initial guess for Ad1_xb n+0.5
    xb[17] = LD -> Aold[0]; //Initial guess for Ad1_xb+0.5 n+0.5

    xb[18] = MD -> Ah[0]; //Initial guess for Ad2_xb n+1
    xb[19] = (MD -> Aold[0] + MD -> Aold[1]) / 2.0; //Initial guess for Ad2_xb n+0.5
    xb[20] = MD -> Aold[0]; //Initial guess for Ad2_xb+0.5 n+0.5

    xb[21] = RD -> Ah[0]; //Initial guess for Ad3_xb n+1
    xb[22] = (RD -> Aold[0] + RD -> Aold[1]) / 2.0; //Initial guess for Ad3_xb n+0.5
    xb[23] = RD -> Aold[0]; //Initial guess for Ad3_xb+0.5 n+0.5

    // This is where a Bernoulli term can be perscribed; set to zero otherwise
    double k7[4];
    double k7h[4];
    k7[0] = 0.0;
    k7[1] = 0.0;
    k7[2] = 0.0;
    k7h[0] = 0.0;
    k7h[1] = 0.0;
    k7h[2] = 0.0;

    // The residuals (fvec), and the Jacobian is determined, and if possible
    // the system of equations is solved.
    while (j <= ntrial && ok == false) // Find the zero
    {
        double fvec[24];

        // The residuals.

        // Characteristic Q residual at n+1
        fvec[0] = k1[0] - xb[0] -
            theta * (sq(xb[2]) / xb[14] + PV -> Bh(N, xb[14])) +
            gamma * (PV -> F(xb[2], xb[14]) + PV -> dBdx1h(N, xb[14]));

        fvec[1] = k1[1] - xb[3] +
            theta * (sq(xb[5]) / xb[17] + LD -> Bh(-1, xb[17])) +
            gamma * (PV -> F(xb[5], xb[17]) + LD -> dBdx1h(-1, xb[17]));

        fvec[2] = k1[2] - xb[6] +
            theta * (sq(xb[8]) / xb[20] + MD -> Bh(-1, xb[20])) +
            gamma * (PV -> F(xb[8], xb[20]) + MD -> dBdx1h(-1, xb[20]));

        fvec[3] = k1[3] - xb[9] +
            theta * (sq(xb[11]) / xb[23] + RD -> Bh(-1, xb[23])) +
            gamma * (PV -> F(xb[11], xb[23]) + RD -> dBdx1h(-1, xb[23]));

        // Characteristic A residual at n+1
        fvec[4] = -theta * xb[2] - xb[12] + k2[0];
        fvec[5] = theta * xb[5] - xb[15] + k2[1];
        fvec[6] = theta * xb[8] - xb[18] + k2[2];
        fvec[7] = theta * xb[11] - xb[21] + k2[3];

        // Flow residuals at n+1/2 (ghost points)
        fvec[8] = -xb[1] + xb[2] / 2.0 + k3[0];
        fvec[9] = -xb[4] + xb[5] / 2.0 + k3[1];
        fvec[10] = -xb[7] + xb[8] / 2.0 + k3[2];
        fvec[11] = -xb[10] + xb[11] / 2.0 + k3[3];

        // Area residuals at n+1/2 (ghost points)
        fvec[12] = -xb[13] + xb[14] / 2.0 + k4[0];
        fvec[13] = -xb[16] + xb[17] / 2.0 + k4[1];
        fvec[14] = -xb[19] + xb[20] / 2.0 + k4[2];
        fvec[15] = -xb[22] + xb[23] / 2.0 + k4[3];

        // Flow conservation residuals (n+1/2 and n+1)
        fvec[16] = -xb[1] + xb[4] + xb[7] + xb[10];
        fvec[17] = -xb[0] + xb[3] + xb[6] + xb[9];

        // Pressure continuity at the n+1/2 time step
        PN = PV -> P(qLnb, N, xb[13]);
        double u_n_half = sq(xb[1] / xb[13]);

        // The if statements here only matter if a minor loss
        // is included.
        fvec[18] = -PN + LD -> P(qLnb, 0, xb[16]) + abs(k7h[0]) * u_n_half;
        fvec[19] = -PN + MD -> P(qLnb, 0, xb[19]) + abs(k7h[1]) * u_n_half;
        fvec[20] = -PN + RD -> P(qLnb, 0, xb[22]) + abs(k7h[2]) * u_n_half;
        // Pressure continuity at the n+1 time step
        PN = PV -> P(qLnb, N, xb[12]);
        double u_n_1 = sq(xb[0] / xb[12]);
        fvec[21] = -PN + LD -> P(qLnb, 0, xb[15]) + abs(k7[0]) * u_n_1;
        fvec[22] = -PN + MD -> P(qLnb, 0, xb[18]) + abs(k7[1]) * u_n_1;
        fvec[23] = -PN + RD -> P(qLnb, 0, xb[21]) + abs(k7[2]) * u_n_1;

        for (int row = 0; row < 4 * 6; row++)
            for (int col = 0; col < 4 * 6; col++)
                ftrif[row][col] = 0.0;

        //        // The Jacobian.
        double chi[16];

        // Here are the residuals for the characteristic matching for flow
        chi[0] = -2.0 * theta * xb[2] / xb[14] + gamma * PV -> dFdQ(xb[14]);
        chi[2] = 2.0 * theta * xb[5] / xb[17] + gamma * PV -> dFdQ(xb[17]);
        chi[4] = 2.0 * theta * xb[8] / xb[20] + gamma * PV -> dFdQ(xb[20]);
        chi[6] = 2.0 * theta * xb[11] / xb[23] + gamma * PV -> dFdQ(xb[23]);

        // Here are the residuals for the area characteristic matching
        chi[1] = theta * (sq(xb[2] / xb[14]) - PV -> dBdAh(N, xb[14])) +
            gamma * (PV -> dFdA(xb[2], xb[14]) + PV -> d2BdAdxh(N, xb[14]));

        chi[3] = theta * (-sq(xb[5] / xb[17]) + LD -> dBdAh(-1, xb[17])) +
            gamma * (PV -> dFdA(xb[5], xb[17]) + LD -> d2BdAdxh(-1, xb[17]));

        chi[5] = theta * (-sq(xb[8] / xb[20]) + MD -> dBdAh(-1, xb[20])) +
            gamma * (PV -> dFdA(xb[8], xb[20]) + MD -> d2BdAdxh(-1, xb[20]));

        chi[7] = theta * (-sq(xb[11] / xb[23]) + RD -> dBdAh(-1, xb[23])) +
            gamma * (PV -> dFdA(xb[11], xb[23]) + RD -> d2BdAdxh(-1, xb[23]));

        // Here is pressure conservation (n+1/2)
        chi[8] = -(PV -> dPdA(N, xb[13])) + sq(xb[1]) / cu(xb[13]) * (-2.0 * abs(k7h[0])); //Loss term
        chi[9] = LD -> dPdA(0, xb[16]);
        chi[10] = MD -> dPdA(0, xb[19]);
        chi[11] = RD -> dPdA(0, xb[22]);

        // Here is pressure conservation (n+1)
        chi[12] = -(PV -> dPdA(N, xb[12])) + sq(xb[1]) / cu(xb[12]) * (-2.0 * abs(k7h[0])); //Loss term
        chi[13] = LD -> dPdA(0, xb[15]);
        chi[14] = MD -> dPdA(0, xb[18]);
        chi[15] = RD -> dPdA(0, xb[21]);

        // The jacobian
        // Order is [row][column]

//        check determinant
        double detFtrif = (chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]) - (chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]) - (chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]) - (chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]) - (chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]) - (chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]) - (chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]) - (chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]) - (chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]) - (chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]) + (chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]) - (chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]) - (chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]) + (chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]) - (chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]) + (chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);

        if (sqrt(sq(detFtrif)) < tol) {
            error("arteries.c", "Singular Trifurcation Jacobian");
            exit(111);
        }

        ftrif[0][0] = -(chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[0][1] = -(chi[1]*chi[9]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[0][2] = -(chi[1]*chi[9]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[0][3] = -(chi[1]*chi[9]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[0][4] = -(chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] + chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] + chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] + chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[0][5] = -(chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] - chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] - chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] + chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[0][6] = -(chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] - chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[0][7] = -(chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] - chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] - chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[0][8] = -(2*(chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] + chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[0][9] = (2*(chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] + chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[0][10] = (2*(chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] + chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[0][11] = (2*(chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] + chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[0][12] = (2*(chi[1]*chi[3]*chi[8]*chi[10]*chi[11] + chi[1]*chi[5]*chi[8]*chi[9]*chi[11] + chi[1]*chi[7]*chi[8]*chi[9]*chi[10]))/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[0][13] = (2*chi[1]*chi[3]*chi[9]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[0][14] = (2*chi[1]*chi[5]*chi[9]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[0][15] = (2*chi[1]*chi[7]*chi[9]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[0][16] = (2*(chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] + chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[0][17] = -(chi[1]*chi[9]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[0][18] = (2*chi[1]*chi[3]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[0][19] = (2*chi[1]*chi[5]*chi[9]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[0][20] = (2*chi[1]*chi[7]*chi[9]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[0][21] = -(chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[14]*chi[15] - chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[15] - chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14] + chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[14]*chi[15] + chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[14]*chi[15] + chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[15] + chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[14]*chi[15] + chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[0][22] = -(chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[15] + chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[15] - chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13] + chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[15] + chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[15] - chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[15] + chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13] + chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[13]*chi[15])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[0][23] = -(chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14] + chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14] + chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14] + chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13] + chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14] - chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13] - chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14] + chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);


        ftrif[1][0] = 0;
        ftrif[1][1] = 0;
        ftrif[1][2] = 0;
        ftrif[1][3] = 0;
        ftrif[1][4] = -(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));
        ftrif[1][5] = -(chi[13]*chi[14]*chi[15])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));
        ftrif[1][6] = -(chi[13]*chi[14]*chi[15])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));
        ftrif[1][7] = -(chi[13]*chi[14]*chi[15])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));
        ftrif[1][8] = -(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[1][9] = (chi[13]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[1][10] = (chi[13]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[1][11] = (chi[13]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[1][12] = 0;
        ftrif[1][13] = 0;
        ftrif[1][14] = 0;
        ftrif[1][15] = 0;
        ftrif[1][16] = (chi[13]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[1][17] = 0;
        ftrif[1][18] = 0;
        ftrif[1][19] = 0;
        ftrif[1][20] = 0;
        ftrif[1][21] = -(chi[14]*chi[15])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));
        ftrif[1][22] = -(chi[13]*chi[15])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));
        ftrif[1][23] = -(chi[13]*chi[14])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));


        ftrif[2][0] = 0;
        ftrif[2][1] = 0;
        ftrif[2][2] = 0;
        ftrif[2][3] = 0;
        ftrif[2][4] = -(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);
        ftrif[2][5] = -(chi[13]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);
        ftrif[2][6] = -(chi[13]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);
        ftrif[2][7] = -(chi[13]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);
        ftrif[2][8] = -(2*chi[13]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[2][9] = (2*chi[13]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[2][10] = (2*chi[13]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[2][11] = (2*chi[13]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[2][12] = 0;
        ftrif[2][13] = 0;
        ftrif[2][14] = 0;
        ftrif[2][15] = 0;
        ftrif[2][16] = (2*chi[13]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[2][17] = 0;
        ftrif[2][18] = 0;
        ftrif[2][19] = 0;
        ftrif[2][20] = 0;
        ftrif[2][21] = -(chi[14]*chi[15])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);
        ftrif[2][22] = -(chi[13]*chi[15])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);
        ftrif[2][23] = -(chi[13]*chi[14])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);


        ftrif[3][0] = -(chi[3]*chi[8]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[3][1] = -(chi[1]*chi[9]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[3][2] = (chi[3]*chi[8]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[3][3] = (chi[3]*chi[8]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[3][4] = -(chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] - chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[3][5] = (chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] - chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] + chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] - chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] - chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[3][6] = -(chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] + chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] - chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[3][7] = -(chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] - chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] + chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[3][8] = -(2*(chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] - chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[3][9] = (2*(chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] - chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[3][10] = (2*(chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] - chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[3][11] = (2*(chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] - chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[3][12] = (2*chi[1]*chi[3]*chi[8]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[3][13] = (2*(chi[1]*chi[3]*chi[9]*chi[10]*chi[11] + chi[3]*chi[5]*chi[8]*chi[9]*chi[11] + chi[3]*chi[7]*chi[8]*chi[9]*chi[10]))/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[3][14] = -(2*chi[3]*chi[5]*chi[8]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[3][15] = -(2*chi[3]*chi[7]*chi[8]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[3][16] = (2*(chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] - chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[3][17] = (chi[3]*chi[8]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[3][18] = (2*(chi[1]*chi[3]*chi[10]*chi[11] + chi[3]*chi[5]*chi[8]*chi[11] + chi[3]*chi[7]*chi[8]*chi[10]))/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[3][19] = -(2*chi[3]*chi[5]*chi[8]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[3][20] = -(2*chi[3]*chi[7]*chi[8]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[3][21] = (chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14] + chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[15] - chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[14]*chi[15] + chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14] - chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[14]*chi[15] + chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[15] + chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14] + chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[15] + chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[15] - chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[14]*chi[15] + chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14] - chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[14]*chi[15])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[3][22] = -(chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[15] + chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[15] + chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13] + chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[15] + chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[15] + chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[15] - chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13] - chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[13]*chi[15])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[3][23] = -(chi[0]*chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14] + chi[1]*chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14] + chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14] - chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13] + chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14] + chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13] + chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14] - chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);


        ftrif[4][0] = 0;
        ftrif[4][1] = 0;
        ftrif[4][2] = 0;
        ftrif[4][3] = 0;
        ftrif[4][4] = -(chi[12]*chi[14]*chi[15])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));
        ftrif[4][5] = (chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] - chi[13]*chi[14]*chi[15])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));
        ftrif[4][6] = -(chi[12]*chi[14]*chi[15])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));
        ftrif[4][7] = -(chi[12]*chi[14]*chi[15])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));
        ftrif[4][8] = -(chi[12]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[4][9] = -(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] - chi[13]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[4][10] = (chi[12]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[4][11] = (chi[12]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[4][12] = 0;
        ftrif[4][13] = 0;
        ftrif[4][14] = 0;
        ftrif[4][15] = 0;
        ftrif[4][16] = (chi[12]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[4][17] = 0;
        ftrif[4][18] = 0;
        ftrif[4][19] = 0;
        ftrif[4][20] = 0;
        ftrif[4][21] = (chi[12]*chi[14] + chi[12]*chi[15] - chi[14]*chi[15])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));
        ftrif[4][22] = -(chi[12]*chi[15])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));
        ftrif[4][23] = -(chi[12]*chi[14])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));


        ftrif[5][0] = 0;
        ftrif[5][1] = 0;
        ftrif[5][2] = 0;
        ftrif[5][3] = 0;
        ftrif[5][4] = -(chi[12]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);
        ftrif[5][5] = (chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] - chi[13]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);
        ftrif[5][6] = -(chi[12]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);
        ftrif[5][7] = -(chi[12]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);
        ftrif[5][8] = -(2*chi[12]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[5][9] = (2*chi[12]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[5][10] = (2*chi[12]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[5][11] = (2*chi[12]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[5][12] = 0;
        ftrif[5][13] = 0;
        ftrif[5][14] = 0;
        ftrif[5][15] = 0;
        ftrif[5][16] = (2*chi[12]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[5][17] = 0;
        ftrif[5][18] = 0;
        ftrif[5][19] = 0;
        ftrif[5][20] = 0;
        ftrif[5][21] = (chi[12]*chi[14] + chi[12]*chi[15] - chi[14]*chi[15])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);
        ftrif[5][22] = -(chi[12]*chi[15])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);
        ftrif[5][23] = -(chi[12]*chi[14])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);


        ftrif[6][0] = -(chi[5]*chi[8]*chi[9]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[6][1] = (chi[5]*chi[8]*chi[9]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[6][2] = -(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[6][3] = (chi[5]*chi[8]*chi[9]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[6][4] = -(chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] + chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] - chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[6][5] = -(chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] + chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] - chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[6][6] = (chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] - chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] - chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] + chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[6][7] = -(chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] - chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] - chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[6][8] = -(2*(chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] - chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[6][9] = (2*(chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] - chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[6][10] = (2*(chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] - chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[6][11] = (2*(chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] - chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[6][12] = (2*chi[1]*chi[5]*chi[8]*chi[9]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[6][13] = -(2*chi[3]*chi[5]*chi[8]*chi[9]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[6][14] = (2*(chi[1]*chi[5]*chi[9]*chi[10]*chi[11] + chi[3]*chi[5]*chi[8]*chi[10]*chi[11] + chi[5]*chi[7]*chi[8]*chi[9]*chi[10]))/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[6][15] = -(2*chi[5]*chi[7]*chi[8]*chi[9]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[6][16] = (2*(chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] - chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[6][17] = (chi[5]*chi[8]*chi[9]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[6][18] = -(2*chi[3]*chi[5]*chi[8]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[6][19] = (2*(chi[1]*chi[5]*chi[9]*chi[11] + chi[3]*chi[5]*chi[8]*chi[11] + chi[5]*chi[7]*chi[8]*chi[9]))/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[6][20] = -(2*chi[5]*chi[7]*chi[8]*chi[9])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[6][21] = -(chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14] + chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[14]*chi[15] + chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[15] + chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[15] + chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[15] - chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[14]*chi[15] + chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[15] - chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[6][22] = (chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13] - chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[15] + chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13] + chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[15] + chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[15] - chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[13]*chi[15] + chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[15] + chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13] - chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[13]*chi[15] + chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13] + chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[15] - chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[15])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[6][23] = -(chi[0]*chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14] + chi[1]*chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13] - chi[2]*chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14] + chi[3]*chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13] + chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13] + chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13] + chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14] - chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);


        ftrif[7][0] = 0;
        ftrif[7][1] = 0;
        ftrif[7][2] = 0;
        ftrif[7][3] = 0;
        ftrif[7][4] = -(chi[12]*chi[13]*chi[15])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));
        ftrif[7][5] = -(chi[12]*chi[13]*chi[15])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));
        ftrif[7][6] = (chi[12]*chi[13]*chi[14] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));
        ftrif[7][7] = -(chi[12]*chi[13]*chi[15])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));
        ftrif[7][8] = -(chi[12]*chi[13]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[7][9] = (chi[12]*chi[13]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[7][10] = -(chi[12]*chi[13]*chi[14] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[7][11] = (chi[12]*chi[13]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[7][12] = 0;
        ftrif[7][13] = 0;
        ftrif[7][14] = 0;
        ftrif[7][15] = 0;
        ftrif[7][16] = (chi[12]*chi[13]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[7][17] = 0;
        ftrif[7][18] = 0;
        ftrif[7][19] = 0;
        ftrif[7][20] = 0;
        ftrif[7][21] = -(chi[12]*chi[15])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));
        ftrif[7][22] = (chi[12]*chi[13] + chi[12]*chi[15] - chi[13]*chi[15])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));
        ftrif[7][23] = -(chi[12]*chi[13])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));


        ftrif[8][0] = 0;
        ftrif[8][1] = 0;
        ftrif[8][2] = 0;
        ftrif[8][3] = 0;
        ftrif[8][4] = -(chi[12]*chi[13]*chi[15])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);
        ftrif[8][5] = -(chi[12]*chi[13]*chi[15])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);
        ftrif[8][6] = (chi[12]*chi[13]*chi[14] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);
        ftrif[8][7] = -(chi[12]*chi[13]*chi[15])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);
        ftrif[8][8] = -(2*chi[12]*chi[13]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[8][9] = (2*chi[12]*chi[13]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[8][10] = (2*chi[12]*chi[13]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[8][11] = (2*chi[12]*chi[13]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[8][12] = 0;
        ftrif[8][13] = 0;
        ftrif[8][14] = 0;
        ftrif[8][15] = 0;
        ftrif[8][16] = (2*chi[12]*chi[13]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[8][17] = 0;
        ftrif[8][18] = 0;
        ftrif[8][19] = 0;
        ftrif[8][20] = 0;
        ftrif[8][21] = -(chi[12]*chi[15])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);
        ftrif[8][22] = (chi[12]*chi[13] + chi[12]*chi[15] - chi[13]*chi[15])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);
        ftrif[8][23] = -(chi[12]*chi[13])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);


        ftrif[9][0] = -(chi[7]*chi[8]*chi[9]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[9][1] = (chi[7]*chi[8]*chi[9]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[9][2] = (chi[7]*chi[8]*chi[9]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[9][3] = -(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[9][4] = -(chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] + chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] + chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] - chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] + chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] - chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[9][5] = -(chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] + chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] + chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] - chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] - chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[9][6] = -(chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] + chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] - chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] + chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] + chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] + chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[9][7] = (chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] - chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] + chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] + chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] - chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] - chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[9][8] = -(2*(chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] + chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] - chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] + chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] - chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[9][9] = (2*(chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] + chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] - chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] + chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] - chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[9][10] = (2*(chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] + chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] - chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] + chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] - chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[9][11] = (2*(chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] + chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] - chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] + chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] - chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[9][12] = (2*chi[1]*chi[7]*chi[8]*chi[9]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[9][13] = -(2*chi[3]*chi[7]*chi[8]*chi[9]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[9][14] = -(2*chi[5]*chi[7]*chi[8]*chi[9]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[9][15] = (2*chi[7]*(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11]))/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[9][16] = (2*(chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] + chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] - chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] + chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] - chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[9][17] = (chi[7]*chi[8]*chi[9]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[9][18] = -(2*chi[3]*chi[7]*chi[8]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[9][19] = -(2*chi[5]*chi[7]*chi[8]*chi[9])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[9][20] = (2*chi[7]*(chi[1]*chi[9]*chi[10] + chi[3]*chi[8]*chi[10] + chi[5]*chi[8]*chi[9]))/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[9][21] = -(chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14] + chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[14]*chi[15] + chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14] + chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[15] + chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14] - chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[14]*chi[15] - chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[15] + chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[9][22] = -(chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[15] + chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13] - chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[15] + chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13] + chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13] + chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13] + chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[15] - chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[15])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[9][23] = (chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13] - chi[0]*chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14] + chi[2]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14] + chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14] + chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13] + chi[4]*chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13] - chi[1]*chi[6]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14] + chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14] + chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13] - chi[3]*chi[6]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14] + chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14] - chi[5]*chi[6]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);


        ftrif[10][0] = 0;
        ftrif[10][1] = 0;
        ftrif[10][2] = 0;
        ftrif[10][3] = 0;
        ftrif[10][4] = -(chi[12]*chi[13]*chi[14])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));
        ftrif[10][5] = -(chi[12]*chi[13]*chi[14])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));
        ftrif[10][6] = -(chi[12]*chi[13]*chi[14])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));
        ftrif[10][7] = (chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));
        ftrif[10][8] = -(chi[12]*chi[13]*chi[14])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[10][9] = (chi[12]*chi[13]*chi[14])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[10][10] = (chi[12]*chi[13]*chi[14])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[10][11] = -(chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[10][12] = 0;
        ftrif[10][13] = 0;
        ftrif[10][14] = 0;
        ftrif[10][15] = 0;
        ftrif[10][16] = (chi[12]*chi[13]*chi[14])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[10][17] = 0;
        ftrif[10][18] = 0;
        ftrif[10][19] = 0;
        ftrif[10][20] = 0;
        ftrif[10][21] = -(chi[12]*chi[14])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));
        ftrif[10][22] = -(chi[12]*chi[13])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));
        ftrif[10][23] = (chi[12]*chi[13] + chi[12]*chi[14] - chi[13]*chi[14])/(2*(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta));


        ftrif[11][0] = 0;
        ftrif[11][1] = 0;
        ftrif[11][2] = 0;
        ftrif[11][3] = 0;
        ftrif[11][4] = -(chi[12]*chi[13]*chi[14])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);
        ftrif[11][5] = -(chi[12]*chi[13]*chi[14])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);
        ftrif[11][6] = -(chi[12]*chi[13]*chi[14])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);
        ftrif[11][7] = (chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);
        ftrif[11][8] = -(2*chi[12]*chi[13]*chi[14])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[11][9] = (2*chi[12]*chi[13]*chi[14])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[11][10] = (2*chi[12]*chi[13]*chi[14])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[11][11] = (2*chi[12]*chi[13]*chi[14])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[11][12] = 0;
        ftrif[11][13] = 0;
        ftrif[11][14] = 0;
        ftrif[11][15] = 0;
        ftrif[11][16] = (2*chi[12]*chi[13]*chi[14])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[11][17] = 0;
        ftrif[11][18] = 0;
        ftrif[11][19] = 0;
        ftrif[11][20] = 0;
        ftrif[11][21] = -(chi[12]*chi[14])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);
        ftrif[11][22] = -(chi[12]*chi[13])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);
        ftrif[11][23] = (chi[12]*chi[13] + chi[12]*chi[14] - chi[13]*chi[14])/(chi[12]*chi[13]*chi[14]*theta + chi[12]*chi[13]*chi[15]*theta + chi[12]*chi[14]*chi[15]*theta - chi[13]*chi[14]*chi[15]*theta);


        ftrif[12][0] = 0;
        ftrif[12][1] = 0;
        ftrif[12][2] = 0;
        ftrif[12][3] = 0;
        ftrif[12][4] = (chi[13]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[12][5] = (chi[13]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[12][6] = (chi[13]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[12][7] = (chi[13]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[12][8] = (2*chi[13]*chi[14]*chi[15]*theta)/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[12][9] = -(2*chi[13]*chi[14]*chi[15]*theta)/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[12][10] = -(2*chi[13]*chi[14]*chi[15]*theta)/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[12][11] = -(2*chi[13]*chi[14]*chi[15]*theta)/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[12][12] = 0;
        ftrif[12][13] = 0;
        ftrif[12][14] = 0;
        ftrif[12][15] = 0;
        ftrif[12][16] = -(2*chi[13]*chi[14]*chi[15]*theta)/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[12][17] = 0;
        ftrif[12][18] = 0;
        ftrif[12][19] = 0;
        ftrif[12][20] = 0;
        ftrif[12][21] = (chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[12][22] = (chi[13]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[12][23] = (chi[13]*chi[14])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);


        ftrif[13][0] = (chi[9]*chi[10]*chi[11])/(2*(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]));
        ftrif[13][1] = -(chi[9]*chi[10]*chi[11])/(2*(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]));
        ftrif[13][2] = -(chi[9]*chi[10]*chi[11])/(2*(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]));
        ftrif[13][3] = -(chi[9]*chi[10]*chi[11])/(2*(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]));
        ftrif[13][4] = (chi[0]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[0]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[0]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] - chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14])/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));
        ftrif[13][5] = (chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[0]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] - chi[2]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] - chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] - chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14])/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));
        ftrif[13][6] = (chi[0]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] - chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] - chi[4]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15])/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));
        ftrif[13][7] = (chi[0]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] - chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[6]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15])/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));
        ftrif[13][8] = -(chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[0]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[13][9] = (chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[0]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[13][10] = (chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[0]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[13][11] = (chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[0]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[13][12] = -(chi[1]*chi[9]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[13][13] = (chi[3]*chi[9]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[13][14] = (chi[5]*chi[9]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[13][15] = (chi[7]*chi[9]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[13][16] = (chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[0]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[13][17] = -(chi[9]*chi[10]*chi[11])/(2*(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]));
        ftrif[13][18] = (chi[3]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[13][19] = (chi[5]*chi[9]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[13][20] = (chi[7]*chi[9]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[13][21] = (chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14] + chi[0]*chi[9]*chi[10]*chi[11]*chi[14]*chi[15] + chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[15] - chi[2]*chi[9]*chi[10]*chi[11]*chi[14]*chi[15] - chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[15] - chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14])/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));
        ftrif[13][22] = (chi[0]*chi[9]*chi[10]*chi[11]*chi[13]*chi[15] - chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[15] + chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13] + chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[15] - chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13] - chi[4]*chi[9]*chi[10]*chi[11]*chi[13]*chi[15])/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));
        ftrif[13][23] = (chi[0]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14] - chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14] - chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13] + chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13] + chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14] - chi[6]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14])/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));


        ftrif[14][0] = (chi[9]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[14][1] = -(chi[9]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[14][2] = -(chi[9]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[14][3] = -(chi[9]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[14][4] = (chi[0]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[0]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[0]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] - chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[14][5] = (chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[0]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] - chi[2]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] - chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] - chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[14][6] = (chi[0]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] - chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] - chi[4]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[14][7] = (chi[0]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] - chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[6]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[14][8] = -(2*(chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[0]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[14][9] = (2*(chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[0]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[14][10] = (2*(chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[0]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[14][11] = (2*(chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[0]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[14][12] = (2*(chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]))/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[14][13] = (2*chi[3]*chi[9]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[14][14] = (2*chi[5]*chi[9]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[14][15] = (2*chi[7]*chi[9]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[14][16] = (2*(chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[0]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[14][17] = -(chi[9]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[14][18] = (2*chi[3]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[14][19] = (2*chi[5]*chi[9]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[14][20] = (2*chi[7]*chi[9]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[14][21] = (chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14] + chi[0]*chi[9]*chi[10]*chi[11]*chi[14]*chi[15] + chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[15] - chi[2]*chi[9]*chi[10]*chi[11]*chi[14]*chi[15] - chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[15] - chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[14][22] = (chi[0]*chi[9]*chi[10]*chi[11]*chi[13]*chi[15] - chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[15] + chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13] + chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[15] - chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13] - chi[4]*chi[9]*chi[10]*chi[11]*chi[13]*chi[15])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[14][23] = (chi[0]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14] - chi[2]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14] - chi[4]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13] + chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13] + chi[6]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14] - chi[6]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);


        ftrif[15][0] = 0;
        ftrif[15][1] = 0;
        ftrif[15][2] = 0;
        ftrif[15][3] = 0;
        ftrif[15][4] = -(chi[12]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[15][5] = -(chi[12]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[15][6] = -(chi[12]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[15][7] = -(chi[12]*chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[15][8] = -(2*chi[12]*chi[14]*chi[15]*theta)/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[15][9] = (2*chi[12]*chi[14]*chi[15]*theta)/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[15][10] = (2*chi[12]*chi[14]*chi[15]*theta)/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[15][11] = (2*chi[12]*chi[14]*chi[15]*theta)/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[15][12] = 0;
        ftrif[15][13] = 0;
        ftrif[15][14] = 0;
        ftrif[15][15] = 0;
        ftrif[15][16] = (2*chi[12]*chi[14]*chi[15]*theta)/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[15][17] = 0;
        ftrif[15][18] = 0;
        ftrif[15][19] = 0;
        ftrif[15][20] = 0;
        ftrif[15][21] = (chi[12]*chi[14] + chi[12]*chi[15] - chi[14]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[15][22] = -(chi[12]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[15][23] = -(chi[12]*chi[14])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);


        ftrif[16][0] = -(chi[8]*chi[10]*chi[11])/(2*(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]));
        ftrif[16][1] = (chi[8]*chi[10]*chi[11])/(2*(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]));
        ftrif[16][2] = (chi[8]*chi[10]*chi[11])/(2*(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]));
        ftrif[16][3] = (chi[8]*chi[10]*chi[11])/(2*(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]));
        ftrif[16][4] = -(chi[0]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[0]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[0]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] - chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14])/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));
        ftrif[16][5] = -(chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[0]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] - chi[2]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] - chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] - chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14])/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));
        ftrif[16][6] = -(chi[0]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] - chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] - chi[4]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15])/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));
        ftrif[16][7] = -(chi[0]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] - chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[6]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15])/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));
        ftrif[16][8] = (chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[16][9] = -(chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[16][10] = -(chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[16][11] = -(chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[16][12] = (chi[1]*chi[8]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[16][13] = -(chi[3]*chi[8]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[16][14] = -(chi[5]*chi[8]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[16][15] = -(chi[7]*chi[8]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[16][16] = -(chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[16][17] = (chi[8]*chi[10]*chi[11])/(2*(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]));
        ftrif[16][18] = (chi[1]*chi[10]*chi[11] + chi[5]*chi[8]*chi[11] + chi[7]*chi[8]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[16][19] = -(chi[5]*chi[8]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[16][20] = -(chi[7]*chi[8]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[16][21] = -(chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14] + chi[0]*chi[8]*chi[10]*chi[11]*chi[14]*chi[15] + chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[15] - chi[2]*chi[8]*chi[10]*chi[11]*chi[14]*chi[15] - chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[15] - chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14])/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));
        ftrif[16][22] = -(chi[0]*chi[8]*chi[10]*chi[11]*chi[13]*chi[15] - chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[15] + chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13] + chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[15] - chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13] - chi[4]*chi[8]*chi[10]*chi[11]*chi[13]*chi[15])/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));
        ftrif[16][23] = -(chi[0]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14] - chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14] - chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13] + chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13] + chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14] - chi[6]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14])/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));


        ftrif[17][0] = -(chi[8]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[17][1] = (chi[8]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[17][2] = (chi[8]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[17][3] = (chi[8]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[17][4] = -(chi[0]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[0]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[0]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] - chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[17][5] = -(chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[0]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] - chi[2]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] - chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] - chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[17][6] = -(chi[0]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] - chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] - chi[4]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[17][7] = -(chi[0]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] - chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[6]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[17][8] = (2*(chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[17][9] = -(2*(chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[17][10] = -(2*(chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[17][11] = -(2*(chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[17][12] = (2*chi[1]*chi[8]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[17][13] = (2*(chi[1]*chi[9]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]))/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[17][14] = -(2*chi[5]*chi[8]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[17][15] = -(2*chi[7]*chi[8]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[17][16] = -(2*(chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[17][17] = (chi[8]*chi[10]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[17][18] = (2*(chi[1]*chi[10]*chi[11] + chi[5]*chi[8]*chi[11] + chi[7]*chi[8]*chi[10]))/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[17][19] = -(2*chi[5]*chi[8]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[17][20] = -(2*chi[7]*chi[8]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[17][21] = -(chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14] + chi[0]*chi[8]*chi[10]*chi[11]*chi[14]*chi[15] + chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[15] - chi[2]*chi[8]*chi[10]*chi[11]*chi[14]*chi[15] - chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[15] - chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[17][22] = -(chi[0]*chi[8]*chi[10]*chi[11]*chi[13]*chi[15] - chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[15] + chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13] + chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[15] - chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13] - chi[4]*chi[8]*chi[10]*chi[11]*chi[13]*chi[15])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[17][23] = -(chi[0]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14] - chi[2]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14] - chi[4]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13] + chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13] + chi[6]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14] - chi[6]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);


        ftrif[18][0] = 0;
        ftrif[18][1] = 0;
        ftrif[18][2] = 0;
        ftrif[18][3] = 0;
        ftrif[18][4] = -(chi[12]*chi[13]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[18][5] = -(chi[12]*chi[13]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[18][6] = -(chi[12]*chi[13]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[18][7] = -(chi[12]*chi[13]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[18][8] = -(2*chi[12]*chi[13]*chi[15]*theta)/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[18][9] = (2*chi[12]*chi[13]*chi[15]*theta)/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[18][10] = (2*chi[12]*chi[13]*chi[15]*theta)/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[18][11] = (2*chi[12]*chi[13]*chi[15]*theta)/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[18][12] = 0;
        ftrif[18][13] = 0;
        ftrif[18][14] = 0;
        ftrif[18][15] = 0;
        ftrif[18][16] = (2*chi[12]*chi[13]*chi[15]*theta)/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[18][17] = 0;
        ftrif[18][18] = 0;
        ftrif[18][19] = 0;
        ftrif[18][20] = 0;
        ftrif[18][21] = -(chi[12]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[18][22] = (chi[12]*chi[13] + chi[12]*chi[15] - chi[13]*chi[15])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[18][23] = -(chi[12]*chi[13])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);


        ftrif[19][0] = -(chi[8]*chi[9]*chi[11])/(2*(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]));
        ftrif[19][1] = (chi[8]*chi[9]*chi[11])/(2*(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]));
        ftrif[19][2] = (chi[8]*chi[9]*chi[11])/(2*(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]));
        ftrif[19][3] = (chi[8]*chi[9]*chi[11])/(2*(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]));
        ftrif[19][4] = -(chi[11]*(chi[0]*chi[8]*chi[9]*chi[12]*chi[13]*chi[14] + chi[0]*chi[8]*chi[9]*chi[12]*chi[13]*chi[15] + chi[0]*chi[8]*chi[9]*chi[12]*chi[14]*chi[15] - chi[2]*chi[8]*chi[9]*chi[12]*chi[14]*chi[15] - chi[4]*chi[8]*chi[9]*chi[12]*chi[13]*chi[15] - chi[6]*chi[8]*chi[9]*chi[12]*chi[13]*chi[14]))/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));
        ftrif[19][5] = -(chi[11]*chi[13]*(chi[2]*chi[8]*chi[9]*chi[12]*chi[14] + chi[0]*chi[8]*chi[9]*chi[14]*chi[15] + chi[2]*chi[8]*chi[9]*chi[12]*chi[15] - chi[2]*chi[8]*chi[9]*chi[14]*chi[15] - chi[4]*chi[8]*chi[9]*chi[12]*chi[15] - chi[6]*chi[8]*chi[9]*chi[12]*chi[14]))/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));
        ftrif[19][6] = -(chi[11]*(chi[0]*chi[8]*chi[9]*chi[13]*chi[14]*chi[15] - chi[2]*chi[8]*chi[9]*chi[12]*chi[14]*chi[15] + chi[4]*chi[8]*chi[9]*chi[12]*chi[13]*chi[14] + chi[4]*chi[8]*chi[9]*chi[12]*chi[14]*chi[15] - chi[6]*chi[8]*chi[9]*chi[12]*chi[13]*chi[14] - chi[4]*chi[8]*chi[9]*chi[13]*chi[14]*chi[15]))/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));
        ftrif[19][7] = -(chi[11]*(chi[0]*chi[8]*chi[9]*chi[13]*chi[14]*chi[15] - chi[2]*chi[8]*chi[9]*chi[12]*chi[14]*chi[15] - chi[4]*chi[8]*chi[9]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[12]*chi[14]*chi[15] - chi[6]*chi[8]*chi[9]*chi[13]*chi[14]*chi[15]))/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));
        ftrif[19][8] = (chi[11]*(chi[2]*chi[8]*chi[9]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[9]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[9]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[19][9] = -(chi[11]*(chi[2]*chi[8]*chi[9]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[9]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[9]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[19][10] = -(chi[11]*(chi[2]*chi[8]*chi[9]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[9]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[9]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[19][11] = -(chi[11]*(chi[2]*chi[8]*chi[9]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[9]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[9]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[19][12] = (chi[1]*chi[8]*chi[9]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[19][13] = -(chi[3]*chi[8]*chi[9]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[19][14] = -(chi[5]*chi[8]*chi[9]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[19][15] = -(chi[7]*chi[8]*chi[9]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[19][16] = -(chi[11]*(chi[2]*chi[8]*chi[9]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[9]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[9]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[19][17] = (chi[8]*chi[9]*chi[11])/(2*(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]));
        ftrif[19][18] = -(chi[3]*chi[8]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[19][19] = (chi[1]*chi[9]*chi[11] + chi[3]*chi[8]*chi[11] + chi[7]*chi[8]*chi[9])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[19][20] = -(chi[7]*chi[8]*chi[9])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[19][21] = -(chi[11]*(chi[2]*chi[8]*chi[9]*chi[12]*chi[14] + chi[0]*chi[8]*chi[9]*chi[14]*chi[15] + chi[2]*chi[8]*chi[9]*chi[12]*chi[15] - chi[2]*chi[8]*chi[9]*chi[14]*chi[15] - chi[4]*chi[8]*chi[9]*chi[12]*chi[15] - chi[6]*chi[8]*chi[9]*chi[12]*chi[14]))/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));
        ftrif[19][22] = -(chi[11]*(chi[0]*chi[8]*chi[9]*chi[13]*chi[15] - chi[2]*chi[8]*chi[9]*chi[12]*chi[15] + chi[4]*chi[8]*chi[9]*chi[12]*chi[13] + chi[4]*chi[8]*chi[9]*chi[12]*chi[15] - chi[6]*chi[8]*chi[9]*chi[12]*chi[13] - chi[4]*chi[8]*chi[9]*chi[13]*chi[15]))/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));
        ftrif[19][23] = -(chi[11]*(chi[0]*chi[8]*chi[9]*chi[13]*chi[14] - chi[2]*chi[8]*chi[9]*chi[12]*chi[14] - chi[4]*chi[8]*chi[9]*chi[12]*chi[13] + chi[6]*chi[8]*chi[9]*chi[12]*chi[13] + chi[6]*chi[8]*chi[9]*chi[12]*chi[14] - chi[6]*chi[8]*chi[9]*chi[13]*chi[14]))/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));


        ftrif[20][0] = -(chi[8]*chi[9]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[20][1] = (chi[8]*chi[9]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[20][2] = (chi[8]*chi[9]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[20][3] = (chi[8]*chi[9]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[20][4] = -(chi[11]*(chi[0]*chi[8]*chi[9]*chi[12]*chi[13]*chi[14] + chi[0]*chi[8]*chi[9]*chi[12]*chi[13]*chi[15] + chi[0]*chi[8]*chi[9]*chi[12]*chi[14]*chi[15] - chi[2]*chi[8]*chi[9]*chi[12]*chi[14]*chi[15] - chi[4]*chi[8]*chi[9]*chi[12]*chi[13]*chi[15] - chi[6]*chi[8]*chi[9]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[20][5] = -(chi[11]*chi[13]*(chi[2]*chi[8]*chi[9]*chi[12]*chi[14] + chi[0]*chi[8]*chi[9]*chi[14]*chi[15] + chi[2]*chi[8]*chi[9]*chi[12]*chi[15] - chi[2]*chi[8]*chi[9]*chi[14]*chi[15] - chi[4]*chi[8]*chi[9]*chi[12]*chi[15] - chi[6]*chi[8]*chi[9]*chi[12]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[20][6] = -(chi[11]*(chi[0]*chi[8]*chi[9]*chi[13]*chi[14]*chi[15] - chi[2]*chi[8]*chi[9]*chi[12]*chi[14]*chi[15] + chi[4]*chi[8]*chi[9]*chi[12]*chi[13]*chi[14] + chi[4]*chi[8]*chi[9]*chi[12]*chi[14]*chi[15] - chi[6]*chi[8]*chi[9]*chi[12]*chi[13]*chi[14] - chi[4]*chi[8]*chi[9]*chi[13]*chi[14]*chi[15]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[20][7] = -(chi[11]*(chi[0]*chi[8]*chi[9]*chi[13]*chi[14]*chi[15] - chi[2]*chi[8]*chi[9]*chi[12]*chi[14]*chi[15] - chi[4]*chi[8]*chi[9]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[12]*chi[14]*chi[15] - chi[6]*chi[8]*chi[9]*chi[13]*chi[14]*chi[15]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[20][8] = (2*chi[11]*(chi[2]*chi[8]*chi[9]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[9]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[9]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[20][9] = -(2*chi[11]*(chi[2]*chi[8]*chi[9]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[9]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[9]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[20][10] = -(2*chi[11]*(chi[2]*chi[8]*chi[9]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[9]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[9]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[20][11] = -(2*chi[11]*(chi[2]*chi[8]*chi[9]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[9]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[9]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[20][12] = (2*chi[1]*chi[8]*chi[9]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[20][13] = -(2*chi[3]*chi[8]*chi[9]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[20][14] = (2*(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]))/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[20][15] = -(2*chi[7]*chi[8]*chi[9]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[20][16] = -(2*chi[11]*(chi[2]*chi[8]*chi[9]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[9]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[9]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[20][17] = (chi[8]*chi[9]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[20][18] = -(2*chi[3]*chi[8]*chi[11])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[20][19] = (2*(chi[1]*chi[9]*chi[11] + chi[3]*chi[8]*chi[11] + chi[7]*chi[8]*chi[9]))/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[20][20] = -(2*chi[7]*chi[8]*chi[9])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[20][21] = -(chi[11]*(chi[2]*chi[8]*chi[9]*chi[12]*chi[14] + chi[0]*chi[8]*chi[9]*chi[14]*chi[15] + chi[2]*chi[8]*chi[9]*chi[12]*chi[15] - chi[2]*chi[8]*chi[9]*chi[14]*chi[15] - chi[4]*chi[8]*chi[9]*chi[12]*chi[15] - chi[6]*chi[8]*chi[9]*chi[12]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[20][22] = -(chi[11]*(chi[0]*chi[8]*chi[9]*chi[13]*chi[15] - chi[2]*chi[8]*chi[9]*chi[12]*chi[15] + chi[4]*chi[8]*chi[9]*chi[12]*chi[13] + chi[4]*chi[8]*chi[9]*chi[12]*chi[15] - chi[6]*chi[8]*chi[9]*chi[12]*chi[13] - chi[4]*chi[8]*chi[9]*chi[13]*chi[15]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[20][23] = -(chi[11]*(chi[0]*chi[8]*chi[9]*chi[13]*chi[14] - chi[2]*chi[8]*chi[9]*chi[12]*chi[14] - chi[4]*chi[8]*chi[9]*chi[12]*chi[13] + chi[6]*chi[8]*chi[9]*chi[12]*chi[13] + chi[6]*chi[8]*chi[9]*chi[12]*chi[14] - chi[6]*chi[8]*chi[9]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);


        ftrif[21][0] = 0;
        ftrif[21][1] = 0;
        ftrif[21][2] = 0;
        ftrif[21][3] = 0;
        ftrif[21][4] = -(chi[12]*chi[13]*chi[14])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[21][5] = -(chi[12]*chi[13]*chi[14])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[21][6] = -(chi[12]*chi[13]*chi[14])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[21][7] = -(chi[12]*chi[13]*chi[14])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[21][8] = -(2*chi[12]*chi[13]*chi[14]*theta)/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[21][9] = (2*chi[12]*chi[13]*chi[14]*theta)/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[21][10] = (2*chi[12]*chi[13]*chi[14]*theta)/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[21][11] = (2*chi[12]*chi[13]*chi[14]*theta)/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[21][12] = 0;
        ftrif[21][13] = 0;
        ftrif[21][14] = 0;
        ftrif[21][15] = 0;
        ftrif[21][16] = (2*chi[12]*chi[13]*chi[14]*theta)/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[21][17] = 0;
        ftrif[21][18] = 0;
        ftrif[21][19] = 0;
        ftrif[21][20] = 0;
        ftrif[21][21] = -(chi[12]*chi[14])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[21][22] = -(chi[12]*chi[13])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);
        ftrif[21][23] = (chi[12]*chi[13] + chi[12]*chi[14] - chi[13]*chi[14])/(chi[12]*chi[13]*chi[14] + chi[12]*chi[13]*chi[15] + chi[12]*chi[14]*chi[15] - chi[13]*chi[14]*chi[15]);


        ftrif[22][0] = -(chi[8]*chi[9]*chi[10])/(2*(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]));
        ftrif[22][1] = (chi[8]*chi[9]*chi[10])/(2*(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]));
        ftrif[22][2] = (chi[8]*chi[9]*chi[10])/(2*(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]));
        ftrif[22][3] = (chi[8]*chi[9]*chi[10])/(2*(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]));
        ftrif[22][4] = -(chi[0]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] + chi[0]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[0]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[2]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[4]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[6]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14])/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));
        ftrif[22][5] = -(chi[10]*chi[13]*(chi[2]*chi[8]*chi[9]*chi[12]*chi[14] + chi[0]*chi[8]*chi[9]*chi[14]*chi[15] + chi[2]*chi[8]*chi[9]*chi[12]*chi[15] - chi[2]*chi[8]*chi[9]*chi[14]*chi[15] - chi[4]*chi[8]*chi[9]*chi[12]*chi[15] - chi[6]*chi[8]*chi[9]*chi[12]*chi[14]))/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));
        ftrif[22][6] = -(chi[0]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] - chi[2]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] + chi[4]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] + chi[4]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[6]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[4]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15])/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));
        ftrif[22][7] = -(chi[0]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] - chi[2]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[4]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[6]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15])/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));
        ftrif[22][8] = (chi[2]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[22][9] = -(chi[2]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[22][10] = -(chi[2]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[22][11] = -(chi[2]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[22][12] = (chi[1]*chi[8]*chi[9]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[22][13] = -(chi[3]*chi[8]*chi[9]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[22][14] = -(chi[5]*chi[8]*chi[9]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[22][15] = -(chi[7]*chi[8]*chi[9]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[22][16] = -(chi[2]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[22][17] = (chi[8]*chi[9]*chi[10])/(2*(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]));
        ftrif[22][18] = -(chi[3]*chi[8]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[22][19] = -(chi[5]*chi[8]*chi[9])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[22][20] = (chi[1]*chi[9]*chi[10] + chi[3]*chi[8]*chi[10] + chi[5]*chi[8]*chi[9])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[22][21] = -(chi[10]*(chi[2]*chi[8]*chi[9]*chi[12]*chi[14] + chi[0]*chi[8]*chi[9]*chi[14]*chi[15] + chi[2]*chi[8]*chi[9]*chi[12]*chi[15] - chi[2]*chi[8]*chi[9]*chi[14]*chi[15] - chi[4]*chi[8]*chi[9]*chi[12]*chi[15] - chi[6]*chi[8]*chi[9]*chi[12]*chi[14]))/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));
        ftrif[22][22] = -(chi[0]*chi[8]*chi[9]*chi[10]*chi[13]*chi[15] - chi[2]*chi[8]*chi[9]*chi[10]*chi[12]*chi[15] + chi[4]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13] + chi[4]*chi[8]*chi[9]*chi[10]*chi[12]*chi[15] - chi[6]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13] - chi[4]*chi[8]*chi[9]*chi[10]*chi[13]*chi[15])/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));
        ftrif[22][23] = -(chi[0]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14] - chi[2]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14] - chi[4]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13] + chi[6]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13] + chi[6]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14] - chi[6]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14])/(2*(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta));


        ftrif[23][0] = -(chi[8]*chi[9]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[23][1] = (chi[8]*chi[9]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[23][2] = (chi[8]*chi[9]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[23][3] = (chi[8]*chi[9]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[23][4] = -(chi[0]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] + chi[0]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[0]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[2]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[4]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[6]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[23][5] = -(chi[10]*chi[13]*(chi[2]*chi[8]*chi[9]*chi[12]*chi[14] + chi[0]*chi[8]*chi[9]*chi[14]*chi[15] + chi[2]*chi[8]*chi[9]*chi[12]*chi[15] - chi[2]*chi[8]*chi[9]*chi[14]*chi[15] - chi[4]*chi[8]*chi[9]*chi[12]*chi[15] - chi[6]*chi[8]*chi[9]*chi[12]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[23][6] = -(chi[0]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] - chi[2]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] + chi[4]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] + chi[4]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[6]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[4]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[23][7] = -(chi[0]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] - chi[2]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[4]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[6]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[23][8] = (2*(chi[2]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[23][9] = -(2*(chi[2]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[23][10] = -(2*(chi[2]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[23][11] = -(2*(chi[2]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[23][12] = (2*chi[1]*chi[8]*chi[9]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[23][13] = -(2*chi[3]*chi[8]*chi[9]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[23][14] = -(2*chi[5]*chi[8]*chi[9]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[23][15] = (2*(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11]))/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[23][16] = -(2*(chi[2]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[0]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15] + chi[4]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] + chi[6]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14] + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14] - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14] - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15] + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15] - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15] + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15] - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]);
        ftrif[23][17] = (chi[8]*chi[9]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[23][18] = -(2*chi[3]*chi[8]*chi[10])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[23][19] = -(2*chi[5]*chi[8]*chi[9])/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[23][20] = (2*(chi[1]*chi[9]*chi[10] + chi[3]*chi[8]*chi[10] + chi[5]*chi[8]*chi[9]))/(chi[1]*chi[9]*chi[10]*chi[11] + chi[3]*chi[8]*chi[10]*chi[11] + chi[5]*chi[8]*chi[9]*chi[11] + chi[7]*chi[8]*chi[9]*chi[10]);
        ftrif[23][21] = -(chi[10]*(chi[2]*chi[8]*chi[9]*chi[12]*chi[14] + chi[0]*chi[8]*chi[9]*chi[14]*chi[15] + chi[2]*chi[8]*chi[9]*chi[12]*chi[15] - chi[2]*chi[8]*chi[9]*chi[14]*chi[15] - chi[4]*chi[8]*chi[9]*chi[12]*chi[15] - chi[6]*chi[8]*chi[9]*chi[12]*chi[14]))/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[23][22] = -(chi[0]*chi[8]*chi[9]*chi[10]*chi[13]*chi[15] - chi[2]*chi[8]*chi[9]*chi[10]*chi[12]*chi[15] + chi[4]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13] + chi[4]*chi[8]*chi[9]*chi[10]*chi[12]*chi[15] - chi[6]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13] - chi[4]*chi[8]*chi[9]*chi[10]*chi[13]*chi[15])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        ftrif[23][23] = -(chi[0]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14] - chi[2]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14] - chi[4]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13] + chi[6]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13] + chi[6]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14] - chi[6]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14])/(chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[14]*theta + chi[1]*chi[9]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[14]*theta - chi[1]*chi[9]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[3]*chi[8]*chi[10]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[13]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[14]*theta - chi[3]*chi[8]*chi[10]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[5]*chi[8]*chi[9]*chi[11]*chi[12]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[13]*chi[15]*theta - chi[5]*chi[8]*chi[9]*chi[11]*chi[13]*chi[14]*chi[15]*theta + chi[7]*chi[8]*chi[9]*chi[10]*chi[12]*chi[14]*chi[15]*theta - chi[7]*chi[8]*chi[9]*chi[10]*chi[13]*chi[14]*chi[15]*theta);
        // Check whether solution is close enough. If not run the loop again.
        int ch = zeroNoDcmp(xb, 24, tol, tol, fvec, ftrif); //
//        int ch = zero(xb, 24, tol, tol, fvec, ftrif);
        if (ch == 1) ok = true;

        j = j + 1;
    }

    // Solutions is applied, and right boundary is updated.
    PV -> Anew[N] = xb[12];
    PV -> Qnew[N] = xb[0];
    LD -> Anew[0] = xb[15];
    LD -> Qnew[0] = xb[3];
    MD -> Anew[0] = xb[18];
    MD -> Qnew[0] = xb[6];
    RD -> Anew[0] = xb[21];
    RD -> Qnew[0] = xb[9];
    //
    if (j >= ntrial) {printf("Error (junction.C)\n Root not found in the trifurcation\n Arteries[%i]\nqLnb = %i", parent, qLnb); exit(222);}
}
