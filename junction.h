// header for junction.c
// J. A. Mackenzie.
// University of Glasgow

#ifndef _JUNCTION_H
#define _JUNCTION_H

#include <cmath>
#include "tools.h"

void junction (int qLnb, double theta, double gamma, Tube *Arteries[], int parent);
void monofurcation (int qLnb, double theta, double gamma, Tube *Arteries[], int parent, int ntrial, double tol, int j, int ok);
void bifurcation (int qLnb, double theta, double gamma, Tube *Arteries[], int parent, int ntrial, double tol, int j, int ok);
void trifurcation (int qLnb, double theta, double gamma, Tube *Arteries[], int parent, int ntrial, double tol, int j, int ok);

#endif
