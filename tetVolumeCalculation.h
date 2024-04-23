#ifndef TETVOLUMECALCULATION_H
#define TETVOLUMECALCULATION_H

#include <libmesh/mesh_tetgen_wrapper.h>

#include "tetgen.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <math.h>



double tetVolumeCalculation(libMesh::TetGenWrapper& tetGenData);
double tetVolumeCalculationByPoints(std::vector< std::vector<double> >all_points);
double tetVolumeTest(); //a test function
#endif
