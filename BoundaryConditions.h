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

#ifndef included_BoundaryConditions
#define included_BoundaryConditions


// IBAMR INCLUDES
#include <ibamr/IBFEMethod.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/libmesh_utilities.h>

// LIBMESH INCLUDES
#include <libmesh/dof_map.h>
#include <libmesh/equation_systems.h>
#include <libmesh/quadrature.h>
#include <libmesh/point.h>

#include <string>

// BoundaryConditions is a static class that provides data and functions
// required to implement loading conditions and boundary constraints for the LV
// model.
class BoundaryConditions
{
public:
    static double P_load, t_load, kappa;
    static BoundaryInfo* boundary_info;
    static double P_current_loading;
    static double P_current_loading_RV;
    static double t_end_diastole;
    static double t_end_systole;
    static double t_isovolumetric_relaxation_end;
    static double t_period;
    static double cycle_no;
    static int current_cycle;
    static double cycle_time_saved;
    static double single_cycle_length;
    static double diastole_systole_difference;
    static int total_node_num;
    static int total_elem_num;
    static int total_endo_nodes;

    static double p_loading_bcast_from_root;
    static double LV_volume;

    static double whole_volume;
    static double wall_volume;

    static std::vector<int>  LV_endo_points_list;
    static std::vector< std::vector<double>  > LV_endo_points;
    static std::vector< std::vector<double>  > LV_all_points;
    static int LV_NoOfEndoNode;

    static double P_old_loading;
    static double dP;
    static double LV_volume_old;
    static double LV_CIsoContract;
    static double LV_CIsoRelax;
    static double dt;

    static int LV_states;
    static double P_AV_open;
    static double P_AV_close;
    static double LV_volume_ED;
    static double LV_volume_ES;
    static double Pwk_AV;
    static double P_Ao;
    static double Pwk_AV_debug;
    static double Pwk_AV_debug_old;
    static double Pwk_LV_debug;
    static double QLV_debug;
    static double pao_initial;
    static double pwk_initial;

    static bool diastolic_filling;
    static bool isovolumetric_contraction;
    static bool systolic_ejection;
    static bool isovolumetric_relaxation;
    static bool early_diastolic_relaxation;

    static double Ta_relax_factor;

    static double p_aortic_hold;

    static double remaining_steps;
    static double p_ao_step_difference;
    static double p_wk_step_difference;
    static double loading_p_step_difference;

    static int remaining_steps_rounded;
    static int val_test;

    static double
    loading_pressure(
        double time);

    static void
    loading_force_function(
        double& P,
        const libMesh::VectorValue<double>& n,
        const libMesh::VectorValue<double>& N,
        const TensorValue<double>& FF,
        const libMesh::Point& x,
        const libMesh::Point& X,
        Elem* const elem,
        unsigned short int side,
        const vector<const vector<double>*>& system_data,
        const vector<const vector<libMesh::VectorValue<double> >*>& system_grad_var_data,
        double time,
        void* ctx);

    static void
    tether_force_function(
        VectorValue<double>& F,
        const libMesh::VectorValue<double>& n,
        const libMesh::VectorValue<double>& N,
        const TensorValue<double>& FF,
        const libMesh::Point& x,
        const libMesh::Point& X,
        Elem* const elem,
        const unsigned short int side,
        const vector<const vector<double>*>& system_data,
        const vector<const vector<libMesh::VectorValue<double> >*>& system_grad_var_data,
        double time,
        void* ctx);


    static void
    readingPoints(MeshBase& mesh);

    static double lump_para_pressure(EquationSystems* equation_systems,
      double time);

    static void phase_determination(double time);

    static void
    readingPointsGeneral(MeshBase& mesh,
                                    std::vector<int>& points_list,
                                    std::vector< std::vector<double> >& points,
                                    std::vector< std::vector<double> >& points_sec,
                                    int& NoOfPoints,
                                    std::string file_name);

    static void
      updatePointsPosition(EquationSystems* equation_systems);

    static void
      updatePointsPositionGeneral(EquationSystems* equation_systems,
                                                    std::vector<int> & points_list,
                                     std::vector< std::vector<double> >& points,
                                     std::vector< std::vector<double> >& points_sec,
                                                    int NoOfPoints);

    static void
    solve_wall_volume(EquationSystems* equation_systems);

    static void P_Ao_compute(double pwk_prev,double q_prev);

    static void print_info(double time);

private:
    BoundaryConditions();
    BoundaryConditions(BoundaryConditions&);
    ~BoundaryConditions();
    BoundaryConditions& operator=(BoundaryConditions&);
};

#endif
