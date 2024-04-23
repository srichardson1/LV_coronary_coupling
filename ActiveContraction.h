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

#ifndef included_ActiveContraction
#define included_ActiveContraction

// IBAMR INCLUDES
#include <ibamr/IBFEMethod.h>
#include <ibamr/app_namespaces.h>

// ActiveContraction is a static class that provides data and functions required
// to implement the active contraction model.
class ActiveContraction
{
public:
    // Activation variable indexing.
    static const unsigned int     CA_I_IDX = 0;
    static const unsigned int     CA_B_IDX = 1;
    static const unsigned int       Q1_IDX = 2;
    static const unsigned int       Q2_IDX = 3;
    static const unsigned int       Q3_IDX = 4;
    static const unsigned int        Z_IDX = 5;
    static const unsigned int NUM_ACT_VARS = 6;

    static double T_current_loading;
    static double dCa_i_dt_hold;
    static double Ca_i_wanted;

    static double remaining_steps_act;
    static double remaining_steps_rounded_act;
    static double difference_act;

    static int val_test_act;
    static int time_to_enter;

    // Equation system numbers.
    static int act_system_num, T_system_num, T_system_num_holder;

    static void
    update_active_tension_model_state_variables(
        EquationSystems* equation_systems,
        double time,
        double dt);


        static void
        save_reduction_value(
            EquationSystems* equation_systems,
            double time,
            double dt);



    static void
    reinitialize_variables(EquationSystems* equation_systems);


    static void
    NHS_RK2_step(
        double& Ca_i,
        double& Ca_b,
        double& Q1,
        double& Q2,
        double& Q3,
        double& z,
        double lambda,
        double dlambda_dt,
        double time,
        double dt);

    static void
    NHS_euler_step(
        double& Ca_i,
        double& Ca_b,
        double& Q1,
        double& Q2,
        double& Q3,
        double& z,
        double lambda,
        double dlambda_dt,
        double time,
        double dt);

     static void
       runOneProcessor(EquationSystems* equation_systems);

private:
    ActiveContraction();
    ActiveContraction(ActiveContraction&);
    ~ActiveContraction();
    ActiveContraction& operator=(ActiveContraction&);
};

#endif //#ifndef included_ActiveContraction
