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

// IBAMR INCLUDES
#include <ibamr/IBFEMethod.h>
#include <ibamr/app_namespaces.h>

// MechanicsModel is a static class that provides data and functions required to
// implement the LV mechanics model.
class MechanicsModel
{
public:
    static int f0_system_num, s0_system_num;
    static bool enable_active_tension, normalize_stress;
    static double T_scale, beta_s;
    static double I1_dev_max, I1_dev_min, I1_dil_max, I1_dil_min;
    static double J_dev_max, J_dev_min, J_dil_max, J_dil_min;

    static int J_system_num, FF_system_num;
    static string J_system_name, FF_system_name;
    static string dJdt_system_name;
    static int dJdt_system_num;
    static string Jold_system_name;
    static int Jold_system_num;


    //static void
    //get_PK1_dev_stress_function_systems(
    //    vector<unsigned int>& systems);

    //static void
    //get_PK1_dev_stress_function_systems(
    //     std::vector<IBTK::SystemData> system_data);

    static void
    PK1_dev_stress_function(
        TensorValue<double>& PP,
        const TensorValue<double>& FF,
        const libMesh::Point& x,
        const libMesh::Point& X,
        Elem* const elem,
        const vector<const vector<double>*>& system_data,
        const vector<const vector<libMesh::VectorValue<double> >*>& system_grad_var_data,
        double data_time,
        void* ctx);

    static void
    PK1_dil_stress_function(
        TensorValue<double>& PP,
        const TensorValue<double>& FF,
        const libMesh::Point& x,
        const libMesh::Point& X,
        Elem* const elem,
        const vector<const vector<double>*>& system_data,
        const vector<const vector<libMesh::VectorValue<double> >*>& system_grad_var_data,
        double data_time,
        void* ctx);


    static void initialize_J_system(EquationSystems* es);

    static void assemble_J_system( EquationSystems& es,
                                         const std::string& system_name);
    static void
           update_J_system(EquationSystems* equation_systems, double dt);

    static void initialize_FF_system(EquationSystems* es);
    static void assemble_FF_system(EquationSystems& es, const std::string& system_name);
    static void update_FF_system(EquationSystems* equations_systems);



private:
    MechanicsModel();
    MechanicsModel(MechanicsModel&);
    ~MechanicsModel();
    MechanicsModel& operator=(MechanicsModel&);
};
