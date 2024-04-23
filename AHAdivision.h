// All rights reserved.
//

#ifndef included_AHAdivision
#define included_AHAdivision

// IBAMR INCLUDES
#include <ibamr/IBFEMethod.h>
#include <ibamr/app_namespaces.h>


#include <string>

extern double   total_aha_volumes[17];
extern int  node_aha_list[16403],element_aha_list[82290];

class AHAdivision
{
public:

    static std::vector<int>  node_aha_list;
    static std::vector<int>  element_aha_list;
    static double total_aha_volumes[17];

    static std::string mesh_division_systemName;
    static std::string mesh_division_variableName;
    static Order mesh_division_variable_order;


             static void
             setup_mesh_division(
                 EquationSystems* equation_systems);

             static void
             initialize_mesh_division(
                     EquationSystems* equation_systems);

            static void
              update_mesh_division(EquationSystems* equation_systems);


              static void
                find_aha_volumes(EquationSystems* equation_systems);

//---------------------------------------------------------


private:
    AHAdivision();
    AHAdivision(AHAdivision&);
    ~AHAdivision();
    AHAdivision& operator=(AHAdivision&);
};


#endif
