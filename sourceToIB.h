// IBAMR INCLUDES
#include <ibamr/IBFEMethod.h>
#include <ibamr/app_namespaces.h>

#include <string>

// DarcySolver is a static class that provides data and functions
// required to sovle the darcy flow.
class sourceToIB
{
public:





    static std::string SstarSystemName;
    static std::string SstarvariableName;
    static Order Sstarvariable_order;



    static double time;
    static double time_step;


//---------------------------------------------------------

    static void
    setup_Sstar(EquationSystems* equation_systems);

    static void
    initialize_Sstar(EquationSystems* equation_systems);


    static void
    assemble_Sstar(EquationSystems& es,const std::string&  system_name);


    static void
    solve_Sstar(EquationSystems* es,const double,const double);


//---------------------------------------------------------


private:
    sourceToIB();
    sourceToIB(sourceToIB&);
    ~sourceToIB();
    sourceToIB& operator=(sourceToIB&);
};
