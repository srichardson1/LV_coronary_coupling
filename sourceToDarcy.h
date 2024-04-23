// IBAMR INCLUDES
#include <ibamr/IBFEMethod.h>
#include <ibamr/app_namespaces.h>

#include <string>

class sourceToDarcy
{
public:

    static BoundaryInfo* boundary_info;

    //all static varialbes need to be initialized outside the class definition
    static std::string sourceToDarcySystemName;
    static std::string sourceToDarcyvariableName;
    static Order sourceToDarcyvariable_order;

    static double p_a;
    static double p_v;
    static double beta_a;
    static double beta_v;

    static double time_std;
    static double time_step_std;

    static void
    setup_std(EquationSystems* equation_systems);

    static void
    initialize_std(EquationSystems* equation_systems);

    static void
    assemble_std(EquationSystems& es,const std::string&  system_name);

    static void
    solve_std(EquationSystems* es);



  private:
      sourceToDarcy();
      sourceToDarcy(sourceToDarcy&);
      ~sourceToDarcy();
      sourceToDarcy& operator=(sourceToDarcy&);


  };
