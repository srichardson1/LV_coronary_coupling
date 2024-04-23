// IBAMR INCLUDES
#include <ibamr/IBFEMethod.h>
#include <ibamr/app_namespaces.h>

#include <string>



class DarcyPressure
{
public:

  static double phi0;
  static double biot_modulus;
  static double kappa0;
  static double cap_K;

  static std::string pSystemName;
  static std::string pvariableName;
  static Order pvariable_order;


  static std::string pstarSystemName;
  static std::string pstarvariableName;
  static Order pstarvariable_order;


  static double time;
  static double time_step;


  static void
  setup_p(EquationSystems* equation_systems);


  static void
  setup_pstar(EquationSystems* equation_systems);


  static void
  initialize_p(EquationSystems* equation_systems);


  static void
  initialize_pstar(EquationSystems* equation_systems);


  static void
  assemble_p(EquationSystems& es,const std::string&  system_name);

  static void
  assemble_pstar(EquationSystems& es,const std::string&  system_name);

  static void
  solve_p(EquationSystems* es,const double,const double);


  static void
  solve_pstar(EquationSystems* es,const double,const double);



  private:
      DarcyPressure();
      DarcyPressure(DarcyPressure&);
      ~DarcyPressure();
      DarcyPressure& operator=(DarcyPressure&);


  };
