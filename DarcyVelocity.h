// IBAMR INCLUDES
#include <ibamr/IBFEMethod.h>
#include <ibamr/app_namespaces.h>

#include <string>


class DarcyVelocity
{
public:


  static std::string wxSystemName;
  static std::string wxvariableName;
  static Order wxvariable_order;

  static std::string wySystemName;
  static std::string wyvariableName;
  static Order wyvariable_order;

  static std::string wzSystemName;
  static std::string wzvariableName;
  static Order wzvariable_order;


  static std::string wxstarSystemName;
  static std::string wxstarvariableName;
  static Order wxstarvariable_order;

  static std::string wystarSystemName;
  static std::string wystarvariableName;
  static Order wystarvariable_order;

  static std::string wzstarSystemName;
  static std::string wzstarvariableName;
  static Order wzstarvariable_order;

  static double time;
  static double time_step;



  static void
  setup_wx(
          EquationSystems* equation_systems);

  static void
  setup_wy(
          EquationSystems* equation_systems);

  static void
  setup_wz(
          EquationSystems* equation_systems);



  static void
  setup_wxstar(EquationSystems* equation_systems);

  static void
  setup_wystar(EquationSystems* equation_systems);

  static void
  setup_wzstar(EquationSystems* equation_systems);

  static void
  initialize_wx(EquationSystems* equation_systems);

  static void
  initialize_wy(EquationSystems* equation_systems);

  static void
  initialize_wz(EquationSystems* equation_systems);


  static void
  initialize_wxstar(EquationSystems* equation_systems);

  static void
  initialize_wystar(EquationSystems* equation_systems);

  static void
  initialize_wzstar(EquationSystems* equation_systems);



  static void
  assemble_wx(EquationSystems& es,const std::string&  system_name);

  static void
  assemble_wy(EquationSystems& es,const std::string&  system_name);

  static void
  assemble_wz(EquationSystems& es,const std::string&  system_name);


  static void
  assemble_wxstar(EquationSystems& es,const std::string&  system_name);

  static void
  assemble_wystar(EquationSystems& es,const std::string&  system_name);

  static void
  assemble_wzstar(EquationSystems& es,const std::string&  system_name);


  static void
  solve_wx(EquationSystems* es);

  static void
  solve_wy(EquationSystems* es);

  static void
  solve_wz(EquationSystems* es);

  static void
  solve_wxstar(EquationSystems* es);

  static void
  solve_wystar(EquationSystems* es);

  static void
  solve_wzstar(EquationSystems* es);



  private:
      DarcyVelocity();
      DarcyVelocity(DarcyVelocity&);
      ~DarcyVelocity();
      DarcyVelocity& operator=(DarcyVelocity&);


  };
