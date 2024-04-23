// IBAMR INCLUDES
#include <ibamr/IBFEMethod.h>
#include <ibamr/app_namespaces.h>

#include <string>



class DarcyMass
{
public:

  static int m_system_num;
  static int mold_system_num;

  static double phi0;
  static double biot_modulus;
  static double kappa0;

  static std::string moldSystemName;
  static std::string moldvariableName;
  static Order moldvariable_order;

  static std::string mSystemName;
  static std::string mvariableName;
  static Order mvariable_order;

  static std::string maverageSystemName;
  static std::string maveragevariableName;
  static Order maveragevariable_order;

  static std::string mstarSystemName;
  static std::string mstarvariableName;
  static Order mstarvariable_order;

  static std::string mextraSystemName;
  static std::string mextravariableName;
  static Order mextravariable_order;

  static double time;
  static double time_step;

  static void
  setup_m(EquationSystems* equation_systems);

  static void
  setup_maverage(EquationSystems* equation_systems);

  static void
  setup_mstar(EquationSystems* equation_systems);

  static void
  setup_mold(EquationSystems* equation_systems);

  static void
  setup_mextra(EquationSystems* equation_systems);

  static void
  initialize_m(EquationSystems* equation_systems);

  static void
  initialize_mstar(EquationSystems* equation_systems);

  static void
  initialize_mold(EquationSystems* equation_systems);

  static void
  initialize_maverage(EquationSystems* equation_systems);

  static void
  initialize_mextra(EquationSystems* equation_systems);


  static void
  assemble_mstar(EquationSystems& es,const std::string&  system_name);

  static void
  assemble_mextra(EquationSystems& es,const std::string&  system_name);

  static void
  assemble_m(EquationSystems& es,const std::string&  system_name);

  static void
  assemble_mold(EquationSystems& es,const std::string&  system_name);

  static void
  assemble_maverage(EquationSystems& es,const std::string&  system_name);

  static void
  solve_mstar(EquationSystems* es,const double,const double);

  static void
  solve_mextra(EquationSystems* es,const double,const double);

  static void
  solve_m(EquationSystems* es,const double,const double);

  static void
  solve_maverage(EquationSystems* es,const double,const double);

  static void
  solve_mold(EquationSystems* es);



  private:
      DarcyMass();
      DarcyMass(DarcyMass&);
      ~DarcyMass();
      DarcyMass& operator=(DarcyMass&);


  };
