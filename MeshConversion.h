// IBAMR INCLUDES
#include <ibamr/IBFEMethod.h>
#include <ibamr/app_namespaces.h>

#include <string>



class MeshConversion
{
public:

  static std::string meshlongSystemName;
  static std::string meshlongvariableName;
  static Order meshlongvariable_order;

  static std::string meshradialSystemName;
  static std::string meshradialvariableName;
  static Order meshradialvariable_order;

  static std::string meshcircumSystemName;
  static std::string meshcircumvariableName;
  static Order meshcircumvariable_order;


  static void
  setup_mesh_long(EquationSystems* equation_systems);

  static void
  initialize_mesh_long(EquationSystems* equation_systems);

  static void
  assemble_mesh_long(EquationSystems& es,const std::string&  system_name);

  static void
  solve_mesh_long(EquationSystems* es);



  static void
  setup_mesh_radial(EquationSystems* equation_systems);

  static void
  initialize_mesh_radial(EquationSystems* equation_systems);

  static void
  assemble_mesh_radial(EquationSystems& es,const std::string&  system_name);

  static void
  solve_mesh_radial(EquationSystems* es);




  static void
  setup_mesh_circum(EquationSystems* equation_systems);

  static void
  initialize_mesh_circum(EquationSystems* equation_systems);

  static void
  assemble_mesh_circum(EquationSystems& es,const std::string&  system_name);

  static void
  solve_mesh_circum(EquationSystems* es);




  private:
      MeshConversion();
      MeshConversion(MeshConversion&);
      ~MeshConversion();
      MeshConversion& operator=(MeshConversion&);


  };
