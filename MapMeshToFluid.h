// All rights reserved.
//

#ifndef included_MapMeshToFluid
#define included_MapMeshToFluid

// IBAMR INCLUDES
#include <ibamr/IBFEMethod.h>
#include <ibamr/app_namespaces.h>
#include <string>

//SAMARI
#include <CartesianGridGeometry.h>
#include <ibtk/HierarchyMathOps.h>



// This class will map the deformed structure mesh into an Eulerian field
class MapMeshToFluid
{
   public:
     //variable in eulerian space
     static std::string phi_mesh_var_name; 
     static std::string phi_mesh_ctx_name;
     
     //variable in solid
     static std::string phiMeshSystemName;
     static std::string phiMeshVariableName;
     
     static void initialize_equation_systems(
        EquationSystems* equation_systems);
        
     static void initialize_phi(EquationSystems *es);
     
     static void spread_phi_solid_fluid(
        SAMRAI::tbox::Pointer<IBFEMethod> ib_method_ops,
        SAMRAI::tbox::Pointer<PatchHierarchy<NDIM> > patch_hierarchy, 
        const int phi_idx,
        const double data_time);
    

  private:
    MapMeshToFluid();
    MapMeshToFluid(MapMeshToFluid&);
    ~MapMeshToFluid();
    MapMeshToFluid& operator=(MapMeshToFluid&);
};

#endif
