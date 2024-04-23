// All rights reserved.
//

#ifndef included_MapmAddedToFluid
#define included_MapmAddedToFluid

// IBAMR INCLUDES
#include <ibamr/IBFEMethod.h>
#include <ibamr/app_namespaces.h>
#include <string>

//SAMARI
#include <CartesianGridGeometry.h>
#include <ibtk/HierarchyMathOps.h>



// This class will map the added m from structural part into an Eulerian field
class MapmAddedToFluid
{
   public:
     //variable in eulerian space
     static std::string m_added_var_name;
     static std::string m_added_ctx_name;

     //variable in solid
     static std::string mAddedSystemName;
     static std::string mAddedVariableName;

     static void initialize_equation_systems(
        EquationSystems* equation_systems);

     static void initialize_m_added(EquationSystems *es);

     static void spread_m_added_solid_fluid(
        SAMRAI::tbox::Pointer<IBFEMethod> ib_method_ops,
        SAMRAI::tbox::Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
        const int idx,
        const double data_time);


  private:
    MapmAddedToFluid();
    MapmAddedToFluid(MapmAddedToFluid&);
    ~MapmAddedToFluid();
    MapmAddedToFluid& operator=(MapmAddedToFluid&);
};

#endif
