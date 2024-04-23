#ifndef included_AdvanceCoronaryCode
#define included_AdvanceCoronaryCode

// IBAMR INCLUDES
#include <ibamr/IBFEMethod.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/libmesh_utilities.h>

// LIBMESH INCLUDES
#include <libmesh/dof_map.h>
#include <libmesh/equation_systems.h>
#include <libmesh/quadrature.h>
#include <libmesh/point.h>

#include <string>


class AdvanceCoronaryCode
{
public:


static void
advance_variables(EquationSystems* equation_systems);





private:
    AdvanceCoronaryCode();
    AdvanceCoronaryCode(AdvanceCoronaryCode&);
    ~AdvanceCoronaryCode();
    AdvanceCoronaryCode& operator=(AdvanceCoronaryCode&);



};
#endif
