// Filename: SetFluidProeprtiesBySolid.h
// Created on 16 May, 2018 by Hao Gao
// based on the file bellow
// Filename: SetFluidProperties.h
// Created on Dec 18, 2017 by Nishant Nangia

/////////////////////// INCLUDE GUARD ////////////////////////////////////

#ifndef included_SetFluidPropertiesBySolid
#define included_SetFluidPropertiesBySolid

///////////////////////////// INCLUDES ///////////////////////////////////

#include <ibamr/AdvDiffHierarchyIntegrator.h>
#include <ibtk/ibtk_utilities.h>
#include <tbox/Pointer.h>
#include <Variable.h>

namespace IBTK
{
class HierarchyMathOps;
}

/*!
 * Pre processing call back function to be hooked into IBAMR::AdvDiffHierarchyIntegrator class.
 *
 * \param rho_idx a patch data index for the current density variable maintained by the integrator.
 * \param ctx is the pointer to SetFluidProperties class object.
 */

void callSetFluidDensityCallbackFunctionBySolid(int rho_idx,
                                         SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > rho_var,
                                         SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                         const int cycle_num,
                                         const double time,
                                         const double current_time,
                                         const double new_time,
                                         void* ctx);

/*!
 * Pre processing call back function to be hooked into IBAMR::AdvDiffHierarchyIntegrator class.
 *
 * \param rho_idx a patch data index for the current density variable maintained by the integrator.
 * \param ctx is the pointer to SetFluidProperties class object.
 */

void callSetFluidViscosityCallbackFunctionBySolid(int mu_idx,
                                           SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > mu_var,
                                           SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                           const int cycle_num,
                                           const double time,
                                           const double current_time,
                                           const double new_time,
                                           void* ctx);

class SetFluidPropertiesBySolid
{
    /*!
     * \brief Class SetFluidProperties is a utility class which sets the fluid and
     * solid Eulerian density based on the current solid structure position
     */
public:
    /*!
     * The only constructor of this class.
     */
    SetFluidPropertiesBySolid(const std::string& object_name,
                       const std::string& ls_name,
                       const unsigned int ls_idx,
                       const std::string& phiMesh_name,
                       const unsigned int phiMesh_idx,
                       const double rho_solid,
                       const double rho_fluid,
                       const double mu_solid,
                       const double mu_fluid);

    /*!
     * Destructor for this class.
     */
    ~SetFluidPropertiesBySolid();

    /*!
     * Set the density based on the current level set information
     */
    void setDensityPatchData(int rho_idx,
                             SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > rho_var,
                             SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                             const int cycle_num,
                             const double time,
                             const double current_time,
                             const double new_time);

    /*!
     * Set the viscosity based on the current level set information
     */
    void setViscosityPatchData(int mu_idx,
                               SAMRAI::tbox::Pointer<SAMRAI::hier::Variable<NDIM> > mu_var,
                               SAMRAI::tbox::Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                               const int cycle_num,
                               const double time,
                               const double current_time,
                               const double new_time);

    //////////////// PRIVATE /////////////////////////////

private:
    /*!
     * Default constructor is not implemented and should not be used.
     */
    SetFluidPropertiesBySolid();

    /*!
     * Default assignment operator is not implemented and should not be used.
     */
    SetFluidPropertiesBySolid& operator=(const SetFluidPropertiesBySolid& that);

    /*!
     * Default copy constructor is not implemented and should not be used.
     */
    SetFluidPropertiesBySolid(const SetFluidPropertiesBySolid& from);

    /*!
     * Name of this object.
     */
    std::string d_object_name;

    /*!
     * Pointer to advection-diffusion solver. no need anymore, the interface is tracked by the solid
     */
    // SAMRAI::tbox::Pointer<IBAMR::AdvDiffHierarchyIntegrator> d_adv_diff_solver;

    /*!
     * Name of the level set.
     */
    std::string d_ls_name;
    /*!
     * idx for added_mass
     */
     unsigned int d_ls_idx;
   /*
   * name and idx for mesh labelled region
   */
   std::string d_phiMesh_name;
   unsigned int d_phiMesh_idx;

    /*!
     * Density of the fluid.
     */
    double d_rho_solid, d_rho_fluid;

    /*!
     * Viscosity of the fluid.
     */
    double d_mu_solid, d_mu_fluid;

    /*!
     * Level set reinitialization interval no need anymore
     */
    //int d_ls_reinit_interval;

    /*!
     * Number of interface cells over which to smooth the material properties
     */
    // double d_num_interface_cells;

}; // SetFluidPropertiesBySolid

#endif // #ifndef included_SetFluidPropertiesBySolid
