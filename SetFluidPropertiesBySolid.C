// Filename: SetFluidPropertiesBySolid.cpp
// Created on 16 May, 2018, by Hao Gao
// based on the file bellow
// Filename: SetFluidProperties.cpp
// Created on Dec 17, 2017 by Nishant Nangia

// APPLICATION INCLUDES
#include "SetFluidPropertiesBySolid.h"

#include <CartesianGridGeometry.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/HierarchyMathOps.h>

// C++ INCLUDES

/////////////////////////////// STATIC ///////////////////////////////////////

// there is no need to change this one, it is the interface called by VCINSHI
// however the interface is different here
void
callSetFluidDensityCallbackFunctionBySolid(int rho_idx,
                                    Pointer<Variable<NDIM> > rho_var,
                                    Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                    const int cycle_num,
                                    const double time,
                                    const double current_time,
                                    const double new_time,
                                    void* ctx)
{
    // Set the density from the level set information
    static SetFluidPropertiesBySolid* ptr_SetFluidProperties = static_cast<SetFluidPropertiesBySolid*>(ctx);
    ptr_SetFluidProperties->setDensityPatchData(rho_idx, rho_var, hier_math_ops, cycle_num, time, current_time, new_time);

    return;

} // callSetFluidDensityCallbackFunction

// their is no need to change for this one, it is the interface called by VCINSHI
void
callSetFluidViscosityCallbackFunctionBySolid(int mu_idx,
                                      Pointer<Variable<NDIM> > mu_var,
                                      Pointer<IBTK::HierarchyMathOps> hier_math_ops,
                                      const int cycle_num,
                                      const double time,
                                      const double current_time,
                                      const double new_time,
                                      void* ctx)
{
    // Set the density from the level set information
    static SetFluidPropertiesBySolid* ptr_SetFluidProperties = static_cast<SetFluidPropertiesBySolid*>(ctx);
    ptr_SetFluidProperties->setViscosityPatchData(mu_idx, mu_var, hier_math_ops, cycle_num, time, current_time, new_time);

    return;

} // callSetFluidViscosityCallbackFunction

/////////////////////////////// PUBLIC //////////////////////////////////////

SetFluidPropertiesBySolid::SetFluidPropertiesBySolid(const std::string& object_name,
                                       const std::string& ls_name,
                                       const unsigned int ls_idx,
                                       const std::string& phiMesh_name,
                                       const unsigned int phiMesh_idx,
                                       const double rho_solid,
                                       const double rho_fluid,
                                       const double mu_solid,
                                       const double mu_fluid)
    : d_object_name(object_name),
      d_ls_name(ls_name),
      d_ls_idx(ls_idx),
      d_phiMesh_name(phiMesh_name),
      d_phiMesh_idx(phiMesh_idx),
      d_rho_solid(rho_solid),
      d_rho_fluid(rho_fluid),
      d_mu_solid(mu_solid),
      d_mu_fluid(mu_fluid)
{
    // intentionally left blank
    return;
} // SetFluidProperties

SetFluidPropertiesBySolid::~SetFluidPropertiesBySolid()
{
    // intentionally left blank
    return;

} //~SetFluidProperties

void
SetFluidPropertiesBySolid::setDensityPatchData(int rho_idx,
                                        Pointer<Variable<NDIM> > rho_var,
                                        SAMRAI::tbox::Pointer<HierarchyMathOps> hier_math_ops,
                                        const int cycle_num,
                                        const double time,
                                        const double /*current_time*/,
                                        const double /*new_time*/)
{
    // Get the current level set information
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    //Pointer<CellVariable<NDIM, double> > ls_var = d_adv_diff_solver->getLevelSetVariable(d_ls_name);
    //const int ls_current_idx = var_db->mapVariableAndContextToIndex(ls_var, d_adv_diff_solver->getCurrentContext());
    const int ls_current_idx = d_ls_idx;

    // Set the density based on the level set
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            //const double* const patch_dx = patch_geom->getDx();
            //double vol_cell = 1.0;
            //for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
            //double alpha = d_num_interface_cells * std::pow(vol_cell, 1.0 / (double)NDIM);

            const Box<NDIM>& patch_box = patch->getBox();
            const Pointer<CellData<NDIM, double> > ls_data = patch->getPatchData(ls_current_idx);
            Pointer<CellData<NDIM, double> > rho_data = patch->getPatchData(rho_idx);

            const Pointer<CellData<NDIM, double> > phiMesh_data = patch->getPatchData(d_phiMesh_idx);

            // Calderer et al, 2014
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());
                const double m_added = (*ls_data)(ci);
                const double phi_mesh_coef = (*phiMesh_data)(ci);

                /* considering spreading a density over the region, then needs to
                * rescale using phi_mesh_coef, which is spreaded from 1 in the solid mesh
                * if phi_Mesh_coef < 0.01, then outside the solid region, that is fluid
                */
                /*double m_added_scaled = m_added;
                if (phi_mesh_coef <= 1.0)
                {
                  m_added_scaled = m_added;
                }
                else
                {
                  m_added_scaled = m_added / phi_mesh_coef;
                }
                */

                (*rho_data)(ci) = d_rho_solid ; // in eulerian form the density is constant due to incompressibility in both the fluid and structure parts 
            }
        }
    }

    return;
} // setDensityPatchData

void
SetFluidPropertiesBySolid::setViscosityPatchData(int mu_idx,
                                          Pointer<Variable<NDIM> > mu_var,
                                          SAMRAI::tbox::Pointer<HierarchyMathOps> hier_math_ops,
                                          const int cycle_num,
                                          const double time,
                                          const double /*current_time*/,
                                          const double /*new_time*/)
{
    // Get the current level set information
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    //Pointer<CellVariable<NDIM, double> > ls_var = d_adv_diff_solver->getLevelSetVariable(d_ls_name);
    //const int ls_current_idx = var_db->mapVariableAndContextToIndex(ls_var, d_adv_diff_solver->getCurrentContext());
    const int ls_current_idx = d_ls_idx;

    // Set the density based on the level set
    Pointer<PatchHierarchy<NDIM> > patch_hierarchy = hier_math_ops->getPatchHierarchy();
    const int coarsest_ln = 0;
    const int finest_ln = patch_hierarchy->getFinestLevelNumber();

    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(ln);
        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            Pointer<CartesianPatchGeometry<NDIM> > patch_geom = patch->getPatchGeometry();
            //const double* const patch_dx = patch_geom->getDx();
            //double vol_cell = 1.0;
            //for (int d = 0; d < NDIM; ++d) vol_cell *= patch_dx[d];
            //double alpha = d_num_interface_cells * std::pow(vol_cell, 1.0 / (double)NDIM);

            const Box<NDIM>& patch_box = patch->getBox();
            //const Pointer<CellData<NDIM, double> > ls_data = patch->getPatchData(ls_current_idx);
            Pointer<CellData<NDIM, double> > mu_data = patch->getPatchData(mu_idx);

            // Calderer et al, 2014
            for (Box<NDIM>::Iterator it(patch_box); it; it++)
            {
                CellIndex<NDIM> ci(it());
                //const double phi = (*ls_data)(ci);
                //double h_phi;
                //if (phi < 0)
                //{
                //    h_phi = 0.0;
                //}
                //else if (phi > 1)
                //{
                //    h_phi = 1.0;
                //}
                //else
                //{
                //    h_phi = phi;
                //}

                (*mu_data)(ci) = d_mu_solid; // mu will be constant
                //(*mu_data)(ci) = d_mu_solid*h_phi + d_mu_fluid*(1.0-h_phi);
            }
        }
    }

    return;
} // setViscosityPatchData

/////////////////////////////// PRIVATE //////////////////////////////////////
