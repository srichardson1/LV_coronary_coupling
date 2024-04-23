// Copyright (c) 2011-2013, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// APPLICATION INCLUDES
#include <ModelInitialization.h>
#include <MechanicsModel.h>
#include <ActiveContraction.h>

// LIBMESH INCLUDES
#include <libmesh/centroid_partitioner.h>
#include <libmesh/equation_systems.h>
#include <libmesh/explicit_system.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_tools.h>

// C++ STDLIB INCLUDES
#include <map>

// CLASS IMPLEMENTATION

void
ModelInitialization::initialize_mesh(
                                     Mesh& mesh,
                                     Pointer<Database> input_db)
{
    pout << "Loading the mesh...\n";
    mesh.allow_renumbering(false);
    mesh.partitioner() = UniquePtr<Partitioner>(new CentroidPartitioner());  // BEG: Needed to run on cardiac (not sure why)
    mesh.read(input_db->getString("MESH_FILENAME"), NULL);
    const string elem_order = input_db->getString("ELEM_ORDER");
    if (elem_order != "FIRST" && elem_order != "SECOND")
    {
        TBOX_ERROR("invalid element order: " << elem_order << "\n");
    }
    if (elem_order == "SECOND") mesh.all_second_order(true);
    pout << "number of elems = " << mesh.n_elem()  << "\n"
    << "number of nodes = " << mesh.n_nodes() << "\n"
    << "max elem ID = " << mesh.max_elem_id() << "\n"
    << "max node ID = " << mesh.max_node_id() << "\n";
    plog << "processor " << SAMRAI_MPI::getRank() << ": number of        local elems = " << distance(mesh.       local_elements_begin(),mesh.       local_elements_end()) << "\n";
    plog << "processor " << SAMRAI_MPI::getRank() << ": number of active local elems = " << distance(mesh.active_local_elements_begin(),mesh.active_local_elements_end()) << "\n";
    plog << "processor " << SAMRAI_MPI::getRank() << ": number of        local nodes = " << distance(mesh.       local_nodes_begin()   ,mesh.       local_nodes_end()   ) << "\n";

    // Rescale the mesh so that length is in cm.
    const string mesh_length_unit = input_db->getString("MESH_LENGTH_UNIT");
    pout << "input mesh length unit is: " << mesh_length_unit << "\n";
    if (strcasecmp(mesh_length_unit.c_str(),"m") == 0)
    {
        MeshTools::Modification::scale(mesh, 100.0);
    }
    else if (strcasecmp(mesh_length_unit.c_str(),"cm") == 0)
    {
        // do nothing
    }
    else if (strcasecmp(mesh_length_unit.c_str(),"mm") == 0)
    {
        MeshTools::Modification::scale(mesh, 0.1);
    }
    else
    {
        TBOX_ERROR("MESH_LENGTH_UNIT = " << mesh_length_unit << " is unrecognized" << endl);
    }

    // Center the mesh in the computational domain.
    //MeshTools::BoundingBox bbox = MeshTools::bounding_box(mesh);

    libMesh::BoundingBox bbox = MeshTools::create_bounding_box(mesh);


    libMesh::Point lower = bbox.min();
    libMesh::Point upper = bbox.max();
    pout << "         mesh bounding box lower = (" << lower(0) << " , " << lower(1) << " , " << lower(2) << ") cm\n"
    << "         mesh bounding box upper = (" << upper(0) << " , " << upper(1) << " , " << upper(2) << ") cm\n";
    //MeshTools::Modification::translate(mesh, -0.5*(upper(0)+lower(0)), -0.5*(upper(1)+lower(1)), -0.5*(upper(2)+lower(2)));
    //bbox = MeshTools::bounding_box(mesh);
    //lower = bbox.min();
    //upper = bbox.max();
    //pout << "centered mesh bounding box lower = (" << lower(0) << " , " << lower(1) << " , " << lower(2) << ") cm\n"
    //     << "centered mesh bounding box upper = (" << upper(0) << " , " << upper(1) << " , " << upper(2) << ") cm\n";


    // Setup ZERO_DISPLACEMENT boundary conditions.
    const MeshBase::const_element_iterator end_el = mesh.elements_end();
    for (MeshBase::const_element_iterator el = mesh.elements_begin(); el != end_el; ++el)
    {
        Elem* const elem = *el;
        for (unsigned int side = 0; side < elem->n_sides(); ++side)
        {
            const bool at_mesh_bdry = elem->neighbor_ptr(side) == nullptr;
            /*if (at_mesh_bdry)
             {
             const short int boundary_id = mesh.boundary_info->boundary_id(elem,side);
             if (boundary_id == 4000)
             {
             mesh.boundary_info->add_side(elem, side, FEDataManager::ZERO_DISPLACEMENT_XYZ_BDRY_ID);
             }
             }*/
        }
    }

    return;
}// initialize_mesh

void
ModelInitialization::initialize_equation_systems(EquationSystems* equation_systems)
{
    System& f0_system = equation_systems->add_system<System>("fiber direction");
     MechanicsModel::f0_system_num = f0_system.number();
     for (unsigned int d = 0; d < NDIM; ++d)
     {
     ostringstream os;
     os << "f0_" << d;
     f0_system.add_variable(os.str(), CONSTANT, MONOMIAL);
     }

     System& s0_system = equation_systems->add_system<System>("sheet direction");
     MechanicsModel::s0_system_num = s0_system.number();
     for (unsigned int d = 0; d < NDIM; ++d)
     {
         ostringstream os;
         os << "s0_" << d;
         s0_system.add_variable(os.str(), CONSTANT, MONOMIAL);
     }

    //adding the J system, determination of the deformation gradient
    // MechanicsModel::J_system_name = "J_detF"; //already initialized in MechanicalModel.C
    System& J_system = equation_systems->add_system<LinearImplicitSystem>(MechanicsModel::J_system_name);
    MechanicsModel::J_system_num = J_system.number();
    J_system.add_variable("J_detF", FIRST, LAGRANGE);

    System& Jold_system = equation_systems->add_system<LinearImplicitSystem>(MechanicsModel::Jold_system_name);
    MechanicsModel::Jold_system_num = Jold_system.number();
    Jold_system.add_variable("Jold", FIRST, LAGRANGE);

    System& dJdt_system = equation_systems->add_system<LinearImplicitSystem>(MechanicsModel::dJdt_system_name);
    MechanicsModel::dJdt_system_num = dJdt_system.number();
    dJdt_system.add_variable("dJ_dt", FIRST, LAGRANGE);

    // add the FF systems
    // MechanicsModel::FF_system_name = "FF"; //already initialized in MechanicalModel.C
    System& FF_system = equation_systems->add_system<LinearImplicitSystem>(MechanicsModel::FF_system_name);
    MechanicsModel::FF_system_num = FF_system.number();
    FF_system.add_variable("FF_00d", FIRST, LAGRANGE);
    FF_system.add_variable("FF_01d", FIRST, LAGRANGE);
    FF_system.add_variable("FF_02d", FIRST, LAGRANGE);
    FF_system.add_variable("FF_10d", FIRST, LAGRANGE);
    FF_system.add_variable("FF_11d", FIRST, LAGRANGE);
    FF_system.add_variable("FF_12d", FIRST, LAGRANGE);
    FF_system.add_variable("FF_20d", FIRST, LAGRANGE);
    FF_system.add_variable("FF_21d", FIRST, LAGRANGE);
    FF_system.add_variable("FF_22d", FIRST, LAGRANGE);

    //active contraction part
   if (!MechanicsModel::enable_active_tension)
   {
     pout << "activation model is not activated!\n";
     return;
   }

   // Build equation systems that store the active tension data and related
   // variables.
   System& T_system = equation_systems->add_system<System>("active tension");
   ActiveContraction::T_system_num = T_system.number();
   T_system.add_variable("T", CONSTANT, MONOMIAL);



   System& T_system_holder = equation_systems->add_system<System>("active_tension_holder");
   ActiveContraction::T_system_num_holder = T_system_holder.number();
   T_system_holder.add_variable("T_holder", CONSTANT, MONOMIAL);



   System& act_system = equation_systems->add_system<System>("activation variables");
   ActiveContraction::act_system_num = act_system.number();
   act_system.add_variable("Ca_i", CONSTANT, MONOMIAL);
   act_system.add_variable("Ca_b", CONSTANT, MONOMIAL);
   act_system.add_variable("Q1"  , CONSTANT, MONOMIAL);
   act_system.add_variable("Q2"  , CONSTANT, MONOMIAL);
   act_system.add_variable("Q3"  , CONSTANT, MONOMIAL);
   act_system.add_variable("z"   , CONSTANT, MONOMIAL);

    return;
}// initialize_equation_systems

void
ModelInitialization::initialize_material_axes(Mesh& mesh,
                                              EquationSystems* equation_systems,
                                              Pointer<Database> input_db)
{
  std::map<int,VectorValue<double> > f0_data;
  ifstream f0_stream(input_db->getString("FIBER_DIRECTION_FILENAME").c_str());
  for (unsigned int e = 0; e < mesh.n_elem(); ++e)
  {
      int e_idx;
      f0_stream >> e_idx;
      TBOX_ASSERT(e_idx == static_cast<int>(e+1));
      VectorValue<double>& f0 = f0_data[e];
      for (unsigned int d = 0; d < NDIM; ++d) f0_stream >> f0(d);
  }
  f0_stream.close();

  std::map<int,VectorValue<double> > s0_data;
  ifstream s0_stream(input_db->getString("SHEET_DIRECTION_FILENAME").c_str());
  for (unsigned int e = 0; e < mesh.n_elem(); ++e)
  {
      int e_idx;
      s0_stream >> e_idx;
      TBOX_ASSERT(e_idx == static_cast<int>(e+1));
      VectorValue<double>& s0 = s0_data[e];
      for (unsigned int d = 0; d < NDIM; ++d) s0_stream >> s0(d);
  }
  s0_stream.close();

  System& f0_system = equation_systems->get_system<System>("fiber direction");
  NumericVector<double>& f0_vec = *f0_system.solution;
  System& s0_system = equation_systems->get_system<System>("sheet direction");
  NumericVector<double>& s0_vec = *s0_system.solution;
  for (MeshBase::element_iterator it = mesh.elements_begin(); it != mesh.elements_end(); ++it)
  {
      Elem* elem = *it;
      const int elem_id = elem->id();

      TBOX_ASSERT(f0_data.find(elem_id) != f0_data.end());
      const VectorValue<double> f0 = f0_data[elem_id].unit();
      for (unsigned int d = 0; d < NDIM; ++d)
      {
          const int dof_index = elem->dof_number(MechanicsModel::f0_system_num,d,0);
          f0_vec.set(dof_index, f0(d));
      }

      TBOX_ASSERT(s0_data.find(elem_id) != s0_data.end());
      const VectorValue<double> s0 = s0_data[elem_id].unit();
      for (unsigned int d = 0; d < NDIM; ++d)
      {
          const int dof_index = elem->dof_number(MechanicsModel::s0_system_num,d,0);
          s0_vec.set(dof_index, s0(d));
      }
  }
  f0_vec.close();
  s0_vec.close();
  f0_vec.localize(*f0_system.current_local_solution);
  s0_vec.localize(*s0_system.current_local_solution);
  return;


    return;
}// initialize_material_axes
