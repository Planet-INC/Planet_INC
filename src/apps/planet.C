//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Planet - An atmospheric code for planetary bodies, adapted to Titan
//
// Copyright (C) 2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

// GRINS
#include "grins/simulation.h"
#include "grins/simulation_builder.h"

// Planet
#include "planet/physics_factory.h"

int main(int argc, char* argv[])
{
  // libMesh input file should be first argument
  std::string libMesh_input_filename = argv[1];

  // Create our GetPot object.
  GetPot libMesh_inputfile( libMesh_input_filename );

  // Initialize libMesh library.
  LibMeshInit libmesh_init(argc, argv);

  GRINS::SimulationBuilder sim_builder;

  // PhysicsFactory handles which GRINS::Physics objects to create
  std::tr1::shared_ptr<GRINS::PhysicsFactory> physics_factory( new Planet::PhysicsFactory );
  sim_builder.attach_physics_factory(physics_factory);

  GRINS::Simulation grins( libMesh_inputfile,
                           sim_builder,
                           libmesh_init.comm() );

  // Do solve here
  grins.run();
  
  return 0;
}
