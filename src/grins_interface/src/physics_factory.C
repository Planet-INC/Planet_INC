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

// This class
#include "planet/physics_factory.h"

// Planet
#include "planet/physics_names.h"
#include "planet/planet_physics.h"

namespace Planet
{

  PhysicsFactory::PhysicsFactory()
    : GRINS::PhysicsFactory()
  {
    return;
  }

  PhysicsFactory::~PhysicsFactory()
  {
    return;
  }

  void PhysicsFactory::add_physics( const GetPot& input,
                                    const std::string& physics_to_add,
                                    GRINS::PhysicsList& physics_list )
  {
    // Convenience
    typedef std::tr1::shared_ptr<GRINS::Physics> PhysicsPtr;

    // Deal with creating PlanetPhysics
    if( physics_to_add == planet_physics )
      {
        /* For now, just using standard types. Can generalize later to other types
           from input. */
        physics_list[physics_to_add] =
          PhysicsPtr( new PlanetPhysics<double,std::vector<double>, std::vector<std::vector<double> > >(physics_to_add,input) );
      }
    // Call base class otherwise
    else
      {
        GRINS::PhysicsFactory::add_physics(input,physics_to_add,physics_list);
      }

    return;
  }

} // end namespace Planet
