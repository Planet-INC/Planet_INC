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

#ifndef PLANET_PLANET_PHYSICS_H
#define PLANET_PLANET_PHYSICS_H

// GRINS
#include "grins/physics.h"
#include "grins/var_typedefs.h"

namespace GRINS
{
  class AssemblyContext;
  class CachedValues;
}

// Planet
#include "planet/planet_physics_helper.h"

// libMesh
class GetPot;
namespace libMesh
{
  class FEMSystem;
}

namespace Planet
{

  template <typename CoeffType, typename VectorCoeffType>
  class PlanetPhysics : public GRINS::Physics
  {
  public:

    PlanetPhysics( const GRINS::PhysicsName& physics_name, const GetPot& input );

    virtual ~PlanetPhysics();

    //! Initialize variables for this physics.
    virtual void init_variables( libMesh::FEMSystem* system );

    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    //! Initialize context for added physics variables
    virtual void init_context( GRINS::AssemblyContext& context );

    //! Time dependent part(s) of physics for element interiors
    virtual void element_time_derivative( bool compute_jacobian,
                                          GRINS::AssemblyContext& context,
                                          GRINS::CachedValues& cache );

    virtual void mass_residual( bool compute_jacobian,
                                GRINS::AssemblyContext& context,
                                GRINS::CachedValues& cache );

  protected:

    //! Number of species
    unsigned int _n_species;

    //! Indices for each (owned) variable;
    std::vector<GRINS::VariableIndex> _species_vars; /* Indicies for species densities */

    //! Names of each (owned) variable in the system
    std::vector<std::string> _species_var_names;

    //! Element type, read from input
    libMeshEnums::FEFamily _species_FE_family;

    //! Element orders, read from input
    libMeshEnums::Order _species_order;

    PlanetPhysicsHelper<CoeffType,VectorCoeffType> _helper;

  private:

    PlanetPhysics();
    

  };

} // end namespace Planet

#endif // PLANET_PLANET_PHYSICS_H
