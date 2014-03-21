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

#ifndef PLANET_PLANET_UPPER_NEUMANN_BC_H
#define PLANET_PLANET_UPPER_NEUMANN_BC_H

// GRINS
#include "grins/neumann_func_obj.h"
#include "grins/assembly_context.h"

// Planet
#include "planet/planet_physics_helper.h"

namespace Planet
{
  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  class PlanetUpperNeumannBC : public GRINS::NeumannFuncObj
  {
  public:

    PlanetUpperNeumannBC(const PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>& physics_helper,
                         const std::vector<GRINS::VariableIndex>& species_vars,
                         unsigned int species);

    virtual ~PlanetUpperNeumannBC();

    virtual libMesh::Real normal_value( const GRINS::AssemblyContext& context, const GRINS::CachedValues& cache,
					const unsigned int qp );

    virtual libMesh::Real normal_derivative( const GRINS::AssemblyContext& context, const GRINS::CachedValues& cache,
					     const unsigned int qp );

  protected:

    const PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>& _physics_helper;

    const std::vector<GRINS::VariableIndex> _species_vars;

    const unsigned int _species;

  private:

    PlanetUpperNeumannBC();

  };

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  PlanetUpperNeumannBC<CoeffType,VectorCoeffType,MatrixCoeffType>::PlanetUpperNeumannBC(const PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>& physics_helper,
                                                                                        const std::vector<GRINS::VariableIndex>& species_vars,
                                                                                        unsigned int species)
    : _physics_helper(physics_helper),
      _species_vars(species_vars),
      _species(species)
  {
    return;
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  PlanetUpperNeumannBC<CoeffType,VectorCoeffType,MatrixCoeffType>::~PlanetUpperNeumannBC()
  {
    return;
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  libMesh::Real PlanetUpperNeumannBC<CoeffType,VectorCoeffType,MatrixCoeffType>::normal_value( const GRINS::AssemblyContext& context,
                                                                                               const GRINS::CachedValues& /*cache*/,
                                                                                               const unsigned int qp )
  {
    unsigned int n_species = _physics_helper.composition().neutral_composition().n_species();

    std::vector<libMesh::Number> molar_densities(n_species, 0);
    for(unsigned int s=0; s < n_species; s++ )
      {
        molar_densities[s] = context.interior_value(this->_species_vars[s],qp);
      }

    return _physics_helper.upper_boundary_neumann(molar_densities,_species);
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  libMesh::Real PlanetUpperNeumannBC<CoeffType,VectorCoeffType,MatrixCoeffType>::normal_derivative( const GRINS::AssemblyContext& /*context*/,
                                                                                                    const GRINS::CachedValues& /*cache*/,
                                                                                                    const unsigned int /*qp*/ )
  {
    libmesh_not_implemented();

    // Dummy value
    return 0.0;
  }

} // end namespace Planet

#endif // PLANET_PLANET_UPPER_NEUMANN_BC_H
