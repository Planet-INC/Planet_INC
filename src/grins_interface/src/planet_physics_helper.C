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

//Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/vector_utils.h"

// This class
#include "planet/planet_physics_helper.h"

//C++

namespace Planet
{

  template<typename CoeffType, typename VectorCoeffType>
  PlanetPhysicsHelper<CoeffType,VectorCoeffType>::PlanetPhysicsHelper(AtmosphericKinetics<CoeffType,VectorCoeffType > *kinetics,
                                                                      DiffusionEvaluator <CoeffType,VectorCoeffType > *diffusion):
        _kinetics(kinetics),
        _diffusion(diffusion)
  {
    _omegas.resize(_kinetics->neutral_kinetics().reaction_set().n_species());
    _omegas_dots.resize(_kinetics->neutral_kinetics().reaction_set().n_species());
    return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  PlanetPhysicsHelper<CoeffType,VectorCoeffType>::~PlanetPhysicsHelper()
  {
    return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  libMesh::Real PlanetPhysicsHelper<CoeffType,VectorCoeffType>::diffusion_term(unsigned int s) const
  {
    return _omegas[s];
  }

  template<typename CoeffType, typename VectorCoeffType>
  libMesh::Real PlanetPhysicsHelper<CoeffType,VectorCoeffType>::chemical_term(unsigned int s) const
  {
    return _omegas_dots[s];
  }

} // end namespace Planet
