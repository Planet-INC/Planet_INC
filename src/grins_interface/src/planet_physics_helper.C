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

// This class
#include "planet/planet_physics_helper.h"

#include "antioch/vector_utils.h"

//C++
#include <limits>

namespace Planet
{

  template<typename CoeffType, typename VectorCoeffType>
  PlanetPhysicsHelper<CoeffType,VectorCoeffType>::PlanetPhysicsHelper(AtmosphericKinetics<CoeffType,VectorCoeffType > *kinetics,
                                                                      DiffusionEvaluator <CoeffType,VectorCoeffType > *diffusion):
        _kinetics(kinetics),
        _diffusion(diffusion),
        _current_z(-1.),
        _eps(std::numeric_limits<double>::epsilon() * 1000.L)
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
  libMesh::Real PlanetPhysicsHelper<CoeffType,VectorCoeffType>::diffusion_term(unsigned int s)
  {
    return _omegas[s];
  }

  template<typename CoeffType, typename VectorCoeffType>
  libMesh::Real PlanetPhysicsHelper<CoeffType,VectorCoeffType>::chemical_term(unsigned int s)
  {
    return _omegas_dots[s];
  }

  template<typename CoeffType, typename VectorCoeffType>
  template<typename StateType, typename VectorStateType, typename MatrixStateType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType>::compute(const VectorStateType & molar_concentrations,
                                                               const VectorStateType & dmolar_concentrations_dz,
                                                               const VectorStateType & other_altitudes,
                                                               const MatrixStateType & other_concentrations,
                                                               const StateType & z)
  {
    if((z - _current_z) > _eps)
    {
//// if other_altitudes from top to bottom, end at z
     std::vector<double> sum_concentration;
     sum_concentration.resize(molar_concentrations.size(),0.L);
     for(unsigned int iz = 0; iz < other_altitudes.size() - 1; iz++)
     {
        for(unsigned int s = 0; s < sum_concentration.size(); s++)
        {
          sum_concentration[s] += other_concentrations[s][iz] * (other_altitudes[iz] - other_altitudes[iz +1]);
        }
     }
     _diffusion->diffusion(molar_concentrations,dmolar_concentrations_dz,z,_omegas);
     _kinetics->chemical_rate(molar_concentrations,sum_concentration,z,_omegas_dots);
     _current_z = z;
    }
    return;
  }

  template <typename CoeffType, typename VectorCoeffType>
  template <typename StateType, typename VectorStateType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType>::set_kinetics(AtmosphericKinetics<StateType,VectorStateType> *kinetics)
  {
     _kinetics = kinetics;
  }

  template <typename CoeffType, typename VectorCoeffType>
  template <typename StateType, typename VectorStateType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType>::set_diffusion(DiffusionEvaluator <StateType,VectorStateType> *diffusion)
  {
     _diffusion = diffusion;
  }

} // end namespace Planet
