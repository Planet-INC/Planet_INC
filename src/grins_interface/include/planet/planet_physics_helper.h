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

#ifndef PLANET_PLANET_PHYSICS_HELPER_H
#define PLANET_PLANET_PHYSICS_HELPER_H

//Planet
#include "planet/diffusion_evaluator.h"
#include "planet/atmospheric_kinetics.h"

// libMesh
#include "libmesh/libmesh_common.h"

namespace Planet
{

  template<typename CoeffType, typename VectorCoeffType>
  class PlanetPhysicsHelper
  {
  public:

    PlanetPhysicsHelper(AtmosphericKinetics<CoeffType,VectorCoeffType > *kinetics = NULL,
                        DiffusionEvaluator <CoeffType,VectorCoeffType > *diffusion = NULL);

    ~PlanetPhysicsHelper();

    template <typename StateType, typename VectorStateType>
    void set_kinetics(AtmosphericKinetics<StateType,VectorStateType> *kinetics);

    template <typename StateType, typename VectorStateType>
    void set_diffusion(DiffusionEvaluator<StateType,VectorStateType> *diffusion);

    libMesh::Real diffusion_term(unsigned int s);

    libMesh::Real chemical_term(unsigned int s);

    template<typename StateType, typename VectorStateType, typename MatrixStateType>
    void compute(const VectorStateType & molar_concentrations,
                 const VectorStateType & dmolar_concentrations_dz,
                 const VectorStateType & other_altitudes,
                 const MatrixStateType & other_concentrations,
                 const StateType & z);

  private:

    AtmosphericKinetics<CoeffType,VectorCoeffType> *_kinetics;
    DiffusionEvaluator <CoeffType,VectorCoeffType> *_diffusion;

    VectorCoeffType _omegas;
    VectorCoeffType _omegas_dots;
    CoeffType _current_z;
    const CoeffType _eps;

  };

} // end namespace Planet

#endif // PLANET_PLANET_PHYSICS_HELPER_H
