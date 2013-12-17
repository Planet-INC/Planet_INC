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

  PlanetPhysicsHelper::PlanetPhysicsHelper(AtmosphericKinetics<double,std::vector<double>,std::vector<std::vector<double> > > *kinetics,
                                           DiffusionEvaluator <double,std::vector<double>,std::vector<std::vector<double> > > *diffusion):
        _kinetics(kinetics),
        _diffusion(diffusion),
        _current_z(-1.),
        _eps(std::numeric_limits<double>::epsilon() * 1000.L)
  {
    _omegas.resize(_kinetics->neutral_kinetics().reaction_set().n_species());
    _omegas_dots.resize(_kinetics->neutral_kinetics().reaction_set().n_species());
    return;
  }

  PlanetPhysicsHelper::~PlanetPhysicsHelper()
  {
    return;
  }

  libMesh::Real PlanetPhysicsHelper::compute_omega(unsigned int s, double z,
                                                   const std::vector<double> & molar_concentrations,
                                                   const std::vector<double> & dmolar_concentrations_dz)
  {
    this->compute(molar_concentrations,dmolar_concentrations_dz,z);
    return _omegas[s];
  }

  libMesh::Real PlanetPhysicsHelper::compute_omega_dot(unsigned int s, double z,
                                                       const std::vector<double> & molar_concentrations,
                                                       const std::vector<double> & dmolar_concentrations_dz)
  {
    this->compute(molar_concentrations,dmolar_concentrations_dz,z);
    return _omegas_dots[s];
  }

  void PlanetPhysicsHelper::compute(const std::vector<double> & molar_concentrations,
                                    const std::vector<double> & dmolar_concentrations_dz,
                                    double z)
  {
    if((z - _current_z) > _eps)
    {
     _diffusion->diffusion(molar_concentrations,dmolar_concentrations_dz,z,_omegas);
     _kinetics->chemical_rate(molar_concentrations,z,_omegas_dots);
     _current_z = z;
    }
    return;
  }

  void PlanetPhysicsHelper::set_kinetics(AtmosphericKinetics<double,std::vector<double>,std::vector<std::vector<double> > > *kinetics)
  {
     _kinetics = kinetics;
  }

  void PlanetPhysicsHelper::set_diffusion(DiffusionEvaluator <double,std::vector<double>,std::vector<std::vector<double> > > *diffusion)
  {
     _diffusion = diffusion;
  }

} // end namespace Planet
