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

  class PlanetPhysicsHelper
  {
  public:

    PlanetPhysicsHelper(AtmosphericKinetics<double,std::vector<double>,std::vector<std::vector<double> > > *kinetics = NULL,
                        DiffusionEvaluator <double,std::vector<double>,std::vector<std::vector<double> > > *diffusion = NULL);

    ~PlanetPhysicsHelper();

    void set_kinetics(AtmosphericKinetics<double,std::vector<double>,std::vector<std::vector<double> > > *kinetics);

    void set_diffusion(DiffusionEvaluator <double,std::vector<double>,std::vector<std::vector<double> > > *diffusion);

    libMesh::Real compute_omega(unsigned int s, double z, const std::vector<double> & molar_concentrations,
                                                          const std::vector<double> & dmolar_concentrations_dz);

    libMesh::Real compute_omega_dot(unsigned int s, double z, const std::vector<double> & molar_concentrations,
                                                              const std::vector<double> & dmolar_concentrations_dz);

  private:

    void compute(const std::vector<double> & molar_concentrations,
                 const std::vector<double> & dmolar_concentrations_dz,
                 double z);

    AtmosphericKinetics<double,std::vector<double>,std::vector<std::vector<double> > > *_kinetics;
    DiffusionEvaluator <double,std::vector<double>,std::vector<std::vector<double> > > *_diffusion;

    std::vector<double> _omegas;
    std::vector<double> _omegas_dots;
    double _current_z;
    const double _eps;

  };

} // end namespace Planet

#endif // PLANET_PLANET_PHYSICS_HELPER_H
