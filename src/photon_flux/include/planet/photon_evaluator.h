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

#ifndef PLANET_PHOTON_EVALUATOR_H
#define PLANET_PHOTON_EVALUATOR_H

//Antioch
#include "antioch/particle_flux.h"
#include "antioch/reaction_set.h"

//Planet
#include "planet/photon_opacity.h"
#include "planet/atmospheric_mixture.h"

//C++
#include <vector>
#include <map>

namespace Planet
{
  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  class PhotonEvaluator
  {
     private:
        //! no default constructor authorized
        PhotonEvaluator(){antioch_error();return;}

//dependencies
        const Antioch::ParticleFlux<VectorCoeffType>                        & _phy_at_top;
        const PhotonOpacity<CoeffType,VectorCoeffType>                      & _hv_tau;
        const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> & _mixture;

     public:
        PhotonEvaluator(const Antioch::ParticleFlux<VectorCoeffType>                        &phy_at_top,
                        const PhotonOpacity<CoeffType,VectorCoeffType>                      &hv_tau, 
                        const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &mix);
        ~PhotonEvaluator();

        //!\return photon flux at top of atmosphere
        const Antioch::ParticleFlux<VectorCoeffType> &photon_flux_at_top() const;

        //!calculate photon flux
        template<typename StateType, typename VectorStateType>
        void update_photon_flux(const VectorStateType &molar_densities, const VectorStateType &sum_dens, 
                                const StateType &z, VectorStateType & flux_at_z) const;

  };

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::PhotonEvaluator(
                                                              const Antioch::ParticleFlux<VectorCoeffType>                        &phy_at_top,
                                                              const PhotonOpacity<CoeffType,VectorCoeffType>                      &hv_tau, 
                                                              const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &mix):
  _phy_at_top(phy_at_top),
  _hv_tau(hv_tau),
  _mixture(mix)
  {
     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::~PhotonEvaluator()
  {
     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::update_photon_flux(const VectorStateType &molar_densities, 
                                                                      const VectorStateType &sum_dens, const StateType &z,
                                                                      VectorStateType & flux_at_z) const
  {
     antioch_assert_equal_to(molar_densities.size(), _mixture.neutral_composition().n_species());
     antioch_assert_equal_to(sum_dens.size(), _mixture.neutral_composition().n_species());
     antioch_assert(!_phy_at_top.abscissa().empty());
     antioch_assert(!_phy_at_top.flux().empty());
     antioch_assert_equal_to(_phy_at_top.flux().size(),flux_at_z.size());

     VectorCoeffType tau;
     _hv_tau.compute_tau(_mixture.a(molar_densities,z),sum_dens,tau);

     antioch_assert_equal_to(tau.size(), _phy_at_top.abscissa().size());

     for(unsigned int ilambda = 0; ilambda < _phy_at_top.abscissa().size(); ilambda++)
     {
       flux_at_z[ilambda] = _phy_at_top.flux()[ilambda] * Antioch::ant_exp(- tau[ilambda]);
     }


     return; 
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const Antioch::ParticleFlux<VectorCoeffType> & PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::photon_flux_at_top() const
  {
     return _phy_at_top;
  }
  

}

#endif
