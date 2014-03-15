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

//parameters & output
        Antioch::ParticleFlux<VectorCoeffType>   _phy_at_top;
        Antioch::ParticleFlux<VectorCoeffType> * _phy;

//dependencies
        PhotonOpacity<CoeffType,VectorCoeffType>      &_hv_tau;
        const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &_mixture;

     public:
        PhotonEvaluator(PhotonOpacity<CoeffType,VectorCoeffType> &hv_tau, 
                        const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &mix);
        ~PhotonEvaluator();

        //!\return const ref photon flux
        const Antioch::ParticleFlux<VectorCoeffType> &photon_flux() const;

        //!\return photon flux pointer at top of atmosphere
        Antioch::ParticleFlux<VectorCoeffType> * photon_flux_ptr();

        //!\return photon flux at top of atmosphere
        const Antioch::ParticleFlux<VectorCoeffType> &photon_flux_at_top() const;

        //!sets the photon flux at the top of the atmosphere
        template<typename StateType, typename VectorStateType>
        void set_photon_flux_at_top(const VectorStateType &lambda, const VectorStateType &hv, const StateType &d = 1.L);

        //!calculate photon flux
        template<typename StateType, typename VectorStateType>
        void update_photon_flux(const VectorStateType &molar_densities, const VectorStateType &sum_dens, const StateType &z);

  };

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::PhotonEvaluator(PhotonOpacity<CoeffType,VectorCoeffType> &hv_tau, 
                                                              const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &mix):
  _phy(NULL),
  _hv_tau(hv_tau),
  _mixture(mix)
  {
     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::~PhotonEvaluator()
  {
     if(_phy)delete _phy;
     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::set_photon_flux_at_top(const VectorStateType &lambda, const VectorStateType &hv, const StateType &d)
  {
     antioch_assert_equal_to(lambda.size(),hv.size());

//phy at top
     _phy_at_top.set_abscissa(lambda);
     VectorCoeffType flux;
     flux.resize(hv.size());
     for(unsigned int i = 0; i < hv.size(); i++)
     {
        flux[i] = hv[i]/(d * d);
     }
     _phy_at_top.set_flux(flux);

//cross sections
     _hv_tau.update_cross_section(lambda);

     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const Antioch::ParticleFlux<VectorCoeffType> &PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::photon_flux_at_top() const
  {
    return _phy_at_top;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const Antioch::ParticleFlux<VectorCoeffType> &PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::photon_flux() const
  {
     return *_phy;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  Antioch::ParticleFlux<VectorCoeffType>  *PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::photon_flux_ptr()
  {
     return _phy;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::update_photon_flux(const VectorStateType &molar_densities, 
                                                                      const VectorStateType &sum_dens, const StateType &z)
  {
     antioch_assert_equal_to(molar_densities.size(), _mixture.neutral_composition().n_species());
     antioch_assert_equal_to(sum_dens.size(), _mixture.neutral_composition().n_species());
     antioch_assert(!_phy_at_top.abscissa().empty());
     antioch_assert(!_phy_at_top.flux().empty());

     if(!_phy)
     {
       _phy = new Antioch::ParticleFlux<VectorCoeffType>;
       _phy->set_abscissa(_phy_at_top.abscissa());
     }

     VectorCoeffType tau;
     _hv_tau.compute_tau(_mixture.a(molar_densities,z),sum_dens,tau);

     antioch_assert_equal_to(tau.size(), _phy_at_top.abscissa().size());

     VectorCoeffType flux;
     flux.resize(_phy_at_top.abscissa().size());
     for(unsigned int ilambda = 0; ilambda < _phy_at_top.abscissa().size(); ilambda++)
     {
       flux[ilambda] = _phy_at_top.flux()[ilambda] * Antioch::ant_exp(- tau[ilambda]);
     }

     _phy->set_flux(flux);

     return; 
  }

}

#endif
