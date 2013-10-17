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

#ifndef PLANET_PHOTON_FLUX_H
#define PLANET_PHOTON_FLUX_H

//Planet
#include "planet/chapman.h"

//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/metaprogramming_decl.h"
#include "antioch/metaprogramming.h"
#include "antioch/cmath_shims.h"

//C++
#include <vector>
#include <map>

namespace Planet
{

//forward declarations
template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
class Atmosphere;

  template <typename CoeffType = double, 
            typename VectorCoeffType = std::vector<CoeffType>, 
            typename MatrixCoeffType = std::vector<std::vector<CoeffType> >
           >
  class PhotonFlux
  {
        private:
          PhotonFlux();
          void make_lambda_map(const VectorCoeffType &lambda);

          VectorCoeffType _phy_at_top;
          MatrixCoeffType _phy; //alt, lambda

          std::map<CoeffType,unsigned int> _map_lambda;
          std::map<unsigned int,CoeffType> _inverse_map_lambda;

          Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType> *_atm;
          Chapman<CoeffType> &_chap;
          
        public:
          PhotonFlux(Chapman<CoeffType> &chapman, Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType> *atm = NULL);
          ~PhotonFlux();

//!
          void set_photon_flux_top_atmosphere(const VectorCoeffType &lambda, const VectorCoeffType& phyAU, const CoeffType& dSunTopAtm);
//! update photon flux at all altitudes
          void update_photon_flux();
//! photon flux at altitude z
          VectorCoeffType phy(const CoeffType& z) const;
//! tau = Chap * sum_species sigma(lambda) int_z^top n_s(z')dz'
          VectorCoeffType tau(const CoeffType& z, const VectorCoeffType &sum_densities) const;

          const std::map<CoeffType,unsigned int> map_lambda()         const {return _map_lambda;}
          const std::map<unsigned int,CoeffType> inverse_map_lambda() const {return _inverse_map_lambda;}

          void set_atmosphere(Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType> *atm);
          void set_chapman(const Chapman<CoeffType> &chap);

          //!lambda vector
          const VectorCoeffType lambda() const;
          //!phy at top vector
          const VectorCoeffType phy_at_top() const;

  };

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  PhotonFlux<CoeffType,VectorCoeffType,MatrixCoeffType>::PhotonFlux(Chapman<CoeffType> &chapman, 
                                                                    Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType> *atm):
  _atm(atm),_chap(chapman)
  {
     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  PhotonFlux<CoeffType,VectorCoeffType,MatrixCoeffType>::PhotonFlux()
  {
     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  PhotonFlux<CoeffType,VectorCoeffType,MatrixCoeffType>::~PhotonFlux()
  {
     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const VectorCoeffType PhotonFlux<CoeffType,VectorCoeffType,MatrixCoeffType>::phy_at_top() const
  {
     return _phy_at_top;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const VectorCoeffType PhotonFlux<CoeffType,VectorCoeffType,MatrixCoeffType>::lambda() const
  {
     VectorCoeffType out_lambda;
     for(unsigned int i = 0; i < _inverse_map_lambda.size(); i++)
     {
        out_lambda.push_back(_inverse_map_lambda.at(i));
     }

     return out_lambda;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  VectorCoeffType PhotonFlux<CoeffType,VectorCoeffType,MatrixCoeffType>::tau(const CoeffType& z, const VectorCoeffType &sum_densities) const
  {
    VectorCoeffType tau;
    tau.resize(_map_lambda.size(),0.L);
    for(unsigned int ilambda = 0; ilambda < _map_lambda.size(); ilambda++) //lambda loop
    {
      for(unsigned int ns = 0; ns < _atm->n_photon_absorbing_species(); ns++)
      {
         tau[ilambda] += sum_densities[ns] * _atm->photon_sigma(ns).y_on_custom()[ilambda]; //filtering
      }
      tau[ilambda] *= _chap.chapman(_atm->a(z));
    }

    return tau;
  }


  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  VectorCoeffType PhotonFlux<CoeffType,VectorCoeffType,MatrixCoeffType>::phy(const CoeffType& z) const
  {
     antioch_assert_less_equal(_phy.size(),_atm->altitude_map().size());
     return _phy[_atm->altitude_map().at(z)];
  }


  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  void PhotonFlux<CoeffType,VectorCoeffType,MatrixCoeffType>::update_photon_flux()
  {

     antioch_assert_equal_to(_map_lambda.size(),_phy_at_top.size());
     _phy.resize(_atm->altitude_map().size());
     VectorCoeffType sum_densities;
     sum_densities.resize(_atm->n_neutral_species(),0.L);
//from top to bottom
     for(unsigned int nalt = 0; nalt < _atm->altitude_map().size(); nalt++)
     {
       _phy[nalt].resize(_map_lambda.size(),0.L);
//add altitude densities
       for(unsigned int nabs = 0; nabs < _atm->n_photon_absorbing_species(); nabs++)
       {
         sum_densities[nabs] += _atm->neutral_molar_density(nabs,nalt);
       }
       VectorCoeffType tau_on_the_fly = this->tau(_atm->altitude(nalt),sum_densities);

       for(unsigned int ilambda = 0; ilambda < _map_lambda.size(); ilambda++)
       {
          _phy[nalt][ilambda] = _phy_at_top[ilambda] * Antioch::ant_exp(- tau_on_the_fly[ilambda]);
       }
     }

     return; 
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  void PhotonFlux<CoeffType,VectorCoeffType,MatrixCoeffType>::set_atmosphere(Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType> *atm)
  {
     _atm = atm;
     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  void PhotonFlux<CoeffType,VectorCoeffType,MatrixCoeffType>::set_chapman(const Chapman<CoeffType> &chap)
  {
     _chap = chap;
     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  void PhotonFlux<CoeffType,VectorCoeffType,MatrixCoeffType>::make_lambda_map(const VectorCoeffType &lambda)
  {
     for(unsigned int i = 0; i < lambda.size(); i++)
     {
        _map_lambda[lambda[i]] = i;
        _inverse_map_lambda[i] = lambda[i];
     }

     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  void PhotonFlux<CoeffType,VectorCoeffType,MatrixCoeffType>::set_photon_flux_top_atmosphere(const VectorCoeffType &lambda, 
                                                                                             const VectorCoeffType& phyAU, 
                                                                                             const CoeffType& dSunTopAtm)
  {
     antioch_assert_equal_to(lambda.size(),phyAU.size());
     this->make_lambda_map(lambda);
     _phy_at_top.resize(phyAU.size(),0.L);
     for(unsigned int i = 0; i < phyAU.size(); i++)
     {
        _phy_at_top[i] = phyAU[i]/(dSunTopAtm * dSunTopAtm);
     }
  }

}

#endif
