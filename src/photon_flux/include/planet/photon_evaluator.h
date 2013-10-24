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
#include "antioch/species_enum.h"
#include "antioch/particle_flux.h"

//Planet
#include "planet/photon_opacity.h"
#include "planet/cross_section.h"
#include "planet/atmospheric_mixture.h"
#include "planet/photon_evaluator.h"

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
        Antioch::ParticleFlux<VectorCoeffType> _phy_at_top;
        std::vector<Antioch::ParticleFlux<VectorCoeffType> > _phy; //alt

//store
        std::map<Antioch::Species, unsigned int> _cross_sections_map;
        std::vector<CrossSection<VectorCoeffType> > _absorbing_species_cs;
        std::vector<Antioch::Species> _absorbing_species;

//dependencies
        PhotonOpacity<CoeffType,VectorCoeffType> &_hv_tau;
        Altitude<CoeffType,VectorCoeffType> &_altitude;
        AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &_mixture;

        //! cross-section on flux grid
        template<typename VectorStateType>
        void update_cross_section(const VectorStateType &custom_grid);

     public:
        PhotonEvaluator(Altitude<CoeffType,VectorCoeffType> &alt,
                        PhotonOpacity<CoeffType,VectorCoeffType> &hv_tau, 
                        AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &mix);
        ~PhotonEvaluator();

        //!sets the photon flux at the top of the atmosphere
        template<typename StateType, typename VectorStateType>
        void set_photon_flux_at_top(const VectorStateType &lambda, const VectorStateType &hv, const StateType &d = 1.L);

        //!\return photon flux
        const std::vector<Antioch::ParticleFlux<VectorCoeffType> > &photon_flux() const;

        //!calculate photon flux
        void update_photon_flux();

        //!\return absorbing species
        const std::vector<Antioch::Species> absorbing_species() const;

        //!\return absorbing species cross-section map
        const std::vector<CrossSection<VectorCoeffType> > &absorbing_species_cs() const;

        //!adds a photon cross-section
        template<typename VectorStateType>
        void add_cross_section(const VectorStateType &lambda, const VectorStateType &cs, const Antioch::Species &sp);

        //!
        void initialize();

  };

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::PhotonEvaluator(Altitude<CoeffType,VectorCoeffType> &alt,
                                                              PhotonOpacity<CoeffType,VectorCoeffType> &hv_tau, 
                                                              AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &mix):
  _hv_tau(hv_tau),
  _altitude(alt),
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
  inline
  void PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::initialize()
  {
    this->update_photon_flux();
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
     for(unsigned int i = 0; i < hv.size(); i++)
     {
        flux.push_back(hv[i]/(d * d));
     }
     _phy_at_top.set_flux(flux);

//phy at all altitudes
     _phy.clear();
     _phy.resize(_altitude.altitudes().size());
     for(unsigned int iz = 0; iz < _phy.size(); iz++)
     {
        _phy[iz].set_abscissa(lambda);
     }
//cross sections
     this->update_cross_section(lambda);

     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename VectorStateType>
  inline
  void PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::add_cross_section(const VectorStateType &lambda, const VectorStateType &cs, const Antioch::Species &sp)
  {
     _absorbing_species.push_back(sp);
     _absorbing_species_cs.push_back(CrossSection<VectorCoeffType>(lambda,cs));
     _cross_sections_map[sp] = _absorbing_species.size() - 1;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const std::vector<Antioch::ParticleFlux<VectorCoeffType> > &PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::photon_flux() const
  {
     return _phy;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const std::vector<Antioch::Species> PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::absorbing_species() const
  {
     return _absorbing_species;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const std::vector<CrossSection<VectorCoeffType> > &PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::absorbing_species_cs() const
  {
      return _absorbing_species_cs;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename VectorStateType>
  inline
  void PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::update_cross_section(const VectorStateType &custom_grid)
  {
     for(unsigned int i = 0; i < _absorbing_species.size(); i++)
     {
        _absorbing_species_cs[i].update_cross_section(custom_grid);
     }

     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  void PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::update_photon_flux()
  {
     antioch_assert(!_phy_at_top.abscissa().empty());
     antioch_assert(!_phy_at_top.flux().empty());

     MatrixCoeffType sigma;
     sigma.resize(_absorbing_species.size());
     for(unsigned int s = 0; s < sigma.size(); s++)
     {
        sigma[s] =  _absorbing_species_cs[s].cross_section_on_custom_grid();
     }

//
     MatrixCoeffType totdens;
     totdens.resize(_absorbing_species.size());
     for(unsigned int s = 0; s < _absorbing_species.size(); s++)
     {
       totdens[s].resize(_altitude.altitudes().size(),0.L);
//from top to bottom
       unsigned int iz = _altitude.altitudes_map().at(_altitude.alt_max());
       totdens[s][iz]  = _mixture.neutral_molar_fraction()[s][iz] * _mixture.total_density()[iz];
       unsigned int izb = iz;
       for(CoeffType z = _altitude.alt_max() - _altitude.alt_step(); z >= _altitude.alt_min(); z -= _altitude.alt_step())
       {
         iz = _altitude.altitudes_map().at(z);
         totdens[s][iz] = totdens[s][izb] + _mixture.neutral_molar_fraction()[s][iz] * _mixture.total_density()[iz];
         izb = iz;
       }
     }

//tau
     _hv_tau.update_tau(_mixture.a_factor(), totdens, sigma);

//finally
     for(unsigned int iz = 0; iz < _altitude.altitudes().size(); iz++)
     {
       VectorCoeffType flux;
       flux.resize(_phy_at_top.abscissa().size());
       for(unsigned int ilambda = 0; ilambda < _phy_at_top.abscissa().size(); ilambda++)
       {
         flux[ilambda] = _phy_at_top.flux()[ilambda] * Antioch::ant_exp(- _hv_tau.tau()[iz][ilambda]);
       }
       _phy[iz].set_flux(flux);
     }

     return; 
  }

}

#endif
