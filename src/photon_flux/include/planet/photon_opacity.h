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

#ifndef PLANET_PHOTON_OPACITY_H
#define PLANET_PHOTON_OPACITY_H

//Planet
#include "planet/chapman.h"
#include "planet/cross_section.h"

//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/cmath_shims.h"

//C++
#include <vector>
#include <map>

namespace Planet
{

  template <typename CoeffType, typename VectorCoeffType>
  class PhotonOpacity
  {
        private:
          PhotonOpacity(){antioch_error();return;}

          //dependencies
          Chapman<CoeffType> &_chapman;

//store
          std::map<unsigned int, unsigned int>    _cross_sections_map;
          std::vector<CrossSection<VectorCoeffType> > _absorbing_species_cs;
          std::vector<unsigned int>               _absorbing_species;
          std::vector<unsigned int>                   _absorbing_species_id;

        public:
          PhotonOpacity(Chapman<CoeffType> &chapman);
          ~PhotonOpacity();

          //! tau = Chap * sum_species sigma(lambda) int_z^top n_s(z')dz'
          template<typename StateType, typename VectorStateType>
          void compute_tau(const StateType &a, const VectorStateType &sum_dens, VectorStateType &tau) const;

          //!\return absorbing species
          const std::vector<unsigned int> absorbing_species() const;

          //!\return absorbing species cross-section
          const std::vector<CrossSection<VectorCoeffType> > &absorbing_species_cs() const;

          //!\return absorbing species cross-section map
          const std::map<unsigned int, unsigned int> cross_sections_map() const;

          //!adds a photon cross-section
          template<typename VectorStateType>
          void add_cross_section(const VectorStateType &lambda, const VectorStateType &cs, const unsigned int &sp, unsigned int id);

          //!update cross-section
          template<typename VectorStateType>
          void update_cross_section(const VectorStateType &custom_grid);


  };

  template<typename CoeffType, typename VectorCoeffType>
  inline
  PhotonOpacity<CoeffType,VectorCoeffType>::PhotonOpacity(Chapman<CoeffType> &chapman):
  _chapman(chapman)
  {
     return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  PhotonOpacity<CoeffType,VectorCoeffType>::~PhotonOpacity()
  {
     return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  template<typename VectorStateType>
  inline
  void PhotonOpacity<CoeffType,VectorCoeffType>::add_cross_section(const VectorStateType &lambda, const VectorStateType &cs, 
                                                                   const unsigned int &sp, unsigned int id)
  {
     _absorbing_species.push_back(sp);
     _absorbing_species_id.push_back(id);
     _absorbing_species_cs.push_back(CrossSection<VectorCoeffType>(lambda,cs));
     _cross_sections_map[sp] = _absorbing_species.size() - 1;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  const std::vector<unsigned int> PhotonOpacity<CoeffType,VectorCoeffType>::absorbing_species() const
  {
     return _absorbing_species;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  const std::vector<CrossSection<VectorCoeffType> > &PhotonOpacity<CoeffType,VectorCoeffType>::absorbing_species_cs() const
  {
      return _absorbing_species_cs;
  }


  template<typename CoeffType, typename VectorCoeffType>
  inline
  const std::map<unsigned int, unsigned int> PhotonOpacity<CoeffType,VectorCoeffType>::cross_sections_map() const
  {
     return _cross_sections_map;
  }

  template<typename CoeffType, typename VectorCoeffType>
  template<typename VectorStateType>
  inline
  void PhotonOpacity<CoeffType,VectorCoeffType>::update_cross_section(const VectorStateType &custom_grid)
  {
     for(unsigned int i = 0; i < _absorbing_species.size(); i++)
     {
        _absorbing_species_cs[i].update_cross_section(custom_grid);
     }

     return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void PhotonOpacity<CoeffType,VectorCoeffType>::compute_tau(const StateType &a, const VectorStateType &sum_dens, VectorStateType &tau) const
  {
      antioch_assert(!_absorbing_species_cs.empty());

      tau.resize(_absorbing_species_cs[0].cross_section_on_custom_grid().size(),0.L);

      for(unsigned int il = 0; il < _absorbing_species_cs[0].cross_section_on_custom_grid().size(); il++) //lambda
      {
          for(unsigned int s = 0; s < _absorbing_species_cs.size(); s++) // neutrals
          {
             tau[il] += _absorbing_species_cs[s].cross_section_on_custom_grid()[il] * sum_dens[_absorbing_species_id[s]]; //cm2 * cm-3.km
          }
          tau[il] *= _chapman(a) * Antioch::constant_clone(a,1e5); //cm-1.km  -> no unit
      }
      return;
  }


}

#endif
