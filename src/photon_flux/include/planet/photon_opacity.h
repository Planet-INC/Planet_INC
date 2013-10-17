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
#include "planet/atmospheric_mixture.h"

//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/cmath_shims.h"
#include "antioch/particle_flux.h"

//C++
#include <vector>
#include <map>

namespace Planet
{

  template <typename CoeffType = double, 
            typename MatrixCoeffType = std::vector<std::vector<CoeffType> >
           >
  class PhotonOpacity
  {
        private:
          PhotonOpacity(){antioch_error();return;}

          MatrixCoeffType _tau;

          //dependencies
          Chapman<CoeffType> &_chap;

        public:
          PhotonOpacity(Chapman<CoeffType> &chapman);
          ~PhotonOpacity();

          //! tau = Chap * sum_species sigma(lambda) int_z^top n_s(z')dz'
          template<typename StateType, typename VectorStateType, typename MatrixStateType>
          void update_tau(const StateType &a, const VectorStateType &totdens, const MatrixStateType &sigma);

          //!\returns photon opacity
          const MatrixCoeffType &tau() const;
  };

  template<typename CoeffType, typename MatrixCoeffType>
  inline
  PhotonOpacity<CoeffType,MatrixCoeffType>::PhotonOpacity(Chapman<CoeffType> &chapman):
  _chap(chapman)
  {
     return;
  }

  template<typename CoeffType, typename MatrixCoeffType>
  inline
  PhotonOpacity<CoeffType,MatrixCoeffType>::~PhotonOpacity()
  {
     return;
  }

  template<typename CoeffType, typename MatrixCoeffType>
  inline
  const MatrixCoeffType &PhotonOpacity<CoeffType,VectorCoeffType>::tau() const
  {
     return _tau;
  }



  template<typename CoeffType, typename MatrixCoeffType>
  template<typename StateType, typename MatrixStateType>
  inline
  void PhotonOpacity<CoeffType,VectorCoefftype>::update_tau(const StateType &a, const MatrixStateType &totdens, const MatrixStateType &sigma)
  {
      _tau.clear();
      _tau.resize(totdens[0].size());
      for(unsigned int iz = 0; iz < totdens[0].size(); iz++)
      {
        _tau[iz].resize(sigma[0].size(),0.L);
        for(unsigned int il = 0; il < sigma[0].size(); il++)
        {
          for(unsigned int s = 0; s < totdens.size(); s++)
          {
             _tau[iz][il] += sigma[s][ilambda] * totdens[s][iz];
          }
          _tau[iz][il] *= _chap.chapman(a);
        }
      }
      return;
  }


}

#endif
