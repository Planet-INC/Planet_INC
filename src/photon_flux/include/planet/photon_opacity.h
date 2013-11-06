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
#include "planet/altitude.h"

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
            typename VectorCoeffType = std::vector<CoeffType>, 
            typename MatrixCoeffType = std::vector<VectorCoeffType> 
           >
  class PhotonOpacity
  {
        private:
          PhotonOpacity(){antioch_error();return;}

          MatrixCoeffType _tau;

          //dependencies
          Altitude<CoeffType,VectorCoeffType> &_altitude;
          Chapman<CoeffType> &_chapman;

        public:
          PhotonOpacity(Altitude<CoeffType,VectorCoeffType> &alt, Chapman<CoeffType> &chapman);
          ~PhotonOpacity();

          //! tau = Chap * sum_species sigma(lambda) int_z^top n_s(z')dz'
          template<typename VectorStateType, typename MatrixStateType>
          void update_tau(const VectorStateType &a, const MatrixStateType &sumdens, const MatrixStateType &sigma);

          //!\returns photon opacity
          const MatrixCoeffType &tau() const;
  };

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  PhotonOpacity<CoeffType,VectorCoeffType,MatrixCoeffType>::PhotonOpacity(Altitude<CoeffType,VectorCoeffType> &alt,
                                                          Chapman<CoeffType> &chapman):
  _altitude(alt),
  _chapman(chapman)
  {
     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  PhotonOpacity<CoeffType,VectorCoeffType,MatrixCoeffType>::~PhotonOpacity()
  {
     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const MatrixCoeffType &PhotonOpacity<CoeffType,VectorCoeffType,MatrixCoeffType>::tau() const
  {
     return _tau;
  }



  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename VectorStateType, typename MatrixStateType>
  inline
  void PhotonOpacity<CoeffType,VectorCoeffType,MatrixCoeffType>::update_tau(const VectorStateType &a, const MatrixStateType &sumdens, const MatrixStateType &sigma)
  {
      antioch_assert_equal_to(_altitude.altitudes().size(),a.size());
      _tau.clear();
      _tau.resize(_altitude.altitudes().size());
      for(unsigned int iz = 0; iz < _altitude.altitudes().size(); iz++) //alt
      {
        _tau[iz].resize(sigma[0].size(),0.L);
        for(unsigned int il = 0; il < sigma[0].size(); il++) //lambda
        {
          for(unsigned int s = 0; s < sumdens.size(); s++) // neutrals
          {
             _tau[iz][il] += sigma[s][il] * sumdens[s][iz];
          }
          _tau[iz][il] *= _chapman(a[iz]) * _altitude.alt_step() * 1e5; //km -> cm
        }
      }
      return;
  }


}

#endif
