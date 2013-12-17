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

#ifndef _PLANET_EDDY_DIFFUSION_
#define _PLANET_EDDY_DIFFUSION_

//Antioch
#include "antioch/cmath_shims.h"

//Planet
#include "planet/altitude.h"
#include "planet/atmospheric_mixture.h"

//C++
#include <vector>

namespace Planet{

  template<typename CoeffType = double, 
           typename VectorCoeffType = std::vector<double>, 
           typename MatrixCoeffType = std::vector<std::vector<double> >
           >
  class EddyDiffusionEvaluator
  {
        private:
         CoeffType _K0;
         CoeffType _ntot_bottom;
         VectorCoeffType _K;

//dependencies
         AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &_mixture;
         Altitude<CoeffType,VectorCoeffType> &_altitude;

         EddyDiffusionEvaluator() {antioch_error();return;}


        public:
         //!\return eddy coefficient vector
         const VectorCoeffType& K() const;

         //!\return K0
         const CoeffType K0() const;

         //!calculate eddy coefficient
         void make_eddy_diffusion();

         //!sets K0
         template<typename StateType>
         void set_K0(const StateType &K0);

         template<typename StateType>
         ANTIOCH_AUTO(StateType)
         K(const StateType &ntot) const
         ANTIOCH_AUTOFUNC(StateType,_K0 * Antioch::ant_sqrt(_ntot_bottom/ntot))

         //!
         EddyDiffusionEvaluator(AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &mix, 
                                Altitude<CoeffType,VectorCoeffType> &alt,
                                const CoeffType &K0 = -1.);

         //!
         ~EddyDiffusionEvaluator();
  };

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
EddyDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::EddyDiffusionEvaluator(AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &mix, 
                                                                          Altitude<CoeffType,VectorCoeffType> &alt,
                                                                          const CoeffType &K0):
  _K0(K0),
  _mixture(mix),
  _altitude(alt)
{
  _ntot_bottom = _mixture.total_density()[_altitude.altitudes_map().at(_altitude.alt_min())];
  return;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
EddyDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::~EddyDiffusionEvaluator()
{
  return;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
template<typename StateType>
inline
void EddyDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::set_K0(const StateType &K0)
{
   _K0 = K0;
   this->make_eddy_diffusion();
   return;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
void EddyDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::make_eddy_diffusion()
{

  antioch_assert_greater(_K0,0.);
  _K.resize(_altitude.altitudes().size(),0.L);
  for(unsigned int ialt = 0; ialt < _altitude.altitudes().size(); ialt++)
  {
     _K[ialt] = this->K(_mixture.total_density()[ialt]);
  }
  return;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
const VectorCoeffType &EddyDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::K() const
{
  return _K;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
const CoeffType EddyDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::K0() const
{
  return _K0;
}

}

#endif
