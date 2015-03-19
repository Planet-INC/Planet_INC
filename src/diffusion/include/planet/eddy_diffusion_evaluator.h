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
#include "planet/atmospheric_mixture.h"

//C++
#include <vector>

namespace Planet{

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  class EddyDiffusionEvaluator
  {
        private:
         CoeffType _K0;

//dependencies
         const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &_mixture;

         EddyDiffusionEvaluator() {antioch_error();return;}


        public:

         //!\return K0
         const CoeffType K0() const;

         //!sets K0
         template<typename StateType>
         void set_K0(const StateType &K0);

         //! \return eddy coefficient in cm2.s-1
         template<typename StateType>
         ANTIOCH_AUTO(StateType)
         K(const StateType &ntot) const
         ANTIOCH_AUTOFUNC(StateType,_K0 * Antioch::ant_sqrt(_mixture.total_bottom_density()/ntot))

         //! \return deddy_dT coefficient in cm2.s-1.K-1
         template<typename StateType>
         ANTIOCH_AUTO(StateType)
         K_deriv_T(const StateType &T) const
         ANTIOCH_AUTOFUNC(StateType,Antioch::zero_clone(T))

         //! \return deddy_dn coefficient in cm2.s-1.K-1.(cm-3)-1
         template<typename StateType>
         ANTIOCH_AUTO(StateType)
         K_deriv_ns(const StateType &ntot) const
         ANTIOCH_AUTOFUNC(StateType,- this->K(ntot) / (ntot * Antioch::constant_clone(ntot,2)))


         //!
         EddyDiffusionEvaluator(const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &mix, 
                                const CoeffType &K0 = -1.);

         //!
         ~EddyDiffusionEvaluator();
  };

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
EddyDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::EddyDiffusionEvaluator(const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &mix, 
                                                                          const CoeffType &K0):
  _K0(K0),
  _mixture(mix)
{
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
   return;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
const CoeffType EddyDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::K0() const
{
  return _K0;
}

}

#endif
