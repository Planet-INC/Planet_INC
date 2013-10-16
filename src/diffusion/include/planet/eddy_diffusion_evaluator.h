//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef _PLANET_EDDY_DIFFUSION_
#define _PLANET_EDDY_DIFFUSION_

//Antioch
#include "antioch/cmath_shims.h"

//Planet

//C++
#include <vector>

namespace Planet{

  template<typename CoeffType = double, typename VectorCoeffType = std::vector<double> >
  class EddyDiffusionEvaluator
  {
        private:
         CoeffType _K0;
         VectorCoeffType &_ntot;
         VectorCoeffType _K;

         EddyDiffusionEvaluator();

         template<typename VectorStateType>
         void make_eddy_diffusion();

         template<typename StateType>
         void reset_K0(const StateType &K0);

         template<typename VectorStateType>
         void reset_ntot(const VectorStateType &ntot);

         const VectorCoeffType K() const;

        public:
         EddyDiffusionEvaluator(const CoeffType K0, const VectorCoeffType &ntot);
         ~EddyDiffusionEvaluator();
  };

template<typename CoeffType, typename VectorCoeffType>
inline
EddyDiffusionEvaluator<CoeffType,VectorCoeffType>::EddyDiffusionEvaluator():
  _K0(-1.L),
{
  return;
}

template<typename CoeffType, typename VectorCoeffType>
inline
EddyDiffusionEvaluator<CoeffType,VectorCoeffType>::EddyDiffusionEvaluator(const CoeffType &K0, const VectorCoeffType &ntot):
  _K0(K0),
  _ntot(ntot)
{
  this->make_eddy_diffusion();
  return;
}

template<typename CoeffType, typename VectorCoeffType>
inline
EddyDiffusionEvaluator<CoeffType,VectorCoeffType>::~EddyDiffusionEvaluator()
{
  return;
}

template<typename CoeffType, typename VectorCoeffType>
template<typename StateType>
inline
void EddyDiffusionEvaluator<CoeffType,VectorCoeffType>::reset_K0(const StateType &K0)
{
   _K0 = K0;
   return;
}

template<typename CoeffType, typename VectorCoeffType>
template<typename VectorStateType>
inline
void EddyDiffusionEvaluator<CoeffType,VectorCoeffType>::reset_K0(const VectorStateType &ntot)
{
   _ntot = ntot;
   return;
}

template<typename CoeffType, typename VectorCoeffType>
inline
void EddyDiffusionEvaluator<CoeffType,VectorCoeffType>::make_eddy_diffusion()
{
//bottom to top
  antioch_assert_greater(_K0,0.)
  _K.resize(tot_dens.size,0.L);
  for(unsigned int ialt = 0; ialt < _ntot.size(); ialt++)
  {
     _K[ialt] = _K0 * Antioch::ant_sqrt(_ntot[ialt]/_ntot[0]);
  }
  return;
}

template<typename CoeffType, typename VectorCoeffType>
inline
const EddyDiffusionEvaluator<CoeffType,VectorCoeffType>::VectorCoeffType K() const
{
  return _K;
}


}

#endif
