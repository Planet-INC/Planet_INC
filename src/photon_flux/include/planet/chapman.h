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

#ifndef PLANET_CHAPMAN_FUNCTION_H
#define PLANET_CHAPMAN_FUNCTION_H

//Antioch
#include "antioch/metaprogramming_decl.h"
#include "antioch/cmath_shims.h"
#include "antioch/antioch_asserts.h"
//
//Planet
#include "planet/math_constants.h"

//C++
#include <cmath>

namespace Planet{

template <typename CoeffType = double>
class Chapman{
   private:

     CoeffType _chi;

     //! degrees to gradians
     template <typename StateType = CoeffType>
     ANTIOCH_AUTO(StateType)
     rad(const StateType &chi) const
     ANTIOCH_AUTOFUNC(StateType,chi*Constants::pi<StateType>()/StateType(180.L))

     //! gradians to degrees
     template <typename StateType = CoeffType>
     ANTIOCH_AUTO(StateType)
     deg(const StateType &chi) const
     ANTIOCH_AUTOFUNC(StateType,chi*StateType(180.L)/Constants::pi<StateType>())

     //! chapman for high angles
     template <typename StateType>
     ANTIOCH_AUTO(StateType)
     chapman_high_angles(const StateType &x) const
     ANTIOCH_AUTOGENFUNC(StateType,antioch_assert_greater(_chi,rad(90.L));antioch_assert_less(_chi,rad(180.L)),
                                Antioch::ant_sqrt(StateType(2.L) * Constants::pi<StateType>() * x) * 
                                (  Antioch::ant_sqrt(Antioch::ant_sin(_chi)) * 
                                   Antioch::ant_exp( x * (StateType(1.L) - Antioch::ant_sin(_chi))) -
                                   StateType(0.5L) *
                                   Antioch::ant_exp( x/StateType(2.L) *
                                            Antioch::ant_pow(Antioch::ant_cos(_chi),2) * 
                                            (StateType(1.L) - this->erf(
                                                                        Antioch::ant_sqrt( x / StateType(2.L)) * 
                                                                        Antioch::ant_abs(Antioch::ant_cos(_chi))
                                                                       )
                                            )
                                                   )
                                )
                        )
     //! chapman for medium angles
     template <typename StateType>
     ANTIOCH_AUTO(StateType)
     chapman_medium_angles(const StateType &x) const
     ANTIOCH_AUTOGENFUNC(StateType,antioch_assert_less(_chi,rad(90.1L));antioch_assert_greater(_chi,rad(75.L)),
                                Antioch::ant_sqrt(Constants::pi<StateType>() * x/StateType(2.L)) * 
                                (StateType(1.L) - this->erf(Antioch::ant_sqrt(x/StateType(2.L))  * Antioch::ant_abs(Antioch::ant_cos(_chi)))) *
                                Antioch::ant_exp(x/StateType(2.L) * Antioch::ant_pow(Antioch::ant_cos(_chi),2))
                        )
     
     //! The approximation from Abramowitz and Stegun, Eq. 7.1.26
     template <typename StateType>
     ANTIOCH_AUTO(StateType)
     erf(StateType x) const
     ANTIOCH_AUTOGENFUNC(StateType, if(x < 0.)x = -x ,
                StateType(1.L) - (StateType( 0.254829592L)  * (StateType(1.L)/(StateType(1.L)                 + StateType(0.3275911L) * x))   +
                                  StateType(-0.284496736L)  * Antioch::ant_pow(StateType(1.L)/(StateType(1.L) + StateType(0.3275911L) * x),2) + 
                                  StateType( 1.421413741L)  * Antioch::ant_pow(StateType(1.L)/(StateType(1.L) + StateType(0.3275911L) * x),3) + 
                                  StateType(-1.453152027L)  * Antioch::ant_pow(StateType(1.L)/(StateType(1.L) + StateType(0.3275911L) * x),4) + 
                                  StateType( 1.061405429L)  * Antioch::ant_pow(StateType(1.L)/(StateType(1.L) + StateType(0.3275911L) * x),5)
                                 ) * Antioch::ant_exp(-Antioch::ant_pow(x,2))
                        )
   public:

     Chapman(){return;}
     Chapman(const CoeffType &chi):_chi(rad(chi)){return;}
     ~Chapman(){return;}

     //!
     template <typename StateType>
     void set_chi(const StateType &chi) {_chi = rad(chi);}

     //!
     CoeffType chapman() const ;
     //!
     template <typename StateType>
     StateType chapman(const StateType &x) const; // x = (R + z)/H

     //! \return the rate evaluated at low angle.
     CoeffType operator()() const;

     //! \return the rate evaluated at \p a.
     template <typename StateType>
     ANTIOCH_AUTO(StateType) 
     operator()(const StateType& a) const
     ANTIOCH_AUTOFUNC(StateType, this->chapman(a))

     //! Approach angle
     CoeffType chi() const;
};

template<typename CoeffType>
inline
CoeffType Chapman<CoeffType>::operator()() const
{
   return this->chapman();
}

template<typename CoeffType>
inline
CoeffType Chapman<CoeffType>::chapman() const
{
   antioch_assert_less(_chi,rad(75.1L));
   return CoeffType(1.L)/Antioch::ant_cos(_chi);
}

template<typename CoeffType>
template<typename StateType>
inline
StateType Chapman<CoeffType>::chapman(const StateType & x) const
{
  if(_chi < rad(75.1L))return this->chapman();
  if(_chi < rad(90.1L))return this->chapman_medium_angles(x);
                       return this->chapman_high_angles(x);
}

template <typename CoeffType>
inline
CoeffType Chapman<CoeffType>::chi() const
{
  return _chi;
}

}

#endif
