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

#ifndef PLANET_ATMOSPHERE_TEMPERATURE_H
#define PLANET_ATMOSPHERE_TEMPERATURE_H

//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/gsl_spliner.h"

//Planet

//C++


/*before it is merged in Antioch*/
// GSL
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

namespace Planet
{

  /*! \class AtmosphericTemperature
      
      This class contains the models for three temperatures:
        - the neutral species temperature
        - the ionic species temperature
        - the electrons temperature

      We use a spline for the neutral temperature, it is read from
      an input file.
      We use the hypothesis \f$T_\text{ion} = T_\text{neu}\f$.
      We use the simple Titan model:
      \f[
           \begin{split}
             T_\text{e} = 180~\text{K}                             & \text{ if }         z  < 900~\text{km}\\
             T_\text{e} = 180 + \frac{z - 900}{500} * (1150 - 180) & \text{ if } 900 \le z < 1400~\text{km} \\
             T_\text{e} = 1150 + 0.1 * (z - 1400)                  & \text{ else}\\
           \end{split}
      \f]
   */
  template<typename CoeffType, typename VectorCoeffType, typename Spliner = Antioch::GSLSpliner>
  class AtmosphericTemperature
  {
      private:
        //!no default constructor
        AtmosphericTemperature(){antioch_error();return;}

        Spliner spline;

      public:
        AtmosphericTemperature(const VectorCoeffType &alt_neu, const VectorCoeffType &T_neu);
        ~AtmosphericTemperature();

        //!\return neutral temperature at custom altitude
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
        neutral_temperature(const StateType &z) const
        ANTIOCH_AUTOFUNC(StateType, spline.interpolated_value(z))

        //!\return derivative neutral temperature at custom altitude
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
        dneutral_temperature_dz(const StateType &z) const
        ANTIOCH_AUTOFUNC(StateType, spline.dinterp_dx(z))

        //!\return ionic temperature at custom altitude
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
        ionic_temperature(const StateType &z)   const
        ANTIOCH_AUTOFUNC(StateType, neutral_temperature(z))

        //!\return derivative ionic temperature at custom altitude
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
        dionic_temperature_dz(const StateType &z)   const
        ANTIOCH_AUTOFUNC(StateType, dneutral_temperature_dz(z))

        //!\return electronic temperature at custom altitude
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
        electronic_temperature(const StateType &z)   const
        ANTIOCH_AUTOFUNC(StateType, (z < Antioch::constant_clone(z,900.))?Antioch::constant_clone(z,180.):
                                         (z < Antioch::constant_clone(z,1400.))?Antioch::constant_clone(z,180.) + Antioch::constant_clone(z,2e-3) * (z - Antioch::constant_clone(z,900.)) * Antioch::constant_clone(z,(1150. - 180.)):
                                                                                Antioch::constant_clone(z,1150.) + Antioch::constant_clone(z,0.1) * (z - Antioch::constant_clone(z,1400.))
                        )

        //!\return derivative electronic temperature at custom altitude
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
        delectronic_temperature_dz(const StateType &z)   const
        ANTIOCH_AUTOFUNC(StateType, (z < Antioch::constant_clone(z,900.))?Antioch::zero_clone(z):
                                         (z < Antioch::constant_clone(z,1400.))?Antioch::constant_clone(z,2e-3):
                                                                                Antioch::constant_clone(z,0.1)
                        )

        //! reset neutral temperature
        void set_neutral_temperature(const VectorCoeffType & alt, const VectorCoeffType &neu);

        //!
        void initialize();

  };

  template<typename CoeffType, typename VectorCoeffType, typename Spliner>
  inline
  AtmosphericTemperature<CoeffType,VectorCoeffType,Spliner>::~AtmosphericTemperature()
  {
    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename Spliner>
  inline
  AtmosphericTemperature<CoeffType,VectorCoeffType,Spliner>::AtmosphericTemperature(const VectorCoeffType &alt_neu, 
                                                                            const VectorCoeffType &T_neu)
        :spline(alt_neu,T_neu)
  {
    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename Spliner>
  inline
  void AtmosphericTemperature<CoeffType,VectorCoeffType,Spliner>::initialize()
  {
     return; //nothing
  }


  template<typename CoeffType, typename VectorCoeffType, typename Spliner>
  inline
  void AtmosphericTemperature<CoeffType,VectorCoeffType,Spliner>::set_neutral_temperature(const VectorCoeffType & alt_neu, const VectorCoeffType & T_neu)
  {
    spline.spline_delete();
    spline.init(alt_neu,T_neu);
  }

}

#endif
