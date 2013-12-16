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

//Planet
#include "planet/altitude.h"

//C++

namespace Planet
{
  template<typename CoeffType, typename VectorCoeffType>
  class AtmosphericTemperature
  {
      private:
        //!no default constructor
        AtmosphericTemperature(){antioch_error();return;}

        VectorCoeffType _neutral_temperature;
        VectorCoeffType _ionic_temperature;
        VectorCoeffType _electronic_temperature;


///dependencie
        Altitude<CoeffType,VectorCoeffType> &_altitude;

      public:
        AtmosphericTemperature(const VectorCoeffType &neu, const VectorCoeffType &ion, Altitude<CoeffType,VectorCoeffType> &alt);
        ~AtmosphericTemperature();

        //!\return neutral temperature
        const VectorCoeffType &neutral_temperature() const;

        //!\return ionic temperature
        const VectorCoeffType &ionic_temperature()   const;

        //!\return electronic temperature
        const VectorCoeffType &electronic_temperature()   const;

        //!
        template<typename VectorStateType>
        void set_neutral_temperature(const VectorStateType &neu);

        //!
        template<typename VectorStateType>
        void set_ionic_temperature(const VectorStateType &ion);

        template<typename VectorStateType>
        void set_electronic_temperature(const VectorStateType &electron);

        //!
        void initialize();

  };

  template<typename CoeffType, typename VectorCoeffType>
  inline
  AtmosphericTemperature<CoeffType,VectorCoeffType>::~AtmosphericTemperature()
  {
    return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  AtmosphericTemperature<CoeffType,VectorCoeffType>::AtmosphericTemperature(const VectorCoeffType &neu, 
                                                                            const VectorCoeffType &ion, 
                                                                            Altitude<CoeffType,VectorCoeffType> &alt):
      _neutral_temperature(neu),
      _ionic_temperature(ion),
      _altitude(alt)
  {
    _electronic_temperature.resize(_altitude.altitudes().size());
    for(unsigned int iz = 0; iz < _altitude.altitudes().size(); iz++)
    {
      if(_altitude.altitudes()[iz] < 900.)
      {
         _electronic_temperature[iz] = 180.L;
      }else if(_altitude.altitudes()[iz] < 1400.)
      {
         _electronic_temperature[iz] = CoeffType(180.L) + (_altitude.altitudes()[iz] - CoeffType(900.L))/CoeffType(500.L) * CoeffType(1150.L - 180.L) ;
      }else
      {
         _electronic_temperature[iz] = CoeffType(1150.L) + CoeffType(0.1L) * (_altitude.altitudes()[iz] - CoeffType(1400.L));
      }
    }
    return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  void AtmosphericTemperature<CoeffType,VectorCoeffType>::initialize()
  {
     return;//nothing right now
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  const VectorCoeffType &AtmosphericTemperature<CoeffType,VectorCoeffType>::ionic_temperature() const
  {
    return _ionic_temperature;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  const VectorCoeffType &AtmosphericTemperature<CoeffType,VectorCoeffType>::neutral_temperature() const
  {
    return _neutral_temperature;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  const VectorCoeffType &AtmosphericTemperature<CoeffType,VectorCoeffType>::electronic_temperature() const
  {
    return _electronic_temperature;
  }


  template<typename CoeffType, typename VectorCoeffType>
  template<typename VectorStateType>
  inline
  void AtmosphericTemperature<CoeffType,VectorCoeffType>::set_neutral_temperature(const VectorStateType &neu)
  {
     _neutral_temperature = neu;
  }

  template<typename CoeffType, typename VectorCoeffType>
  template<typename VectorStateType>
  inline
  void AtmosphericTemperature<CoeffType,VectorCoeffType>::set_ionic_temperature(const VectorStateType &ion)
  {
     _ionic_temperature = ion;
  }

  template<typename CoeffType, typename VectorCoeffType>
  template<typename VectorStateType>
  inline
  void AtmosphericTemperature<CoeffType,VectorCoeffType>::set_electronic_temperature(const VectorStateType &electron)
  {
     _electronic_temperature = electron;
  }


}

#endif
