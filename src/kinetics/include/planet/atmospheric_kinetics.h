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

#ifndef PLANET_ATMOSPHERIC_KINETICS_H
#define PLANET_ATMOSPHERIC_KINETICS_H

//Antioch
#include "antioch/kinetics_evaluator.h"

//Planet
#include "planet/atmospheric_temperature.h"
#include "planet/atmospheric_mixture.h"
#include "planet/photon_evaluator.h"

//C++

namespace Planet
{
  template<typename CoeffType, typename VectorCoeffType>
  class AtmosphericKinetics
  {
      private:
        //! no default constructor
        AtmosphericKinetics() {antioch_error();return;}

        Antioch::KineticsEvaluator<CoeffType> &_neutral_reactions;
        Antioch::KineticsEvaluator<CoeffType> &_ionic_reactions;

//
        AtmosphericTemperature<CoeffType,VectorCoeffType> &_temperature;
        PhotonEvaluator<CoeffType,VectorCoeffType>        &_photon;
        AtmosphericMixture<CoeffType,VectorCoeffType>     &_composition;
      public:
        //!
        AtmosphericKinetics(Antioch::KineticsEvaluator<CoeffType>             &neu,
                            Antioch::KineticsEvaluator<CoeffType>             &ion,
                            AtmosphericTemperature<CoeffType,VectorCoeffType> &temperature,
                            PhotonEvaluator<CoeffType,VectorCoeffType>        &photon,
                            AtmosphericMixture<CoeffType,VectorCoeffType>     &composition);
        //!
        ~AtmosphericKinetics();

        //!\return neutral kinetics system, writable reference
        Antioch::KineticsEvaluator<CoeffType> &neutral_kinetics();

        //!\return ionic kinetics system, writable reference
        Antioch::KineticsEvaluator<CoeffType> &ionic_kinetics();

        //! compute chemical net rate and provide them in kin_rates
        template<typename StateType, typename VectorStateType>
        void chemical_rate(const VectorStateType &molar_concentrations, const VectorStateType &sum_concentrations, 
                           const StateType &z, VectorStateType &kin_rates) const;
  };


  template<typename CoeffType, typename VectorCoeffType>
  inline
  AtmosphericKinetics<CoeffType,VectorCoeffType>::AtmosphericKinetics(Antioch::KineticsEvaluator<CoeffType>             &neu,
                                                                      Antioch::KineticsEvaluator<CoeffType>             &ion,
                                                                      AtmosphericTemperature<CoeffType,VectorCoeffType> &temperature,
                                                                      PhotonEvaluator<CoeffType,VectorCoeffType>        &photon,
                                                                      AtmosphericMixture<CoeffType,VectorCoeffType>     &composition):
   _neutral_reactions(neu),
   _ionic_reactions(ion),
   _temperature(temperature),
   _photon(photon),
   _composition(composition)
  {
    return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  AtmosphericKinetics<CoeffType,VectorCoeffType>::~AtmosphericKinetics()
  {
    return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  Antioch::KineticsEvaluator<CoeffType> &AtmosphericKinetics<CoeffType,VectorCoeffType>::neutral_kinetics()
  {
     return _neutral_reactions;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  Antioch::KineticsEvaluator<CoeffType> &AtmosphericKinetics<CoeffType,VectorCoeffType>::ionic_kinetics()
  {
     return _ionic_reactions;
  }

  template<typename CoeffType, typename VectorCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void AtmosphericKinetics<CoeffType,VectorCoeffType>::chemical_rate(const VectorStateType &molar_concentrations, 
                                                                     const VectorStateType &sum_concentrations, 
                                                                     const StateType &z,
                                                                     VectorStateType &kin_rates) const
  {
     kin_rates.resize(_composition.neutral_composition().n_species(),0.L);
     VectorCoeffType dummy;
     dummy.resize(_composition.neutral_composition().n_species(),0.L); //everything is irreversible
     _photon.update_photon_flux(molar_concentrations, sum_concentrations, z);
     _neutral_reactions.compute_mole_sources(_temperature.neutral_temperature(z),
                                             molar_concentrations,dummy,kin_rates);
     return;
  }
}

#endif
