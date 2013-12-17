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
#include "planet/altitude.h"
#include "planet/atmospheric_temperature.h"
#include "planet/atmospheric_mixture.h"
#include "planet/photon_evaluator.h"

//C++

namespace Planet
{
  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  class AtmosphericKinetics
  {
      private:
        //! no default constructor
        AtmosphericKinetics() {antioch_error();return;}

        Antioch::KineticsEvaluator<CoeffType> &_neutral_reactions;
        Antioch::KineticsEvaluator<CoeffType> &_ionic_reactions;

        MatrixCoeffType _P_minus_nL;
        MatrixCoeffType _chemical_loss; 
        MatrixCoeffType _chemical_production;

//
        Altitude<CoeffType,VectorCoeffType> &_altitude;
        AtmosphericTemperature<CoeffType,VectorCoeffType> &_temperature;
        PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> &_photon;
        AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &_composition;
      public:
        //!
        AtmosphericKinetics(Antioch::KineticsEvaluator<CoeffType> &neu, 
                            Antioch::KineticsEvaluator<CoeffType> &ion,
                            Altitude<CoeffType,VectorCoeffType> &altitude,
                            AtmosphericTemperature<CoeffType,VectorCoeffType> &temperature,
                            PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> &photon,
                            AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &composition);
        //!
        ~AtmosphericKinetics();

        //!\return neutral kinetics system, writable reference
        Antioch::KineticsEvaluator<CoeffType> &neutral_kinetics();

        //!\return ionic kinetics system, writable reference
        Antioch::KineticsEvaluator<CoeffType> &ionic_kinetics();

        //!
        void initialize();

        //! computes P and L
        void compute_production_and_loss_rates();

        //! computes the total rates, r = P - n L
        void compute_P_minus_nL();

        //!\return chemical production
        const MatrixCoeffType &chemical_production() const;

        //!\return chemical loss
        const MatrixCoeffType &chemical_loss() const;

        //!\return chemical net rate
        const MatrixCoeffType &chemical_net_rate() const;

        //! compute chemical net rate and provide them in kin_rates
        template<typename StateType, typename VectorStateType>
        void chemical_rate(const VectorStateType &molar_concentrations, const StateType &z,
                           VectorStateType &kin_rates) const;
  };


  //!please be aware that here, the convention is [iz][s]
  //for chemical production and loss, for portability to
  //lower level calculations, the total rate is with the
  //normal convention [s][iz]
  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  AtmosphericKinetics<CoeffType,VectorCoeffType,MatrixCoeffType>::AtmosphericKinetics(Antioch::KineticsEvaluator<CoeffType> &neu,
                                                      Antioch::KineticsEvaluator<CoeffType> &ion,
                                                      Altitude<CoeffType,VectorCoeffType> &altitude,
                                                      AtmosphericTemperature<CoeffType,VectorCoeffType> &temperature,
                                                      PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> &photon,
                                                      AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &composition):
   _neutral_reactions(neu),
   _ionic_reactions(ion),
   _altitude(altitude),
   _temperature(temperature),
   _photon(photon),
   _composition(composition)
  {
    _chemical_loss.resize(_altitude.altitudes().size());
    _chemical_production.resize(_altitude.altitudes().size());
    for(unsigned int iz = 0; iz < _altitude.altitudes().size(); iz++)
    {
       _chemical_loss[iz].resize(_composition.neutral_composition().n_species());
       _chemical_production[iz].resize(_composition.neutral_composition().n_species());
    }

    _P_minus_nL.resize(_composition.neutral_composition().n_species());
    for(unsigned int s = 0; s < _composition.neutral_composition().n_species(); s++)
    {
      _P_minus_nL[s].resize(_altitude.altitudes().size());
    }


    this->initialize();
    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  AtmosphericKinetics<CoeffType,VectorCoeffType,MatrixCoeffType>::~AtmosphericKinetics()
  {
    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  void AtmosphericKinetics<CoeffType,VectorCoeffType,MatrixCoeffType>::initialize()
  {
    this->compute_production_and_loss_rates();
    this->compute_P_minus_nL();
  }


  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  void AtmosphericKinetics<CoeffType,VectorCoeffType,MatrixCoeffType>::compute_production_and_loss_rates()
  {
    VectorCoeffType dummy_rev;
    dummy_rev.resize(_composition.neutral_composition().n_species());
    for(unsigned int iz = 0; iz < _altitude.altitudes().size(); iz++)
    {
       VectorCoeffType molar_densities;
       molar_densities.resize(_composition.neutral_composition().n_species());
       for(unsigned int s = 0; s < _composition.neutral_composition().n_species(); s++)
       {
          molar_densities[s] = _composition.total_density()[iz] * _composition.neutral_molar_fraction()[s][iz];
       }
       VectorCoeffType net_rate,kfwd_const,kbkwd_const,kfwd,kbkwd,fwd_conc,bkwd_conc;

       //setting the photon flux
       _photon.set_photon_flux(iz);
        //getting reaction infos
       _neutral_reactions.reaction_set().get_reactive_scheme(_temperature.neutral_temperature()[iz],
                                                             molar_densities,
                                                             dummy_rev, //all reactions are irreversible
                                                             net_rate,kfwd_const,kbkwd_const,kfwd,kbkwd,fwd_conc,bkwd_conc);

      //loss & prod rates (P and nL)
      for(unsigned int rxn = 0; rxn < _neutral_reactions.reaction_set().n_reactions(); rxn++)
      {
        const Antioch::Reaction<CoeffType>& reaction = _neutral_reactions.reaction_set().reaction(rxn);
        for (unsigned int r = 0; r < reaction.n_reactants(); r++)
        {
           _chemical_loss[iz][reaction.reactant_id(r)] += static_cast<CoeffType>(reaction.reactant_stoichiometric_coefficient(r)) * kfwd[rxn];
        }
        for (unsigned int p = 0; p < reaction.n_products(); p++)
        {
           _chemical_production[iz][reaction.product_id(p)] += static_cast<CoeffType>(reaction.product_stoichiometric_coefficient(p)) * kfwd[rxn];
        }
      }

    }
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  void AtmosphericKinetics<CoeffType,VectorCoeffType,MatrixCoeffType>::compute_P_minus_nL()
  {
     for(unsigned int iz = 0; iz < _altitude.altitudes().size(); iz++)
     {
       for(unsigned int s = 0; s < _composition.neutral_composition().n_species(); s++)
       {
         _P_minus_nL[s][iz] = _chemical_production[iz][s] - _chemical_loss[iz][s];
       }
     }
     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  Antioch::KineticsEvaluator<CoeffType> &AtmosphericKinetics<CoeffType,VectorCoeffType,MatrixCoeffType>::neutral_kinetics()
  {
     return _neutral_reactions;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  Antioch::KineticsEvaluator<CoeffType> &AtmosphericKinetics<CoeffType,VectorCoeffType,MatrixCoeffType>::ionic_kinetics()
  {
     return _ionic_reactions;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void AtmosphericKinetics<CoeffType,VectorCoeffType,MatrixCoeffType>::chemical_rate(const VectorStateType &molar_concentrations, const StateType &z,
                                                                                     VectorStateType &kin_rates) const
  {
     kin_rates.resize(_composition.neutral_composition().n_species(),0.L);
     VectorCoeffType dummy;
     dummy.resize(_composition.neutral_composition().n_species(),0.L); //everything is irreversible
     _photon.set_photon_flux(z);
     _neutral_reactions.compute_mole_sources(_temperature.neutral_temperature(z),
                                             molar_concentrations,dummy,kin_rates);
     return;
  }
}

#endif
