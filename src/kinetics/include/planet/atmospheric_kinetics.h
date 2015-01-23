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
#include "antioch/kinetics_conditions.h"
#include "antioch/kinetics_evaluator.h"

//Planet
#include "planet/atmospheric_temperature.h"
#include "planet/atmospheric_mixture.h"
#include "planet/photon_evaluator.h"
#include "planet/atmospheric_steady_state.h"

//eigen
#include <Eigen/Dense>

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

        bool _ionic_coupling;
        std::vector<Antioch::Species>                     _ions_species;
        AtmosphericSteadyState<CoeffType,VectorCoeffType> _newton_solver;
        
//
        const AtmosphericTemperature<CoeffType,VectorCoeffType>             & _temperature;
        PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>          & _photon;
        const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> & _composition;
      public:
        //!
        AtmosphericKinetics(Antioch::KineticsEvaluator<CoeffType>                               &neu,
                            Antioch::KineticsEvaluator<CoeffType>                               &ion,
                            const AtmosphericTemperature<CoeffType,VectorCoeffType>             &temperature,
                            PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>          &photon,
                            const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &composition,
                            std::vector<Antioch::Species>                                       ionic_species = std::vector<Antioch::Species>());
        //!
        ~AtmosphericKinetics();

        //!\return neutral kinetics system, writable reference
        Antioch::KineticsEvaluator<CoeffType> &neutral_kinetics();

        //!\return ionic kinetics system, writable reference
        Antioch::KineticsEvaluator<CoeffType> &ionic_kinetics();

        //! compute chemical net rate and provide them in kin_rates
        template<typename StateType, typename VectorStateType>
        void chemical_rate(const VectorStateType &molar_concentrations, 
                           const Antioch::KineticsConditions<StateType> &KC, 
                           const StateType & z, VectorStateType &kin_rates);

        template<typename StateType, typename VectorStateType, typename MatrixStateType>
        void chemical_rate_and_derivs(const VectorStateType &molar_concentrations,
                                      const Antioch::KineticsConditions<StateType> &KC, 
                                      const StateType & z,
                                      VectorStateType &kin_rates, MatrixStateType &dkin_rates_dn);

        //! Newton solver for the ionic system
        template<typename StateType, typename VectorStateType>
        void add_ionic_contribution(const VectorStateType &molar_concentrations, const Antioch::KineticsConditions<StateType> &KC, 
                                    const StateType & z, VectorStateType &kin_rates);

        //! Newton solver for the ionic system
        template<typename StateType, typename VectorStateType, typename MatrixStateType>
        void add_ionic_contribution_and_derivs(const VectorStateType &neutral_concentrations, 
                                               const Antioch::KineticsConditions<StateType> & KC, const StateType &z, 
                                               VectorStateType &kin_rates, MatrixStateType &dkin_rates_dn);
  };


  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  AtmosphericKinetics<CoeffType,VectorCoeffType,MatrixCoeffType>::AtmosphericKinetics(Antioch::KineticsEvaluator<CoeffType>               &neu,
                                                                      Antioch::KineticsEvaluator<CoeffType>                               &ion,
                                                                      const AtmosphericTemperature<CoeffType,VectorCoeffType>             &temperature,
                                                                      PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>          &photon,
                                                                      const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &composition,
                                                                      std::vector<Antioch::Species>                                       ionic_species ):
   _neutral_reactions(neu),
   _ionic_reactions(ion),
   _ions_species(ionic_species),
   _newton_solver(_ions_species,ion),
   _temperature(temperature),
   _photon(photon),
   _composition(composition)
  {
    _ionic_coupling = !ionic_species.empty();
    if(_ionic_coupling)
    {
      _newton_solver.build_map();        //ionospheric solver maps
    }
    
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
  void AtmosphericKinetics<CoeffType,VectorCoeffType,MatrixCoeffType>::chemical_rate(const VectorStateType &molar_concentrations, 
                                                                     const Antioch::KineticsConditions<StateType> &KC,
                                                                     const StateType & z,
                                                                     VectorStateType &kin_rates)
  {
     antioch_assert_equal_to(kin_rates.size(),_composition.neutral_composition().n_species());
     VectorStateType dummy;
     dummy.resize(_composition.neutral_composition().n_species(),0.L); //everything is irreversible
     _neutral_reactions.compute_mole_sources(KC,
                                             molar_concentrations,dummy,kin_rates);

     if(_ionic_coupling)this->add_ionic_contribution(molar_concentrations,KC,z,kin_rates);

     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType, typename MatrixStateType>
  inline
  void AtmosphericKinetics<CoeffType,VectorCoeffType,MatrixCoeffType>::chemical_rate_and_derivs(const VectorStateType &molar_concentrations,
                                      const Antioch::KineticsConditions<StateType> & kinetics_conditions, 
                                      const StateType & z, VectorStateType &kin_rates, MatrixStateType &dkin_rates_dn)
  {
     antioch_assert_equal_to(kin_rates.size(),_composition.neutral_composition().n_species());
     antioch_assert_equal_to(dkin_rates_dn.size(),_composition.neutral_composition().n_species());
#ifdef NDEBUG
#else
     for(unsigned int s = 0; s < _composition.neutral_composition().n_species(); s++)
     {
       antioch_assert_equal_to(dkin_rates_dn[s].size(),_composition.neutral_composition().n_species());
     }
#endif
     VectorStateType dummy;
     VectorStateType ddummy_dT;
     VectorStateType dkin_dT;
     dummy.resize(_composition.neutral_composition().n_species(),0.L); //everything is irreversible
     ddummy_dT.resize(_composition.neutral_composition().n_species(),0.L); //everything is irreversible
     dkin_dT.resize(_composition.neutral_composition().n_species(),0.L); //no temp
     _neutral_reactions.compute_mole_sources_and_derivs(kinetics_conditions,molar_concentrations,dummy,ddummy_dT,
                                                        kin_rates,dkin_dT,dkin_rates_dn);

     if(_ionic_coupling)this->add_ionic_contribution_and_derivs(molar_concentrations,kinetics_conditions,z,kin_rates,dkin_rates_dn);
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void AtmosphericKinetics<CoeffType,VectorCoeffType,MatrixCoeffType>::add_ionic_contribution(const VectorStateType &neutral_concentrations, 
                                                                                              const Antioch::KineticsConditions<StateType> &KC, 
                                                                                              const StateType & z,
                                                                                              VectorStateType &kin_rates)
  {

//    if(z < 800. || z > 1200.)return;

 // neutrals resized
    VectorStateType full_concentrations;
    full_concentrations.resize(_composition.ionic_composition().n_species(),0.L);
    for(unsigned int s = 0; s < neutral_concentrations.size(); s++)
    {
       unsigned int i = _composition.ionic_composition().species_list()[_composition.neutral_composition().species_list()[s]];
       full_concentrations[i] = neutral_concentrations[s];
    }

//solve for ions
    VectorCoeffType source_ions;
    source_ions.resize(_ionic_reactions.n_species(),0.L);
//all temperature conditions, solver deal with it

    _newton_solver.precompute_rates(full_concentrations,KC, KC.T(), _temperature.electronic_temperature(z)); 
    if(_newton_solver.steady_state(source_ions))
    {

//     std::cout << "Ionospheric activity at " << z << " km" << std::endl;
// update sources
      for(unsigned int s = 0; s < _composition.neutral_composition().n_species(); s++)
      {
        unsigned int i_neu = _composition.ionic_composition().species_list()[_composition.neutral_composition().species_list()[s]];
        kin_rates[s] += source_ions[i_neu];
      }
    }
  }


  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType, typename MatrixStateType>
  inline
  void AtmosphericKinetics<CoeffType,VectorCoeffType,MatrixCoeffType>::add_ionic_contribution_and_derivs(const VectorStateType &neutral_concentrations, 
                                                                                                         const Antioch::KineticsConditions<StateType> & KC, const StateType &z, 
                                                                                              VectorStateType &kin_rates, MatrixStateType &dkin_rates_dn)
  {

//    if(z < 800. || z > 1200.)return;

 // neutrals resized
    VectorStateType full_concentrations;
    full_concentrations.resize(_composition.ionic_composition().n_species(),0.L);
    for(unsigned int s = 0; s < neutral_concentrations.size(); s++)
    {
       unsigned int i = _composition.ionic_composition().species_list()[_composition.neutral_composition().species_list()[s]];
       full_concentrations[i] = neutral_concentrations[s];
    }

//solve for ions
    VectorCoeffType source_ions;
    MatrixCoeffType drate_dn;
    source_ions.resize(_ionic_reactions.n_species(),0.L);
    drate_dn.resize(_ionic_reactions.n_species());
    for(unsigned int s = 0; s < _ionic_reactions.n_species(); s++)
    {
      drate_dn[s].resize(_ionic_reactions.n_species(),0.L);
    }
//all temperature conditions, solver deal with it
    _newton_solver.precompute_rates(full_concentrations,KC, KC.T(), _temperature.electronic_temperature(z)); 
    if(_newton_solver.steady_state_and_derivs(source_ions,drate_dn))
    {

//       std::cout << "Ionospheric activity at " << z << " km" << std::endl;
// update sources and derivs
      for(unsigned int s = 0; s < _composition.neutral_composition().n_species(); s++)
      {
        unsigned int i_neu = _composition.ionic_composition().species_list()[_composition.neutral_composition().species_list()[s]];
        kin_rates[s] += source_ions[i_neu];
        for(unsigned int q = 0; q < _composition.neutral_composition().n_species(); q++)
        {
           unsigned int j_neu = _composition.ionic_composition().species_list()[_composition.neutral_composition().species_list()[q]];
           dkin_rates_dn[s][q] += drate_dn[i_neu][j_neu];
        }
      }
    }
  }
}

#endif
