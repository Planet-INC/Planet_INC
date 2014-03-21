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
        Antioch::ReactionSet<CoeffType>       &_neutral_reactions_set;
        Antioch::KineticsEvaluator<CoeffType> &_ionic_reactions;

        bool _ionic_coupling;
        std::vector<Antioch::Species> _ions_species;
        VectorCoeffType _cache_concentrations;
        
//
        const AtmosphericTemperature<CoeffType,VectorCoeffType>             &_temperature;
        PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>    &_photon;
        const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &_composition;
      public:
        //!
        AtmosphericKinetics(Antioch::KineticsEvaluator<CoeffType>                         &neu,
                            Antioch::ReactionSet<CoeffType>                               &neu_set,
                            Antioch::KineticsEvaluator<CoeffType>                         &ion,
                            const AtmosphericTemperature<CoeffType,VectorCoeffType>             &temperature,
                            PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>    &photon,
                            const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &composition);
        //!
        ~AtmosphericKinetics();

        //!\return neutral kinetics system, writable reference
        Antioch::KineticsEvaluator<CoeffType> &neutral_kinetics();

        //!\return ionic kinetics system, writable reference
        Antioch::KineticsEvaluator<CoeffType> &ionic_kinetics();

        //! compute chemical net rate and provide them in kin_rates
        template<typename StateType, typename VectorStateType>
        void chemical_rate(const VectorStateType &molar_concentrations, const VectorStateType &sum_concentrations, 
                           const StateType &z, VectorStateType &kin_rates);

        template<typename StateType, typename VectorStateType, typename MatrixStateType>
        void chemical_rate_and_derivs(const VectorStateType &molar_concentrations, const VectorStateType &sum_concentrations, 
                                      const StateType &z, VectorStateType &kin_rates, MatrixStateType &dkin_rates_dn) const;

        //! Newton solver for the ionic system
        template<typename StateType, typename VectorStateType>
        void add_ionic_contribution(const VectorStateType &molar_concentrations, const StateType &z, VectorStateType &kin_rates);
  };


  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  AtmosphericKinetics<CoeffType,VectorCoeffType,MatrixCoeffType>::AtmosphericKinetics(Antioch::KineticsEvaluator<CoeffType>         &neu,
                                                                      Antioch::ReactionSet<CoeffType>                               &neu_set,
                                                                      Antioch::KineticsEvaluator<CoeffType>                         &ion,
                                                                      const AtmosphericTemperature<CoeffType,VectorCoeffType>             &temperature,
                                                                      PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>    &photon,
                                                                      const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &composition):
   _neutral_reactions(neu),
   _neutral_reactions_set(neu_set),
   _ionic_reactions(ion),
   _temperature(temperature),
   _photon(photon),
   _composition(composition)
  {
    _ionic_coupling = (_ionic_reactions.n_reactions() != 0);
    if(_ionic_coupling)
    {
       for(unsigned int s = 0; s < _composition.ionic_composition().n_species(); s++)
       {
           if(!_composition.neutral_composition().species_list_map().count(_composition.ionic_composition().species_list()[s])) //if not in the neutral system
                        _ions_species.push_back(_composition.ionic_composition().species_list()[s]); // then adds here
       }
       if(_ions_species.empty())_ionic_coupling = false;
      _cache_concentrations.resize(_composition.ionic_composition().n_species(),-1.L);
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
                                                                     const VectorStateType &sum_concentrations, 
                                                                     const StateType &z,
                                                                     VectorStateType &kin_rates)
  {
     antioch_assert_equal_to(kin_rates.size(),_composition.neutral_composition().n_species());
     VectorStateType dummy;
     dummy.resize(_composition.neutral_composition().n_species(),0.L); //everything is irreversible
     _photon.update_photon_flux(molar_concentrations, sum_concentrations, z);
     _neutral_reactions_set.update_particle_flux_chemistry();
     _neutral_reactions.compute_mole_sources(_temperature.neutral_temperature(z),
                                             molar_concentrations,dummy,kin_rates);

     this->add_ionic_contribution(molar_concentrations,z,kin_rates);

     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType, typename MatrixStateType>
  inline
  void AtmosphericKinetics<CoeffType,VectorCoeffType,MatrixCoeffType>::chemical_rate_and_derivs(const VectorStateType &molar_concentrations, const VectorStateType &sum_concentrations, 
                                      const StateType &z, VectorStateType &kin_rates, MatrixStateType &dkin_rates_dn) const
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
     _photon.update_photon_flux(molar_concentrations, sum_concentrations, z);
     _neutral_reactions.compute_mole_sources_and_derivs(_temperature.neutral_temperature(z),molar_concentrations,dummy,ddummy_dT,
                                                        kin_rates,dkin_dT,dkin_rates_dn);

     this->add_ionic_contribution(molar_concentrations,z,kin_rates);
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void AtmosphericKinetics<CoeffType,VectorCoeffType,MatrixCoeffType>::add_ionic_contribution(const VectorStateType &neutral_concentrations, const StateType &z, 
                                                                                              VectorStateType &kin_rates)
  {
    if(!_ionic_coupling)return;
    antioch_assert(!_cache_concentrations.empty());


 // update neutrals
    for(unsigned int s = 0; s < neutral_concentrations.size(); s++)
    {
       unsigned int i = _composition.ionic_composition().species_list_map().at(_composition.neutral_composition().species_list()[s]);
       _cache_concentrations[i] = neutral_concentrations[s];
    }

//solve for ions
    AtmosphericSteadyState solver; //ionospheric solver
    VectorCoeffType source_ions;
    source_ions.resize(_ionic_reactions.n_species(),0.L);
    solver.steady_state(_ionic_reactions,_ions_species,_composition.ionic_composition(),_temperature.neutral_temperature(z),_cache_concentrations,source_ions);

// update sources
    for(unsigned int s = 0; s < _composition.neutral_composition().n_species(); s++)
    {
        unsigned int i_neu = _composition.ionic_composition().species_list_map().at(_composition.neutral_composition().species_list()[s]);
        kin_rates[s] += source_ions[i_neu];
    }
  }
}

#endif
