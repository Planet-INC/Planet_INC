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
        std::vector<Antioch::Species> _ions_species;
        
//
        const AtmosphericTemperature<CoeffType,VectorCoeffType>             &_temperature;
        PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>    &_photon;
        const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &_composition;
      public:
        //!
        AtmosphericKinetics(Antioch::KineticsEvaluator<CoeffType>                         &neu,
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
                           const StateType &z, VectorStateType &kin_rates) const;

        template<typename StateType, typename VectorStateType, typename MatrixStateType>
        void chemical_rate_and_derivs(const VectorStateType &molar_concentrations, const VectorStateType &sum_concentrations, 
                                      const StateType &z, VectorStateType &kin_rates, MatrixStateType &dkin_rates_dn) const;

        //! Newton solver for the ionic system
        template<typename StateType, typename VectorStateType>
        void add_ionic_contribution(const VectorStateType &molar_concentrations, const StateType &z, VectorStateType &kin_rates) const;
  };


  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  AtmosphericKinetics<CoeffType,VectorCoeffType,MatrixCoeffType>::AtmosphericKinetics(Antioch::KineticsEvaluator<CoeffType>         &neu,
                                                                      Antioch::KineticsEvaluator<CoeffType>                         &ion,
                                                                      const AtmosphericTemperature<CoeffType,VectorCoeffType>             &temperature,
                                                                      PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>    &photon,
                                                                      const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &composition):
   _neutral_reactions(neu),
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
                                                                     VectorStateType &kin_rates) const
  {
     antioch_assert_equal_to(kin_rates.size(),_composition.neutral_composition().n_species());
     VectorStateType dummy;
     dummy.resize(_composition.neutral_composition().n_species(),0.L); //everything is irreversible
     _photon.update_photon_flux(molar_concentrations, sum_concentrations, z);
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
                                                                                              VectorStateType &kin_rates) const
  {
    if(!_ionic_coupling)return;
    
// Newton solver here
// Ax + b = 0
// A is jacobian, b is what goes to 0 (dc/dt here)
    Eigen::Matrix<CoeffType,Eigen::Dynamic,Eigen::Dynamic> A;
    Eigen::Matrix<CoeffType,Eigen::Dynamic,1> b;


    VectorCoeffType molar_concentrations;
    molar_concentrations.resize(_ionic_reactions.n_species(),0.L); //full system
    for(unsigned int s = 0; s < neutral_concentrations.size(); s++)// full system
    {
       unsigned int i = _composition.ionic_composition().species_list_map().at(_composition.neutral_composition().species_list()[s]);
       molar_concentrations[i] = neutral_concentrations[s]; // add it
    }


    VectorCoeffType h_RT_minus_s_R;
    VectorCoeffType dh_RT_minus_s_R_dT;
    VectorCoeffType mole_sources;
    VectorCoeffType dmole_dT;
    MatrixCoeffType dmole_dX_s;
    h_RT_minus_s_R.resize(_ionic_reactions.n_species(),0.L); //irreversible
    dh_RT_minus_s_R_dT.resize(_ionic_reactions.n_species());
    mole_sources.resize(_ionic_reactions.n_species());
    dmole_dT.resize(_ionic_reactions.n_species());
    dmole_dX_s.resize(_ionic_reactions.n_species());

    CoeffType lim(1.L);
    CoeffType thresh = std::numeric_limits<CoeffType>::epsilon();
    if(thresh < 1e-10)thresh = 1e-10; // physically this precision is ridiculous, which is nice
    unsigned int loop_max(50);
    unsigned int nloop(0);
    while(lim > thresh)
    {

      _ionic_reactions.compute_mole_sources_and_derivs(_temperature.neutral_temperature(z), molar_concentrations,
                                                       h_RT_minus_s_R, dh_RT_minus_s_R_dT,
                                                       mole_sources, dmole_dT, dmole_dX_s );


      for(unsigned int i = 0; i < _ions_species.size(); i++)
      {
        unsigned int i_ion = _composition.ionic_composition().species_list_map().at(_ions_species[i]);
        for(unsigned int j = 0; j < _ions_species.size(); j++)
        {
           unsigned int j_ion = _composition.ionic_composition().species_list_map().at(_ions_species[j]);
           A(i,j) = dmole_dX_s[i_ion][j_ion];
        }
        b(i) = - mole_sources[i_ion];
      }
    
      Eigen::PartialPivLU<Eigen::Matrix<CoeffType,Eigen::Dynamic,Eigen::Dynamic> > mypartialPivLu(A);
      Eigen::Matrix<CoeffType,Eigen::Dynamic,1> x(_ions_species.size());
      x = mypartialPivLu.solve(b);

      Antioch::set_zero(lim);
      for(unsigned int s = 0; s < _ions_species.size(); s++)
      {
        unsigned int i_ion = _composition.ionic_composition().species_list_map().at(_ions_species[s]);
        molar_concentrations[i_ion] += x(s);
        lim += (x(s) < 0.)?-x(s):x(s);
      }

      nloop++;
      if(nloop > loop_max)antioch_error();

    }

    for(unsigned int s = 0; s < _composition.neutral_composition().n_species(); s++)
    {
        unsigned int i_neu = _composition.ionic_composition().species_list_map().at(_composition.neutral_composition().species_list()[s]);
        kin_rates[s] += mole_sources[i_neu];
    }
      
  }
}

#endif
