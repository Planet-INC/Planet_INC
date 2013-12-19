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

#ifndef PLANET_ATMOSPHERIC_MIXTURE_H
#define PLANET_ATMOSPHERIC_MIXTURE_H

//Antioch
#include "antioch/chemical_mixture.h"
#include "antioch/cmath_shims.h"

//Planet
#include "planet/atmospheric_temperature.h"
#include "planet/planet_constants.h"
#include "planet/math_constants.h"

//C++
#include <vector>

namespace Planet
{
   template<typename CoeffType, typename VectorCoeffType>
   class AtmosphericMixture
   {
      private:
        //!no default constructor
        AtmosphericMixture(){antioch_error();return;}

        Antioch::ChemicalMixture<CoeffType> &_neutral_composition;
        Antioch::ChemicalMixture<CoeffType> &_ionic_composition;

        CoeffType _total_bottom_density;
        VectorCoeffType _neutral_molar_fraction_bottom;

// neutral only
        VectorCoeffType _thermal_coefficient;
        VectorCoeffType _hard_sphere_radius;

////dependencies
        AtmosphericTemperature<CoeffType,VectorCoeffType> &_temperature;


        //! \return scale height of species s at altitude z, H = kb*T/(g*Ms)
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        H(const StateType &Ms, const StateType &temp, const StateType &alt) const
        ANTIOCH_AUTOFUNC(StateType, Constants::Universal::kb<StateType>() * Antioch::Constants::Avogadro<StateType>() * temp / 
                                    (StateType(1e-3L) * Ms // to kg
                                     * Constants::g(Constants::Titan::radius<StateType>(), alt, Constants::Titan::mass<StateType>()))
                                        
                        )

      public:
        AtmosphericMixture(Antioch::ChemicalMixture<CoeffType> &neutral, Antioch::ChemicalMixture<CoeffType> &ion,
                           AtmosphericTemperature<CoeffType,VectorCoeffType> &temp);
        ~AtmosphericMixture();

        //!
        template<typename StateType, typename VectorStateType>
        void init_composition(const VectorStateType &bot_compo,const StateType &dens_tot_bot);

        //!sets the hard sphere radius
        template <typename VectorStateType>
        void set_hard_sphere_radius(const VectorStateType &hsr);

        //!sets the thermal coefficient
        template <typename VectorStateType>
        void set_thermal_coefficient(const VectorStateType &tc);

/////

        //!\return thermal coefficient
        const VectorCoeffType &thermal_coefficient() const;

        //!\return hard sphere radius
        const VectorCoeffType &hard_sphere_radius() const;

        //!\return const reference to neutral composition
        const Antioch::ChemicalMixture<CoeffType> &neutral_composition() const;

        //!\return const reference to ionic composition
        const Antioch::ChemicalMixture<CoeffType> &ionic_composition() const;

        //! \return Jeans' escape flux
        //
        // \param ms: mass of molecule (kg)
        // \param ns: molecular density (any density unit, reported in the flux)
        // \param T: temperature (K)
        // \param z: altitude (km)
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        Jeans_flux(const StateType &ms, const StateType &ns, const StateType &T, const StateType &z) const
        ANTIOCH_AUTOFUNC(StateType, ns * Antioch::ant_sqrt(Constants::Universal::kb<StateType>() * T / (StateType(2.L) * ms * Constants::pi<StateType>())) 
                                       * Antioch::ant_exp(- ms * Constants::Universal::G<StateType>() * Constants::Titan::mass<StateType>() 
                                                         / (StateType(1e3L) * (Constants::Titan::radius<StateType>() + z) * Constants::Universal::kb<StateType>() * T)
                                                         )
                                       * (StateType(1.L) + 
                                           ((ms * Constants::Universal::G<StateType>() * Constants::Titan::mass<StateType>())
                                                / (StateType(1e3L) * (Constants::Titan::radius<StateType>() + z) * Constants::Universal::kb<StateType>() * T))
                                         ))


        //!
        template<typename StateType, typename VectorStateType>
        void scale_heights(const StateType &z, VectorStateType &Hs) const;

        //!
        template<typename StateType,typename VectorStateType>
        const CoeffType a(const VectorStateType &molar_densities,const StateType &z) const;

        //!
        template<typename StateType, typename VectorStateType>
        const CoeffType atmospheric_scale_height(const VectorStateType &molar_densities,const StateType &z) const;

        //!
        const CoeffType total_bottom_density() const;

   };


  template<typename CoeffType, typename VectorCoeffType>
  inline
  AtmosphericMixture<CoeffType,VectorCoeffType>::AtmosphericMixture(Antioch::ChemicalMixture<CoeffType> &neutral,
                                                                    Antioch::ChemicalMixture<CoeffType> &ion,
                                                                    AtmosphericTemperature<CoeffType,VectorCoeffType> &temp):
  _neutral_composition(neutral),
  _ionic_composition(ion),
  _temperature(temp)
  {
    return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  AtmosphericMixture<CoeffType,VectorCoeffType>::~AtmosphericMixture()
  {
     return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  const Antioch::ChemicalMixture<CoeffType> &AtmosphericMixture<CoeffType,VectorCoeffType>::neutral_composition() const
  {
     return _neutral_composition;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  const Antioch::ChemicalMixture<CoeffType> &AtmosphericMixture<CoeffType,VectorCoeffType>::ionic_composition() const
  {
     return _ionic_composition;
  }

  template<typename CoeffType, typename VectorCoeffType>
  template <typename VectorStateType>
  inline
  void AtmosphericMixture<CoeffType,VectorCoeffType>::set_hard_sphere_radius(const VectorStateType &hsr)
  {
      antioch_assert_equal_to(hsr.size(), _neutral_composition.n_species());
      _hard_sphere_radius = hsr;
     return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  template <typename VectorStateType>
  inline
  void AtmosphericMixture<CoeffType,VectorCoeffType>::set_thermal_coefficient(const VectorStateType &tc)
  {
      antioch_assert_equal_to(tc.size(), _neutral_composition.n_species());
      _thermal_coefficient = tc;
     return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  const VectorCoeffType &AtmosphericMixture<CoeffType,VectorCoeffType>::thermal_coefficient() const
  {
     return _thermal_coefficient;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  const VectorCoeffType &AtmosphericMixture<CoeffType,VectorCoeffType>::hard_sphere_radius() const
  {
     return _hard_sphere_radius;
  }

  template<typename CoeffType, typename VectorCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void AtmosphericMixture<CoeffType,VectorCoeffType>::init_composition(const VectorStateType &bot_compo,
                                                                       const StateType &dens_tot_bot)
  {
    _total_bottom_density = dens_tot_bot;
    _neutral_molar_fraction_bottom = bot_compo;

  }

  template<typename CoeffType, typename VectorCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void AtmosphericMixture<CoeffType,VectorCoeffType>::scale_heights(const StateType &z, VectorStateType &Hs) const
  {
      Hs.resize(_neutral_composition.n_species(),0.L);
      CoeffType T = _temperature.neutral_temperature(z);
      for(unsigned int s = 0; s < _neutral_composition.n_species(); s++)
      {
         Hs[s] = this->H(_neutral_composition.M(s),T,z);
      }
      return;
  }

  
  template<typename CoeffType, typename VectorCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  const CoeffType AtmosphericMixture<CoeffType,VectorCoeffType>::atmospheric_scale_height(const VectorStateType &molar_densities, 
                                                                                          const StateType &z) const
  {
    antioch_assert_equal_to(molar_densities.size(),_neutral_composition.n_species());

    CoeffType Mm;
    Antioch::set_zero(Mm);
    CoeffType nTot;
    Antioch::set_zero(nTot);
    for(unsigned int s = 0; s < _neutral_composition.n_species(); s++)
    {
      Mm   += molar_densities[s] * _neutral_composition.M(s);
      nTot += molar_densities[s];
    }
    Mm /= nTot;

    return (this->H(Mm,_temperature.neutral_temperature(z),z));
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  const CoeffType AtmosphericMixture<CoeffType,VectorCoeffType>::total_bottom_density() const
  {
     return _total_bottom_density;
  }

  template<typename CoeffType, typename VectorCoeffType>
  template<typename StateType,typename VectorStateType>
  inline
  const CoeffType AtmosphericMixture<CoeffType,VectorCoeffType>::a(const VectorStateType &molar_densities,const StateType &z) const
  {
     return (Constants::Titan::radius<CoeffType>() + z) / this->atmospheric_scale_height(molar_densities,z) * CoeffType(1e3); // to m
  }

}

#endif
