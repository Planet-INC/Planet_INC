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
   template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
   class AtmosphericMixture
   {
      private:
        //!no default constructor
        AtmosphericMixture(){antioch_error();return;}
        Antioch::Species iH;  // unsigned int
        Antioch::Species iH2; // unsigned int
//boundaries
        CoeffType _zmin;
        CoeffType _zmax;

        Antioch::ChemicalMixture<CoeffType> &_neutral_composition;
        Antioch::ChemicalMixture<CoeffType> &_ionic_composition;

        CoeffType _total_bottom_density;
        VectorCoeffType _neutral_molar_fraction_bottom;

// neutral only
        VectorCoeffType _thermal_coefficient;
        VectorCoeffType _hard_sphere_radius;

// neutral only, precomputations
        MatrixCoeffType _mean_free_path_precompute;

////dependencies
        AtmosphericTemperature<CoeffType,VectorCoeffType> &_temperature;


        void precompute_mean_free_path();

        /*! \return scale height of species s at altitude z, H = R * T / (g * Ms)
         *
         * \param  Ms: molar mass in kg/mol
         * \param  T: temperature in K
         * \param  alt: altitude in km
         * \return scale height in km
         */
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        H(const CoeffType &Ms, const StateType &T, const StateType &alt) const
        ANTIOCH_AUTOFUNC(StateType, Antioch::constant_clone(T,1e-3) * // m -> km
                                    Antioch::Constants::R_universal<StateType>() * T / 
                                    ( Ms * Constants::g(Constants::Titan::radius<StateType>(), alt, Constants::Titan::mass<StateType>()))
                        )

        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        barometry_density(const StateType &z, const StateType &zmin, const StateType &zmin_dens, const StateType &T, const StateType &Mm) const
        ANTIOCH_AUTOFUNC(StateType, zmin_dens * Antioch::ant_exp( - (z - zmin) /
                                                                 ( (Planet::Constants::Titan::radius<StateType>() + z) * 
                                                                   (Planet::Constants::Titan::radius<StateType>() + zmin) * Antioch::constant_clone(z,1e3) * //to SI (km -> m)
                                                                    Antioch::Constants::Avogadro<StateType>() * Planet::Constants::Universal::kb<StateType>() * T / 
                                                                   (Planet::Constants::Universal::G<StateType>() * Planet::Constants::Titan::mass<StateType>() * Mm)
                                                                 )
                                                               )
                        )


      public:
        AtmosphericMixture(Antioch::ChemicalMixture<CoeffType> &neutral, Antioch::ChemicalMixture<CoeffType> &ion,
                           AtmosphericTemperature<CoeffType,VectorCoeffType> &temp);
        ~AtmosphericMixture();

        //!
        template<typename StateType, typename VectorStateType>
        void init_composition(const VectorStateType &bot_compo,const StateType &dens_tot_bot,
                              const StateType &zmin, const StateType &zmax);

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

        //!\return the mean free path in km
        //
        // \param densities: cm-3
        // \return \f$l_s = 10^{-5} / \sum_i n_i  \sigma_{is} \sqrt{1 + \frac{m_s}{m_i}}\f$
        template <typename VectorStateType>
        void mean_free_path(const VectorStateType &densities, VectorStateType &mean_free_path) const;

        //! \return Jeans' escape flux (cm-3.km.s-1)
        //
        // \param ms: mass of molecule (kg/mol)
        // \param ns: molecular density (density)
        // \param T: temperature (K)
        // \param z: altitude (km)
        // \return Jeans' escape flux in density.km/s
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        Jeans_flux(const StateType &ms, const StateType &ns, const StateType &T, const StateType &z) const
        ANTIOCH_AUTOFUNC(StateType, Antioch::constant_clone(T,1e-3) * // m -> km
                                    ns * Antioch::ant_sqrt(Antioch::Constants::R_universal<StateType>() * T / (Antioch::constant_clone(T,2.) * ms * Constants::pi<StateType>())) 
                                       * Antioch::ant_exp(- ms * Constants::Universal::G<StateType>() * Constants::Titan::mass<StateType>() 
                                                         / (Antioch::constant_clone(T,1e3) * (Constants::Titan::radius<StateType>() + z) * Antioch::Constants::R_universal<StateType>() * T)
                                                         )
                                       * (Antioch::constant_clone(T,1.) + 
                                           ((ms * Constants::Universal::G<StateType>() * Constants::Titan::mass<StateType>())
                                                / (Antioch::constant_clone(T,1e3) * (Constants::Titan::radius<StateType>() + z) * Antioch::Constants::R_universal<StateType>() * T))
                                         ))

        //! \return the derivative of Jeans' escape flux (km.s-1)
        //
        // \param ms: mass of molecule (kg/mol)
        // \param T: temperature (K)
        // \param z: altitude (km)
        // \return Jeans' escape speed in km/s
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        Jeans_velocity(const CoeffType &ms, const StateType &T, const CoeffType &z) const
        ANTIOCH_AUTOFUNC(StateType, //Antioch::constant_clone(T,1e-3) * //ns m -> km
                                         Antioch::ant_sqrt(Antioch::Constants::R_universal<StateType>() * T / (Antioch::constant_clone(T,2.) * ms * Constants::pi<StateType>())) 
                                       * Antioch::ant_exp(- ms * Constants::Universal::G<StateType>() * Constants::Titan::mass<StateType>() 
                                                         / (Antioch::constant_clone(T,1e3) * (Constants::Titan::radius<StateType>() + z) * Antioch::Constants::R_universal<StateType>() * T)
                                                         )
                                       * (Antioch::constant_clone(T,1.) + 
                                           ((ms * Constants::Universal::G<StateType>() * Constants::Titan::mass<StateType>())
                                                / (Antioch::constant_clone(T,1e3) * (Constants::Titan::radius<StateType>() + z) * Antioch::Constants::R_universal<StateType>() * T))
                                         ))

        //!use isobaric equation and molar fractions at bottom
        template <typename StateType, typename VectorStateType>
        void first_guess_densities(const StateType &z, VectorStateType &densities) const;

        template <typename StateType>
        StateType first_guess_density(const StateType &z, unsigned int species) const;

        //!use isobaric equation, molar fractions at bottom and zstep = 10
        template <typename StateType, typename VectorStateType>
        void first_guess_densities_sum(const StateType &z, VectorStateType &sum_densities) const;

        //!lower boundary concentrations
        template <typename VectorStateType>
        void lower_boundary_concentrations(VectorStateType &low_densities) const;

        //!lower boundary concentrations
        template <typename VectorStateType>
        void upper_boundary_fluxes(VectorStateType &upper_fluxes, const VectorStateType &molar_concentrations) const;

        //!upper boundary condition
        template <typename VectorStateType>
        typename Antioch::value_type<VectorStateType>::type
          upper_boundary_flux(const VectorStateType &molar_concentrations, unsigned int s) const;

        //!upper boundary derived condition
        template <typename StateType>
        StateType upper_boundary_velocity(unsigned int s, const StateType & /*ex*/) const;

        //!
        template<typename StateType, typename VectorStateType>
        void scale_heights(const StateType &z, VectorStateType &Hs) const;

        //! \return a factor (see model documentation Eq. 2.7, no dimension
        //
        //  a = (Rtitan + z) / H
        template<typename StateType,typename VectorStateType>
        ANTIOCH_AUTO(StateType)
        a(const VectorStateType &molar_densities,const StateType &z) const
        ANTIOCH_AUTOFUNC(StateType, (Constants::Titan::radius<StateType>() + z) / 
                                       this->atmospheric_scale_height(molar_densities,z))

        //!
        template<typename StateType, typename VectorStateType>
        const CoeffType atmospheric_scale_height(const VectorStateType &molar_densities,const StateType &z) const;

        //!
        template<typename StateType, typename VectorStateType>
        void datmospheric_scale_height_dn_i(const VectorStateType &molar_densities, const StateType &z, StateType & Ha, VectorStateType & dHa_dn_i) const;

        //!
        const CoeffType total_bottom_density() const;

        //!
        const VectorCoeffType neutral_molar_fraction_bottom() const;

   };


  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::AtmosphericMixture(Antioch::ChemicalMixture<CoeffType> &neutral,
                                                                    Antioch::ChemicalMixture<CoeffType> &ion,
                                                                    AtmosphericTemperature<CoeffType,VectorCoeffType> &temp):
  iH(neutral.n_species()+1),
  iH2(neutral.n_species()+1),
  _zmin(0.L),
  _zmax(0.L),
  _neutral_composition(neutral),
  _ionic_composition(ion),
  _temperature(temp)
  {
    if(neutral.species_name_map().count("H"))iH = neutral.species_name_map().at("H");
    if(neutral.species_name_map().count("H2"))iH2 = neutral.species_name_map().at("H2");

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::~AtmosphericMixture()
  {
     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const Antioch::ChemicalMixture<CoeffType> &AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::neutral_composition() const
  {
     return _neutral_composition;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const Antioch::ChemicalMixture<CoeffType> &AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::ionic_composition() const
  {
     return _ionic_composition;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template <typename VectorStateType>
  inline
  void AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::set_hard_sphere_radius(const VectorStateType &hsr)
  {
      antioch_assert_equal_to(hsr.size(), _neutral_composition.n_species());
      _hard_sphere_radius = hsr;
      this->precompute_mean_free_path();
     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template <typename VectorStateType>
  inline
  void AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::set_thermal_coefficient(const VectorStateType &tc)
  {
      antioch_assert_equal_to(tc.size(), _neutral_composition.n_species());
      _thermal_coefficient = tc;
     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const VectorCoeffType &AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::thermal_coefficient() const
  {
     return _thermal_coefficient;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const VectorCoeffType &AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::hard_sphere_radius() const
  {
     return _hard_sphere_radius;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::init_composition(const VectorStateType &bot_compo,
                                                                       const StateType &dens_tot_bot,
                                                                       const StateType &zmin, const StateType &zmax)
  {
    _total_bottom_density = dens_tot_bot;
    _neutral_molar_fraction_bottom.resize(_neutral_composition.n_species(),0.L);
    for(unsigned int s = 0; s < _neutral_composition.n_species(); s++)
    {
      _neutral_molar_fraction_bottom[s] = bot_compo[s];
    }
    _zmin = zmin;
    _zmax = zmax;   

  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::scale_heights(const StateType &z, VectorStateType &Hs) const
  {
      Hs.resize(_neutral_composition.n_species(),0.L);
      CoeffType T = _temperature.neutral_temperature(z);
      for(unsigned int s = 0; s < _neutral_composition.n_species(); s++)
      {
         Hs[s] = this->H(_neutral_composition.M(s),T,z);
      }
      return;
  }

  
  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  const CoeffType AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::atmospheric_scale_height(const VectorStateType &molar_densities, 
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

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::datmospheric_scale_height_dn_i(const VectorStateType &molar_densities, 
                                                                                                    const StateType &z,
                                                                                                    StateType & Ha, VectorStateType & dHa_dn_i) const
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

    Ha = this->H(Mm / nTot,_temperature.neutral_temperature(z),z); 

    for(unsigned int s = 0; s < _neutral_composition.n_species(); s++)
    {
      dHa_dn_i[s] = Ha / nTot;
    }
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const CoeffType AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::total_bottom_density() const
  {
     return _total_bottom_density;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const VectorCoeffType AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::neutral_molar_fraction_bottom() const
  {
     return _neutral_molar_fraction_bottom;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename VectorStateType>
  inline
  void AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::mean_free_path(const VectorStateType &densities, VectorStateType &mean_free_path) const
  {
     antioch_assert_equal_to(densities.size(),_neutral_composition.n_species());
     antioch_assert(!_mean_free_path_precompute.empty());

     mean_free_path.resize(densities.size(),0.L);
     CoeffType out(0.L);
     for(unsigned int s = 0; s < densities.size(); s++)
     {
       for(unsigned int n = 0; n < densities.size(); n++)
       {
          out += densities[n] * _mean_free_path_precompute[s][n];
       }
       mean_free_path[s] = Antioch::constant_clone(densities[0],1e5) / out;
     }
  }

  
  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  void AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::precompute_mean_free_path()
  {
     antioch_assert(!_hard_sphere_radius.empty());
     _mean_free_path_precompute.resize(_hard_sphere_radius.size());

     for(unsigned int s = 0; s < _hard_sphere_radius.size(); s++)
     {
        _mean_free_path_precompute[s].resize(_hard_sphere_radius.size());
        for(unsigned int n = 0; n < _hard_sphere_radius.size(); n++)
        {
            _mean_free_path_precompute[s][n] = Constants::pi<CoeffType>() * (_hard_sphere_radius[s] + _hard_sphere_radius[n]) 
                                                                          * (_hard_sphere_radius[s] + _hard_sphere_radius[n])
                                               * Antioch::ant_sqrt(CoeffType(1.L) + _neutral_composition.M(s)/_neutral_composition.M(n));
        }
     }
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  void AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::first_guess_densities(const StateType &z, VectorStateType &densities) const
  {
      antioch_assert_equal_to(densities.size(),_neutral_composition.n_species());
      antioch_assert(!_neutral_molar_fraction_bottom.empty());

      CoeffType Mm(0.L);
      for(unsigned int s = 0; s < _neutral_composition.n_species(); s++)
      {
         Mm += _neutral_molar_fraction_bottom[s] * _neutral_composition.M(s);
      }

      CoeffType nTot = this->barometry_density(z, _zmin, _total_bottom_density, _temperature.neutral_temperature(z), Mm);

      for(unsigned int s = 0; s < _neutral_composition.n_species(); s++)
      {
         densities[s] = _neutral_molar_fraction_bottom[s] * nTot;
      }
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template <typename StateType>
  inline
  StateType AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::first_guess_density(const StateType& z, unsigned int species) const
  {
      antioch_assert_less(species,_neutral_composition.n_species());
      antioch_assert(!_neutral_molar_fraction_bottom.empty());

      CoeffType Mm(0.L);
      for(unsigned int s = 0; s < _neutral_composition.n_species(); s++)
      {
         Mm += _neutral_molar_fraction_bottom[s] * _neutral_composition.M(s);
      }

      return _neutral_molar_fraction_bottom[species] * this->barometry_density(z, _zmin, _total_bottom_density, _temperature.neutral_temperature(z), Mm);
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  void AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::first_guess_densities_sum(const StateType &z, VectorStateType &sum_densities) const
  {
      antioch_assert_equal_to(sum_densities.size(),_neutral_composition.n_species());
      antioch_assert(!_neutral_molar_fraction_bottom.empty());

      CoeffType z_tmp_step(10.L);
      Antioch::set_zero(sum_densities);
      VectorCoeffType densities_z = Antioch::zero_clone(sum_densities);

      for(CoeffType z_tmp = z; z_tmp <= _zmax - z_tmp_step; z_tmp += z_tmp_step) //integration is ns_{i} = n_{i} * (z_{i+1} - z_{i}), i from bottom to top
      {
         this->first_guess_densities(z_tmp,densities_z);
         for(unsigned int s = 0; s < sum_densities.size(); s++)
         {
            sum_densities[s] += densities_z[s]* z_tmp_step;;
         }
      }
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template <typename VectorStateType>
  inline
  void AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::lower_boundary_concentrations(VectorStateType &low_densities) const
  {
      antioch_assert_equal_to(low_densities.size(),_neutral_composition.n_species());

      for(unsigned int s = 0; s < _neutral_composition.n_species(); s++)
      {
          low_densities[s] = _neutral_molar_fraction_bottom[s] * _total_bottom_density;
      }

      return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template <typename VectorStateType>
  inline
  void AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::upper_boundary_fluxes(VectorStateType &upper_fluxes, const VectorStateType &molar_concentrations) const
  {
      antioch_assert_equal_to(upper_fluxes.size(),_neutral_composition.n_species());
      antioch_assert_equal_to(upper_fluxes.size(),molar_concentrations.size());
      std::fill(upper_fluxes.begin(),upper_fluxes.end(),0.);
      for(unsigned int s = 0; s < _neutral_composition.n_species(); s++)
      {
          if(_neutral_composition.species_list()[s] != iH &&
             _neutral_composition.species_list()[s] != iH2)continue;
          upper_fluxes[s] = - this->Jeans_flux(_neutral_composition.M(s),molar_concentrations[s],_temperature.neutral_temperature(_zmax),_zmax); // cm-3.km/s, escaping flux, term < 0
      }
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template <typename VectorStateType>
  inline
  typename Antioch::value_type<VectorStateType>::type
    AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::upper_boundary_flux(const VectorStateType &molar_concentrations, unsigned int s) const
  {
      antioch_assert_less(s,_neutral_composition.n_species());
      antioch_assert_less(s,molar_concentrations.size());

      return (_neutral_composition.species_list()[s] != iH && _neutral_composition.species_list()[s] != iH2)?
                        0.:
                        - this->Jeans_flux(_neutral_composition.M(s), molar_concentrations[s],_temperature.neutral_temperature(_zmax),_zmax); // cm-3.km.s-1, escaping flux, term < 0;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template <typename StateType>
  inline
  StateType AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::upper_boundary_velocity(unsigned int s, const StateType & /*ex*/) const
  {
      antioch_assert_less(s,_neutral_composition.n_species());

      return (_neutral_composition.species_list()[s] != iH && _neutral_composition.species_list()[s] != iH2)?
                        0.:
                        - this->Jeans_velocity(_neutral_composition.M(s), _temperature.neutral_temperature(_zmax),_zmax); // cm-3.km.s-1, escaping flux, term < 0
  }


}

#endif
