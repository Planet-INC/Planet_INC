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
#include "planet/altitude.h"
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

        Antioch::ChemicalMixture<CoeffType> &_neutral_composition;
        Antioch::ChemicalMixture<CoeffType> &_ionic_composition;

        VectorCoeffType _total_density;
        MatrixCoeffType _neutral_molar_fraction; //nneus,alt
        MatrixCoeffType _ionic_molar_fraction;
// neutral only
        VectorCoeffType _thermal_coefficient;
        VectorCoeffType _hard_sphere_radius;
        MatrixCoeffType _scale_height;
        MatrixCoeffType _free_path;
        std::vector<unsigned int> _exobase; //index of exobase altitudes
// mean atmosphere
        VectorCoeffType _mean_free_path;
        VectorCoeffType _mean_scale_height;
        VectorCoeffType _a_factor;
        unsigned int    _mean_exobase; 

////dependencies
        Altitude<CoeffType,VectorCoeffType> &_altitude;
        AtmosphericTemperature<CoeffType,VectorCoeffType> &_temperature;


        //! \return scale height of species s at altitude z, H = kb*T/(g*Ms)
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        H(const StateType &Ms, const StateType &temp, const StateType &alt) const
        ANTIOCH_AUTOFUNC(StateType, Constants::Universal::kb<StateType>() * temp / 
                                    (StateType(1e-3L) * Constants::g(Constants::Titan::radius<StateType>(), alt, Constants::Titan::mass<StateType>()) *
                                        Ms / Antioch::Constants::Avogadro<StateType>())
                        )

        //! initializes exobases, defaults values here
        void initialize_exobases();

      public:
        AtmosphericMixture(Antioch::ChemicalMixture<CoeffType> &neutral, Antioch::ChemicalMixture<CoeffType> &ion,
                           Altitude<CoeffType,VectorCoeffType> &alt, AtmosphericTemperature<CoeffType,VectorCoeffType> &temp);
        ~AtmosphericMixture();

        //!
        template<typename StateType, typename VectorStateType>
        void init_composition(const VectorStateType &bot_compo,const StateType &dens_tot_bot);

        //\return total density
        const VectorCoeffType &total_density() const;

        //!\return neutral molar fraction
        const MatrixCoeffType &neutral_molar_fraction() const;

        //!\return ionic molar fraction
        const MatrixCoeffType &ionic_molar_fraction() const;

        //!\return thermal coefficient
        const VectorCoeffType &thermal_coefficient() const;

        //!\return scales height
        const MatrixCoeffType &scale_height() const;

        //!\return hard sphere radius
        const VectorCoeffType &hard_sphere_radius() const;

        //!\return scales height
        const MatrixCoeffType &free_path() const;

        //!\return exobase altitudes
        const std::vector<unsigned int> &exobase() const;


        //!\return mean free path
        const VectorCoeffType &atmosphere_free_path() const;

        //!\return mean scale height
        const VectorCoeffType &atmosphere_scale_height() const;

        //!\return a factor (R_Titan + z)/H(z)
        const VectorCoeffType &a_factor() const;

        //!\return mean exobase altitude
        unsigned int atmosphere_exobase() const;


        //! \return mean neutral molar mass at altitude z
        //
        //Only neutral species
        CoeffType mean_neutral_molar_mass(unsigned int iz) const;

        //!\return const reference to neutral composition
        const Antioch::ChemicalMixture<CoeffType> &neutral_composition() const;

        //!\return const reference to ionic composition
        const Antioch::ChemicalMixture<CoeffType> &ionic_composition() const;

        //!\return the total number of species
        unsigned int n_total_species() const;


        //!sets the hard sphere radius
        template <typename VectorStateType>
        void set_hard_sphere_radius(const VectorStateType &hsr);

        //!sets the thermal coefficient
        template <typename VectorStateType>
        void set_thermal_coefficient(const VectorStateType &tc);

        //! initializes/updates free pathes
        void make_free_pathes();

        //! initializes/updates scale height and a factor
        void make_scale_height();

        //! initializes/updates exobases
        //
        // Beware about the INT hypothesises, not accounted for:
        //   - all species exobase = N2 exobase except CH4, H and H2
        //   - if no coupling with neutral, exobase = 1590 km by default
        void make_exobase();

        //!
        void initialize();


   };


  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::AtmosphericMixture(Antioch::ChemicalMixture<CoeffType> &neutral,
                                                                                    Antioch::ChemicalMixture<CoeffType> &ion,
                                                                                    Altitude<CoeffType,VectorCoeffType> &alt, 
                                                                                    AtmosphericTemperature<CoeffType,VectorCoeffType> &temp):
  _neutral_composition(neutral),
  _ionic_composition(ion),
  _altitude(alt),
  _temperature(temp)
  {
    this->initialize_exobases();
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
  void AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::initialize()
  {
     this->make_free_pathes();
     this->make_scale_height();
     this->make_exobase();
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
  inline
  void AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::make_free_pathes()
  {
    antioch_assert_equal_to(_hard_sphere_radius.size(),_neutral_composition.n_species());

    _free_path.clear();
    _free_path.resize(_neutral_composition.n_species());
    _mean_free_path.resize(_altitude.altitudes().size(),0.L);

    for(unsigned int ineu = 0; ineu < _neutral_composition.n_species(); ineu++)
    {
       _free_path[ineu].resize(_altitude.altitudes().size(),1.L);
       for(unsigned int iz = 0; iz < _altitude.altitudes().size(); iz++)
       {
         CoeffType sum(0.L);
         for(unsigned int jneu = 0; jneu <_neutral_composition.n_species(); jneu++)
         {
           sum += _neutral_molar_fraction[jneu][iz] * _total_density[iz] * 1e6 * //cm-3 -> m-3
                  Constants::pi<CoeffType>() *
                  Antioch::ant_pow(_hard_sphere_radius[jneu] + _hard_sphere_radius[ineu],2) * 
                  Antioch::ant_sqrt(CoeffType(1.L) + _neutral_composition.M(ineu)/_neutral_composition.M(jneu));
         }
         _free_path[ineu][iz] /= sum;
         _mean_free_path[iz] += _free_path[ineu][iz] * _neutral_molar_fraction[ineu][iz];
      }
    }

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  void AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::make_exobase()
  {
//they are necessary for this determination
     antioch_assert_equal_to(_scale_height.size(),_neutral_composition.n_species());
     antioch_assert_equal_to(_free_path.size(),_neutral_composition.n_species());
     for(unsigned int s = 0; s < _neutral_composition.n_species(); s++)
     {
        antioch_assert_equal_to(_scale_height[s].size(),_altitude.altitudes().size());
        antioch_assert_equal_to(_free_path[s].size(),_altitude.altitudes().size());
     }
     antioch_assert_equal_to(_mean_scale_height.size(),_altitude.altitudes().size());
     antioch_assert_equal_to(_mean_free_path.size(),_altitude.altitudes().size());

//initialized in the constructor, should not be resized
     antioch_assert_equal_to(_exobase.size(),_neutral_composition.n_species());

     VectorCoeffType min;
     min.resize(_neutral_composition.n_species());
     CoeffType minatm;
     unsigned int limit = _altitude.altitudes().size() * 9 / 10; //by hypothese
     if(_altitude.altitudes()[limit] < _altitude.altitudes()[_altitude.altitudes().size() / 10])
                limit = _altitude.altitudes().size() / 10;
//init
     for(unsigned int s = 0; s < _neutral_composition.n_species(); s++)
     {
        _exobase[s] = 0;
        min[s] = Antioch::ant_abs(_scale_height[s][0] - _free_path[s][0]);
     }
     _mean_exobase = 0;
      minatm = Antioch::ant_abs(_mean_scale_height[0] - _mean_free_path[0]);
     

     for(unsigned int iz = 1; iz < _altitude.altitudes().size(); iz++)
     {
        for(unsigned int s = 0; s < _neutral_composition.n_species(); s++)
        {
           if(Antioch::ant_abs(_scale_height[s][iz] - _free_path[s][iz]) < min[s])
           {
             _exobase[s] = iz;
             min[s] = Antioch::ant_abs(_scale_height[s][iz] - _free_path[s][iz]);
           }
        }
        if(Antioch::ant_abs(_mean_scale_height[iz] - _mean_free_path[iz]) < minatm)
        {
          _mean_exobase = iz;
          minatm = Antioch::ant_abs(_mean_scale_height[iz] - _mean_free_path[iz]);
        }
     }

//no more than the limit
    for(unsigned int s = 0; s < _neutral_composition.n_species(); s++)
    {
        if(_altitude.altitudes()[_exobase[s]] > _altitude.altitudes()[limit])
                _exobase[s] = limit;
    }
    if(_altitude.altitudes()[_mean_exobase] > _altitude.altitudes()[limit])
                _mean_exobase = limit;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template <typename VectorStateType>
  inline
  void AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::set_hard_sphere_radius(const VectorStateType &hsr)
  {
      antioch_assert_equal_to(hsr.size(), _neutral_composition.n_species());
      _hard_sphere_radius = hsr;
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
  void AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::make_scale_height()
  {

     _scale_height.resize(_neutral_composition.n_species());
     for(unsigned int ineu = 0; ineu < _neutral_composition.n_species(); ineu++)
     {
        _scale_height[ineu].resize(_altitude.altitudes().size(),0.L);
        for(unsigned int iz = 0; iz < _altitude.altitudes().size(); iz++)
        {
           _scale_height[ineu][iz] = this->H(_neutral_composition.M(ineu),
                                             _temperature.neutral_temperature()[iz],
                                             _altitude.altitudes()[iz]);
        }
     }

     antioch_assert_equal_to(_neutral_molar_fraction.size(),_neutral_composition.n_species());
     for(unsigned int s = 0; s < _neutral_molar_fraction.size(); s++)
     {
        antioch_assert_equal_to(_neutral_molar_fraction[s].size(),_altitude.altitudes().size());
     }

        //mean and a
     _mean_scale_height.resize(_altitude.altitudes().size(),0.L);
     _a_factor.resize(_altitude.altitudes().size(),0.L);
     for(unsigned int iz = 0; iz < _altitude.altitudes().size(); iz++)
     {
        CoeffType Mm;
        Antioch::set_zero(Mm);
        for(unsigned int s = 0; s < _neutral_composition.n_species(); s++)
        {
           Mm += _neutral_molar_fraction[s][iz] * _neutral_composition.M(s);
        }
        _mean_scale_height[iz] = this->H(Mm,_temperature.neutral_temperature()[iz],_altitude.altitudes()[iz]);
        _a_factor[iz] = (Constants::Titan::radius<CoeffType>() + _altitude.altitudes()[iz]) / _mean_scale_height[iz];
     }
     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const MatrixCoeffType &AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::scale_height() const
  {
     return _scale_height;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const VectorCoeffType &AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::total_density() const
  {
     return _total_density;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const MatrixCoeffType &AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::neutral_molar_fraction() const
  {
     return _neutral_molar_fraction;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const MatrixCoeffType &AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::ionic_molar_fraction() const
  {
     return _ionic_molar_fraction;
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
  inline
  const MatrixCoeffType &AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::free_path() const
  {
     return _free_path;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const std::vector<unsigned int> &AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::exobase() const
  {
     return _exobase;
  }



  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const VectorCoeffType &AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::atmosphere_free_path() const
  {
     return _mean_free_path;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const VectorCoeffType &AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::atmosphere_scale_height() const
  {
     return _mean_scale_height;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const VectorCoeffType &AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::a_factor() const
  {
     return _a_factor;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  unsigned int AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::atmosphere_exobase() const
  {
     return _mean_exobase;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  CoeffType AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::mean_neutral_molar_mass(unsigned int iz) const
  {
    CoeffType meanmass;
    Antioch::set_zero(meanmass);
    unsigned int n(_neutral_composition.n_species());
    for(unsigned int ispec = 0; ispec < n; ispec++)
    {
        meanmass += _neutral_molar_fraction[ispec][iz] * _neutral_composition.M(ispec);
    }

    return meanmass;
  }


  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::init_composition(const VectorStateType &bot_compo,
                                                                                       const StateType &dens_tot_bot)
  {

    _total_density.resize(_altitude.altitudes().size(),dens_tot_bot);
    _neutral_molar_fraction.resize(_neutral_composition.n_species());
    _ionic_molar_fraction.resize(bot_compo.size() - _neutral_composition.n_species());
  // molar fractions kept constant
  // first neutral, then ionic
    unsigned int bc(0);
    for(unsigned int n = 0; n < _neutral_composition.n_species(); n++)
    {
       _neutral_molar_fraction[n].resize(_altitude.altitudes().size(),bot_compo[bc]);
       bc++;
    }
    for(unsigned int n = 0; n < bot_compo.size() - _neutral_composition.n_species(); n++)
    {
       _ionic_molar_fraction[n].resize(_altitude.altitudes().size(),bot_compo[bc]);
       bc++;
    }

    for(unsigned int i = 0 ; i < _altitude.altitudes().size(); i++)
    {
  //mean
      StateType meanmass = this->mean_neutral_molar_mass(i);

      _total_density[i] = dens_tot_bot * Antioch::ant_exp(-(_altitude.altitudes()[i] - _altitude.alt_min())/   //n_tot * exp(-(z - z0) / ( 
                         (
                           (Constants::Titan::radius<StateType>() + _altitude.altitudes()[i]) *                //(r_Titan + z) *
                           (Constants::Titan::radius<StateType>() + _altitude.alt_min()) * 1e3L * //to m       //(r_Titan + z0) *
                           (Antioch::Constants::Avogadro<StateType>() * Constants::Universal::kb<StateType>() 
                                                                      * _temperature.neutral_temperature()[i] /   //(kb * T /
                                (Constants::Universal::G<StateType>() * Constants::Titan::mass<StateType>() * meanmass * 1e-3L))//g/mol to kg/mol //(G * m_Titan * <M>))))
                         )                              );

    }

    return;
  }


  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  void AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::initialize_exobases()
  {
     /// default is 9/10 alt max, this is the max allowed
     _mean_exobase = _altitude.altitudes().size() * 9 / 10;
     if(_altitude.altitudes()[_mean_exobase] < _altitude.altitudes()[_altitude.altitudes().size() / 10])
                _mean_exobase = _altitude.altitudes().size() / 10;
     _exobase.resize(_neutral_composition.n_species(),_mean_exobase);
     return;
  }

}

#endif
