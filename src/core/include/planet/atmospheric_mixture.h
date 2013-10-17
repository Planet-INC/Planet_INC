//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef PLANET_ATMOSPHERIC_MIXTURE_H
#define PLANET_ATMOSPHERIC_MIXTURE_H

//Antioch
#include "antioch/chemical_mixture.h"
#include "antioch/cmath_shims.h"

//Planet
#include "planet/altitude.h"

//C++
#include <vector>

namespace Planet
{
   template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
   AtmosphericMixture
   {
      private:
        //!no default constructor
        ChemicalAtmosphere(){antioch_error();return;}

        Antioch::ChemicalMixture<CoeffType> &_neutral_composition;
        Antioch::ChemicalMixture<CoeffType> &_ionic_composition;

        VectorCoeffType _total_density;
        MatrixCoeffType _neutral_molar_fraction; //alt,nneus
        MatrixCoeffType _ionic_molar_fraction;
// neutral only
        VectorCoeffType _thermal_coefficient;
// neutral only
        VectorCoeffType _hard_sphere_radius;
        MatrixCoeffType _mean_free_path;
        MatrixCoeffType _scale_height;
        VectorCoeffType _mean_scale_height;

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
      public:
        ChemicalAtmosphere(Altitude<CoeffType,VectorCoeffType> &alt, AtmosphericTemperature<CoeffType,VectorCoeffType> &temp);
        ~ChemicalAtmosphere();

        //!
        template<typename VectorStateType>
        void init_composition(const VectorStateType &bot_compo,const StateType &dens_tot_bot);

        //!make scale heights
        void update_scale_height();

        //!\return scales height
        const MatrixCoeffType scale_height();

        //! \return mean neutral molar mass at altitude z
        //
        //Only neutral species
        CoeffType mean_neutral_molar_mass(unsigned int iz) const;

        //!sets the hard sphere radius
        template <typename VectorStateType>
        void set_hard_sphere_radius(const VectorStateType &hsr);

        //!sets the thermal coefficient
        template <typename VectorStateType>
        void set_thermal_coefficient(const VectorStateType &tc);

   };


  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  ChemicalAtmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::ChemicalAtmosphere(Altitude<CoeffType,VectorCoeffType> &alt, 
                                                                                    AtmosphericTemperature<CoeffType,VectorCoeffType> &temp):
  _altitude(alt),
  _temperature(temp)
  {
    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  void AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::make_free_pathes()
  {
    _mean_free_path.clear();
    _mean_free_path.resize(_neutral_composition.n_species());

    for(unsigned int ineu = 0; ineu < _neutral_composition.n_species(); ineu++)
    {
       _mean_free_path[ineu].resize(_altitude.altitudes().size(),1.L);
       for(unsigned int iz = 0; iz < _altitude.altitudes().size(); iz++)
       {
         CoeffType sum(0.L);
         for(unsigned int jneu = 0; jneu <_neutral_composition.n_species(); jneu++)
         {
           sum += _neutral_molar_density[jneu][iz] * Planet::Constants::pi<CoeffType>() * 
                  Antioch::ant_pow(_hard_sphere_radius[jneu] + _hard_sphere_radius[ineu],2) * 
                  Antioch::ant_sqrt(CoeffType(1.L) + _neutral_composition.M(ineu)/_neutral_composition.M(jneu));
         }
         _mean_free_path[ineu][iz] /= sum;
      }
    }

    return;
  }


  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template <typename VectorStateType>
  inline
  void AtmosphericMixture<Coefftype,VectorCoeffType,MatrixCoeffType>::set_hard_sphere_radius(const VectorStateType &hsr)
  {
      antioch_assert_equal_to(hsr.size(), _neutral_composition.n_species());

      _hard_sphere_radius = hsr;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template <typename VectorStateType>
  inline
  void AtmosphericMixture<Coefftype,VectorCoeffType,MatrixCoeffType>::set_thermal_coefficient(const VectorStateType &tc);
  {
      antioch_assert_equal_to(tc.size(), _neutral_composition.n_species());
      _thermal_coefficient = tc;
  }


  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  void AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::make_scale_height()
  {
     antioch_assert_not_equal_to(_temperature,NULL);
     antioch_assert_not_equal_to(_altitude,NULL);

     if(_scale_height.empty())_scale_height.resize(_neutral_composition.n_species());
     if(_mean_scale_height.empty())_mean_scale_height.resize(_neutral_composition.n_species());
     for(unsigned int ineu = 0; ineu < _neutral_composition.n_species(); ineu++)
     {
        if(_scale_height[ineu].empty())_scale_height[ineu].resize(_altitude.altitudes().size(),0.L);
        for(unsigned int iz = 0; iz < _altitude.altitudes().size(); iz++)
        {
           _scale_height[ineu][iz] = this->H(_neutral_composition.M(ineu), 
                                             _temperature->neutral_temperature()[iz],
                                             _altitude.altitudes()[iz]);
        }
     }

     antioch_assert_not_equal_to(_molar_fraction.size(),0);
        //mean
     for(unsigned int iz = 0; iz < _altitude.altitudes().size(); iz++)
     {
        CoeffType Mm(0.L);
        for(unsigned int s = 0; s < _neutral_composition.n_species(); s++)
        {
           Mm += _molar_fraction[s][iz] * _neutral_composition.M(s);
        }
        _mean_scale_height[iz] = this->H(Mm,_temperature->neutral_temperature()[iz],_altitude.altitudes[iz]);
     }
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const MatrixCoeffType AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>::scale_height()
  {
     return _scale_height;
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
    unsigned int n_species = _neutral_composition.n_species() + _ionic_composition.n_species();

    antioch_assert_not_equal_to(_altitude.altitudes().size(),0);
    antioch_assert_equal_to(_temperature.neutral_temperature().size(),_altitude.altitudes().size());
    antioch_assert_equal_to(bot_compo.size(),n_species);


    _total_density.resize(_altitude.altitudes().size(),dens_tot_bot);
    _neutral_molar_fraction.resize(_neutral_composition.n_species());
    _ionic_molar_fraction.resize(_ionic_composition.n_species());
  // molar fractions kept constant
  // first neutral, then ionic
    unsigned int bc(0);
    for(unsigned int n = 0; n < _neutral_composition.n_species(); n++)
    {
       _neutral_molar_fraction[n].resize(_altitude.altitudes().size(),bot_compo[bc]);
       bc++;
    }
    for(unsigned int n = 0; n < _ionic_composition.n_species(); n++)
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


}

#endif
