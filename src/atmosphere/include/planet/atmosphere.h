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

#ifndef PLANET_ATMOSPHERE_H
#define PLANET_ATMOSPHERE_H

//Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/metaprogramming_decl.h"
#include "antioch/vector_utils.h"
#include "antioch/metaprogramming.h"
#include "antioch/cmath_shims.h"
#include "antioch/antioch_asserts.h"
#include "antioch/chemical_mixture.h"
#include "antioch/physical_constants.h"
#include "antioch/kinetics_evaluator.h"

//Planet
#include "planet/planet_constants.h"
#include "planet/photon_flux.h"
#include "planet/diffusion_evaluator.h"
#include "planet/absorption_grid.h"

//C++
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <iostream>

namespace Planet{

template <typename CoeffType = double, 
          typename VectorCoeffType = std::vector<CoeffType>, 
          typename MatrixCoeffType = std::vector<std::vector<CoeffType> > 
         >
class Atmosphere{
      private:

        Atmosphere(){antioch_error();return;}

//composition
        Antioch::ChemicalMixture<CoeffType> _neutral_composition;
        Antioch::ChemicalMixture<CoeffType> _ionic_composition;
//
        VectorCoeffType _total_density;
        MatrixCoeffType _neutral_molar_fraction; //alt,nneus
        MatrixCoeffType _ionic_molar_fraction;
// neutral only
        VectorCoeffType _thermal_coefficient;
// neutral only
        VectorCoeffType _hard_sphere_radius;
        VectorCoeffType _mean_free_path;
// neutral only
        MatrixCoeffType _scale_height;
        VectorCoeffType _mean_scale_height;

//kinetics, photo/electrochemistry to be swallowed
        Antioch::KineticsEvaluator<CoeffType> _neutral_reactions;
        Antioch::KineticsEvaluator<CoeffType> _ionic_reactions;
//photon
        PhotonFlux<CoeffType,VectorCoeffType,MatrixCoeffType> &_hv_flux;
        std::map<Antioch::Species, AbsorptionGrid<CoeffType,VectorCoeffType> > _photon_absorbing_species_sigma;
        std::vector<Antioch::Species >   _photon_absorbing_species;
//electron

//temperature
        VectorCoeffType _temperature;
        VectorCoeffType _ionic_temperature;

//altitude grid
        std::map<CoeffType,unsigned int> _altitudes;
        std::vector<CoeffType>           _altitudes_list;
        CoeffType _min_alt;
        CoeffType _max_alt;
        CoeffType _step_alt;

//diffusion, bimolecular + eddy
        DiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> _diffusion;

// intern functions
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        scale_height(const StateType &temp, const StateType &g, const StateType &Mm)
        ANTIOCH_AUTOFUNC(StateType, Constants::Universal::kb<StateType>() * temp / 
                                    (g * StateType(1e-3L) * Mm / Antioch::Constants::Avogadro<StateType>())
                        )

      public:

        Atmosphere(const std::vector<std::string> &neutral_spec,
                   const std::vector<std::string> &ionic_spec,
                   PhotonFlux<CoeffType,VectorCoeffType,MatrixCoeffType> &hv_flux);
        Atmosphere(const std::vector<std::string> &neutral_spec,
                   const std::vector<std::string> &ionic_spec,
                   PhotonFlux<CoeffType,VectorCoeffType,MatrixCoeffType> &hv_flux,
                   const VectorCoeffType & T,      const VectorCoeffType &bot_compo, 
                   const CoeffType & dens_tot_bot, const VectorCoeffType &alti);
        ~Atmosphere();

        //!From bottom up to top
        template<typename StateType,typename VectorStateType>
        void init_composition(const VectorStateType &bot_compo,const StateType &dens_tot_bot);

/* mean free path*/
        //! set hard sphere radius and calculate mean free pathes
        template<typename VectorStateType>
        void set_hard_sphere_radius(const VectorStateType &hsp);

        //! computes mean free pathes
        void make_free_pathes();

        //! mean free path of species s
        VectorCoeffType mean_free_path(unsigned int s) const;

/*exobase (altitude where mean free path is equal to scale height) */
        //! \return the exobase altitude for one species
        const VectorCoeffType exobase_altitude(unsigned int s) const;

        //! fills exo_altitudes with all the exobases
        void exobase_altitude(VectorCoeffType &exo_altitudes) const;

/* diffusion */
        //! \return the diffusion evaluator object
        DiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> diffusion();

        //!
        void make_diffusion();

// eddy diffusion
        //! \return eddy diffusion factor on altitude column
        const VectorCoeffType K() const;

// bimolecular
        const MatrixCoeffType Dtilde() const;

/* photon */
        //!
        template<typename VectorStateType>
        void add_photoabsorption(const std::string &name,const VectorStateType &lambda, const VectorStateType &sigma);

        //!any change in composition changes photon flux at altitude z
        void update_photochemistry();

        //!
        template<typename StateType, typename VectorStateType, typename MatrixStateType>
        void set_photon_flux(const PhotonFlux<StateType,VectorStateType,MatrixStateType> &hv_flux);

        //!
        void set_photon_absorption();

        //!
        const AbsorptionGrid<CoeffType,VectorCoeffType> photon_sigma(unsigned int is) const;

        //!
        template<typename StateType = CoeffType,typename VectorStateType = VectorCoeffType, typename MatrixStateType = MatrixCoeffType>
        const PhotonFlux<StateType,VectorStateType,MatrixStateType> hv_flux() const;


/* altitude */

        //! Compute altitude grid with parameters
        template<typename StateType>
        void make_altitude_grid(const StateType &zlow, const StateType &zhigh, const StateType &zstep);

        //! Copy given altitude grid
        template<typename VectorStateType>
        void set_altitude_grid(const VectorStateType &grid);

        //!
        const std::map<CoeffType,unsigned int> altitude_map() const;

        //!
        const VectorCoeffType altitude() const;

        //!
        CoeffType min_alt() const;

        //!
        CoeffType max_alt() const;

        //!
        CoeffType step_alt() const;

        //!
        unsigned int n_altitudes() const;

/*temperature*/
        //!
        template<typename VectorStateType>
        void set_temperature(const VectorStateType &T, const VectorStateType &altT);

/* thermal coefficient */
        //! returns thermal coefficients
        const VectorCoeffType neutral_thermal_coefficient() const;

        //!sets a thermal coefficient
        template<typename VectorStateType>
        void set_neutral_thermal_coefficient(const VectorStateType &alpha);

/* species */

        //!
        Antioch::ChemicalMixture<CoeffType> neutral_composition() const {return _neutral_composition;}

        //!
        Antioch::ChemicalMixture<CoeffType> ionic_composition()   const {return _ionic_composition;}

        //!
        template<typename StateType>
        void set_composition(Antioch::ChemicalMixture<StateType> &comp);

        //!
        template<typename VectorStateType, typename MatrixStateType>
        void set_molar_fraction(const VectorStateType &alt, const VectorStateType &dens_tot, const MatrixStateType &frac_mol);

        //!
        const std::map<Antioch::Species, AbsorptionGrid<CoeffType,VectorCoeffType> > absorbing_species_sigma() const;

//global characteristics
        //!\return temperature 
        const VectorCoeffType temperature() const;

        //! \return temperature at given altitude
        template <typename StateType>
        ANTIOCH_AUTO(StateType)
        temperature(const StateType &z) const 
        ANTIOCH_AUTOFUNC(StateType,_temperature[_altitudes.at(z)])

        //! \return pressure at given altitude
        template <typename StateType>
        ANTIOCH_AUTO(StateType)
        pressure(const StateType &z) const 
        ANTIOCH_AUTOFUNC(StateType,_total_density[_altitudes.at(z)] * 1e-6 / Antioch::Constants::Avogadro<StateType>() * //cm-3 -> m-3 -> mol/m-3
                                   Antioch::Constants::R_universal<StateType>() * 1e-3 * // J/kmol/K -> J/mol/K
                                   _temperature[_altitudes.at(z)]) // K


        //!\return total density top to bottom
        const VectorCoeffType total_density()    const; 

        //! return total density at altitude z
        template <typename StateType>
        ANTIOCH_AUTO(StateType)
        total_density(const StateType &z) const
        ANTIOCH_AUTOFUNC(StateType,_total_density[_altitudes.at(z)])


        //! \return gravity
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
         g(const StateType &z) const 
        ANTIOCH_AUTOFUNC(StateType, Constants::Universal::G<StateType>() * Constants::Titan::mass<StateType>() / 
                                    Antioch::ant_pow(StateType(1e3L)*(Constants::Titan::radius<StateType>() + z),2))

        //! \return scale height of species s at altitude z, H = kb*T/(g*Ms)
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        H(unsigned int s, const StateType &z) const
        ANTIOCH_AUTOFUNC(StateType, Constants::Universal::kb<StateType>() * this->temperature(z) / 
                                    (this->g(z) * StateType(1e-3L) * _neutral_composition.M(s) / Antioch::Constants::Avogadro<StateType>())
                        )

        //! \return scale height of species s H = kb*T/(g*Ms)
        void H(unsigned int s, VectorCoeffType &H_top_bottom) const;

        //! \return mean scale height for all altitudes
        void H_top_to_bottom(VectorCoeffType &H_top_bottom) const;

        //! \return mean scale height at altitude z, H = kb*T/(g*<m>)
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        H(const StateType &z) const
        ANTIOCH_AUTOFUNC(StateType, Constants::Universal::kb<StateType>() * this->temperature(z) / 
                                    (this->g(z) * StateType(1e-3L) * this->M(z) / Antioch::Constants::Avogadro<StateType>())
                        )

        //! \return a factor at altitude z, a = (R + z)/H
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        a(const StateType &z) const 
        ANTIOCH_AUTOFUNC(StateType,StateType(1e3L)*(Constants::Titan::radius<StateType>() + z)/this->H(z))

        //! \return molar mass of atmosphere at given altitude
        template <typename StateType>
        ANTIOCH_AUTO(StateType)
        M(const StateType &z) const
        ANTIOCH_AUTOFUNC(StateType, this->mean_neutral_molar_mass(z))

// scale height
        void make_scale_height();

// molecule
        //! \return neutral mass fraction
        template <typename StateType, typename VectorStateType>
        void neutral_mass_fraction(const StateType &z, VectorStateType &neutral_y) const;

        //! \return mean neutral molar mass at altitude z
        template <typename StateType>
        StateType mean_neutral_molar_mass(const StateType &z) const;

        //!
        template <typename StateType>
        StateType neutral_molar_fraction(const unsigned int ispec, const StateType &z) const;

        //!
        template <typename StateType>
        StateType ionic_molar_fraction(const unsigned int ispec, const StateType &z) const;

        //!
        CoeffType neutral_molar_density(unsigned int nneu, unsigned int alt) const; 

        //!
        void neutral_molar_density_bottom_to_top(unsigned int nneu, VectorCoeffType &neutral_n_bottom_top) const;

        //!
        void neutral_molar_density_top_to_bottom(unsigned int nneu, VectorCoeffType &neutral_n_top_bottom) const;

        //!
        template <typename StateType>
        ANTIOCH_AUTO(StateType)
        neutral_molar_density(unsigned int nneu, const StateType &alt) const 
        ANTIOCH_AUTOFUNC(StateType,neutral_molar_density(nneu,_altitudes.at(alt)))

        //!
        template <typename StateType>
        ANTIOCH_AUTO(StateType)
        ionic_molar_density(unsigned int nion, const StateType &alt) const 
        ANTIOCH_AUTOFUNC(StateType,this->ionic_molar_fraction(nion,alt) * this->total_density(alt))

        //!
        unsigned int n_neutral_species()          const;
        //!
        unsigned int n_ionic_species()            const;
        //!
        unsigned int n_photon_absorbing_species() const;

        //!printing methods
        void print_composition(std::ostream &out) const;
        void print_temperature(std::ostream &out) const;
        void print_photon_flux(std::ostream &out) const;

};


template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::make_scale_height()
{
    _scale_height.resize(_neutral_composition.n_species(),0.L);
    _mean_scale_height.resize(_altitudes_list.size(),0.L);

    VectorCoeffType mass_fraction;
    mass_fraction.resize(_neutral_composition.n_species(),0.L);  
    for(unsigned int s = 0; s < _neutral_composition.n_species(); s++)
    {
      
       _scale_height[s].resize(_altitudes_list.size(),0.L);
       for(unsigned int iz = 0; iz < _altitudes_list.size(); iz++)
       {
          _scale_height[s][iz]   = H(_temperature[iz],this->g(_altitudes_list[iz]),_neutral_composition.M(s));
          for(unsigned int p = 0; p < _neutral_composition.n_species(); p++)
          {
             mtot += _molar_fraction[] * _neutral_composition.M(s);
          }
          _mean_scale_height[iz] = H(_temperature[iz],this->g(_altitudes_list[iz]),_neutral_composition.M(mass_fraction));
       }
    }
    return;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
const VectorCoeffType Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::exobase_altitude(unsigned int s) const
{
  antioch_assert_lower(s,_neutral_composition.n_species());
  CoeffType alt(0.L);
  CoeffType tmp = 1e303;
  for(CoeffType z = _min_alt; z <= _max_alt; z += _step_alt)
  {
     if(Antioch::ant_abs(_mean_free_path[s][_altitudes[z]] - H(s,z)) < tmp)
     {
        alt = z;
     }
  }

  return alt;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::exobase_altitude(VectorCoeffType &exo_altitudes) const
{
  exo_altitudes.resize(_neutral_composition.n_species(),0.L);
  for(unsigned int s = 0; s < _neutral_composition.n_species(); s++)
  {
      exo_altitudes[s] = exobase_altitude(s);
  }

  return;  
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
const VectorCoeffType Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::altitude() const
{
  return _altitudes_list;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
template<typename VectorStateType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::set_hard_sphere_radius(const VectorStateType &hsp)
{
  antioch_assert_equal(hsp.size(),_neutral_composition.n_species());
         
  _hard_sphere_radius.clear();
  _hard_sphere_radius = hsp;

  this->make_free_pathes();

}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::make_free_pathes()
{
  _mean_free_path.clear();
  _mean_free_path.resize(this->n_neutral_species());

  for(unsigned int ineu = 0; ineu < _neutral_composition.n_species(); ineu++)
  {
     _mean_free_path[ineu].resize(this->n_altitudes(),1.L);
     for(unsigned int iz = 0; iz < this->n_altitudes(); iz++)
     {
       CoeffType sum(0.L);
       for(unsigned int jneu = 0; jneu <_neutral_composition.n_species(); jneu++)
       {
        sum += this->neutral_molar_density(jneu,iz) * Planet::Constants::pi<CoeffType>() * 
               Antioch::ant_pow(_hard_sphere_radius[jneu] + _hard_sphere_radius[ineu],2) * 
               Antioch::ant_sqrt(CoeffType(1.L) + _neutral_composition.M(ineu)/_neutral_composition.M(jneu));

       }
       _mean_free_path[ineu][iz] /= sum;
    }
  }

  return;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
VectorCoeffType Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::mean_free_path(unsigned int s) const
{
  return _mean_free_path[s];
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
const VectorCoeffType Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::neutral_thermal_coefficient() const
{
  return _thermal_coefficient;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
template<typename VectorStateType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::set_neutral_thermal_coefficient(const VectorStateType &alpha)
{
  antioch_assert_equal(alpha.size(),_neutral_composition.n_species());
  _thermal_coefficient.clear();
  _thermal_coefficient = alpha;
  return;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
const VectorCoeffType Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::K() const
{
  return _diffusion.eddy_coefficient_coeff();
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
const VectorCoeffType Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::Dtilde() const
{
  return _diffusion.molecular_diffusion();
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::make_diffusion()
{
  _diffusion.make_molecular_diffusion(*this);
  _diffusion.make_eddy_diffusion();
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
MolecularDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::diffusion()
{
  return _diffusion;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::neutral_molar_density_bottom_to_top(unsigned int nneu, VectorCoeffType &neutral_n_bottom_top) const
{
   neutral_n_bottom_top.resize(_altitudes.size(),0.L);
   unsigned int i(0);
   for(CoeffType z = _min_alt; z <= _max_alt; z += _step_alt)
   {
      neutral_n_bottom_top[i] = neutral_molar_density(nneu,z);
      i++;
   }

   return;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::neutral_molar_density_top_to_bottom(unsigned int nneu, VectorCoeffType &neutral_n_top_bottom) const
{
   neutral_n_top_bottom.resize(_altitudes.size(),0.L);
   unsigned int i(0);
   for(CoeffType z = _max_alt; z >= _min_alt; z -= _step_alt)
   {
      neutral_n_top_bottom[i] = neutral_molar_density(nneu,z);
      i++;
   }

   return;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::H(VectorCoeffType &H_bottom_top) const
{
  H_bottom_top.resize(_altitudes_list.size(),0.L);

  unsigned int i(0);
  for(CoeffType z = _min_alt; z <= _max_alt; z += _step_alt)
  {
     H_bottom_top[i] = H(z);
  }

  return;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::H(unsigned int s, VectorCoeffType &H_bottom_top) const
{
  H_bottom_top.resize(_altitudes_list.size(),0.L);
  unsigned int i(0);
  for(CoeffType z = _min_alt; z <= _max_alt; z += _step_alt)
  {
     H_bottom_top[i] = H(s,z);
  }

  return;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
CoeffType Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::min_alt() const
{
  return _min_alt;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
CoeffType Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::max_alt() const
{
  return _max_alt;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
CoeffType Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::step_alt() const
{
  return _step_alt;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
const std::map<CoeffType,unsigned int> Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::altitude_map() const 
{
  return _altitudes;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
const VectorCoeffType Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::temperature() const
{
  return _temperature;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
const VectorCoeffType Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::total_density() const
{
  return _total_density;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
unsigned int Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::n_altitudes() const
{
  return _altitudes_list.size();
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
CoeffType Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::altitude(unsigned int nalt) const
{
  antioch_assert_less(nalt,_altitudes_list.size());
  return _altitudes_list[nalt];
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
const std::map<Antioch::Species, AbsorptionGrid<CoeffType,VectorCoeffType> > 
  Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::absorbing_species_sigma() const
{
  return _photon_absorbing_species_sigma;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
const AbsorptionGrid<CoeffType,VectorCoeffType> 
  Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::photon_sigma(unsigned int is) const
{
  return _photon_absorbing_species_sigma.at(_photon_absorbing_species.at(is));
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
template<typename StateType, typename VectorStateType, typename MatrixStateType>
inline
const PhotonFlux<StateType,VectorStateType,MatrixStateType> 
  Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::hv_flux() const
{
  return _hv_flux;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
CoeffType Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::neutral_molar_density(unsigned int nneu, unsigned int alt) const
{
  return (_total_density[alt] * _neutral_molar_fraction[nneu][alt]);
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::Atmosphere(const std::vector<std::string> &neutral_spec,
                                  const std::vector<std::string> &ionic_spec,
                                  PhotonFlux<CoeffType,VectorCoeffType,MatrixCoeffType> &hv_flux):
  _neutral_composition(neutral_spec),
  _ionic_composition(ionic_spec),
  _hv_flux(hv_flux),
  _mol_diffusion(_neutral_composition)
{
  _hv_flux.set_atmosphere(this);
  return;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::Atmosphere(const std::vector<std::string> &neutral_spec,
                                  const std::vector<std::string> &ionic_spec,
                                  PhotonFlux<CoeffType,VectorCoeffType,MatrixCoeffType> &hv_flux,
                                  const VectorCoeffType & temp, const VectorCoeffType &bot_compo, const CoeffType &dens_tot_bot, const VectorCoeffType &alti):
_neutral_composition(neutral_spec),
_ionic_composition(ionic_spec),
_hv_flux(hv_flux),
 _mol_diffusion(_neutral_composition)
{
  _hv_flux.set_atmosphere(this);
  set_altitude_grid(alti);
  set_temperature(temp,alti);
  init_composition(bot_compo,dens_tot_bot);
  return;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::~Atmosphere()
{
  return;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
template<typename VectorStateType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::set_temperature(const VectorStateType &T, const VectorStateType &altT)
{
  antioch_assert_greater(_altitudes.size(),0);

  _temperature.clear();
  _temperature.resize(_altitudes_list.size(),0.L);
  unsigned int j(0);
  for(int ialt = (int)_altitudes_list.size() - 1; ialt >= 0 ; ialt--)//from bottom to top
  {
     while(altT[j] < _altitudes_list[ialt])
     {
        j++;
        if(j >= T.size() - 1)
        {
           j = T.size() - 1;
           break;
        }
     }
     CoeffType yup = T[j];
     CoeffType ydown = T[j-1];
     CoeffType xup = altT[j];
     CoeffType xdown = altT[j-1];
     CoeffType lin_extr = (xup == xdown)?T[j]:(yup - ydown)/(xup - xdown) * _altitudes_list[ialt] + yup - (yup - ydown)/(xup - xdown) * xup;
     _temperature[ialt] = lin_extr;
  }

  return;

}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
template<typename VectorStateType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::add_photoabsorption(const std::string &name,const VectorStateType &lambda, const VectorStateType &sigma)
{
  if(!(_neutral_composition.species_name_map().count(name) || _ionic_composition.species_name_map().count(name)))antioch_error();
  Antioch::Species species = (_neutral_composition.species_name_map().count(name))?
                                _neutral_composition.species_name_map().at(name):
                                _ionic_composition.species_name_map().at(name);

  unsigned int s = _photon_absorbing_species_sigma.size();
  _photon_absorbing_species_sigma.insert(std::make_pair(species,AbsorptionGrid<CoeffType,VectorCoeffType>(lambda,sigma)));
  _photon_absorbing_species.push_back(species);
  _photon_absorbing_species_sigma[_photon_absorbing_species[s]].y_on_x_grid(this->hv_flux<CoeffType,VectorCoeffType,MatrixCoeffType>().lambda());

  return;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
template<typename StateType, typename VectorStateType, typename MatrixStateType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::set_photon_flux(const PhotonFlux<StateType,VectorStateType,MatrixStateType> &hv_flux)
{
  _hv_flux = hv_flux;
  _hv_flux.set_atmosphere(this);
  return;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
unsigned int Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::n_neutral_species() const
{
  return _neutral_composition.n_species();
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
unsigned int Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::n_photon_absorbing_species() const
{
  return _photon_absorbing_species_sigma.size();
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
unsigned int Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::n_ionic_species() const
{
  return _ionic_composition.n_species();
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
template<typename StateType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::make_altitude_grid(const StateType &zlow, const StateType &zhigh, const StateType &zstep)
{
  _altitudes.clear();
  _altitudes_list.clear();
  unsigned int i(0);
  for(StateType z = zlow; z <= zhigh; z += zstep)//bottom to top
  {
     _altitudes[z] = i;
     _altitudes_list.push_back(z);
     i++;
  }
  _min_alt  = zlow;
  _max_alt  = zhigh;
  _step_alt = zstep;

  return;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::set_photon_absorption()
{
  for(unsigned int s = 0; s < _photon_absorbing_species_sigma.size(); s++)
  {
     _photon_absorbing_species_sigma[_photon_absorbing_species[s]].y_on_x_grid(this->hv_flux().lambda());
  }
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
template<typename VectorStateType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::set_altitude_grid(const VectorStateType &grid)
{
  _altitudes.clear();
  _altitudes_list.clear();
  _altitudes_list.resize(grid.size(),0.L);
  bool reverse = (grid.back() < grid.front());
  for(unsigned int i = 0; i < grid.size(); i++)
  {
     unsigned int j(i);
     if(reverse)j = grid.size() -1 - i;
     _altitudes[grid[j]] = j;
     _altitudes_list[j] = grid[j];
  }

  _min_alt = _altitudes_list.front();
  _max_alt = _altitudes_list.back();
  _step_alt = Antioch::ant_abs(grid[1] - grid[0]);

  return;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
template<typename StateType>
inline
StateType Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::neutral_molar_fraction(const unsigned int ispec, const StateType &z) const
{
  antioch_assert_less(ispec,this->n_neutral_species());

  return _neutral_molar_fraction[ispec][_altitudes.at(z)];
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
template<typename StateType>
inline
StateType Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::ionic_molar_fraction(const unsigned int ispec, const StateType &z) const
{
  antioch_assert_less_than(ispec,this->n_ionic_species());

  return _ionic_molar_fraction[ispec][_altitudes[z]];
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
template<typename StateType>
inline
StateType Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::mean_neutral_molar_mass(const StateType &z) const
{
  StateType meanmass;
  Antioch::set_zero(meanmass);
  unsigned int n(this->n_neutral_species());
  for(unsigned int ispec = 0; ispec < n; ispec++)
  {
      meanmass += this->neutral_molar_fraction(ispec,z) * _neutral_composition.M(ispec);
  }

  return meanmass;
}


template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
template<typename VectorStateType, typename MatrixStateType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::set_molar_fraction(const VectorStateType &alt, 
                                                                               const VectorStateType &dens_tot, 
                                                                               const MatrixStateType &frac_mol)
{
  if(_altitudes_list.empty())
  {
//sets everything
     this->set_altitude_grid(alt);
     _total_density.resize(_altitudes.size());
     _neutral_molar_fraction.resize(_neutral_composition.n_species());
     _ionic_molar_fraction.resize(_ionic_composition.n_species());
     bool reverse = (alt.front() > alt.back());
     for(unsigned int i = 0; i < alt.size(); i++)
     {
       unsigned int j = (reverse)?alt.size() - 1 - i:i;
       _total_density[i] = dens_tot[j];
       unsigned int bc(0);
       for(unsigned int n = 0; n < this->n_neutral_species(); n++)
       {
          _neutral_molar_fraction[n].push_back(frac_mol[bc][j]);
          bc++;
       }
       for(unsigned int n = 0; n < this->n_ionic_species(); n++)
       {
          _ionic_molar_fraction[n].push_back(frac_mol[bc][j]);
          bc++;
       }
     }
  }else
  {
//scale to stored altitude scale
     for(unsigned int iz = 0; iz < _altitudes_list.size(); iz++)
     {
        unsigned int jz(0);
        while( (_altitudes_list[iz] - alt[jz]) * (_altitudes_list[iz] - alt[jz+1]) > 0. &&
               (jz < alt.size() - 1) )jz++;
       
        CoeffType rel_dist = (_altitudes_list[iz] - alt[jz])/(alt[jz+1] - alt[jz]);
        _total_density[iz] = dens_tot[jz] + rel_dist * (dens_tot[jz+1] - dens_tot[jz]);
       unsigned int bc(0);
       for(unsigned int n = 0; n < this->n_neutral_species(); n++)
       {
          _neutral_molar_fraction[n][iz] = frac_mol[bc][jz] + rel_dist * (frac_mol[bc][jz+1] - frac_mol[bc][jz]);
          bc++;
       }
       for(unsigned int n = 0; n < this->n_ionic_species(); n++)
       {
          _ionic_molar_fraction[n][iz] =  frac_mol[bc][jz] + rel_dist * (frac_mol[bc][jz+1] - frac_mol[bc][jz]);
          bc++;
       }
     }
  }


}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
template<typename StateType, typename VectorStateType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::init_composition(const VectorStateType &bot_compo,const StateType &dens_tot_bot)
{
  using std::exp;
  unsigned int n_species = this->n_neutral_species() + this->n_ionic_species();

  antioch_assert_not_equal_to(_altitudes.size(),0);
  antioch_assert_equal_to(_temperature.size(),_altitudes.size());
  antioch_assert_equal_to(bot_compo.size(),n_species);


  _total_density.resize(_altitudes.size(),dens_tot_bot);
  _neutral_molar_fraction.resize(_neutral_composition.n_species());
  _ionic_molar_fraction.resize(_ionic_composition.n_species());
// molar fractions kept constant
// first neutral, then ionic
  unsigned int bc(0);
  for(unsigned int n = 0; n < this->n_neutral_species(); n++)
  {
      _neutral_molar_fraction[n].resize(_altitudes_list.size(),bot_compo[bc]);
      bc++;
  }
  for(unsigned int n = 0; n < this->n_ionic_species(); n++)
  {
      _ionic_molar_fraction[n].resize(_altitudes_list.size(),bot_compo[bc]);
      bc++;
  }

  for(unsigned int i = 0 ; i < _altitudes_list.size(); i++)//top to bottom
  {
//mean
    StateType meanmass = this->mean_neutral_molar_mass(_altitudes_list[i]);

    _total_density[i] = dens_tot_bot * exp(-(_altitudes_list[i] - _altitudes_list.back())/ //n_tot * exp(-(z - z0) / ( 
                        (
                         (Constants::Titan::radius<StateType>() + _altitudes_list[i]) *                         //(r_Titan + z) *
                         (Constants::Titan::radius<StateType>() + _altitudes_list.back()) * 1e3L * //to m       //(r_Titan + z0) *
                         (Antioch::Constants::Avogadro<StateType>() * Constants::Universal::kb<StateType>() * _temperature[i] /   //(kb * T /
                                (Constants::Universal::G<StateType>() * Constants::Titan::mass<StateType>() * meanmass * 1e-3L))//g/mol to kg/mol //(G * m_Titan * <M>))))
                        )                 );

  }

  int i = system("date \"+%S,%N\"");
  _hv_flux.update_photon_flux();
  i = system("date \"+%S,%N\"");

  return;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::print_composition(std::ostream &out) const
{
  out << "altitude dens ";
//neutral
  for(unsigned int n = 0; n < this->n_neutral_species(); n++)
  {
      out << _neutral_composition.species_inverse_name_map().at(_neutral_composition.species_list()[n]) << " ";
  }
//ionic
  for(unsigned int n = 0; n < this->n_ionic_species(); n++)
  {
      out << _ionic_composition.species_inverse_name_map().at(_ionic_composition.species_list()[n]) << " ";
  }
  out << "(km cm-3 -)" << std::endl;
  for(unsigned int nalt = 0; nalt < _altitudes.size(); nalt++)
  {
     out << _altitudes_list[nalt] << " " << _total_density[nalt] << " ";
//neutral
    for(unsigned int n = 0; n < this->n_neutral_species(); n++)
    {
        out << _neutral_molar_fraction[n][nalt] << " ";
    }
//ionic
    for(unsigned int n = 0; n < this->n_ionic_species(); n++)
    {
        out << _ionic_molar_fraction[n][nalt] << " ";
    }
    out << std::endl;
  }
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::print_temperature(std::ostream &out) const
{
  out << "Temperature altitudes (K km)" << std::endl; 
  for(unsigned int nalt = 0; nalt < _altitudes.size(); nalt++)
  {
    out << _temperature[nalt] << " " << _altitudes_list[nalt] << std::endl;
  }
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::print_photon_flux(std::ostream &out) const
{
  out << "Altitude lambda photon_flux (K nm W/m2/nm)" << std::endl; 
  _hv_flux.update_photon_flux();
  VectorCoeffType lambda = _hv_flux.lambda();
  for(unsigned int nalt = 0; nalt < _altitudes.size(); nalt++)
  {
    VectorCoeffType phy = _hv_flux.phy(_altitudes_list[nalt]);
    for(unsigned int il = 0; il < lambda.size(); il++)
    {
      out << _altitudes_list[nalt] << " " << lambda[il] << " " << phy[il] << std::endl;
    }
  }
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
template<typename StateType, typename VectorStateType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::neutral_mass_fraction(const StateType &z, VectorStateType &neutral_y) const
{
  neutral_y.resize(_neutral_molar_fraction.at(_altitudes.at(z)).size(),0.L);
  StateType mass_sum;
  Antioch::set_zero(mass_sum);
  for(unsigned int s = 0; s < this->n_neutral_species(); s++)
  {
    neutral_y[s] = _neutral_molar_fraction.at(_altitudes.at(z)).at(s) * this->_neutral_composition.M(s);
    mass_sum += neutral_y[s];
  }
  for(unsigned int s = 0; s < this->n_neutral_species(); s++)
  {
    neutral_y[s] /= mass_sum;
  }

  return;
}

} //namespace Planet

#endif
