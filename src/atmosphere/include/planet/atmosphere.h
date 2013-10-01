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

//Planet
#include "planet/planet_constants.h"
#include "planet/photon_flux.h"
#include "planet/molecular_diffusion_evaluator.h"
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

        VectorCoeffType _temperature;
        VectorCoeffType _ionic_temperature;
        Antioch::ChemicalMixture<CoeffType> _neutral_composition;
        Antioch::ChemicalMixture<CoeffType> _ionic_composition;
//
        std::map<CoeffType,unsigned int> _altitudes;
        std::vector<CoeffType>           _altitudes_list;
        CoeffType _min_alt;
        CoeffType _max_alt;
        CoeffType _step_alt;
//
        VectorCoeffType _total_density;
        MatrixCoeffType _neutral_molar_fraction; //alt,nneus
        MatrixCoeffType _ionic_molar_fraction;
//photon
        PhotonFlux<CoeffType,VectorCoeffType,MatrixCoeffType> &_hv_flux;
        std::map<Antioch::Species, AbsorptionGrid<CoeffType,VectorCoeffType> > _photon_absorbing_species_sigma;
        std::vector<Antioch::Species >   _photon_absorbing_species;
//electron
//diffusion
        MolecularDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> _diffusion; //molecular diffusion
        CoeffType _K0;          //eddy diffusion
        VectorCoeffType _K;     //eddy diffusion
//thermal
        VectorCoeffType _thermal_coefficient;
//mean free path
        VectorCoeffType _hard_sphere_radius;
        VectorCoeffType _mean_free_path;

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

        //! mean free path of species s at altitude z
        CoeffType mean_free_path(unsigned int s, const CoeffType &z) const;

        //! mean free path of species s
        VectorCoeffType mean_free_path_top_to_bottom(unsigned int s) const;

        //! mean free path of species s
        void mean_free_path_bottom_to_top(unsigned int s, VectorCoeffType &mfp_top_bottom) const;

        //! calculate mean free pathes
        void make_free_pathes();

/*exobase (altitude where mean free path is equal to scale height) */
        //!
        const VectorCoeffType exobase_altitude(unsigned int s) const;

        //!
        void exobase_altitude(VectorCoeffType &exo_altitudes) const;

/* diffusion */
        //! return the diffusion evaluator object
        MolecularDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> diffusion();

        //!
        void make_diffusion();

        //! diffusion of species ineu on an altitude column
        void diffusion_species_bottom_to_top(unsigned int ineu, VectorCoeffType &diff_bottom_top) const;

        //! diffusion of species ineu on an altitude column
        const VectorCoeffType diffusion_species_top_to_bottom(unsigned int ineu) const;

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

/* altitude */

        //!
        template<typename StateType>
        void make_altitude_grid(const StateType &zlow, const StateType &zhigh, const StateType &zstep);

        //!
        template<typename VectorStateType>
        void set_altitude_grid(const VectorStateType &grid);

        //!
        template<typename VectorStateType>
        void set_temperature(const VectorStateType &T, const VectorStateType &altT);

        //!
        const std::map<CoeffType,unsigned int> altitude_map() const;

        //!
        void altitude_bottom_to_top(VectorCoeffType &alt_bottom_top) const;

        //!
        const VectorCoeffType altitude_top_to_bottom() const;

        //!
        CoeffType altitude(unsigned int nalt) const;

        //!
        CoeffType min_alt() const;

        //!
        CoeffType max_alt() const;

        //!
        CoeffType step_alt() const;

        //!
        unsigned int n_altitudes() const;

/* thermal coefficient */
        const VectorCoeffType neutral_thermal_coefficient() const;

        template<typename VectorStateType>
        void set_neutral_thermal_coefficient(const VectorStateType &alpha);

/* eddy diffusion*/
        template<typename StateType>
        void set_K0(const StateType &K0);

        void make_thermal_coefficient();

        //! \return thermal factor on altitude column
        void K_bottom_to_top(VectorCoeffType &K_bottom_top) const;

        //! \return thermal factor on altitude column
        const VectorCoeffType K_top_to_bottom() const;

        //! \return thermal factor at altitude z, K = K0 * sqrt(ntot(z)/ntot_bottom)
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        K(const StateType &z) const 
        ANTIOCH_AUTOFUNC(StateType,_K[_altitudes.at(z)])

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

        //!
        const AbsorptionGrid<CoeffType,VectorCoeffType> photon_sigma(unsigned int is) const;

        //!
        template<typename StateType = CoeffType,typename VectorStateType = VectorCoeffType, typename MatrixStateType = MatrixCoeffType>
        const PhotonFlux<StateType,VectorStateType,MatrixStateType> hv_flux() const;

//global characteristics
        //!\return temperature top to bottom
        const VectorCoeffType temperature_top_to_bottom()      const;

        //!\return temperature bottom to top
        void temperature_bottom_to_top(VectorCoeffType &temp_bottom_top)      const;

        //!\return total density top to bottom
        const VectorCoeffType total_density_top_to_bottom()    const; 

        //!\return total density bottom to top
        void total_density_bottom_to_bottom(VectorCoeffType &ntot_bottom_top) const; 

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
        void H_top_to_bottom(unsigned int s, VectorCoeffType &H_top_bottom) const;

        //! \return scale height of species s at altitude z, H = kb*T/(g*Ms)
        void H_bottom_to_top(unsigned int s, VectorCoeffType &H_bottom_top) const;

        //! \return mean scale height for all altitudes
        void H_bottom_to_top(VectorCoeffType &H_bottom_top) const;

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

        //! \return temperature at given altitude
        template <typename StateType>
        ANTIOCH_AUTO(StateType)
        temperature(const StateType &z) const 
        ANTIOCH_AUTOFUNC(StateType,_temperature[_altitudes.at(z)])

        //! \return temperature at given altitude
        template <typename StateType>
        ANTIOCH_AUTO(StateType)
        pressure(const StateType &z) const 
        ANTIOCH_AUTOFUNC(StateType,_total_density[_altitudes.at(z)] * 1e-6 / Antioch::Constants::Avogadro<StateType>() * //cm-3 -> m-3 -> mol/m-3
                                   Antioch::Constants::R_universal<StateType>() * 1e-3 * // J/kmol/K -> J/mol/K
                                   _temperature[_altitudes.at(z)]) // K

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
       template <typename StateType>
        ANTIOCH_AUTO(StateType)
        total_density(const StateType &z) const
        ANTIOCH_AUTOFUNC(StateType,_total_density[_altitudes.at(z)])

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
const VectorCoeffType Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::diffusion_species_top_to_bottom(unsigned int ineu) const
{
  return _diffusion.Dtilde_top_to_bottom(ineu);
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::diffusion_species_bottom_to_top(unsigned int ineu, VectorCoeffType &diff_bottom_top) const
{
  return _diffusion.Dtilde_bottom_to_top(ineu,diff_bottom_top);
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
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::altitude_bottom_to_top(VectorCoeffType &alt_bottom_top) const
{
  alt_bottom_top.resize(this->n_altitudes(),0.L);

  for(unsigned int iz = 0; iz < _altitudes_list.size(); iz++)
  {
     alt_bottom_top[_altitudes_list.size() - 1 - iz] = _altitudes_list[iz];
  }

  return;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
const VectorCoeffType Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::altitude_top_to_bottom() const
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
CoeffType Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::mean_free_path(unsigned int s, const CoeffType &z) const
{
  return _mean_free_path[s][_altitudes[z]];
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
VectorCoeffType Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::mean_free_path_top_to_bottom(unsigned int s) const
{
  return _mean_free_path[s];
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::mean_free_path_bottom_to_top(unsigned int s, VectorCoeffType &mfp_bottom_top) const
{
  mfp_bottom_top.resize(this->n_altitudes(),0.L);

  for(unsigned int iz = 0; iz < _mean_free_path[s].size(); iz++)
  {
     mfp_bottom_top[_mean_free_path[s].size() - 1 - iz] = _mean_free_path[s][iz];
  }

  return;
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
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::K_bottom_to_top(VectorCoeffType &K_bottom_top) const
{
  K_bottom_top.resize(_K.size(),0.L);
  for(unsigned int iz = 0; iz < _K.size(); iz++)
  {
     K_bottom_top[_K.size() - 1 - iz] = _K[iz];
  }

  return;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
const VectorCoeffType Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::K_top_to_bottom() const
{
  return _K;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::make_diffusion()
{
  _diffusion.make_diffusion(*this);
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
template<typename StateType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::set_K0(const StateType &K0)
{
  _K0 = K0;
  return;
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
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::H_bottom_to_top(VectorCoeffType &H_bottom_top) const
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
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::H_top_to_bottom(VectorCoeffType &H_top_bottom) const
{
  H_top_bottom.resize(_altitudes_list.size(),0.L);
  unsigned int i(0);
  for(CoeffType z = _max_alt; z >= _min_alt; z -= _step_alt)
  {
     H_top_bottom[i] = H(z);
  }

  return;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::H_bottom_to_top(unsigned int s, VectorCoeffType &H_bottom_top) const
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
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::H_top_to_bottom(unsigned int s, VectorCoeffType &H_top_bottom) const
{
  H_top_bottom.resize(_altitudes_list.size(),0.L);
  unsigned int i(0);
  for(CoeffType z = _max_alt; z >= _min_alt; z -= _step_alt)
  {
     H_top_bottom[i] = H(s,z);
  }

  return;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::make_thermal_coefficient()
{
  _K.resize(this->n_altitudes());
  for(CoeffType z = this->min_alt(); z <= this->max_alt(); z += this->step_alt())
  {
     _K[_altitudes[z]] = _K0 * Antioch::ant_sqrt(this->total_density(z)/this->total_density(_min_alt));
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
const VectorCoeffType Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::temperature_top_to_bottom() const
{
  return _temperature;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::temperature_bottom_to_top(VectorCoeffType &temp_bottom_top) const
{
  temp_bottom_top.resize(_temperature.size(),0.L);
  for(int ialt = _temperature.size() - 1; ialt >= 0; ialt--)
  {
     temp_bottom_top[ialt] = _temperature[ialt];
  }

  return;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
const VectorCoeffType Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::total_density_top_to_bottom() const
{
  return _total_density;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::total_density_bottom_to_bottom(VectorCoeffType &ntot_bottom_top) const
{
  ntot_bottom_top.resize(_total_density.size(),0.L);
  for(int ialt = _total_density.size() - 1; ialt >= 0; ialt--)
  {
     ntot_bottom_top[ialt] = _total_density[ialt];
  }

  return;
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
  return (_total_density[alt] * _neutral_molar_fraction[alt][nneu]);
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::Atmosphere(const std::vector<std::string> &neutral_spec,
                                  const std::vector<std::string> &ionic_spec,
                                  PhotonFlux<CoeffType,VectorCoeffType,MatrixCoeffType> &hv_flux):
  _neutral_composition(neutral_spec),
  _ionic_composition(ionic_spec),
  _hv_flux(hv_flux),
  _diffusion(_neutral_composition)
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
 _diffusion(_neutral_composition)
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

/*template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
template<typename VectorStateType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::set_ionic_temperature(VectorStateType &temp)
{
  antioch_assert_equal_to(temp.size(),_ionic_altitudes.size());
  _ionic_temperature = temp;
}*/

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
  for(StateType z = zhigh; z >= zlow; z -= zstep)//top to bottom
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
  bool reverse = (grid.back() > grid.front());
  for(unsigned int i = 0; i < grid.size(); i++)
  {
     unsigned int j(i);
     if(reverse)j = grid.size() -1 - i;
     _altitudes[grid[j]] = j;
     _altitudes_list[j] = grid[j];
  }

  _min_alt = _altitudes_list.back();
  _max_alt = _altitudes_list.front();
  _step_alt = Antioch::ant_abs(grid[1] - grid[0]);

  return;
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
template<typename StateType>
inline
StateType Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::neutral_molar_fraction(const unsigned int ispec, const StateType &z) const
{
  antioch_assert_less(ispec,this->n_neutral_species());

  return _neutral_molar_fraction[_altitudes.at(z)][ispec];
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
template<typename StateType>
inline
StateType Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::ionic_molar_fraction(const unsigned int ispec, const StateType &z) const
{
  antioch_assert_less_than(ispec,this->n_ionic_species());
  return _ionic_molar_fraction[_altitudes[z]][ispec];
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
     _neutral_molar_fraction.resize(_altitudes.size());
     _ionic_molar_fraction.resize(_altitudes.size());
     bool reverse = (alt.front() > alt.back());
     for(unsigned int i = 0; i < alt.size(); i++)
     {
       unsigned int j = (reverse)?alt.size() - 1 - i:i;
       _total_density[i] = dens_tot[j];
       unsigned int bc(0);
       for(unsigned int n = 0; n < this->n_neutral_species(); n++)
       {
          _neutral_molar_fraction[i].push_back(frac_mol[bc][j]);
          bc++;
       }
       for(unsigned int n = 0; n < this->n_ionic_species(); n++)
       {
          _ionic_molar_fraction[i].push_back(frac_mol[bc][j]);
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
          _neutral_molar_fraction[iz][n] = frac_mol[bc][jz] + rel_dist * (frac_mol[bc][jz+1] - frac_mol[bc][jz]);
          bc++;
       }
       for(unsigned int n = 0; n < this->n_ionic_species(); n++)
       {
          _ionic_molar_fraction[iz][n] =  frac_mol[bc][jz] + rel_dist * (frac_mol[bc][jz+1] - frac_mol[bc][jz]);
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
  _neutral_molar_fraction.resize(_altitudes.size());
  _ionic_molar_fraction.resize(_altitudes.size());
  for(unsigned int i = 0 ; i < _altitudes.size(); i++)//top to bottom
  {
// molar fractions kept constant
// first neutral, then ionic
    unsigned int bc(0);
    for(unsigned int n = 0; n < this->n_neutral_species(); n++)
    {
        _neutral_molar_fraction[i].push_back(bot_compo[bc]);
        bc++;
    }
    for(unsigned int n = 0; n < this->n_ionic_species(); n++)
    {
        _ionic_molar_fraction[i].push_back(bot_compo[bc]);
        bc++;
    }

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
        out << _neutral_molar_fraction[nalt][n] << " ";
    }
//ionic
    for(unsigned int n = 0; n < this->n_ionic_species(); n++)
    {
        out << _ionic_molar_fraction[nalt][n] << " ";
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
