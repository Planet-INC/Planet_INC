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
#include "antioch/metaprogramming.h"
#include "antioch/cmath_shims.h"
#include "antioch/antioch_asserts.h"
#include "antioch/physical_constants.h"

//Planet
#include "planet/altitude.h"
#include "planet/planet_constants.h"
#include "planet/photon_evaluator.h"
#include "planet/diffusion_evaluator.h"
#include "planet/atmospheric_mixture.h"
#include "planet/atmospheric_kinetics.h"
#include "planet/atmospheric_temperature.h"

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

//caracteristics
        VectorCoeffType _exobase;


//dependencies

//altitude grid
        Altitude<CoeffType,VectorCoeffType> &_altitude;

//composition
        AtmosphericMixture<CoeffType,VectorCoeffType> &_composition;

//kinetics, photo/electrochemistry swallowed
        AtmosphericKinetics<CoeffType> &_kinetics;

//temperature
        AtmosphericTemperature<CoeffType,VectorCoeffType> &_temperature;

//photon
        PhotonEvaluator<CoeffType,VectorCoeffType> &_photon_eval;
//electron

//diffusion, bimolecular + eddy
        DiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> &_diffusion;


      public:

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

        //! \return a factor at altitude z, a = (R + z)/H
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        a(const StateType &z) const 
        ANTIOCH_AUTOFUNC(StateType,StateType(1e3L)*(Constants::Titan::radius<StateType>() + z)/this->H(z))

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
Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::~Atmosphere()
{
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
