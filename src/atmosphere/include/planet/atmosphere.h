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

//dependencies

//altitude grid
        Altitude<CoeffType,VectorCoeffType> &_altitude;
//temperature
        AtmosphericTemperature<CoeffType,VectorCoeffType> &_temperature;
//composition
        AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &_composition;
//kinetics, photo/electrochemistry swallowed in Antioch
        AtmosphericKinetics<CoeffType> &_kinetics;
//photon
        PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> &_photon_eval;
//electron
//        ElectronEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> &_electron_eval;
//diffusion, bimolecular + eddy
        DiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> &_diffusion;

        //!initialize everything
        void initialize_atmosphere();

      public:
        //!
        Atmosphere(Altitude<CoeffType,VectorCoeffType> &alt,
                   AtmosphericTemperature<CoeffType,VectorCoeffType> &temp,
                   AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &comp,
                   AtmosphericKinetics<CoeffType> &kin,
                   PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> &hv_flux,
                   DiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> &diff);

        //!
        ~Atmosphere();

        //! \return pressure at altitude iz, p = nRT
        CoeffType pressure(unsigned int iz) const;

        //! ideal gas equation
        //
        //BEWARE, once Antioch unit merged, R in SI!!!
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        p_ideal_gas(const StateType &n, const StateType &T) const
        ANTIOCH_AUTOFUNC(StateType, n * Antioch::Constants::R_universal<CoeffType>() * CoeffType(1e-3) * T) // J/kmol/K -> J/mol/K

        //!printing composition
        void print_composition(std::ostream &out) const;

        //!printing temperature
        void print_temperature(std::ostream &out) const;

        //!printing photon
        void print_photon_flux(std::ostream &out) const;
        
};

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
CoeffType Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::pressure(unsigned int iz) const
{
  return p_ideal_gas(_composition.total_density()[iz] * CoeffType(1e-6) / Antioch::Constants::Avogadro<CoeffType>(), * //cm-3 -> m-3 -> mol/m-3
                     _temperature.altitudes()[iz]);
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::Atmosphere(Altitude<CoeffType,VectorCoeffType> &alt,
                                  AtmosphericTemperature<CoeffType,VectorCoeffType> &temp,
                                  AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &comp,
                                  AtmosphericKinetics<CoeffType> &kin,
                                  PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> &hv_flux,
                                  DiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> &diff):
        _altitude(alt),
        _temperature(temp),
        _composition(comp),
        _kinetics(kin),
        _photon_eval(hv_flux),
        _diffusion(diff)
{
  initialize_atmosphere();
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
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::initialize_atmosphere()
{
   _temperature.initialize();
   _composition.initialize();
   _kinetics.initialize();
   _photon_eval.initialize();
   _diffusion.initialize();
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::print_composition(std::ostream &out) const
{
  out << "altitude dens ";
//neutral
  for(unsigned int n = 0; n < _composition.neutral_composition().n_species(); n++)
  {
      out << _composition.neutral_composition().species_inverse_name_map().at(_composition.neutral_composition().species_list()[n]) << " ";
  }
//ionic
  for(unsigned int n = 0; n < _composition.ionic_composition().n_species(); n++)
  {
      out << _composition.ionic_composition().species_inverse_name_map().at(_composition.ionic_composition().species_list()[n]) << " ";
  }
  out << "(km cm-3 -)" << std::endl;
  for(unsigned int nalt = 0; nalt < _altitude.altitudes().size(); nalt++)
  {
     out << _altitude.altitudes()[nalt] << " " << _composition.total_density()[nalt] << " ";
//neutral
    for(unsigned int n = 0; n < _composition.neutral_composition().n_species(); n++)
    {
        out << _composition.neutral_molar_fraction()[n][nalt] << " ";
    }
//ionic
    for(unsigned int n = 0; n < _composition.ionic_composition().n_species(); n++)
    {
        out << _composition.ionic_molar_fraction()[n][nalt] << " ";
    }
    out << std::endl;
  }
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::print_temperature(std::ostream &out) const
{
  out << "Temperature altitudes (K km)" << std::endl; 
  for(unsigned int nalt = 0; nalt < _altitude.altitudes().size(); nalt++)
  {
    out << _temperature.neutral_temperature()[nalt] << " " << _altitude.altitudes()[nalt] << std::endl;
  }
}

template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
inline
void Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType>::print_photon_flux(std::ostream &out) const
{
  out << "Altitude lambda photon_flux (K nm W/m2/nm)" << std::endl; 
  for(unsigned int nalt = 0; nalt < _altitude.altitudes().size(); nalt++)
  {
    for(unsigned int il = 0; il < _photon_eval.photon_flux()[nalt].abscissa().size(); il++)
    {
      out << _altitude.altitudes()[nalt] << " " 
          << _photon_eval.photon_flux()[nalt].abscissa()[il] << " " 
          << _photon_eval.photon_flux()[nalt].flux()[il] << std::endl;
    }
  }
}

} //namespace Planet

#endif
