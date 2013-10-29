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

#ifndef PLANET_DIFFUSION_EVALUATOR_H
#define PLANET_DIFFUSION_EVALUATOR_H

//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/chemical_mixture.h"

//Planet
#include "planet/molecular_diffusion_evaluator.h"
#include "planet/eddy_diffusion_evaluator.h"

//C++
#include <string>

namespace Planet{

  /*!\class DiffusionEvaluator
 * Stores all kind of diffusions
 *
 */
  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  class DiffusionEvaluator
  {
      private:
       DiffusionEvaluator() {antioch_error();return;}

       MatrixCoeffType _omega;

//dependencies
       MolecularDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> &_molecular_diffusion;
       EddyDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>      &_eddy_diffusion;
       AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>          &_mixture;
       Altitude<CoeffType,VectorCoeffType>                                    &_altitude;
       AtmosphericTemperature<CoeffType,VectorCoeffType>                      &_temperature;

       //! \returns \f$\frac{\partial n_s}{\partial z}\f$ at a given altitude for a given species
       CoeffType dn_dz(unsigned int iz, unsigned int s) const;

       //! \returns \f$\frac{\partial T}{\partial z}\f$ at a given altitude
       CoeffType dT_dz(unsigned int iz) const;

      public:
       //!
       DiffusionEvaluator(MolecularDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> &mol_diff,
                          EddyDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>      &eddy_diff,
                          AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>          &mix,
                          Altitude<CoeffType,VectorCoeffType>                                    &alt,
                          AtmosphericTemperature<CoeffType,VectorCoeffType>                      &temp);

        //!
       ~DiffusionEvaluator();

       //!
       void make_diffusion();

       //!
       void initialize();

       //!
       const MatrixCoeffType &diffusion() const;

  };

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  DiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::DiffusionEvaluator(
                          MolecularDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> &mol_diff,
                          EddyDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>      &eddy_diff,
                          AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>          &mix,
                          Altitude<CoeffType,VectorCoeffType>                                    &alt,
                          AtmosphericTemperature<CoeffType,VectorCoeffType>                      &temp):
    _molecular_diffusion(mol_diff),
    _eddy_diffusion(eddy_diff),
    _mixture(mix),
    _altitude(alt),
    _temperature(temp)
  {
     return;
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  DiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::~DiffusionEvaluator()
  {
     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const MatrixCoeffType &DiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::diffusion() const
  {
     return _omega;
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  void DiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::initialize()
  {
     _molecular_diffusion->make_molecular_diffusion();
     _eddy_diffusion->make_eddy_diffusion();
     this->make_diffusion();
     return;
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  void DiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::make_diffusion()
  {

/*     
 *     omega = - Dtilde * (1/ns * dns_dz + 1/Hs + 1/T * dT_dz * (1 + (1 - xs) * alphas)) 
 *             - K      * (1/ns * dns_dz + 1/Ha + 1/T * dT_dz) 
 */
     _omega.clear();
     _omega.resize(_mixture.neutral_composition().n_species());
     for(unsigned int s = 0; s < _mixture.neutral_composition().n_species(); s++)
     {
        _omega[s].resize(_altitude.altitudes().size(),0.L);
        for(unsigned int iz = 1; iz < _altitude.altitudes().size() - 1; iz++)
        {
            _omega[s][iz] =  - _molecular_diffusion.Dtilde()[s][iz] * // - Dtilde * (
            (
              CoeffType(1.)/(_mixture.total_density()[iz] * _mixture.neutral_molar_fraction()[s][iz]) * this->dn_dz(iz,s) // 1/ns * dns_dz
            + CoeffType(1.)/_mixture.scale_height()[s][iz]  // + 1/Hs
            + CoeffType(1.)/_temperature.neutral_temperature()[iz] * this->dT_dz(iz) // + 1/T * dT_dz * (
                * (CoeffType(1.) + (CoeffType(1.) - _mixture.neutral_molar_fraction()[s][iz]) * _mixture.thermal_coefficient()[s]) //1 + (1 - xs)*alphas ) )
            )
             - _eddy_diffusion.K()[iz] * // - K * (
            ( 
              CoeffType(1.)/(_mixture.total_density()[iz] * _mixture.neutral_molar_fraction()[s][iz]) * this->dn_dz(iz,s) // 1/ns * dns_dz
            + CoeffType(1.)/_mixture.atmosphere_scale_height()[iz]  // + 1/Ha
            + CoeffType(1.)/_temperature.neutral_temperature()[iz] * this->dT_dz(iz) //+1/T * dT_dz )
            );
        }
     }
  }
  
  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  CoeffType DiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::dn_dz(unsigned int iz, unsigned int s) const
  {
     antioch_assert_greater(iz,0);
     antioch_assert_less(iz,_altitude.altitudes().size()-1);
     antioch_assert_less(s,_mixture.neutral_composition().n_species());
     return  (_mixture.total_density()[iz+1] * _mixture.neutral_molar_fraction()[s][iz+1] - 
                                        _mixture.total_density()[iz-1] * _mixture.neutral_molar_fraction()[s][iz-1])  / 
                        (_altitude.altitudes()[iz+1] - _altitude.altitudes()[iz-1]);
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  CoeffType DiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::dT_dz(unsigned int iz) const
  {
     antioch_assert_greater(iz,0);
     antioch_assert_less(iz,_altitude.altitudes().size()-1);
     return  (_temperature.neutral_temperature()[iz+1] - _temperature.neutral_temperature()[iz-1])  
                                / (_altitude.altitudes()[iz+1] - _altitude.altitudes()[iz-1]);
  }


}

#endif
