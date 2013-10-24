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


       
       void derive(VectorCoeffType &deriv, const VectorCoeffType &vec, const CoeffType &dx);

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

/*
  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  CoeffType DiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::derivative(const VectorCoeffType &vec, const VectorCoeffType &dx, unsigned int iz)
  {
     return (deriv[ix] = (vec[iz+1] - vec[iz-1]) / (dx[iz+1] - dx[iz-1]));
  }
*/

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

     _omega.clear();
     _omega.resize(_mixture.neutral_composition().n_species());
     for(unsigned int s = 0; s < _mixture.neutral_composition().n_species(); s++)
     {
        _omega[s].resize(_altitude.altitudes().size(),0.L);
        for(unsigned int iz = 1; iz < _altitude.altitudes().size() - 1; iz++)
        {
            _omega[s][iz] = - _molecular_diffusion.Dtilte()[s][iz] * (
            CoeffType(1.)/(_mixture.total_density()[iz] * _mixture.neutral_molar_fraction()[s][iz]) * 
( /// dn_dz 
   _mixture.total_density()[iz+1] * _mixture.neutral_molar_fraction()[s][iz+1] - _mixture.total_density()[iz-1] * _mixture.neutral_molar_fraction()[s][iz-1] 
)  / (_altitude.altitudes()[iz+1] - _altitude.altitudes()[iz-1]) + 
        CoeffType(1.)/_mixture.scale_height()[s][iz] + CoeffType(1.)/_temperature.neutral_temperature()[iz] * 
(_temperature.neutral_temperature()[iz+1] - _temperature.neutral_temperature()[iz-1])  / (_altitude.altitudes()[iz+1] - _altitude.altitudes()[iz-1])
        * (CoeffType(1.) + (CoeffType(1.) - _mixture.neutral_molar_fraction()[s][iz]) * _mixture.thermal_coefficient()[s])
                                                ) - _eddy_diffusion.K()[iz] * (
            CoeffType(1.)/(_mixture.total_density()[iz] * _mixture.neutral_molar_fraction()[s][iz]) * 
( /// dn_dz 
   _mixture.total_density()[iz+1] * _mixture.neutral_molar_fraction()[s][iz+1] - _mixture.total_density()[iz-1] * _mixture.neutral_molar_fraction()[s][iz-1] 
)  / (_altitude.altitudes()[iz+1] - _altitude.altitudes()[iz-1]) + 
        CoeffType(1.)/_mixture.mean_scale_height()[iz] + 
CoeffType(1.)/_temperature.neutral_temperature()[iz] * 
(_temperature.neutral_temperature()[iz+1] - _temperature.neutral_temperature()[iz-1])  / (_altitude.altitudes()[iz+1] - _altitude.altitudes()[iz-1]));
        }
     }
  }
  
}

#endif
