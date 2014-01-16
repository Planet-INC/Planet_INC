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

//dependencies
       MolecularDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> &_molecular_diffusion;
       EddyDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>      &_eddy_diffusion;
       AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>          &_mixture;
       AtmosphericTemperature<CoeffType,VectorCoeffType>                      &_temperature;

      public:
       //!
       DiffusionEvaluator(MolecularDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> &mol_diff,
                          EddyDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>      &eddy_diff,
                          AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>          &mix,
                          AtmosphericTemperature<CoeffType,VectorCoeffType>                      &temp);

        //!
       ~DiffusionEvaluator();

       //!
       template<typename StateType, typename VectorStateType>
       void diffusion(const VectorStateType &molar_concentrations,
                      const VectorStateType &dmolar_concentrations_dz,
                      const StateType &z, VectorStateType &omegas) const;

  };

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  DiffusionEvaluator<CoeffType, VectorCoeffType,MatrixCoeffType>::DiffusionEvaluator(
                          MolecularDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> &mol_diff,
                          EddyDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>      &eddy_diff,
                          AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>          &mix,
                          AtmosphericTemperature<CoeffType,VectorCoeffType>                      &temp):
    _molecular_diffusion(mol_diff),
    _eddy_diffusion(eddy_diff),
    _mixture(mix),
    _temperature(temp)
  {
     return;
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  DiffusionEvaluator<CoeffType, VectorCoeffType,MatrixCoeffType>::~DiffusionEvaluator()
  {
     return;
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void DiffusionEvaluator<CoeffType, VectorCoeffType,MatrixCoeffType>::diffusion(const VectorStateType &molar_concentrations,
                                                                 const VectorStateType &dmolar_concentrations_dz,
                                                                 const StateType &z, VectorStateType &omegas) const
  {

     antioch_assert_equal_to(molar_concentrations.size(),_mixture.neutral_composition().n_species());
     antioch_assert_equal_to(dmolar_concentrations_dz.size(),_mixture.neutral_composition().n_species());

// Dtilde
     VectorCoeffType molecular;
     _molecular_diffusion.Dtilde(molar_concentrations,z,molecular);// Dtilde

// nTot
     CoeffType nTot(0.L);
     for(unsigned int s = 0; s < molar_concentrations.size(); s++)
     {
        nTot += molar_concentrations[s];
     }

// scale heights
     VectorCoeffType Hs;
     _mixture.scale_heights(z,Hs);
     //mean
     CoeffType Ha = _mixture.atmospheric_scale_height(molar_concentrations,z);

// temperature
     CoeffType T = _temperature.neutral_temperature(z);
     CoeffType dT_dz = _temperature.dneutral_temperature_dz(z);

     omegas.resize(_mixture.neutral_composition().n_species(),0.L);
// eddy diff
     CoeffType eddy_K = _eddy_diffusion.K(nTot);

     for(unsigned int s = 0; s < _mixture.neutral_composition().n_species(); s++)
     {
            omegas[s] =  //omega = Dtilde * [
            - molecular[s] * 
            (
                dmolar_concentrations_dz[s]/molar_concentrations[s] // 1/ns * dns_dz
              + CoeffType(1.L)/Hs[s]  // + 1/Hs
              + dT_dz/T // + 1/T * dT_dz * (
                * (CoeffType(1.L) + ((nTot - molar_concentrations[s])/nTot) * _mixture.thermal_coefficient()[s]) //1 + (1 - xs)*alphas ) ]
            )
             - eddy_K * // - K * (
            ( 
                dmolar_concentrations_dz[s]/molar_concentrations[s] // 1/ns * dns_dz
              + CoeffType(1.L)/Ha // + 1/Ha
              + dT_dz/T //+1/T * dT_dz )
            );
     }
     return;
  }


}

#endif
