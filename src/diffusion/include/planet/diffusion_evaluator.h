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
       const MolecularDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> & _molecular_diffusion;
       const EddyDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>      & _eddy_diffusion;
       const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>          & _mixture;
       const AtmosphericTemperature<CoeffType,VectorCoeffType>                      & _temperature;

      public:
       //!
       DiffusionEvaluator(const MolecularDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> &mol_diff,
                          const EddyDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>      &eddy_diff,
                          const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>          &mix,
                          const AtmosphericTemperature<CoeffType,VectorCoeffType>                      &temp);

        //!
       ~DiffusionEvaluator();

       //!
       template<typename StateType, typename VectorStateType>
       void diffusion(const VectorStateType &molar_concentrations,
                      const VectorStateType &dmolar_concentrations_dz,
                      const StateType &z, VectorStateType &omegas) const;

       //!
       template<typename StateType, typename VectorStateType, typename MatrixStateType>
       void diffusion_and_derivs(const VectorStateType &molar_concentrations,
                                 const VectorStateType &dmolar_concentrations_dz,
                                 const StateType &z, 
                                 VectorStateType &omegas_A_term,
                                 VectorStateType &omegas_B_term,
                                 MatrixStateType &domegas_dn_i_A_TERM,
                                 MatrixStateType &domegas_dn_i_B_TERM) const;

  };

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  DiffusionEvaluator<CoeffType, VectorCoeffType,MatrixCoeffType>::DiffusionEvaluator(
                          const MolecularDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> &mol_diff,
                          const EddyDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>      &eddy_diff,
                          const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>          &mix,
                          const AtmosphericTemperature<CoeffType,VectorCoeffType>                      &temp):
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

// in cm-3.km.s-1
     for(unsigned int s = 0; s < _mixture.neutral_composition().n_species(); s++)
     {
            omegas[s] = Antioch::constant_clone(T,1e-10) * (//omega = - ns * Dtilde * [
            - molecular[s] * 
            (
                dmolar_concentrations_dz[s] // 1/ns * dns_dz
              + molar_concentrations[s]/Hs[s]  // + 1/Hs
              + molar_concentrations[s] * dT_dz/T // + 1/T * dT_dz * (
                * (Antioch::constant_clone(T,1.) + ((nTot - molar_concentrations[s])/nTot) * _mixture.thermal_coefficient()[s]) //1 + (1 - xs)*alphas ) ]
            )
             - eddy_K * // - ns * K * (
            ( 
                dmolar_concentrations_dz[s] // 1/ns * dns_dz
              + molar_concentrations[s]/Ha // + 1/Ha
              + molar_concentrations[s] * dT_dz/T //+1/T * dT_dz )
            ));
     }
     return;
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template <typename StateType, typename VectorStateType, typename MatrixStateType>
  inline
  void DiffusionEvaluator<CoeffType, VectorCoeffType,MatrixCoeffType>::diffusion_and_derivs(const VectorStateType &molar_concentrations,
                                                                                            const VectorStateType &dmolar_concentrations_dz,
                                                                                            const StateType &z, 
                                                                                            VectorStateType &omegas_A_TERM,
                                                                                            VectorStateType &omegas_B_TERM,
                                                                                            MatrixStateType &domegas_dn_i_A_TERM,
                                                                                            MatrixStateType &domegas_dn_i_B_TERM) const
  {

     antioch_assert_equal_to(_mixture.neutral_composition().n_species(),molar_concentrations.size());
     antioch_assert_equal_to(_mixture.neutral_composition().n_species(),dmolar_concentrations_dz.size());
//params
     StateType nTot;
     Antioch::set_zero(nTot);
     for(unsigned int s = 0; s < molar_concentrations.size();s++)
     {
        nTot += molar_concentrations[s];
     }
     StateType T       = _temperature.neutral_temperature(z);
     StateType dT_dz   = _temperature.dneutral_temperature_dz(z);
     StateType dT_dz_T = dT_dz / T;

//eddy
     CoeffType dK_dn = _eddy_diffusion.K_deriv_ns(nTot);
     CoeffType K     = _eddy_diffusion.K(nTot);

//molecular
     VectorStateType Dtilde;
     Dtilde.resize(_mixture.neutral_composition().n_species(),0.);
     MatrixCoeffType dDtilde_dn;
     dDtilde_dn.resize(_mixture.neutral_composition().n_species());
     for(unsigned int s = 0; s < _mixture.neutral_composition().n_species(); s++)
     {
        dDtilde_dn[s].resize(_mixture.neutral_composition().n_species());
     }

     _molecular_diffusion.Dtilde_and_derivs_dn(molar_concentrations,T,nTot,Dtilde,dDtilde_dn);

//scale heights
     CoeffType Ha;
     VectorCoeffType dHa_dn_i;
     VectorCoeffType Hs;
     Hs.resize(_mixture.neutral_composition().n_species(),0.);
     dHa_dn_i.resize(_mixture.neutral_composition().n_species(),0.);

     _mixture.datmospheric_scale_height_dn_i(molar_concentrations,z,Ha,dHa_dn_i);
     _mixture.scale_heights(z,Hs);

     for(unsigned int s = 0; s < _mixture.neutral_composition().n_species(); s++)
     {
        omegas_A_TERM[s] = - Dtilde[s] - K ;
        omegas_B_TERM[s] = 
         - Dtilde[s]/Hs[s]  // + D/Hs
         - Dtilde[s] * dT_dz_T // + 1/T * dT_dz * (
                * (Antioch::constant_clone(T,1.) + ((nTot - molar_concentrations[s])/nTot) * _mixture.thermal_coefficient()[s]) //1 + (1 - xs)*alphas ) ]
        - K/Ha // + 1/Ha
        - K  * dT_dz_T; //+1/T * dT_dz )

       domegas_dn_i_A_TERM[s].clear();
       domegas_dn_i_A_TERM[s].resize(_mixture.neutral_composition().n_species());
       for(unsigned int i = 0; i < _mixture.neutral_composition().n_species(); i++)
       {
          domegas_dn_i_A_TERM[s][i] = - (dDtilde_dn[s][i] + dK_dn) * Antioch::constant_clone(T,1e-10); //to cm-3.km.s-1
          domegas_dn_i_B_TERM[s][i] = -  dDtilde_dn[s][i] * (   Antioch::constant_clone(Hs[s],1.)/Hs[s]
                                                              + dT_dz_T * (Antioch::constant_clone(T,1.) + 
                                                                             ( (nTot - molar_concentrations[s]) / nTot ) * _mixture.thermal_coefficient()[s]
                                                                          )
                                                            )
                                      - Dtilde[s] * dT_dz_T * _mixture.thermal_coefficient()[s] /nTot * (Antioch::constant_clone(nTot,1) - molar_concentrations[s] / nTot)
                                      - dK_dn * (Antioch::constant_clone(Ha,1) / Ha + dT_dz / T)
                                      - K * dHa_dn_i[i] / (Ha * Ha);
// in cm-3.km.s-1
         domegas_dn_i_A_TERM[s][i] *= Antioch::constant_clone(T,1e-10);
       }
// in cm-3.km.s-1
       omegas_A_TERM[s] *= Antioch::constant_clone(T,1e-10);
       omegas_B_TERM[s] *= Antioch::constant_clone(T,1e-10);
     }

     return;

  }


}

#endif
