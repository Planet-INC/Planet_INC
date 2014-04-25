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

#ifndef _PLANET_MOLECULAR_DIFFUSION_EVALUATOR_
#define _PLANET_MOLECULAR_DIFFUSION_EVALUATOR_

//Antioch
#include "antioch/cmath_shims.h"
#include "antioch/chemical_mixture.h"

//Planet
#include "planet/binary_diffusion.h"
#include "planet/atmospheric_mixture.h"
#include "planet/atmospheric_temperature.h"

//C++


namespace Planet
{
  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  class MolecularDiffusionEvaluator
  {
     private:
        //!
        MolecularDiffusionEvaluator() {antioch_error();return;}

        unsigned int _n_medium;
        std::vector<unsigned int> _i_medium;

    //dependencies
        //! stock of binary diffusion coefficients
        const std::vector<std::vector<BinaryDiffusion<CoeffType> > >          & _diffusion;
        const AtmosphericMixture<CoeffType, VectorCoeffType,MatrixCoeffType>  &_mixture;
        const AtmosphericTemperature<CoeffType, VectorCoeffType>              &_temperature;

        //! The coefficients are known
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        binary_coefficient_known(unsigned int i, unsigned int j, const StateType &T, const StateType &P) const
        ANTIOCH_AUTOFUNC(StateType,_diffusion[i][j].binary_coefficient(T,P))

        //! The coefficients are unknown, i heavier than j
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        binary_coefficient_unknown_ji(unsigned int i, unsigned int j, const StateType &T, const StateType &P) const
        ANTIOCH_AUTOFUNC(StateType,_diffusion[i][i].binary_coefficient(T,P) * 
                                   Antioch::ant_sqrt( (_mixture.neutral_composition().M(j)/_mixture.neutral_composition().M(i) + StateType(1.L)) / StateType(2.L) )
                                   )
        //! The coefficients are unknown, j heavier than i
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        binary_coefficient_unknown_ij(unsigned int i, unsigned int j, const StateType &T, const StateType &P) const
        ANTIOCH_AUTOFUNC(StateType,_diffusion[i][i].binary_coefficient(T,P) * Antioch::ant_sqrt(_mixture.neutral_composition().M(j)/_mixture.neutral_composition().M(i)))

        void set_medium_species(const std::vector<std::string> &medium_species);


     public:
        //!
        MolecularDiffusionEvaluator( const std::vector<std::vector<BinaryDiffusion<CoeffType> > > &diff,
                                     const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &comp,
                                     const AtmosphericTemperature<CoeffType,VectorCoeffType> &temp,
                                     const std::vector<std::string> & medium);
        //!
        ~MolecularDiffusionEvaluator();

        //!
        template<typename StateType, typename VectorStateType>
        void Dtilde(const VectorStateType &molar_concentrations, const StateType &z,
                    VectorStateType &Dtilde) const;// Dtilde

        //!
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        binary_coefficient(unsigned int i, unsigned int j, const StateType &T, const StateType &P) const
        ANTIOCH_AUTOFUNC(StateType,(_diffusion[i][j].diffusion_model() != DiffusionType::NoData)?
                                             this->binary_coefficient_known(i,j,T,P):
                                                (_mixture.neutral_composition().M(j) < _mixture.neutral_composition().M(i))?
                                                        this->binary_coefficient_unknown_ji(i,j,T,P):
                                                        this->binary_coefficient_unknown_ij(i,j,T,P))

        template<typename StateType, typename VectorStateType, typename MatrixStateType>
        void Dtilde_and_derivs_dn(const VectorStateType &molar_concentrations, const StateType &T, const StateType &nTot, 
                                  VectorStateType &Dtilde, MatrixStateType &dD_dns) const;

        template<typename StateType, typename VectorStateType>
        void Dtilde_s_dn_i(unsigned int s, const VectorStateType &molar_concentrations, const StateType &T,const StateType &nTot,
                                                                                        StateType & Dtilde_s, VectorStateType &dDtilde_s_dn_i) const;

  };

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  MolecularDiffusionEvaluator<CoeffType, VectorCoeffType,MatrixCoeffType>::~MolecularDiffusionEvaluator()
  {
     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  MolecularDiffusionEvaluator<CoeffType, VectorCoeffType,MatrixCoeffType>::MolecularDiffusionEvaluator
                       (const std::vector<std::vector<BinaryDiffusion<CoeffType> > > &diff,
                        const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &comp,
                        const AtmosphericTemperature<CoeffType,VectorCoeffType> &temp,
                        const std::vector<std::string> & medium
                       ):
       _n_medium(diff.size()),
       _diffusion(diff),
       _mixture(comp),
       _temperature(temp)
  {
     this->set_medium_species(medium);
     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  void MolecularDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::set_medium_species(const std::vector<std::string> &medium_species)
  {
    antioch_assert_equal_to(_n_medium,medium_species.size());
    _i_medium.resize(_n_medium);
    for(unsigned int i = 0; i < _n_medium; i++)
    {
      _i_medium[i] = _mixture.neutral_composition().active_species_name_map().at(medium_species[i]);
    }
  
    return; 
   }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void MolecularDiffusionEvaluator<CoeffType, VectorCoeffType,MatrixCoeffType>::Dtilde(const VectorStateType &molar_concentrations, 
                                                                       const StateType &z,VectorStateType &Dtilde) const
  {
     antioch_assert_equal_to(molar_concentrations.size(),_mixture.neutral_composition().n_species());

     Dtilde.resize(molar_concentrations.size(),0.L);
     CoeffType nTot;
     Antioch::set_zero(nTot);
     for(unsigned int s = 0; s < molar_concentrations.size(); s++)
     {
        nTot += molar_concentrations[s];
     }
// p = n * kb * T  (Pa)
     CoeffType T = _temperature.neutral_temperature(z);
     CoeffType p = nTot * Antioch::constant_clone(nTot,1e6) //cm-3 -> m-3
                   * Constants::Universal::kb<CoeffType>() * T;
     for(unsigned int s = 0; s < _mixture.neutral_composition().n_species(); s++)
     {
//M_{/=}
        CoeffType meanM;
        Antioch::set_zero(meanM);
        CoeffType ntot_s = nTot - molar_concentrations[s]; //ntot - ns
        for(unsigned int i = 0; i < _mixture.neutral_composition().n_species(); i++)
        {
          if(i == s)continue;
          meanM += _mixture.neutral_composition().M(i) * molar_concentrations[i] / ntot_s; //x_i without s: ni/(ntot - ns)
        }
//Ds denominator : sum_{j_m} n_{j_m}/D_{s,j_m}
        CoeffType n_D;
        Antioch::set_zero(n_D);
        for(unsigned int i = 0; i < _n_medium; i++)
        {
          if(_i_medium[i] == s)continue;
          n_D += molar_concentrations[_i_medium[i]] / this->binary_coefficient(i,s,T,p);
        }
//Dtilde = Ds numerator (ntot - n_s) / Ds denom ... 
// cm2.s-1
        Dtilde[s] = (nTot - molar_concentrations[s])
                            / ( n_D * (CoeffType(1.L) - molar_concentrations[s]/nTot * 
                                      (CoeffType(1.L) - _mixture.neutral_composition().M(s) / meanM))
                              );
       }
       return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType, typename MatrixStateType>
  inline
  void MolecularDiffusionEvaluator<CoeffType, VectorCoeffType,MatrixCoeffType>::Dtilde_and_derivs_dn(const VectorStateType &molar_concentrations, const StateType &T, const StateType &nTot, 
                                                                                                VectorStateType &Dtilde, MatrixStateType &dD_dns) const
  {
     antioch_assert_equal_to(molar_concentrations.size(),_mixture.neutral_composition().n_species());
     antioch_assert_equal_to(Dtilde.size(),_mixture.neutral_composition().n_species());
     antioch_assert_equal_to(dD_dns.size(),_mixture.neutral_composition().n_species());
#ifdef NDEBUG
#else
     for(unsigned int s = 0; s < dD_dns.size(); s++)
     {
        antioch_assert_equal_to(dD_dns[s].size(),_mixture.neutral_composition().n_species());
     }
#endif

     for(unsigned int s = 0; s < dD_dns.size(); s++)
     {
        this->Dtilde_s_dn_i(s,molar_concentrations,T,nTot,Dtilde[s],dD_dns[s]);
     }
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void MolecularDiffusionEvaluator<CoeffType, VectorCoeffType,MatrixCoeffType>::Dtilde_s_dn_i(unsigned int s,
                                                                                                   const VectorStateType &molar_concentrations,
                                                                                                   const StateType &T,const StateType &nTot,
                                                                                                   StateType & Dtilde_s, VectorStateType &dDtilde_s_dn) const
  {

        antioch_assert_equal_to(molar_concentrations.size(),_mixture.neutral_composition().n_species());
        antioch_assert_equal_to(dDtilde_s_dn.size(),_mixture.neutral_composition().n_species());

        //D_s = (nT - ns) / (sum_{medium} n_{medium}/D_{medium,s})
        CoeffType Ds = (nTot - molar_concentrations[s]); //cm-3
//Ds denominator : sum_{j_m} n_{j_m}/D_{s,j_m}
        CoeffType n_D;
        Antioch::set_zero(n_D);
        CoeffType p = nTot * Antioch::constant_clone(nTot,1e6) //cm-3 -> m-3
                           * Constants::Universal::kb<CoeffType>() * T;

        for(unsigned int m = 0; m < _n_medium; m++)
        {
          if(_i_medium[m] == s)continue;
          n_D += molar_concentrations[_i_medium[m]] / this->binary_coefficient(m,s,T,p);
        }
        Ds /= n_D;

        CoeffType meanM;
        Antioch::set_zero(meanM);
        CoeffType ntot_s = nTot - molar_concentrations[s]; //ntot - ns
        for(unsigned int j = 0; j < _mixture.neutral_composition().n_species(); j++)
        {
          if(j == s)continue;
          meanM += _mixture.neutral_composition().M(j) * molar_concentrations[j] / ntot_s; //x_i without s: ni/(ntot - ns)
        }

        CoeffType sum = (Antioch::constant_clone(nTot,1) - molar_concentrations[s]/nTot * 
                        (Antioch::constant_clone(nTot,1) - _mixture.neutral_composition().M(s)/meanM));

        Dtilde_s = Ds / sum;

        for(unsigned int i = 0; i < molar_concentrations.size(); i++)
        {

          CoeffType first = (i == s)?Antioch::constant_clone(nTot,1.):Antioch::zero_clone(nTot);

          CoeffType second;
          Antioch::set_zero(second);
          for(unsigned int m = 0; m < _n_medium; m++)
          {
            if(_i_medium[m] == s)continue;
            if(_i_medium[m] == i)second = Antioch::constant_clone(nTot,1)/ this->binary_coefficient(m,s,T,p);
          }

          CoeffType dDs_dn_i = Ds * ( first - (second + n_D/nTot)/n_D);

          CoeffType third = _mixture.neutral_composition().M(i) * molar_concentrations[i];
          if(i == s)Antioch::set_zero(third);
          CoeffType dmeanM_dn_i = (third - meanM * first) / (nTot - molar_concentrations[s]);

          CoeffType dsum_dn_i = molar_concentrations[s]/(nTot * nTot) + _mixture.neutral_composition().M(s)/(meanM * meanM) * dmeanM_dn_i;
          if(i == s)dsum_dn_i -= Antioch::constant_clone(nTot,1) / nTot;


          dDtilde_s_dn[i] = Dtilde_s * (dDs_dn_i/Ds - dsum_dn_i/sum);
        }

        return;
  }

}

#endif
