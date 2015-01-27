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
        const AtmosphericMixture<CoeffType, VectorCoeffType,MatrixCoeffType>  & _mixture;
        const AtmosphericTemperature<CoeffType, VectorCoeffType>              & _temperature;

        //! The coefficients are unknown, i heavier than j
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        binary_coefficient_unknown_ji(unsigned int m, unsigned int j, const StateType &T, const StateType &P) const
        ANTIOCH_AUTOFUNC(StateType,_diffusion[m][m].binary_coefficient(T,P) * 
                                   Antioch::ant_sqrt( (_mixture.neutral_composition().M(j)/_mixture.neutral_composition().M(_i_medium[m]) + StateType(1.L)) / StateType(2.L) )
                                   )
        //! The coefficients are unknown, j heavier than i
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        binary_coefficient_unknown_ij(unsigned int m, unsigned int j, const StateType &T, const StateType &P) const
        ANTIOCH_AUTOFUNC(StateType,_diffusion[m][m].binary_coefficient(T,P) * Antioch::ant_sqrt(_mixture.neutral_composition().M(j)/_mixture.neutral_composition().M(_i_medium[m])))

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
        ANTIOCH_AUTO(StateType) 
        Dtilde(unsigned int s, const StateType & nTot, const StateType & T, const StateType & p, const VectorStateType & molar_concentrations) const;

        //!
        template<typename StateType, typename VectorStateType>
        void Dtilde(const VectorStateType &molar_concentrations, const StateType &T,
                    VectorStateType &Dtilde) const;// Dtilde

        //!
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        binary_coefficient(unsigned int m, unsigned int j, const StateType &T, const StateType &P) const
        ANTIOCH_AUTOFUNC(StateType,(_diffusion[m][j].diffusion_model() != DiffusionType::NoData)?
                                             _diffusion[m][j].binary_coefficient(T,P):
                                                (_mixture.neutral_composition().M(j) < _mixture.neutral_composition().M(_i_medium[m]))?
                                                        this->binary_coefficient_unknown_ji(m,j,T,P):
                                                        this->binary_coefficient_unknown_ij(m,j,T,P))

        //! \f$\frac{\partial D_{m,s}}{\partial n_i}\f$
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        binary_coefficient_deriv_n(unsigned int m, unsigned int s, unsigned int i, const StateType & T, const StateType & P, const StateType & nTot, const StateType & ni) const
        ANTIOCH_AUTOFUNC(StateType,(_diffusion[m][s].diffusion_model() != DiffusionType::NoData)?
                                             _diffusion[m][s].binary_coefficient_deriv_n(T,P,nTot,ni):
                                                (_mixture.neutral_composition().M(_i_medium[m]) < _mixture.neutral_composition().M(s))?
                                                    Antioch::ant_sqrt((_mixture.neutral_composition().M(_i_medium[m]) / _mixture.neutral_composition().M(s) + 
                                                                        Antioch::constant_clone(T,1))/Antioch::constant_clone(T,2)) * 
                                                                        _diffusion[m][m].binary_coefficient_deriv_n(T,P,nTot,ni):
                                                    Antioch::ant_sqrt(_mixture.neutral_composition().M(_i_medium[m]) / _mixture.neutral_composition().M(s)) * 
                                                                        _diffusion[m][m].binary_coefficient_deriv_n(T,P,nTot,ni)
                        )

        template<typename StateType, typename VectorStateType, typename MatrixStateType>
        void Dtilde_and_derivs_dn(const VectorStateType &molar_concentrations, const StateType &T, const StateType &nTot, 
                                  VectorStateType &Dtilde, MatrixStateType &dD_dns) const;

        template<typename StateType, typename VectorStateType>
        void Dtilde_s_dn_i(unsigned int s, unsigned int i, const VectorStateType &molar_concentrations, const StateType &T, const StateType &nTot,
                           StateType &dDtilde_s_dn_i) const;

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
      _i_medium[i] = _mixture.neutral_composition().species_name_map().at(medium_species[i]);
    }
  
    return; 
   }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void MolecularDiffusionEvaluator<CoeffType, VectorCoeffType,MatrixCoeffType>::Dtilde(const VectorStateType &molar_concentrations, 
                                                                       const StateType &T,VectorStateType &Dtilde) const
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
     CoeffType p = nTot * Antioch::constant_clone(nTot,1e6) //cm-3 -> m-3
                   * Constants::Universal::kb<CoeffType>() * T;
     for(unsigned int s = 0; s < _mixture.neutral_composition().n_species(); s++)
     {
        Dtilde[s] = this->Dtilde(s,nTot,T,p,molar_concentrations);
     }
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  ANTIOCH_AUTO(StateType) 
    MolecularDiffusionEvaluator<CoeffType, VectorCoeffType,MatrixCoeffType>::Dtilde(unsigned int s, const StateType & nTot, const StateType & T,
                                                                                    const StateType & p, const VectorStateType & molar_concentrations) const
  {
//M_{/=}
        StateType meanM;
        Antioch::set_zero(meanM);
        StateType ntot_s = nTot - molar_concentrations[s]; //ntot - ns
        for(unsigned int i = 0; i < _mixture.neutral_composition().n_species(); i++)
        {
          if(i == s)continue;
          meanM += _mixture.neutral_composition().M(i) * molar_concentrations[i]; //x_i without s: ni/(ntot - ns)
        }
        meanM /= ntot_s;
//Ds denominator : sum_{j_m} n_{j_m}/D_{s,j_m}
        StateType n_D;
        Antioch::set_zero(n_D);
        for(unsigned int m = 0; m < _n_medium; m++)
        {
          if(_i_medium[m] == s)continue;
          n_D += molar_concentrations[_i_medium[m]] / this->binary_coefficient(m,s,T,p);
        }
//Dtilde = Ds numerator (ntot - n_s) / Ds denom ... 
// cm2.s-1
        return (nTot - molar_concentrations[s])
                            / ( n_D * ( CoeffType(1.L) - molar_concentrations[s]/nTot * 
                                         (CoeffType(1.L) - _mixture.neutral_composition().M(s) / meanM)
                                      )
                              );
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

     this->Dtilde(molar_concentrations,T,Dtilde);
     for(unsigned int s = 0; s < dD_dns.size(); s++)
     {
       for(unsigned int i = 0; i < dD_dns.size(); i++)
       {
          this->Dtilde_s_dn_i(s,i,molar_concentrations,T,nTot,dD_dns[s][i]);
       }
     }
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void MolecularDiffusionEvaluator<CoeffType, VectorCoeffType,MatrixCoeffType>::Dtilde_s_dn_i(unsigned int s,unsigned int i,
                                                                                              const VectorStateType &molar_concentrations,
                                                                                              const StateType &T,const StateType &nTot,
                                                                                              StateType &dDtilde_s_dni) const
  {

        antioch_assert_equal_to(molar_concentrations.size(),_mixture.neutral_composition().n_species());

/////////////
// Wilke part
/////////////


        //D_s = (nT - ns) / (sum_{medium} n_{medium}/D_{medium,s})
        StateType Ds = nTot - molar_concentrations[s];
        StateType dDs_dni = (s == i)?Antioch::zero_clone(T):Antioch::constant_clone(T,1);
//Ds denominator : sum_{j_m} n_{j_m}/D_{s,j_m}
        StateType n_D;
        StateType dn_D_dni;
        Antioch::set_zero(n_D);
        Antioch::set_zero(dn_D_dni);
        StateType p = nTot * Antioch::constant_clone(nTot,1e6) //cm-3 -> m-3
                           * Constants::Universal::kb<CoeffType>() * T;

        for(unsigned int m = 0; m < _n_medium; m++)
        {
          if(_i_medium[m] == s)continue;

          StateType Dsjm = this->binary_coefficient(m,s,T,p);
          StateType dDsjm_dni = this->binary_coefficient_deriv_n(m,s,i,T, p, nTot, molar_concentrations[i]);

          n_D += molar_concentrations[_i_medium[m]] / Dsjm;
          dn_D_dni -= molar_concentrations[_i_medium[m]] / (Dsjm * Dsjm) * dDsjm_dni;
          if(_i_medium[m] == i)dn_D_dni += StateType(1.L) / Dsjm;
        }
        Ds /= n_D;
        dDs_dni /= n_D; // D(nT - ns)/Dni part

        dDs_dni -= (nTot - molar_concentrations[i]) * dn_D_dni / (n_D * n_D); // sum part

///////////
// Dtilde, De La Haye modification
//////

        StateType meanM;
        StateType dmeanM_dni = (i == s)?Antioch::zero_clone(meanM):_mixture.neutral_composition().M(i);
        Antioch::set_zero(meanM);
        Antioch::set_zero(dmeanM_dni);
        for(unsigned int j = 0; j < _mixture.neutral_composition().n_species(); j++)
        {
          if(j == s)continue;
          meanM += _mixture.neutral_composition().M(j) * molar_concentrations[j]; 
        }
        StateType ntot_s = nTot - molar_concentrations[s]; //ntot - ns
        StateType dntot_s_dni = (s == i)?Antioch::zero_clone(ntot_s):Antioch::constant_clone(ntot_s,1); //ntot - ns
        dmeanM_dni /= ntot_s;
        dmeanM_dni -= meanM * dntot_s_dni / (ntot_s * ntot_s);

        meanM /=  ntot_s; //x_i without s: ni/(ntot - ns)

        StateType denom = (Antioch::constant_clone(nTot,1) - molar_concentrations[s]/nTot * 
                        (Antioch::constant_clone(nTot,1) - _mixture.neutral_composition().M(s)/meanM));
        StateType ddenom_dni = (i == s)?StateType(1.)/nTot:Antioch::zero_clone(nTot);
        ddenom_dni -= molar_concentrations[s] / (nTot / nTot);
        ddenom_dni *= (StateType(1.) - _mixture.neutral_composition().M(s)/meanM);
        ddenom_dni += molar_concentrations[s]/nTot * _mixture.neutral_composition().M(s)/(meanM * meanM) * dmeanM_dni;

        dDtilde_s_dni = dDs_dni / denom - Ds * ddenom_dni / (denom * denom);

        return;
  }

}

#endif
