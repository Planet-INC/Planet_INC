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
                                   Antioch::ant_sqrt( (_mixture.neutral_composition().M(j)/_mixture.neutral_composition().M(_i_medium[m]) + Antioch::constant_clone(T,1)) / Antioch::constant_clone(T,2) )
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
        binary_coefficient_deriv_n(unsigned int m, unsigned int s, unsigned int /*i*/, const StateType & T, const StateType & P, const StateType & nTot) const
        ANTIOCH_AUTOFUNC(StateType,(_diffusion[m][s].diffusion_model() != DiffusionType::NoData)?
                                             _diffusion[m][s].binary_coefficient_deriv_n(T,P,nTot):
                                                (_mixture.neutral_composition().M(_i_medium[m]) < _mixture.neutral_composition().M(s))?
                                                    Antioch::ant_sqrt(_mixture.neutral_composition().M(s) / _mixture.neutral_composition().M(_i_medium[m])) * 
                                                                        _diffusion[m][_i_medium[m]].binary_coefficient_deriv_n(T,P,nTot)
                                                        :
                                                    Antioch::ant_sqrt( _mixture.neutral_composition().M(s) /( Antioch::constant_clone(T,2) * _mixture.neutral_composition().M(_i_medium[m])) + 
                                                                        Antioch::constant_clone(T,0.5) ) * 
                                                                        _diffusion[m][_i_medium[m]].binary_coefficient_deriv_n(T,P,nTot)
                        )

        //! \f$\frac{\partial D_{m,s}}{\partial T}\f$
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        binary_coefficient_deriv_T(unsigned int m, unsigned int s, const StateType & T, const StateType & P) const
        ANTIOCH_AUTOFUNC(StateType,(_diffusion[m][s].diffusion_model() != DiffusionType::NoData)?
                                             _diffusion[m][s].binary_coefficient_deriv_T(T,P):
                                                (_mixture.neutral_composition().M(_i_medium[m]) < _mixture.neutral_composition().M(s))?
                                                    Antioch::ant_sqrt(_mixture.neutral_composition().M(s) / _mixture.neutral_composition().M(_i_medium[m])) * 
                                                                        _diffusion[m][_i_medium[m]].binary_coefficient_deriv_T(T,P)
                                                    :
                                                    Antioch::ant_sqrt( _mixture.neutral_composition().M(s) /( Antioch::constant_clone(T,2) * _mixture.neutral_composition().M(_i_medium[m])) + 
                                                                        Antioch::constant_clone(T,0.5) ) * 
                                                                        _diffusion[m][_i_medium[m]].binary_coefficient_deriv_T(T,P)
                        )

        //! \f$\tilde{D}\f$ and all the derivatives with respect to concentrations
        template<typename StateType, typename VectorStateType, typename MatrixStateType>
        void Dtilde_and_derivs_dn(const VectorStateType &molar_concentrations, const StateType &T, const StateType &nTot, 
                                  VectorStateType &Dtilde, MatrixStateType &dD_dns) const;

        //! \f$\tilde{D}_s\f$ and all its derivatives with respect to concentrations
        template<typename StateType, typename VectorStateType>
        void Dtilde_and_derivative_n(unsigned int s, const VectorStateType &molar_concentrations, const StateType &T, const StateType &nTot,
                                     StateType & Dtilde, VectorStateType & dDtilde_dn) const;

        //! derivative with respect to \f$T\f$, kept separated at the moment for efficiency
        template<typename StateType, typename VectorStateType>
        void dDtilde_dT(const VectorStateType &molar_concentrations, const StateType &T, VectorStateType & dDtilde_dT) const;


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

     for(unsigned int s = 0; s < dD_dns.size(); s++)
     {
       this->Dtilde_and_derivative_n(s,molar_concentrations,T,nTot,Dtilde[s],dD_dns[s]);
     }
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void MolecularDiffusionEvaluator<CoeffType, VectorCoeffType,MatrixCoeffType>::Dtilde_and_derivative_n(unsigned int s,
                                                                                                        const VectorStateType &molar_concentrations,
                                                                                                        const StateType &T,const StateType &nTot,
                                                                                                        StateType & Dtilde,
                                                                                                        VectorStateType & dDtilde_dn) const
  {

        antioch_assert_equal_to(molar_concentrations.size(),_mixture.neutral_composition().n_species());
        antioch_assert_equal_to(dDtilde_dn.size(),_mixture.neutral_composition().n_species());

/////////////
// Wilke part
/////////////

///// temp
        // nt - ns
        StateType nTot_diff = nTot - molar_concentrations[s];
        // 1/(nt - ns)
        StateType one_over_nTot_diff = Antioch::constant_clone(T,1) / nTot_diff;

//// M_{\neq}
        StateType meanM;
        Antioch::set_zero(meanM);
        for(unsigned int j = 0; j < _mixture.neutral_composition().n_species(); j++)
        {
          if(j == s)continue;
          meanM += _mixture.neutral_composition().M(j) * molar_concentrations[j]; 
        }
        meanM /= nTot_diff;
// 1 / nt
        StateType one_over_nTot = Antioch::constant_clone(T,1) / nTot;
// Ms / M_{\neq}
        StateType Ms_over_Mmean = _mixture.neutral_composition().M(s) / meanM;

// ns / nt * Ms / M_{\neq}^2 / (ntot - ns)
        StateType molar_mass_term = molar_concentrations[s] * one_over_nTot * Ms_over_Mmean / meanM * one_over_nTot_diff;
// ns / nt^2
        StateType ns_over_ntot_square = molar_concentrations[s] * one_over_nTot * one_over_nTot;
// 1 - Ms / M_{\neq}
        StateType one_minus_Ms_over_Mdiff = Antioch::constant_clone(T,1) - Ms_over_Mmean;

// bimolecular diffusion, Wilke rule
        StateType sum_bimol;
        Antioch::set_zero(sum_bimol);
        StateType p = nTot * Antioch::constant_clone(nTot,1e6) //cm-3 -> m-3
                           * Constants::Universal::kb<CoeffType>() * T;

        for(unsigned int m = 0; m < _n_medium; m++)
        {
          if(_i_medium[m] == s)continue;
          sum_bimol += molar_concentrations[_i_medium[m]] / this->binary_coefficient(m,s,T,p);
        }
// Ds
        StateType one_over_Ds = sum_bimol / nTot_diff;

// term for derivative
        std::map<unsigned int, StateType> Ds_medium_term;
        for(unsigned int m = 0; m < _n_medium; m++)
        {
          if(_i_medium[m] == s)continue;
          Ds_medium_term[_i_medium[m]] = Antioch::constant_clone(T,1) / (this->binary_coefficient(m,s,T,p) * sum_bimol);          
        }

///////////
// Dtilde, De La Haye modification
//////
        
        StateType denom = Antioch::constant_clone(T,1) - molar_concentrations[s] * one_over_nTot * one_minus_Ms_over_Mdiff;

        Dtilde = nTot_diff / (denom * sum_bimol);
//      dDtilde_dT = dDs_dT * Dtilde / Ds

        for(unsigned int k = 0; k < _mixture.neutral_composition().n_species(); k++)
        {
          dDtilde_dn[k] = - ns_over_ntot_square;
          if(k == s)dDtilde_dn[k] += one_over_nTot;
          dDtilde_dn[k] *= one_minus_Ms_over_Mdiff;
          if(k != s)dDtilde_dn[k] += molar_mass_term * (_mixture.neutral_composition().M(k) - meanM);
          dDtilde_dn[k] *= Dtilde * one_over_Ds;
          dDtilde_dn[k] -= one_over_nTot;
          if(Ds_medium_term.count(k))dDtilde_dn[k] -= Ds_medium_term[k];
          if(k != s)dDtilde_dn[k] += one_over_nTot_diff;
          dDtilde_dn[k] *= Dtilde;

        }

        return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void MolecularDiffusionEvaluator<CoeffType, VectorCoeffType,MatrixCoeffType>::dDtilde_dT(const VectorStateType &molar_concentrations, 
                                                                                           const StateType &T, VectorStateType & dDtilde_dT) const
  {
     antioch_assert_equal_to(dDtilde_dT.size(),_mixture.neutral_composition().n_species());
     antioch_assert_equal_to(molar_concentrations.size(),_mixture.neutral_composition().n_species());


// global
     StateType nTot;
     Antioch::set_zero(nTot);
     for(unsigned int j = 0; j < _mixture.neutral_composition().n_species(); j++)
     {
        nTot += molar_concentrations[j]; 
     }

     StateType p = nTot * Antioch::constant_clone(nTot,1e6) //cm-3 -> m-3
                        * Constants::Universal::kb<CoeffType>() * T;


     for(unsigned int s = 0; s < molar_concentrations.size(); s++)
     {

// dDs_dT
// Ds denominator : sum_{j_m} n_{j_m}/D_{s,j_m}
// D_s = (nT - ns) / (sum_{medium} n_{medium}/D_{medium,s})
        StateType nTot_diff = nTot - molar_concentrations[s];
        StateType sum_bimol;
        Antioch::set_zero(sum_bimol);
        StateType numerator;
        Antioch::set_zero(numerator);
        for(unsigned int m = 0; m < _n_medium; m++)
        {
          if(_i_medium[m] == s)continue;

          StateType Dsjm = Antioch::constant_clone(T,1) / this->binary_coefficient(m,s,T,p);
          StateType dDsjm_dT = this->binary_coefficient_deriv_T(m,s,T, p);

          sum_bimol += molar_concentrations[_i_medium[m]] * Dsjm;
          numerator += molar_concentrations[_i_medium[m]] * Dsjm * Dsjm * dDsjm_dT;
        }
// Ds
        StateType dDs_dT = nTot_diff /(sum_bimol * sum_bimol) * numerator;

// DTilde, De La Haye version
        StateType meanM;
        Antioch::set_zero(meanM);
        for(unsigned int j = 0; j < _mixture.neutral_composition().n_species(); j++)
        {
          if(j == s)continue;
          meanM += _mixture.neutral_composition().M(j) * molar_concentrations[j]; 
        }
        meanM /= nTot_diff;

        dDtilde_dT[s] = dDs_dT / (Antioch::constant_clone(T,1) - molar_concentrations[s] / nTot * 
                                                                 (Antioch::constant_clone(T,1) - _mixture.neutral_composition().M(s) / meanM) );
     }

  }

}

#endif
