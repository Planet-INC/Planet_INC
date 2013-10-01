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
#include "antioch/metaprogramming_decl.h"
#include "antioch/cmath_shims.h"
#include "antioch/chemical_mixture.h"

//Planet
#include "planet/binary_diffusion.h"

//C++


namespace Planet
{
  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  class MolecularDiffusionEvaluator
  {
     private:
        //!
        MolecularDiffusionEvaluator();

        MatrixCoeffType _Dtilde;
        std::vector<std::vector<BinaryDiffusion<CoeffType> > > _diffusion;
        Antioch::ChemicalMixture<CoeffType> &_composition;
        const unsigned int _n_medium;


     public:
        //!
        MolecularDiffusionEvaluator(Antioch::ChemicalMixture<CoeffType> &comp);
        //!
        ~MolecularDiffusionEvaluator();

        //! update Dtilde
        template<typename StateType, typename VectorStateType, typename MatrixStateType>
        void make_diffusion(const Atmosphere<StateType,VectorStateType,MatrixStateType> &atm);

        //!
        CoeffType Dtilde(unsigned int s, unsigned int iz) const;

        //!
        VectorCoeffType Dtilde_top_to_bottom(unsigned int ineu) const;

        //!
        void Dtilde_bottom_to_top(unsigned int ineu, VectorCoeffType &diff_bottom_top) const;

        //!
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        binary_coefficient(unsigned int i, unsigned int j, const StateType &T, const StateType &P) const
        ANTIOCH_AUTOFUNC(StateType,(_diffusion[i][j].diffusion_model() != DiffusionType::NoData)?
                                             this->binary_coefficient_known(i,j,T,P):
                                                (_composition.M(j) < _composition.M(i))?
                                                        this->binary_coefficient_unknown_ji(i,j,T,P):
                                                        this->binary_coefficient_unknown_ij(i,j,T,P))

        //!
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        binary_coefficient_known(unsigned int i, unsigned int j, const StateType &T, const StateType &P) const
        ANTIOCH_AUTOFUNC(StateType,_diffusion[i][j].binary_coefficient(T,P))

        //!
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        binary_coefficient_unknown_ji(unsigned int i, unsigned int j, const StateType &T, const StateType &P) const
        ANTIOCH_AUTOFUNC(StateType,_diffusion[i][i].binary_coefficient(T,P) * 
                                   Antioch::ant_sqrt( (_composition.M(j)/_composition.M(i) + StateType(1.L)) / StateType(2.L) )
                                   )
        //!
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        binary_coefficient_unknown_ij(unsigned int i, unsigned int j, const StateType &T, const StateType &P) const
        ANTIOCH_AUTOFUNC(StateType,_diffusion[i][i].binary_coefficient(T,P) * Antioch::ant_sqrt(_composition.M(j)/_composition.M(i)))

        //!
        void set_chemical_mixture(Antioch::ChemicalMixture<CoeffType> &comp);

        //!
        void set_binary_coefficient(unsigned int i, unsigned int j, const BinaryDiffusion<CoeffType> &bin_coef);

  };


  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  CoeffType MolecularDiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::Dtilde(unsigned int s, unsigned int iz) const
  {
      return _Dtilde.at(iz).at(s);
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  VectorCoeffType MolecularDiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::Dtilde_top_to_bottom(unsigned int ineu) const
  {
     return _Dtilde[ineu];
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  void MolecularDiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::Dtilde_bottom_to_top(unsigned int ineu, VectorCoeffType &diff_bottom_top) const
  {
     diff_bottom_top.resize(_Dtilde[ineu].size(),0.L);
     for(unsigned int iz = 0; iz < _Dtilde[ineu].size(); iz++)
     {
        diff_bottom_top[iz] = _Dtilde[ineu][_Dtilde[ineu].size() - 1 - iz];
     }

     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType, typename MatrixStateType>
  inline
  void MolecularDiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::make_diffusion(const Atmosphere<StateType,VectorStateType,MatrixStateType> &atm)
  {
    _Dtilde.resize(_composition.n_species()); //from top to bottom
    std::vector<CoeffType> meanM;

    for(unsigned int s = 0; s < _composition.n_species(); s++)
    {
      meanM.clear();
      meanM.resize(atm.n_altitudes(),0.L);
      for(StateType z = atm.min_alt(); z <= atm.min_alt(); z += atm.step_alt())
      {
        unsigned int iz = atm.altitude_map().at(z);
        for(unsigned int i = 0; i < _composition.n_species(); i++)
        {
          if(i == s)continue;
          meanM[iz] += _composition.M(i) * atm.neutral_molar_fraction(i,z);
        }
        meanM[iz] /= StateType(_composition.n_species() - 1);
      }
       _Dtilde[s].resize(atm.n_altitudes(),0.L);
       for(StateType z = atm.min_alt(); z <= atm.min_alt(); z += atm.step_alt())
       {
          unsigned int iz = atm.altitude_map().at(z);
          CoeffType n_D;
          Antioch::set_zero(n_D);
          for(unsigned int i = 0; i < _n_medium; i++)
          {
            if(i == s)continue;
            n_D += atm.neutral_molar_density(i,z) / this->binary_coefficient(i,s,atm.temperature(z),atm.pressure(z));
          }
          CoeffType Ds = (atm.total_density(z) - atm.neutral_molar_density(s,z))/n_D;
          _Dtilde[s][iz] = Ds / (CoeffType(1.L) - atm.neutral_molar_density(s,z) / atm.total_density(z) * 
                                  (CoeffType(1.L) - _composition.M(s) / meanM[iz])
                                );
       }
    }
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  MolecularDiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::MolecularDiffusionEvaluator()
  {
     antioch_error();
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  MolecularDiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::~MolecularDiffusionEvaluator()
  {
     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  MolecularDiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::MolecularDiffusionEvaluator(Antioch::ChemicalMixture<CoeffType> &comp):
       _composition(comp),
       _n_medium(2)
  {
     _diffusion.resize(_n_medium); //N2,CH4
     for(unsigned int i = 0; i < _diffusion.size(); i++)
     {
        _diffusion[i].resize(_composition.n_species());
     }
     return; 
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  void MolecularDiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::set_binary_coefficient(unsigned int i, unsigned int j, const BinaryDiffusion<CoeffType> &bin_coef)
  {
      _diffusion[i][j] = bin_coef;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  void MolecularDiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::set_chemical_mixture(Antioch::ChemicalMixture<CoeffType> &comp)
  {
     _composition = comp;
     return;
  }

}

#endif
