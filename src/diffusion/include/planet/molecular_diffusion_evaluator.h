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
#include "planet/altitude.h"
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

        MatrixCoeffType _Dtilde;
        std::vector<std::vector<BinaryDiffusion<CoeffType> > > _diffusion;
        const unsigned int _n_medium;

//dependencies
        AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &_mixture;
        Altitude<CoeffType,VectorCoeffType> &_altitude;
        AtmosphericTemperature<CoeffType,VectorCoeffType> &_temperature;


     public:
        //!
        MolecularDiffusionEvaluator(AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &comp,
                                    Altitude<CoeffType,VectorCoeffType> &alt,
                                    AtmosphericTemperature<CoeffType,VectorCoeffType> &temp,
                                    const std::vector<std::vector<BinaryDiffusion<CoeffType> > > &diff);
        //!
        ~MolecularDiffusionEvaluator();

        //!
        template<typename StateType>
        void set_binary_coefficient(unsigned int i, unsigned int j, const BinaryDiffusion<StateType> &bin_coef);

        //! update Dtilde
        void make_molecular_diffusion();

        //! returns the Dtilde matrix (species,altitudes) bottom to top
        const MatrixCoeffType Dtilde() const;

        //!
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        binary_coefficient(unsigned int i, unsigned int j, const StateType &T, const StateType &P) const
        ANTIOCH_AUTOFUNC(StateType,(_diffusion[i][j].diffusion_model() != DiffusionType::NoData)?
                                             this->binary_coefficient_known(i,j,T,P):
                                                (_mixture.neutral_composition().M(j) < _mixture.neutral_composition().M(i))?
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
                                   Antioch::ant_sqrt( (_mixture.neutral_composition().M(j)/_mixture.neutral_composition().M(i) + StateType(1.L)) / StateType(2.L) )
                                   )
        //!
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        binary_coefficient_unknown_ij(unsigned int i, unsigned int j, const StateType &T, const StateType &P) const
        ANTIOCH_AUTOFUNC(StateType,_diffusion[i][i].binary_coefficient(T,P) * Antioch::ant_sqrt(_mixture.neutral_composition().M(j)/_mixture.neutral_composition().M(i)))

  };


  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const MatrixCoeffType MolecularDiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::Dtilde() const
  {
      return _Dtilde;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  void MolecularDiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::make_molecular_diffusion()
  {
    _Dtilde.resize(_mixture.neutral_composition().n_species()); //from bottom to top
    std::vector<CoeffType> meanM;

    for(unsigned int s = 0; s < _mixture.neutral_composition().n_species(); s++)
    {
      meanM.clear();
      meanM.resize(_altitude.altitudes().size(),0.L);
      for(unsigned int iz = 0; iz < _altitude.altitudes().size(); iz++)
      {
        for(unsigned int i = 0; i < _mixture.neutral_composition().n_species(); i++)
        {
          if(i == s)continue;
          meanM[iz] += _mixture.neutral_composition().M(i) * _mixture.neutral_molar_fraction[i][iz];
        }
        meanM[iz] /= CoeffType(_mixture.neutral_composition().n_species() - 1);
      }

       _Dtilde[s].resize(_altitude.altitudes().size(),0.L);
       for(unsigned int iz = 0; iz < _altitude.altitudes().size(); iz++)
       {
          CoeffType n_D;
          Antioch::set_zero(n_D);
          for(unsigned int i = 0; i < _n_medium; i++)
          {
            if(i == s)continue;
            CoeffType p = _mixture.total_density()[iz] * Constants::Universal::kb<CoeffType>() * _temperature.neutral_temperature()[iz];
            n_D += _mixture.total_density()[iz] * _mixture.neutral_molar_fraction[i][iz] / this->binary_coefficient(i,s,_temperature.neutral_temperature()[iz],p);
          }
          CoeffType Ds = (_mixture.total_density()[iz] - _mixture.neutral_molar_fraction[s][iz] * _mixture.total_density()[iz])/n_D;
          _Dtilde[s][iz] = Ds / (CoeffType(1.L) - _mixture.neutral_molar_fraction[s][iz] * 
                                (CoeffType(1.L) - _mixture.neutral_composition().M(s) / meanM[iz])
                                );
       }
    }
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  MolecularDiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::~MolecularDiffusionEvaluator()
  {
     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  MolecularDiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::MolecularDiffusionEvaluator
                       (AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> &comp,
                        Altitude<CoeffType,VectorCoeffType> &alt,
                        AtmosphericTemperature<CoeffType,VectorCoeffType> &temp,
                        const std::vector<std::vector<BinaryDiffusion<CoeffType> > > &diff
                       ):
       _temperature(temp),
       _diffusion(diff),
       _n_medium(2),
       _mixture(comp),
       _altitude(alt)
  {
     _diffusion.resize(_n_medium); //N2,CH4
     for(unsigned int i = 0; i < _diffusion.size(); i++)
     {
        _diffusion[i].resize(_mixture.neutral_composition().n_species());
     }
     return; 
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType>
  inline
  void MolecularDiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::set_binary_coefficient(unsigned int i, 
                                                                                                        unsigned int j, 
                                                                                                        const BinaryDiffusion<StateType> &bin_coef)
  {
     antioch_assert_less(i,_n_medium);
     antioch_assert_less(j,_mixture.neutral_composition().n_species());
      _diffusion[i][j] = bin_coef;
  }

}

#endif
