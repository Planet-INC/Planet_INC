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

#ifndef PLANET_DIFFUSION_H
#define PLANET_DIFFUSION_H

//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/species_enum.h"
#include "antioch/cmath_shims.h"

//Planet
#include "planet/diffusion_enum.h"
#include "planet/planet_constants.h"

//C++

namespace Planet{

/*!
 * Massman model saved
 * Translation from Wilson and Wakeham model
 */
template <typename CoeffType = double>
class BinaryDiffusion{

      private:
        CoeffType _D01;
        CoeffType _beta;
        DiffusionType _diffusion_model;

        Antioch::Species _mol1;
        Antioch::Species _mol2;

      public:
        //!
        BinaryDiffusion();
        //!
        BinaryDiffusion(const BinaryDiffusion &rhs);
        //! no data
        BinaryDiffusion(const Antioch::Species &mol1, const Antioch::Species & mol2);
        //!
        BinaryDiffusion(const Antioch::Species &mol1, const Antioch::Species & mol2, const CoeffType &par1, const CoeffType &par2, const DiffusionType &model);
        //!
        ~BinaryDiffusion();

        //!
        void set_molecules(const Antioch::Species &mol1, const Antioch::Species &mol2);
        //!
        void set_diffusion_model(const DiffusionType &model);
        //!
        template<typename StateType>
        void set_parameters(const StateType &par1, const StateType &par2);
        //!
        template<typename StateType>
        void set_binary_diffusion(const Antioch::Species &mol1, const Antioch::Species & mol2, 
                                  const StateType &par1, const StateType &par2, const DiffusionType &model);

        //!
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        binary_coefficient(const StateType &T, const StateType &P) const
        ANTIOCH_AUTOFUNC(StateType,_D01 * Constants::Convention::P_normal<StateType>() / P * 
                                    Antioch::ant_pow((T/Constants::Convention::T_standard<StateType>()),_beta))

        //!
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        operator()(const StateType &T, const StateType &P) const
        ANTIOCH_AUTOFUNC(StateType,this->binary_coefficient(T,P))

        //!
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        binary_coefficient_deriv_T(const StateType &T, const StateType &P) const
        ANTIOCH_AUTOFUNC(StateType,this->binary_coefficient(T,P) / T * _beta)

        //!
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        binary_coefficient_deriv_n(const StateType &T, const StateType &P, const StateType &nTot) const
        ANTIOCH_AUTOFUNC(StateType, - this->binary_coefficient(T,P) / nTot)

        //!
        template<typename StateType>
        void binary_coefficient_and_derivatives(const StateType &T, const StateType &P, const StateType &nTot,
                                                StateType &Dij, StateType &Dij_dT, StateType &Dij_dP) const;
       
        //!
        CoeffType D01()  const;
        //!
        CoeffType beta() const;

        //!
        DiffusionType diffusion_model() const;
        //!
        Antioch::Species mol1() const {return _mol1;}
        //!
        Antioch::Species mol2() const {return _mol2;}

};

template<typename CoeffType>
inline
BinaryDiffusion<CoeffType>::BinaryDiffusion(const BinaryDiffusion &rhs)
{
  _D01 = rhs.D01();
  _beta = rhs.beta();
  _mol1 = rhs.mol1();
  _mol2 = rhs.mol2();
  _diffusion_model = rhs.diffusion_model();
  return;
}

template<typename CoeffType>
inline
CoeffType BinaryDiffusion<CoeffType>::D01() const
{
  return _D01;
}

template<typename CoeffType>
inline
CoeffType BinaryDiffusion<CoeffType>::beta() const
{
  return _beta;
}

template<typename CoeffType>
inline
DiffusionType BinaryDiffusion<CoeffType>::diffusion_model() const
{
  return _diffusion_model;
}

template<typename CoeffType>
inline
BinaryDiffusion<CoeffType>::BinaryDiffusion():
_D01(-1.L),
_beta(0.L),
_diffusion_model(DiffusionType::NoData)
{
  return;
}

template<typename CoeffType>
inline
BinaryDiffusion<CoeffType>::BinaryDiffusion(const Antioch::Species &mol1, 
                                            const Antioch::Species &mol2):
_D01(-1.L),
_beta(0.L),
_diffusion_model(DiffusionType::NoData),
_mol1(mol1),
_mol2(mol2)
{
  return;
}

template<typename CoeffType>
inline
BinaryDiffusion<CoeffType>::BinaryDiffusion(const Antioch::Species &mol1, 
                                            const Antioch::Species &mol2, 
                                            const CoeffType &par1, 
                                            const CoeffType &par2, 
                                            const DiffusionType &model):
_diffusion_model(model),
_mol1(mol1),
_mol2(mol2)
{
  set_parameters(par1,par2);
  return;
}

template<typename CoeffType>
inline
BinaryDiffusion<CoeffType>::~BinaryDiffusion()
{
  return;
}

template<typename CoeffType>
inline
void BinaryDiffusion<CoeffType>::set_molecules(const Antioch::Species &mol1, const Antioch::Species &mol2)
{
  _mol1 = mol1;
  _mol2 = mol2;
  return;
}


template<typename CoeffType>
inline
void BinaryDiffusion<CoeffType>::set_diffusion_model(const DiffusionType &model)
{
  _diffusion_model = model;
  return;
}

template<typename CoeffType>
template<typename StateType>
inline
void BinaryDiffusion<CoeffType>::set_parameters(const StateType &par1, const StateType &par2)
{

  switch(_diffusion_model)
  {
    case DiffusionType::Massman:
    {
      _D01  = par1;
      _beta = par2;
      return;
    }
    case DiffusionType::Wilson:
    {
      _D01 = par1 * Antioch::ant_pow(Constants::Convention::T_standard<CoeffType>(),par2 + CoeffType(1.L)) * 
                    Constants::Universal::kb<CoeffType>() /
                    Constants::Convention::P_normal<CoeffType>();
      _beta = par2 + CoeffType(1.L);
      return;
    } 
    case DiffusionType::Wakeham:
    {
       _D01 = par1 * Antioch::ant_pow(Constants::Convention::T_standard<CoeffType>(),par2);
       _beta = par2;
       return;
    }
    case DiffusionType::NoData:
    {
       _D01 = -1.L;
       _beta = -1.L;
       return;
    }
    default:
    {
       antioch_not_implemented();
       return;
    }
  }

  return;
}

template<typename CoeffType>
template<typename StateType>
inline
void BinaryDiffusion<CoeffType>::set_binary_diffusion(const Antioch::Species &mol1, const Antioch::Species & mol2, 
                                                      const StateType &par1, const StateType &par2, const DiffusionType &model)
{
  this->set_diffusion_model(model);
  this->set_molecules(mol1,mol2);
  this->set_parameters(par1,par2);
}

template<typename CoeffType>
template<typename StateType>
inline
void BinaryDiffusion<CoeffType>::binary_coefficient_and_derivatives(const StateType &T, const StateType &P, const StateType &nTot,
                                                StateType &Dij, StateType &Dij_dT, StateType &Dij_dns) const
{
  Dij = this->binary_coefficient(T,P);
  Dij_dT  = this->binary_coefficient_deriv_T(T,P);
  Dij_dns = this->binary_coefficient_deriv_P(T,P,nTot);
  return;
}

} //namespace Planet

#endif
