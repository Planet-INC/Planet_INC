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

#ifndef PLANET_NORT_PDF_H
#define PLANET_NORT_PDF_H

//Antioch
#include "antioch/antioch_asserts.h"

//Planet

//C++

namespace Planet
{
  template <typename CoeffType>
  class NorTPdf: public  BasePdf<CoeffType>
  {
      private:
        CoeffType _mu;
        CoeffType _sigma;
        CoeffType _min;
        CoeffType _max;

      public:
        NorTPdf();
        NorTPdf(const CoeffType &mu, const CoeffType &sigma, const CoeffType &min, const CoeffType &max);
        ~NorTPdf();

        template <typename StateType>
        void set_parameters(const std::vector<StateType> &pars);

        template <typename StateType>
        void set_mu(const StateType &mu);

        template <typename StateType>
        void set_sigma(const StateType &sigma);

        template <typename StateType>
        void set_min(const StateType &min);

        template <typename StateType>
        void set_max(const StateType &max);

        const CoeffType mu()    const;

        const CoeffType sigma() const;

        const CoeffType min()   const;

        const CoeffType max()   const;

        const CoeffType value(unsigned int ip = 0) const;
  };

  template <typename CoeffType>
  inline
  NorTPdf<CoeffType>::NorTPdf():
      _mu(0.L),
      _sigma(0.L),
      _min(0.L),
      _max(0.L)
  {
     return;
  }

  template <typename CoeffType>
  inline
  NorTPdf<CoeffType>::NorTPdf(const CoeffType &mu, const CoeffType &sigma, const CoeffType &min, const CoeffType &max):
      _mu(mu),
      _sigma(sigma),
      _min(min),
      _max(max)
  {
     return;
  }

  template <typename CoeffType>
  inline
  NorTPdf<CoeffType>::~NorTPdf()
  {
     return;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void NorTPdf<CoeffType>::set_mu(const StateType &mu)
  {
     _mu = mu;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void NorTPdf<CoeffType>::set_sigma(const StateType &sigma)
  {
     _sigma = sigma;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void NorTPdf<CoeffType>::set_min(const StateType &min)
  {
     _min = min;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void NorTPdf<CoeffType>::set_max(const StateType &max)
  {
     _max = max;
  }

  template <typename CoeffType>
  inline
  const CoeffType NorTPdf<CoeffType>::mu() const
  {
      return _mu;
  }

  template <typename CoeffType>
  inline
  const CoeffType NorTPdf<CoeffType>::sigma() const
  {
      return _sigma;
  }

  template <typename CoeffType>
  inline
  const CoeffType NorTPdf<CoeffType>::min() const
  {
      return _min;
  }

  template <typename CoeffType>
  inline
  const CoeffType NorTPdf<CoeffType>::max() const
  {
      return _max;
  }

  template <typename CoeffType>
  inline
  const CoeffType NorTPdf<CoeffType>::value(unsigned int ip) const
  {
      return this->mu();
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void NorTPdf<CoeffType>::set_parameters(const std::vector<StateType> &pars)
  {
      antioch_assert_equal_to(pars.size(),4);
      _mu    = pars[0];
      _sigma = pars[1];
      _min   = pars[2];
      _max   = pars[3];
  }
}

#endif
