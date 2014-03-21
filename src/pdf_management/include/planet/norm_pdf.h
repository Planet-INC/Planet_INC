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

#ifndef PLANET_NORM_PDF_H
#define PLANET_NORM_PDF_H

//Antioch
#include "antioch/antioch_asserts.h"

//Planet
#include "planet/base_pdf.h"

//C++

namespace Planet
{
  template <typename CoeffType>
  class NormPdf: public BasePdf<CoeffType>
  {
      private:
        CoeffType _mu;
        CoeffType _sigma;

      public:
        NormPdf();
        NormPdf(const CoeffType &mu, const CoeffType &sigma);
        ~NormPdf();

        template <typename StateType>
        void set_parameters(const std::vector<StateType> &pars);

        template <typename StateType>
        void get_parameters(std::vector<StateType> &pars) const;

        template <typename StateType>
        void set_mu(const StateType &mu);

        template <typename StateType>
        void set_sigma(const StateType &sigma);

        const CoeffType mu() const;

        const CoeffType sigma() const;

        const CoeffType value(unsigned int ip = 0) const;

        void print(std::ostream &out = std::cout)  const;

  };

  template <typename CoeffType>
  inline
  NormPdf<CoeffType>::NormPdf():
      BasePdf<CoeffType>(PDFName::Norm),
      _mu(0.L),
      _sigma(0.L)
  {
     return;
  }

  template <typename CoeffType>
  inline
  NormPdf<CoeffType>::NormPdf(const CoeffType &mu, const CoeffType &sigma):
      BasePdf<CoeffType>(PDFName::Norm),
      _mu(mu),
      _sigma(sigma)
  {
     return;
  }

  template <typename CoeffType>
  inline
  NormPdf<CoeffType>::~NormPdf()
  {
     return;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void NormPdf<CoeffType>::set_mu(const StateType &mu)
  {
     _mu = mu;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void NormPdf<CoeffType>::set_sigma(const StateType &sigma)
  {
     _sigma = sigma;
  }

  template <typename CoeffType>
  inline
  const CoeffType NormPdf<CoeffType>::mu() const
  {
      return _mu;
  }

  template <typename CoeffType>
  inline
  const CoeffType NormPdf<CoeffType>::sigma() const
  {
      return _sigma;
  }

  template <typename CoeffType>
  inline
  const CoeffType NormPdf<CoeffType>::value(unsigned int ip) const
  {
      return this->mu();
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void NormPdf<CoeffType>::set_parameters(const std::vector<StateType> &pars)
  {
      antioch_assert_equal_to(pars.size(),2);
      _mu    = pars[0];
      _sigma = pars[1];
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void NormPdf<CoeffType>::get_parameters(std::vector<StateType> &pars) const
  {
     pars.resize(2,0.L);
     pars[0] = _mu;
     pars[1] = _sigma;
  }

  template <typename CoeffType>
  inline
  void NormPdf<CoeffType>::print(std::ostream &out)  const
  {
     out << "Norm(" << _mu << "," << _sigma << ")";
  }

}//end namespace

#endif
