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

#ifndef PLANET_LOGN_PDF_H
#define PLANET_LOGN_PDF_H

//Antioch
#include "antioch/antioch_asserts.h"

//Planet
#include "planet/base_pdf.h"

//C++

namespace Planet
{
  template <typename CoeffType>
  class LogNPdf: public  BasePdf<CoeffType>
  {
      private:
        CoeffType _mu;
        CoeffType _f;

      public:
        LogNPdf();
        LogNPdf(const CoeffType &mu, const CoeffType &f);
        ~LogNPdf();

        template <typename StateType>
        void set_parameters(const std::vector<StateType> &pars);

        template <typename StateType>
        void set_mu(const StateType &mu);

        template <typename StateType>
        void set_f(const StateType &f);

        const CoeffType mu() const;

        const CoeffType f() const;

        const CoeffType value(unsigned int ip = 0) const;

        void print(std::ostream &out = std::cout)  const;
  };

  template <typename CoeffType>
  inline
  LogNPdf<CoeffType>::LogNPdf():
      BasePdf<CoeffType>(PDFName::LogN),
      _mu(0.L),
      _f(0.L)
  {
     return;
  }

  template <typename CoeffType>
  inline
  LogNPdf<CoeffType>::LogNPdf(const CoeffType &mu, const CoeffType &f):
      BasePdf<CoeffType>(PDFName::LogN),
      _mu(mu),
      _f(f)
  {
     return;
  }

  template <typename CoeffType>
  inline
  LogNPdf<CoeffType>::~LogNPdf()
  {
     return;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void LogNPdf<CoeffType>::set_mu(const StateType &mu)
  {
     _mu = mu;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void LogNPdf<CoeffType>::set_f(const StateType &f)
  {
     _f = f;
  }

  template <typename CoeffType>
  inline
  const CoeffType LogNPdf<CoeffType>::mu() const
  {
      return _mu;
  }

  template <typename CoeffType>
  inline
  const CoeffType LogNPdf<CoeffType>::f() const
  {
      return _f;
  }

  template <typename CoeffType>
  inline
  const CoeffType LogNPdf<CoeffType>::value(unsigned int ip) const
  {
      return this->mu();
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void LogNPdf<CoeffType>::set_parameters(const std::vector<StateType> &pars)
  {
      antioch_assert_equal_to(pars.size(),2);
      _mu = pars[0];
      _f  = pars[1];
  }

  template <typename CoeffType>
  inline
  void LogNPdf<CoeffType>::print(std::ostream &out)  const
  {
     out << "LogN(" << _mu  << "," << _f
                    << ")";
  }
}

#endif
