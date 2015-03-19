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

#ifndef PLANET_DIUN_PDF_H
#define PLANET_DIUN_PDF_H

//Antioch
#include "antioch/antioch_asserts.h"

//Planet
#include "planet/base_pdf.h"

//C++
#include <vector>

namespace Planet
{
  template <typename CoeffType>
  class DiUnPdf: public BasePdf<CoeffType>
  {
      private:
        unsigned int _n;

      public:
        DiUnPdf();
        DiUnPdf(unsigned int n);
        ~DiUnPdf();

        template <typename StateType>
        void set_parameters(const std::vector<StateType> &pars);

        template <typename StateType>
        void get_parameters(std::vector<StateType> &pars) const;

        void set_n(unsigned int n);

        unsigned int n() const;

        const CoeffType value(unsigned int ip = 0) const;

        void print(std::ostream &out = std::cout)  const;
  };

  template <typename CoeffType>
  inline
  DiUnPdf<CoeffType>::DiUnPdf():
      BasePdf<CoeffType>(PDFName::DiUn),
      _n(0)
  {
     return;
  }

  template <typename CoeffType>
  inline
  DiUnPdf<CoeffType>::DiUnPdf(unsigned int n):
      BasePdf<CoeffType>(PDFName::DiUn),
      _n(n)
  {
     return;
  }

  template <typename CoeffType>
  inline
  DiUnPdf<CoeffType>::~DiUnPdf()
  {
     return;
  }

  template <typename CoeffType>
  inline
  void DiUnPdf<CoeffType>::set_n(unsigned int n)
  {
     _n = n;
  }

  template <typename CoeffType>
  inline
  unsigned int DiUnPdf<CoeffType>::n() const
  {
      return _n;
  }

  template <typename CoeffType>
  inline
  const CoeffType DiUnPdf<CoeffType>::value(unsigned int /*ip*/) const
  {
      return CoeffType(1.L)/(CoeffType)_n;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void DiUnPdf<CoeffType>::set_parameters(const std::vector<StateType> &pars)
  {
      antioch_assert_equal_to(pars.size(),1);
      _n = (unsigned int)(pars[0]);
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void DiUnPdf<CoeffType>::get_parameters(std::vector<StateType> &pars) const
  {
     pars.resize(1,0.);
     pars[0] = (StateType)_n;
  }

  template <typename CoeffType>
  inline
  void DiUnPdf<CoeffType>::print(std::ostream &out)  const
  {
     out << "DiUn(" 
         << _n
         << ")";
  }
}

#endif
