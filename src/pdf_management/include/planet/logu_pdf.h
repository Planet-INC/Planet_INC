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

#ifndef PLANET_LOGU_PDF_H
#define PLANET_LOGU_PDF_H

//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/cmath_shims.h"

//Planet
#include "planet/base_pdf.h"

//C++
#include <cmath>

namespace Planet
{
  template <typename CoeffType>
  class LogUPdf: public  BasePdf<CoeffType>
  {
      private:
        CoeffType _min;
        CoeffType _max;

      public:
        LogUPdf();
        LogUPdf(const CoeffType &min, const CoeffType &max);
        ~LogUPdf();

        template <typename StateType>
        void set_parameters(const std::vector<StateType> &pars);

        template <typename StateType>
        void get_parameters(std::vector<StateType> &pars) const;

        template <typename StateType>
        void set_min(const StateType &min);

        template <typename StateType>
        void set_max(const StateType &max);

        const CoeffType min()   const;

        const CoeffType max()   const;

        const CoeffType value(unsigned int ip = 0) const;

        void print(std::ostream &out = std::cout)  const;
  };

  template <typename CoeffType>
  inline
  LogUPdf<CoeffType>::LogUPdf():
      BasePdf<CoeffType>(PDFName::LogU),
      _min(0.L),
      _max(0.L)
  {
     return;
  }

  template <typename CoeffType>
  inline
  LogUPdf<CoeffType>::LogUPdf(const CoeffType &min, const CoeffType &max):
      BasePdf<CoeffType>(PDFName::LogU),
      _min(min),
      _max(max)
  {
     return;
  }

  template <typename CoeffType>
  inline
  LogUPdf<CoeffType>::~LogUPdf()
  {
     return;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void LogUPdf<CoeffType>::set_min(const StateType &min)
  {
     _min = min;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void LogUPdf<CoeffType>::set_max(const StateType &max)
  {
     _max = max;
  }

  template <typename CoeffType>
  inline
  const CoeffType LogUPdf<CoeffType>::min() const
  {
      return _min;
  }

  template <typename CoeffType>
  inline
  const CoeffType LogUPdf<CoeffType>::max() const
  {
      return _max;
  }

  template <typename CoeffType>
  inline
  const CoeffType LogUPdf<CoeffType>::value(unsigned int /* ip */) const
  {
      return Antioch::ant_pow(this->min() * this->max(),CoeffType(0.5));
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void LogUPdf<CoeffType>::set_parameters(const std::vector<StateType> &pars)
  {
      antioch_assert_equal_to(pars.size(),2);
      _min = pars[0];
      _max = pars[1];
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void LogUPdf<CoeffType>::get_parameters(std::vector<StateType> &pars) const
  {
     pars.resize(2,0.L);
     pars[0] = _min;
     pars[1] = _max;
  }

  template <typename CoeffType>
  inline
  void LogUPdf<CoeffType>::print(std::ostream &out)  const
  {
     out << "LogU(" 
         << _min  << "," << _max
         << ")";
  }
}

#endif
