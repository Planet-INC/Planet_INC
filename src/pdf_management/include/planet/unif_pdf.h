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

#ifndef PLANET_UNIF_PDF_H
#define PLANET_UNIF_PDF_H

//Antioch
#include "antioch/antioch_asserts.h"

//Planet
#include "planet/base_pdf.h"

//C++

namespace Planet
{
  template <typename CoeffType>
  class UnifPdf: public BasePdf<CoeffType>
  {
      private:
        CoeffType _min;
        CoeffType _max;

      public:
        UnifPdf();
        UnifPdf(const CoeffType &min, const CoeffType &max);
        ~UnifPdf();

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
  UnifPdf<CoeffType>::UnifPdf():
      BasePdf<CoeffType>(PDFName::Unif),
      _min(0.L),
      _max(0.L)
  {
     return;
  }

  template <typename CoeffType>
  inline
  UnifPdf<CoeffType>::UnifPdf(const CoeffType &min, const CoeffType &max):
      BasePdf<CoeffType>(PDFName::Unif),
      _min(min),
      _max(max)
  {
     return;
  }

  template <typename CoeffType>
  inline
  UnifPdf<CoeffType>::~UnifPdf()
  {
     return;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void UnifPdf<CoeffType>::set_min(const StateType &min)
  {
     _min = min;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void UnifPdf<CoeffType>::set_max(const StateType &max)
  {
     _max = max;
  }

  template <typename CoeffType>
  inline
  const CoeffType UnifPdf<CoeffType>::min() const
  {
      return _min;
  }

  template <typename CoeffType>
  inline
  const CoeffType UnifPdf<CoeffType>::max() const
  {
      return _max;
  }

  template <typename CoeffType>
  inline
  const CoeffType UnifPdf<CoeffType>::value(unsigned int ip) const
  {
      return (this->min() + this->max())/2.L;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void UnifPdf<CoeffType>::set_parameters(const std::vector<StateType> &pars)
  {
      antioch_assert_equal_to(pars.size(),2);
      _min = pars[0];
      _max = pars[1];
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void UnifPdf<CoeffType>::get_parameters(std::vector<StateType> &pars) const
  {
     pars.resize(2,0.L);
     pars[0] = _min;
     pars[1] = _max;
  }

  template <typename CoeffType>
  inline
  void UnifPdf<CoeffType>::print(std::ostream &out)  const
  {
     out << "Unif(" 
         << _min  << "," << _max
         << ")";
  }
}

#endif
