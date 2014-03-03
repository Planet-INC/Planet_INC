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

#ifndef PLANET_DIRI_PDF_H
#define PLANET_DIRI_PDF_H

//Antioch
#include "antioch/antioch_asserts.h"

//Planet

//C++
#include <vector>

namespace Planet
{
  template <typename CoeffType>
  class DiriPdf: public  BasePdf<CoeffType>
  {
      private:
        std::vector<CoeffType> _v;
        CoeffType _r;

      public:
        DiriPdf();
        DiriPdf(const std::vector<CoeffType> &v, const CoeffType &r);
        ~DiriPdf();

        template <typename StateType>
        void set_parameters(const std::vector<StateType> &pars);

        template <typename StateType>
        void set_v(const std::vector<StateType> &v);

        template <typename StateType>
        void set_r(const StateType &r);

        const std::vector<CoeffType> v()   const;

        const CoeffType r()   const;

        const CoeffType value(unsigned int ip = 0) const;
  };

  template <typename CoeffType>
  inline
  DiriPdf<CoeffType>::DiriPdf():
      _r(0.L)
  {
     return;
  }

  template <typename CoeffType>
  inline
  DiriPdf<CoeffType>::DiriPdf(const std::vector<CoeffType> &v, const CoeffType &r):
      _v(v),
      _r(r)
  {
     return;
  }

  template <typename CoeffType>
  inline
  DiriPdf<CoeffType>::~DiriPdf()
  {
     return;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void DiriPdf<CoeffType>::set_v(const std::vector<StateType> &v)
  {
     _v.resize(v.size(),0.L);
     for(unsigned int i = 0; i < v.size(); i++)
     {
       _v[i] = v[i];
     }
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void DiriPdf<CoeffType>::set_r(const StateType &r)
  {
     _r = r;
  }

  template <typename CoeffType>
  inline
  const std::vector<CoeffType> DiriPdf<CoeffType>::v() const
  {
      return _v;
  }

  template <typename CoeffType>
  inline
  const CoeffType DiriPdf<CoeffType>::r() const
  {
      return _r;
  }

  template <typename CoeffType>
  inline
  const CoeffType DiriPdf<CoeffType>::value(unsigned int ip) const
  {
      return _v[ip];
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void DiriPdf<CoeffType>::set_parameters(const std::vector<StateType> &pars)
  {
     _v.resize(pars.size() - 1,0.L);
     for(unsigned int i = 0; i < pars.size() - 1; i++)
     {
       _v[i] = pars[i];
     }
      _r = pars.back();
  }
}

#endif
