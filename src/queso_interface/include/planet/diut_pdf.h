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

#ifndef PLANET_DIUT_PDF_H
#define PLANET_DIUT_PDF_H

//Antioch
#include "antioch/antioch_asserts.h"

//Planet
#include "planet/base_pdf.h"

//C++
#include <vector>

namespace Planet
{
  template <typename CoeffType>
  class DiUTPdf: public BasePdf<CoeffType>
  {
      private:
        std::vector<CoeffType> _min;
        std::vector<CoeffType> _max;

      public:
        DiUTPdf();
        DiUTPdf(const std::vector<CoeffType> &min, const std::vector<CoeffType> &max);
        ~DiUTPdf();

        template <typename StateType>
        void set_parameters(const std::vector<StateType> &pars);

        template <typename StateType>
        void set_min(const std::vector<StateType> &min);

        template <typename StateType>
        void set_max(const std::vector<StateType> &max);

        const std::vector<CoeffType> min()  const;

        const std::vector<CoeffType> max() const;

        const CoeffType value(unsigned int ip = 0) const;

        void print(std::ostream &out = std::cout)  const;
  };

  template <typename CoeffType>
  inline
  DiUTPdf<CoeffType>::DiUTPdf():
      BasePdf<CoeffType>(PDFName::DiUT)
  {
     return;
  }

  template <typename CoeffType>
  inline
  DiUTPdf<CoeffType>::DiUTPdf(const std::vector<CoeffType> &min, const std::vector<CoeffType> &max):
      BasePdf<CoeffType>(PDFName::DiUT),
      _min(min),
      _max(max)
  {
     return;
  }

  template <typename CoeffType>
  inline
  DiUTPdf<CoeffType>::~DiUTPdf()
  {
     return;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void DiUTPdf<CoeffType>::set_min(const std::vector<StateType> &min)
  {
     _min.resize(min.size(),0.L);
     for(unsigned int i = 0; i < min.size(); i++)
     {
       _min[i] = min[i];
     }
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void DiUTPdf<CoeffType>::set_max(const std::vector<StateType> &max)
  {
     _max.resize(max.size(),0.L);
     for(unsigned int i = 0; i < max.size(); i++)
     {
       _max[i] = max[i];
     }
  }

  template <typename CoeffType>
  inline
  const std::vector<CoeffType> DiUTPdf<CoeffType>::min() const
  {
      return _min;
  }

  template <typename CoeffType>
  inline
  const std::vector<CoeffType> DiUTPdf<CoeffType>::max() const
  {
      return _max;
  }


  template <typename CoeffType>
  inline
  const CoeffType DiUTPdf<CoeffType>::value(unsigned int ip) const
  {
      std::vector<CoeffType> out;
      out.resize(_min.size(),0.L);
      CoeffType sum(0.L);
      for(unsigned int i = 0; i < _min.size(); i++)
      {
          out[i] = (_min[i] + _max[i])/2.L;
          sum += out[i];
      }
      for(unsigned int i = 0; i < _min.size(); i++)
      {
          out[i] /= sum;
      }

      return out[ip];
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void DiUTPdf<CoeffType>::set_parameters(const std::vector<StateType> &pars)
  {
      antioch_assert_equal_to(pars.size()%2,0);

     _min.resize(pars.size()/2,0.L);
     _max.resize(pars.size()/2,0.L);
     for(unsigned int i = 0; i < pars.size()/2; i++)
     {
       _min[i] = pars[i];
       _max[i] = pars[i + pars.size()/2];
     }
  }

  template <typename CoeffType>
  inline
  void DiUTPdf<CoeffType>::print(std::ostream &out)  const
  {
     out << "DiUT(" ;
     out << _min[0];

     for(unsigned int ibr = 1; ibr < _min.size(); ibr++)
     {
        out << "," << _min[ibr];
     }
     out << "; " << _max[0];
     for(unsigned int ibr = 1; ibr < _min.size(); ibr++)
     {
        out << "," << _max[ibr];
     }

     out << ")";
  }
}

#endif
