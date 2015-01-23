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

#ifndef PLANET_CROSS_SECTION_H
#define PLANET_CROSS_SECTION_H

//Antioch
#include "antioch/sigma_bin_converter.h"

//Planet

//C++

namespace Planet
{
  template<typename VectorCoeffType>
  class CrossSection
  {
      private:
        VectorCoeffType _abscissa;
        VectorCoeffType _cross_section;
        VectorCoeffType _cross_section_on_custom_grid;

        Antioch::SigmaBinConverter<VectorCoeffType> _converter;

      public:

        CrossSection();
        CrossSection(const CrossSection<VectorCoeffType> &rhs);
        CrossSection(const VectorCoeffType &x, const VectorCoeffType &y);
        ~CrossSection();

        //!sets the abscissa
        template<typename VectorStateType>
        void set_abscissa(const VectorStateType &x);

        //!sets the cross section
        template<typename VectorStateType>
        void set_cross_section(const VectorStateType &cs);

        template<typename VectorStateType>
        void update_cross_section(const VectorStateType &custom_x);

        //!\return the abscissa
        const VectorCoeffType abscissa()      const;

        //!\return the cross-section
        const VectorCoeffType cross_section() const;

        //!\return the cross-section on custom grid
        const VectorCoeffType &cross_section_on_custom_grid() const;
  };

  template<typename VectorCoeffType>
  inline
  CrossSection<VectorCoeffType>::CrossSection(const VectorCoeffType &x, const VectorCoeffType &y):
  _abscissa(x),
  _cross_section(y)
  {
    return;
  }

  template<typename VectorCoeffType>
  inline
  CrossSection<VectorCoeffType>::~CrossSection()
  {
    return;
  }

  template<typename VectorCoeffType>
  inline
  CrossSection<VectorCoeffType>::CrossSection(const CrossSection<VectorCoeffType> &rhs):
  _abscissa(rhs.abscissa()),
  _cross_section(rhs.cross_section())
  {
    return;
  }

  template<typename VectorCoeffType>
  template<typename VectorStateType>
  inline
  void CrossSection<VectorCoeffType>::set_abscissa(const VectorStateType &x)
  {
      _abscissa = x;
  }

  template<typename VectorCoeffType>
  template<typename VectorStateType>
  inline
  void CrossSection<VectorCoeffType>::set_cross_section(const VectorStateType &cs)
  {
     _cross_section = cs;
  }

  template<typename VectorCoeffType>
  template<typename VectorStateType>
  inline
  void CrossSection<VectorCoeffType>::update_cross_section(const VectorStateType &custom_x)
  {
    _cross_section_on_custom_grid.resize(custom_x.size());
    _converter.y_on_custom_grid(_abscissa, _cross_section,
                                 custom_x, _cross_section_on_custom_grid);

    for(unsigned int i = 0; i < _cross_section_on_custom_grid.size(); i++)
    {
       if(_cross_section_on_custom_grid[i] < typename Antioch::value_type<VectorStateType>::type(0.L))
                _cross_section_on_custom_grid[i] *= typename Antioch::value_type<VectorStateType>::type(-1.L);
    }

     return;
  }

  template<typename VectorCoeffType>
  inline
  const VectorCoeffType CrossSection<VectorCoeffType>::abscissa() const
  {
     return _abscissa;
  }

  template<typename VectorCoeffType>
  inline
  const VectorCoeffType CrossSection<VectorCoeffType>::cross_section() const
  {
     return _cross_section;
  }

  template<typename VectorCoeffType>
  inline
  const VectorCoeffType &CrossSection<VectorCoeffType>::cross_section_on_custom_grid() const
  {
     return _cross_section_on_custom_grid;
  }

}

#endif
