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

#ifndef PLANET_MATH_FUNCTIONS_H
#define PLANET_MATH_CUNCTIONS_H

namespace Planet 
{
  namespace Functions
  {
     /*!
      * looking for index
      */
    template<typename CoeffType, typename VectorCoeffType>
    inline
    unsigned int find_floor_index(const VectorCoeffType &alt, const CoeffType &value)
    {
      unsigned int iz;
      for(iz = 0; iz < alt.size() - 1; iz++)
      {
         if((alt[iz] <= value && alt[iz + 1] >  value) ||
            (alt[iz] >  value && alt[iz + 1] <= value))break;
      }
      return iz;
    }

    /*!
     * linear evaluation
     */
    template<typename CoeffType, typename VectorCoeffType>
    inline
    CoeffType linear_evaluation(const VectorCoeffType &alt, const VectorCoeffType &data, const CoeffType &value)
    {
      CoeffType interpolation;
      if(value <= alt[0])
      {
        interpolation = data[0];
      }else if(value >= alt.back())
      {
        interpolation = data.back();
      }else
      {
        unsigned int iz = find_floor_index(alt,value);
        CoeffType a = (data[iz + 1] - data[iz]) / (alt[iz + 1] - alt[iz]);
        CoeffType b = data[iz] - a * alt[iz];
        interpolation = a * value + b;
      }

      return interpolation;
    }

    /*!
     * linear evaluation of derivative
     */
    template<typename CoeffType, typename VectorCoeffType>
    inline
    CoeffType linear_evaluation_dz(const VectorCoeffType &alt, const VectorCoeffType &data, const CoeffType &value)
    {
      unsigned int iz = find_floor_index(alt,value);
      return (data[iz + 1] - data[iz]) / (alt[iz] - alt[iz + 1]);
    }

  } // end namespace Functions
} // end namespace Planet

#endif //PLANET_MATH_FUNCTIONS_H
