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

#ifndef PLANET_CONSTANTS_H
#define PLANET_CONSTANTS_H

namespace Planet 
{
  namespace Constants
  {

    namespace Universal
    {
      /*!
       * Universal gravity constant in m3.kg-1.s-2
       * error is 0.00080e-11 m3.kg-1.s-2
       */
    template<typename CoeffType>
    inline
    CoeffType G()
    {
      return 6.67384e-11;
    }

    /*!
     * Boltzmann constant in J.K-1
     * error is 0.0000013e-23 J.K-1
     */
    template<typename CoeffType>
    inline
    CoeffType kb()
    {
      return 1.3806488e-23;
    }

    } //end namespace Universal

    /*!
     * Gravitationnal constant, giving
     * radius of body (km),
     * altitude (in km) and 
     * mass (kg)
     */
    template<typename CoeffType>
    inline
    CoeffType g(const CoeffType &radius, const CoeffType &alt, const CoeffType &mass)
    {
       return Constants::Universal::G<CoeffType>() * mass / (CoeffType(1e6L) * (radius + alt) * (radius + alt));
    }

    namespace Saturn
    {
    /*!
     * Sun-Saturn distance in AU
     *  error unknown
     *  from www.solarviews.com
     */
    template<typename CoeffType>
    inline
    CoeffType d_Sun()
    {
      return 9.5388;
    }
    } //end namespace Saturn


    namespace Titan
    {
    /*!
     * mTitan in kg,
     * error is 2.0155e22 kg
     * from Jacobson et al., AJ, 132:2520--2526, 2006
     */
    template<typename CoeffType>
    inline
    CoeffType mass()
    {
      return 1.34520029e23L;
    }

    /*!
     * RTitan in km,
     * error is 2 km
     * from Jacobson et al., AJ, 132:2520--2526, 2006
     *          from Thomas et al., 1995
     */
    template<typename CoeffType>
    inline
    CoeffType radius()
    {
      return 2575.5L;
    }
    } //end namespace Titan


    namespace Convention
    {
    /*!
     * Normal Pressure: 101,325 Pa (1 atm)
     */
    template<typename CoeffType>
    inline
    CoeffType P_normal()
    {
      return 1.01325e5L;
    }

    /*!
     * Normal Temperature : 300 K
     */
    template<typename CoeffType>
    inline
    CoeffType T_normal()
    {
      return 300.L;
    }

    /*!
     * Standard Pressure: 1e5 Pa (1 bar)
     */
    template<typename CoeffType>
    inline
    CoeffType P_standard()
    {
      return 1e5L;
    }

    /*!
     * Standard Temperature: 273.15 K (0 °C)
     */
    template<typename CoeffType>
    inline
    CoeffType T_standard()
    {
      return 273.15L;
    }
    } //end namespace Convention

  } // end namespace Constants
} // end namespace Planet

#endif //PLANET_MATH_CONSTANTS_H
