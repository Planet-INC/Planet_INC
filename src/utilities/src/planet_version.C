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

// planet
#include "planet/planet_version.h"

// C++
#include <iostream>

namespace Planet
{

  void planet_version_stdout()
  {
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "planet Package: Version = " << PLANET_LIB_VERSION;
    std::cout << " (" << get_planet_version() << ")" << std::endl << std::endl;
  
    std::cout << PLANET_LIB_RELEASE << std::endl << std::endl;
  
    std::cout << "Build Date   = " << PLANET_BUILD_DATE     << std::endl;
    std::cout << "Build Host   = " << PLANET_BUILD_HOST     << std::endl;
    std::cout << "Build User   = " << PLANET_BUILD_USER     << std::endl;
    std::cout << "Build Arch   = " << PLANET_BUILD_ARCH     << std::endl;
    std::cout << "Build Rev    = " << PLANET_BUILD_VERSION  << std::endl << std::endl;
  
    std::cout << "C++ Config   = " << PLANET_CXX << " " << PLANET_CXXFLAGS << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
  
    return;
  }

  int get_planet_version()
  {
    /* Note: return format follows the versioning convention xx.yy.zz where
   
       xx = major version number
       yy = minor version number
       zz = micro version number
     
       For example:
       v.   0.23  -> 002300 = 2300
       v   0.23.1 -> 002301 = 2301
       v. 10.23.2 -> 102302         */

    int major_version = 0;
    int minor_version = 0;
    int micro_version = 0;

#ifdef PLANET_MAJOR_VERSION
    major_version = PLANET_MAJOR_VERSION;
#endif

#ifdef PLANET_MINOR_VERSION
    minor_version = PLANET_MINOR_VERSION;
#endif

#ifdef PLANET_MICRO_VERSION
    micro_version = PLANET_MICRO_VERSION;
#endif
      
    return major_version*10000 + minor_version*100 + micro_version;
  }

} // end namespace Planet
