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

//Antioch
#include "antioch/cmath_shims.h"
//Planet
#include "planet/altitude.h"
//C++
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

template<typename Scalar>
int check(const Scalar &test, const Scalar &ref, const Scalar &tol, const std::string &model)
{
  if(Antioch::ant_abs(test - ref)/ref > tol)
  {
     std::cout << std::scientific << std::setprecision(20)
               << "Error in binary diffusion calculations" << std::endl
               << "model is " << model << std::endl
               << "calculated coefficient = " << test << std::endl
               << "solution = " << ref << std::endl
               << "relative error = " << Antioch::ant_abs(test - ref)/ref << std::endl
               << "tolerance = " << tol << std::endl;
     return 1;
  }
  return 0;
}

template<typename Scalar>
int tester()
{
  Planet::Altitude<Scalar,std::vector<Scalar> > altitude(600,1400,10);

  Scalar min(600.L),max(1400.L),step(10.L);  

  int return_flag(0);
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100.;
  return_flag = check(altitude.alt_min(),min,tol,"minimum altitude") ||
                check(altitude.alt_max(),max,tol,"maximum altitude") ||
                check(altitude.alt_step(),step,tol,"step altitude");

  unsigned int iz(0);
  for(Scalar alt = min; alt <= max; alt += step)
  {
      return_flag = return_flag || check(altitude.altitudes()[iz],alt,tol,"altitude");
      iz++;
  }

  return return_flag;
}


int main()
{
  return (tester<float>()  ||
          tester<double>() ||
          tester<long double>());
}
