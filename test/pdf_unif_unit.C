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
#include "antioch/physical_constants.h"
#include "antioch/cmath_shims.h"
//Planet
#include "planet/unif_pdf.h"
//C++
#include <limits>
#include <string>
#include <vector>
#include <iomanip>


template<typename Scalar>
int check(const Scalar &test, const Scalar &ref, const std::string &blabla)
{
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100.;
  if(Antioch::ant_abs(test - ref)/ref > tol)
  {
     std::cout << std::scientific << std::setprecision(20)
               << "Error in pdf unif" << std::endl
               << "test is "          << blabla << std::endl
               << "given value = "    << test << std::endl
               << "exact value = "    << ref << std::endl
               << "relative error = " << Antioch::ant_abs(test - ref)/ref << std::endl
               << "tolerance = "      << tol << std::endl;
     return 1;
  }
  return 0;
}

template<typename Scalar>
int tester()
{
  Scalar min(1.5),max(3.);
  std::vector<Scalar> pars;
  pars.push_back(min);
  pars.push_back(max);
  Planet::UnifPdf<Scalar> unif1(min,max);
  Planet::UnifPdf<Scalar> unif2,unif3;
  unif2.set_min(min);
  unif2.set_max(max);
  unif3.set_parameters(pars);
  Scalar value = (min + max)/2.L;

  int return_flag(0);
  return_flag = check(unif1.min(),min,"unif1 min") ||
                check(unif1.max(),max,"unif1 max") ||
                check(unif1.value(),value,"unif1 value") ||
                check(unif2.min(),min,"unif2 min") ||
                check(unif2.max(),max,"unif2 max") ||
                check(unif2.value(),value,"unif2 value") ||
                check(unif3.min(),min,"unif3 min") ||
                check(unif3.max(),max,"unif3 max") ||
                check(unif3.value(),value,"unif3 value");

  return return_flag;
}


int main()
{

  return (tester<float>()  || 
          tester<double>() || 
          tester<long double>());
}
