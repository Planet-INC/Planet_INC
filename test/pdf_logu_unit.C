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
#include "planet/logu_pdf.h"
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
               << "Error in pdf logu" << std::endl
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
  Scalar min(1.5e3),max(3.2e8);
  std::vector<Scalar> pars;
  pars.push_back(min);
  pars.push_back(max);
  Planet::LogUPdf<Scalar> logu1(min,max);
  Planet::LogUPdf<Scalar> logu2,logu3;
  logu2.set_min(min);
  logu2.set_max(max);
  logu3.set_parameters(pars);
  Scalar value = Antioch::ant_sqrt(min*max);

  int return_flag(0);
  return_flag = check(logu1.min(),min,"logu1 min") ||
                check(logu1.max(),max,"logu1 max") ||
                check(logu1.value(),value,"logu1 value") ||
                check(logu2.min(),min,"logu2 min") ||
                check(logu2.max(),max,"logu2 max") ||
                check(logu2.value(),value,"logu2 value") ||
                check(logu3.min(),min,"logu3 min") ||
                check(logu3.max(),max,"logu3 max") ||
                check(logu3.value(),value,"logu3 value");

  return return_flag;
}


int main()
{

  return (tester<float>()  || 
          tester<double>() || 
          tester<long double>());
}
