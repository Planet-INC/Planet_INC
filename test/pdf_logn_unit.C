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
#include "planet/logn_pdf.h"
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
               << "Error in pdf logn" << std::endl
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
  Scalar mu(1.5),f(3.);
  std::vector<Scalar> pars;
  pars.push_back(mu);
  pars.push_back(f);
  Planet::LogNPdf<Scalar> logn1(mu,f);
  Planet::LogNPdf<Scalar> logn2,logn3;
  logn2.set_mu(mu);
  logn2.set_f(f);
  logn3.set_parameters(pars);


  int return_flag(0);
  return_flag = check(logn1.mu(),mu,"logn1 mu") ||
                check(logn1.f(),f,"logn1 f") ||
                check(logn1.value(),mu,"logn1 value") ||
                check(logn2.mu(),mu,"logn2 mu") ||
                check(logn2.f(),f,"logn2 f") ||
                check(logn2.value(),mu,"logn2 value") ||
                check(logn3.mu(),mu,"logn3 mu") ||
                check(logn3.f(),f,"logn3 f") ||
                check(logn3.value(),mu,"logn3 value");

  return return_flag;
}


int main()
{

  return (tester<float>()  || 
          tester<double>() || 
          tester<long double>());
}
