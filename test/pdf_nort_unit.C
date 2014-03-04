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
#include "planet/nort_pdf.h"
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
               << "Error in pdf nort" << std::endl
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
  Scalar mu(1.5),sigma(3.),min(-1.),max(2.);
  std::vector<Scalar> pars;
  pars.push_back(mu);
  pars.push_back(sigma);
  pars.push_back(min);
  pars.push_back(max);
  Planet::NorTPdf<Scalar> nort1(mu,sigma,min,max);
  Planet::NorTPdf<Scalar> nort2,nort3;
  nort2.set_mu(mu);
  nort2.set_sigma(sigma);
  nort2.set_min(min);
  nort2.set_max(max);
  nort3.set_parameters(pars);


  int return_flag(0);
  return_flag = check(nort1.mu(),mu,"nort1 mu") ||
                check(nort1.sigma(),sigma,"nort1 sigma") ||
                check(nort1.value(),mu,"nort1 value") ||
                check(nort1.min(),min,"nort1 min") ||
                check(nort1.max(),max,"nort1 max") ||
                check(nort2.mu(),mu,"nort2 mu") ||
                check(nort2.sigma(),sigma,"nort2 sigma") ||
                check(nort2.value(),mu,"nort2 value") ||
                check(nort2.min(),min,"nort1 min") ||
                check(nort2.max(),max,"nort1 max") ||
                check(nort3.mu(),mu,"nort3 mu") ||
                check(nort3.sigma(),sigma,"nort3 sigma") ||
                check(nort3.value(),mu,"nort3 value") ||
                check(nort3.min(),min,"nort1 min") ||
                check(nort3.max(),max,"nort1 max");

  return return_flag;
}


int main()
{

  return (tester<float>()  || 
          tester<double>() || 
          tester<long double>());
}
