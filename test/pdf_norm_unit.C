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
#include "planet/norm_pdf.h"
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
               << "Error in pdf norm" << std::endl
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
  Scalar mu(1.5),sigma(3.);
  std::vector<Scalar> pars;
  pars.push_back(mu);
  pars.push_back(sigma);
  Planet::NormPdf<Scalar> norm1(mu,sigma);
  Planet::NormPdf<Scalar> norm2,norm3;
  norm2.set_mu(mu);
  norm2.set_sigma(sigma);
  norm3.set_parameters(pars);


  int return_flag(0);
  return_flag = check(norm1.mu(),mu,"norm1 mu") ||
                check(norm1.sigma(),sigma,"norm1 sigma") ||
                check(norm1.value(),mu,"norm1 value") ||
                check(norm2.mu(),mu,"norm2 mu") ||
                check(norm2.sigma(),sigma,"norm2 sigma") ||
                check(norm2.value(),mu,"norm2 value") ||
                check(norm3.mu(),mu,"norm3 mu") ||
                check(norm3.sigma(),sigma,"norm3 sigma") ||
                check(norm3.value(),mu,"norm3 value");

  return return_flag;
}


int main()
{

  return (tester<float>()  || 
          tester<double>() || 
          tester<long double>());
}
