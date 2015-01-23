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
#include "planet/diun_pdf.h"
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
               << "Error in pdf diun" << std::endl
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
  unsigned int n(3);
  std::vector<Scalar> v;
  std::vector<Scalar> pars;
  v.push_back(1.L/3.);
  v.push_back(1.L/3.);
  v.push_back(1.L/3.);
  pars.push_back((Scalar)n);
  Planet::DiUnPdf<Scalar> diun1(n);
  Planet::DiUnPdf<Scalar> diun2,diun3;
  diun2.set_n(n);
  diun3.set_parameters(pars);


  int return_flag(0);
  return_flag = check(diun1.n(),n,"diun1 n") ||
                check(diun1.value(0),v[0],"diun1 value(0)") ||
                check(diun1.value(1),v[1],"diun1 value(1)") ||
                check(diun1.value(2),v[2],"diun1 value(2)") ||
                check(diun2.n(),n,"diun2 n") ||
                check(diun2.value(0),v[0],"diun2 value(0)") ||
                check(diun2.value(1),v[1],"diun2 value(1)") ||
                check(diun2.value(2),v[2],"diun2 value(2)") ||
                check(diun3.n(),n,"diun3 n") ||
                check(diun3.value(0),v[0],"diun3 value(0)") ||
                check(diun3.value(1),v[1],"diun3 value(1)") ||
                check(diun3.value(2),v[2],"diun3 value(2)");

  return return_flag;
}


int main()
{

  return (tester<float>()  || 
          tester<double>() || 
          tester<long double>());
}
