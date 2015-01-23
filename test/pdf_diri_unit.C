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
#include "planet/diri_pdf.h"
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
               << "Error in pdf diri" << std::endl
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
  Scalar r(3.);
  std::vector<Scalar> v;
  std::vector<Scalar> pars;
  v.push_back(0.5);
  v.push_back(0.2);
  v.push_back(0.3);
  for(unsigned int i = 0; i < v.size(); i++)
  {
    pars.push_back(v[i]);
  }
  pars.push_back(r);
  Planet::DiriPdf<Scalar> diri1(v,r);
  Planet::DiriPdf<Scalar> diri2,diri3;
  diri2.set_v(v);
  diri2.set_r(r);
  diri3.set_parameters(pars);


  int return_flag(0);
  return_flag = check(diri1.v()[0],v[0],"diri1 v[0]") ||
                check(diri1.v()[1],v[1],"diri1 v[1]") ||
                check(diri1.v()[2],v[2],"diri1 v[2]") ||
                check(diri1.r(),r,"diri1 r") ||
                check(diri1.value(0),v[0],"diri1 value(0)") ||
                check(diri1.value(1),v[1],"diri1 value(1)") ||
                check(diri1.value(2),v[2],"diri1 value(2)") ||
                check(diri2.v()[0],v[0],"diri2 v[0]") ||
                check(diri2.v()[1],v[1],"diri2 v[1]") ||
                check(diri2.v()[2],v[2],"diri2 v[2]") ||
                check(diri2.r(),r,"diri2 r") ||
                check(diri2.value(0),v[0],"diri2 value(0)") ||
                check(diri2.value(1),v[1],"diri2 value(1)") ||
                check(diri2.value(2),v[2],"diri2 value(2)") ||
                check(diri3.v()[0],v[0],"diri3 v[0]") ||
                check(diri3.v()[1],v[1],"diri3 v[1]") ||
                check(diri3.v()[2],v[2],"diri3 v[2]") ||
                check(diri3.r(),r,"diri3 r") ||
                check(diri3.value(0),v[0],"diri3 value(0)") ||
                check(diri3.value(1),v[1],"diri3 value(1)") ||
                check(diri3.value(2),v[2],"diri3 value(2)");

  return return_flag;
}


int main()
{

  return (tester<float>()  || 
          tester<double>() || 
          tester<long double>());
}
