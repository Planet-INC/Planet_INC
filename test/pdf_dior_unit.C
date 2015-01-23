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
#include "planet/dior_pdf.h"
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
               << "Error in pdf dior" << std::endl
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
  std::vector<unsigned int> n;
  std::vector<Scalar> pars;
  n.push_back(2);
  n.push_back(1);
  n.push_back(3);
  for(unsigned int i = 0; i < n.size(); i++)
  {
    pars.push_back((Scalar)n[i]);
  }
  Planet::DiOrPdf<Scalar> dior1(n);
  Planet::DiOrPdf<Scalar> dior2,dior3;
  dior2.set_n(n);
  dior3.set_parameters(pars);
  std::vector<Scalar> v;
  for(unsigned int i = 0; i < n.size(); i++)
  {
      v.push_back( 2. * ((Scalar)n.size() + 1.L - (Scalar)n[i]) / ((Scalar)n.size() * (Scalar)(n.size() + 1)) );
  }


  int return_flag(0);
  return_flag = check(dior1.n()[0],n[0],"dior1 n[0]") ||
                check(dior1.n()[1],n[1],"dior1 n[1]") ||
                check(dior1.n()[2],n[2],"dior1 n[2]") ||
                check(dior1.value(0),v[0],"dior1 value(0)") ||
                check(dior1.value(1),v[1],"dior1 value(1)") ||
                check(dior1.value(2),v[2],"dior1 value(2)") ||
                check(dior2.n()[0],n[0],"dior2 n[0]") ||
                check(dior2.n()[1],n[1],"dior2 n[1]") ||
                check(dior2.n()[2],n[2],"dior2 n[2]") ||
                check(dior2.value(0),v[0],"dior2 value(0)") ||
                check(dior2.value(1),v[1],"dior2 value(1)") ||
                check(dior2.value(2),v[2],"dior2 value(2)") ||
                check(dior3.n()[0],n[0],"dior3 n[0]") ||
                check(dior3.n()[1],n[1],"dior3 n[1]") ||
                check(dior3.n()[2],n[2],"dior3 n[2]") ||
                check(dior3.value(0),v[0],"dior3 value(0)") ||
                check(dior3.value(1),v[1],"dior3 value(1)") ||
                check(dior3.value(2),v[2],"dior3 value(2)");

  return return_flag;
}


int main()
{

  return (tester<float>()  || 
          tester<double>() || 
          tester<long double>());
}
