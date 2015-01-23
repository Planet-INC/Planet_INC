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
#include "planet/dirg_pdf.h"
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
               << "Error in pdf dirg" << std::endl
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
  std::vector<Scalar> v;
  std::vector<Scalar> dv;
  std::vector<Scalar> pars;
  v.push_back(0.5);
  v.push_back(0.2);
  v.push_back(0.3);
  dv.push_back(0.3);
  dv.push_back(0.2);
  dv.push_back(0.2);
  for(unsigned int i = 0; i < v.size(); i++)
  {
    pars.push_back(v[i]);
  }
  for(unsigned int i = 0; i < dv.size(); i++)
  {
    pars.push_back(dv[i]);
  }
  Planet::DirGPdf<Scalar> dirg1(v,dv);
  Planet::DirGPdf<Scalar> dirg2,dirg3;
  dirg2.set_v(v);
  dirg2.set_dv(dv);
  dirg3.set_parameters(pars);


  int return_flag(0);
  return_flag = check(dirg1.v()[0],v[0],"dirg1 v[0]") ||
                check(dirg1.v()[1],v[1],"dirg1 v[1]") ||
                check(dirg1.v()[2],v[2],"dirg1 v[2]") ||
                check(dirg1.dv()[0],dv[0],"dirg1 dv[0]") ||
                check(dirg1.dv()[1],dv[1],"dirg1 dv[1]") ||
                check(dirg1.dv()[2],dv[2],"dirg1 dv[2]") ||
                check(dirg1.value(0),v[0],"dirg1 value(0)") ||
                check(dirg1.value(1),v[1],"dirg1 value(1)") ||
                check(dirg1.value(2),v[2],"dirg1 value(2)") ||
                check(dirg2.v()[0],v[0],"dirg2 v[0]") ||
                check(dirg2.v()[1],v[1],"dirg2 v[1]") ||
                check(dirg2.v()[2],v[2],"dirg2 v[2]") ||
                check(dirg2.dv()[0],dv[0],"dirg2 dv[0]") ||
                check(dirg2.dv()[1],dv[1],"dirg2 dv[1]") ||
                check(dirg2.dv()[2],dv[2],"dirg2 dv[2]") ||
                check(dirg2.value(0),v[0],"dirg2 value(0)") ||
                check(dirg2.value(1),v[1],"dirg2 value(1)") ||
                check(dirg2.value(2),v[2],"dirg2 value(2)") ||
                check(dirg3.v()[0],v[0],"dirg3 v[0]") ||
                check(dirg3.v()[1],v[1],"dirg3 v[1]") ||
                check(dirg3.v()[2],v[2],"dirg3 v[2]") ||
                check(dirg3.dv()[0],dv[0],"dirg3 dv[0]") ||
                check(dirg3.dv()[1],dv[1],"dirg3 dv[1]") ||
                check(dirg3.dv()[2],dv[2],"dirg3 dv[2]") ||
                check(dirg3.value(0),v[0],"dirg3 value(0)") ||
                check(dirg3.value(1),v[1],"dirg3 value(1)") ||
                check(dirg3.value(2),v[2],"dirg3 value(2)");

  return return_flag;
}


int main()
{

  return (tester<float>()  || 
          tester<double>() || 
          tester<long double>());
}
