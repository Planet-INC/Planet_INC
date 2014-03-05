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
#include "planet/diut_pdf.h"
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
               << "Error in pdf diut" << std::endl
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
  std::vector<Scalar> min;
  std::vector<Scalar> max;
  std::vector<Scalar> pars;
  min.push_back(0.1);
  min.push_back(0.2);
  min.push_back(0.3);
  max.push_back(0.5);
  max.push_back(0.8);
  max.push_back(0.4);
  for(unsigned int i = 0; i < min.size(); i++)
  {
    pars.push_back(min[i]);
  }
  for(unsigned int i = 0; i < max.size(); i++)
  {
    pars.push_back(max[i]);
  }
  Planet::DiUTPdf<Scalar> diut1(min,max);
  Planet::DiUTPdf<Scalar> diut2,diut3;
  diut2.set_min(min);
  diut2.set_max(max);
  diut3.set_parameters(pars);

  std::vector<Scalar> v;
  Scalar sum(0.L);
  for(unsigned int i = 0; i < min.size(); i++)
  {
    v.push_back((min[i] + max[i])/2.L);
    sum += v[i];
  }
  for(unsigned int i = 0; i < v.size(); i++)
  {
    v[i] /= sum;
  }
  

  int return_flag(0);
  return_flag = check(diut1.min()[0],min[0],"diut1 min[0]") ||
                check(diut1.min()[1],min[1],"diut1 min[1]") ||
                check(diut1.min()[2],min[2],"diut1 min[2]") ||
                check(diut1.max()[0],max[0],"diut1 max[0]") ||
                check(diut1.max()[1],max[1],"diut1 max[1]") ||
                check(diut1.max()[2],max[2],"diut1 max[2]") ||
                check(diut1.value(0),v[0],"diut1 value(0)") ||
                check(diut1.value(1),v[1],"diut1 value(1)") ||
                check(diut1.value(2),v[2],"diut1 value(2)") ||
                check(diut2.min()[0],min[0],"diut2 min[0]") ||
                check(diut2.min()[1],min[1],"diut2 min[1]") ||
                check(diut2.min()[2],min[2],"diut2 min[2]") ||
                check(diut2.max()[0],max[0],"diut2 max[0]") ||
                check(diut2.max()[1],max[1],"diut2 max[1]") ||
                check(diut2.max()[2],max[2],"diut2 max[2]") ||
                check(diut2.value(0),v[0],"diut2 value(0)") ||
                check(diut2.value(1),v[1],"diut2 value(1)") ||
                check(diut2.value(2),v[2],"diut2 value(2)") ||
                check(diut3.min()[0],min[0],"diut3 min[0]") ||
                check(diut3.min()[1],min[1],"diut3 min[1]") ||
                check(diut3.min()[2],min[2],"diut3 min[2]") ||
                check(diut3.max()[0],max[0],"diut3 max[0]") ||
                check(diut3.max()[1],max[1],"diut3 max[1]") ||
                check(diut3.max()[2],max[2],"diut3 max[2]") ||
                check(diut3.value(0),v[0],"diut3 value(0)") ||
                check(diut3.value(1),v[1],"diut3 value(1)") ||
                check(diut3.value(2),v[2],"diut3 value(2)");

  return return_flag;
}


int main()
{

  return (tester<float>()  || 
          tester<double>() || 
          tester<long double>());
}
