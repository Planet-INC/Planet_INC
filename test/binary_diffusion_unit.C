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
#include "planet/binary_diffusion.h"
#include "planet/planet_constants.h"
//C++
#include <limits>
#include <string>
#include <iomanip>


template<typename Scalar>
int check(const Scalar &test, const Scalar &ref, const Scalar &tol, const std::string &model)
{
  if(Antioch::ant_abs(test - ref)/ref > tol)
  {
     std::cout << std::scientific << std::setprecision(20)
               << "Error in binary diffusion calculations" << std::endl
               << "model is " << model << std::endl
               << "calculated coefficient = " << test << std::endl
               << "solution = " << ref << std::endl
               << "relative error = " << Antioch::ant_abs(test - ref)/ref << std::endl
               << "tolerance = " << tol << std::endl;
     return 1;
  }
  return 0;
}

template<typename Scalar>
int tester()
{

  Scalar p11(0.1783),p12(1.81);
  Scalar p21(1.04e-5),p22(1.76);
  Scalar p31(5.73e16),p32(0.5);

  Planet::BinaryDiffusion<Scalar> N2N2(   0,  0, p11, p12, Planet::DiffusionType::Massman);
  Planet::BinaryDiffusion<Scalar> N2CH4(  0,  1, p21, p22, Planet::DiffusionType::Wakeham);
  Planet::BinaryDiffusion<Scalar> CH4CH4( 1,  1, p31, p32, Planet::DiffusionType::Wilson);

  Scalar T(1500.),P(1e5);
  Scalar n = P / (T * Planet::Constants::Universal::kb<Scalar>());

  Scalar n2n2 = p11 * Planet::Constants::Convention::P_normal<Scalar>() / P * Antioch::ant_pow(T/Planet::Constants::Convention::T_standard<Scalar>(),p12);
  Scalar n2ch4 = p21 * Antioch::ant_pow(T,p22)* Planet::Constants::Convention::P_normal<Scalar>() / P;
  Scalar ch4ch4 = p31 * Antioch::ant_pow(T,p32)/n;

  Scalar nn = N2N2.binary_coefficient(T,P);
  Scalar nc = N2CH4.binary_coefficient(T,P);
  Scalar cc = CH4CH4.binary_coefficient(T,P);

  int return_flag(0);
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100.;
  return_flag = check(nn,n2n2,tol,"Massman")  || 
                check(nc,n2ch4,tol,"Wilson") ||
                check(cc,ch4ch4,tol,"Wakeham");

  return return_flag;
}


int main()
{

  return (tester<float>()  || 
          tester<double>() || 
          tester<long double>());
}
