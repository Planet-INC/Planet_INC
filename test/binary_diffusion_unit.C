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
#include "antioch/units.h"

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
  Scalar dist = (std::abs(ref) < tol)?Antioch::ant_abs(test - ref):
                                      Antioch::ant_abs(test - ref)/ref;
  if(dist > tol)
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

template <typename Scalar>
Scalar Dij(const Scalar & T, const Scalar & P, const Scalar & D01, const Scalar & beta)
{
  return D01 * Planet::Constants::Convention::P_normal<Scalar>() / P * Antioch::ant_pow(T/Planet::Constants::Convention::T_standard<Scalar>(),beta);
}

template <typename Scalar>
Scalar dDij_dT(const Scalar & T, const Scalar & P, const Scalar & D01, const Scalar & beta)
{
  return Dij(T,P,D01,beta) * (beta - 1.L)/T;
}

template <typename Scalar>
Scalar dDij_dn(const Scalar & T, const Scalar & P, const Scalar & D01, const Scalar & beta)
{
  return -Dij(T,P,D01,beta) * (T*Planet::Constants::Universal::kb<Scalar>())/(P*P);
}

template<typename Scalar>
int tester()
{

  Scalar p11(0.1783L),p12(1.81L);
  Scalar p21(1.04e-5L),p22(1.76L);
  Scalar p31(5.73e16L),p32(0.5L);

  Planet::BinaryDiffusion<Scalar> N2N2(   0,  0, p11, p12, Planet::DiffusionType::Massman);
  Planet::BinaryDiffusion<Scalar> N2CH4(  0,  1, p21, p22, Planet::DiffusionType::Wakeham);
  Planet::BinaryDiffusion<Scalar> CH4CH4( 1,  1, p31, p32, Planet::DiffusionType::Wilson);

// on range of temperatures, common pressure of
// measurements is 1 atm
  Scalar P(1.L); // in atm
  P *= Antioch::Units<Scalar>("atm").get_SI_factor(); // in Pa

  int return_flag(0);
  const Scalar tol = (std::numeric_limits<Scalar>::epsilon() < 1e-17)?std::numeric_limits<Scalar>::epsilon() * 20.:
                                                                      std::numeric_limits<Scalar>::epsilon() * 10.;
  for(Scalar T = 100.; T < 3500.; T += 100.)
  {

    Scalar n = P / (T * Planet::Constants::Universal::kb<Scalar>());

    Scalar n2n2   = Dij(T,P,p11,p12);
    Scalar n2ch4  = p21 * Antioch::ant_pow(T,p22);
    Scalar ch4ch4 = p31 * Antioch::ant_pow(T,p32) / n;

    Scalar dn2n2_dT   = dDij_dT(T,P,p11,p12);
    Scalar dn2n2_dn   = dDij_dn(T,P,p11,p12);
    Scalar dch4ch4_dT = p32 * p31 / n * Antioch::ant_pow(T,p32 - 1.L);
    Scalar dch4ch4_dn = - p31 * Antioch::ant_pow(T,p32) / (n * n);

    return_flag = check(N2N2.binary_coefficient(T,P),n2n2,tol,"Massman binary coefficient")    ||
                  check(N2CH4.binary_coefficient(T,P),n2ch4,tol,"Wakeham binary coefficient")  ||
                  check(CH4CH4.binary_coefficient(T,P),ch4ch4,tol,"Wilson binary coefficient") ||
                  check(N2N2.binary_coefficient_deriv_T(T,P),dn2n2_dT,tol,"Massman binary coefficient derived to T")    ||
                  check(CH4CH4.binary_coefficient_deriv_T(T,P),dch4ch4_dT,tol,"Wilson binary coefficient derived to T") ||
                  check(N2N2.binary_coefficient_deriv_n(T,P,n),dn2n2_dn,tol,"Massman binary coefficient derived to n")    ||
                  check(CH4CH4.binary_coefficient_deriv_n(T,P,n),dch4ch4_dn,tol,"Wilson binary coefficient derived to n") ||
                  return_flag;
  }

  // golden values
  Scalar T(201.L);
  P = 0.53L;
  Scalar n = P / (T * Planet::Constants::Universal::kb<Scalar>());
  Scalar truth     = 1.95654918135100954675984465985829467e4L;
  Scalar dtruth_dT = 7.88460117857869518843519489793641134e1L;
  Scalar dtruth_dn = -1.0244580436868377275736362395882496e-16L;

  return_flag = check(N2N2.binary_coefficient(T,P),truth,tol,"Massman binary coefficient") ||
                check(N2N2.binary_coefficient_deriv_T(T,P),dtruth_dT,tol,"Massman binary coefficient derived to T") ||
                check(N2N2.binary_coefficient_deriv_n(T,P,n),dtruth_dn,tol,"Massman binary coefficient derived to n") ||
                return_flag;

  return return_flag;
}


int main()
{

  return (tester<float>()  || 
          tester<double>() || 
          tester<long double>());
}
