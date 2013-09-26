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

//Planet
#include "planet/chapman.h"

//C++
#include <cmath>
#include <limits>
#include <iomanip>


template<typename Scalar>
int check_test(Scalar theory, Scalar cal, const std::string &words)
{
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100.;
  if(std::abs((theory-cal)/cal) < tol)return 0;
  std::cout << "failed test:\n"
            << "theory: " << theory
            << "\ncalculated: " << cal
            << "\ndifference: " << std::abs((theory-cal)/cal)
            << "\ntolerance: " << tol << std::endl;
  std::cout << words << std::endl;
  return 1;
}

template<typename Scalar>
Scalar erf(Scalar x)
{
  if(x < 0.)x = -x;
  return Scalar(1.L) - (Scalar(0.254829592)   * (Scalar(1.L)/(Scalar(1.L)         + Scalar(0.3275911L) * x)) +
                        Scalar(-0.284496736L) * std::pow(Scalar(1.L)/(Scalar(1.L) + Scalar(0.3275911L) * x),2) + 
                        Scalar(1.421413741L)  * std::pow(Scalar(1.L)/(Scalar(1.L) + Scalar(0.3275911L) * x),3) + 
                        Scalar(-1.453152027L) * std::pow(Scalar(1.L)/(Scalar(1.L) + Scalar(0.3275911L) * x),4) + 
                        Scalar(1.061405429L)  * std::pow(Scalar(1.L)/(Scalar(1.L) + Scalar(0.3275911L) * x),5)
                       ) * std::exp(-std::pow(x,2));

}

template <typename Scalar>
int test()
{
  Planet::Chapman<Scalar> chap;

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100.;
  int return_flag(0);
  Scalar x(70e3L);// (R + z) / H
  for(Scalar chi = 30.; chi <= 75.; chi += 5.)// Chap = 1/cos(chi)
  {
    Scalar rchi = Planet::Constants::pi<Scalar>()/180.L * chi;
    chap.set_chi(chi);
    Scalar chap_one = chap.chapman();
    Scalar chap_two = chap.chapman(x);
    Scalar chap_theo = 1.L/std::cos(rchi);
    if(std::abs((chap_theo - chap_one)/chap_theo) > tol)
    {
      std::cout << "failed test: Chapman low angle\n"
                << "chi: " << chi
                << std::setprecision(25)
                << "\ntheory: " << chap_theo
                << "\ncalculated: " << chap_one
                << "\ndifference: " << std::abs((chap_theo-chap_one)/chap_theo)
                << "\ntolerance: " << tol << std::endl;
      return_flag = 1;
    }
    if(std::abs((chap_theo - chap_two)/chap_theo) > tol)
    {
      std::cout << "failed test: Chapman low angle\n"
                << "chi: " << chi
                << std::setprecision(25)
                << "theory: " << chap_theo
                << "\ncalculated: " << chap_two
                << "\ndifference: " << std::abs((chap_theo-chap_two)/chap_two)
                << "\ntolerance: " << tol << std::endl;
      return_flag = 1;
    }
  }

  for(Scalar chi = 80.L; chi <= 89.L; chi += 2.L)// Chap = sqrt(pi*x/2)...
  {
    Scalar rchi = Planet::Constants::pi<Scalar>()/180.L * chi;
    chap.set_chi(chi);
    Scalar chap_cal = chap.chapman(x);
    Scalar chap_theo = std::sqrt(Planet::Constants::pi<Scalar>() * x /2.L)        *
                       (1.L - erf(std::sqrt(x/2.L) * std::abs(std::cos(rchi)) ) ) *
                       std::exp(x/2.L * std::pow(std::cos(rchi),2));

    if(std::abs((chap_theo - chap_cal)/chap_theo) > tol)
    {
      std::cout << "failed test: Chapman med angle\n"
                << "chi: " << chi
                << std::setprecision(25)
                << "\ntheory: " << chap_theo
                << "\ncalculated: " << chap_cal
                << "\ndifference: " << std::abs((chap_theo-chap_cal)/chap_theo)
                << "\ntolerance: " << tol << std::endl;
      return_flag = 1;
    }
  }

  for(Scalar chi = 92.; chi <= 180.; chi += 2.)// Chap = sqrt(2*pi*x)...
  {
    Scalar rchi = Planet::Constants::pi<Scalar>()/180.L * chi;
    chap.set_chi(chi);
    Scalar chap_cal = chap.chapman(x);
    Scalar chap_theo = std::sqrt(Planet::Constants::pi<Scalar>() * x * 2.L) * 
                (
                   std::sqrt(std::sin(rchi)) * std::exp(x * (1.L - std::sin(rchi)) ) -
                   0.5L * std::exp(
                                     x/2.L * std::cos(rchi) * std::cos(rchi) * 
                                     (1.L - erf(std::sqrt(x/2.L) * std::abs(std::cos(rchi))))
                                  )
                );
    if(std::abs((chap_theo - chap_cal)/chap_theo) > tol)
    {
      std::cout << "failed test: Chapman high angle\n"
                << "chi: " << chi
                << std::scientific << std::setprecision(25)
                << "\ntheory: " << chap_theo
                << "\ncalculated: " << chap_cal
                << "\ndifference: " << std::abs((chap_theo-chap_cal)/chap_theo)
                << "\ntolerance: " << tol << std::endl;
      return_flag = 1;
    }
  }

  return return_flag;
}

int main()
{

  return (test<float>() ||
          test<double>());/* ||
          test<long double>());*/
}
