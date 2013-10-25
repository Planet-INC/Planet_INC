//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
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

  Scalar p11(5.09e16),p12(0.81);
  Scalar p21(7.34e16),p22(0.75);
  Scalar p31(1.04e-5),p32(1.76);

  Planet::BinaryDiffusion<Scalar> N2N2(   Antioch::Species::N2,  Antioch::Species::N2 , p11, p12, Planet::DiffusionType::Massman);
  Planet::BinaryDiffusion<Scalar> N2CH4(  Antioch::Species::N2,  Antioch::Species::CH4, p21, p22, Planet::DiffusionType::Wilson);
  Planet::BinaryDiffusion<Scalar> CH4CH4( Antioch::Species::CH4, Antioch::Species::CH4, p31, p32, Planet::DiffusionType::Wakeham);

  Scalar T(1500.),P(1e5);
  Scalar n = P / (T * Antioch::Constants::R_universal<Scalar>()/1000.L);

  Scalar n2n2 = p11 * Planet::Constants::Convention::P_normal<Scalar>() / P * Antioch::ant_pow(T/Planet::Constants::Convention::T_standard<Scalar>(),p12);
  Scalar n2ch4 = p21 * Antioch::ant_pow(T,p22)/n;
  Scalar ch4ch4 = p31 * Antioch::ant_pow(T,p32) * Planet::Constants::Convention::P_normal<Scalar>() / P;

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
