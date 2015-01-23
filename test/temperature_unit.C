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
#include "antioch/vector_utils_decl.h"
#include "antioch/cmath_shims.h"
#include "antioch/vector_utils.h"
//Planet
#include "planet/atmospheric_temperature.h"
//C++
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

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
Scalar linear_interpolation(const Scalar &x0, const Scalar &x1,
                            const Scalar &y0, const Scalar &y1,
                            const Scalar &x)
{

   Scalar a = (y1 - y0)/(x1 - x0);
   Scalar b = y0 - a * x0;
   return (a * x + b);
}

template<typename Scalar>
Scalar dlinear_interpolation(const Scalar &x0, const Scalar &x1,
                             const Scalar &y0, const Scalar &y1)
{

   return (y1 - y0)/(x1 - x0);
}

template<typename Scalar>
Scalar electron_temperature(const Scalar &alt)
{
  if(alt < 900.L)
  {
     return 180.L;
  }else if(alt < 1400.L)
  {
     return 180.L + (alt - 900.L)/500.L * (1150.L - 180.L);
  }else
  {
     return 1150.L + (alt - 1400.L)/10.L;
  }

  return 0;
}

template<typename Scalar, typename VectorScalar = std::vector<Scalar> >
void read_temperature(VectorScalar &T0, VectorScalar &Tz, const std::string &file)
{
  T0.clear();
  Tz.clear();
  std::string line;
  std::ifstream temp(file);
  getline(temp,line);
  while(!temp.eof())
  {
     Scalar t(0.),tz(0.),dt(0.),dtz(0.);
     temp >> t >> tz >> dt >> dtz;
     T0.push_back(t);
     Tz.push_back(tz);
  }
  temp.close();
  return;
}

template<typename Scalar>
int tester(const std::string & input_file)
{
//temperature
  std::vector<Scalar> T0,Tz;
  read_temperature<Scalar>(T0,Tz,input_file);
  Planet::AtmosphericTemperature<Scalar,std::vector<Scalar> > temperature(T0,T0,Tz,Tz); //neutral, ionic, altitude

  int return_flag(0);
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100.;

  for(unsigned int iz = 0; iz < T0.size() - 1; iz++)
  {
      Scalar z = (T0[iz] + T0[iz + 1]) / Scalar(2.L);
      Scalar neu_temp  = linear_interpolation(Tz[iz],Tz[iz+1],T0[iz],T0[iz+1],z);
      Scalar neu_dtemp = dlinear_interpolation(Tz[iz],Tz[iz+1],T0[iz],T0[iz+1]);
      Scalar e_temp    = electron_temperature(z);
      return_flag = return_flag && 
                    check(temperature.neutral_temperature(z)   ,neu_temp,tol,"neutral temperature") &&
                    check(temperature.dneutral_temperature_dz(z)   ,neu_dtemp,tol,"neutral temperature differentiate") &&
                    check(temperature.ionic_temperature(z)     ,neu_temp,tol,"ionic temperature")   &&
                    check(temperature.electronic_temperature(z),e_temp  ,tol,"electron temperature");
  }

  return return_flag;
}


int main(int argc, char** argv)
{
  // Check command line count.
  if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify input file." << std::endl;
      antioch_error();
    }

  return (tester<float>(std::string(argv[1])) ||
          tester<double>(std::string(argv[1])) ||
          tester<long double>(std::string(argv[1])));
}
