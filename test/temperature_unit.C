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

template<typename VectorScalar>
void linear_interpolation(const VectorScalar &temp0, const VectorScalar &alt0,
                          const VectorScalar &alt1, VectorScalar &temp1)
{
  unsigned int j(0);
  typename Antioch::value_type<VectorScalar>::type a;
  typename Antioch::value_type<VectorScalar>::type b;

  temp1.resize(alt1.size(),0.L);
  for(unsigned int iz = 0; iz < alt1.size(); iz++)
  {
     while(alt0[j] < alt1[iz])
     {
        j++;
        if(!(j < alt0.size()))break;
     }
     if(j == 0)
     {
        Antioch::set_zero(a);
        b = temp0[j];
     }else if(j < alt0.size() - 1)
     {
        a = (temp0[j] - temp0[j-1])/(alt0[j] - alt0[j-1]);
        b = temp0[j] - a * alt0[j];
     }else
     {
        Antioch::set_zero(a);
        b = temp0.back();
     }
     temp1[iz] = a * alt1[iz] + b;
  }
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
     Scalar t,tz,dt,dtz;
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
  Planet::Altitude<Scalar,std::vector<Scalar>> altitude(600,1400,10);
//temperature
  std::vector<Scalar> T0,Tz;
  read_temperature<Scalar>(T0,Tz,input_file);
  std::vector<Scalar> neutral_temperature;
  linear_interpolation(T0,Tz,altitude.altitudes(),neutral_temperature);
  Planet::AtmosphericTemperature<Scalar,std::vector<Scalar> > temperature(neutral_temperature,neutral_temperature,altitude);

  int return_flag(0);
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100.;

  for(unsigned int iz = 0; iz < altitude.altitudes().size(); iz++)
  {
      Scalar neu_temp(neutral_temperature[iz]);
      Scalar e_temp = electron_temperature(altitude.altitudes()[iz]);
      return_flag = return_flag && 
                    check(temperature.neutral_temperature()[iz],neu_temp,tol,"neutral temperature") &&
                    check(temperature.ionic_temperature()[iz],neu_temp,tol,"ionic temperature")     &&
                    check(temperature.electronic_temperature()[iz],e_temp,tol,"electron temperature");
  }

  return return_flag;
}


int main(int argc, char** argv)
{
  // Check command line count.
  if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify reaction set XML input file." << std::endl;
      antioch_error();
    }

  return (tester<float>(std::string(argv[1])) ||
          tester<double>(std::string(argv[1])) ||
          tester<long double>(std::string(argv[1])));
}
