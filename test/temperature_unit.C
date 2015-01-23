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

template<typename Scalar>
Scalar delectron_temperature(const Scalar &alt)
{
  if(alt < 900.L)
  {
     return 0.L;
  }else if(alt < 1400.L)
  {
     return 2e-3L;
  }else
  {
     return 0.1L;
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
     Scalar t(0.),tz(0.);
     temp >> t >> tz;
     if(temp.good())
     {
       T0.push_back(t);
       Tz.push_back(tz);
     }
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

  gsl_interp_accel * acc = gsl_interp_accel_alloc();
  gsl_spline       * spline = gsl_spline_alloc(gsl_interp_cspline, T0.size());

  // GLS takes only double, raaaaahhhh
  std::vector<double> gsl_x_point(T0.size(),0);
  std::vector<double> gsl_y_point(T0.size(),0);
  for(unsigned int i = 0; i < T0.size(); i++)
  {
    gsl_x_point[i] = (const double)Tz[i];
    gsl_y_point[i] = (const double)T0[i];
  }

  const double * x = &gsl_x_point[0];
  const double * y = &gsl_y_point[0];

  gsl_spline_init(spline, x, y, T0.size());

  Planet::AtmosphericTemperature<Scalar,std::vector<Scalar> > temperature(Tz,T0); //neutral, ionic, altitude

  int return_flag(0);
  const Scalar tol = (std::numeric_limits<Scalar>::epsilon() * 100. < 6e-17)?6e-17:
                      std::numeric_limits<Scalar>::epsilon() * 100.;

  for(Scalar z = 600.; z < 1401.; z += 10)
  {
      Scalar neu_temp  = gsl_spline_eval(spline,z,acc);
      Scalar neu_dtemp = gsl_spline_eval_deriv(spline,z,acc);
      Scalar ion_temp  = gsl_spline_eval(spline,z,acc);
      Scalar ion_dtemp = gsl_spline_eval_deriv(spline,z,acc);
      Scalar e_temp    = electron_temperature(z);
      Scalar de_dtemp  = delectron_temperature(z);

      return_flag = check(temperature.neutral_temperature(z),        neu_temp,  tol, "neutral temperature")               ||
                    check(temperature.dneutral_temperature_dz(z),    neu_dtemp, tol, "neutral temperature differentiate") ||
                    check(temperature.ionic_temperature(z),          ion_temp,  tol, "ionic temperature")                 ||
                    check(temperature.dionic_temperature_dz(z),      ion_dtemp, tol, "ionic temperature differentiate")   ||
                    check(temperature.electronic_temperature(z),     e_temp,    tol, "electron temperature")              ||
                    check(temperature.delectronic_temperature_dz(z), de_dtemp,  tol, "electron temperature differentiate") ||
                    return_flag;
  }

  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);

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
