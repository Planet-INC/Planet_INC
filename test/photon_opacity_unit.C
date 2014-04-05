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
#include "antioch/sigma_bin_converter.h"
#include "antioch/vector_utils.h"
#include "antioch/physical_constants.h"

//Planet
#include "planet/photon_opacity.h"
#include "planet/planet_constants.h"
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
  if(Antioch::ant_abs(test - ref)/ref > tol && Antioch::ant_abs(test - ref) > tol)
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

template<typename Scalar, typename VectorScalar = std::vector<Scalar> >
void read_crossSection(VectorScalar & lambda, const std::string &file, VectorScalar &sigma)
{
  std::string line;
  std::ifstream sig_f(file);
  getline(sig_f,line);
  while(!sig_f.eof())
  {
     Scalar wv(-1),sigt;
     sig_f >> wv >> sigt;
     if(!getline(sig_f,line))break;
     if(wv < 0.)break;
     lambda.push_back(wv/10.);//A -> nm
     sigma.push_back(sigt*10.);//cm-2/A -> m-2/nm
  }
  sig_f.close();

  return;
}

template<typename Scalar>
Scalar barometry(const Scalar &zmin, const Scalar &z, const Scalar &T, const Scalar &Mm, const Scalar &botdens)
{
   return botdens * Antioch::ant_exp(-(z - zmin)/((Planet::Constants::Titan::radius<Scalar>() + z) * (Planet::Constants::Titan::radius<Scalar>() + zmin) * 1e3 *
                                             Antioch::Constants::Avogadro<Scalar>() * Planet::Constants::Universal::kb<Scalar>() * T / 
                                                        (Planet::Constants::Universal::G<Scalar>() * Planet::Constants::Titan::mass<Scalar>() * Mm))
                              );
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

template<typename Scalar, typename VectorScalar>
void densities(const Scalar &z, const Scalar &zmin, const Scalar &zmax, const Scalar &dens_tot,
               const Planet::AtmosphericTemperature<Scalar, std::vector<Scalar> > &temperature,
               const Scalar &Mmean, const VectorScalar &molar_frac, VectorScalar &sum_dens)
{
  Scalar zstep(1.L);
  sum_dens.resize(molar_frac.size(),0.L);
  for(Scalar ztmp = zmax; ztmp > z; ztmp -= zstep)
  {
     Scalar T     = temperature.neutral_temperature(z);
     Scalar nTot  = barometry(zmin,ztmp,T,Mmean,dens_tot);
     for(unsigned int s = 0; s < molar_frac.size(); s++)
     {
       sum_dens[s] += nTot * molar_frac[s] * zstep;
     }
  }
}

template<typename Scalar>
Scalar a(const Scalar &T, const Scalar &Mmean, const Scalar &z)
{
  return Planet::Constants::g<Scalar>(Planet::Constants::Titan::radius<Scalar>(),z,Planet::Constants::Titan::mass<Scalar>()) * Mmean * 
            (Planet::Constants::Titan::radius<Scalar>() + z) * Scalar(1e3) / (Antioch::Constants::Avogadro<Scalar>() * Planet::Constants::Universal::kb<Scalar>() * T);
}

template<typename Scalar>
int tester(const std::string &input_T, const std::string &input_N2, const std::string &input_CH4)
{
  Scalar chi(120);
//
  Scalar zmin(600.);
  Scalar zmax(1400.);
  Scalar zstep(10.);
//
  Scalar MN(14.008L), MC(12.011), MH(1.008L);
  Scalar MN2 = 2.L*MN , MCH4 = MC + 4.L*MH;
  std::vector<Scalar> Mm;
  Mm.push_back(MN2);
  Mm.push_back(MCH4);

  std::vector<Scalar> molar_frac;
  molar_frac.push_back(0.96L);
  molar_frac.push_back(0.04L);
  molar_frac.push_back(0.L);
  Scalar dens_tot(1e12L);

  Scalar Mmean(0.L);
  for(unsigned int s = 0; s < Mm.size(); s++)
  {
     Mmean += Mm[s] * molar_frac[s];
  }
  Mmean *= 1e-3; // to kg

////cross-section
  std::vector<Scalar> lambda;
  for(Scalar l = 0.1; l < 300.; l += 1.)
  {
      lambda.push_back(l);
  }
  std::vector<std::vector<Scalar> > sigmas;
  std::vector<std::vector<Scalar> > lambdas;
  sigmas.resize(2);
  lambdas.resize(2);
  read_crossSection<Scalar>(lambdas[0],input_N2,sigmas[0]);
  read_crossSection<Scalar>(lambdas[1],input_CH4,sigmas[1]);
  std::vector<std::vector<Scalar> > sigma_ref;
  sigma_ref.resize(2);
  Antioch::SigmaBinConverter<std::vector<Scalar> > binconv;
  for(unsigned int i = 0; i < 2; i++)
  {
    binconv.y_on_custom_grid(lambdas[i],sigmas[i],lambda,sigma_ref[i]);
  }

//////
  Planet::Chapman<Scalar> chapman(chi);
/////
  Planet::PhotonOpacity<Scalar,std::vector<Scalar> > tau(chapman);
  tau.add_cross_section(lambdas[0],sigmas[0],Antioch::Species::N2, 0);
  tau.add_cross_section(lambdas[1],sigmas[1],Antioch::Species::CH4, 1);
  tau.update_cross_section(lambda);

//temperature
  std::vector<Scalar> T0,Tz;
  read_temperature<Scalar>(T0,Tz,input_T);
  Planet::AtmosphericTemperature<Scalar, std::vector<Scalar> > temperature(T0, T0, Tz, Tz);


////////////////////:
  molar_frac.pop_back();

  int return_flag(0);
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100.;
  
  for(Scalar z = zmin; z <= zmax; z += zstep)
  {
     Scalar T = temperature.neutral_temperature(z);

     std::vector<Scalar> sum_dens;
     densities(z,zmin,zmax,dens_tot,temperature,Mmean,molar_frac,sum_dens);

     Scalar x = a(T,Mmean,z);
     std::vector<Scalar> tau_cal;
     tau.compute_tau(x,sum_dens,tau_cal);

     for(unsigned int il = 0; il < lambda.size(); il++)
     {
        Scalar tau_exact(0.L);
        for(unsigned int s = 0; s < 2; s++)
        {
           tau_exact += sigma_ref[s][il] * sum_dens[s]; // cm-3.km
           
           return_flag = check(tau.absorbing_species_cs()[s].cross_section_on_custom_grid()[il],
                               sigma_ref[s][il],tol,"sigma ref of species at altitude and wavelength") ||
                         return_flag;
        }
        tau_exact *= chapman(x) * 1e5; //cm-1.km to no unit
        return_flag = check(tau_cal[il],tau_exact,tol,"tau at altitude and wavelength") ||
                      return_flag;
                      
     }
  }

  return return_flag;
}


int main(int argc, char** argv)
{
  // Check command line count.
  if( argc < 4 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify input file." << std::endl;
      antioch_error();
    }
  return (tester<float>(std::string(argv[1]), std::string(argv[2]), std::string(argv[3])) ||
          tester<double>(std::string(argv[1]), std::string(argv[2]), std::string(argv[3])) ||
          tester<long double>(std::string(argv[1]), std::string(argv[2]), std::string(argv[3])));
}
