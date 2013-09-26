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

//Planet
#include "planet/atmosphere.h"
#include "planet/planet_constants.h"

//C++
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>


template<typename Scalar>
int check_test(Scalar theory, Scalar cal, const std::string &words)
{
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100.;
  if(std::abs((theory-cal)/theory) < tol)return 0;
  std::cout << "failed test:\n"
            << "theory: " << theory
            << "\ncalculated: " << cal
            << "\ndifference: " << std::abs((theory-cal)/cal)
            << "\ntolerance: " << tol << std::endl;
  std::cout << words << std::endl;
  return 1;
}

template<typename Scalar, typename VectorScalar = std::vector<Scalar> >
void read_hv_flux(VectorScalar &lambda, VectorScalar &phy1AU, const std::string &file)
{
  std::string line;
  std::ifstream flux_1AU(file);
  getline(flux_1AU,line);
  while(!flux_1AU.eof())
  {
     Scalar wv,ir,dirr;
     flux_1AU >> wv >> ir >> dirr;
     if(!lambda.empty() && wv == lambda.back())continue;
     lambda.push_back(wv);//nm
     phy1AU.push_back(ir);//W/m2/nm
  }
  flux_1AU.close();
  return;
}

template<typename Scalar, typename VectorScalar = std::vector<Scalar> >
void read_crossSection(VectorScalar & lambda, VectorScalar & sigma, const std::string &file, unsigned int nbr)
{
  lambda.clear();
  sigma.clear();
  std::string line;
  std::ifstream sig_f(file);
  getline(sig_f,line);
  while(!sig_f.eof())
  {
     Scalar wv,sigt,sigbr;
     sig_f >> wv >> sigt;
     for(unsigned int i = 0; i < nbr; i++)sig_f >> sigbr;
     lambda.push_back(wv/10.);//A -> nm
     sigma.push_back(sigt*10.);//cm-2/A -> m-2/nm
  }
  sig_f.close();
  return;
}

template<typename Scalar, typename VectorScalar = std::vector<Scalar> >
void  read_temperature(VectorScalar &T0, VectorScalar &Tz, const std::string &file)
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

template <typename Scalar>
int test()
{
  std::vector<std::string> neutrals;
  std::vector<std::string> ions;
  neutrals.push_back("N2");
  neutrals.push_back("CH4");
  ions.push_back("N2+");
  ions.push_back("e");

//photons
  std::vector<Scalar> lambda,phy1AU;
  read_hv_flux<Scalar>(lambda,phy1AU,"./input/hv_SSI.dat");
  Scalar chi(100.);
  Planet::Chapman<Scalar> chapman(chi);
  Planet::PhotonFlux<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > photon(chapman);
  photon.set_photon_flux_top_atmosphere(lambda,phy1AU,Planet::Constants::Saturn::d_Sun<Scalar>());

//densities
  std::vector<Scalar> molar_frac = {0.96,0.04,0.,0.};
  Scalar dens_tot(1e12);

//diffusion
  Planet::BinaryDiffusion<Scalar> N2N2(   Antioch::Species::N2,  Antioch::Species::N2 , 5.09e+16,0.81, Planet::DiffusionType::Massman);
  Planet::BinaryDiffusion<Scalar> N2CH4(  Antioch::Species::N2,  Antioch::Species::CH4, 7.34e+16,0.75, Planet::DiffusionType::Massman);
  Planet::BinaryDiffusion<Scalar> CH4CH4( Antioch::Species::CH4, Antioch::Species::CH4, 5.73e+16,0.50, Planet::DiffusionType::Massman);

//atmosphere
  Planet::Atmosphere<Scalar, std::vector<Scalar>, std::vector<std::vector<Scalar> > > atm(neutrals,ions,photon);
  atm.make_altitude_grid(600.,1400.,10.);

//temperature
  std::vector<Scalar> T0,Tz;
  read_temperature<Scalar>(T0,Tz,"input/temperature.dat");
  atm.set_temperature(T0,Tz);

//cross-section
  std::vector<Scalar> sigma;
  read_crossSection<Scalar>(lambda,sigma,"./input/N2_hv_cross-sections.dat",3);
  atm.add_photoabsorption("N2",lambda,sigma);

  read_crossSection<Scalar>(lambda,sigma,"./input/CH4_hv_cross-sections.dat",9);
  atm.add_photoabsorption("N2",lambda,sigma);
  atm.add_photoabsorption("CH4",lambda,sigma);

//init compo
  atm.init_composition(molar_frac,dens_tot);

//diffusion
  atm.diffusion().set_binary_coefficient(0,0,N2N2);
  atm.diffusion().set_binary_coefficient(0,1,N2CH4);
  atm.diffusion().set_binary_coefficient(1,1,CH4CH4);
  atm.make_diffusion();

//thermal coefficient
  Scalar K0(4.3e6L);
  atm.set_K0(K0);
  atm.make_thermal_coefficient();

  std::vector<Scalar> T = atm.temperature_top_to_bottom();

  unsigned int i(0);
  Scalar MN(14.008L), MC(12.011), MH(1.008L);
  Scalar MN2 = 2.L*MN , MCH4 = MC + 4.L*MH;
  int return_flag(0);
  std::ofstream out("data.dat");
  out << "z M g H a n q" << std::endl;
  for(Scalar z = 1400.; z >= 600.; z -= 10. )
  {
    Scalar M_the = molar_frac[0] * MN2 + molar_frac[1] * MCH4;
    if(check_test<Scalar>(M_the,atm.M(z),"Atmosphere mean atomic mass"))return_flag = 1;

    Scalar g_the = Planet::Constants::Universal::G<Scalar>() * Planet::Constants::Titan::mass<Scalar>()/
                   (1e6L * (Planet::Constants::Titan::radius<Scalar>() + z) * (Planet::Constants::Titan::radius<Scalar>() + z));
    if(check_test<Scalar>(g_the,atm.g(z),"Titan gravitationnal fields"))return_flag = 1;

    Scalar n_tot_the = dens_tot * std::exp(-(z - 600.) / 
                       ( (z + Planet::Constants::Titan::radius<Scalar>())    *
                         (600. + Planet::Constants::Titan::radius<Scalar>()) * 1e3L * //to m
                         ( (Planet::Constants::Universal::kb<Scalar>() * Antioch::Constants::Avogadro<Scalar>() * T[i]) /
                           (Planet::Constants::Universal::G<Scalar>() * Planet::Constants::Titan::mass<Scalar>() * M_the * 1e-3L) //to kg/mol
                         )
                       ));
    if(check_test<Scalar>(n_tot_the,atm.total_density(z),"Titan's atmosphere's total density"))return_flag = 1;

    Scalar K_theo = K0 * std::sqrt(n_tot_the/dens_tot);
    if(check_test<Scalar>(K_theo,atm.K(z),"Thermal coefficient at altitude z"))return_flag = 1;

    Scalar H_the = Planet::Constants::Universal::kb<Scalar>() * T[i] /
                   (1e-3L * g_the *  M_the / Antioch::Constants::Avogadro<Scalar>()); //kb*T/(g(z) * M/Navo * n)
    if(check_test<Scalar>(H_the,atm.H(z),"Atmosphere scale height"))return_flag = 1;

    Scalar a_the = Scalar(1e3) * (Planet::Constants::Titan::radius<Scalar>() + z) / H_the;
    if(check_test<Scalar>(a_the,atm.a(z),"Titan's atmosphere a"))return_flag = 1;

    out << z << " " << M_the     << " " << g_the << " " 
                    << H_the     << " " << a_the << " " 
                    << n_tot_the << " " << Planet::Constants::Universal::kb<Scalar>() * T[i] / (g_the) //kb*T/(g(z) * M/Navo * n)
                    << " "       << std::endl;

    i++;
  }

  out.close();

  return return_flag;
}

int main()
{

  return (test<float>() ||
          test<double>());// ||
//          test<long double>());
}
