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
#include "antioch/physical_constants.h"
#include "antioch/sigma_bin_converter.h"
#include "antioch/vector_utils.h"

//Planet
#include "planet/photon_evaluator.h"
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
  Scalar coeff = (std::numeric_limits<Scalar>::epsilon() < 1e-17)?5e3:100.;
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * coeff;
  if(std::abs((theory-cal)/theory) < tol)return 0;
  std::cout << std::scientific << std::setprecision(20)
            << "failed test: " << words << "\n"
            << "theory: " << theory
            << "\ncalculated: " << cal
            << "\ndifference: " << std::abs((theory-cal)/cal)
            << "\ntolerance: " << tol << std::endl;
  return 1;
}

template<typename VectorScalar>
void linear_interpolation(const VectorScalar &temp0, const VectorScalar &alt0,
                          const VectorScalar &alt1, VectorScalar &temp1)
{
  unsigned int j(0);
  typename Antioch::value_type<VectorScalar>::type a;
  typename Antioch::value_type<VectorScalar>::type b;
  temp1.resize(alt1.size());
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

template<typename Scalar, typename VectorScalar = std::vector<Scalar> >
void read_crossSection(const std::string &file, unsigned int nbr, VectorScalar &lambda, VectorScalar &sigma)
{
  std::string line;
  std::ifstream sig_f(file);
  getline(sig_f,line);
  while(!sig_f.eof())
  {
     Scalar wv,sigt,sigbr;
     sig_f >> wv >> sigt;
     for(unsigned int i = 0; i < nbr; i++)sig_f >> sigbr;
     lambda.push_back(wv);//A
     sigma.push_back(sigt);//cm-2/A
  }
  sig_f.close();

  return;
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
     lambda.push_back(wv * 10.L);//nm -> A
     phy1AU.push_back(ir * 1e3L * (wv*1e-9L) / (Antioch::Constants::Planck_constant<Scalar>() * 
                                        Antioch::Constants::light_celerity<Scalar>()));//W/m2/nm -> J/s/cm2/A -> s-1/cm-2/A
  }
  flux_1AU.close();
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

template<typename Scalar, typename VectorScalar, typename MatrixScalar>
void calculate_densities(MatrixScalar &densities, const Scalar &tot_dens, const VectorScalar &molar_frac, 
                        const Scalar &zmin,const Scalar &zmax,const Scalar &zstep, 
                        const VectorScalar &T, const VectorScalar &mm)
{
   unsigned int iz(0);
   Scalar Mm;
   Antioch::set_zero(Mm);
   for(unsigned int s = 0; s < molar_frac.size(); s++)
   {
      Mm += molar_frac[s] * mm[s];
   }
   Mm *= 1e-3;//to kg
   densities.clear();
   densities.resize(molar_frac.size());
   for(Scalar z = zmin; z <= zmax; z += zstep)
   {
      for(unsigned int s = 0; s < molar_frac.size(); s++)
      {
        densities[s].push_back(molar_frac[s] * barometry(zmin,z,T[iz],Mm,tot_dens));
      }
      iz++;
   }

   return;
}

template<typename Scalar, typename VectorScalar, typename MatrixScalar>
void calculate_tau(MatrixScalar &opacity, const Planet::Chapman<Scalar> &chapman, 
                   const std::vector<VectorScalar*> &cs, const VectorScalar &lambda_ref,
                   const MatrixScalar &sum_dens, const VectorScalar &molar_frac,
                   const VectorScalar &mm, const VectorScalar &T,
                   const Scalar &zmin, const Scalar &zmax, const Scalar &zstep)
{
  Scalar Mm;
  Antioch::set_zero(Mm);
  for(unsigned int s = 0; s < molar_frac.size(); s++)
  {
     Mm += molar_frac[s] * mm[s];
  }
  Mm *= 1e-3; //to kg
  opacity.resize(T.size());
  Antioch::SigmaBinConverter<VectorScalar> bin_converter;

  unsigned int iz(0);
  for(Scalar z = zmin; z <= zmax; z += zstep)
  {
    Scalar g = Planet::Constants::g(Planet::Constants::Titan::radius<Scalar>(),z,Planet::Constants::Titan::mass<Scalar>());
    Scalar H = Planet::Constants::Universal::kb<Scalar>() * T[iz] / (g * Mm ) * Antioch::Constants::Avogadro<Scalar>(); //to m-3
    Scalar a = (Planet::Constants::Titan::radius<Scalar>() + z) / H; //to m

    opacity[iz].resize(lambda_ref.size(),0.L);
    for(unsigned int s = 0; s < sum_dens.size(); s++)
    {
      VectorScalar sigma_process;
      bin_converter.y_on_custom_grid(*(cs[2*s]),*(cs[2*s+1]),lambda_ref,sigma_process);
      for(unsigned int il = 0; il < lambda_ref.size(); il++)
      {
        opacity[iz][il] += sigma_process[il] * sum_dens[s][iz];
      }
    }
    for(unsigned int il = 0; il < lambda_ref.size(); il++)
    {
      opacity[iz][il] *= chapman(a) * zstep * 1e5;
    }
    iz++;
  }

}


template <typename Scalar>
int tester(const std::string &input_T, const std::string &input_hv, 
           const std::string &input_N2, const std::string &input_CH4)
{
//description
  std::vector<std::string> neutrals;
  std::vector<std::string> ions;
  neutrals.push_back("N2");
  neutrals.push_back("CH4");
//ionic system contains neutral system
  ions = neutrals;
  ions.push_back("N2+");
  Scalar MN(14.008L), MC(12.011), MH(1.008L);
  Scalar MN2 = 2.L*MN , MCH4 = MC + 4.L*MH;
  std::vector<Scalar> Mm;
  Mm.push_back(MN2);
  Mm.push_back(MCH4);

//densities
  std::vector<Scalar> molar_frac;
  molar_frac.push_back(0.96L);
  molar_frac.push_back(0.04L);
  molar_frac.push_back(0.L);
  Scalar dens_tot(1e12L);

//hard sphere radius
  std::vector<Scalar> hard_sphere_radius;
  hard_sphere_radius.push_back(2.0675e-8L * 1e-2L); //N2  in cm -> m
  hard_sphere_radius.push_back(2.3482e-8L * 1e-2L); //CH4 in cm -> m

//zenith angle
  Scalar chi(120);

//photon flux
  std::vector<Scalar> lambda_hv,phy1AU;
  read_hv_flux<Scalar>(lambda_hv,phy1AU,input_hv);

////cross-section
  std::vector<Scalar> lambda_N2,sigma_N2;
  std::vector<Scalar> lambda_CH4, sigma_CH4;
  read_crossSection<Scalar>(input_N2,3,lambda_N2,sigma_N2);
  read_crossSection<Scalar>(input_CH4,9,lambda_CH4,sigma_CH4);

//altitudes
  Scalar zmin(600.),zmax(1400.),zstep(10.);

//temperature
  std::vector<Scalar> T0,Tz;
  read_temperature<Scalar>(T0,Tz,input_T);
  std::vector<Scalar> neutral_temperature;

/************************
 * first level
 ************************/

//altitude
  Planet::Altitude<Scalar,std::vector<Scalar> > altitude(zmin,zmax,zstep);

//neutrals
  Antioch::ChemicalMixture<Scalar> neutral_species(neutrals); 

//ions
  Antioch::ChemicalMixture<Scalar> ionic_species(ions); 

//chapman
  Planet::Chapman<Scalar> chapman(chi);

//binary diffusion
//not needed

/************************
 * second level
 ************************/

//temperature
  linear_interpolation(T0,Tz,altitude.altitudes(),neutral_temperature);
  Planet::AtmosphericTemperature<Scalar, std::vector<Scalar> > temperature(neutral_temperature, neutral_temperature, altitude);

//photon opacity
  Planet::PhotonOpacity<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > tau(altitude,chapman);

//reaction sets
//not needed

/************************
 * third level
 ************************/

//atmospheric mixture
  Planet::AtmosphericMixture<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > composition(neutral_species, ionic_species, altitude, temperature);
  composition.init_composition(molar_frac,dens_tot);
  composition.set_hard_sphere_radius(hard_sphere_radius);
  composition.initialize();

//kinetics evaluators
//not needed

/************************
 * fourth level
 ************************/

//photon evaluator
  Planet::PhotonEvaluator<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > photon(altitude,tau,composition);
  photon.add_cross_section(lambda_N2,  sigma_N2,  Antioch::Species::N2);
  photon.add_cross_section(lambda_CH4, sigma_CH4, Antioch::Species::CH4);
  photon.set_photon_flux_at_top(lambda_hv, phy1AU, Planet::Constants::Saturn::d_Sun<Scalar>());
  photon.update_photon_flux();

//molecular diffusion
//not needed

//eddy diffusion
//not needed

/************************
 * checks
 ************************/

  molar_frac.pop_back();//get the ion outta here

  std::vector<std::vector<Scalar> > densities;
  std::vector<Scalar> mm;
  mm.push_back(MN2);
  mm.push_back(MCH4);
  calculate_densities(densities,dens_tot,molar_frac,zmin,zmax,zstep,neutral_temperature,mm);

  for(Scalar z = zmax-zstep; z >= zmin; z -= zstep)
  {
     unsigned int iz = altitude.altitudes_map().at(z);
     unsigned int jz = altitude.altitudes_map().at(z+zstep);
     for(unsigned int s = 0; s < densities.size(); s++)
     {
        densities[s][iz] += densities[s][jz];
     }
  }

  std::vector<std::vector<Scalar> > opacity;
  std::vector<std::vector<Scalar>*> cs;
  cs.push_back(&lambda_N2);
  cs.push_back(&sigma_N2);
  cs.push_back(&lambda_CH4);
  cs.push_back(&sigma_CH4);
  calculate_tau(opacity,chapman,cs,lambda_hv,
                densities,molar_frac,mm,
                neutral_temperature,
                zmin,zmax,zstep);


  std::vector<std::vector<Scalar> > phy_theo;
  phy_theo.resize(altitude.altitudes().size());
  for(unsigned int iz = 0; iz < altitude.altitudes().size(); iz++)
  {
    phy_theo[iz].resize(lambda_hv.size());
    for(unsigned int il = 0; il < lambda_hv.size(); il++)  
    {
      phy_theo[iz][il] = phy1AU[il] / (Planet::Constants::Saturn::d_Sun<Scalar>() * Planet::Constants::Saturn::d_Sun<Scalar>()) 
                        * Antioch::ant_exp(-opacity[iz][il]);
    }
  }


  int return_flag(0);

  for(unsigned int iz = 0; iz < altitude.altitudes().size(); iz++)
  {
    for(unsigned int il = 0; il < lambda_hv.size(); il++)
    {
        return_flag = return_flag ||
                      check_test(phy_theo[iz][il], photon.photon_flux()[iz].flux()[il], "phy at altitude and wavelength");
        if(return_flag)
        {
          std::cout << altitude.altitudes()[iz] << " " << lambda_hv[il] << std::endl;
          return return_flag;
        }
    }
  }

  return return_flag;
}

int main(int argc, char** argv)
{
  // Check command line count.
  if( argc < 5 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify reaction set XML input file." << std::endl;
      antioch_error();
    }

  return (tester<float>(std::string(argv[1]),std::string(argv[2]),std::string(argv[3]),std::string(argv[4])) ||
          tester<double>(std::string(argv[1]),std::string(argv[2]),std::string(argv[3]),std::string(argv[4])) ||
          tester<long double>(std::string(argv[1]),std::string(argv[2]),std::string(argv[3]),std::string(argv[4])));
}
