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
#include "planet/diffusion_evaluator.h"
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
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 6000.L;
  if(std::abs((theory-cal)/theory) < tol)return 0;
  std::cout << std::scientific << std::setprecision(20)
            << "\nfailed test: " << words << "\n"
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
Scalar barometry(const Scalar &zmin, const Scalar &z, const Scalar &T, const Scalar &Mm, const Scalar &botdens)
{
   return botdens * Antioch::ant_exp(-(z - zmin)/((Planet::Constants::Titan::radius<Scalar>() + z) * (Planet::Constants::Titan::radius<Scalar>() + zmin) * 1e3 *
                                             Antioch::Constants::Avogadro<Scalar>() * Planet::Constants::Universal::kb<Scalar>() * T / 
                                                        (Planet::Constants::Universal::G<Scalar>() * Planet::Constants::Titan::mass<Scalar>() * Mm * 1e-3))
                              );
}

template<typename Scalar>
Scalar dbarometry_dz(const Scalar &zmin, const Scalar &z, const Scalar &T, const Scalar &Mm, const Scalar &botdens)
{
   return barometry(zmin, z, T, Mm, botdens) / ((Planet::Constants::Titan::radius<Scalar>() + z) * (Planet::Constants::Titan::radius<Scalar>() + zmin) * 1e6 *
                                             Antioch::Constants::Avogadro<Scalar>() * Planet::Constants::Universal::kb<Scalar>() * T / 
                                                        (Planet::Constants::Universal::G<Scalar>() * Planet::Constants::Titan::mass<Scalar>() * Mm * 1e-3))
          * (Scalar(1.L) - Scalar(2.L) * (z - zmin)/(Planet::Constants::Titan::radius<Scalar>() + zmin));
}

template<typename Scalar, typename VectorScalar>
void calculate_densities(VectorScalar &densities, VectorScalar &dn_dz, const VectorScalar &molar_frac, const Scalar &nTot, const Scalar &dnTot_dz)
{
   densities.resize(molar_frac.size(),0.L);
   dn_dz.resize(molar_frac.size(),0.L);
   for(unsigned int s = 0; s < molar_frac.size(); s++)
   {
     densities[s] = molar_frac[s] * nTot;
     dn_dz[s] = molar_frac[s] * dnTot_dz;
   }

   return;
}

template<typename Scalar>
Scalar binary_coefficient(const Scalar &T, const Scalar &P, const Scalar &D01, const Scalar &beta)
{
   return D01 * Planet::Constants::Convention::P_normal<Scalar>() / P * Antioch::ant_pow(T/Planet::Constants::Convention::T_standard<Scalar>(),beta);
}

template<typename Scalar>
Scalar binary_coefficient(const Scalar &Dii, const Scalar &Mi, const Scalar &Mj)
{
   return (Mj < Mi)?Dii * Antioch::ant_sqrt((Mj/Mi + Scalar(1.L))/Scalar(2.L))
                   :
                    Dii * Antioch::ant_sqrt((Mj/Mi));
}

template<typename Scalar>
Scalar pressure(const Scalar &n, const Scalar &T)
{
   return n * 1e6L * Planet::Constants::Universal::kb<Scalar>() * T;
}

template<typename Scalar>
Scalar scale_height(const Scalar &T, const Scalar &z, const Scalar &Mm)
{
  return Antioch::Constants::R_universal<Scalar>() * T / 
         (Planet::Constants::g<Scalar>(Planet::Constants::Titan::radius<Scalar>(),z,Planet::Constants::Titan::mass<Scalar>()) * Mm);
}


template <typename Scalar>
int tester(const std::string &input_T)
{
//description
  std::vector<std::string> neutrals;
  std::vector<std::string> ions;
  neutrals.push_back("N2");
  neutrals.push_back("CH4");
//ionic system contains neutral system
  ions = neutrals;
  Scalar MN(14.008L), MC(12.011), MH(1.008L);
  Scalar MN2 = 2.L*MN , MCH4 = MC + 4.L*MH;
  std::vector<Scalar> Mm;
  Mm.push_back(MN2);
  Mm.push_back(MCH4);
  std::vector<std::string> medium;
  medium.push_back("N2");
  medium.push_back("CH4");

//densities
  std::vector<Scalar> molar_frac;
  molar_frac.push_back(0.96L);
  molar_frac.push_back(0.04L);
  Scalar dens_tot(1.7e13L);

//zenith angle
//not necessary

//photon flux
//not necessary

////cross-section
//not necessary

//altitudes
  Scalar zmin(600.),zmax(1400.),zstep(10.);

//binary diffusion
  Scalar bCN1(1.04e-5),bCN2(1.76); //cm2.s-1
  Planet::DiffusionType CN_model(Planet::DiffusionType::Wakeham);
  Scalar bCC1(5.73e16),bCC2(0.5); //cm2.s-1
  Planet::DiffusionType CC_model(Planet::DiffusionType::Wilson);
  Scalar bNN1(0.1783),bNN2(1.81); //cm2.s-1
  Planet::DiffusionType NN_model(Planet::DiffusionType::Massman);

//thermal coefficient
  std::vector<Scalar> tc;
  tc.push_back(0.L); //N2
  tc.push_back(0.L); //CH4

//eddy
  Scalar K0(4.3e6L);//cm2.s-1

/************************
 * first level
 ************************/

//neutrals
  Antioch::ChemicalMixture<Scalar> neutral_species(neutrals); 

//ions
  Antioch::ChemicalMixture<Scalar> ionic_species(ions); 

//chapman
//not needed

//binary diffusion
  Planet::BinaryDiffusion<Scalar> N2N2(   0, 0, bNN1, bNN2, NN_model);
  Planet::BinaryDiffusion<Scalar> N2CH4(  0, 1, bCN1, bCN2, CN_model);
  Planet::BinaryDiffusion<Scalar> CH4CH4( 1, 1, bCC1, bCC2, CC_model);
  std::vector<std::vector<Planet::BinaryDiffusion<Scalar> > > bin_diff_coeff;
  bin_diff_coeff.resize(2);
  for(unsigned int n = 0; n < 2; n++)
  {
    bin_diff_coeff[n].resize(2);
  }
  bin_diff_coeff[0][0] = N2N2;
  bin_diff_coeff[0][1] = N2CH4;
  bin_diff_coeff[1][0] = N2CH4;
  bin_diff_coeff[1][1] = CH4CH4;


/************************
 * second level
 ************************/

//temperature
  std::vector<Scalar> T0,Tz;
  read_temperature<Scalar>(T0,Tz,input_T);
  Planet::AtmosphericTemperature<Scalar, std::vector<Scalar> > temperature(Tz, T0);

//photon opacity
//not needed

//reaction sets
//not needed

/************************
 * third level
 ************************/

//atmospheric mixture
  Planet::AtmosphericMixture<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > composition(neutral_species, ionic_species, temperature);
  composition.init_composition(molar_frac,dens_tot,zmin,zmax);
  composition.set_thermal_coefficient(tc);

//kinetics evaluators
//not needed

/************************
 * fourth level
 ************************/

//photon evaluator
//not needed

//molecular diffusion
  Planet::MolecularDiffusionEvaluator<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > molecular_diffusion(bin_diff_coeff,composition,temperature, medium);

//eddy diffusion
  Planet::EddyDiffusionEvaluator<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > eddy_diffusion(composition,K0);

/**************************
 * fifth level
 **************************/

  Planet::DiffusionEvaluator<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > diffusion(molecular_diffusion,eddy_diffusion,composition,temperature);

/************************
 * checks
 ************************/

  Scalar mean_M;
  Antioch::set_zero(mean_M);
  for(unsigned int s = 0; s < molar_frac.size(); s++)
  {
    mean_M += molar_frac[s] * Mm[s];
  }

  std::vector<std::vector<Scalar> > Dij;
  Dij.resize(2);
  for(unsigned int s = 0; s < molar_frac.size(); s++)
  {
    Dij[s].resize(2,0.L);
  }

  std::vector<Scalar> Dtilde;
  Dtilde.resize(molar_frac.size(),0.L);

  int return_flag(0);
  for(Scalar z = zmin; z <= zmax; z += zstep)
  {
     std::stringstream walt;
     walt << z;

     Scalar T        = temperature.neutral_temperature(z);
     Scalar dT_dz    = temperature.dneutral_temperature_dz(z);
     Scalar nTot     = barometry(zmin,z,T,mean_M,dens_tot);
     Scalar dnTot_dz = dbarometry_dz(zmin,z,T,mean_M,dens_tot);
     Scalar P        = pressure(nTot,T);
     Scalar Ha       = scale_height(T,z,mean_M);

     std::vector<Scalar> dns_dz;
     std::vector<Scalar> densities;
     calculate_densities(densities,dns_dz,molar_frac,nTot,dnTot_dz);

//eddy
     Scalar K = K0 * Antioch::ant_sqrt(dens_tot/nTot);

//mol
     Dij[0][0] = binary_coefficient(T,P,bNN1,bNN2); //N2 N2
     Dij[1][1] = binary_coefficient(T,P,bCC1 * Antioch::ant_pow(Planet::Constants::Convention::T_standard<Scalar>(),bCC2 + Scalar(1.L)) 
                                             * Planet::Constants::Universal::kb<Scalar>()
                                             / Planet::Constants::Convention::P_normal<Scalar>(),bCC2 + Scalar(1.L)); //CH4 CH4
     Dij[0][1] = binary_coefficient(T,P,bCN1 * Antioch::ant_pow(Planet::Constants::Convention::T_standard<Scalar>(),bCN2),bCN2); //N2 CH4
     Dij[1][0] = Dij[0][1]; //CH4 N2

     std::vector<Scalar> total_diffusion;
     diffusion.diffusion(densities,dns_dz,z,total_diffusion);

     for(unsigned int s = 0; s < molar_frac.size(); s++)
     {
       Scalar tmp(0.L);
       for(unsigned int medium = 0; medium < 2; medium++)
       {
          if(s == medium)continue;
          tmp += densities[medium]/Dij[medium][s];
       }
       Scalar Ds = (nTot - densities[s]) / tmp;

       Scalar M_diff(0.L);
       Scalar totdens_diff(0.L);
       for(unsigned int j = 0; j < molar_frac.size(); j++)
       {
          if(s == j)continue;
          M_diff += densities[j] * Mm[j];
          totdens_diff += densities[j];
       }
       M_diff /= totdens_diff;

       Scalar Dtilde = Ds / (Scalar(1.L) - molar_frac[s] * 
                            (Scalar(1.L) - Mm[s] / M_diff)
                            );
//
       Scalar Hs = scale_height(T,z,Mm[s]); //km
       Scalar omega_theo = - Dtilde * ( dns_dz[s] /densities[s]
                                      + Scalar(1.L)/Hs 
                                      + dT_dz /T * (Scalar(1.L) + (Scalar(1.L) - molar_frac[s]) * tc[s]))
                           - K      * ( dns_dz[s] /densities[s]
                                      + Scalar(1.L)/Ha
                                      + dT_dz/T);
       omega_theo *= Scalar(1e-10) * densities[s];

       return_flag = check_test(omega_theo,total_diffusion[s],"omega of species " + neutrals[s] + " at altitude " + walt.str()) || return_flag;
                        
     }
  }

  Scalar z(1.00126188021535153894e+03);
  Scalar nN2(8.43005804625429992676e+10);
  Scalar dnN2_dz(-1.48999877435783243179e+09);
  Scalar nCH4(1.20988959053107666969e+09);
  Scalar dnCH4_dz(-2.13845977940878272057e+07);

  std::vector<Scalar> densities(2,0.);
  densities[0] = nN2;
  densities[1] = nCH4;
  std::vector<Scalar> ddensities_dz(2,0.);
  ddensities_dz[0] = dnN2_dz;
  ddensities_dz[1] = dnCH4_dz;

  std::vector<Scalar> A_lib(2,0.),B_lib(2,0.);
  std::vector<std::vector<Scalar> > dA_lib(2,std::vector<Scalar>(2,0.)),dB_lib(2,std::vector<Scalar>(2,0.));

  diffusion.diffusion_and_derivs(densities,ddensities_dz,z,A_lib,B_lib,dA_lib,dB_lib);

  Scalar omega_A_N2(-8.91745150154882304580961160994058011150693141477e-3);
  Scalar omega_A_CH4(-1.104779285772000424689969814217345776876057607747e-2);
  Scalar omega_B_N2(-8.890663229272582684868201945551452637241904819e-5);
  Scalar omega_B_CH4(-8.472520943252638403954699093459571430051180243e-5);

  Scalar domega_A_N2_dCH4(5.468351719402458916575928694020385445e-14);
  Scalar domega_A_CH4_dN2(9.41013244365896810745955626496851621e-14);
  Scalar domega_B_N2_dCH4(9.0325080649590294607336826345205794e-16);
  Scalar domega_B_CH4_dN2(6.3478962472226290094290797587472414e-16);
  Scalar domega_A_N2_dN2(6.90365411299394253965738848931126499e-14);
  Scalar domega_B_N2_dN2(6.8398665016645709138999378530991927e-16);
  Scalar domega_A_CH4_dCH4(6.90365411299394253965738848931126499e-14);
  Scalar domega_B_CH4_dCH4(8.7507876626170793617393231111165119e-16);


  return_flag = check_test(omega_A_N2,A_lib[0], "golden value omega_A N2")  || return_flag;
  return_flag = check_test(omega_B_N2,B_lib[0], "golden value omega_B N2")  || return_flag;
  return_flag = check_test(omega_A_CH4,A_lib[1],"golden value omega_A CH4") || return_flag;
  return_flag = check_test(omega_B_CH4,B_lib[1],"golden value omega_B CH4") || return_flag;

  return_flag = check_test(domega_A_N2_dCH4, dA_lib[0][1],"golden value domega_A N2 CH4")  || return_flag;
  return_flag = check_test(domega_B_N2_dCH4, dB_lib[0][1],"golden value domega_B N2 CH4")  || return_flag;
  return_flag = check_test(domega_A_CH4_dN2, dA_lib[1][0],"golden value domega_A CH4 N2")  || return_flag;
  return_flag = check_test(domega_B_CH4_dN2, dB_lib[1][0],"golden value domega_B CH4 N2")  || return_flag;
  return_flag = check_test(domega_A_N2_dN2,  dA_lib[0][0],"golden value domega_A N2 N2")   || return_flag;
  return_flag = check_test(domega_B_N2_dN2,  dB_lib[0][0],"golden value domega_B N2 N2")   || return_flag;
  return_flag = check_test(domega_A_CH4_dCH4,dA_lib[1][1],"golden value domega_A CH4 CH4") || return_flag;
  return_flag = check_test(domega_B_CH4_dCH4,dB_lib[1][1],"golden value domega_B CH4 CH4") || return_flag;

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
