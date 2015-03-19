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
int check_test(Scalar theory, Scalar cal, const std::string &words, const Scalar & tol, Scalar & max_diff)
{
  const Scalar diff = std::abs((theory-cal)/theory);
  if(diff > max_diff)max_diff = diff;

  if(diff < tol)return 0;
  std::cout << std::scientific << std::setprecision(20)
            << "\nfailed test: " << words << "\n"
            << "theory: " << theory
            << "\ncalculated: " << cal
            << "\ndifference: " << std::abs((theory-cal)/theory)
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
     Scalar t(0),tz(0);
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
   return botdens * Antioch::ant_exp(-(z - zmin)/((Planet::Constants::Titan::radius<Scalar>() + z) * (Planet::Constants::Titan::radius<Scalar>() + zmin) * Scalar(1e3) *
                                             Antioch::Constants::Avogadro<Scalar>() * Planet::Constants::Universal::kb<Scalar>() * T / 
                                                        (Planet::Constants::Universal::G<Scalar>() * Planet::Constants::Titan::mass<Scalar>() * Mm * Scalar(1e-3)))
                              );
}

template<typename Scalar, typename VectorScalar>
void calculate_densities(VectorScalar &densities, const VectorScalar &molar_frac, const Scalar &nTot )
{
   densities.resize(molar_frac.size(),0);
   for(unsigned int s = 0; s < molar_frac.size(); s++)
   {
     densities[s] = molar_frac[s] * nTot;
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
   return (Mj < Mi)?Dii * Antioch::ant_sqrt((Mj/Mi + Scalar(1))/Scalar(2))
                   :
                    Dii * Antioch::ant_sqrt((Mj/Mi));
}

template<typename Scalar>
Scalar pressure(const Scalar &n, const Scalar &T)
{
   return n * Scalar(1e6) * Planet::Constants::Universal::kb<Scalar>() * T;
}

template<typename Scalar>
Scalar scale_height(const Scalar &T, const Scalar &z, const Scalar &Mm)
{
  return Scalar(1e-3) * Antioch::Constants::R_universal<Scalar>() * T / 
         (Planet::Constants::g<Scalar>(Planet::Constants::Titan::radius<Scalar>(),z,Planet::Constants::Titan::mass<Scalar>()) * Mm);
}


template <typename Scalar>
int tester(const std::string &input_T, const std::string & type)
{
//description
  std::vector<std::string> neutrals;
  std::vector<std::string> ions;
  neutrals.push_back("N2");
  neutrals.push_back("CH4");
  neutrals.push_back("H2");
//ionic system contains neutral system
  ions = neutrals;
  Scalar MN(14.008e-3L), MC(12.011e-3L), MH(1.008e-3L);
  Scalar MN2 = 2 * MN , MCH4 = MC + 4 * MH, MH2 = 2 * MH;
  std::vector<Scalar> Mm;
  Mm.push_back(MN2);
  Mm.push_back(MCH4);
  Mm.push_back(MH2);
  std::vector<std::string> medium;
  medium.push_back("N2");
  medium.push_back("CH4");

//densities
  std::vector<Scalar> molar_frac;
  molar_frac.push_back(0.98525004881495660129211947533421479L);
  molar_frac.push_back(0.01414039822253783945793469946368639L);
  molar_frac.push_back(0.00060955296250555924994582520209881L);
  Scalar dens_tot(1.7e13L);

//zenith angle
//not necessary

//photon flux
//not necessary

////cross-section
//not necessary

//altitudes
  Scalar zmin(600),zmax(1400),zstep(50);

//binary diffusion
//N2 - N2
  Scalar bNN1(5.09e16),bNN2(0.81L); //cm2
  Planet::DiffusionType NN_model(Planet::DiffusionType::Wilson);

//N2 - CH4
  Scalar bCN1(7.34e16L),bCN2(0.75L); //cm2
  Planet::DiffusionType CN_model(Planet::DiffusionType::Wilson);

//CH4 - CH4
  Scalar bCC1(5.73e16L),bCC2(0.5); //cm2
  Planet::DiffusionType CC_model(Planet::DiffusionType::Wilson);

//N2 - H2
  Scalar bNH1(1.88e17L),bNH2(0.82L); //cm2
  Planet::DiffusionType NH_model(Planet::DiffusionType::Wilson);

//CH4 - H2
  Scalar bCH1(2.3e17L),bCH2(0.765L); //cm2
  Planet::DiffusionType CH_model(Planet::DiffusionType::Wilson);

  std::vector<std::vector<Scalar> > Massman(5,std::vector<Scalar>(2,0) );
// N2 - N2
  Massman[0][0] = bNN1 * Antioch::ant_pow(Planet::Constants::Convention::T_standard<Scalar>(),bNN2 + 1) *  
                              Planet::Constants::Universal::kb<Scalar>() / Planet::Constants::Convention::P_normal<Scalar>();
  Massman[0][1] = bNN2 + 1;
// N2 - CH4
  Massman[1][0] = bCN1 * Antioch::ant_pow(Planet::Constants::Convention::T_standard<Scalar>(),bCN2 + 1) * 
                              Planet::Constants::Universal::kb<Scalar>() / Planet::Constants::Convention::P_normal<Scalar>();
  Massman[1][1] = bCN2 + 1;
// CH4 - CH4
  Massman[2][0] = bCC1 * Antioch::ant_pow(Planet::Constants::Convention::T_standard<Scalar>(),bCC2 + 1) * 
                              Planet::Constants::Universal::kb<Scalar>() / Planet::Constants::Convention::P_normal<Scalar>();
  Massman[2][1] = bCC2 + 1;
// N2 - H2
  Massman[3][0] = bNH1 * Antioch::ant_pow(Planet::Constants::Convention::T_standard<Scalar>(),bNH2 + 1) * 
                              Planet::Constants::Universal::kb<Scalar>() / Planet::Constants::Convention::P_normal<Scalar>();
  Massman[3][1] = bNH2 + 1;
// CH4 - H2
  Massman[4][0] = bCH1 * Antioch::ant_pow(Planet::Constants::Convention::T_standard<Scalar>(),bCH2 + 1) * 
                              Planet::Constants::Universal::kb<Scalar>() / Planet::Constants::Convention::P_normal<Scalar>();
  Massman[4][1] = bCH2 + 1;


//thermal coefficient
  std::vector<Scalar> tc;
  tc.push_back(0);   // N2
  tc.push_back(0);   // CH4
  tc.push_back(-0.38L); // H2

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
  Planet::BinaryDiffusion<Scalar> N2H2(   0, 2, bNH1, bNH2, NH_model);
  Planet::BinaryDiffusion<Scalar> CH4H2(  1, 2, bCH1, bCH2, CH_model);
  std::vector<std::vector<Planet::BinaryDiffusion<Scalar> > > bin_diff_coeff;
  bin_diff_coeff.resize(2);
  bin_diff_coeff[0].push_back(N2N2);
  bin_diff_coeff[0].push_back(N2CH4);
  bin_diff_coeff[0].push_back(N2H2);
  bin_diff_coeff[1].push_back(N2CH4);
  bin_diff_coeff[1].push_back(CH4CH4);
  bin_diff_coeff[1].push_back(CH4H2);


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
  for(unsigned int s = 0; s < Dij.size(); s++)
  {
    Dij[s].resize(3,0);
  }

  std::vector<Scalar> Dtilde;
  Dtilde.resize(molar_frac.size(),0);

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 6000;
  Scalar max_diff(-1);
  std::cout << "Type: " << type << ", tolerance = " << tol << " ... ";

  int return_flag(0);
  for(Scalar z = zmin; z <= zmax; z += zstep)
  {
     std::stringstream walt;
     walt << z;

     Scalar T        = temperature.neutral_temperature(z);
     Scalar dT_dz    = temperature.dneutral_temperature_dz(z);
     Scalar nTot     = barometry(zmin,z,T,mean_M,dens_tot);
     Scalar P        = pressure(nTot,T);
     Scalar Ha       = scale_height(T,z,mean_M);

     std::vector<Scalar> densities;
     calculate_densities(densities,molar_frac,nTot);

//eddy
     Scalar K = K0 * Antioch::ant_sqrt(dens_tot/nTot);

     Dij[0][0] = binary_coefficient(T,P,Massman[0][0],Massman[0][1]); //N2 N2
     Dij[0][1] = binary_coefficient(T,P,Massman[1][0],Massman[1][1]); //N2 CH4
     Dij[0][2] = binary_coefficient(T,P,Massman[3][0],Massman[3][1]); //N2 H2
     Dij[1][0] = Dij[0][1]; //CH4 N2
     Dij[1][1] = binary_coefficient(T,P,Massman[2][0],Massman[2][1]); //CH4 CH4
     Dij[1][2] = binary_coefficient(T,P,Massman[4][0],Massman[4][1]); //CH4 H2

     std::vector<Scalar> omega_a(densities.size()),omega_b(densities.size());
     diffusion.diffusion(densities,z,omega_a,omega_b);

     for(unsigned int s = 0; s < molar_frac.size(); s++)
     {
       Scalar tmp(0);
       for(unsigned int medium = 0; medium < 2; medium++)
       {
          if(s == medium)continue;
          tmp += densities[medium]/Dij[medium][s];
       }
       Scalar Ds = (nTot - densities[s]) / tmp;

       Scalar M_diff(0);
       Scalar totdens_diff(0);
       for(unsigned int j = 0; j < molar_frac.size(); j++)
       {
          if(s == j)continue;
          M_diff += densities[j] * Mm[j];
          totdens_diff += densities[j];
       }
       M_diff /= totdens_diff;

       Scalar Dtilde = Ds / (1 - molar_frac[s] * 
                              (1 - Mm[s] / M_diff)
                            );
//
       Scalar Hs = scale_height(T,z,Mm[s]); //km
       Scalar omega_a_theo = - (Dtilde + K) * Scalar(1e-10);
       Scalar omega_b_theo = - Scalar(1e-10) *
                               (   Dtilde / Hs 
                                 + Dtilde * dT_dz / T * (1 + (1 - molar_frac[s]) * tc[s])
                                 + K / Ha
                                 + K * dT_dz / T
                               );

       return_flag = check_test(omega_a_theo,omega_a[s],"omega A of species " + neutrals[s] + " at altitude " + walt.str(), tol, max_diff) || return_flag;
       return_flag = check_test(omega_b_theo,omega_b[s],"omega B of species " + neutrals[s] + " at altitude " + walt.str(), tol, max_diff) || return_flag;
     }
  }

  const Scalar z(1.00126188021535153894e3L);
  const Scalar nN2(1.67492504780800000000e13L);
  const Scalar nCH4(2.40386768896000000000e11L);
  const Scalar nH2(1.03623997440000000000e10L);

  std::vector<Scalar> densities(3,0.);
  densities[0] = nN2;
  densities[1] = nCH4;
  densities[2] = nH2;

  std::vector<Scalar> omega_A(3,0.);
  std::vector<Scalar> omega_B(3,0.);
  std::vector<std::vector<Scalar> > domega_A_dn(3,std::vector<Scalar>(3,0.));
  std::vector<std::vector<Scalar> > domega_B_dn(3,std::vector<Scalar>(3,0.));

  diffusion.diffusion_and_derivs(densities,z,omega_A,omega_B,domega_A_dn,domega_B_dn);

  std::vector<Scalar> A_theory(3,0.),B_theory(3,0.);
  std::vector<std::vector<Scalar> > dA_theory(3,std::vector<Scalar>(3,0.)),
                                    dB_theory(3,std::vector<Scalar>(3,0.));
  A_theory[0] = -0.0004300000176435734125008421080247927451079;
  A_theory[1] = -0.0004300000273536820303058928815120652201946;
  A_theory[2] = -0.0004300000886980048406228781609380718626931;

  B_theory[0] = -4.27421456269783728448831552830770374850214667e-6;
  B_theory[1] = -4.27421454279530772749581406930371947621281983e-6;
  B_theory[2] = -4.27421440077198098665129413126672792396876733e-6;

  dA_theory[0][0] = 1.26470599979531896797220504230e-17;
  dA_theory[0][1] = 1.264705995983158054989764521749298e-17;
  dA_theory[0][2] = 1.264705242344990091473083187836913e-17;
  dA_theory[1][0] = 1.2647060572996819308235210972765e-17;
  dA_theory[1][1] = 1.26470599937842429008274340768931e-17;
  dA_theory[1][2] = 1.26470592169591530184670015105645e-17;
  dA_theory[2][0] = 1.26470641687479427915876977803256e-17;
  dA_theory[2][1] = 1.2647064629892737622485275516851780e-17;
  dA_theory[2][2] = 1.264705957666781620405559513885556e-17;

  dB_theory[0][0] = 1.2370924919932135783887209825042e-19;
  dB_theory[0][1] = 2.5323205556404673293013943589522e-19;
  dB_theory[0][2] = 4.0497476884979917624445017703953e-19;
  dB_theory[1][0] = 1.2370924802109365910395142572773e-19;
  dB_theory[1][1] = 2.5323205193046541487278097547526e-19;
  dB_theory[1][2] = 4.0497483691956866256085514421392e-19;
  dB_theory[2][0] = 1.2370923962618499256633961848765e-19;
  dB_theory[2][1] = 2.5323204620907203470505917733277e-19;
  dB_theory[2][2] = 4.0497483310274272097250604001168e-19;

  for(unsigned int s = 0; s < densities.size(); s++)
  {
    return_flag = check_test(A_theory[s],omega_A[s], "golden value omega_A " + neutrals[s], tol, max_diff)  || return_flag;
    return_flag = check_test(B_theory[s],omega_B[s], "golden value omega_B " + neutrals[s], tol, max_diff)  || return_flag;
    for(unsigned int k = 0; k < densities.size(); k++)
    {
       return_flag = check_test(dA_theory[s][k], domega_A_dn[s][k], "golden value domega_A " + neutrals[s] + " " + neutrals[k], tol, max_diff)  || return_flag;
       return_flag = check_test(dB_theory[s][k], domega_B_dn[s][k], "golden value domega_B " + neutrals[s] + " " + neutrals[k], tol, max_diff)  || return_flag;
    }
  }

  std::cout << "max diff = " << max_diff << std::endl;

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

  return (tester<float>(std::string(argv[1]), "float") ||
          tester<double>(std::string(argv[1]), "double"));// ||
//          tester<long double>(std::string(argv[1]), "long double"));
}
