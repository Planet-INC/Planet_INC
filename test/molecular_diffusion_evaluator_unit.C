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
#include "planet/molecular_diffusion_evaluator.h"
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
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 80.L;
  if(std::abs((theory-cal)/theory) < tol)return 0;
  std::cout << std::scientific << std::setprecision(20)
            << "\nfailed test: " << words << "\n"
            << "theory: " << theory
            << "\ncalculated: " << cal
            << "\ndifference: " << std::abs((theory-cal)/cal)
            << "\ntolerance: " << tol << std::endl;
  return 1;
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
     Scalar t,tz;
     temp >> t >> tz;
     if(!temp.good())break;
     T0.push_back(t);
     Tz.push_back(tz);
  }
  temp.close();
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

template<typename Scalar, typename VectorScalar>
void calculate_densities(VectorScalar &densities, const Scalar &tot_dens, const VectorScalar &molar_frac, 
                        const Scalar &zmin,const Scalar &z,
                        const Scalar &T, const VectorScalar &mm)
{
   Scalar Mm;
   Antioch::set_zero(Mm);
   for(unsigned int s = 0; s < molar_frac.size(); s++)
   {
      Mm += molar_frac[s] * mm[s];
   }

   densities.clear();
   densities.resize(molar_frac.size());
   for(unsigned int s = 0; s < molar_frac.size(); s++)
   {
     densities[s] = molar_frac[s] * barometry(zmin,z,T,Mm,tot_dens);
   }

   return;
}

template<typename Scalar>
Scalar binary_coefficient(const Scalar &T, const Scalar &P, const Scalar &D01, const Scalar &beta)
{
   return D01 * Planet::Constants::Convention::P_normal<Scalar>() / P * Antioch::ant_pow(T/Planet::Constants::Convention::T_standard<Scalar>(),beta);
}

template<typename Scalar>
Scalar dbinary_coefficient_dn(const Scalar &T, const Scalar &P, const Scalar &D01, const Scalar &beta)
{
   return - binary_coefficient(T,P,D01,beta) / P * Planet::Constants::Universal::kb<Scalar>() * T * 1e6L;
}

template<typename Scalar>
Scalar dbinary_coefficient_dT(const Scalar &T, const Scalar &P, const Scalar &D01, const Scalar &beta)
{
   return binary_coefficient(T,P,D01,beta) / T * (beta - Scalar(1.L));
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
   return n * 1e6L * Planet::Constants::Universal::kb<Scalar>() * T; //cm-3 -> m-3
}

template <typename Scalar>
int tester(const std::string &input_T)
{
//description
  std::vector<std::string> neutrals;
  std::vector<std::string> ions;
  neutrals.push_back("N2");
  neutrals.push_back("CH4");
  neutrals.push_back("C2H");
//ionic system contains neutral system
  ions = neutrals;
  Scalar MN(14.008e-3L), MC(12.011e-3L), MH(1.008e-3L);
  Scalar MN2 = 2.L*MN , MCH4 = MC + 4.L*MH, MC2H = 2.L * MC + MH;
  std::vector<Scalar> Mm;
  Mm.push_back(MN2);
  Mm.push_back(MCH4);
  Mm.push_back(MC2H);
  std::vector<std::string> medium;
  medium.push_back("N2");
  medium.push_back("CH4");

//densities
  std::vector<Scalar> molar_frac;
  molar_frac.push_back(0.50000L);
  molar_frac.push_back(0.04000L);
  molar_frac.push_back(0.46000L);
  Scalar dens_tot(1e12L); // cm^-3

//altitudes
  Scalar zmin(600.),zmax(1400.),zstep(10.);

//binary diffusion
  Scalar bCN1(1.04e-5L),bCN2(1.76L); //cm2
  Planet::DiffusionType CN_model(Planet::DiffusionType::Wakeham);
  Scalar bCC1(5.73e16L),bCC2(0.5L); //cm2
  Planet::DiffusionType CC_model(Planet::DiffusionType::Wilson);
  Scalar bNN1(0.1783L),bNN2(1.81L); //cm2
  Planet::DiffusionType NN_model(Planet::DiffusionType::Massman);

  std::vector<std::vector<Scalar> > Massman(3,std::vector<Scalar>(2,0.) );
// N2 - N2
  Massman[0][0] = bNN1; 
  Massman[0][1] = bNN2;
// N2 - CH4
  Massman[1][0] = bCN1 * Antioch::ant_pow(Planet::Constants::Convention::T_standard<Scalar>(),bCN2);
  Massman[1][1] = bCN2;
// CH4 - CH4
  Massman[2][0] = bCC1 * Antioch::ant_pow(Planet::Constants::Convention::T_standard<Scalar>(),bCC2 + Scalar(1.L)) * 
                              Planet::Constants::Universal::kb<Scalar>() / Planet::Constants::Convention::P_normal<Scalar>();
  Massman[2][1] = bCC2 + Scalar(1.L);

//neutrals
  Antioch::ChemicalMixture<Scalar> neutral_species(neutrals); 

//ions
  Antioch::ChemicalMixture<Scalar> ionic_species(ions); 

//binary diffusion
  Planet::BinaryDiffusion<Scalar> N2N2(   0, 0, bNN1, bNN2, NN_model);
  Planet::BinaryDiffusion<Scalar> N2CH4(  0, 1, bCN1, bCN2, CN_model);
  Planet::BinaryDiffusion<Scalar> CH4CH4( 1, 1, bCC1, bCC2, CC_model);
  Planet::BinaryDiffusion<Scalar> N2C2H(  0, 2);
  Planet::BinaryDiffusion<Scalar> CH4C2H( 1, 2);
  std::vector<std::vector<Planet::BinaryDiffusion<Scalar> > > bin_diff_coeff;
  bin_diff_coeff.resize(2);
  bin_diff_coeff[0].push_back(N2N2);
  bin_diff_coeff[0].push_back(N2CH4);
  bin_diff_coeff[0].push_back(N2C2H);
  bin_diff_coeff[1].push_back(N2CH4);
  bin_diff_coeff[1].push_back(CH4CH4);
  bin_diff_coeff[1].push_back(CH4C2H);

//temperature
  std::vector<Scalar> T0,Tz;
  read_temperature<Scalar>(T0,Tz,input_T);
  Planet::AtmosphericTemperature<Scalar, std::vector<Scalar> > temperature(Tz, T0);

//atmospheric mixture
  Planet::AtmosphericMixture<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > composition(neutral_species, ionic_species, temperature);
  composition.init_composition(molar_frac,dens_tot,zmin,zmax);

//molecular diffusion
  Planet::MolecularDiffusionEvaluator<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > molecular_diffusion(bin_diff_coeff,composition,temperature,medium);

/************************
 * checks
 ************************/

  Scalar Matm(0.L);
  for(unsigned int s = 0; s < molar_frac.size(); s++)
  {
     Matm += molar_frac[s] * composition.neutral_composition().M(s);
  }

//N2, CH4, C2H
  std::vector<std::vector<Scalar> > Dij;
  Dij.resize(2);
  Dij[0].resize(3,0.L);
  Dij[1].resize(3,0.L);

  std::vector<std::vector<std::vector<Scalar> > > dDij;
  dDij.resize(2);
  dDij[0].resize(3,std::vector<Scalar>(2,0.));
  dDij[1].resize(3,std::vector<Scalar>(2,0.));

  int return_flag(0);
  for(Scalar z = zmin; z <= zmax; z += zstep)
  {
      std::stringstream walt;
      walt << z;
      Scalar T    = temperature.neutral_temperature(z);
      Scalar nTot = barometry(zmin,z,T,Matm,dens_tot);
      Scalar P    = pressure(nTot,T);

      std::vector<Scalar> densities;
      calculate_densities(densities, dens_tot, molar_frac, zmin, z, T, Mm);

      std::vector<Scalar> molecular_diffusion_Dtilde(densities.size(),0.);
      std::vector<Scalar> molecular_diffusion_Dtilde_2(densities.size(),0.);
      std::vector<Scalar> molecular_diffusion_dDtilde_dT(densities.size(),0.);
      std::vector<std::vector<Scalar> > molecular_diffusion_dDtilde_dn(densities.size(), std::vector<Scalar>(densities.size(),0.));

      molecular_diffusion.Dtilde(densities,T,molecular_diffusion_Dtilde);
      molecular_diffusion.Dtilde_and_derivs_dn(densities,T,nTot,molecular_diffusion_Dtilde_2,molecular_diffusion_dDtilde_dn);
      molecular_diffusion.dDtilde_dT(densities, T, molecular_diffusion_dDtilde_dT);

      Dij[0][0] = binary_coefficient(T,P,Massman[0][0],Massman[0][1]); //N2 N2
      Dij[0][1] = binary_coefficient(T,P,Massman[1][0],Massman[1][1]); //N2 CH4
      Dij[0][2] = binary_coefficient(Dij[0][0],Mm[0],Mm[2]); //N2 C2H
      Dij[1][0] = Dij[0][1]; //CH4 N2
      Dij[1][1] = binary_coefficient(T,P,Massman[2][0],Massman[2][1]); //CH4 CH4
      Dij[1][2] = binary_coefficient(Dij[1][1],Mm[1],Mm[2]); //CH4 C2H

      dDij[0][0][0] = dbinary_coefficient_dn(T,P,Massman[0][0],Massman[0][1]); //N2 N2 dn
      dDij[0][0][1] = dbinary_coefficient_dT(T,P,Massman[0][0],Massman[0][1]); //N2 N2 dT
      dDij[0][1][0] = dbinary_coefficient_dn(T,P,Massman[1][0],Massman[1][1]); //N2 CH4 dn
      dDij[0][1][1] = dbinary_coefficient_dT(T,P,Massman[1][0],Massman[1][1]); //N2 CH4 dT
      dDij[0][2][0] = binary_coefficient(dbinary_coefficient_dn(T,P,Massman[0][0],Massman[0][1]),Mm[0],Mm[2]); //N2 C2H dn
      dDij[0][2][1] = binary_coefficient(dbinary_coefficient_dT(T,P,Massman[0][0],Massman[0][1]),Mm[0],Mm[2]); //N2 C2H dT
      dDij[1][0][0] = dDij[0][1][0];
      dDij[1][0][1] = dDij[0][1][1];
      dDij[1][1][0] = dbinary_coefficient_dn(T,P,Massman[2][0],Massman[2][1]); //CH4 CH4 dn
      dDij[1][1][1] = dbinary_coefficient_dT(T,P,Massman[2][0],Massman[2][1]); //CH4 CH4 dT
      dDij[1][2][0] = binary_coefficient(dbinary_coefficient_dn(T,P,Massman[2][0],Massman[2][1]),Mm[1],Mm[2]); //CH4 C2H dn
      dDij[1][2][1] = binary_coefficient(dbinary_coefficient_dT(T,P,Massman[2][0],Massman[2][1]),Mm[1],Mm[2]); //CH4 C2H dT

      return_flag = check_test(Dij[0][0],molecular_diffusion.binary_coefficient(0,0,T,P),"binary molecular coefficient N2 N2 at altitude " + walt.str())   || 
                    check_test(Dij[0][1],molecular_diffusion.binary_coefficient(0,1,T,P),"binary molecular coefficient N2 CH4 at altitude " + walt.str())  || 
                    check_test(Dij[0][2],molecular_diffusion.binary_coefficient(0,2,T,P),"binary molecular coefficient N2 C2H at altitude " + walt.str())  || 
                    check_test(Dij[1][0],molecular_diffusion.binary_coefficient(1,0,T,P),"binary molecular coefficient CH4 N2 at altitude " + walt.str())  || 
                    check_test(Dij[1][1],molecular_diffusion.binary_coefficient(1,1,T,P),"binary molecular coefficient CH4 CH4 at altitude " + walt.str()) || 
                    check_test(Dij[1][2],molecular_diffusion.binary_coefficient(1,2,T,P),"binary molecular coefficient CH4 C2H at altitude " + walt.str()) ||
                    check_test(dDij[0][0][0],molecular_diffusion.binary_coefficient_deriv_n(0,0,0,T,P,nTot),"binary molecular coefficient N2 N2 derivative with respect to n at altitude " + walt.str()) ||
                    check_test(dDij[0][0][1],molecular_diffusion.binary_coefficient_deriv_T(0,0,T,P),"binary molecular coefficient N2 N2 derivative with respect to T at altitude " + walt.str()) ||
                    check_test(dDij[0][1][0],molecular_diffusion.binary_coefficient_deriv_n(0,1,0,T,P,nTot),"binary molecular coefficient N2 CH4 derivative with respect to n at altitude " + walt.str()) ||
                    check_test(dDij[0][1][1],molecular_diffusion.binary_coefficient_deriv_T(0,1,T,P),"binary molecular coefficient N2 CH4 derivative with respect to T at altitude " + walt.str()) ||
                    check_test(dDij[0][2][0],molecular_diffusion.binary_coefficient_deriv_n(0,2,0,T,P,nTot),"binary molecular coefficient N2 C2H derivative with respect to n at altitude " + walt.str()) ||
                    check_test(dDij[0][2][1],molecular_diffusion.binary_coefficient_deriv_T(0,2,T,P),"binary molecular coefficient N2 C2H derivative with respect to T at altitude " + walt.str()) ||
                    check_test(dDij[1][0][0],molecular_diffusion.binary_coefficient_deriv_n(1,0,0,T,P,nTot),"binary molecular coefficient CH4 N2 derivative with respect to n at altitude " + walt.str()) ||
                    check_test(dDij[1][0][1],molecular_diffusion.binary_coefficient_deriv_T(1,0,T,P),"binary molecular coefficient CH4 N2 derivative with respect to T at altitude " + walt.str()) ||
                    check_test(dDij[1][1][0],molecular_diffusion.binary_coefficient_deriv_n(1,1,0,T,P,nTot),"binary molecular coefficient CH4 CH4 derivative with respect to n at altitude " + walt.str()) ||
                    check_test(dDij[1][1][1],molecular_diffusion.binary_coefficient_deriv_T(1,1,T,P),"binary molecular coefficient CH4 CH4 derivative with respect to T at altitude " + walt.str()) ||
                    check_test(dDij[1][2][0],molecular_diffusion.binary_coefficient_deriv_n(1,2,0,T,P,nTot),"binary molecular coefficient CH4 C2H derivative with respect to n at altitude " + walt.str()) ||
                    check_test(dDij[1][2][1],molecular_diffusion.binary_coefficient_deriv_T(1,2,T,P),"binary molecular coefficient CH4 C2H derivative with respect to T at altitude " + walt.str()) ||
                    return_flag;
// per species
      for(unsigned int s = 0; s < molar_frac.size(); s++)
      {

        Scalar tmp(0.L);
        Scalar dDs_dT(0.L);
        for(unsigned int imedium = 0; imedium < medium.size(); imedium++)
        {
           if(s == imedium)continue;
           tmp    += densities[imedium] / Dij[imedium][s];
           dDs_dT += densities[imedium] / (Dij[imedium][s] * Dij[imedium][s]) * dDij[imedium][s][1];
        }
        Scalar Ds = (nTot - densities[s]) / tmp;
        dDs_dT *= Ds / tmp;

        Scalar M_diff(0.L);
        Scalar totdens_diff = nTot - densities[s];
        for(unsigned int j = 0; j < molar_frac.size(); j++)
        {
           if(s == j)continue;
           M_diff += densities[j] * Mm[j];
        }
        M_diff /= totdens_diff;

        Scalar Dtilde = Ds / (Scalar(1.L) - molar_frac[s] * (Scalar(1.L) - Mm[s]/M_diff));
        dDs_dT /= (Scalar(1.L) - molar_frac[s] * (Scalar(1.L) - Mm[s]/M_diff));

        Scalar Dtilde_3 = molecular_diffusion.Dtilde(s, nTot, T, P, densities);

        return_flag = check_test(Dtilde,molecular_diffusion_Dtilde[s],     "Dtilde of species "   + neutrals[s] + " at altitude " + walt.str()) || 
                      check_test(Dtilde,molecular_diffusion_Dtilde_2[s],   "Dtilde 2 of species " + neutrals[s] + " at altitude " + walt.str()) || 
                      check_test(Dtilde,Dtilde_3,                          "Dtilde 3 of species " + neutrals[s] + " at altitude " + walt.str()) || 
                      check_test(dDs_dT,molecular_diffusion_dDtilde_dT[s], "Dtilde derived with respect to T of species " + neutrals[s] + " at altitude " + walt.str()) || 
                      return_flag;

// now the derivatives with n
        for(unsigned int k = 0; k < molar_frac.size(); k++)
        {
           Scalar dDs_dn;
           Antioch::set_zero(dDs_dn);
           Scalar dDs_dn_common;
           Antioch::set_zero(dDs_dn_common);
           for(unsigned int m = 0; m < medium.size(); m++)
           {
              if(m == s)continue;
              dDs_dn += densities[m] / (Dij[m][s] * Dij[m][s]) * dDij[m][s][0];
              dDs_dn_common += densities[m] / (Dij[m][s] * Dij[m][s]) * dDij[m][s][0];
              if(m == k)dDs_dn -= Antioch::constant_clone(T,1) / Dij[m][s];
           }
           dDs_dn *= Ds / tmp;
           if(s != k) dDs_dn += Ds / (nTot - densities[s]);


           Scalar dDtilde_dn;
           Antioch::set_zero(dDtilde_dn);
           if(k == s)dDtilde_dn = -Antioch::constant_clone(T,1) / nTot;
           dDtilde_dn += densities[s] / (nTot * nTot);
           dDtilde_dn *= Antioch::constant_clone(T,1) - Mm[s] / M_diff;
           if(k != s)dDtilde_dn += (Mm[k] - M_diff) / (nTot - densities[s]) * Mm[s] / (M_diff * M_diff) * densities[s] / nTot;
           dDtilde_dn *= (- Antioch::constant_clone(T,1) / (Antioch::constant_clone(T,1) - molar_frac[s] * (Antioch::constant_clone(T,1) - Mm[s]/M_diff) ));
           dDtilde_dn += dDs_dn / Ds;
           dDtilde_dn *= Dtilde;

           return_flag = check_test(dDtilde_dn,molecular_diffusion_dDtilde_dn[s][k],"Dtilde derivative with respect to " + neutrals[k] + " of species " + neutrals[s] + " at altitude " + walt.str()) || 
                         return_flag;
        }

      }
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
