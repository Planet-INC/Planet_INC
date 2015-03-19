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
int check_test(Scalar theory, Scalar cal, const std::string &words, const Scalar & tol, Scalar & max_diff)
{
  

  Scalar diff = std::abs((theory-cal)/theory);
  if(max_diff < diff)max_diff = diff;

  if(diff < tol)return 0;
  std::cout << std::scientific << std::setprecision(20)
            << "\nfailed test: " << words << "\n"
            << "theory: " << theory
            << "\ncalculated: " << cal
            << "\ndifference: " << diff
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
   return - binary_coefficient(T,P,D01,beta) / P * Planet::Constants::Universal::kb<Scalar>() * T * Scalar(1e6);
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
   return n * Scalar(1e6) * Planet::Constants::Universal::kb<Scalar>() * T; //cm-3 -> m-3
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
  Scalar dens_tot(1.7e13L); // cm^-3

//altitudes
  Scalar zmin(600.),zmax(1400.),zstep(50.);

//binary diffusion
//N2 - N2
  Scalar bNN1(5.09e16),bNN2(0.81L); //cm2
  Planet::DiffusionType NN_model(Planet::DiffusionType::Wilson);

//N2 - CH4
  Scalar bCN1(7.34e16L),bCN2(0.75L); //cm2
  Planet::DiffusionType CN_model(Planet::DiffusionType::Wilson);

//CH4 - CH4
  Scalar bCC1(5.73e16L),bCC2(0.5L); //cm2
  Planet::DiffusionType CC_model(Planet::DiffusionType::Wilson);

//N2 - H2
  Scalar bNH1(1.88e17L),bNH2(0.82L); //cm2
  Planet::DiffusionType NH_model(Planet::DiffusionType::Wilson);

//CH4 - H2
  Scalar bCH1(2.3e17L),bCH2(0.765L); //cm2
  Planet::DiffusionType CH_model(Planet::DiffusionType::Wilson);

  std::vector<std::vector<Scalar> > Massman(5,std::vector<Scalar>(2,0.) );
// N2 - N2
  Massman[0][0] = bNN1 * Antioch::ant_pow(Planet::Constants::Convention::T_standard<Scalar>(),bNN2 + Scalar(1.L)) *  
                              Planet::Constants::Universal::kb<Scalar>() / Planet::Constants::Convention::P_normal<Scalar>();
  Massman[0][1] = bNN2 + Scalar(1.L);
// N2 - CH4
  Massman[1][0] = bCN1 * Antioch::ant_pow(Planet::Constants::Convention::T_standard<Scalar>(),bCN2 + Scalar(1.L)) * 
                              Planet::Constants::Universal::kb<Scalar>() / Planet::Constants::Convention::P_normal<Scalar>();
  Massman[1][1] = bCN2 + Scalar(1.L);
// CH4 - CH4
  Massman[2][0] = bCC1 * Antioch::ant_pow(Planet::Constants::Convention::T_standard<Scalar>(),bCC2 + Scalar(1.L)) * 
                              Planet::Constants::Universal::kb<Scalar>() / Planet::Constants::Convention::P_normal<Scalar>();
  Massman[2][1] = bCC2 + Scalar(1.L);
// N2 - H2
  Massman[3][0] = bNH1 * Antioch::ant_pow(Planet::Constants::Convention::T_standard<Scalar>(),bNH2 + Scalar(1.L)) * 
                              Planet::Constants::Universal::kb<Scalar>() / Planet::Constants::Convention::P_normal<Scalar>();
  Massman[3][1] = bNH2 + Scalar(1.L);
// CH4 - H2
  Massman[4][0] = bCH1 * Antioch::ant_pow(Planet::Constants::Convention::T_standard<Scalar>(),bCH2 + Scalar(1.L)) * 
                              Planet::Constants::Universal::kb<Scalar>() / Planet::Constants::Convention::P_normal<Scalar>();
  Massman[4][1] = bCH2 + Scalar(1.L);

//neutrals
  Antioch::ChemicalMixture<Scalar> neutral_species(neutrals); 

//ions
  Antioch::ChemicalMixture<Scalar> ionic_species(ions); 

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

//N2, CH4, H2
  std::vector<std::vector<Scalar> > Dij;
  Dij.resize(2);
  Dij[0].resize(3,0);
  Dij[1].resize(3,0);

  std::vector<std::vector<std::vector<Scalar> > > dDij;
  dDij.resize(2);
  dDij[0].resize(3,std::vector<Scalar>(2,0));
  dDij[1].resize(3,std::vector<Scalar>(2,0));

  int return_flag(0);

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 4000;
  Scalar max_diff(-1);

  std::cout << type << " tolerance: " << tol << " and max diff ... ";

  for(Scalar z = zmin; z <= zmax; z += zstep)
  {
      std::stringstream walt;
      walt << z;
      Scalar T    = temperature.neutral_temperature(z);
      Scalar nTot = barometry(zmin,z,T,Matm,dens_tot);
      Scalar P    = pressure(nTot,T);

      std::vector<Scalar> densities;
      calculate_densities(densities, dens_tot, molar_frac, zmin, z, T, Mm);

      std::vector<Scalar> molecular_diffusion_Dtilde(densities.size(),0);
      std::vector<Scalar> molecular_diffusion_Dtilde_2(densities.size(),0);
      std::vector<Scalar> molecular_diffusion_dDtilde_dT(densities.size(),0);
      std::vector<std::vector<Scalar> > molecular_diffusion_dDtilde_dn(densities.size(), std::vector<Scalar>(densities.size(),0.));

      molecular_diffusion.Dtilde(densities,T,molecular_diffusion_Dtilde);
      molecular_diffusion.Dtilde_and_derivs_dn(densities,T,nTot,molecular_diffusion_Dtilde_2,molecular_diffusion_dDtilde_dn);
      molecular_diffusion.dDtilde_dT(densities, T, molecular_diffusion_dDtilde_dT);

      Dij[0][0] = binary_coefficient(T,P,Massman[0][0],Massman[0][1]); //N2 N2
      Dij[0][1] = binary_coefficient(T,P,Massman[1][0],Massman[1][1]); //N2 CH4
      Dij[0][2] = binary_coefficient(T,P,Massman[3][0],Massman[3][1]); //N2 H2
      Dij[1][0] = Dij[0][1]; //CH4 N2
      Dij[1][1] = binary_coefficient(T,P,Massman[2][0],Massman[2][1]); //CH4 CH4
      Dij[1][2] = binary_coefficient(T,P,Massman[4][0],Massman[4][1]); //CH4 H2

      dDij[0][0][0] = dbinary_coefficient_dn(T,P,Massman[0][0],Massman[0][1]); //N2 N2 dn
      dDij[0][0][1] = dbinary_coefficient_dT(T,P,Massman[0][0],Massman[0][1]); //N2 N2 dT
      dDij[0][1][0] = dbinary_coefficient_dn(T,P,Massman[1][0],Massman[1][1]); //N2 CH4 dn
      dDij[0][1][1] = dbinary_coefficient_dT(T,P,Massman[1][0],Massman[1][1]); //N2 CH4 dT
      dDij[0][2][0] = dbinary_coefficient_dn(T,P,Massman[3][0],Massman[3][1]); //N2 H2 dn
      dDij[0][2][1] = dbinary_coefficient_dT(T,P,Massman[3][0],Massman[3][1]); //N2 H2 dT
      dDij[1][0][0] = dDij[0][1][0];
      dDij[1][0][1] = dDij[0][1][1];
      dDij[1][1][0] = dbinary_coefficient_dn(T,P,Massman[2][0],Massman[2][1]); //CH4 CH4 dn
      dDij[1][1][1] = dbinary_coefficient_dT(T,P,Massman[2][0],Massman[2][1]); //CH4 CH4 dT
      dDij[1][2][0] = dbinary_coefficient_dn(T,P,Massman[4][0],Massman[4][1]); //CH4 H2 dn
      dDij[1][2][1] = dbinary_coefficient_dT(T,P,Massman[4][0],Massman[4][1]); //CH4 H2 dT

      return_flag = check_test(Dij[0][0],molecular_diffusion.binary_coefficient(0,0,T,P),"binary molecular coefficient N2 N2 at altitude " + walt.str(),tol,max_diff)   || 
                    check_test(Dij[0][1],molecular_diffusion.binary_coefficient(0,1,T,P),"binary molecular coefficient N2 CH4 at altitude " + walt.str(),tol,max_diff)  || 
                    check_test(Dij[0][2],molecular_diffusion.binary_coefficient(0,2,T,P),"binary molecular coefficient N2 H2 at altitude " + walt.str(),tol,max_diff)  || 
                    check_test(Dij[1][0],molecular_diffusion.binary_coefficient(1,0,T,P),"binary molecular coefficient CH4 N2 at altitude " + walt.str(),tol,max_diff)  || 
                    check_test(Dij[1][1],molecular_diffusion.binary_coefficient(1,1,T,P),"binary molecular coefficient CH4 CH4 at altitude " + walt.str(),tol,max_diff) || 
                    check_test(Dij[1][2],molecular_diffusion.binary_coefficient(1,2,T,P),"binary molecular coefficient CH4 H2 at altitude " + walt.str(),tol,max_diff) ||
                    check_test(dDij[0][0][0],molecular_diffusion.binary_coefficient_deriv_n(0,0,0,T,P,nTot),"binary molecular coefficient N2 N2 derivative with respect to n at altitude " + walt.str(),tol,max_diff) ||
                    check_test(dDij[0][0][1],molecular_diffusion.binary_coefficient_deriv_T(0,0,T,P),"binary molecular coefficient N2 N2 derivative with respect to T at altitude " + walt.str(),tol,max_diff) ||
                    check_test(dDij[0][1][0],molecular_diffusion.binary_coefficient_deriv_n(0,1,0,T,P,nTot),"binary molecular coefficient N2 CH4 derivative with respect to n at altitude " + walt.str(),tol,max_diff) ||
                    check_test(dDij[0][1][1],molecular_diffusion.binary_coefficient_deriv_T(0,1,T,P),"binary molecular coefficient N2 CH4 derivative with respect to T at altitude " + walt.str(),tol,max_diff) ||
                    check_test(dDij[0][2][0],molecular_diffusion.binary_coefficient_deriv_n(0,2,0,T,P,nTot),"binary molecular coefficient N2 H2 derivative with respect to n at altitude " + walt.str(),tol,max_diff) ||
                    check_test(dDij[0][2][1],molecular_diffusion.binary_coefficient_deriv_T(0,2,T,P),"binary molecular coefficient N2 H2 derivative with respect to T at altitude " + walt.str(),tol,max_diff) ||
                    check_test(dDij[1][0][0],molecular_diffusion.binary_coefficient_deriv_n(1,0,0,T,P,nTot),"binary molecular coefficient CH4 N2 derivative with respect to n at altitude " + walt.str(),tol,max_diff) ||
                    check_test(dDij[1][0][1],molecular_diffusion.binary_coefficient_deriv_T(1,0,T,P),"binary molecular coefficient CH4 N2 derivative with respect to T at altitude " + walt.str(),tol,max_diff) ||
                    check_test(dDij[1][1][0],molecular_diffusion.binary_coefficient_deriv_n(1,1,0,T,P,nTot),"binary molecular coefficient CH4 CH4 derivative with respect to n at altitude " + walt.str(),tol,max_diff) ||
                    check_test(dDij[1][1][1],molecular_diffusion.binary_coefficient_deriv_T(1,1,T,P),"binary molecular coefficient CH4 CH4 derivative with respect to T at altitude " + walt.str(),tol,max_diff) ||
                    check_test(dDij[1][2][0],molecular_diffusion.binary_coefficient_deriv_n(1,2,0,T,P,nTot),"binary molecular coefficient CH4 H2 derivative with respect to n at altitude " + walt.str(),tol,max_diff) ||
                    check_test(dDij[1][2][1],molecular_diffusion.binary_coefficient_deriv_T(1,2,T,P),"binary molecular coefficient CH4 H2 derivative with respect to T at altitude " + walt.str(),tol,max_diff) ||
                    return_flag;
// per species
      for(unsigned int s = 0; s < molar_frac.size(); s++)
      {

        Scalar tmp(0);
        Scalar dDs_dT(0);
        for(unsigned int imedium = 0; imedium < medium.size(); imedium++)
        {
           if(s == imedium)continue;
           tmp    += densities[imedium] / Dij[imedium][s];
           dDs_dT += densities[imedium] / (Dij[imedium][s] * Dij[imedium][s]) * dDij[imedium][s][1];
        }
        Scalar Ds = (nTot - densities[s]) / tmp;
        dDs_dT *= Ds / tmp;

        Scalar M_diff(0);
        Scalar totdens_diff = nTot - densities[s];
        for(unsigned int j = 0; j < molar_frac.size(); j++)
        {
           if(s == j)continue;
           M_diff += densities[j] * Mm[j];
        }
        M_diff /= totdens_diff;

        Scalar Dtilde = Ds / (1 - molar_frac[s] * (1 - Mm[s]/M_diff));
        dDs_dT /= (1 - molar_frac[s] * (1 - Mm[s]/M_diff));

        Scalar Dtilde_3 = molecular_diffusion.Dtilde(s, nTot, T, P, densities);

        return_flag = check_test(Dtilde,molecular_diffusion_Dtilde[s],     "Dtilde of species "   + neutrals[s] + " at altitude " + walt.str(),tol,max_diff) || 
                      check_test(Dtilde,molecular_diffusion_Dtilde_2[s],   "Dtilde 2 of species " + neutrals[s] + " at altitude " + walt.str(),tol,max_diff) || 
                      check_test(Dtilde,Dtilde_3,                          "Dtilde 3 of species " + neutrals[s] + " at altitude " + walt.str(),tol,max_diff) || 
                      check_test(dDs_dT,molecular_diffusion_dDtilde_dT[s], "Dtilde derived with respect to T of species " + neutrals[s] + " at altitude " + walt.str(),tol,max_diff) || 
                      return_flag;

// now the derivatives with n
        for(unsigned int k = 0; k < molar_frac.size(); k++)
        {
           Scalar term_bimol;
           Antioch::set_zero(term_bimol);
           if(s != k && k < medium.size())
           {
             Scalar sum(0);
             for(unsigned int m = 0; m < medium.size(); m++)
             {
                if(m == s)continue;
                sum += densities[m] / Dij[m][s];
             }
              term_bimol = Antioch::constant_clone(T,1) / (sum * Dij[k][s]);
           }

           Scalar dDtilde_dn(-densities[s]/(nTot*nTot));
           if(k == s)dDtilde_dn += Antioch::constant_clone(T,1) / nTot;
           dDtilde_dn *= (Antioch::constant_clone(T,1) - Mm[s] / M_diff);
           if(k != s)dDtilde_dn += (Mm[k] - M_diff) / (nTot - densities[s]) * Mm[s] / (M_diff * M_diff) * densities[s] / nTot;
           dDtilde_dn *= Dtilde / Ds;
           dDtilde_dn -= Antioch::constant_clone(T,1) / nTot;
           if(s != k && k < medium.size())dDtilde_dn -= term_bimol;
           if(s != k)dDtilde_dn += Antioch::constant_clone(T,1) / (nTot - densities[s]);
           dDtilde_dn *= Dtilde;

           return_flag = check_test(dDtilde_dn,molecular_diffusion_dDtilde_dn[s][k],"Dtilde derivative with respect to " + neutrals[k] + " of species " + neutrals[s] + " at altitude " + walt.str(),tol,max_diff) || 
                         return_flag;
        }

      }
  }
  std::vector<Scalar> densities(3,0.);
  densities[0] = 1.67492504780800000000e13L;
  densities[1] = 2.40386768896000000000e11L;
  densities[2] = 1.03623997440000000000e10L;
  const Scalar nTot = densities[0] + densities[1] + densities[2];
  const Scalar T = 197.3444020114407067012507468461990356;
  std::vector<Scalar> test(3,0.);
  std::vector<std::vector<Scalar> > dtest(3,std::vector<Scalar>(3,0.));

  molecular_diffusion.Dtilde_and_derivs_dn(densities,T,nTot,test,dtest);

  const Scalar exact_N2(0.13175620959620829365933043902850430735644492531321);
  const Scalar dexact_N2_N2(-7.80192658843725743470610638359192186e-15);
  const Scalar dexact_N2_CH4(-7.42071048099607652288172861518391368e-15);
  const Scalar dexact_N2_H2(6.794310950670392087834409003544762293e-14);

  const Scalar exact_CH4(0.22885729988609127375665445647101399645508956757757);
  const Scalar dexact_CH4_N2(-1.355236312822976495111972234499004441e-14);
  const Scalar dexact_CH4_CH4(-7.76023711888293341444280146924910631e-15);
  const Scalar dexact_CH4_H2(8.01410889421362104507059311510769e-18);

  const Scalar exact_H2(0.84230055725467281158924055356031866834452738609897);
  const Scalar dexact_H2_N2(-4.950987607664367816185865630216608636e-14);
  const Scalar dexact_H2_CH4(-5.41213244060823733480112116153913187e-14);
  const Scalar dexact_H2_H2(-3.58907268929640668494924664343440539e-15);

  return_flag = check_test(exact_N2,test[0],"Dtilde N2 compared to bc",tol,max_diff) || return_flag;
  return_flag = check_test(dexact_N2_N2, dtest[0][0],"Dtilde N2 derived with respect to N2 compared to bc",tol,max_diff)  || return_flag;
  return_flag = check_test(dexact_N2_CH4,dtest[0][1],"Dtilde N2 derived with respect to CH4 compared to bc",tol,max_diff) || return_flag;
  return_flag = check_test(dexact_N2_H2, dtest[0][2],"Dtilde N2 derived with respect to H2 compared to bc",tol,max_diff)  || return_flag;

  return_flag = check_test(exact_CH4,test[1],"Dtilde CH4 compared to bc",tol,max_diff) || return_flag;
  return_flag = check_test(dexact_CH4_N2, dtest[1][0],"Dtilde CH4 derived with respect to N2 compared to bc",tol,max_diff)  || return_flag;
  return_flag = check_test(dexact_CH4_CH4,dtest[1][1],"Dtilde CH4 derived with respect to CH4 compared to bc",tol,max_diff) || return_flag;
  return_flag = check_test(dexact_CH4_H2, dtest[1][2],"Dtilde CH4 derived with respect to H2 compared to bc",tol,max_diff)  || return_flag;

  return_flag = check_test(exact_H2,test[2],"Dtilde H2 compared to bc",tol,max_diff) || return_flag;
  return_flag = check_test(dexact_H2_N2, dtest[2][0],"Dtilde H2 derived with respect to N2 compared to bc",tol,max_diff)  || return_flag;
  return_flag = check_test(dexact_H2_CH4,dtest[2][1],"Dtilde H2 derived with respect to CH4 compared to bc",tol,max_diff) || return_flag;
  return_flag = check_test(dexact_H2_H2, dtest[2][2],"Dtilde H2 derived with respect to H2 compared to bc",tol,max_diff)  || return_flag;

  std::cout << max_diff << std::endl;

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
          tester<double>(std::string(argv[1]), "double") ||
          tester<long double>(std::string(argv[1]), "long double"));
}
