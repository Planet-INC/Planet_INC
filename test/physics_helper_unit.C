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
#include "planet/planet_physics_helper.h"

//C++
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <cstdlib> //atof


template<typename Scalar>
int check_test(Scalar theory, Scalar cal, const std::string &words)
{
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 6000.L;
  if(std::abs((theory-cal)/theory) < tol)return 0;
  std::cout << std::scientific << std::setprecision(20)
            << "failed test: " << words << "\n"
            << "theory: " << theory
            << "\ncalculated: " << cal
            << "\ndifference: " << std::abs((theory-cal)/cal)
            << "\ntolerance: " << tol << std::endl;
  return 1;
}


void shave_string(std::string &str)
{
  while(str[0] == ' ')str.erase(0,1);
  while(str[str.length() - 1] == ' ')str.erase(str.length() - 1,1);
}

void shave_strings(std::vector<std::string> &stock)
{
  for(unsigned int i = 0; i < stock.size(); i++)
  {
     shave_string(stock[i]);
  }
}

void condense_molecule(std::vector<unsigned int> &stoi, std::vector<std::string> &mol)
{
      for(unsigned int ir = 1; ir < mol.size(); ir++)
      {
         for(unsigned int jt = 0, jr < ir; jr++)
         {
            if(mol[jr] == mol[ir])
            {
               stoi[jr]++;
               stoi.erase(ir,1);
               mol.erase(ir,1);
               break;
            }
         }
      }
}

template<typename Scalar>
void read_photochemistry_reac(const std::string &hv_file, const std::string &reac,
                              Antioch::ReactionSet<Scalar> &neutral_reaction_set, 
                              Antioch::ReactionSet<Scalar> &neut_reac_theo)
{
   Antioch::KineticsModel::KineticsModel kineticsModel(Antioch::KineticsModel::PHOTOCHEMISTRY);
   Antioch::ReactionType::ReactionType reactionType(Antioch::ReactionType::ELEMENTARY);

   const Antioch::ChemicalMixture<Scalar>& chem_mixture = neutral_reaction_set.chemical_mixture();

   std::ifstream data(hv_file.c_str());
   std::string line;
   getline(data,line);
   std::vector<std::string> out;
   unsigned int nbr = splitString(line,' ',out,false);
   if(nbr == 0)antioch_error();
   std::vector<std::vector<Scalar> > datas;
   std::vector<std::vector<std::string> > produc;
   product.resize(nbr - 2);
   for(unsigned int ibr = 2; ibr < nbr - 2; ibr++)
   {
      splitString(out[ibr],"/",produc[ibr - 2],false);
   }
   std::vector<std::vector<unsigned int> > stoi_prod;
   stoi_prod.resize(nbr - 2);

   for(unsigned int ibr = 0; ibr < nbr - 2; ibr++)
   {
      condense_molecule(stoi_prod[ibr],produc[ibr]);
   }

   datas.resize(nbr - 1);
   while(!data.eof())
   {
      Scalar lambda, total;
      std::vector<Scalar> sigmas;
      sigmas.resize(nbr - 2,0.L)
      data >> lambda >> total;
      for(unsigned int ibr = 0; ibr < nbr - 2; ibr++)data >> sigmas[ibr];
      datas[0].push_back(lambda);
      for(unsigned int ibr = 0; ibr < nbr - 2; ibr++)datas[ibr].push_back(sigmas[ibr]);
   }
   data.close();
   for(unsigned int ibr = 0; ibr < nbr - 2; ibr++)
   {
     if(skip[ibr])continue;
     std::string equation(reac + " -> ");
     for(unsigned int ip = 0; ip < iprod[ibr]; ip++)
     {
        equation += produc[ibr][ip] + " + ";
     }
     equation.erase(equation.size() - 4, 3);
     Antioch::Reaction<Scalar> * reaction = build_reaction(chem_mixture.n_species(), equation, false, reactionType, kineticsModel);
     Antioch::Reaction<Scalar> * reaction2 = build_reaction(chem_mixture.n_species(), equation, false, reactionType, kineticsModel);

     reaction->add_reactant( reac,chem_mixture.active_species_name_map().find(reac)->second,1);
     reaction2->add_reactant( reac,chem_mixture.active_species_name_map().find(reac)->second,1);
     for(unsigned int ip = 0; ip < iprod[ibr]; ip++)
     {
        reaction->add_product( produc[ibr][ip],chem_mixture.active_species_name_map().find(produc[ibr][ip])->second,stoi_prod[ibr][ip]);
        reaction2->add_product( produc[ibr][ip],chem_mixture.active_species_name_map().find(produc[ibr][ip])->second,stoi_prod[ibr][ip]);
     }

     std::vector<Scalar> dataf = datas[0];
     for(unsigned int i = 0; i < datas[ibr].size(); i++)
     {
        dataf.push_back(datas[ibr][i]);
     }
     Antioch::KineticsType<Scalar, VectorScalar>* rate = build_rate(dataf,kineticsModel); //kinetics rate
     Antioch::KineticsType<Scalar, VectorScalar>* rate2 = build_rate(dataf,kineticsModel); //kinetics rate
     reaction.add_forward_rate(rate);
     reaction2.add_forward_rate(rate2);

     neutral_reaction_set.add_reaction(reaction);
     neut_reac_theo.add_reaction(reaction2);
   }
}

template<typename Scalar, typename VectorScalar>
void fill_neutral_reactions(const std::string &neutral_reactions_file,
                            const std::string &N2_hv_file,
                            const std::string &CH4_hv_file,
                            Antioch::ReactionSet<Scalar> &neutral_reaction_set,
                            Antioch::ReactionSet<Scalar> &neut_reac_theo)
{
//here only simple ones: bimol Kooij/Arrhenius model
   std::ifstream data(neutral_reactions_file.c_str());
   std::string line;
   getline(data,line); //title
   Antioch::KineticsModel::KineticsModel kineticsModel(Antioch::KineticsModel::KOOIJ);
   Antioch::ReactionType::ReactionType reactionType(Antioch::ReactionType::ELEMENTARY);
   const Antioch::ChemicalMixture<Scalar>& chem_mixture = neutral_reaction_set.chemical_mixture();
   while(!data.eof())
   {
      if(!getline(data,line))break;
      std::vector<std::string> out;
      SplitString(line,";",out,false);
      if(out.size() != 2)antioch_error();
      std::string equation(out[0]);
      std::string parameters(out[1]);
///// equation
      std::vector<std::string> molecules;
      SplitString(equation,"->",molecules,false);
      if(molecules.size() != 2)antioch_error();

      std::vector<std::string> reactants;
      std::vector<std::string> products;
      SplitString(molecules[0],"+",reactants,false);
      SplitString(molecules[1],"+",products,false);
      shave_strings(reactants);
      shave_strings(products);
      std::vector<unsigned int> stoi_reac;
      std::vector<unsigned int> stoi_prod;
      stoi_reac.resize(reactants.size(),1);
      stoi_prod.resize(products.size(),1);

      condense_molecule(stoi_reac,reactants);
      condense_molecule(stoi_prod,products);

      bool skip(false);
      for(unsigned int ir = 0; ir < reactants.size(); ir++)
      {
        if( !chem_mixture.active_species_name_map().count(reactants[ir]))
        {
           skip = true;
           break;
        }
      }
      if(skip)continue;
      for(unsigned int ip = 0; ip < products.size(); ip++)
      {
        if( !chem_mixture.active_species_name_map().count(products[ip]))
        {
           skip = true;
           break;
        }
      }
      if(skip)continue;
      std::vector<Scalar> dataf;
      std::vector<std::string> str_data;
      SplitString(out[1],' ',str_data,false);
      if(str_data.size() != 4)antioch_error();
      dataf.push_back(std::atof(str_data)); //Cf
      if(dataf[1] == 0.) //Arrhenius
      {
         kineticsModel = Antioch::KineticsModel::ARRHENIUS;
      }else
      {
        dataf.push_back(std::atof(str_data[1])); //beta
      }
      dataf.push_back(std::atof(srt_data[2]));//Ea
      if(kineticsModel == Antioch::KineticsModel::KOOIJ)data.push_back(1.); //Tref
      dataf.push_back(Antioch::Constants::R_universal<Scalar>()*1e-3); //scale (R in J/mol/K)

      Antioch::KineticsType<Scalar, VectorScalar>* rate = build_rate(dataf,kineticsModel); //kinetics rate
      Antioch::Reaction<Scalar> * reaction = build_reaction(chem_mixture.n_species(), equation, false, reactionType, kineticsModel);
      reaction->add_forward_rate(rate);

      Antioch::KineticsType<Scalar, VectorScalar>* rate2 = build_rate(dataf,kineticsModel); //kinetics rate
      Antioch::Reaction<Scalar> * reaction2 = build_reaction(chem_mixture.n_species(), equation, false, reactionType, kineticsModel);
      reaction2->add_forward_rate(rate2);

      for(unsigned int ir = 0; ir < reactants.size(); ir++)
      {
        reaction->add_reactant( reactants[ir],chem_mixture.active_species_name_map().find( reactants[ir] )->second,stoi_reac[ir]);
        reaction2->add_reactant( reactants[ir],chem_mixture.active_species_name_map().find( reactants[ir] )->second,stoi_reac[ir]);
      }
      for(unsigned int ip = 0; ip < products.size(); ip++)
      {
        reaction->add_product( product[ip],chem_mixture.active_species_name_map().find( product[ip] )->second,stoi_prod[ip]);
        reaction2->add_product( product[ip],chem_mixture.active_species_name_map().find( product[ip] )->second,stoi_prod[ip]);
      }
      neutral_reaction_set.add_reaction(reaction);
      neut_reac_theo.add_reaction(reaction2);
   }
   data.close();
//now the photochemical ones
   read_photochemistry_reac(N2_photochem,neutral_reaction_set, neut_reac_theo);
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
Scalar barometry(const Scalar &zmin, const Scalar &z, const Scalar &T, const Scalar &Mm, const Scalar &botdens)
{
   return botdens * Antioch::ant_exp(-(z - zmin)/((Planet::Constants::Titan::radius<Scalar>() + z) * (Planet::Constants::Titan::radius<Scalar>() + zmin) * 1e3 *
                                             Antioch::Constants::Avogadro<Scalar>() * Planet::Constants::Universal::kb<Scalar>() * T / 
                                                        (Planet::Constants::Universal::G<Scalar>() * Planet::Constants::Titan::mass<Scalar>() * Mm))
                              );
}

template<typename Scalar>
Scalar dbarometry_dz(const Scalar &zmin, const Scalar &z, const Scalar &T, const Scalar &Mm, const Scalar &botdens)
{
   return barometry(zmin, z, T, Mm, botdens) / ((Planet::Constants::Titan::radius<Scalar>() + z) * (Planet::Constants::Titan::radius<Scalar>() + zmin) * 1e6 *
                                             Antioch::Constants::Avogadro<Scalar>() * Planet::Constants::Universal::kb<Scalar>() * T / 
                                                        (Planet::Constants::Universal::G<Scalar>() * Planet::Constants::Titan::mass<Scalar>() * Mm))
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
Scalar binary_coefficient_Massman(const Scalar &T, const Scalar &P, const Scalar &D01, const Scalar &beta)
{
   return D01 * Planet::Constants::Convention::P_normal<Scalar>() / P * Antioch::ant_pow(T/Planet::Constants::Convention::T_standard<Scalar>(),beta);
}

template<typename Scalar>
Scalar binary_coefficient_unknown(const Scalar &Dii, const Scalar &Mi, const Scalar &Mj)
{
   return (Mj < Mi)?Dii * Antioch::ant_sqrt((Mj/Mi + Scalar(1.L))/Scalar(2.L))
                   :
                    Dii * Antioch::ant_sqrt((Mj/Mi));
}

template<typename Scalar>
Scalar binary_coefficient(const Antioch::DiffusionType &model ,const Scalar &T, const Scalar &P, const Scalar &par1, const Scalar &par2, const Scalar &par3)
{
   
  switch(model)
  {
     case Antioch::DiffusionType::Massman:
     {
        return binary_coefficient_massman(T,P,par1,par2);
        break;
     }
     case Antioch::DiffusionType::Wakeham:
     {
        Scalar D01 = par1 * Antioch::ant_pow(Constants::Convention::T_standard<Scalar>(),par2);
        Scalar beta = par2;
        return binary_coefficient_massman(T,P,D01,beta);
        break;
     }
     case Antioch::DiffusionType::Wilson:
     {
        Scalar D01 = par1 * Antioch::ant_pow(Constants::Convention::T_standard<Scalar>(),par2 + Scalar(1.L)) * 
                      Constants::Universal::kb<Scalar>() /
                      Constants::Convention::P_normal<Scalar>();
        Scalar beta = par2 + Scalar(1.L);
        return binary_coefficient_massman(T,P,D01,beta);
        break;
     }
     case Antioch::DiffusionType::NoData:
     {
        return binary_coefficient_unknown(par1,par2,par3);
        break;
     }
  }

   return 0.L;
}

template<typename Scalar>
Scalar pressure(const Scalar &n, const Scalar &T)
{
   return n * 1e6L * Planet::Constants::Universal::kb<Scalar>() * T;
}

template<typename Scalar>
Scalar scale_height(const Scalar &T, const Scalar &z, const Scalar &Mm)
{
  return Planet::Constants::Universal::kb<Scalar>() * T / 
         (Planet::Constants::g<Scalar>(Planet::Constants::Titan::radius<Scalar>(),z,Planet::Constants::Titan::mass<Scalar>()) *
          Mm/Antioch::Constants::Avogadro<Scalar>());
}

template<typename Scalar, typename VectorScalar, typename TensorScalar>
void compute_diffusion(const Scalar &K, const VectorScalar &densities,
                       const VectorScalar &dn_dz, const Scalar &nTot,
                       const VectorScalar &molar_frac, const VectorScalar &Mm,
                       const Scalar &T, const Scalar &P, const Scalar &z,
                       const TensorScalar &bin_coeff_data,
                       const std::vector<std::vector<Planet::DiffusionType> > &bin_coeff_model,
                       VectorScalar &omega_theo)
{

    omega_theo.resize(densities.size());

    std::vector<std::vector<Scalar> > Dij;
    Dij.resize(2);
    Dij[0].resize(densities.size());
    Dij[1].resize(densities.size());
//mol
     for(unsigned int m = 0; m < 2; m++)
     {
       for(unsigned int s = 0; s < densities.size(); s++)
       {
         Dij[m][s] = binary_coefficient(bin_coeff_model[m][s],T,P,bin_coeff_data[m][s][0],bin_coeff_data[m][s][1],Mm[s]);
       }
     }

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
       Scalar Hs = scale_height(T,z,Mm[s] * Scalar(1e-3));

       Scalar omega_theo[s] = - Dtilde * ( dns_dz[s] /densities[s]
                                      + Scalar(1.L)/Hs 
                                      + dT_dz /T * (Scalar(1.L) + (Scalar(1.L) - molar_frac[s]) * tc[s]))
                              - K      * ( dns_dz[s] /densities[s]
                                      + Scalar(1.L)/Ha
                                      + dT_dz/T);
    }

    return;
}


template <typename Scalar>
int tester(const std::string &input_T,const std::string & input_hv, const std::string &input_reactions,
           const std::string &input_N2,const std::string &input_CH4)
{
//description
  std::vector<std::string> neutrals;
  std::vector<std::string> ions;
  neutrals.push_back("N2");
  neutrals.push_back("CH4");
  neutrals.push_back("N");
  neutrals.push_back("(1)CH2");
  neutrals.push_back("H");
//ionic system contains neutral system
  ions = neutrals;
  ions.push_back("N2+");
  Scalar MN(14.008L), MC(12.011), MH(1.008L);
  Scalar MN2 = 2.L*MN , MCH4 = MC + 4.L*MH, MCH2 = 2.L * MH + MC;
  std::vector<Scalar> Mm;
  Mm.push_back(MN2);
  Mm.push_back(MCH4);
  Mm.push_back(MN);
  Mm.push_back(MCH2);
  Mm.push_back(MH);
  std::vector<std::string> medium;
  medium.push_back("N2");
  medium.push_back("CH4");

//densities
  std::vector<Scalar> molar_frac;
  molar_frac.push_back(0.95L);
  molar_frac.push_back(0.04L);
  molar_frac.push_back(0.0006L);
  molar_frac.push_back(0.0004L);
  molar_frac.push_back(0.009L);
  molar_frac.push_back(0.L); //ion
  Scalar dens_tot(1e12L);

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

//binary diffusion * 1e-4 to m2 from cm2
  std::vector<std::vector<std::vector<Scalar> > > bin_diff_datas;
  std::vector<std::vector< Planet::DiffusionType> > bin_diff_model;

  bin_diff_datas.resize(2);
  bin_diff_model.resize(2);

  std::vector<Antioch::Species> spec;
  spec.push_back(Antioch::Species::N2);
  spec.push_back(Antioch::Species::CH4);
  spec.push_back(Antioch::Species::N);
  spec.push_back(Antioch::Species::1CH2);
  spec.push_back(Antioch::Species::H);

//N2 with ...
  bin_diff_datas[0].resize(Mm.size())
  bin_diff_data[0][0].push_back(0.1783 * 1e-4); //N2, A
  bin_diff_data[0][0].push_back(0.5); //N2, s
  bin_diff_data[0][1].push_back(0.04e-5 * 1e-4); //CH4, A
  bin_diff_data[0][1].push_back(1.76); //CH4, s
  bin_diff_data[0][2].push_back(6.234e16 * 1e-4); //N, A
  bin_diff_data[0][2].push_back(0.81); //N, s
  bin_diff_data[0][3].push_back(6.234e16 * 1e-4); //(1)CH2, A
  bin_diff_data[0][3].push_back(0.81); //(1)CH2, s
  bin_diff_data[0][4].push_back(4.87e17 * 1e-4); //H, A
  bin_diff_data[0][4].push_back(0.698); //H, s
//CH4 with..
  bin_diff_data[1][0].push_back(0.04e-5 * 1e-4); //N2, A
  bin_diff_data[1][0].push_back(1.76); //N2, s
  bin_diff_data[1][1].push_back(5.73e16 * 1e-4); //CH4, A
  bin_diff_data[1][1].push_back(0.5); //CH4, s
  bin_diff_data[1][2].push_back(1.6706e17 * 1e-4); //N, A
  bin_diff_data[1][2].push_back(0.5); //N, s
  bin_diff_data[1][3].push_back(5.9311e16 * 1e-4); //(1)CH2, A
  bin_diff_data[1][3].push_back(0.5); //(1)CH2, s
  bin_diff_data[1][4].push_back(1.6706e17 * 1e-4); //H, A
  bin_diff_data[1][4].push_back(0.5); //H, s
/// models
//N2 with ...
  bin_diff_model[0].push_back(Planet::DiffusionType::Massman); //N2
  bin_diff_model[0].push_back(Planet::DiffusionType::Wilson); //CH4
  bin_diff_model[0].push_back(Planet::DiffusionType::Wilson); //N
  bin_diff_model[0].push_back(Planet::DiffusionType::Wilson); //(1)CH2
  bin_diff_model[0].push_back(Planet::DiffusionType::Wilson); //H
//CH4 with ...
  bin_diff_model[1].push_back(Planet::DiffusionType::Wakeham); //N2
  bin_diff_model[1].push_back(Planet::DiffusionType::Wakeham); //CH4
  bin_diff_model[1].push_back(Planet::DiffusionType::Wilson); //N
  bin_diff_model[1].push_back(Planet::DiffusionType::Wilson); //(1)CH2
  bin_diff_model[1].push_back(Planet::DiffusionType::Wilson); //H

  std::vector<std::vector<Planet::BinaryDiffusion<Scalar> > > bin_diff_coeff;
  bin_diff_coeff.resize(2);
  for(unsigned int n = 0; n < 2; n++)
  {
    bin_diff_coeff[n].resize(Mm.size());
    for(unsigned int s = 0; s < Mm.size(); s++)
    {
      bin_diff_coeff[n][s] = Planet::BinaryDiffusion<Scalar>( spec[n], spec[s], bin_diff_data[n][s][0], bin_diff_data[n][s][1], bin_diff_model[n][s]);
    }
  }

//thermal coefficient
  std::vector<Scalar> tc;
  tc.push_back(0.L); //N2
  tc.push_back(0.L); //CH4
  tc.push_back(0.L); //N
  tc.push_back(0.L); //(1)CH2
  tc.push_back(-0.38L); //H

//eddy
  Scalar K0(4.3e6L * 1e-4);//cm2 -> m2

/************************
 * first level
 ************************/

//neutrals
  Antioch::ChemicalMixture<Scalar> neutral_species(neutrals); 

//ions
  Antioch::ChemicalMixture<Scalar> ionic_species(ions); 

//chapman
  Planet::Chapman<Scalar> chapman(chi);

/************************
 * second level
 ************************/

//temperature
  std::vector<Scalar> T0,Tz;
  read_temperature<Scalar>(T0,Tz,input_T);
  Planet::AtmosphericTemperature<Scalar, std::vector<Scalar> > temperature(T0, T0, Tz);

//photon opacity
  Planet::PhotonOpacity<Scalar,std::vector<Scalar> > tau(chapman);
  tau.add_cross_section(lambda_N2,  sigma_N2,  Antioch::Species::N2, neutral_species.active_species_name_map().at("N2"));
  tau.add_cross_section(lambda_CH4, sigma_CH4, Antioch::Species::CH4, neutral_species.active_species_name_map().at("CH4"));
  tau.update_cross_section(lambda_hv);

//reaction sets
  Antioch::ReactionSet<Scalar> neutral_reaction_set(neutral_species);
  Antioch::ReactionSet<Scalar> ionic_reaction_set(ionic_species);

  Antioch::ReactionSet<Scalar> neut_reac_theo(neutral_species);

  fill_neutral_reactions<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > >
                (input_reactions,input_N2,input_CH4,neutral_reaction_set,neut_reac_theo); //here only simple ones

/************************
 * third level
 ************************/

//atmospheric mixture
  Planet::AtmosphericMixture<Scalar,std::vector<Scalar> > composition(neutral_species, ionic_species, temperature);
  composition.init_composition(molar_frac,dens_tot);
  composition.set_thermal_coefficient(tc);

//kinetics evaluators
  Antioch::KineticsEvaluator<Scalar> neutral_kinetics( neutral_reaction_set, 0 );
  Antioch::KineticsEvaluator<Scalar> ionic_kinetics( ionic_reaction_set, 0 );
//theo, we have confidence in Antioch
  Antioch::KineticsEvaluator<Scalar> neutral_theo( neut_reac_theo, 0 );

/************************
 * fourth level
 ************************/

//photon evaluator
  Planet::PhotonEvaluator<Scalar,std::vector<Scalar> > photon(tau,composition);
  photon.set_photon_flux_at_top(lambda_hv, phy1AU, Planet::Constants::Saturn::d_Sun<Scalar>());
  neutral_reaction_set.set_particle_flux(photon.photon_flux_ptr());
  neut_reac_theo.set_particle_flux(photon.photon_flux_ptr());

//molecular diffusion
  Planet::MolecularDiffusionEvaluator<Scalar,std::vector<Scalar> > molecular_diffusion(bin_diff_coeff,composition,temperature);
  molecular_diffusion.set_medium_species(medium);

//eddy diffusion
  Planet::EddyDiffusionEvaluator<Scalar,std::vector<Scalar> > eddy_diffusion(composition,K0);

/**************************
 * fifth level
 **************************/

//full diffusion
  Planet::DiffusionEvaluator<Scalar,std::vector<Scalar> > diffusion(molecular_diffusion,eddy_diffusion,composition,temperature);

//full chemistry
  Planet::AtmosphericKinetics<Scalar,std::vector<Scalar> > kinetics(neutral_kinetics, ionic_kinetics, temperature, photon, composition);

/**************************
 * sixth level
 **************************/

  Planet::PlanetPhysicsHelper<Scalar,std::vector<Scalar> > helper(&kinetics,&diffusion);

/************************
 * checks
 ************************/

  molar_frac.pop_back();//get the ion outta here
  Scalar mean_M;
  Antioch::set_zero(mean_M);
  for(unsigned int s = 0; s < molar_frac.size(); s++)
  {
    mean_M += molar_frac[s] * Mm[s];
  }
  mean_M *= 1e-3;//to kg


  std::vector<Scalar> other_z;
  std::vector<std::vector<Scalar> > other_n;

  int return_flag(0);
  for(Scalar z = zmax; z >= zmin; z -= zstep)
  {

     Scalar T        = temperature.neutral_temperature(z);
     Scalar dT_dz    = temperature.dneutral_temperature_dz(z);
     Scalar nTot     = barometry(zmin,z,T,mean_M,dens_tot);
     Scalar dnTot_dz = dbarometry_dz(zmin,z,T,mean_M,dens_tot);
     Scalar P        = pressure(nTot,T);
     Scalar Ha       = scale_height(T,z,mean_M);

     std::vector<Scalar> dns_dz;
     std::vector<Scalar> densities;
     calculate_densities(densities,dns_dz,molar_frac,nTot,dnTot_dz);

     std::vector<Scalar> omega_theo;
     std::vector<Scalar> chemical_theo;
     std::vector<Scalar> dummy;
     dummy.resize(densities.size());
     chemical_theo.resize(densities.size(),0.L);

     helper.compute(densities, dns_dz, other_z, other_n, z); //to update phy for everyone

     compute_diffusion(K0 * Antioch::ant_sqrt(dens_tot/nTot), densities,
                       dn_dz, nTot, molar_frac, Mm, T, P, z, omega_theo);

     neutral_theo.compute_mole_sources(T, densities, dummy, chemical_theo);
     
     for(unsigned int s = 0; s < molar_frac.size(); s++)
     {

       return_flag = return_flag ||
                        check_test(omega_theo[s],helper.diffusion_term(s),"diffusion term of species at altitude") ||
                        check_test(chemical_theo[s],helper.chemical_term(s),"chemical term of species at altitude");
     }
     other_z.push_back(z);
     other_n.push_back(densities);
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

  return (tester<float>(std::string(argv[1]),std::string(argv[2]),std::string(argv[3])) ||
          tester<double>(std::string(argv[1]),std::string(argv[2]),std::string(argv[3])));//||
          //tester<long double>(std::string(argv[1])));
}
