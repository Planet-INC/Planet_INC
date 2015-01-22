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
#include "antioch/string_utils.h"
#include "antioch/physical_constants.h"
#include "antioch/sigma_bin_converter.h"
#include "antioch/kinetics_parsing.h"
#include "antioch/reaction_parsing.h"
#include "antioch/vector_utils.h"

//Planet
#include "planet/diffusion_evaluator.h"
#include "planet/planet_constants.h"
#include "planet/planet_physics_helper.h"
#include "planet/planet_physics_evaluator.h"

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
  Scalar test = (theory-cal);
  if(theory != 0.)test = std::abs(test/cal);
  if(test < tol)return 0;
  std::cout << std::scientific << std::setprecision(20)
            << "failed test: " << words << "\n"
            << "theory: " << theory
            << "\ncalculated: " << cal
            << "\ndifference: " << std::abs((theory-cal)/cal)
            << "\ntolerance: " << tol << std::endl;
  return 1;
}


void shave_string(std::string& str)
{
  // Trim from the left
  str.erase(str.begin(), std::find_if(str.begin(), str.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));

  // Trim from the right
  str.erase(std::find_if(str.rbegin(), str.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), str.end());

  return;
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
//first we suppress well species
    std::vector<std::string> mol_cpy(mol);
    for(unsigned int imol = 0; imol < mol_cpy.size(); imol++)
    {
        //if we find it in the copy
        if(mol_cpy[imol] == "well" ||
           mol_cpy[imol] == "CxHyNz+")
        {
          //look for and delete it in the original
           for(unsigned int jmol = 0; jmol < mol.size(); jmol++)
           {
              if(mol[jmol] == "well" || mol[jmol] == "CxHyNz+")
              {
                mol.erase(mol.begin() + jmol);
                break;
              }
           }
        }
    }

//now count similar molecules
    bool again(true);
    while(again)
    {
// for every time we find a couple, we start over again,
// to be sure nothing weird with the indices happens
      again = false;
      for(unsigned int ir = 1; ir < mol.size(); ir++)//end of vector
        {
          for(unsigned int jr = 0; jr < ir; jr++)//compared to beginning
            {
              if(mol[jr] == mol[ir])
                {
                  stoi[jr]++;
                  stoi.erase(stoi.begin() + ir);
                  mol.erase(mol.begin() + ir);
                  again = true;
                  break;
                }
            }
            if(again)break;
        }
    }

    return;
}

template <typename Scalar>
void parse_equation(std::vector<std::string> &reactants, std::vector<std::string> &products, std::string &line,
                     bool &skip, const Antioch::ChemicalMixture<Scalar>& chem_mixture, std::string &equation,
                     std::vector<unsigned int> &stoi_reac, std::vector<unsigned int> &stoi_prod)
{
   std::vector<std::string> out;
   Antioch::SplitString(line,";",out,false);
   if(out.size() != 2)antioch_error();
   equation = out[0];
   std::string parameters(out[1]);
///// equation
   std::vector<std::string> molecules;
   Antioch::SplitString(equation,"->",molecules,false);
   if(molecules.size() != 2)antioch_error();

   Antioch::SplitString(molecules[0],"+",reactants,false);
   if(reactants.empty())reactants.push_back(molecules[0]);
   Antioch::SplitString(molecules[1],"+",products,false);
   if(products.empty())products.push_back(molecules[1]);
   shave_strings(reactants);
   shave_strings(products);
   stoi_reac.resize(reactants.size(),1);
   stoi_prod.resize(products.size(),1);

   condense_molecule(stoi_reac,reactants);
   condense_molecule(stoi_prod,products);

   for(unsigned int ir = 0; ir < reactants.size(); ir++)
   {
     if( !chem_mixture.species_name_map().count(reactants[ir]))
     {
       skip = true;
       break;
     }
   }
   if(!skip)
   {
     for(unsigned int ip = 0; ip < products.size(); ip++)
     {
       if( !chem_mixture.species_name_map().count(products[ip]))
       {
         skip = true;
         break;
       }
     }
   }
   line.erase(0,line.find(';') + 1);
}

template<typename Scalar>
void read_photochemistry_reac(const std::string &hv_file, const std::string &reac,
                              Antioch::ReactionSet<Scalar> &neutral_reaction_set, 
                              Antioch::ReactionSet<Scalar> &neut_reac_theo)
{
   Antioch::KineticsModel::KineticsModel kineticsModel(Antioch::KineticsModel::PHOTOCHEM);
   Antioch::ReactionType::ReactionType reactionType(Antioch::ReactionType::ELEMENTARY);

   const Antioch::ChemicalMixture<Scalar>& chem_mixture = neutral_reaction_set.chemical_mixture();

   std::ifstream data(hv_file.c_str());
   if( !data )
     {
       std::cerr << "Could not open file " << hv_file << std::endl;
       antioch_error();
     }
   std::string line;
   getline(data,line);
   std::vector<std::string> out;
   Antioch::SplitString(line," ",out,false);
   unsigned int nbr = out.size();
   if(nbr == 0)antioch_error();

   std::vector<std::vector<Scalar> > datas;
   std::vector<std::vector<std::string> > produc;
   std::vector<bool > skip;
   skip.resize(nbr,false);
   produc.resize(nbr - 2);
   for(unsigned int ibr = 2; ibr < nbr; ibr++)
   {
      Antioch::SplitString(out[ibr],"/",produc[ibr - 2],false);
      if(produc[ibr - 2].empty())produc[ibr - 2].push_back(out[ibr]);
   }
   std::vector<std::vector<unsigned int> > stoi_prod;
   stoi_prod.resize(nbr - 2);

   for(unsigned int ibr = 0; ibr < nbr - 2; ibr++)
   {
      stoi_prod[ibr].resize(produc[ibr].size(),1);
      condense_molecule(stoi_prod[ibr],produc[ibr]);
      for(unsigned int ip = 0; ip < produc[ibr].size(); ip++)
      {
        if( !chem_mixture.species_name_map().count(produc[ibr][ip]))
        {
           skip[ibr] = true;
           break;
        }
      }
   }

   datas.resize(nbr - 1);
   while(!data.eof())
   {
      Scalar lambda(-1), total;
      std::vector<Scalar> sigmas;
      sigmas.resize(nbr - 2,0.L);
      data >> lambda >> total;
      if(!data.good())break;
      for(unsigned int ibr = 0; ibr < nbr - 2; ibr++)data >> sigmas[ibr]; //only br here
      datas[0].push_back(lambda);
      for(unsigned int ibr = 0; ibr < nbr - 2; ibr++)datas[ibr + 1].push_back(sigmas[ibr]);// + \lambda
   }
   data.close();

   for(unsigned int ibr = 0; ibr < nbr - 2; ibr++)
   {
     if(skip[ibr])continue;
     std::string equation(reac + " -> ");
     for(unsigned int ip = 0; ip < produc[ibr].size(); ip++)
     {
        for(unsigned int i = 0; i < stoi_prod[ibr][ip]; i++)equation += produc[ibr][ip] + " + ";
     }
     equation.erase(equation.size() - 3, 3);
     Antioch::Reaction<Scalar> * reaction  = Antioch::build_reaction<Scalar>(chem_mixture.n_species(), equation, false, reactionType, kineticsModel);
     Antioch::Reaction<Scalar> * reaction2 = Antioch::build_reaction<Scalar>(chem_mixture.n_species(), equation, false, reactionType, kineticsModel);

     reaction->add_reactant( reac,chem_mixture.species_name_map().find(reac)->second,1);
     reaction2->add_reactant( reac,chem_mixture.species_name_map().find(reac)->second,1);
     for(unsigned int ip = 0; ip < produc[ibr].size(); ip++)
     {
        reaction->add_product( produc[ibr][ip],chem_mixture.species_name_map().find(produc[ibr][ip])->second,stoi_prod[ibr][ip]);
        reaction2->add_product( produc[ibr][ip],chem_mixture.species_name_map().find(produc[ibr][ip])->second,stoi_prod[ibr][ip]);
     }

     std::vector<Scalar> dataf;
     dataf.resize(datas[0].size() + datas[ibr].size(),0.);
     int istep(1);
     int start(0);
     if(datas[0].back() < datas[0].front())
     {
       istep = -1;
       start = datas[0].size() - 1;
     }
     unsigned int j(0);
     for(int i = start; i < (int)datas[0].size() && i > -1; i += istep)
     {
       dataf[j] = datas[ibr + 1][i];
       j++;
     }
     for(int i = start; i < (int)datas[0].size() && i > -1; i += istep)
     {
       dataf[j] = datas[0][i];
       j++;
     }

     Antioch::KineticsType<Scalar, std::vector<Scalar> > * rate  = Antioch::build_rate<Scalar,std::vector<Scalar> >(dataf,kineticsModel); //kinetics rate
     Antioch::KineticsType<Scalar, std::vector<Scalar> > * rate2 = Antioch::build_rate<Scalar,std::vector<Scalar> >(dataf,kineticsModel); //kinetics rate
     reaction->add_forward_rate(rate);
     reaction2->add_forward_rate(rate2);

     neutral_reaction_set.add_reaction(reaction);
     neut_reac_theo.add_reaction(reaction2);
   }

}

template<typename Scalar, typename VectorScalar>
void fill_neutral_reactions_falloff(const std::string &neutral_reactions_file,
                                    Antioch::ReactionSet<Scalar> &neutral_reaction_set,
                                    Antioch::ReactionSet<Scalar> &theo_neutral_reaction_set)
{
//Lindemann
   std::ifstream data(neutral_reactions_file.c_str());
   if( !data )
     {
       std::cerr << "Could not open file " << neutral_reactions_file << std::endl;
       antioch_error();
     }
   std::string line;
   getline(data,line); //title
   getline(data,line); //title
   getline(data,line); //title
   Antioch::KineticsModel::KineticsModel kineticsModel(Antioch::KineticsModel::KOOIJ);
   Antioch::ReactionType::ReactionType reactionType(Antioch::ReactionType::LINDEMANN_FALLOFF);
   const Antioch::ChemicalMixture<Scalar>& chem_mixture = neutral_reaction_set.chemical_mixture();
   while(!data.eof())
   {
      if(!getline(data,line))break;
      if(!data.good())break;
      std::vector<std::string> reactants;
      std::vector<std::string> products;

      bool skip(false);
      std::string equation;
      std::vector<unsigned int> stoi_reac; 
      std::vector<unsigned int> stoi_prod;
      parse_equation(reactants,products,line,skip,chem_mixture,equation,stoi_reac,stoi_prod);

      kineticsModel = Antioch::KineticsModel::KOOIJ;

      if(skip)continue;
      VectorScalar dataf1,dataf2;
      std::vector<std::string> str_data;
      Antioch::SplitString(line," ",str_data,false);
      if(str_data.size() != 7)
      {
        std::cerr << "data are badly shaped, need 7 numbers in this line\n"
                  << line << std::endl;
        antioch_error();
      }

      dataf1.push_back(std::atof(str_data[0].c_str()));
      dataf2.push_back(std::atof(str_data[3].c_str()));
      Scalar temp1 = std::atof(str_data[1].c_str());
      Scalar temp2 = std::atof(str_data[4].c_str());
      if(temp1 == 0 && temp2 == 0)
      {
         kineticsModel = Antioch::KineticsModel::ARRHENIUS;
      }else
      {
         dataf1.push_back(std::atof(str_data[1].c_str()));
         dataf2.push_back(std::atof(str_data[4].c_str()));
      }
      dataf1.push_back(std::atof(str_data[2].c_str()));
      dataf2.push_back(std::atof(str_data[5].c_str()));
      if(kineticsModel == Antioch::KineticsModel::KOOIJ)
      {
        dataf1.push_back(1.); //Tref
        dataf2.push_back(1.); //Tref
      }
      dataf1.push_back(Antioch::Constants::R_universal<Scalar>() * Scalar(1e-3L)); //scale (R in kJ/mol/K)
      dataf2.push_back(Antioch::Constants::R_universal<Scalar>() * Scalar(1e-3L)); //scale (R in kJ/mol/K)

      Antioch::KineticsType<Scalar, VectorScalar>* rate1 = Antioch::build_rate<Scalar,VectorScalar>(dataf1,kineticsModel); //kinetics rate
      Antioch::KineticsType<Scalar, VectorScalar>* rate2 = Antioch::build_rate<Scalar,VectorScalar>(dataf2,kineticsModel); //kinetics rate
      Antioch::Reaction<Scalar> * reaction = Antioch::build_reaction<Scalar>(chem_mixture.n_species(), equation, false, reactionType, kineticsModel);
      reaction->add_forward_rate(rate1);
      reaction->add_forward_rate(rate2);
      
      for(unsigned int ir = 0; ir < reactants.size(); ir++)
      {
        reaction->add_reactant( reactants[ir],chem_mixture.species_name_map().find( reactants[ir] )->second,stoi_reac[ir]);
      }
      for(unsigned int ip = 0; ip < products.size(); ip++)
      {
        reaction->add_product( products[ip],chem_mixture.species_name_map().find( products[ip] )->second,stoi_prod[ip]);
      }
      neutral_reaction_set.add_reaction(reaction);
/// theo
      Antioch::KineticsType<Scalar, VectorScalar>* rate21 = Antioch::build_rate<Scalar,VectorScalar>(dataf1,kineticsModel); //kinetics rate
      Antioch::KineticsType<Scalar, VectorScalar>* rate22 = Antioch::build_rate<Scalar,VectorScalar>(dataf2,kineticsModel); //kinetics rate
      Antioch::Reaction<Scalar> * reaction2 = Antioch::build_reaction<Scalar>(chem_mixture.n_species(), equation, false, reactionType, kineticsModel);
      reaction2->add_forward_rate(rate21);
      reaction2->add_forward_rate(rate22);
      
      for(unsigned int ir = 0; ir < reactants.size(); ir++)
      {
        reaction2->add_reactant( reactants[ir],chem_mixture.species_name_map().find( reactants[ir] )->second,stoi_reac[ir]);
      }
      for(unsigned int ip = 0; ip < products.size(); ip++)
      {
        reaction2->add_product( products[ip],chem_mixture.species_name_map().find( products[ip] )->second,stoi_prod[ip]);
      }
      theo_neutral_reaction_set.add_reaction(reaction2);

   }
   data.close();
}

template<typename Scalar, typename VectorScalar>
void fill_neutral_reactions(const std::string &neutral_reactions_elem,
                            const std::string &neutral_reactions_fall,
                            const std::string &N2_hv_file,
                            const std::string &CH4_hv_file,
                            Antioch::ReactionSet<Scalar> &neutral_reaction_set,
                            Antioch::ReactionSet<Scalar> &neut_reac_theo)
{
//here only simple ones: bimol Kooij/Arrhenius model
   std::ifstream data(neutral_reactions_elem.c_str());
   if( !data )
     {
       std::cerr << "Could not open file " << neutral_reactions_elem << std::endl;
       antioch_error();
     }
   std::string line;
   getline(data,line); //title
   Antioch::KineticsModel::KineticsModel kineticsModel(Antioch::KineticsModel::KOOIJ);
   Antioch::ReactionType::ReactionType reactionType(Antioch::ReactionType::ELEMENTARY);
   const Antioch::ChemicalMixture<Scalar>& chem_mixture = neutral_reaction_set.chemical_mixture();
   while(!data.eof())
   {
      if(!getline(data,line))break;
      if(!data.good())break;
      std::vector<std::string> reactants;
      std::vector<std::string> products;

      bool skip(false);

      std::string equation;
      std::vector<unsigned int> stoi_reac; 
      std::vector<unsigned int> stoi_prod;
      parse_equation(reactants,products,line,skip,chem_mixture,equation,stoi_reac,stoi_prod);

      kineticsModel = Antioch::KineticsModel::KOOIJ;

      if(skip)continue;

      VectorScalar dataf;
      std::vector<std::string> str_data;
      Antioch::SplitString(line," ",str_data,false);
      if(str_data.size() != 4)
      {
        std::cerr << "data are badly shaped, need 4 numbers in this line\n"
                  << line << std::endl;
        antioch_error();
      }
      dataf.push_back(std::atof(str_data[0].c_str())); //Cf
      Scalar temp = std::atof(str_data[1].c_str());
      Scalar temp2 = std::atof(str_data[2].c_str());
      if(temp == 0 && temp2 == 0) //constant
      {
         kineticsModel = Antioch::KineticsModel::CONSTANT;
      }else if(temp == 0) //Arrhenius
      {
         kineticsModel = Antioch::KineticsModel::ARRHENIUS;
         dataf.push_back(std::atof(str_data[2].c_str()));//Ea
         dataf.push_back(Antioch::Constants::R_universal<Scalar>() * Scalar(1e-3L)); //scale (R in kJ/mol/K)
      }else
      {
        dataf.push_back(std::atof(str_data[1].c_str())); //beta
        dataf.push_back(std::atof(str_data[2].c_str()));//Ea
        dataf.push_back(1.); //Tref
        dataf.push_back(Antioch::Constants::R_universal<Scalar>() * Scalar(1e-3L)); //scale (R in kJ/mol/K)
      }

      Antioch::KineticsType<Scalar, VectorScalar>* rate = Antioch::build_rate<Scalar,VectorScalar>(dataf,kineticsModel); //kinetics rate
      Antioch::Reaction<Scalar> * reaction = Antioch::build_reaction<Scalar>(chem_mixture.n_species(), equation, false, reactionType, kineticsModel);
      reaction->add_forward_rate(rate);

      Antioch::KineticsType<Scalar, VectorScalar>* rate2 = Antioch::build_rate<Scalar,VectorScalar>(dataf,kineticsModel); //kinetics rate
      Antioch::Reaction<Scalar> * reaction2 = Antioch::build_reaction<Scalar>(chem_mixture.n_species(), equation, false, reactionType, kineticsModel);
      reaction2->add_forward_rate(rate2);

      for(unsigned int ir = 0; ir < reactants.size(); ir++)
      {
        reaction->add_reactant( reactants[ir],chem_mixture.species_name_map().find( reactants[ir] )->second,stoi_reac[ir]);
        reaction2->add_reactant( reactants[ir],chem_mixture.species_name_map().find( reactants[ir] )->second,stoi_reac[ir]);
      }
      for(unsigned int ip = 0; ip < products.size(); ip++)
      {
        reaction->add_product( products[ip],chem_mixture.species_name_map().find( products[ip] )->second,stoi_prod[ip]);
        reaction2->add_product( products[ip],chem_mixture.species_name_map().find( products[ip] )->second,stoi_prod[ip]);
      }
      neutral_reaction_set.add_reaction(reaction);
      neut_reac_theo.add_reaction(reaction2);
   }
   data.close();
//now falloff
   fill_neutral_reactions_falloff<Scalar,VectorScalar>(neutral_reactions_fall, neutral_reaction_set, neut_reac_theo);
//now the photochemical ones
   read_photochemistry_reac(N2_hv_file, "N2", neutral_reaction_set, neut_reac_theo);
   read_photochemistry_reac(CH4_hv_file, "CH4", neutral_reaction_set, neut_reac_theo);
}

template<typename Scalar, typename VectorScalar = std::vector<Scalar> >
void read_temperature(VectorScalar &T0, VectorScalar &Tz, const std::string &file)
{
  T0.clear();
  Tz.clear();
  std::string line;
  std::ifstream temp(file);
  if( !temp )
    {
      std::cerr << "Could not open file " << file << std::endl;
      antioch_error();
    }
  getline(temp,line);
  while(!temp.eof())
  {
     Scalar t(-1.),tz;
     temp >> t >> tz;
     if(!temp.good())break;
     T0.push_back(t);
     Tz.push_back(tz);
  }
  temp.close();
  return;
}

template<typename Scalar, typename VectorScalar = std::vector<Scalar> >
void read_hv_flux(VectorScalar &lambda, VectorScalar &phy1AU, const std::string &file)
{
  std::string line;
  std::ifstream flux_1AU(file);
  if( !flux_1AU )
    {
      std::cerr << "Could not open file " << file << std::endl;
      antioch_error();
    }
  getline(flux_1AU,line);
  while(!flux_1AU.eof())
  {
     Scalar wv(-1.),ir,dirr;
     flux_1AU >> wv >> ir >> dirr;
     if(!flux_1AU.good())break;
     lambda.push_back(wv);//A * 10.L);//nm -> A
     phy1AU.push_back(ir);// * 1e3L * (wv*1e-9L) / (Antioch::Constants::Planck_constant<Scalar>() * 
                            //            Antioch::Constants::light_celerity<Scalar>()));//W/m2/nm -> J/s/cm2/A -> s-1/cm-2/A
  }
  flux_1AU.close();

  if(lambda.back() < lambda.front())
  {
      VectorScalar tmp_l(lambda.size());
      VectorScalar tmp_p(phy1AU.size());
      for(unsigned int i = 0; i < lambda.size(); i++)
      {
         tmp_l[lambda.size() - 1 - i] = lambda[i];
         tmp_p[lambda.size() - 1 - i] = phy1AU[i];
      }
      lambda.clear();
      phy1AU.clear();
      lambda = tmp_l;
      phy1AU = tmp_p;
  }

  return;
}

template<typename Scalar, typename VectorScalar = std::vector<Scalar> >
void read_crossSection(const std::string &file, VectorScalar &lambda, VectorScalar &sigma)
{
  std::string line;
  std::ifstream sig_f(file);
  if( !sig_f )
    {
      std::cerr << "Could not open file " << file << std::endl;
      antioch_error();
    }
  getline(sig_f,line);
  while(!sig_f.eof())
  {
     Scalar wv(-1.),sigt;
     sig_f >> wv >> sigt;
     if(!sig_f.good())break;
     if(!getline(sig_f,line))break;
     lambda.push_back(wv);//A
     sigma.push_back(sigt);//cm-2/A
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

template<typename Scalar>
Scalar dbarometry_dz(const Scalar &zmin, const Scalar &z, const Scalar &T, const Scalar &Mm, const Scalar &botdens)
{
   return barometry(zmin, z, T, Mm, botdens) / ((Planet::Constants::Titan::radius<Scalar>() + z) * (Planet::Constants::Titan::radius<Scalar>() + zmin) * 1e6 *
                                             Antioch::Constants::Avogadro<Scalar>() * Planet::Constants::Universal::kb<Scalar>() * T / 
                                                        (Planet::Constants::Universal::G<Scalar>() * Planet::Constants::Titan::mass<Scalar>() * Mm))
          * (Scalar(1.L) - Scalar(2.L) * (z - zmin)/(Planet::Constants::Titan::radius<Scalar>() + zmin));
}

template<typename Scalar, typename VectorScalar>
void calculate_densities(VectorScalar &densities, VectorScalar &sum_densities, VectorScalar &dn_dz, 
                         const VectorScalar &molar_frac, const Scalar &nTot, const Scalar &dnTot_dz, const Scalar &dz,
                         const bool compute)
{
   if(sum_densities.size() != molar_frac.size())sum_densities.resize(molar_frac.size(),0.L);
   densities.resize(molar_frac.size(),0.L);
   dn_dz.resize(molar_frac.size(),0.L);
   for(unsigned int s = 0; s < molar_frac.size(); s++)
   {
     densities[s] = molar_frac[s] * nTot;
     dn_dz[s] = molar_frac[s] * dnTot_dz;
   }

   if(compute)
   {
     for(unsigned int s = 0; s < molar_frac.size(); s++)
     {
       sum_densities[s] += densities[s] * dz;
     }
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
Scalar binary_coefficient(const Planet::DiffusionType &model ,
                          const Scalar &T, const Scalar &P, const Scalar &par1, const Scalar &par2, const Scalar &par3)
{
   
  switch(model)
  {
     case Planet::DiffusionType::Massman:
     {
        return binary_coefficient_Massman(T,P,par1,par2);
        break;
     }
     case Planet::DiffusionType::Wakeham:
     {
        Scalar D01 = par1 * Antioch::ant_pow(Planet::Constants::Convention::T_standard<Scalar>(),par2);
        Scalar beta = par2;
        return binary_coefficient_Massman(T,P,D01,beta);
        break;
     }
     case Planet::DiffusionType::Wilson:
     {
        Scalar D01 = par1 * Antioch::ant_pow(Planet::Constants::Convention::T_standard<Scalar>(),par2 + Scalar(1.L)) * 
                      Planet::Constants::Universal::kb<Scalar>() /
                      Planet::Constants::Convention::P_normal<Scalar>();
        Scalar beta = par2 + Scalar(1.L);
        return binary_coefficient_Massman(T,P,D01,beta);
        break;
     }
     case Planet::DiffusionType::NoData:
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
  return Scalar(1e-3) * Planet::Constants::Universal::kb<Scalar>() * T /  // in km
         (Planet::Constants::g<Scalar>(Planet::Constants::Titan::radius<Scalar>(),z,Planet::Constants::Titan::mass<Scalar>()) *
          Mm/Antioch::Constants::Avogadro<Scalar>());
}

template<typename Scalar, typename VectorScalar, typename TensorScalar>
void compute_diffusion(const Scalar &K, const VectorScalar &densities,
                       const VectorScalar &dns_dz, const Scalar &nTot,
                       const VectorScalar &molar_frac, const VectorScalar &Mm,
                       const Scalar &T, const Scalar &dT_dz, const Scalar &P, const Scalar &z,
                       const TensorScalar &bin_coeff_data,
                       const std::vector<std::vector<Planet::DiffusionType> > &bin_coeff_model,
                       const VectorScalar &tc, const Scalar &Ha,
                       VectorScalar &omega_theo_A, VectorScalar &omega_theo_B)
{

    omega_theo_A.resize(molar_frac.size(),0.L);
    omega_theo_B.resize(molar_frac.size(),0.L);

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

       Scalar Dt_theo = Ds / (  Scalar(1.L) - molar_frac[s] * 
                               (Scalar(1.L) - Mm[s] / M_diff)
                             );
//
       Scalar Hs = scale_height(T,z,Mm[s]);

       omega_theo_A[s] = - Dt_theo - K;
       omega_theo_B[s] = - Dt_theo * 
                                    ( Scalar(1.L)/Hs 
                                     + dT_dz /T * (Scalar(1.L) + (Scalar(1.L) - molar_frac[s]) * tc[s])
                                    )
                         - K * 
                                    ( Scalar(1.L)/Ha
                                     + dT_dz/T
                                    );
       omega_theo_A[s] *= Scalar(1e-10);
       omega_theo_B[s] *= Scalar(1e-10);
    }

    return;
}


template <typename Scalar>
int tester(const std::string &input_T,const std::string & input_hv, 
           const std::string &input_reactions_elem, const std::string &input_reactions_fall, 
           const std::string &input_reactions_photochem_N2, const std::string &input_reactions_photochem_CH4, 
           const std::string &input_N2,const std::string &input_CH4, const std::string& input_filename,
           const std::string &species_file )
{

//description
  std::vector<std::string> neutrals;
  std::vector<std::string> ions;
  neutrals.push_back("N2");
  neutrals.push_back("CH4");
  neutrals.push_back("N(4S)");
  neutrals.push_back("CH3");
  neutrals.push_back("(1)CH2");
  neutrals.push_back("(3)CH2");
  neutrals.push_back("H");
  neutrals.push_back("H2");
//ionic system contains neutral system
  ions = neutrals;
  Scalar MN(14.008e-3L), MC(12.011e-3L), MH(1.008e-3L);
  Scalar MN2  = 2.L * MN , 
         MCH4 = MC + 4.L * MH, 
         MCH2 = 2.L * MH + MC, 
         MH2  = 2.L * MH,
         MCH3 = MC + 3.L * MH;
  std::vector<Scalar> Mm;
  Mm.push_back(MN2);
  Mm.push_back(MCH4);
  Mm.push_back(MN);
  Mm.push_back(MCH3);
  Mm.push_back(MCH2);
  Mm.push_back(MCH2);
  Mm.push_back(MH);
  Mm.push_back(MH2);
  std::vector<std::string> medium;
  medium.push_back("N2");
  medium.push_back("CH4");

//densities
  std::vector<Scalar> molar_frac;
  molar_frac.push_back(0.95L);    // N2
  molar_frac.push_back(0.04L);    // CH4
  molar_frac.push_back(0.00025L); // N
  molar_frac.push_back(0.00025L); // CH3
  molar_frac.push_back(0.00025L); // (1)CH2
  molar_frac.push_back(0.00025L); // (3)CH2
  molar_frac.push_back(0.004L);   // H
  molar_frac.push_back(0.005L);   // H2
  Scalar dens_tot(1e12L);

//zenith angle
  Scalar chi(120);

//photon flux
  std::vector<Scalar> lambda_hv,phy1AU;
  read_hv_flux<Scalar>(lambda_hv,phy1AU,input_hv);
  Antioch::ParticleFlux<std::vector<Scalar> > phy_at_top;
  phy_at_top.set_abscissa(lambda_hv);
  std::vector<Scalar> flux(phy1AU.size(),0.);
  for(unsigned int il = 0; il < phy1AU.size(); il++)
  {
     flux[il] = phy1AU[il] / (Planet::Constants::Saturn::d_Sun<Scalar>() * Planet::Constants::Saturn::d_Sun<Scalar>());
  }
  phy_at_top.set_flux(flux);

  Antioch::ParticleFlux<std::vector<Scalar> > phy_at_z;
  phy_at_z.set_abscissa(lambda_hv);

////cross-section
  std::vector<Scalar> lambda_N2,sigma_N2;
  std::vector<Scalar> lambda_CH4, sigma_CH4;
  read_crossSection<Scalar>(input_N2,lambda_N2,sigma_N2);
  read_crossSection<Scalar>(input_CH4,lambda_CH4,sigma_CH4);

//altitudes
  Scalar zmin(600.),zmax(1400.),zstep(10.);

//binary diffusion * 1e-4 to m2 from cm2
  std::vector<std::vector<std::vector<Scalar> > > bin_diff_data;
  std::vector<std::vector< Planet::DiffusionType> > bin_diff_model;

  bin_diff_data.resize(2);
  bin_diff_model.resize(2);

  std::vector<unsigned int> spec;
  spec.push_back(0);
  spec.push_back(1);
  spec.push_back(2);
  spec.push_back(3);
  spec.push_back(4);
  spec.push_back(5);
  spec.push_back(6);
  spec.push_back(7);

//N2 with ...
  bin_diff_data.resize(medium.size());
  bin_diff_data[0].resize(Mm.size());
  bin_diff_data[0][0].push_back(5.9e16);   // N2, A
  bin_diff_data[0][0].push_back(0.81);             // N2, s
  bin_diff_data[0][1].push_back(7.34e16);  // CH4, A
  bin_diff_data[0][1].push_back(0.75);            // CH4, s
  bin_diff_data[0][2].push_back(6.234e16); // N, A
  bin_diff_data[0][2].push_back(0.81);            // N, s
  bin_diff_data[0][3].push_back(6.094e16); // CH3, A
  bin_diff_data[0][3].push_back(0.81);            // CH3, s
  bin_diff_data[0][4].push_back(6.234e16); // (1)CH2, A
  bin_diff_data[0][4].push_back(0.81);            // (1)CH2, s
  bin_diff_data[0][5].push_back(6.234e16); // (3)CH2, A
  bin_diff_data[0][5].push_back(0.81);            // (3)CH2, s
  bin_diff_data[0][6].push_back(4.87e17);  // H, A
  bin_diff_data[0][6].push_back(0.698);           // H, s
  bin_diff_data[0][7].push_back(1.88e17);  // H2, A
  bin_diff_data[0][7].push_back(0.82);            // H2, s
//CH4 with..
  bin_diff_data[1].resize(Mm.size());
  bin_diff_data[1][0].push_back(7.34e16);   // N2, A
  bin_diff_data[1][0].push_back(0.75);             // N2, s
  bin_diff_data[1][1].push_back(5.73e16);   // CH4, A
  bin_diff_data[1][1].push_back(0.5);              // CH4, s
  bin_diff_data[1][2].push_back(5.9311e16); // N, A
  bin_diff_data[1][2].push_back(0.5);              // N, s
  bin_diff_data[1][3].push_back(5.8247e16); // CH3, A
  bin_diff_data[1][3].push_back(0.5);              // CH3, s
  bin_diff_data[1][4].push_back(5.9311e16); // (1)CH2, A
  bin_diff_data[1][4].push_back(0.5);              // (1)CH2, s
  bin_diff_data[1][5].push_back(5.9311e16); // (3)CH2, A
  bin_diff_data[1][5].push_back(0.5);              // (3)CH2, s
  bin_diff_data[1][6].push_back(1.6706e17); // H, A
  bin_diff_data[1][6].push_back(0.5);              // H, s
  bin_diff_data[1][7].push_back(2.3e17);    // H2, A
  bin_diff_data[1][7].push_back(0.765);            // H2, s

/// models
//N2 with ...
  bin_diff_model.resize(medium.size());
  bin_diff_model[0].push_back(Planet::DiffusionType::Wilson); // N2
  bin_diff_model[0].push_back(Planet::DiffusionType::Wilson); // CH4
  bin_diff_model[0].push_back(Planet::DiffusionType::Wilson);  // N
  bin_diff_model[0].push_back(Planet::DiffusionType::Wilson);  // CH3
  bin_diff_model[0].push_back(Planet::DiffusionType::Wilson);  // (1)CH2
  bin_diff_model[0].push_back(Planet::DiffusionType::Wilson);  // (3)CH2
  bin_diff_model[0].push_back(Planet::DiffusionType::Wilson);  // H
  bin_diff_model[0].push_back(Planet::DiffusionType::Wilson);  // H2
//CH4 with ...
  bin_diff_model[1].push_back(Planet::DiffusionType::Wilson); // N2
  bin_diff_model[1].push_back(Planet::DiffusionType::Wilson);  // CH4
  bin_diff_model[1].push_back(Planet::DiffusionType::Wilson);  // N
  bin_diff_model[1].push_back(Planet::DiffusionType::Wilson);  // CH3
  bin_diff_model[1].push_back(Planet::DiffusionType::Wilson);  // (1)CH2
  bin_diff_model[1].push_back(Planet::DiffusionType::Wilson);  // (3)CH2
  bin_diff_model[1].push_back(Planet::DiffusionType::Wilson);  // H
  bin_diff_model[1].push_back(Planet::DiffusionType::Wilson);  // H2


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
  tc.push_back(0.L);    // N2
  tc.push_back(0.L);    // CH4
  tc.push_back(0.L);    // N
  tc.push_back(0.L);    // CH3
  tc.push_back(0.L);    // (1)CH2
  tc.push_back(0.L);    // (3)CH2
  tc.push_back(-0.38L); // H
  tc.push_back(-0.38L); // H2

//eddy
  Scalar K0(4.3e6L);//cm2.s-1


/************************
 * first level
 ************************/

//neutrals
  Antioch::ChemicalMixture<Scalar> neutral_species(neutrals,true,species_file); 

//ions
  Antioch::ChemicalMixture<Scalar> ionic_species(ions,true,species_file); 

//chapman
  Planet::Chapman<Scalar> chapman(chi);

//temperature
  std::vector<Scalar> T0,Tz;
  read_temperature<Scalar>(T0,Tz,input_T);
  Planet::AtmosphericTemperature<Scalar, std::vector<Scalar> > temperature(T0, T0, Tz, Tz);

/************************
 * second level
 ************************/

//photon opacity
  Planet::PhotonOpacity<Scalar,std::vector<Scalar> > tau(chapman);
  tau.add_cross_section(lambda_N2,  sigma_N2,  0,  neutral_species.species_name_map().at("N2"));
  tau.add_cross_section(lambda_CH4, sigma_CH4, 1, neutral_species.species_name_map().at("CH4"));
  tau.update_cross_section(lambda_hv);

//reaction sets
  Antioch::ReactionSet<Scalar> neutral_reaction_set(neutral_species);
  Antioch::ReactionSet<Scalar> ionic_reaction_set(ionic_species);

  Antioch::ReactionSet<Scalar> neut_reac_theo(neutral_species);

  fill_neutral_reactions<Scalar,std::vector<Scalar> >
                (input_reactions_elem,input_reactions_fall,input_reactions_photochem_N2,input_reactions_photochem_CH4,neutral_reaction_set,neut_reac_theo); //here only simple ones

//atmospheric mixture
  Planet::AtmosphericMixture<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > composition(neutral_species, ionic_species, temperature);
  composition.init_composition(molar_frac,dens_tot,zmin,zmax);
  composition.set_thermal_coefficient(tc);

/************************
 * third level
 ************************/

//kinetics evaluators
  Antioch::KineticsEvaluator<Scalar> neutral_kinetics( neutral_reaction_set, 0 );
  Antioch::KineticsEvaluator<Scalar> ionic_kinetics( ionic_reaction_set, 0 );
//theo, we have confidence in Antioch
  Antioch::KineticsEvaluator<Scalar> neutral_theo( neut_reac_theo, 0 );

//photon evaluator
  Planet::PhotonEvaluator<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > photon(phy_at_top,tau,composition);


/**************************
 * fifth level
 **************************/
  GetPot input(input_filename);

  Planet::PlanetPhysicsHelper<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > helper(input);

  Planet::PlanetPhysicsEvaluator<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > evaluator(helper);

/************************
 * checks
 ************************/

  Scalar mean_M;
  Antioch::set_zero(mean_M);
  for(unsigned int s = 0; s < molar_frac.size(); s++)
  {
    mean_M += molar_frac[s] * Mm[s];
  }

  int return_flag(0);
  const Scalar nTot_virtual(dens_tot);
  std::vector<Scalar> flux_at_z(lambda_hv.size(),0.);
  std::vector<Scalar> sum_densities;

  for(Scalar z = zmax; z >= zmax; z -= zstep)
  {

     Scalar T        = temperature.neutral_temperature(z);
     Scalar dT_dz    = temperature.dneutral_temperature_dz(z);
     Scalar nTot     = barometry(zmin,z,T,mean_M,dens_tot);
     Scalar dnTot_dz = dbarometry_dz(zmin,z,T,mean_M,dens_tot);
     Scalar P        = pressure(nTot,T);
     Scalar Ha       = scale_height(T,z,mean_M);

     std::vector<Scalar> dns_dz;
     std::vector<Scalar> densities;

     calculate_densities(densities,sum_densities,dns_dz,molar_frac,nTot,dnTot_dz, zstep, (z != zmax));
 
     photon.update_photon_flux(densities,sum_densities,z,flux_at_z);
     phy_at_z.set_flux(flux_at_z);

     std::vector<Scalar> omega_theo_A;
     std::vector<Scalar> omega_theo_B;
     std::vector<Scalar> chemical_theo;
     std::vector<Scalar> dummy;
     std::vector<Scalar> virtual_densities;
     std::vector<Scalar> virtual_dns_dz;

     dummy.resize(densities.size());
     chemical_theo.resize(densities.size(),0.L);

     Antioch::KineticsConditions<Scalar> KC(T);
     for(unsigned int hv = 0; hv < helper.index_photochemistry().size(); hv++)
     {
        KC.add_particle_flux(phy_at_z,helper.index_photochemistry()[hv]);
     }

     neutral_theo.compute_mole_sources(KC, densities, dummy, chemical_theo);

     Scalar K = K0 * Antioch::ant_sqrt(dens_tot/nTot);

     compute_diffusion(K, densities,
                       dns_dz, nTot, molar_frac, Mm, T, dT_dz, P, z, 
                       bin_diff_data, bin_diff_model, 
                       tc, Ha, omega_theo_A, omega_theo_B);

     virtual_densities.resize(densities.size());
     virtual_dns_dz.resize(densities.size());
//
// The solver use rescaled densities, so we
// descale here before giving it to the evaluator
// who's rescaling
     for(unsigned int s = 0; s < molar_frac.size(); s++)
     {
         virtual_densities[s] = densities[s] / nTot_virtual;
         virtual_dns_dz[s]    = dns_dz[s] / nTot_virtual;
     }
     evaluator.compute(virtual_densities, virtual_dns_dz, z);//compute with real densities

     std::stringstream alt;
     alt << z;
     for(unsigned int s = 0; s < molar_frac.size(); s++)
     {
       std::stringstream spec;
       spec << s;
       Scalar diff_A = evaluator.diffusion_A_term(s);
       Scalar diff_B = evaluator.diffusion_B_term(s);
       Scalar chem = evaluator.chemical_term(s) * nTot_virtual;
       return_flag =  check_test(omega_theo_A[s],   diff_A, "diffusion term A of species " + spec.str() + " at altitude " + alt.str()) ||
         return_flag;
       return_flag =  check_test(omega_theo_B[s],   diff_B, "diffusion term B of species " + spec.str() + " at altitude " + alt.str()) ||
         return_flag;
       return_flag =  check_test(chemical_theo[s],chem, "chemical term of species "  + spec.str() + " at altitude " + alt.str())  ||
         return_flag;
     }

  }

  return return_flag;
}

int main(int argc, char** argv)
{
  // Check command line count.
  if( argc < 10 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify input files." << std::endl;
      antioch_error();
    }

  return (/*tester<float>(std::string(argv[1]),std::string(argv[2]),std::string(argv[3]),std::string(argv[4]),std::string(argv[5]),
                        std::string(argv[6]),std::string(argv[7]),std::string(argv[8]),std::string(argv[9]),std::string(argv[10])) ||*/ //float can't take some rate constants
          tester<double>(std::string(argv[1]),std::string(argv[2]),std::string(argv[3]),std::string(argv[4]),
                         std::string(argv[5]),std::string(argv[6]),std::string(argv[7]),std::string(argv[8]),std::string(argv[9]),std::string(argv[10]))||
          tester<long double>(std::string(argv[1]),std::string(argv[2]),std::string(argv[3]),std::string(argv[4]),
                         std::string(argv[5]),std::string(argv[6]),std::string(argv[7]),std::string(argv[8]),std::string(argv[9]),std::string(argv[10])));
}
