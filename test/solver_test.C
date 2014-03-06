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


void shave_string(std::string &str)
{
  while(str[0] == ' ')str.erase(0,1);
  while(str[str.size() - 1] == ' ')str.erase(str.size() - 1,1);
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
         for(unsigned int jr = 0; jr < ir; jr++)
         {
            if(mol[jr] == mol[ir])
            {
               stoi[jr]++;
               stoi.erase(stoi.begin() + ir);
               mol.erase(mol.begin() + ir);
               break;
            }
         }
      }
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
   Antioch::SplitString(molecules[1],"+",products,false);
   shave_strings(reactants);
   shave_strings(products);
   stoi_reac.resize(reactants.size(),1);
   stoi_prod.resize(products.size(),1);

   condense_molecule(stoi_reac,reactants);
   condense_molecule(stoi_prod,products);

   for(unsigned int ir = 0; ir < reactants.size(); ir++)
   {
     if( !chem_mixture.active_species_name_map().count(reactants[ir]))
     {
       skip = true;
       break;
     }
   }
   if(!skip)
   {
     for(unsigned int ip = 0; ip < products.size(); ip++)
     {
       if( !chem_mixture.active_species_name_map().count(products[ip]))
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
                              Antioch::ReactionSet<Scalar> &neutral_reaction_set)
{
   Antioch::KineticsModel::KineticsModel kineticsModel(Antioch::KineticsModel::PHOTOCHEM);
   Antioch::ReactionType::ReactionType reactionType(Antioch::ReactionType::ELEMENTARY);

   const Antioch::ChemicalMixture<Scalar>& chem_mixture = neutral_reaction_set.chemical_mixture();

   std::ifstream data(hv_file.c_str());
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
        if( !chem_mixture.active_species_name_map().count(produc[ibr][ip]))
        {
           skip[ibr] = true;
           break;
        }
      }
   }

   datas.resize(nbr - 1);
   while(!data.eof())
   {
      Scalar lambda, total;
      std::vector<Scalar> sigmas;
      sigmas.resize(nbr - 2,0.L);
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
     for(unsigned int ip = 0; ip < produc[ibr].size(); ip++)
     {
        for(unsigned int i = 0; i < stoi_prod[ibr][ip]; i++)equation += produc[ibr][ip] + " + ";
     }
     equation.erase(equation.size() - 3, 3);
     Antioch::Reaction<Scalar> * reaction = Antioch::build_reaction<Scalar>(chem_mixture.n_species(), equation, false, reactionType, kineticsModel);

     reaction->add_reactant( reac,chem_mixture.active_species_name_map().find(reac)->second,1);
     for(unsigned int ip = 0; ip < produc[ibr].size(); ip++)
     {
        reaction->add_product( produc[ibr][ip],chem_mixture.active_species_name_map().find(produc[ibr][ip])->second,stoi_prod[ibr][ip]);
     }

     std::vector<Scalar> dataf = datas[0];
     for(unsigned int i = 0; i < datas[ibr].size(); i++)
     {
        dataf.push_back(datas[ibr][i]);
     }
     Antioch::KineticsType<Scalar, std::vector<Scalar> > * rate  = Antioch::build_rate<Scalar,std::vector<Scalar> >(dataf,kineticsModel); //kinetics rate
     reaction->add_forward_rate(rate);

     neutral_reaction_set.add_reaction(reaction);
   }

}

template<typename Scalar, typename VectorScalar>
void fill_neutral_reactions_elementary(const std::string &neutral_reactions_file,
                                       const std::string &N2_hv_file,
                                       const std::string &CH4_hv_file,
                                       Antioch::ReactionSet<Scalar> &neutral_reaction_set)
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
      std::vector<std::string> reactants;
      std::vector<std::string> products;

      bool skip(false);
      std::string equation;
      std::vector<unsigned int> stoi_reac; 
      std::vector<unsigned int> stoi_prod;
      parse_equation(reactants,products,line,skip,chem_mixture,equation,stoi_reac,stoi_prod);

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
      if(dataf[1] == 0.) //Arrhenius
      {
         kineticsModel = Antioch::KineticsModel::ARRHENIUS;
      }else
      {
        dataf.push_back(std::atof(str_data[1].c_str())); //beta
      }
      dataf.push_back(std::atof(str_data[2].c_str()));//Ea

      if(kineticsModel == Antioch::KineticsModel::KOOIJ)dataf.push_back(1.); //Tref
      dataf.push_back(Antioch::Constants::R_universal<Scalar>()*1e-3); //scale (R in J/mol/K)


      Antioch::KineticsType<Scalar, VectorScalar>* rate = Antioch::build_rate<Scalar,VectorScalar>(dataf,kineticsModel); //kinetics rate
      Antioch::Reaction<Scalar> * reaction = Antioch::build_reaction<Scalar>(chem_mixture.n_species(), equation, false, reactionType, kineticsModel);
      reaction->add_forward_rate(rate);

      for(unsigned int ir = 0; ir < reactants.size(); ir++)
      {
        reaction->add_reactant( reactants[ir],chem_mixture.active_species_name_map().find( reactants[ir] )->second,stoi_reac[ir]);
      }
      for(unsigned int ip = 0; ip < products.size(); ip++)
      {
        reaction->add_product( products[ip],chem_mixture.active_species_name_map().find( products[ip] )->second,stoi_prod[ip]);
      }
      neutral_reaction_set.add_reaction(reaction);
   }
   data.close();
//now the photochemical ones
   read_photochemistry_reac(N2_hv_file, "N2", neutral_reaction_set);
   read_photochemistry_reac(CH4_hv_file, "CH4", neutral_reaction_set);
}

template<typename Scalar, typename VectorScalar>
void fill_neutral_reactions_falloff(const std::string &neutral_reactions_file,
                                    Antioch::ReactionSet<Scalar> &neutral_reaction_set)
{
//Lindemann
   std::ifstream data(neutral_reactions_file.c_str());
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
      std::vector<std::string> reactants;
      std::vector<std::string> products;

      bool skip(false);
      std::string equation;
      std::vector<unsigned int> stoi_reac; 
      std::vector<unsigned int> stoi_prod;
      parse_equation(reactants,products,line,skip,chem_mixture,equation,stoi_reac,stoi_prod);
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
      if(str_data[1] == "0" && str_data[4] == "0")
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
      dataf1.push_back(Antioch::Constants::R_universal<Scalar>()*1e-3); //scale (R in J/mol/K)
      dataf2.push_back(Antioch::Constants::R_universal<Scalar>()*1e-3); //scale (R in J/mol/K)

      Antioch::KineticsType<Scalar, VectorScalar>* rate1 = Antioch::build_rate<Scalar,VectorScalar>(dataf1,kineticsModel); //kinetics rate
      Antioch::KineticsType<Scalar, VectorScalar>* rate2 = Antioch::build_rate<Scalar,VectorScalar>(dataf2,kineticsModel); //kinetics rate
      Antioch::Reaction<Scalar> * reaction = Antioch::build_reaction<Scalar>(chem_mixture.n_species(), equation, false, reactionType, kineticsModel);
      reaction->add_forward_rate(rate1);
      reaction->add_forward_rate(rate2);
      
      for(unsigned int ir = 0; ir < reactants.size(); ir++)
      {
        reaction->add_reactant( reactants[ir],chem_mixture.active_species_name_map().find( reactants[ir] )->second,stoi_reac[ir]);
      }
      for(unsigned int ip = 0; ip < products.size(); ip++)
      {
        reaction->add_product( products[ip],chem_mixture.active_species_name_map().find( products[ip] )->second,stoi_prod[ip]);
      }
      neutral_reaction_set.add_reaction(reaction);

   }
   data.close();
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

template <typename Scalar>
void fill_molar_frac(const std::vector<std::string> &neutrals, std::vector<Scalar> &molar_frac, const std::string &file, const std::string &root_input)
{
   std::ifstream frac((root_input + file).c_str());
   std::string line;
   getline(frac,line);//title
   molar_frac.resize(neutrals.size(),0.L);
   while(!frac.eof())
   {
      std::string neutral;
      Scalar frac_600, frac_1050, Dfrac_1050;
      frac >> neutral >> frac_600 >> frac_1050 >> Dfrac_1050;
      for(unsigned int s = 0; s < neutrals.size(); s++)
      {
         if(neutral == neutrals[s])
         {
            molar_frac[s] = frac_600;
            break;
         }
      }
   }
   frac.close(); 
}

template<typename Scalar>
void read_flyby_infos(const std::vector<std::string> &neutrals, Scalar &dens_tot, std::vector<Scalar> &molar_frac, 
                      Scalar &chi, Scalar &K0, const std::string &file_flyby, const std::string &root_input)
{
   std::ifstream flyby(file_flyby.c_str());
   std::string line;
   while(getline(flyby,line))
   {
      if(line.substr(0,line.find(':')) == "zenith angle")
      {
         line = line.substr(line.find(':') + 1,std::string::npos);
         line = line.substr(0,line.rfind(' '));
         shave_string(line);
         chi = std::atof(line.c_str()); // deg, alright
      }else if(line.substr(0,line.find(':')) == "lower boundary total density")
      {
         line = line.substr(line.find(':') + 1,std::string::npos);
         line = line.substr(0,line.rfind(' '));
         shave_string(line);
         dens_tot = std::atof(line.c_str()); // cm-3
      }else if(line.substr(0,line.find(':')) == "K0")
      {
         line = line.substr(line.find(':') + 1,std::string::npos);
         line = line.substr(0,line.rfind(' '));
         shave_string(line);
         K0 = std::atof(line.c_str()) * 1e-4; //cm2/s -> m2/s
      }else if(line.substr(0,line.find(':')) == "file mix")
      {
         line = line.substr(line.find(':') + 1,std::string::npos);
         shave_string(line);
         fill_molar_frac(neutrals,molar_frac,line,root_input);
      }
   }

  flyby.close();
}

void read_neutral_system(std::vector<std::string> &neutrals, std::vector<std::string> &ions, const std::string &file_neutrals)
{
  std::ifstream neu_spec(file_neutrals.c_str());
  while(!neu_spec.eof())
  {
     std::string tmp;
     neu_spec >> tmp;
     if(tmp.empty())continue;
     neutrals.push_back(tmp);
  }
  neu_spec.close();
  shave_strings(neutrals);
  ions = neutrals;
}

template <typename Scalar>
void read_neutral_characteristics(const std::vector<std::string> & neutrals,
                                 std::vector<Scalar> &tc, 
                                 std::vector<std::vector<std::vector<Scalar> > > & bin_diff_data,
                                 std::vector<std::vector< Planet::DiffusionType> > & bin_diff_model, 
                                 const std::string & file_neutral_charac)
{
  for(unsigned int m = 0; m < bin_diff_data.size(); m++) //N2, then CH4
  {
      bin_diff_data[m].resize(neutrals.size());
      bin_diff_model[m].resize(neutrals.size(),Planet::DiffusionType::NoData); // no data by default
  }

  tc.resize(neutrals.size(),0.L);

  std::ifstream neu(file_neutral_charac.c_str());
  std::string line;
  getline(neu,line); //first line
  while(!neu.eof())
  {
    std::string name; 
    Scalar mass,AN2,sN2,ACH4,sCH4,alpha,enth,hsr,mr;
    neu >> name >> mass >> AN2 >> sN2 >> ACH4 >> sCH4 >> alpha >> enth >> hsr >> mr;
    for(unsigned int s = 0; s < neutrals.size(); s++)
    {
       if(name == neutrals[s])
       {
          tc[s] = alpha; // no unit
          bin_diff_model[0][s] = Planet::DiffusionType::Wilson; //N2
          bin_diff_model[1][s] = Planet::DiffusionType::Wilson; //CH4
          bin_diff_data[0][s].push_back(AN2 * 1e-4);  //N2 - A  -- cm2/s -> m2/s
          bin_diff_data[0][s].push_back(sN2);  //N2 - s  -- no unit
          bin_diff_data[1][s].push_back(ACH4 * 1e-4); //CH4 - A  -- cm2/s -> m2/s
          bin_diff_data[1][s].push_back(sCH4); //CH4 - s  -- no unit
          break;
       }
    }
  }
  neu.close();
  
}

template <typename Scalar>
int tester(const std::string &input_T,const std::string & input_hv, 
           const std::string &input_reactions_elem, const std::string &input_reactions_fall, 
           const std::string &input_N2, const std::string &input_CH4,
           const std::string &file_neutrals, const std::string &file_flyby,
           const std::string &file_neutral_charac,
           const std::string &root_input)
{

///////////////////////////////////////////////////
// reading part, will evolve:
//  get a nice one for diffusion & solar flux,
//  Antioch will take care of chemistry
//////////////////////////////////////////////////

//description
  std::vector<std::string> neutrals;
  std::vector<std::string> ions;

// medium for the diffusion
// decided by user, should have consistent input file
  std::vector<std::string> medium;
  medium.push_back("N2");
  medium.push_back("CH4");

//densities
  std::vector<Scalar> molar_frac;
  Scalar dens_tot;

//zenith angle
  Scalar chi;

//photon flux
  std::vector<Scalar> lambda_hv,phy1AU;

////cross-section
  std::vector<Scalar> lambda_N2,sigma_N2;
  std::vector<Scalar> lambda_CH4, sigma_CH4;

//binary diffusion
  std::vector<std::vector<std::vector<Scalar> > > bin_diff_data;    // (medium,species)[]
  std::vector<std::vector< Planet::DiffusionType> > bin_diff_model; // (medium,species)

  bin_diff_data.resize(2);
  bin_diff_model.resize(2);

//thermal coefficient
  std::vector<Scalar> tc;

//eddy
  Scalar K0;

//temperature
  std::vector<Scalar> T0,Tz;

////////////////////////////////////
// now filling with input files
/////////////////////////////////////

  read_neutral_system(neutrals, ions, file_neutrals); //ionic system in neutral, here ions = neutrals, another method is necessary to add ions (not here)

/* Flyby defines:
 * - zenith angle (chi)
 * - eddy coefficient (K0)
 * - lower boundary densities (dens_tot)
 * - file where to find mixing ratios (molar_frac)
 */
  read_flyby_infos(neutrals,dens_tot,molar_frac,chi,K0,file_flyby, root_input);

/* neutral characteristics:
 *  - thermal coeff
 *  - bimolecular diffusion
 */
  read_neutral_characteristics(neutrals,tc,bin_diff_data,bin_diff_model, file_neutral_charac);

  std::vector<Antioch::Species> spec;
  std::vector<std::vector<Planet::BinaryDiffusion<Scalar> > > bin_diff_coeff;
  bin_diff_coeff.resize(2);

/* solar flux
 * here only N2 and CH4 absorb
 */
  read_hv_flux<Scalar>(lambda_hv,phy1AU,input_hv);
  read_crossSection<Scalar>(input_N2, 3,lambda_N2, sigma_N2);
  read_crossSection<Scalar>(input_CH4,9,lambda_CH4,sigma_CH4);

/* temperature */
  read_temperature<Scalar>(T0,Tz,input_T);


//altitudes, boundaries def, either input file or user-defined, here is fine
  Scalar zmin(600.),zmax(1400.);

/************************
 * first level
 ************************/

//neutrals
  Antioch::ChemicalMixture<Scalar> neutral_species(neutrals); 
 // for diffusion, need the neutral species map
  for(unsigned int s = 0; s < neutrals.size(); s++)
  {
     spec.push_back(neutral_species.species_name_map().at(neutrals[s]));
  }
  for(unsigned int n = 0; n < 2; n++)
  {
    bin_diff_coeff[n].resize(neutrals.size());
    for(unsigned int s = 0; s < neutrals.size(); s++)
    {
      bin_diff_coeff[n][s] = Planet::BinaryDiffusion<Scalar>( spec[n], spec[s], bin_diff_data[n][s][0], bin_diff_data[n][s][1], bin_diff_model[n][s]);
    }
  }

//ions
  Antioch::ChemicalMixture<Scalar> ionic_species(ions); 

//chapman
  Planet::Chapman<Scalar> chapman(chi);

//temperature
  Planet::AtmosphericTemperature<Scalar, std::vector<Scalar> > temperature(T0, T0, Tz, Tz);

/************************
 * second level
 ************************/

//photon opacity
  Planet::PhotonOpacity<Scalar,std::vector<Scalar> > tau(chapman);
/* here only N2 and CH4 absorb */
  tau.add_cross_section(lambda_N2,  sigma_N2,  Antioch::Species::N2,  neutral_species.active_species_name_map().at("N2"));
  tau.add_cross_section(lambda_CH4, sigma_CH4, Antioch::Species::CH4, neutral_species.active_species_name_map().at("CH4"));
  tau.update_cross_section(lambda_hv);

//reaction sets
  Antioch::ReactionSet<Scalar> neutral_reaction_set(neutral_species);
  Antioch::ReactionSet<Scalar> ionic_reaction_set(ionic_species);

// here read the reactions, Antioch will take care of it once hdf5, no ionic reactions there
  fill_neutral_reactions_elementary<Scalar,std::vector<Scalar> > //Kooij / Arrhenius + photochem
                (input_reactions_elem,input_N2,input_CH4,neutral_reaction_set);

  fill_neutral_reactions_falloff<Scalar,std::vector<Scalar> >
                (input_reactions_fall,neutral_reaction_set);

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

//photon evaluator
  Planet::PhotonEvaluator<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > photon(tau,composition);
  photon.set_photon_flux_at_top(lambda_hv, phy1AU, Planet::Constants::Saturn::d_Sun<Scalar>());
  neutral_reaction_set.set_particle_flux(photon.photon_flux_ptr()); // reactions know the solar flux

//molecular diffusion
  Planet::MolecularDiffusionEvaluator<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > molecular_diffusion(bin_diff_coeff,composition,temperature);
  molecular_diffusion.set_medium_species(medium);

//eddy diffusion
  Planet::EddyDiffusionEvaluator<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > eddy_diffusion(composition,K0);

/************************
 * fourth level
 ************************/

//full diffusion
  Planet::DiffusionEvaluator<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > diffusion(molecular_diffusion,eddy_diffusion,composition,temperature);

//full chemistry
  Planet::AtmosphericKinetics<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > kinetics(neutral_kinetics, ionic_kinetics, temperature, photon, composition);

/**************************
 * fifth level
 **************************/

  Planet::PlanetPhysicsHelper<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > helper(&composition,&kinetics,&diffusion);

/************************
 * solver here, Paul work here
 ************************/

//  Planet::PlanetPhysics<Scalar,std::vector<Scalar>,std::vector<std::vector<Scalar> > > solver();

  int return_flag(0);

/* how the helper works:
 *  - void helper.first_guess(VectorStateType & molar_concentrations, const StateType &z) const;     first guess at altitude z
 *  - void helper.lower_boundary_dirichlet(VectorStateType &molar_concentrations) const;             lower boundary, constrain on molar concentrations
 *  - void helper.upper_boundary_neumann(VectorStateType & upper_boundary, const VectorStateType & molar_concentrations) const;   time dependant upper boundary on fluxes
 *
 *  - void helper.compute(const VectorStateType &molar_concentrations, const VectorStatType &dmolar_dz, const StateType &z);  ...
 *                                    ...      takes care of everything, in PlanetPhysics loop
 *  - libMesh::Real helper.diffusion_term(s)   for species s, in PlanetPhysics loop
 *  - libMesh::Real helper.chemical_term(s)    for species s, in PlanetPhysics loop
 */

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

  return (tester<float>(std::string(argv[1]),std::string(argv[2]),std::string(argv[3]),std::string(argv[4]),std::string(argv[5]),
                        std::string(argv[6]),std::string(argv[7]),std::string(argv[8]),std::string(argv[9]),std::string(argv[10])) ||
         tester<double>(std::string(argv[1]),std::string(argv[2]),std::string(argv[3]),std::string(argv[4]),std::string(argv[5]),
                        std::string(argv[6]),std::string(argv[7]),std::string(argv[8]),std::string(argv[9]),std::string(argv[10])));
          //tester<long double>(std::string(argv[1])));
}
