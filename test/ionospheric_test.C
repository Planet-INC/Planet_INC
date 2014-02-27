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
#include "antioch/vector_utils.h"
#include "antioch/chemical_mixture.h"
#include "antioch/reaction_set.h"
#include "antioch/kinetics_evaluator.h"
#include "antioch/string_utils.h"

//Planet
#include "planet/atmospheric_steady_state.h"

//C++
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>

struct LocalReaction
{
  std::string name;
  std::vector<std::string> channels;
  std::vector<std::vector<unsigned int> > br_path;
};

void read_the_species(const std::string &file_spec, std::vector<std::string> & all_species)
{
  all_species.clear();
  std::ifstream data(file_spec.c_str());
  while(!data.eof())
  {
      std::string spec;
      data >> spec;
      if(spec.empty())continue;
      all_species.push_back(spec);
  }
  data.close();
}

template <typename Scalar>
void get_the_ions(const Antioch::ChemicalMixture<Scalar> &mixture, std::vector<Antioch::Species> & ss_species)
{
  ss_species.clear();
  for(unsigned int s = 0; s < mixture.n_species(); s++)
  {
      if( mixture.species_inverse_name_map().at(mixture.species_list()[s]).find('+') != std::string::npos)
                        ss_species.push_back(mixture.species_list()[s]);
  }
}

template <typename Scalar>
void treat_reaction(LocalReaction &reac, Antioch::ReactionSet<Scalar> &reaction_set)
{
   for(unsigned int ibr = 0; ibr < reac.channels.size(); ibr++)
   {
       std::string & line = reac.channels[ibr];
       
   }
}

template <typename Scalar>
void read_reactions(const std::string &file_reac, Antioch::ReactionSet<Scalar> &reaction_set)
{
  std::ifstream data(file_reac.c_str());
  std::string line;
  getline(data,line);
  LocalReaction cur_reac;
  while(!data.eof())
  {
     if(!getline(data,line))break;
     std::vector<std::string> out;
     Antioch::SplitString(line,";",out,false);
     std::string name = out[0].substr(0,out[0].find(' '));
     std::vector<std::string> branch;
     Antioch::SplitString(name,".",branch,false);
     std::vector<unsigned int> br;
     for(unsigned int i = 2; i < branch.size(); i++)
     {
        std::stringstream b(branch[i]);
        unsigned int inode;
        b >> inode;
        br.push_back(inode);
     }
     if(cur_reac.name.empty())
     {
        cur_reac.name = branch[1];
        cur_reac.br_path.push_back(br);
     }else if(cur_reac.name == branch[1])
     {
        cur_reac.channels.push_back(line);
        cur_reac.br_path.push_back(br);
     }else
     {
        treat_reaction(cur_reac,reaction_set);

        cur_reac.br_path.clear();
        cur_reac.channels.clear();

        cur_reac.name = branch[1];
        cur_reac.channels.push_back(line);
        cur_reac.br_path.push_back(br);
     }
  }
  data.close();
}

template <typename Scalar>
void prepare_the_ionosphere(const std::string &file_neu_conc, std::vector<Scalar> &concentrations,
                            const Antioch::ChemicalMixture<Scalar> &mixture)
{
  concentrations.resize(mixture.n_species(),-1.L);
  std::ifstream data(file_neu_conc.c_str());
  Scalar ntot;
  std::string name;
  data >> name >> ntot;
  while(!data.eof())
  {
     Scalar xmol;
     data >> name >> xmol;
     concentrations[mixture.active_species_name_map().at(name)] = ntot * xmol;
  }
  data.close();
}

template <typename Scalar>
int tester(const std::string & file_spec, const std::string &file_reac, const std::string &file_neu_conc)
{
 // first, the species
  std::vector<std::string> all_species;
  read_the_species(file_spec,all_species);

  Antioch::ChemicalMixture<Scalar> mixture(all_species);

// then, the ions
  std::vector<Antioch::Species>  ss_species;
  get_the_ions(mixture,ss_species);

  Antioch::ReactionSet<Scalar> reaction_set(mixture);

// then, the ionospheric reactions
  read_reactions(file_reac,reaction_set);
  Antioch::KineticsEvaluator<Scalar> reactions_system(reaction_set,0);
  std::vector<Scalar> molar_concentrations, molar_sources;

// now the concentrations
  prepare_the_ionosphere(file_neu_conc,molar_concentrations,mixture);

//solve
  Scalar T(200.);
  Planet::AtmosphericSteadyState solver;
  solver.steady_state(reactions_system, ss_species, mixture, T, molar_concentrations, molar_sources);

  int return_flag(1);

  return return_flag;

}

int main(int argc, char** argv)
{
  // Check command line count.
  if( argc < 4 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify input files." << std::endl;
      antioch_error();
    }

  return (tester<float>(std::string(argv[1]),std::string(argv[2]),std::string(argv[3])) ||
          tester<double>(std::string(argv[1]),std::string(argv[2]),std::string(argv[3])) ||
          tester<long double>(std::string(argv[1]),std::string(argv[2]),std::string(argv[3])));

}
