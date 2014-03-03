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
#include "antioch/kinetics_parsing.h"
#include "antioch/reaction_parsing.h"
#include "antioch/physical_constants.h"
#include "antioch/string_utils.h"

//Planet
#include "planet/kinetics_branching_structure.h"
#include "planet/branching_ratio_node.h"
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

template <typename Scalar>
struct LocalPdf
{
  std::string name;
  std::vector<Scalar> params;
};


void shave_string(std::string &str)
{
  while(str[0] == ' ')str.erase(0,1);
  while(str[str.size() - 1] == ' ')str.erase(str.size() - 1,1);
}

void shave_string(std::vector<std::string> &stock)
{
  for(unsigned int i = 0; i < stock.size(); i++)
  {
     shave_string(stock[i]);
  }
}


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
void fill_local_pdf(LocalPdf<Scalar> & pdf,const std::string &line)
{
   std::string pdf_name = line.substr(0,line.find("("));
   if(pdf.name.empty())pdf.name = pdf_name;
   if(pdf.name != pdf_name)antioch_error();

   std::string param = line.substr(line.find("(")+1,line.find(")") - line.find("(") - 1);
   std::vector<std::string> params;
   Antioch::SplitString(param,",",params,false);
   if(params.empty())params.push_back(param);
   for(unsigned int i = 0; i < params.size(); i++)
   {
      pdf.params.push_back(std::atof(params[i].c_str()));
   }
   if(pdf_name == "DiUn")
   {
      if(pdf.params.size() > 1)
      {
         pdf.params[0]++;
         pdf.params.resize(1);
      }
   }
   return;
}

template <typename Scalar>
void sort_pdf(LocalPdf<Scalar> &pdf)
{
    if(pdf.params.size()%2 != 0)return;

    std::vector<Scalar> data;
    data.resize(pdf.params.size());
    for(unsigned int i = 0; i < pdf.params.size()/2; i++)
    {
       data[i] = pdf.params[2*i];
       data[i+pdf.params.size()/2] = pdf.params[2*i+1];
    }

}

template <typename Scalar>
void treat_reaction(LocalReaction &reac, Antioch::ReactionSet<Scalar> &reaction_set, unsigned int n_species)
{
   Antioch::ReactionType::ReactionType typeReaction(Antioch::ReactionType::ELEMENTARY);
   Antioch::KineticsModel::KineticsModel kineticsModel(Antioch::KineticsModel::HERCOURT_ESSEN);

   Planet::KineticsBranchingStructure<Scalar> prob_reac;
   prob_reac.set_n_channels(reac.channels.size());

   std::vector<std::string> out;
   Antioch::SplitString(reac.channels[0],";",out,false);
   shave_string(out);
   std::map<std::string, LocalPdf<Scalar> > params_diri;
   std::vector<std::string> node_names;

/// first, fill prob_reac
// first br level, in out[3]
   Planet::BranchingRatioNode<Scalar> master_node;
   master_node.set_id("master");
   Planet::PDFName::PDFName pdf_type;
   pdf_type = prob_reac.pdf_map().at(out[3].substr(0,out[3].find('(')));
   master_node.set_pdf(pdf_type);
   std::vector<std::vector<Scalar> > master_data;
   master_data.resize(reac.channels.size());
   for(unsigned int ibr = 0; ibr < reac.channels.size(); ibr++)
   {
       std::string data = out[3].substr(out[3].find('(') + 1, out[3].find(')') - out[3].find('(') - 1);
       std::vector<std::string> datas;
       Antioch::SplitString(data,",",datas,false);
       master_data[ibr].resize(datas.size());
       for(unsigned int i = 0; i < datas.size(); i++)
       {
          Scalar dat;
          std::stringstream oss(datas[i]);
          oss >> dat;
          master_data[ibr][i] = dat;
       }
   }
   std::vector<Scalar> data;
   for(unsigned int n = 0; n < master_data.front().size(); n++)
   {
     for(unsigned int i = 0; i < master_data.size(); i++)
     {
        data.push_back(master_data[i][n]);
     }
   }
   master_node.pdf_object_ptr()->set_parameters(data);
   prob_reac.add_node(master_node);

   for(unsigned int ibr = 0; ibr < reac.channels.size(); ibr++)
   {

       prob_reac.add_to_br_path(ibr,"master",ibr);

       std::string & line = reac.channels[ibr];
       Antioch::SplitString(line,";",out,false);
/// equation
       std::vector<std::string> reacProd,reactants,products;
       Antioch::SplitString(out[0],"->",reacProd,false);
       Antioch::SplitString(reacProd[0]," ",reactants,false);
       Antioch::SplitString(reacProd[1]," ",products,false);


      if(ibr == 0)
      {
// k
// create necessary pdf objects
         std::vector<Planet::PDFName::PDFName> pdf_type;
         std::string loc_pdf = out[2].substr(0,out[2].find('('));
         shave_string(loc_pdf);
         pdf_type.push_back(prob_reac.pdf_map().at(loc_pdf));
         if(line.find("beta") != std::string::npos)
         {
            std::string str = out.back().substr(0,out.back().find('('));
            shave_string(str);
            str = str.substr(str.find(" ") + 1);
            shave_string(str);
            pdf_type.push_back(prob_reac.pdf_map().at(str));

         }
         prob_reac.set_k_pdf(pdf_type);

/// fill 'em
         std::string str_param = out[2].substr(out[2].find('('), out[2].find(')') - out[2].find('('));
         std::vector<std::string> k_pdf_params_str;
         Antioch::SplitString(str_param,",",k_pdf_params_str,false);

         std::vector<std::vector<Scalar> > k_pdf_params;
         k_pdf_params.push_back(std::vector<Scalar>());
         for(unsigned int i = 0; i < k_pdf_params_str.size(); i++)
         {
           std::stringstream oss(k_pdf_params_str[i]);
           Scalar p;
           oss >> p;
           k_pdf_params[0].push_back(p);
         }
         if(line.find("beta") != std::string::npos)
         {
            str_param = out.back().substr(out.back().find('('), out.back().find(')') - out.back().find('('));
            Antioch::SplitString(str_param,",",k_pdf_params_str,false);
            k_pdf_params.push_back(std::vector<Scalar>());
            for(unsigned int i = 0; i < k_pdf_params_str.size(); i++)
            {
              std::stringstream oss(k_pdf_params_str[i]);
              Scalar p;
              oss >> p;
              k_pdf_params[1].push_back(p);
            }
           
          }

          for(unsigned int i = 0; i < prob_reac.pdf_k_object().size(); i++)
          {
             prob_reac.pdf_k_object()[i]->set_parameters(k_pdf_params[i]);
          }

        } // if(ibr == 0)

        std::stringstream oss;
        oss << "master";
// levels
        unsigned int ilev(4);
        for(unsigned int ino = 0; ino < reac.br_path[ibr].size(); ino++)
        {
           oss << "." << reac.br_path[ibr][ino]; 
           std::string node_name = oss.str();
           if(!prob_reac.nodes_map().count(node_name))
           {
                node_names.push_back(node_name);
                prob_reac.create_node(node_name);
                prob_reac.node(node_name).set_n_channels(0);
           }
           prob_reac.node(node_name).set_n_channels(prob_reac.node(node_name).n_channels() + 1);
           prob_reac.add_to_br_path(ibr,node_name, prob_reac.node(node_name).n_channels() - 1);
           fill_local_pdf(params_diri[node_name],out[ilev]);
           ilev++;
        }
   }//end fill by br

   for(unsigned int inode = 0; inode < params_diri.size(); inode++) // sort nodes
   {
      sort_pdf(params_diri[node_names[inode]]);
   }

   if(out[1].find("DR") == std::string::npos)kineticsModel = Antioch::KineticsModel::CONSTANT; //type reaction

// then, create the reaction

   std::cout << "here's reaction " << reac.name << std::endl;
   std::cout << "there are " << reac.channels.size() << " channels:" << std::endl;
   for(unsigned int ibr = 0; ibr < reac.channels.size(); ibr++)
   {
std::cout << "1." << ibr << " / " << reac.channels.size() - 1 << ": " << out[0] << std::endl;
       Antioch::Reaction<Scalar>* my_rxn = Antioch::build_reaction<Scalar>(n_species, out[0], false, typeReaction, kineticsModel);

       Antioch::KineticsType<Scalar>* rate = prob_reac.create_rate_constant(ibr,kineticsModel);

       my_rxn->add_forward_rate(rate);

       reaction_set.add_reaction(my_rxn);     
std::cout << "added" << std::endl;

   }
}

template <typename Scalar>
void read_reactions(const std::string &file_reac, Antioch::ReactionSet<Scalar> &reaction_set, unsigned int n_species)
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
     line.erase(0,name.size());
     std::vector<std::string> branch;
     Antioch::SplitString(name,".",branch,false);

     std::vector<unsigned int> br;
     for(unsigned int i = 3; i < branch.size(); i++)
     {
        if(branch[i] == "0")break;
        std::stringstream b(branch[i]);
        unsigned int inode;
        b >> inode;
        br.push_back(inode);
     }
     if(cur_reac.name.empty())
     {
        cur_reac.name = branch[1];
        cur_reac.br_path.push_back(br);
        cur_reac.channels.push_back(line);
     }else if(cur_reac.name == branch[1])
     {
        cur_reac.channels.push_back(line);
        cur_reac.br_path.push_back(br);
     }else
     {
        treat_reaction(cur_reac,reaction_set, n_species);

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
std::cout << "0" << std::endl;
  read_reactions(file_reac,reaction_set, mixture.n_species());
std::cout << "1" << std::endl;
  Antioch::KineticsEvaluator<Scalar> reactions_system(reaction_set,0);
  std::vector<Scalar> molar_concentrations, molar_sources;

// now the concentrations
  prepare_the_ionosphere(file_neu_conc,molar_concentrations,mixture);

//solve
  Scalar T(200.);
  Planet::AtmosphericSteadyState solver;
//  solver.steady_state(reactions_system, ss_species, mixture, T, molar_concentrations, molar_sources);

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
