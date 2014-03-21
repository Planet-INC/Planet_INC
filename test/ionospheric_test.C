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

struct BrPath 
{
  std::vector<std::string> name;
  std::vector<unsigned int > id_channel;
};


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
    // Trim from the left
  str.erase(str.begin(), std::find_if(str.begin(), str.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));

    // Trim from the right
 str.erase(std::find_if(str.rbegin(), str.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), str.end());

 return;
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


void condense_molecule(std::vector<std::string> &molecule, std::vector<unsigned int> &stoi)
{
   for(unsigned int imol = 0; imol < molecule.size(); imol++)
   {
      for(unsigned int jmol = imol + 1; jmol < molecule.size(); jmol++)
      {
         if(molecule[imol] != molecule[jmol])continue;
         stoi[imol]++;
         stoi.erase(stoi.begin() + jmol);
         molecule.erase(molecule.begin() + jmol);
      }
   }
}


void sanity_check_reaction(LocalReaction &cur_reac)
{
  unsigned int idmax(0);
  for(unsigned int ibr = 0; ibr < cur_reac.channels.size(); ibr++)
  {
      if(cur_reac.br_path[ibr].size() > idmax)idmax = cur_reac.br_path[ibr].size();
  }
   
   for(unsigned int ideep = 0; ideep < idmax; ideep++)
   {
     unsigned int cur_node_id(0);
     unsigned int pre_node_id(0);
     for(unsigned int i = 0; i < cur_reac.channels.size(); i++)
     {
         if(cur_reac.br_path[i].size() <= ideep)
         {
           if(cur_node_id != 0)pre_node_id++;
           cur_node_id = 0;
           continue;
         }
        unsigned int tmp = cur_reac.br_path[i][ideep];    
         if(tmp != cur_node_id)
         {
           if(cur_node_id != 0)
           {
              pre_node_id++;
           }
           if(tmp != pre_node_id + 1)
           {
               std::cerr << "duplicate node\n\t" << cur_reac.channels[i]
                         << "\nin reaction\n\t" << cur_reac.name << std::endl;
               antioch_error();
           }
           cur_node_id = tmp;
        }
     }
   }
   return;
}


template <typename Scalar>
void treat_reaction(LocalReaction &reac, Antioch::ReactionSet<Scalar> &reaction_set, unsigned int n_species)
{
   Antioch::ReactionType::ReactionType typeReaction(Antioch::ReactionType::ELEMENTARY);
   Antioch::KineticsModel::KineticsModel kineticsModel(Antioch::KineticsModel::CONSTANT);

   Planet::KineticsBranchingStructure<Scalar> prob_reac;
   prob_reac.set_n_channels(reac.channels.size());

   std::vector<std::string> out;
   Antioch::SplitString(reac.channels.front(),";",out,false);
   shave_string(out);


// k
// create necessary pdf objects
   std::string line = reac.channels.front();
   std::vector<Planet::PDFName::PDFName> k_pdf_type;
   std::string loc_pdf = out[2].substr(0,out[2].find('('));
   shave_string(loc_pdf);
   k_pdf_type.push_back(prob_reac.pdf_map().at(loc_pdf));
   if(line.find("beta") != std::string::npos) //DR then
   {
     if(out[1].find("RD") == std::string::npos)antioch_error();
     kineticsModel = Antioch::KineticsModel::HERCOURT_ESSEN; //type reaction
     std::string str = out.back().substr(0,out.back().find('('));
     shave_string(str);
     str = str.substr(str.find(" ") + 1);
     shave_string(str);
     k_pdf_type.push_back(prob_reac.pdf_map().at(str));

   }
   prob_reac.set_k_pdf(k_pdf_type);


/// fill 'em
   std::string str_param = out[2].substr(out[2].find('(')+1, out[2].find(')') - out[2].find('(')-1);
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
     k_pdf_params_str.clear();
     str_param = out.back().substr(out.back().find('(')+1, out.back().find(')') - out.back().find('(')-1);
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

//// dirichlet

   std::vector<BrPath> br_path;
   br_path.resize(reac.channels.size());
   out.clear();
   Antioch::SplitString(reac.channels[0],";",out,false);
   shave_string(out);
   std::vector<std::string> node_names;


/// first, fill prob_reac
// first br level, in out[3]
   Planet::PDFName::PDFName pdf_type;
   pdf_type = prob_reac.pdf_map().at(out[3].substr(0,out[3].find('(')));

   std::vector<std::vector<Scalar> >  master_data;
   Planet::BranchingRatioNode<Scalar> *master_node = new Planet::BranchingRatioNode<Scalar>();
   master_node->set_id("master");
   master_node->set_pdf(pdf_type);

   unsigned int master_chan(0);
   unsigned int current_node(0);

   for(unsigned int ibr = 0; ibr < reac.channels.size(); ibr++)
   {
//scan conditions to skip
      if(!reac.br_path[ibr].empty()) 
      {
        if(current_node == reac.br_path[ibr].front()) //if already done this one
        {
           //we just add master in his path
           br_path[ibr].name.push_back("master");
           br_path[ibr].id_channel.push_back(master_chan - 1); //index, not size !!!
           continue;
        }else
        {
          current_node = reac.br_path[ibr].front();
        }
      }

      out.clear();
      Antioch::SplitString(reac.channels[ibr],";",out,false);
      shave_string(out);
      std::string data = out[3].substr(out[3].find('(') + 1, out[3].find(')') - out[3].find('(') - 1);
      std::vector<std::string> datas;
      master_data.push_back(std::vector<Scalar>());
      if(!Antioch::SplitString(data,",",datas,false))datas.push_back(data);
      master_data.back().resize(datas.size(),0.L);
      for(unsigned int id = 0; id < datas.size(); id++)
      {
        master_data.back()[id] = std::atof(datas[id].c_str());
      }
      br_path[ibr].name.push_back("master");
      br_path[ibr].id_channel.push_back(master_chan);
      master_chan++;
   }
   master_node->set_n_channels(master_chan);
   std::vector<Scalar> data;
   for(unsigned int n = 0; n < master_data.front().size(); n++)
   {
     for(unsigned int i = 0; i < master_data.size(); i++)
     {
        data.push_back(master_data[i][n]);
     }
   }
   if(pdf_type == Planet::PDFName::DiUn)
   {
       data[0] = data.size();
       data.resize(1);
   }
   master_node->set_pdf_parameters(data);
   prob_reac.add_node(master_node);


   unsigned int nlevel_max(0); //size
   for(unsigned int ibr = 0; ibr < reac.channels.size(); ibr++)
   {
       if(reac.br_path[ibr].size() > nlevel_max)nlevel_max = reac.br_path[ibr].size();
   }
   for(unsigned int ilevel = 0; ilevel < nlevel_max; ilevel++)//scanning levels
   {
     unsigned int max_node(0);


// finding the number of nodes at that deepness
// nodes are simply counter, 1,2,...,3
     for(unsigned int ibr = 0; ibr < reac.channels.size(); ibr++)
     {
        if(reac.br_path[ibr].size()  < ilevel + 1)continue;
        if(reac.br_path[ibr][ilevel] > max_node)max_node = reac.br_path[ibr][ilevel]; //this is a size, not an index
     }

     for(unsigned int inode = 1; inode <= max_node; inode++) //scanning the inner br nodes by name
     {
       unsigned int child_node(0);
       Planet::BranchingRatioNode<Scalar> * node = new Planet::BranchingRatioNode<Scalar>();
       std::string node_id;
       std::string node_pdf;
       std::vector<Scalar> node_pdf_params;

       unsigned int nchan(0);
       std::vector<std::vector<Scalar> > data_chan;

       for(unsigned int ibr = 0; ibr < reac.channels.size(); ibr++)
       {
         if(reac.br_path[ibr].size() < ilevel + 1)continue; //if current deepness is not enough
         if(reac.br_path[ibr][ilevel] != inode)continue;    //if that deepness but not that node

        if(reac.br_path[ibr].size() > ilevel + 1) //if deeper, check we're not scanning a node
        {
          if(child_node == reac.br_path[ibr][ilevel + 1]) //if already done this one
          {
           //we just add the name in his path
           br_path[ibr].name.push_back(node_id);
           br_path[ibr].id_channel.push_back(nchan - 1); //index, not size !!!
           continue;
          }
          child_node = reac.br_path[ibr][ilevel + 1];
        }

         // name is simply the point-separated lineage i.e. master."node at level 0"."node at level 1".[...]."current node" 
         if(node_id.empty())
         {
           std::stringstream oss;
           oss << "master";
           for(unsigned int il = 0; il <= ilevel;il++)oss << "." << reac.br_path[ibr][il];
           node_id = oss.str();
         }
         data_chan.push_back(std::vector<Scalar>());

         std::vector<std::string> tmp;
         Antioch::SplitString(reac.channels[ibr],";",tmp,false);//
         if(tmp.size() < 4)antioch_error();
//pdf
         shave_string(tmp);
         if(node_pdf.empty())
         {
           node_pdf = tmp[4 + ilevel].substr(0,tmp[4 + ilevel].find('(')); //level 0 starts at 4
           shave_string(node_pdf);
         }

//data
         std::string data = tmp[4 + ilevel].substr(tmp[4 + ilevel].find('(') + 1, 
                                                   tmp[4 + ilevel].find(')') - tmp[4 + ilevel].find('(') - 1);
         std::vector<std::string> datas_chan;
         if(!Antioch::SplitString(data,",",datas_chan,false))datas_chan.push_back(data);
         data_chan.back().resize(datas_chan.size(),0.L);
         for(unsigned int id = 0; id < datas_chan.size(); id++)
         {
            data_chan.back()[id] = std::atof(datas_chan[id].c_str());
         }

         br_path[ibr].name.push_back(node_id);     //name of node
         br_path[ibr].id_channel.push_back(nchan); //channel in node

         nchan++;

       }/// end ibr loop

//reformatting the data in the right order
       for(unsigned int n = 0; n < data_chan.front().size(); n++)
       {
         for(unsigned int i = 0; i < data_chan.size(); i++)
         {
           node_pdf_params.push_back(data_chan[i][n]);
         }
       }
       if(prob_reac.pdf_map().at(node_pdf) == Planet::PDFName::DiUn)
       {
         node_pdf_params[0] = node_pdf_params.size();
         node_pdf_params.resize(1);
       }

       node->set_id(node_id);
       node->set_n_channels(nchan);
       node->set_pdf(prob_reac.pdf_map().at(node_pdf));
       node->set_pdf_parameters(node_pdf_params);
       prob_reac.add_node(node);
     } //end node scan
   } // end level scan

//filling the reaction

   for(unsigned int ibr = 0; ibr < reac.channels.size(); ibr++)
   {
     for(unsigned int i = 0; i < br_path[ibr].name.size(); i++)
     {
       prob_reac.add_to_br_path(ibr,br_path[ibr].name[i],br_path[ibr].id_channel[i]);
     }
   }

// then, create the reaction
   for(unsigned int ibr = 0; ibr < reac.channels.size(); ibr++)
   {
       std::string equation = reac.channels[ibr].substr(0,reac.channels[ibr].find(";"));
       shave_string(equation);
//reaction
       Antioch::Reaction<Scalar>* my_rxn = Antioch::build_reaction<Scalar>(n_species, equation, false, typeReaction, kineticsModel);

       Antioch::KineticsType<Scalar>* rate = prob_reac.create_rate_constant(ibr,kineticsModel);
       my_rxn->add_forward_rate(rate);

       std::vector<std::string> reactants;
       size_t loc(0);
       while(loc < equation.find("->"))
       {
          size_t locnext = equation.find(" ",loc);
          if(locnext == loc)
          {
             loc++;
             continue;
          }
          reactants.push_back(equation.substr(loc,locnext-loc));
          loc = locnext;
       }
       shave_string(reactants);

       std::vector<unsigned int> stoi;
       stoi.resize(reactants.size(),1);
       condense_molecule(reactants,stoi);

       for(unsigned int ir = 0; ir < reactants.size(); ir++)
       {
          my_rxn->add_reactant(reactants[ir],reaction_set.chemical_mixture().active_species_name_map().find( reactants[ir] )->second,stoi[ir]);
       }

       reactants.clear();
       loc = equation.find("->") + 2;
       while(loc < equation.size())
       {
          size_t locnext = equation.find(" ",loc);
          if(locnext == loc)
          {
             loc++;
             continue;
          }
          reactants.push_back(equation.substr(loc,locnext-loc));
          loc = locnext;
       }
       shave_string(reactants);
       stoi.clear();
       stoi.resize(reactants.size(),1);
       condense_molecule(reactants,stoi);

       for(unsigned int ip = 0; ip < reactants.size(); ip++)
       {
          my_rxn->add_product(reactants[ip],reaction_set.chemical_mixture().active_species_name_map().find( reactants[ip])->second,stoi[ip]);
       }
      
       reaction_set.add_reaction(my_rxn);     
   }
}

template <typename Scalar>
void read_reactions(const std::string &file_reac, Antioch::ReactionSet<Scalar> &reaction_set, unsigned int n_species)
{
  std::ifstream data(file_reac.c_str());
  std::string line;
  LocalReaction cur_reac;
  while(!data.eof())
  {
     if(!getline(data,line))break;
     if(line[0] == '#')continue;
     std::vector<std::string> out;
     Antioch::SplitString(line,";",out,false);
     if(out.empty())antioch_error();
     std::string name = out[0].substr(0,out[0].find(' '));
     line.erase(0,name.size());
     std::vector<std::string> branch;
     Antioch::SplitString(name,".",branch,false);
     if(branch.empty())antioch_error();

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

        sanity_check_reaction(cur_reac);


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
  read_reactions(file_reac,reaction_set, mixture.n_species());
  Antioch::KineticsEvaluator<Scalar> reactions_system(reaction_set,0);
  std::vector<Scalar> molar_concentrations, molar_sources;

// now the concentrations
  prepare_the_ionosphere(file_neu_conc,molar_concentrations,mixture);

//solve
  Scalar T(200.);
  Planet::AtmosphericSteadyState solver;
  solver.steady_state(reactions_system, ss_species, mixture, T, molar_concentrations, molar_sources);

  int return_flag(1);

  for(unsigned int s = 0; s < ss_species.size(); s++)
  {
      std::cout << ss_species[s] << " " << molar_sources[s] << std::endl;
  }


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
