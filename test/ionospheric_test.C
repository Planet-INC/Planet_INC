//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

//Antioch
#include "antioch/chemical_mixture.h"
#include "antioch/kinetics_evaluator.h"

//Planet
#include "planet/atmospheric_steady_state.h"

//C++
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>

void read_the_species(const std::string &file_spec, std::vector<std::string> & all_species)
{
  all_species.clear();
  std::ifstream data(file_spec.c_str());
  while(!data.eof())
  {
      std::string spec;
      data >> spec;
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
void read_reactions(const std::string &file_reac, Antioch::ReactionSet<Scalar> &reaction_set)
{
  std::ifstream data(file_reac.c_str());
  data.close();
}

template <typename Scalar>
void prepare_the_ionosphere(const std::string &file_neu_conc, std::vector<Scalar> &concentrations,
                            const Antioch::ChemicalMixture<Scalar> &mixture)
{
  concentrations.resize(mixture.n_species(),-1.L);
  std::ifstream data(file_neu.c_str());
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
  std::vector<Scalar> concentrations, sources;

// now the concentrations
  prepare_the_ionosphere(file_neu_conc,concentrations,mixture);


//solve
  Scalar T(200.);
  Planet::AtmosphericSteadyState solver;
  solver.steady_state(reactions_system, ss_species, mixture, T, molar_concentrations, molar_sources);

  int return_flag(0);

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
