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

#ifndef PLANET_PLANET_PHYSICS_HELPER_H
#define PLANET_PLANET_PHYSICS_HELPER_H

//Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/metaprogramming.h"
#include "antioch/vector_utils.h"
#include "antioch/physical_constants.h"
#include "antioch/string_utils.h"
#include "antioch/kinetics_parsing.h"
#include "antioch/reaction_parsing.h"

//Planet
#include "planet/diffusion_evaluator.h"
#include "planet/atmospheric_kinetics.h"
#include "planet/kinetics_branching_structure.h"
#include "planet/branching_ratio_node.h"

// libMesh
#include "libmesh/libmesh_common.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"

namespace Planet
{

  //! small structure to help parsing reactions
  struct BrPath 
  {
    std::vector<std::string> name;
    std::vector<unsigned int > id_channel;
  };

  //! small structure to help parsing reactions
  struct LocalReaction
  {
    std::string name;
    std::vector<std::string> channels;
    std::vector<std::vector<unsigned int> > br_path;
  };

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  class PlanetPhysicsHelper
  {
  public:

    PlanetPhysicsHelper( const GetPot& input );

    ~PlanetPhysicsHelper();

    //!fills molar_concentrations_first_guess with barometric equation
    template<typename StateType, typename VectorStateType>
    void first_guess(VectorStateType & molar_concentrations_first_guess, const StateType z) const;

    template<typename StateType>
    StateType first_guess(unsigned int s, const StateType z) const;

    //!fills lower boundary conditions
    template<typename VectorStateType>
    void lower_boundary_dirichlet(VectorStateType & lower_boundary) const;

    //!fills upper boundary conditions
    template<typename VectorStateType>
    void upper_boundary_neumann(VectorStateType & upper_boundary, const VectorStateType &molar_densities) const;

    //!return upper boundary derivatives of s with respect to i
    CoeffType dupper_boundary_neumann_s_dn_i(unsigned int s, unsigned int i) const;

    CoeffType upper_boundary_neumann(const VectorCoeffType &molar_densities, unsigned int s) const;

    const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>& composition() const;

    const Antioch::ReactionSet<CoeffType>& neutral_reaction_set() const;

    const Antioch::ReactionSet<CoeffType>& ionic_reaction_set() const;

    const PhotonOpacity<CoeffType,VectorCoeffType>& tau() const;

    const std::vector<std::vector<BinaryDiffusion<CoeffType> > >& bin_diff_coeff() const;

    const AtmosphericTemperature<CoeffType,VectorCoeffType>& temperature() const;

    const VectorCoeffType& lambda_hv() const;

    const VectorCoeffType& phy1AU() const;

    const Antioch::ParticleFlux<VectorCoeffType> & phy_at_top() const;

    const std::vector<std::string>& medium() const;

    CoeffType K0() const;

    CoeffType scaling_factor() const;

    const std::vector<Antioch::Species> ss_species() const;

  private:

    AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>*  _composition; //for first guess

    // Additional data structures that need to be cached
// temperatures
    AtmosphericTemperature<CoeffType,VectorCoeffType>* _temperature;

// mixtures
    Antioch::ChemicalMixture<CoeffType>* _neutral_species;
    Antioch::ChemicalMixture<CoeffType>* _ionic_species;

// chemistry
    Antioch::ReactionSet<CoeffType>* _neutral_reaction_set;
    Antioch::ReactionSet<CoeffType>* _ionic_reaction_set;

// photons related
    Chapman<CoeffType>* _chapman;

    Antioch::ParticleFlux<VectorCoeffType> _phy1AU;
    Antioch::ParticleFlux<VectorCoeffType> _phy_at_top;

    std::vector<CrossSection<VectorCoeffType> > _hv_cross_section;

    PhotonOpacity<CoeffType,VectorCoeffType>* _tau;

//diffusion
//  molecular

    std::vector<std::vector<BinaryDiffusion<CoeffType> > > _bin_diff_coeff;

    std::vector<std::string>      _medium;
    std::vector<Antioch::Species> _ss_species;

//  eddy
    CoeffType _K0;


    CoeffType _scaling_factor;

    /*! Convenience method to hide all the construction code for
        composition, kinetics, and diffusion */
    void build( const GetPot& input );

    /*! Convenience method within a convenience method */
    void build_temperature( const GetPot& input );

    /*! Convenience method within a convenience method */
    void build_species( const GetPot& input, std::vector<std::string>& neutrals, std::vector<std::string>& ionic_species );

    /*! Convenience method within a convenience method */
    void build_opacity( const GetPot& input );

    /*! Convenience method within a convenience method */
    void build_reaction_sets( const GetPot& input );

    /*! Convenience method within a convenience method */
    void build_composition( const GetPot& input, VectorCoeffType& molar_frac,
                            CoeffType dens_tot, VectorCoeffType& tc, VectorCoeffType& hard_sphere_radius);

    void build_diffusion( std::vector<std::vector<std::vector<CoeffType> > >& bin_diff_data,
                          std::vector<std::vector<DiffusionType> >& bin_diff_model,
                          const std::vector<std::string>& neutrals);

    // Helper functions for parsing data
    void read_temperature(VectorCoeffType& T0, VectorCoeffType& Tz, const std::string& file) const;

    void fill_neutral_reactions_elementary(const std::string &neutral_reactions_file,
                                           Antioch::ReactionSet<CoeffType>& neutral_reaction_set ) const;

    void fill_neutral_reactions_falloff(const std::string &neutral_reactions_file,
                                        Antioch::ReactionSet<CoeffType>& neutral_reaction_set) const;

    void fill_ionic_system_reaction(const std::string &file_ions_reac,
                                    Antioch::ReactionSet<CoeffType> &ionic_reaction_set) const;

    void sanity_check_reaction(LocalReaction &cur_reac) const;

    void treat_reaction(LocalReaction &cur_reac, Antioch::ReactionSet<CoeffType> &reaction_set) const;

    void read_flyby_infos(const std::vector<std::string>& neutrals,
                          CoeffType& dens_tot, std::vector<CoeffType>& molar_frac,
                          CoeffType& chi, CoeffType& K0,
                          const std::string& file_flyby,
                          const std::string& root_input) const;

    void read_cross_section( const std::string &file,
                            VectorCoeffType &lambda, VectorCoeffType &sigma ) const;

    void read_hv_flux(Antioch::ParticleFlux<VectorCoeffType> & phy1AU, 
                      Antioch::ParticleFlux<VectorCoeffType> & phy_at_top, 
                      const std::string &file) const;

    void read_neutral_characteristics(const std::vector<std::string>& neutrals,
                                      VectorCoeffType& tc,VectorCoeffType& hard_sphere_radius,
                                      std::vector<std::vector<std::vector<CoeffType> > >& bin_diff_data,
                                      std::vector<std::vector<DiffusionType> >& bin_diff_model,
                                      const std::string & file_neutral_charac) const;

    void fill_molar_frac(const std::vector<std::string> &neutrals, VectorCoeffType& molar_frac,
                         const std::string &file, const std::string &root_input) const;

    void parse_equation(std::vector<std::string> &reactants, std::vector<std::string> &products,
                        std::string &line, bool &skip, const Antioch::ChemicalMixture<CoeffType>& chem_mixture,
                        std::string &equation, std::vector<unsigned int> &stoi_reac,
                        std::vector<unsigned int> &stoi_prod) const;

    void condense_molecule(std::vector<unsigned int> &stoi, std::vector<std::string> &mol) const;

    void read_photochemistry_reac(const std::string &hv_file, const std::string &reac,
                                  Antioch::ReactionSet<CoeffType> &neutral_reaction_set) const;

    void shave_string(std::string &str) const;

    void shave_strings(std::vector<std::string> &stock) const;

    void sanity_check_chemical_system() const;

    bool check_chemical_balance(const Antioch::ReactionSet<CoeffType> &reaction_set,const std::vector<Antioch::Species> &species, const std::string & help) const;

  };

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::PlanetPhysicsHelper( const GetPot& input )
    : _composition(NULL),
      _temperature(NULL),
      _neutral_species(NULL),
      _ionic_species(NULL),
      _neutral_reaction_set(NULL),
      _ionic_reaction_set(NULL),
      _chapman(NULL),
      _tau(NULL),
      _scaling_factor(1e13) // sensible default
  {
    this->build(input);
    this->sanity_check_chemical_system();

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::~PlanetPhysicsHelper()
  {
    delete _tau;
    delete _chapman;
    delete _ionic_reaction_set;
    delete _neutral_reaction_set;
    delete _ionic_species;
    delete _neutral_species;
    delete _temperature;
    delete _composition;

    return;
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::first_guess(VectorStateType & molar_concentrations_first_guess, const StateType z) const
  {
      _composition->first_guess_densities(z,molar_concentrations_first_guess);
      for(unsigned int i = 0; i < molar_concentrations_first_guess.size(); i++)
      {
         molar_concentrations_first_guess[i] /= _scaling_factor;
      }
      return;
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType>
  StateType PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::first_guess(unsigned int s, const StateType z) const
  {
    return _composition->first_guess_density(z,s) / _scaling_factor;
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename VectorStateType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::lower_boundary_dirichlet(VectorStateType & lower_boundary) const
  {
      _composition->lower_boundary_concentrations(lower_boundary);
      for(unsigned int i = 0; i < lower_boundary.size(); i++)
      {
         lower_boundary[i] /= _scaling_factor;
      }
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename VectorStateType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::upper_boundary_neumann(VectorStateType & upper_boundary, const VectorStateType &molar_densities) const
  {
      _composition->upper_boundary_fluxes(upper_boundary, molar_densities);
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  CoeffType PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::upper_boundary_neumann(const VectorCoeffType &molar_densities, unsigned int s) const
  {
    return _composition->upper_boundary_flux(molar_densities,s);
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::build(const GetPot& input)
  {
    // Parse medium
    if( !input.have_variable("Planet/medium") )
      {
        std::cerr << "Error: Could not find medium!" << std::endl;
        antioch_error();
      }
    unsigned int n_medium = input.vector_variable_size("Planet/medium");
    _medium.resize(n_medium);
    for( unsigned int s = 0; s < n_medium; s++ )
      {
        _medium[s] = input("Planet/medium", "DIE!", s);
      }

    // Parse neutrals, ions
    std::vector<std::string> neutrals;
    std::vector<std::string> ions;

    this->build_species(input, neutrals, ions);

    // Parse stuff
    if( !input.have_variable("Planet/file_flyby") )
      {
        std::cerr << "Error: Could not find file_flyby filename!" << std::endl;
        antioch_error();
      }
    if( !input.have_variable("Planet/root_input") )
      {
        std::cerr << "Error: Could not find root_input filename!" << std::endl;
        antioch_error();
      }

    std::string file_flyby = input("Planet/file_flyby", "DIE!");
    std::string root_input = input("Planet/root_input", "DIE!");

    VectorCoeffType molar_frac;
    CoeffType dens_tot;
    CoeffType chi;

    this->read_flyby_infos( neutrals, dens_tot, molar_frac, chi, _K0,
                            file_flyby, root_input );

    // Parse more stuff
    if( !input.have_variable("Planet/file_neutral_charac") )
      {
        std::cerr << "Error: Could not find file_neutral_charac filename!" << std::endl;
        antioch_error();
      }

    std::string file_neutral_charac = input("Planet/file_neutral_charac", "DIE!");

    //binary diffusion
    std::vector<std::vector<std::vector<CoeffType> > > bin_diff_data;    // (medium,species)[]
    std::vector<std::vector<DiffusionType> > bin_diff_model; // (medium,species)

    bin_diff_data.resize(_medium.size());
    bin_diff_model.resize(_medium.size());
    _bin_diff_coeff.resize(_medium.size());

    VectorCoeffType tc;
    VectorCoeffType hard_sphere_radius;

    this->read_neutral_characteristics( neutrals, tc, hard_sphere_radius, bin_diff_data, bin_diff_model,
                                        file_neutral_charac);

    _chapman = new Chapman<CoeffType>(chi);

    this->build_diffusion(bin_diff_data, bin_diff_model, neutrals);

    // Must be called after: build_species, chapman
    this->build_opacity(input);

    this->build_temperature(input);

    // Must be called after: build_species
    this->build_reaction_sets(input);

    // Must be called after: build_temperature, build_species
    this->build_composition(input, molar_frac, dens_tot, tc, hard_sphere_radius);

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::build_temperature( const GetPot& input )
  {
    if( !input.have_variable("Planet/temperature_file") )
      {
        std::cerr << "Error: temperature_file not found in input file!" << std::endl;
        antioch_error();
      }

    std::string input_T = input( "Planet/temperature_file", "DIE!" );
    std::vector<CoeffType> T0,Tz;
    this->read_temperature(T0,Tz,input_T);
    _temperature = new AtmosphericTemperature<CoeffType,VectorCoeffType>(T0,T0,Tz,Tz);

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::build_species( const GetPot& input, std::vector<std::string>& neutrals, std::vector<std::string>& ions )
  {
    if( !input.have_variable("Planet/neutral_species") )
      {
        std::cerr << "Error: neutral_species not found in input file!" << std::endl;
        antioch_error();
      }

    unsigned int n_ionic(0);
    if( input.have_variable("Planet/ionic_species") )
      {
        n_ionic = input.vector_variable_size("Planet/ionic_species");
      }

    // Read neutral and ionic species from input
    unsigned int n_neutral = input.vector_variable_size("Planet/neutral_species");

    neutrals.resize(n_neutral);       //neutral system
    ions.resize(n_neutral + n_ionic); //ionic system
    _ss_species.resize(n_ionic);      //sub system to steady state

    for( unsigned int s = 0; s < n_neutral; s++ )
      {
        neutrals[s] = input("Planet/neutral_species", "DIE!", s);
        ions[s]     = input("Planet/neutral_species", "DIE!", s);
      }

    for( unsigned int s = 0; s < n_ionic; s++ )
      {
        ions[n_neutral + s] = input("Planet/ionic_species", "DIE!", s);
      }

    _neutral_species = new Antioch::ChemicalMixture<CoeffType>(neutrals);
    _ionic_species   = new Antioch::ChemicalMixture<CoeffType>(ions);

    for( unsigned int s = 0; s < n_ionic; s++ )
      {
        _ss_species[s] = _ionic_species->species_name_map().at(ions[s + n_neutral]);
      }

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::build_opacity( const GetPot& input )
  {
    _tau = new PhotonOpacity<CoeffType,VectorCoeffType>(*_chapman);


//hv
    if( !input.have_variable("Planet/input_hv") )
      {
        std::cerr << "Error: input_hv not found in input file!" << std::endl;
        antioch_error();
      }

    std::string input_hv  = input("Planet/input_hv", "DIE!" );

    this->read_hv_flux(_phy1AU, _phy_at_top, input_hv);//photon.angstrom-1.s-1

//cross-section

    if( !input.have_variable("Planet/input_cross_section_root") )
      {
        std::cerr << "Error: input_cross_section_root not found in input file!" << std::endl;
        antioch_error();
      }

    if( !input.have_variable("Planet/absorbing_species") )
      {
        std::cerr << "Error: absorbing_species not found in input file!" << std::endl;
        antioch_error();
      }


    unsigned int n_absorbing = input.vector_variable_size("Planet/absorbing_species");

    std::vector<std::string> abs_file(n_absorbing);

    for( unsigned int s = 0; s < n_absorbing; s++ )
      {
        abs_file[s] = std::string(input("Planet/input_cross_section_root","DIE!")) + std::string(input("Planet/absorbing_species", "DIE!", s));
      }


   for(unsigned int s = 0; s < n_absorbing; s++)
   {

      std::string species = input("Planet/absorbing_species", "DIE!", s);

      if(!_neutral_species->active_species_name_map().count(species))
      {
         std::cerr << "Unknown species \"" << species << "\".  Forgot to add it to the neutral_species entry?" << std::endl;
         antioch_error();
      }

      std::vector<CoeffType> lambda, sigma;

      this->read_cross_section(abs_file[s], lambda, sigma); //cm2.angstrom-1

      _tau->add_cross_section( lambda, sigma, _neutral_species->species_name_map().at(species),          // Antioch::Species
                                              _neutral_species->active_species_name_map().at(species) ); // id

    }

    _tau->update_cross_section(_phy1AU.abscissa());

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::build_reaction_sets( const GetPot& input )
  {
    // Build up reaction sets
    _neutral_reaction_set = new Antioch::ReactionSet<CoeffType>(*_neutral_species);
    _ionic_reaction_set = new Antioch::ReactionSet<CoeffType>(*_ionic_species);

    if( !input.have_variable("Planet/input_reactions_elem") )
      {
        std::cerr << "Error: could not find input_reactions_elem filename!" << std::endl;
        antioch_error();
      }

    if( !input.have_variable("Planet/input_reactions_fall") )
      {
        std::cerr << "Error: could not find input_reactions_fall filename!" << std::endl;
        antioch_error();
      }

    if( !input.have_variable("Planet/input_photoreactions_root") )
      {
        std::cerr << "Error: could not find input_photoreactions_root path name!" << std::endl;
        antioch_error();
      }

    if( !input.have_variable("Planet/photo_reacting_species") )
      {
        std::cerr << "Error: photo_reacting_species not found!" << std::endl;
        antioch_error();
      }

    std::string input_reactions_elem = input( "Planet/input_reactions_elem", "DIE!" );
    std::string input_reactions_fall = input( "Planet/input_reactions_fall", "DIE!" );

    // here read the reactions, Antioch will take care of it once hdf5, no ionic reactions there
    //Kooij / Arrhenius + photochem
    this->fill_neutral_reactions_elementary(input_reactions_elem, *_neutral_reaction_set);

    this->fill_neutral_reactions_falloff(input_reactions_fall, *_neutral_reaction_set);

    //now the photochemical ones
    unsigned int n_hv_reacting = input.vector_variable_size("Planet/photo_reacting_species");

    std::vector<std::string> hv_file(n_hv_reacting);

    for( unsigned int s = 0; s < n_hv_reacting; s++ )
      {
        hv_file[s] = std::string(input("Planet/input_photoreactions_root","DIE!")) + std::string(input("Planet/photo_reacting_species", "DIE!", s));
      }


   for(unsigned int s = 0; s < n_hv_reacting; s++)
   {

      std::string species = input("Planet/photo_reacting_species", "DIE!", s);

      if(!_neutral_species->active_species_name_map().count(species))
      {
         std::cerr << "Unknown species \"" << species << "\".  Forgot to add it to the neutral_species entry?" << std::endl;
         antioch_error();
      }

      this->read_photochemistry_reac(hv_file[s], species, *_neutral_reaction_set);

    }


    if( input.have_variable("Planet/input_ions_reactions") )
      {
        std::string file_ions_reac = input( "Planet/input_ions_reactions", "DIE!" );
        this->fill_ionic_system_reaction(file_ions_reac,*_ionic_reaction_set);
      }

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::build_composition( const GetPot& input, VectorCoeffType& molar_frac, CoeffType dens_tot, VectorCoeffType& tc, VectorCoeffType& hard_sphere_radius)
  {

    // Build AtmosphericMixture
    _composition = new AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>( *_neutral_species, *_ionic_species, *_temperature );

    if( !input.have_variable("Planet/zmin") )
      {
        std::cerr << "Error: zmin not found in input file!" << std::endl;
        antioch_error();
      }

    if( !input.have_variable("Planet/zmax") )
      {
        std::cerr << "Error: zmax not found in input file!" << std::endl;
        antioch_error();
      }

    CoeffType zmin = input("Planet/zmin", 0.0 );
    CoeffType zmax = input("Planet/zmax", 0.0 );

    _composition->init_composition(molar_frac, dens_tot, zmin, zmax);
    _composition->set_thermal_coefficient(tc);
    _composition->set_hard_sphere_radius(hard_sphere_radius);

    _scaling_factor = dens_tot;

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::build_diffusion( std::vector<std::vector<std::vector<CoeffType> > >& bin_diff_data,
                                                                                        std::vector<std::vector<DiffusionType> >& bin_diff_model,
                                                                                        const std::vector<std::string>& neutrals)
  {
    std::vector<Antioch::Species> spec;

    for(unsigned int s = 0; s < neutrals.size(); s++)
      {
        spec.push_back(_neutral_species->species_name_map().at(neutrals[s]));
      }

    for(unsigned int n = 0; n < _medium.size(); n++)
      {
        _bin_diff_coeff[n].resize(neutrals.size());

        for(unsigned int s = 0; s < neutrals.size(); s++)
          {
            _bin_diff_coeff[n][s] = BinaryDiffusion<CoeffType>( spec[n], spec[s], bin_diff_data[n][s][0], bin_diff_data[n][s][1], bin_diff_model[n][s]);
          }
      }

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::read_temperature(VectorCoeffType& T0, VectorCoeffType& Tz, const std::string& file) const
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
        CoeffType t(0.),tz(0.),dt(0.),dtz(0.);
        temp >> t >> tz >> dt >> dtz;
        if(tz < 1.)continue; //altitude < 1. possible only if nothing is read
        T0.push_back(t);
        Tz.push_back(tz);
      }
    temp.close();

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::fill_neutral_reactions_elementary(const std::string &neutral_reactions_file,
                                                                                                         Antioch::ReactionSet<CoeffType> &neutral_reaction_set) const
  {
    //here only simple ones: bimol Kooij/Arrhenius model
    std::ifstream data(neutral_reactions_file.c_str());
    if( !data )
      {
        std::cerr << "Could not open file " << neutral_reactions_file << std::endl;
        antioch_error();
      }
    std::string line;
    Antioch::KineticsModel::KineticsModel kineticsModel(Antioch::KineticsModel::KOOIJ);
    Antioch::ReactionType::ReactionType reactionType(Antioch::ReactionType::ELEMENTARY);
    const Antioch::ChemicalMixture<CoeffType>& chem_mixture = neutral_reaction_set.chemical_mixture();
    while(!data.eof())
      {
        if(!getline(data,line))break;
        shave_string(line);
        if(line[0] == '#' || line.empty())continue;
        std::vector<std::string> reactants;
        std::vector<std::string> products;

        bool skip(false);
        std::string equation;
        std::vector<unsigned int> stoi_reac;
        std::vector<unsigned int> stoi_prod;
        parse_equation(reactants,products,line,skip,chem_mixture,equation,stoi_reac,stoi_prod);

        kineticsModel = Antioch::KineticsModel::KOOIJ;

        if(skip)continue;

        VectorCoeffType dataf;
        std::vector<std::string> str_data;
        Antioch::SplitString(line," ",str_data,false);
        if(str_data.size() != 4)
          {
            std::cerr << "data are badly shaped, need 4 numbers in this line\n"
                      << line << std::endl;
            antioch_error();
          }

        dataf.push_back(std::atof(str_data[0].c_str())); //Cf
        CoeffType temp = std::atof(str_data[1].c_str());
        CoeffType temp2 = std::atof(str_data[2].c_str());
        if(temp == 0 && temp2 == 0) //constant
          {
            kineticsModel = Antioch::KineticsModel::CONSTANT;
          }else if(temp == 0)
          {
            kineticsModel = Antioch::KineticsModel::ARRHENIUS;
            dataf.push_back(std::atof(str_data[2].c_str()));//Ea
            dataf.push_back(Antioch::Constants::R_universal<CoeffType>()*1e-3); //scale (R in J/mol/K)
          }else
          {
            dataf.push_back(std::atof(str_data[1].c_str())); //beta
            dataf.push_back(std::atof(str_data[2].c_str()));//Ea
            dataf.push_back(1.); //Tref
            dataf.push_back(Antioch::Constants::R_universal<CoeffType>()*1e-3); //scale (R in J/mol/K)
          }


        Antioch::KineticsType<CoeffType, VectorCoeffType>* rate = Antioch::build_rate<CoeffType,VectorCoeffType>(dataf,kineticsModel); //kinetics rate
        Antioch::Reaction<CoeffType> * reaction = Antioch::build_reaction<CoeffType>(chem_mixture.n_species(), equation, false, reactionType, kineticsModel);
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
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::fill_neutral_reactions_falloff(const std::string &neutral_reactions_file,
                                                                                                      Antioch::ReactionSet<CoeffType> &neutral_reaction_set) const
  {
    //Lindemann
    std::ifstream data(neutral_reactions_file.c_str());
    if( !data )
      {
        std::cerr << "Could not open file " << neutral_reactions_file << std::endl;
        antioch_error();
      }
    std::string line;
    Antioch::KineticsModel::KineticsModel kineticsModel(Antioch::KineticsModel::KOOIJ);
    Antioch::ReactionType::ReactionType reactionType(Antioch::ReactionType::LINDEMANN_FALLOFF);
    const Antioch::ChemicalMixture<CoeffType>& chem_mixture = neutral_reaction_set.chemical_mixture();
    while(!data.eof())
      {
        if(!getline(data,line))break;
        shave_string(line);
        if(line[0] == '#' || line.empty())continue;
        std::vector<std::string> reactants;
        std::vector<std::string> products;

        bool skip(false);
        std::string equation;
        std::vector<unsigned int> stoi_reac;
        std::vector<unsigned int> stoi_prod;
        parse_equation(reactants,products,line,skip,chem_mixture,equation,stoi_reac,stoi_prod);

        kineticsModel = Antioch::KineticsModel::KOOIJ;

        if(skip)continue;
        VectorCoeffType dataf1,dataf2;
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
        CoeffType temp1 = std::atof(str_data[1].c_str());
        CoeffType temp2 = std::atof(str_data[4].c_str());
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
        dataf1.push_back(Antioch::Constants::R_universal<CoeffType>()*1e-3); //scale (R in J/mol/K)
        dataf2.push_back(Antioch::Constants::R_universal<CoeffType>()*1e-3); //scale (R in J/mol/K)

        Antioch::KineticsType<CoeffType, VectorCoeffType>* rate1 = Antioch::build_rate<CoeffType,VectorCoeffType>(dataf1,kineticsModel); //kinetics rate
        Antioch::KineticsType<CoeffType, VectorCoeffType>* rate2 = Antioch::build_rate<CoeffType,VectorCoeffType>(dataf2,kineticsModel); //kinetics rate
        Antioch::Reaction<CoeffType> * reaction = Antioch::build_reaction<CoeffType>(chem_mixture.n_species(), equation, false, reactionType, kineticsModel);
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

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::fill_ionic_system_reaction(const std::string &file_ions_reac,
                                                                                                  Antioch::ReactionSet<CoeffType> &ionic_reaction_set) const
  {
    std::ifstream data(file_ions_reac.c_str());
    std::string line;
    LocalReaction cur_reac;
    while(!data.eof())
    {
       if(!getline(data,line))break;
       shave_string(line);
       if(line[0] == '#' || line.empty())continue;
       std::vector<std::string> out;
       Antioch::SplitString(line,";",out,false);
       if(out.empty())
       {
           std::cerr << "Error: the line\n\t" << line << "\nis badly shaped" << std::endl;
           antioch_error();
       }
       std::string name = out[0].substr(0,out[0].find(' '));
       line.erase(0,name.size());
       std::vector<std::string> branch;
       Antioch::SplitString(name,".",branch,false);
       if(branch.empty())
       {
           std::cerr << "Error: the reaction name of line\n\t" << line << "\nis badly shaped: " << name << std::endl;
          antioch_error();
       }
  
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
  
          this->sanity_check_reaction(cur_reac);
  
          this->treat_reaction(cur_reac,ionic_reaction_set);
  
          cur_reac.br_path.clear();
          cur_reac.channels.clear();
  
          cur_reac.name = branch[1];
          cur_reac.channels.push_back(line);
          cur_reac.br_path.push_back(br);
       }
    }
// last one
    treat_reaction(cur_reac,ionic_reaction_set);

    data.close();

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::sanity_check_reaction(LocalReaction &cur_reac) const
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

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::treat_reaction(LocalReaction &reac, Antioch::ReactionSet<CoeffType> &reaction_set) const
  {
     Antioch::ReactionType::ReactionType typeReaction(Antioch::ReactionType::ELEMENTARY);
     Antioch::KineticsModel::KineticsModel kineticsModel(Antioch::KineticsModel::CONSTANT);
  
     Planet::KineticsBranchingStructure<CoeffType> prob_reac;
     prob_reac.set_n_channels(reac.channels.size());
  
     std::vector<std::string> out;
     Antioch::SplitString(reac.channels.front(),";",out,false);
     shave_strings(out);
  
  
  // k
  // create necessary pdf objects
     std::string line = reac.channels.front();
     std::vector<Planet::PDFName::PDFName> k_pdf_type;
     std::string loc_pdf = out[2].substr(0,out[2].find('('));
     shave_string(loc_pdf);
     k_pdf_type.push_back(prob_reac.pdf_map().at(loc_pdf));
     if(line.find("beta") != std::string::npos) //DR then
     {
       if(out[1].find("RD") == std::string::npos)
       {
           std::cerr << "Error: found beta parameter for not a DR reaction (keyword looked for is \"RD\") in reaction description\n\t" << line << std::endl;
           antioch_error();
       }
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
  
     std::vector<std::vector<CoeffType> > k_pdf_params;
     k_pdf_params.push_back(std::vector<CoeffType>());
     for(unsigned int i = 0; i < k_pdf_params_str.size(); i++)
     {
       std::stringstream oss(k_pdf_params_str[i]);
       CoeffType p;
       oss >> p;
       k_pdf_params[0].push_back(p);
     }
     if(line.find("beta") != std::string::npos)
     {
       k_pdf_params_str.clear();
       str_param = out.back().substr(out.back().find('(')+1, out.back().find(')') - out.back().find('(')-1);
       Antioch::SplitString(str_param,",",k_pdf_params_str,false);
       k_pdf_params.push_back(std::vector<CoeffType>());
       for(unsigned int i = 0; i < k_pdf_params_str.size(); i++)
       {
         std::stringstream oss(k_pdf_params_str[i]);
         CoeffType p;
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
     shave_strings(out);
     std::vector<std::string> node_names;
  
  
  /// first, fill prob_reac
  // first br level, in out[3]
     Planet::PDFName::PDFName pdf_type;
     pdf_type = prob_reac.pdf_map().at(out[3].substr(0,out[3].find('(')));
  
     std::vector<std::vector<CoeffType> >  master_data;
     Planet::BranchingRatioNode<CoeffType> *master_node = new Planet::BranchingRatioNode<CoeffType>();
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
        shave_strings(out);
        std::string data = out[3].substr(out[3].find('(') + 1, out[3].find(')') - out[3].find('(') - 1);
        std::vector<std::string> datas;
        master_data.push_back(std::vector<CoeffType>());
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
     std::vector<CoeffType> data;
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
         Planet::BranchingRatioNode<CoeffType> * node = new Planet::BranchingRatioNode<CoeffType>();
         std::string node_id;
         std::string node_pdf;
         std::vector<CoeffType> node_pdf_params;
  
         unsigned int nchan(0);
         std::vector<std::vector<CoeffType> > data_chan;
  
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
           data_chan.push_back(std::vector<CoeffType>());
  
           std::vector<std::string> tmp;
           Antioch::SplitString(reac.channels[ibr],";",tmp,false);//
           if(tmp.size() < 4)
           {
                std::cerr << "Error: 4 parameters minimum are required in this description:\n\t" << reac.channels[ibr] << std::endl;
                antioch_error();
           }
  //pdf
           shave_strings(tmp);
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
     const unsigned int n_species = reaction_set.n_species();
     for(unsigned int ibr = 0; ibr < reac.channels.size(); ibr++)
     {
         std::string equation = reac.channels[ibr].substr(0,reac.channels[ibr].find(";"));
         shave_string(equation);
  //reaction
         Antioch::Reaction<CoeffType>* my_rxn = Antioch::build_reaction<CoeffType>(n_species, equation, false, typeReaction, kineticsModel);
  
         Antioch::KineticsType<CoeffType>* rate = prob_reac.create_rate_constant(ibr,kineticsModel);
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
         shave_strings(reactants);
  
         std::vector<unsigned int> stoi;
         stoi.resize(reactants.size(),1);
         condense_molecule(stoi,reactants);
  
//checking reactants
         bool count_reaction(true);
         for(unsigned int ir = 0; ir < reactants.size(); ir++)
         {
            if(!reaction_set.chemical_mixture().active_species_name_map().count(reactants[ir]))
                count_reaction = false;
         }

         if(!count_reaction)
         {
             delete my_rxn;
             continue;
         }

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
         shave_strings(reactants);
         stoi.clear();
         stoi.resize(reactants.size(),1);
         condense_molecule(stoi,reactants);

         for(unsigned int ip = 0; ip < reactants.size(); ip++)
         {
            if(!reaction_set.chemical_mixture().active_species_name_map().count( reactants[ip]))
                      count_reaction = false;
         }
  
         if(!count_reaction)
         {
             delete my_rxn;
             continue;
         }

         for(unsigned int ip = 0; ip < reactants.size(); ip++)
         {
            my_rxn->add_product(reactants[ip],reaction_set.chemical_mixture().active_species_name_map().find( reactants[ip])->second,stoi[ip]);
         }
        
         reaction_set.add_reaction(my_rxn);
     }
     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::read_flyby_infos(const std::vector<std::string>& neutrals,
                                                                                        CoeffType& dens_tot, std::vector<CoeffType>& molar_frac,
                                                                                        CoeffType& chi, CoeffType& K0,
                                                                                        const std::string& file_flyby,
                                                                                        const std::string& root_input) const
  {
//TODO, change to a GetPot
    std::ifstream flyby(file_flyby.c_str());
    if( !flyby )
      {
        std::cerr << "Could not open file " << file_flyby << std::endl;
        antioch_error();
      }

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
            K0 = std::atof(line.c_str()); //cm2/s
          }else if(line.substr(0,line.find(':')) == "file mix")
          {
            line = line.substr(line.find(':') + 1,std::string::npos);
            shave_string(line);
            fill_molar_frac(neutrals,molar_frac,line,root_input);
          }
      }

    flyby.close();
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::read_cross_section( const std::string &file,
                            VectorCoeffType &lambda, VectorCoeffType &sigma ) const
  {
    std::string line;
    std::ifstream sig_f(file);
    if( !sig_f)
      {
        std::cerr << "Could not open file " << file << std::endl;
        antioch_error();
      }
    getline(sig_f,line);
    while(!sig_f.eof())
      {
        CoeffType wv(-1.),sigt(-1.);
        sig_f >> wv >> sigt;
        if(!getline(sig_f,line))break;
        if(wv < 0.)break;
        lambda.push_back(wv);//A
        sigma.push_back(sigt);//cm-2/A
      }
    sig_f.close();

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::read_hv_flux(Antioch::ParticleFlux<VectorCoeffType> &phy1AU, 
                                                                                    Antioch::ParticleFlux<VectorCoeffType> &phy_at_top, 
                                                                                    const std::string &file) const
  {
    std::string line;
    std::ifstream flux_1AU(file);
    if( !flux_1AU)
      {
        std::cerr << "Could not open file " << file << std::endl;
        antioch_error();
      }
    getline(flux_1AU,line);
    VectorCoeffType lambda,flux;
//TODO unit management !! use SwRI only
    while(!flux_1AU.eof())
      {
        CoeffType wv(-1),ir(-1),dirr(-1);
        flux_1AU >> wv >> ir >> dirr;
        if(wv < 0.)break;
        lambda.push_back(wv);// * 10.L);//nm -> A
        flux.push_back(ir);/* * 1e3L * (wv*1e-9L) / (Antioch::Constants::Planck_constant<CoeffType>() *
                                                   Antioch::Constants::light_celerity<CoeffType>()));//W/m2/nm -> J/s/cm2/A -> s-1/cm-2/A*/
      }
    flux_1AU.close();

//if in reverse order
   if(lambda.back() < lambda.front())
   {
      VectorCoeffType tmp_l(lambda.size());
      VectorCoeffType tmp_p(flux.size());
      for(unsigned int i = 0; i < lambda.size(); i++)
      {
         tmp_l[lambda.size() - 1 - i] = lambda[i];
         tmp_p[lambda.size() - 1 - i] = flux[i];
      }
      lambda.clear();
      flux.clear();
      lambda = tmp_l;
      flux = tmp_p;
   }

    phy1AU.set_abscissa(lambda);
    phy1AU.set_flux(flux);

//TODO: generalize it with input file
   CoeffType d = Constants::Saturn::d_Sun<CoeffType>();
   for(unsigned int i = 0; i < lambda.size(); i++)
   {
      flux[i] /= (d * d);
   }
   phy_at_top.set_abscissa(lambda);
   phy_at_top.set_flux(flux);

   return;

  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::read_neutral_characteristics(const std::vector<std::string>& neutrals,
                                      VectorCoeffType& tc,
                                      VectorCoeffType& hard_sphere_radius,
                                      std::vector<std::vector<std::vector<CoeffType> > >& bin_diff_data,
                                      std::vector<std::vector<DiffusionType> >& bin_diff_model,
                                      const std::string & file_neutral_charac) const
  {
    for(unsigned int m = 0; m < bin_diff_data.size(); m++) //medium species
      {
        bin_diff_data[m].resize(neutrals.size());
        for(unsigned int s = 0; s < neutrals.size(); s++)
        {
            bin_diff_data[m][s].resize(2,0.);
        }
        bin_diff_model[m].resize(neutrals.size(),Planet::DiffusionType::NoData); // no data by default
      }

    tc.resize(neutrals.size(),0.L);
    hard_sphere_radius.resize(neutrals.size(),0.L);

    std::ifstream neu(file_neutral_charac.c_str());
    if( !neu)
      {
        std::cerr << "Could not open file " << file_neutral_charac << std::endl;
        antioch_error();
      }
    std::string line;
    getline(neu,line); //first line
    while(!neu.eof())
      {
        std::string name;
        CoeffType mass,AN2,sN2,ACH4,sCH4,alpha,enth,hsr,mr;
        neu >> name >> mass >> AN2 >> sN2 >> ACH4 >> sCH4 >> alpha >> enth >> hsr >> mr;
        for(unsigned int s = 0; s < neutrals.size(); s++)
          {
            if(name == neutrals[s])
              {
                tc[s] = alpha; // no unit
                bin_diff_model[0][s] = Planet::DiffusionType::Wilson; //N2
                bin_diff_model[1][s] = Planet::DiffusionType::Wilson; //CH4
                bin_diff_data[0][s][0] = AN2;  //N2 - A  -- cm2/s
                bin_diff_data[0][s][1] = sN2;  //N2 - s  -- no unit
                bin_diff_data[1][s][0] = ACH4; //CH4 - A  -- cm2/s
                bin_diff_data[1][s][1] = sCH4; //CH4 - s  -- no unit
                hard_sphere_radius[s] = hsr;
                break;
              }
          }
      }
    neu.close();

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::fill_molar_frac(const std::vector<std::string> &neutrals, VectorCoeffType& molar_frac,
                                                                                       const std::string &file, const std::string &root_input) const
  {
    std::ifstream frac((root_input + file).c_str());
    if( !frac )
      {
        std::cerr << "Could not open file " << root_input + file << std::endl;
        antioch_error();
      }
    std::string line;
    getline(frac,line);//title
    molar_frac.resize(neutrals.size(),0.L);
    while(!frac.eof())
      {
        std::string neutral;
        CoeffType frac_600, frac_1050, Dfrac_1050;
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

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::parse_equation(std::vector<std::string> &reactants, std::vector<std::string> &products,
                                                                                      std::string &line, bool &skip, const Antioch::ChemicalMixture<CoeffType>& chem_mixture,
                                                                                      std::string &equation, std::vector<unsigned int> &stoi_reac,
                                                                                      std::vector<unsigned int> &stoi_prod) const
  {
    std::vector<std::string> out;
    Antioch::SplitString(line,";",out,false);
    if(out.size() != 2)
    {
        std::cerr << "Error: badly shaped line\n\t" << line << std::endl;
        antioch_error();
    }
    equation = out[0];
    std::string parameters(out[1]);
    ///// equation
    std::vector<std::string> molecules;
    Antioch::SplitString(equation,"->",molecules,false);
    if(molecules.size() != 2)
    {
        std::cerr << "Error: badly shaped equation\n\t" << equation << std::endl;
        antioch_error();
    }

    Antioch::SplitString(molecules[0],"+",reactants,false);
    if(reactants.empty())reactants.push_back(molecules[0]);
    Antioch::SplitString(molecules[1],"+",products,false);
    if(products.empty())products.push_back(molecules[1]);
    this->shave_strings(reactants);
    this->shave_strings(products);
    stoi_reac.resize(reactants.size(),1);
    stoi_prod.resize(products.size(),1);

    this->condense_molecule(stoi_reac,reactants);
    this->condense_molecule(stoi_prod,products);

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

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::condense_molecule(std::vector<unsigned int> &stoi, std::vector<std::string> &mol) const
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

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::read_photochemistry_reac(const std::string &hv_file, const std::string &reac,
                                                                                                Antioch::ReactionSet<CoeffType> &neutral_reaction_set) const
  {
    Antioch::KineticsModel::KineticsModel kineticsModel(Antioch::KineticsModel::PHOTOCHEM);
    Antioch::ReactionType::ReactionType reactionType(Antioch::ReactionType::ELEMENTARY);

    const Antioch::ChemicalMixture<CoeffType>& chem_mixture = neutral_reaction_set.chemical_mixture();

    std::ifstream data(hv_file.c_str());
    std::string line;
    getline(data,line);
    std::vector<std::string> out;
    Antioch::SplitString(line," ",out,false);
    unsigned int nbr = out.size();
    if(nbr == 0)
    {
        std::cerr << "Error: badly shaped line\n\t" << line << std::endl;
        antioch_error();
    }

    MatrixCoeffType datas;
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
        CoeffType lambda(-1), total(-1);
        VectorCoeffType sigmas;
        sigmas.resize(nbr - 2,0.L);
        data >> lambda >> total;
        if(lambda < 0.)break;
        for(unsigned int ibr = 0; ibr < nbr - 2; ibr++)data >> sigmas[ibr]; //only br here
        datas[0].push_back(lambda);
        for(unsigned int ibr = 0; ibr < nbr - 2; ibr++)datas[ibr + 1].push_back(sigmas[ibr]); // + \lambda
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
        Antioch::Reaction<CoeffType> * reaction = Antioch::build_reaction<CoeffType>(chem_mixture.n_species(), equation, false, reactionType, kineticsModel);

        reaction->add_reactant( reac,chem_mixture.active_species_name_map().find(reac)->second,1);
        for(unsigned int ip = 0; ip < produc[ibr].size(); ip++)
          {
            reaction->add_product( produc[ibr][ip],chem_mixture.active_species_name_map().find(produc[ibr][ip])->second,stoi_prod[ibr][ip]);
          }

        VectorCoeffType dataf;
        dataf.resize(datas[0].size() + datas[ibr].size(),0.);
        int istep(1);
        int start(0);
        if(datas[0].back() < datas[0].front())
        {
           istep = -1;
           start = datas[0].size() - 1;
        }
        unsigned int j(0);
        for(int i = start; i < (int)datas[0].size() && i > -1; i += istep) //cs first
        {
          dataf[j] = datas[ibr + 1][i];
          j++;
        }
        for(int i = start; i < (int)datas[0].size() && i > -1; i += istep) //\lambda then
        {
           dataf[j] = datas[0][i];
           j++;
        }

        Antioch::KineticsType<CoeffType,VectorCoeffType> * rate  = Antioch::build_rate<CoeffType,VectorCoeffType>(dataf,kineticsModel); //kinetics rate
        reaction->add_forward_rate(rate);

        neutral_reaction_set.add_reaction(reaction);
      }

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::shave_string(std::string &str) const
  {
    // Trim from the left
    str.erase(str.begin(), std::find_if(str.begin(), str.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));

    // Trim from the right
    str.erase(std::find_if(str.rbegin(), str.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), str.end());

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::shave_strings(std::vector<std::string> &stock) const
  {
    for(unsigned int i = 0; i < stock.size(); i++)
      {
        this->shave_string(stock[i]);
      }

    return;
  }


  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::sanity_check_chemical_system() const
  {
    // checks neutral and ionic system are chemically balanced, if not, 
    // write to std::cerr the names of unbalanced molecules and sends an antioch_error()
    bool balanced(true);
    balanced = balanced && check_chemical_balance(*_neutral_reaction_set,_neutral_reaction_set->chemical_mixture().species_list(),"in photochemical model");
    balanced = balanced && check_chemical_balance(*_ionic_reaction_set,_ss_species,"in ionospheric model");

   // if(!balanced)antioch_error();
    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  bool PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::check_chemical_balance(const Antioch::ReactionSet<CoeffType> &reaction_set,
                                                                                              const std::vector<Antioch::Species> &species,
                                                                                              const std::string & help) const
  {
      std::vector<unsigned int> prod(species.size(),0);
      std::vector<unsigned int> loss(species.size(),0);

      std::map<unsigned int, unsigned int> species_map;

      for(unsigned int s = 0; s < reaction_set.chemical_mixture().species_list().size(); s++)
      {
          for(unsigned int j = 0; j < species.size(); j++)
          {
              if(reaction_set.chemical_mixture().species_list()[s] == species[j])
              {
                 species_map[s] = j;
                 break;
              }
          }
      }

      for(unsigned int rxn = 0; rxn < reaction_set.n_reactions(); rxn++)
      {
        //reactants
          for(unsigned int r = 0; r < reaction_set.reaction(rxn).n_reactants(); r++)
          {
              if(species_map.count(reaction_set.reaction(rxn).reactant_id(r)))
                             loss[species_map.at(reaction_set.reaction(rxn).reactant_id(r))]++;
          }

        //products
          for(unsigned int p = 0; p < reaction_set.reaction(rxn).n_products(); p++)
          {
              if(species_map.count(reaction_set.reaction(rxn).product_id(p)))
                             prod[species_map.at(reaction_set.reaction(rxn).product_id(p))]++;
          }
      }

      bool balanced(true);
      for(unsigned int s = 0; s < species.size(); s++)
      {
          if(prod[s] == 0 && loss[s] == 0)
          {
              std::cerr << "Species " << reaction_set.chemical_mixture().species_inverse_name_map().at(species[s])
                        << " is unreactive " << help << std::endl;
          }
          if((prod[s] == 0 && loss[s] != 0) || (prod[s] != 0 && loss[s] == 0))
          {
              std::cerr << "Species " << reaction_set.chemical_mixture().species_inverse_name_map().at(species[s])
                        << " is not balanced, it is only ";
              (prod[s] == 0)?std::cerr << "consumed ":std::cerr << "produced ";
              std::cerr << help << std::endl;
              balanced = false;
          }
      }

      return balanced;
  }


  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>& PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::composition() const
  {
    return *_composition;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  const Antioch::ReactionSet<CoeffType>& PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::neutral_reaction_set() const
  {
    return *_neutral_reaction_set;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  const Antioch::ReactionSet<CoeffType>& PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::ionic_reaction_set() const
  {
    return *_ionic_reaction_set;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  const PhotonOpacity<CoeffType,VectorCoeffType>& PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::tau() const
  {
    return *_tau;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  const std::vector<std::vector<BinaryDiffusion<CoeffType> > >& PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::bin_diff_coeff() const
  {
    return _bin_diff_coeff;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  const AtmosphericTemperature<CoeffType,VectorCoeffType>& PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::temperature() const
  {
    return *_temperature;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  const VectorCoeffType& PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::lambda_hv() const
  {
    return _phy1AU.abscissa();
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  const VectorCoeffType& PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::phy1AU() const
  {
    return _phy1AU.flux();
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  const Antioch::ParticleFlux<VectorCoeffType> & PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::phy_at_top() const
  {
    return _phy_at_top;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  const std::vector<std::string>& PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::medium() const
  {
    return _medium;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  CoeffType PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::K0() const
  {
    return _K0;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  CoeffType PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::scaling_factor() const
  {
    return _scaling_factor;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  CoeffType PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::dupper_boundary_neumann_s_dn_i(unsigned int s, unsigned int i) const
  {
    return (i == s)?_composition->upper_boundary_velocity(s):Antioch::zero_clone(_K0);
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  const std::vector<Antioch::Species> PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::ss_species() const
  {
    return _ss_species;
  }

} // end namespace Planet

#endif // PLANET_PLANET_PHYSICS_HELPER_H
