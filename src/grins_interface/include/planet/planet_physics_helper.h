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

// libMesh
#include "libmesh/libmesh_common.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"

namespace Planet
{

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  class PlanetPhysicsHelper
  {
  public:

    PlanetPhysicsHelper( const GetPot& input );

    PlanetPhysicsHelper(AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> *compo,
                        AtmosphericKinetics<CoeffType,VectorCoeffType,MatrixCoeffType > *kinetics = NULL,
                        DiffusionEvaluator <CoeffType,VectorCoeffType,MatrixCoeffType > *diffusion = NULL);

    ~PlanetPhysicsHelper();

    libMesh::Real diffusion_term(unsigned int s) const;

    libMesh::Real chemical_term(unsigned int s)  const;

    //computes omega_dot and omega
    template<typename StateType, typename VectorStateType>
    void compute(const VectorStateType & molar_concentrations,
                 const VectorStateType & dmolar_concentrations_dz,
                 const StateType & z);

    //!fills molar_concentrations_first_guess with barometric equation
    template<typename StateType, typename VectorStateType>
    void first_guess(VectorStateType & molar_concentrations_first_guess, const StateType z) const;

    //!fills lower boundary conditions
    template<typename VectorStateType>
    void lower_boundary_dirichlet(VectorStateType & lower_boundary) const;

    //!fills upper boundary conditions
    template<typename VectorStateType>
    void upper_boundary_neumann(VectorStateType & upper_boundary, const VectorStateType &molar_densities) const;

    const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>& composition() const;

    /*! \todo This should really be const. Need to fix up ParticleFlux stuff. */
    Antioch::ReactionSet<CoeffType>& neutral_reaction_set();

    const Antioch::ReactionSet<CoeffType>& ionic_reaction_set() const;

    PhotonOpacity<CoeffType,VectorCoeffType>& tau();

    const std::vector<std::vector<BinaryDiffusion<CoeffType> > >& bin_diff_coeff() const;

    const AtmosphericTemperature<CoeffType,VectorCoeffType>& temperature() const;

    const VectorCoeffType& lambda_hv() const;

    const VectorCoeffType& phy1AU() const;

    const std::vector<std::string>& medium() const;

    CoeffType K0() const;

  private:

    AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>*  _composition; //for first guess
    AtmosphericKinetics<CoeffType,VectorCoeffType,MatrixCoeffType>* _kinetics;
    DiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>* _diffusion;

    template<typename VectorStateType, typename StateType>
    void update_cache(const VectorStateType &molar_concentrations, const StateType &z);

    void cache_recompute();

    //! uses compo.barometric_density(z);
    template <typename StateType>
    const VectorCoeffType get_cache(const StateType &z) const;

    VectorCoeffType _omegas;
    VectorCoeffType _omegas_dots;
    MatrixCoeffType _cache_composition;
    VectorCoeffType _cache_altitudes;
    std::map<CoeffType,VectorCoeffType> _cache;

    // Additional data structures that need to be cached
    AtmosphericTemperature<CoeffType,VectorCoeffType>* _temperature;
    Antioch::ChemicalMixture<CoeffType>* _neutral_species;
    Antioch::ChemicalMixture<CoeffType>* _ionic_species;

    Antioch::ReactionSet<CoeffType>* _neutral_reaction_set;
    Antioch::ReactionSet<CoeffType>* _ionic_reaction_set;

    Chapman<CoeffType>* _chapman;

    std::vector<std::vector<BinaryDiffusion<CoeffType> > > _bin_diff_coeff;

    //! Parameter for diffusion
    /*! Needs to be cached because some Evaluators depend on this value */
    CoeffType _K0;

    PhotonOpacity<CoeffType,VectorCoeffType>* _tau;

    VectorCoeffType _lambda_hv;
    VectorCoeffType _phy1AU;

    std::vector<std::string> _medium;

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
                            CoeffType dens_tot, VectorCoeffType& tc);

    void build_diffusion( std::vector<std::vector<std::vector<CoeffType> > >& bin_diff_data,
                          std::vector<std::vector<DiffusionType> >& bin_diff_model,
                          const std::vector<std::string>& neutrals);

    // Helper functions for parsing data
    void read_temperature(VectorCoeffType& T0, VectorCoeffType& Tz, const std::string& file) const;

    void fill_neutral_reactions_elementary(const std::string &neutral_reactions_file,
                                           const std::string &N2_hv_file,
                                           const std::string &CH4_hv_file,
                                           Antioch::ReactionSet<CoeffType>& neutral_reaction_set ) const;

    void fill_neutral_reactions_falloff(const std::string &neutral_reactions_file,
                                        Antioch::ReactionSet<CoeffType>& neutral_reaction_set) const;

    void read_flyby_infos(const std::vector<std::string>& neutrals,
                          CoeffType& dens_tot, std::vector<CoeffType>& molar_frac,
                          CoeffType& chi, CoeffType& K0,
                          const std::string& file_flyby,
                          const std::string& root_input) const;

    void read_crossSection( const std::string &file, unsigned int nbr,
                            VectorCoeffType &lambda, VectorCoeffType &sigma ) const;

    void read_hv_flux(VectorCoeffType& lambda, VectorCoeffType& phy1AU, const std::string &file) const;

    void read_neutral_characteristics(const std::vector<std::string>& neutrals,
                                      VectorCoeffType& tc,
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
  };

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::PlanetPhysicsHelper( const GetPot& input )
    : _composition(NULL),
      _kinetics(NULL),
      _diffusion(NULL),
      _temperature(NULL),
      _neutral_species(NULL),
      _ionic_species(NULL),
      _neutral_reaction_set(NULL),
      _ionic_reaction_set(NULL),
      _chapman(NULL),
      _tau(NULL)
  {
    this->build(input);

    //_omegas.resize(_kinetics->neutral_kinetics().reaction_set().n_species());
    //_omegas_dots.resize(_kinetics->neutral_kinetics().reaction_set().n_species());

    return;
  }


  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::compute(const VectorStateType & molar_concentrations,
                                                               const VectorStateType & dmolar_concentrations_dz,
                                                               const StateType & z)
  {
   _diffusion->diffusion(molar_concentrations,dmolar_concentrations_dz,z,_omegas);
   _kinetics->chemical_rate(molar_concentrations,this->get_cache(z),z,_omegas_dots);

   this->update_cache(molar_concentrations,z);

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::PlanetPhysicsHelper(AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType> *compo,
                                                        AtmosphericKinetics<CoeffType,VectorCoeffType,MatrixCoeffType > *kinetics,
                                                        DiffusionEvaluator <CoeffType,VectorCoeffType,MatrixCoeffType > *diffusion):
        _kinetics(kinetics),
        _diffusion(diffusion),
        _composition(compo)
  {
    _omegas.resize(_kinetics->neutral_kinetics().reaction_set().n_species());
    _omegas_dots.resize(_kinetics->neutral_kinetics().reaction_set().n_species());
    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::~PlanetPhysicsHelper()
  {
    /*
      delete _tau;
    delete _chapman;
    delete _ionic_reaction_set;
    delete _neutral_reaction_set;
    delete _ionic_species;
    delete _neutral_species;
    delete _temperature;
    delete _composition;
    */

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  libMesh::Real PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::diffusion_term(unsigned int s) const
  {
    return _omegas[s];
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  libMesh::Real PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::chemical_term(unsigned int s) const
  {
    return _omegas_dots[s];
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template <typename StateType>
  const VectorCoeffType PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::get_cache(const StateType &z) const
  {
     if(!_cache.count(z))
     {
        VectorCoeffType first_sum_guess = Antioch::zero_clone(_composition->neutral_molar_fraction_bottom());
        _composition->first_guess_densities_sum(z,first_sum_guess);
        return first_sum_guess;
     }else
     {
        return _cache.at(z);
     }
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template <typename VectorStateType, typename StateType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::update_cache(const VectorStateType &molar_concentrations, const StateType &z)
  {
    bool recompute(true);
    if(!_cache.count(z))
    {
        recompute = false;
       _cache[z] = get_cache(z);
    }

     _cache_composition.push_back(molar_concentrations);
     _cache_altitudes.push_back(z);
     if(recompute && _cache_composition.size() == _cache.size())this->cache_recompute();
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::cache_recompute()
  {
   //from highest altitude to lowest altitude
   unsigned int istart(0);
   int istep(1);
   if(_cache_altitudes.back() > _cache_altitudes.front())
   {
      istep = -1;
      istart = _cache_altitudes.size() - 1;
   }
 

   //sum densities are sdens_{i} = n(z_{i+1}) * (z_{i+1} - z_{i}), top composition is useless
   for(unsigned int i = 1; i < _cache_altitudes.size(); i++)
   {
      unsigned int j = istart + istep * i;
      unsigned int jbottom = istart + istep * (i - 1);
      for(unsigned int s = 0; s < _cache_composition[j].size(); s++)
      {
        _cache.at(_cache_altitudes[j])[s] = _cache.at(_cache_altitudes[jbottom])[s] + 
                                           _cache_composition[j][s] * (_cache_altitudes[j] - _cache_altitudes[jbottom]);
      }
   }

    _cache_composition.clear();
    _cache_altitudes.clear();
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::first_guess(VectorStateType & molar_concentrations_first_guess, const StateType z) const
  {
      _composition->first_guess_densities(z,molar_concentrations_first_guess);
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename VectorStateType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::lower_boundary_dirichlet(VectorStateType & lower_boundary) const
  {
      _composition->lower_boundary_concentrations(lower_boundary);
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename VectorStateType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::upper_boundary_neumann(VectorStateType & upper_boundary, const VectorStateType &molar_densities) const
  {
      _composition->upper_boundary_fluxes(upper_boundary, molar_densities);
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

    /*! \todo What are these magic numbers? */
    bin_diff_data.resize(_medium.size());
    bin_diff_model.resize(_medium.size());
    _bin_diff_coeff.resize(_medium.size());

    VectorCoeffType tc;

    this->read_neutral_characteristics( neutrals, tc, bin_diff_data, bin_diff_model,
                                        file_neutral_charac);

    _chapman = new Chapman<CoeffType>(chi);

    this->build_diffusion(bin_diff_data, bin_diff_model, neutrals);

    // Must be called after: build_species, chapman
    this->build_opacity(input);

    this->build_temperature(input);

    // Must be called after: build_species
    this->build_reaction_sets(input);

    // Must be called after: build_temperature, build_species
    this->build_composition(input, molar_frac, dens_tot, tc);

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

    if( !input.have_variable("Planet/ionic_species") )
      {
        std::cerr << "Error: ionic_species not found in input file!" << std::endl;
        antioch_error();
      }

    // Read neutral and ionic species from input
    unsigned int n_neutral = input.vector_variable_size("Planet/neutral_species");
    unsigned int n_ionic = input.vector_variable_size("Planet/ionic_species");

    neutrals.resize(n_neutral);
    ions.resize(n_ionic);

    for( unsigned int s = 0; s < n_neutral; s++ )
      {
        neutrals[s] = input("Planet/neutral_species", "DIE!", s);
      }

    for( unsigned int s = 0; s < n_neutral; s++ )
      {
        ions[s] = input("Planet/ionic_species", "DIE!", s);
      }

    _neutral_species = new Antioch::ChemicalMixture<CoeffType>(neutrals);
    _ionic_species = new Antioch::ChemicalMixture<CoeffType>(ions);

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::build_opacity( const GetPot& input )
  {
    _tau = new PhotonOpacity<CoeffType,VectorCoeffType>(*_chapman);

    std::vector<CoeffType> lambda_N2, sigma_N2;
    std::vector<CoeffType> lambda_CH4, sigma_CH4;

    if( !input.have_variable("Planet/input_N2") )
      {
        std::cerr << "Error: input_N2 not found in input file!" << std::endl;
        antioch_error();
      }

    if( !input.have_variable("Planet/input_CH4") )
      {
        std::cerr << "Error: input_CH4 not found in input file!" << std::endl;
        antioch_error();
      }

    if( !input.have_variable("Planet/input_hv") )
      {
        std::cerr << "Error: input_hv not found in input file!" << std::endl;
        antioch_error();
      }

    std::string input_hv  = input("Planet/input_hv", "DIE!" );
    std::string input_N2  = input("Planet/input_N2", "DIE!" );
    std::string input_CH4 = input("Planet/input_CH4", "DIE!" );

    this->read_hv_flux(_lambda_hv, _phy1AU, input_hv);

    /*! \todo What are these magic numbers? */
    this->read_crossSection(input_N2,  3, lambda_N2,  sigma_N2);
    this->read_crossSection(input_CH4, 9, lambda_CH4, sigma_CH4);

    /* here only N2 and CH4 absorb */
    _tau->add_cross_section( lambda_N2, sigma_N2, Antioch::Species::N2,
                             _neutral_species->active_species_name_map().at("N2") );

    _tau->add_cross_section( lambda_CH4, sigma_CH4, Antioch::Species::CH4,
                             _neutral_species->active_species_name_map().at("CH4") );

    _tau->update_cross_section(_lambda_hv);

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

    if( !input.have_variable("Planet/input_N2") )
      {
        std::cerr << "Error: could not find input_N2 filename!" << std::endl;
        antioch_error();
      }

    if( !input.have_variable("Planet/input_CH4") )
      {
        std::cerr << "Error: could not find input_CH4 filename!" << std::endl;
        antioch_error();
      }

    std::string input_reactions_elem = input( "Planet/input_reactions_elem", "DIE!" );
    std::string input_reactions_fall = input( "Planet/input_reactions_fall", "DIE!" );
    std::string input_N2  = input( "Planet/input_N2", "DIE!" );
    std::string input_CH4 = input( "Planet/input_CH4", "DIE!" );

    // here read the reactions, Antioch will take care of it once hdf5, no ionic reactions there
    //Kooij / Arrhenius + photochem
    this->fill_neutral_reactions_elementary(input_reactions_elem, input_N2, input_CH4, *_neutral_reaction_set);

    this->fill_neutral_reactions_falloff(input_reactions_fall, *_neutral_reaction_set);

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::build_composition( const GetPot& input, VectorCoeffType& molar_frac, CoeffType dens_tot, VectorCoeffType& tc )
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

    for(unsigned int n = 0; n < 2; n++)
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
    getline(temp,line);
    while(!temp.eof())
      {
        CoeffType t,tz,dt,dtz;
        temp >> t >> tz >> dt >> dtz;
        T0.push_back(t);
        Tz.push_back(tz);
      }
    temp.close();

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::fill_neutral_reactions_elementary(const std::string &neutral_reactions_file,
                                                                                                         const std::string &N2_hv_file,
                                                                                                         const std::string &CH4_hv_file,
                                                                                                         Antioch::ReactionSet<CoeffType> &neutral_reaction_set) const
  {
    //here only simple ones: bimol Kooij/Arrhenius model
    std::ifstream data(neutral_reactions_file.c_str());
    std::string line;
    getline(data,line); //title
    Antioch::KineticsModel::KineticsModel kineticsModel(Antioch::KineticsModel::KOOIJ);
    Antioch::ReactionType::ReactionType reactionType(Antioch::ReactionType::ELEMENTARY);
    const Antioch::ChemicalMixture<CoeffType>& chem_mixture = neutral_reaction_set.chemical_mixture();
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
        if(dataf[1] == 0.) //Arrhenius
          {
            kineticsModel = Antioch::KineticsModel::ARRHENIUS;
          }else
          {
            dataf.push_back(std::atof(str_data[1].c_str())); //beta
          }
        dataf.push_back(std::atof(str_data[2].c_str()));//Ea

        if(kineticsModel == Antioch::KineticsModel::KOOIJ)dataf.push_back(1.); //Tref
        dataf.push_back(Antioch::Constants::R_universal<CoeffType>()*1e-3); //scale (R in J/mol/K)


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
    //now the photochemical ones
    read_photochemistry_reac(N2_hv_file, "N2", neutral_reaction_set);
    read_photochemistry_reac(CH4_hv_file, "CH4", neutral_reaction_set);
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::fill_neutral_reactions_falloff(const std::string &neutral_reactions_file,
                                                                                                      Antioch::ReactionSet<CoeffType> &neutral_reaction_set) const
  {
    //Lindemann
    std::ifstream data(neutral_reactions_file.c_str());
    std::string line;
    getline(data,line); //title
    getline(data,line); //title
    getline(data,line); //title
    Antioch::KineticsModel::KineticsModel kineticsModel(Antioch::KineticsModel::KOOIJ);
    Antioch::ReactionType::ReactionType reactionType(Antioch::ReactionType::LINDEMANN_FALLOFF);
    const Antioch::ChemicalMixture<CoeffType>& chem_mixture = neutral_reaction_set.chemical_mixture();
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
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::read_flyby_infos(const std::vector<std::string>& neutrals,
                                                                                        CoeffType& dens_tot, std::vector<CoeffType>& molar_frac,
                                                                                        CoeffType& chi, CoeffType& K0,
                                                                                        const std::string& file_flyby,
                                                                                        const std::string& root_input) const
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

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::read_crossSection( const std::string &file, unsigned int nbr,
                            VectorCoeffType &lambda, VectorCoeffType &sigma ) const
  {
    std::string line;
    std::ifstream sig_f(file);
    getline(sig_f,line);
    while(!sig_f.eof())
      {
        CoeffType wv,sigt,sigbr;
        sig_f >> wv >> sigt;
        for(unsigned int i = 0; i < nbr; i++)sig_f >> sigbr;
        lambda.push_back(wv);//A
        sigma.push_back(sigt);//cm-2/A
      }
    sig_f.close();

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::read_hv_flux(VectorCoeffType &lambda, VectorCoeffType &phy1AU, const std::string &file) const
  {
    std::string line;
    std::ifstream flux_1AU(file);
    getline(flux_1AU,line);
    while(!flux_1AU.eof())
      {
        CoeffType wv,ir,dirr;
        flux_1AU >> wv >> ir >> dirr;
        if(!lambda.empty() && wv == lambda.back())continue;
        lambda.push_back(wv);// * 10.L);//nm -> A
        phy1AU.push_back(ir);/* * 1e3L * (wv*1e-9L) / (Antioch::Constants::Planck_constant<CoeffType>() *
                                                   Antioch::Constants::light_celerity<CoeffType>()));//W/m2/nm -> J/s/cm2/A -> s-1/cm-2/A*/
      }
    flux_1AU.close();
    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::read_neutral_characteristics(const std::vector<std::string>& neutrals,
                                      VectorCoeffType& tc,
                                      std::vector<std::vector<std::vector<CoeffType> > >& bin_diff_data,
                                      std::vector<std::vector<DiffusionType> >& bin_diff_model,
                                      const std::string & file_neutral_charac) const
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
        CoeffType mass,AN2,sN2,ACH4,sCH4,alpha,enth,hsr,mr;
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

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::fill_molar_frac(const std::vector<std::string> &neutrals, VectorCoeffType& molar_frac,
                                                                                       const std::string &file, const std::string &root_input) const
  {
    std::ifstream frac((root_input + file).c_str());
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
    if(out.size() != 2)antioch_error();
    equation = out[0];
    std::string parameters(out[1]);
    ///// equation
    std::vector<std::string> molecules;
    Antioch::SplitString(equation,"->",molecules,false);
    if(molecules.size() != 2)antioch_error();

    Antioch::SplitString(molecules[0],"+",reactants,false);
    Antioch::SplitString(molecules[1],"+",products,false);
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
    if(nbr == 0)antioch_error();

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
        CoeffType lambda, total;
        VectorCoeffType sigmas;
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
        Antioch::Reaction<CoeffType> * reaction = Antioch::build_reaction<CoeffType>(chem_mixture.n_species(), equation, false, reactionType, kineticsModel);

        reaction->add_reactant( reac,chem_mixture.active_species_name_map().find(reac)->second,1);
        for(unsigned int ip = 0; ip < produc[ibr].size(); ip++)
          {
            reaction->add_product( produc[ibr][ip],chem_mixture.active_species_name_map().find(produc[ibr][ip])->second,stoi_prod[ibr][ip]);
          }

        VectorCoeffType dataf = datas[0];
        for(unsigned int i = 0; i < datas[ibr].size(); i++)
          {
            dataf.push_back(datas[ibr][i]);
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
    while(str[0] == ' ')str.erase(0,1);
    while(str[str.size() - 1] == ' ')str.erase(str.size() - 1,1);

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
  const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>& PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::composition() const
  {
    return *_composition;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  Antioch::ReactionSet<CoeffType>& PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::neutral_reaction_set()
  {
    return *_neutral_reaction_set;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  const Antioch::ReactionSet<CoeffType>& PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::ionic_reaction_set() const
  {
    return *_ionic_reaction_set;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  PhotonOpacity<CoeffType,VectorCoeffType>& PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::tau()
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
    return _lambda_hv;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  const VectorCoeffType& PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::phy1AU() const
  {
    return _phy1AU;
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

} // end namespace Planet

#endif // PLANET_PLANET_PHYSICS_HELPER_H
