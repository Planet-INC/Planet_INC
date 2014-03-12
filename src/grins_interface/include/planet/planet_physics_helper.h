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

    template <typename StateType, typename VectorStateType, typename MatrixStateType>
    void set_kinetics(AtmosphericKinetics<StateType,VectorStateType,MatrixStateType> *kinetics);

    template <typename StateType, typename VectorStateType, typename MatrixStateType>
    void set_diffusion(DiffusionEvaluator<StateType,VectorStateType,MatrixStateType> *diffusion);

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

    /*! Convenience method to hide all the construction code for
        composition, kinetics, and diffusion */
    void build_helper( const GetPot& input );

    // Additional data structures that need to be cached
    AtmosphericTemperature<CoeffType,VectorCoeffType>* _temperature;
    Antioch::ChemicalMixture<CoeffType>* _neutral_species;
    Antioch::ChemicalMixture<CoeffType>* _ionic_species;

    Antioch::ReactionSet<CoeffType>* _neutral_reaction_set;
    Antioch::ReactionSet<CoeffType>* _ionic_reaction_set;
    Antioch::ReactionSet<CoeffType>* _neut_reac_theo;

    // Helper functions for build_helper
    void read_temperature(VectorCoeffType& T0, VectorCoeffType& Tz, const std::string& file) const;

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
      _neut_reac_theo(NULL)
  {
    this->build_helper(input);

    _omegas.resize(_kinetics->neutral_kinetics().reaction_set().n_species());
    _omegas_dots.resize(_kinetics->neutral_kinetics().reaction_set().n_species());

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

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template <typename StateType, typename VectorStateType, typename MatrixStateType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::set_kinetics(AtmosphericKinetics<StateType,VectorStateType,MatrixStateType> *kinetics)
  {
     _kinetics = kinetics;
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template <typename StateType, typename VectorStateType, typename MatrixStateType>
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::set_diffusion(DiffusionEvaluator <StateType,VectorStateType,MatrixStateType> *diffusion)
  {
     _diffusion = diffusion;
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
  void PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::build_helper(const GetPot& input)
  {
    // Read temperature profile
    std::string input_T = input( "Planet/temperature_file", "DIE!" );
    std::vector<CoeffType> T0,Tz;
    this->read_temperature(T0,Tz,input_T);
    _temperature = new AtmosphericTemperature<CoeffType,VectorCoeffType>(T0,T0,Tz,Tz);

    // Read neutral and ionic species from input
    unsigned int n_neutral = input.vector_variable_size("Planet/neutral_species");
    unsigned int n_ionic = input.vector_variable_size("Planet/ionic_species");

    std::vector<std::string> neutrals(n_neutral);
    std::vector<std::string> ions(n_ionic);

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

    // Build up reaction sets
    _neutral_reaction_set = new Antioch::ReactionSet<CoeffType>(*_neutral_species);
    _ionic_reaction_set = new Antioch::ReactionSet<CoeffType>(*_ionic_species);
    _neut_reac_theo = new Antioch::ReactionSet<CoeffType>(*_neutral_species);

    _composition = new AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>( _neutral_species, _ionic_species, _temperature );

    /*
    _diffusion = new DiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>( _molecular_diffusion, _eddy_diffusion, _composition, _temperature );

    _kinetics = new AtmosphericKinetics<CoeffType,VectorCoeffType,MatrixCoeffType>( _neutral_kinetics, _ionic_kinetics, _temperature, _photon, *_composition );
    */

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
  const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>& PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>::composition() const
  {
    return *_composition;
  }

} // end namespace Planet

#endif // PLANET_PLANET_PHYSICS_HELPER_H
