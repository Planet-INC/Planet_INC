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

#ifndef PLANET_PLANET_PHYSICS_EVALUATOR_H
#define PLANET_PLANET_PHYSICS_EVALUATOR_H

namespace Planet
{

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  class PlanetPhysicsEvaluator
  {
  public:

    PlanetPhysicsEvaluator( PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>& helper );
    ~PlanetPhysicsEvaluator();

    //computes omega_dot and omega
    template<typename StateType, typename VectorStateType>
    void compute(const VectorStateType & molar_concentrations,
                 const VectorStateType & dmolar_concentrations_dz,
                 const StateType & z);

    libMesh::Real diffusion_term(unsigned int s) const;

    libMesh::Real chemical_term(unsigned int s)  const;

    //!fills molar_concentrations_first_guess with barometric equation
    template<typename StateType, typename VectorStateType>
    void first_guess(VectorStateType & molar_concentrations_first_guess, const StateType z) const;

  protected:

    template<typename VectorStateType, typename StateType>
    void update_cache(const VectorStateType &molar_concentrations, const StateType &z);

    void cache_recompute();

    //! uses compo.barometric_density(z);
    template <typename StateType>
    const VectorCoeffType get_cache(const StateType &z) const;

    const AtmosphericMixture<CoeffType,VectorCoeffType,MatrixCoeffType>& _composition;

    Antioch::KineticsEvaluator<CoeffType> _neutral_kinetics;
    Antioch::KineticsEvaluator<CoeffType> _ionic_kinetics;

    Planet::PhotonEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> _photon;

    AtmosphericKinetics<CoeffType,VectorCoeffType,MatrixCoeffType> _kinetics;
    DiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> _diffusion;

    Planet::MolecularDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> _molecular_diffusion;
    Planet::EddyDiffusionEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType> _eddy_diffusion;

    VectorCoeffType _omegas;
    VectorCoeffType _omegas_dots;

    MatrixCoeffType _cache_composition;
    VectorCoeffType _cache_altitudes;
    std::map<CoeffType,VectorCoeffType> _cache;

  private:

    PlanetPhysicsEvaluator();

  };

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  PlanetPhysicsEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::PlanetPhysicsEvaluator( PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>& helper )
    : _composition(helper.composition()),
      _neutral_kinetics(helper.neutral_reaction_set(),0), /*! \todo generalize 0 for other types*/
      _ionic_kinetics(helper.ionic_reaction_set(),0), /*! \todo generalize 0 for other types*/
      _photon(helper.tau(),_composition),
      _molecular_diffusion(helper.bin_diff_coeff(),_composition,helper.temperature()),
      _eddy_diffusion(_composition,helper.K0()),
      _kinetics(_neutral_kinetics,_ionic_kinetics,helper.temperature(),_photon,_composition),
      _diffusion(_molecular_diffusion,_eddy_diffusion,_composition,helper.temperature())
  {
    _omegas.resize(_kinetics->neutral_kinetics().reaction_set().n_species());
    _omegas_dots.resize(_kinetics->neutral_kinetics().reaction_set().n_species());

    _photon.set_photon_flux_at_top(helper.lambda_hv(), helper.phy1AU(), Planet::Constants::Saturn::d_Sun<CoeffType>());

    /*! \todo This call to set_particle_flux is going to kill thread safety because
              it's resetting stuff in the ReactionSet */
    helper.neutral_reaction_set().set_particle_flux(photon.photon_flux_ptr()); // reactions know the solar flux

    molecular_diffusion.set_medium_species(helper.medium());

    return;
  }


  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType>
  void PlanetPhysicsEvaluator<CoeffType,VectorCoeffType,MatrixCoeffType>::compute(const VectorStateType & molar_concentrations,
                                                                                  const VectorStateType & dmolar_concentrations_dz,
                                                                                  const StateType & z)
  {
   _diffusion->diffusion(molar_concentrations,dmolar_concentrations_dz,z,_omegas);
   _kinetics->chemical_rate(molar_concentrations,this->get_cache(z),z,_omegas_dots);

   this->update_cache(molar_concentrations,z);

    return;
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

} // end namespace Planet

#endif // PLANET_PLANET_PHYSICS_EVALUATOR_H
