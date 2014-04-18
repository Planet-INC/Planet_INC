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
#ifndef PLANET_ATMOSPHERIC_STEADY_STATE_H
#define PLANET_ATMOSPHERIC_STEADY_STATE_H

//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/kinetics_evaluator.h"
#include "antioch/cmath_shims.h"

//Planet
#include "planet/atmospheric_mixture.h"

//eigen
#include <Eigen/Dense>

//C++
#include <map>

namespace Planet
{
//name convention:
// mole_*  -> full system, neutral inside
// molar_* -> ionic system, only ions
  template <typename CoeffType, typename VectorCoeffType = std::vector<CoeffType> >
  class AtmosphericSteadyState
  {
      private:
        //don't use it
        AtmosphericSteadyState();

        //from AtmosphericKinetics
        //!contains all the system
        Antioch::KineticsEvaluator<CoeffType> &_reactions_system; 
        //!targets only the species at steady state
        const std::vector<Antioch::Species>   & _ss_species;

        //internal stuff
        std::map<unsigned int,unsigned int> _ionic_map;         //from n_species to ss_species
        std::map<unsigned int,unsigned int> _ionic_inverse_map; //to ss_species to n_species
        VectorCoeffType _rates;
        VectorCoeffType _mole_concentrations;
        VectorCoeffType _updated_rates;
        VectorCoeffType _molar_concentrations;


        //solver library, for the Ax = b solve
        Eigen::Matrix<CoeffType,Eigen::Dynamic,Eigen::Dynamic> A;
        Eigen::Matrix<CoeffType,Eigen::Dynamic,1>              b;


        //! to avoid antioch long calculations
        template <typename VectorStateType, typename MatrixStateType>
        void compute_sources_and_jacob(VectorStateType & mole_sources, MatrixStateType & dmole_dX_s) const;

        //! closes the system, ionic system with electric neutrality
        template <typename VectorSolveType, typename MatrixSolveType>
        void bring_me_closure(VectorSolveType & deriv, MatrixSolveType & jacob) const;

        //! gives a first approximation if needed
        void first_approximation();

        //! final calculations for output, only sources here
        template <typename VectorStateType>
        void compute_full_sources(VectorStateType & mole_sources) const;

        //! final calculations for output, sources and derivs here
        // not const as it updates the ionic values in the full system
        template <typename VectorStateType, typename MatrixStateType>
        void compute_full_sources_and_derivs(VectorStateType &mole_sources, MatrixStateType &dmole_dX_s);

        //! Newton solver here
        void solve();

      public:

        AtmosphericSteadyState(const std::vector<Antioch::Species> &ss_species, Antioch::KineticsEvaluator<CoeffType> &reactions_system);
        ~AtmosphericSteadyState();

        //! 
        void build_map();

        //! caches updated rate constants
        template <typename StateType, typename VectorStateType>
        void precompute_rates(const VectorStateType & mole_concentrations, 
                              const StateType & T_neutral,
                              const StateType & T_ions,
                              const StateType & T_electrons);

        //! Newton solver
        template <typename VectorStateType>
        void steady_state(VectorStateType & mole_sources);

        //! Newton solver
        template <typename VectorStateType, typename MatrixStateType>
        void steady_state_and_derivs(VectorStateType & mole_sources, MatrixStateType & drate_dn);

  };


  template <typename CoeffType, typename VectorCoeffType>
  inline
  void AtmosphericSteadyState<CoeffType,VectorCoeffType>::build_map()
  {
     const Antioch::ChemicalMixture<CoeffType> & mixture = _reactions_system.reaction_set().chemical_mixture();
     for(unsigned int s = 0; s < _ss_species.size(); s++)
     {
         _ionic_inverse_map[s] = mixture.species_list_map().at(_ss_species[s]);
         _ionic_map[mixture.species_list_map().at(_ss_species[s])] = s;
     }
  }

// this will calculate the rate constant and multiply by the neutral concentrations
// these values won't change, so no recomputing during the loop
// TODO how can we generalize the temperature specialization?
  template <typename CoeffType, typename VectorCoeffType>
  template <typename StateType, typename VectorStateType>
  inline
  void AtmosphericSteadyState<CoeffType,VectorCoeffType>::precompute_rates(const VectorStateType & mole_concentrations,
                                                                           const StateType & T_neutral,
                                                                           const StateType & T_ions,
                                                                           const StateType & T_electrons)
  {
    Antioch::set_zero(_updated_rates);
    Antioch::set_zero(_rates);

    _mole_concentrations = mole_concentrations;

    // compute the requisite reaction rates
    for(unsigned int rxn = 0; rxn < _reactions_system.n_reactions(); rxn++)
    {
       const Antioch::Reaction<CoeffType> & reaction = _reactions_system.reaction_set().reaction(rxn);
       StateType kfwd;
       Antioch::set_zero(kfwd);
       if(reaction.kinetics_model() == Antioch::KineticsModel::PHOTOCHEM)
       {
          kfwd = reaction.compute_photo_rate_of_progress(mole_concentrations, (*_reactions_system.particle_flux()[_reactions_system.reaction_particle_flux_map().at(rxn)] ),T_neutral);
       }else
       {
//TODO separate the different kind of reactions in a better way
          kfwd = (reaction.kinetics_model() == Antioch::KineticsModel::HERCOURT_ESSEN)?
                    reaction.compute_forward_rate_coefficient( mole_concentrations, T_neutral):   // any reaction that is not a DR
                    reaction.compute_forward_rate_coefficient( mole_concentrations, T_electrons); // DR
       }

       _updated_rates[rxn] = kfwd;
       _rates[rxn] = kfwd;
       // adding neutral contribution, updating
       // update with neutral concentrations
       for (unsigned int r=0; r< reaction.n_reactants(); r++)
       {
          if(!_ionic_map.count(reaction.reactant_id(r)))
          {
             _updated_rates[rxn] *= Antioch::ant_pow( mole_concentrations[reaction.reactant_id(r)],
                                              static_cast<int>(reaction.reactant_stoichiometric_coefficient(r)) );
          }
        }
    }
    return;
  }



  template <typename CoeffType, typename VectorCoeffType>
  inline
  void AtmosphericSteadyState<CoeffType, VectorCoeffType>::first_approximation()
  {

   const Antioch::ChemicalMixture<CoeffType> & mixture = _reactions_system.reaction_set().chemical_mixture();

//all ions to 0
    Antioch::set_zero(_molar_concentrations);

//first approx C_i = prod_i / (dloss_dCi) = kfwd_const * fwd_conc / (kfwd_const * conc_{no ss species}) ~ updated_rates_fwd / updated_rate_bkwd
    VectorCoeffType sum_forward;
    VectorCoeffType sum_backward;
    sum_forward.resize(_ss_species.size());
    sum_backward.resize(_ss_species.size());
    Antioch::set_zero(sum_forward);
    Antioch::set_zero(sum_backward);
    for(unsigned int rxn = 0; rxn < _updated_rates.size(); rxn++)
    {
//prod
       const Antioch::Reaction<CoeffType>& reaction = _reactions_system.reaction_set().reaction(rxn);
       for (unsigned int p=0; p<reaction.n_products(); p++)
       {
         if(!_ionic_map.count(reaction.product_id(p)))continue;
         sum_forward[_ionic_map[reaction.product_id(p)]] += _updated_rates[rxn];
       }

// loss
       for (unsigned int r=0; r<reaction.n_reactants(); r++)
       {
         if(!_ionic_map.count(reaction.reactant_id(r)))continue;
         sum_backward[_ionic_map[reaction.reactant_id(r)]] += _updated_rates[rxn];
        }
     }

     CoeffType sum;
     Antioch::set_zero(sum);
     unsigned int s_electron(_ionic_map.at(mixture.species_list_map().at(Antioch::Species::e)));
     for(unsigned int ss = 0; ss < _ss_species.size(); ss++)
     {
       if(ss == s_electron)continue;
        _molar_concentrations[ss] = Antioch::ant_sqrt(sum_forward[ss] / sum_backward[ss]);
        sum += _molar_concentrations[ss];
     }
     _molar_concentrations[s_electron] = sum;

//beurk, comparison to zero
//TODO find something better than comparing a real number to zero
     const CoeffType tol = std::numeric_limits<CoeffType>::epsilon() * 100.L;
     if(sum < tol)
     {
        std::cerr << "Error: A first approximation is global 0" << std::endl;
        antioch_error();
     }
     return;
  }


  template <typename CoeffType, typename VectorCoeffType>
  template <typename VectorSolveType, typename MatrixSolveType>
  inline
  void AtmosphericSteadyState<CoeffType,VectorCoeffType>::bring_me_closure(VectorSolveType & deriv,
                                                                           MatrixSolveType & jacob) const
  {
        const Antioch::ChemicalMixture<CoeffType> & mixture = _reactions_system.reaction_set().chemical_mixture();

// neutral atmosphere approximation: [e] = sum [ions]
        CoeffType sum;
        Antioch::set_zero(sum);
        unsigned int s_electron(_ionic_map.at(mixture.species_list_map().at(Antioch::Species::e)));
        for(unsigned int s = 0; s < _ss_species.size(); s++)
        {
           if(s == s_electron)continue;
           jacob[s_electron][s] = 1.L;
           sum += _molar_concentrations[s];
        }
        jacob[s_electron][s_electron] = -1.L;
        deriv[s_electron] = sum - _molar_concentrations[s_electron];

        return;
  }


//this will do the second par of Antioch::KineticsEvaluator<CoeffType,StateType>::compute_mole_sources_and_derivs()
  template <typename CoeffType, typename VectorCoeffType>
  template <typename VectorStateType, typename MatrixStateType>
  inline
  void AtmosphericSteadyState<CoeffType,VectorCoeffType>::compute_sources_and_jacob(VectorStateType & molar_sources,
                                                                                    MatrixStateType & dmolar_dX_s) const
  {

//initialization
    Antioch::set_zero(molar_sources);
    for(unsigned int ss = 0; ss < _ss_species.size(); ss++)
    {
       Antioch::set_zero(dmolar_dX_s[ss]);
    }

    for(unsigned int rxn = 0; rxn < _reactions_system.reaction_set().n_reactions(); rxn++)
    {

        const Antioch::Reaction<CoeffType> & reac = _reactions_system.reaction_set().reaction(rxn);

        VectorStateType dRfwd_dX_s(reac.n_species(), 0.L);
        CoeffType facfwd(1.L);
    
        // pre-fill the participating species partials with the updated rates
        for (unsigned int r=0; r< reac.n_reactants(); r++)
        {
             if(!_ionic_map.count(reac.reactant_id(r)))continue;

                dRfwd_dX_s[_ionic_map.at(reac.reactant_id(r))] = _updated_rates[rxn];
        }
           // product of concentrations and derivatives term
        for (unsigned int ro=0; ro < reac.n_reactants(); ro++)
        {
//we consider only the ions
           if(!_ionic_map.count(reac.reactant_id(ro)))continue;

           const CoeffType val = 
                Antioch::ant_pow(_molar_concentrations[_ionic_map.at(reac.reactant_id(ro))],
                                 static_cast<int>(reac.reactant_stoichiometric_coefficient(ro)) );

           facfwd *= val;

           const CoeffType dval = 
                 ( static_cast<CoeffType>(reac.reactant_stoichiometric_coefficient(ro))*
                Antioch::ant_pow(_molar_concentrations[_ionic_map.at(reac.reactant_id(ro))],
                                 static_cast<int>(reac.reactant_stoichiometric_coefficient(ro))-1 ) 
                );

           for (unsigned int ri=0; ri<reac.n_reactants(); ri++)
           {
              if(!_ionic_map.count(reac.reactant_id(ri)))continue;

              dRfwd_dX_s[_ionic_map.at(reac.reactant_id(ri))] *= (ri == ro) ? dval : val;
           }
         }

        /// calculate now sources and derivatives

        // reactants contributions
        for (unsigned int r = 0; r < reac.n_reactants(); r++)
          {
            if(!_ionic_map.count(reac.reactant_id(r)))continue;

            const unsigned int s_id = _ionic_map.at(reac.reactant_id(r));
            const unsigned int r_stoich = reac.reactant_stoichiometric_coefficient(r);
            
            // d/dX_s rate contributions, no need to consider other species than reactants
            for (unsigned int rr=0; rr < reac.n_reactants(); rr++)
              {
                if(!_ionic_map.count(reac.reactant_id(rr)))continue;
                unsigned int s = _ionic_map.at(reac.reactant_id(rr));
                dmolar_dX_s[s_id][s] -= (static_cast<CoeffType>(r_stoich) * dRfwd_dX_s[s]);
              }
            // sources, loss term
            molar_sources[_ionic_map.at(reac.reactant_id(r))] -= static_cast<CoeffType>(r_stoich) * facfwd * _updated_rates[rxn];
          }
        
        // product contributions
        for (unsigned int p=0; p < reac.n_products(); p++)
          {
            if(!_ionic_map.count(reac.product_id(p)))continue;
            const unsigned int s_id = _ionic_map.at(reac.product_id(p));
            const unsigned int p_stoich = reac.product_stoichiometric_coefficient(p);
            
            // d/dX_s rate contributions, no need to consider other species than reactants
            for (unsigned int r=0; r < reac.n_reactants(); r++)
              {
                if(!_ionic_map.count(reac.reactant_id(r)))continue;
                unsigned int s = _ionic_map.at(reac.reactant_id(r));
                dmolar_dX_s[s_id][s] += (static_cast<CoeffType>(p_stoich) * dRfwd_dX_s[s]);
              }
            // sources, prod term
            molar_sources[_ionic_map.at(reac.product_id(p))] += static_cast<CoeffType>(p_stoich) * facfwd * _updated_rates[rxn];
          }
    }
        
  }

//this will do the second par of Antioch::KineticsEvaluator<CoeffType,StateType>::compute_mole_sources_and_derivs()
  template <typename CoeffType, typename VectorCoeffType>
  template <typename VectorStateType>
  inline
  void AtmosphericSteadyState<CoeffType,VectorCoeffType>::compute_full_sources(VectorStateType & mole_sources) const
  {
    for(unsigned int rxn = 0; rxn < _reactions_system.reaction_set().n_reactions(); rxn++)
    {
        const Antioch::Reaction<CoeffType> & reac = _reactions_system.reaction_set().reaction(rxn);

        CoeffType facfwd(1.L);

        // product of concentrations, only ions are missing
        for (unsigned int r = 0; r < reac.n_reactants(); r++)
          {
            if(_ionic_map.count(reac.reactant_id(r)))
                                facfwd *= Antioch::ant_pow(_molar_concentrations[_ionic_map.at(reac.reactant_id(r))],
                                          static_cast<int>(reac.reactant_stoichiometric_coefficient(r)) );
          }


        // reactants contributions
        for (unsigned int r = 0; r < reac.n_reactants(); r++)
          {
            const unsigned int r_stoich = reac.reactant_stoichiometric_coefficient(r);
            // sources, loss term
            mole_sources[reac.reactant_id(r)] -= static_cast<CoeffType>(r_stoich) * facfwd * _updated_rates[rxn];
          }
        
        // product contributions
        for (unsigned int p=0; p < reac.n_products(); p++)
          {
            const unsigned int p_stoich = reac.product_stoichiometric_coefficient(p);
            // sources, prod term
            mole_sources[reac.product_id(p)] += static_cast<CoeffType>(p_stoich) * facfwd * _updated_rates[rxn];
          }
    }
  }

  template <typename CoeffType, typename VectorCoeffType>
  template <typename VectorStateType, typename MatrixStateType>
  inline
  void AtmosphericSteadyState<CoeffType,VectorCoeffType>::compute_full_sources_and_derivs(VectorStateType &mole_sources, MatrixStateType &dmole_dX_s)
  {
//initialization
    Antioch::set_zero(mole_sources);
    for(unsigned int ss = 0; ss < _ss_species.size(); ss++)
    {
       Antioch::set_zero(dmole_dX_s[ss]);
       _mole_concentrations[_ionic_inverse_map.at(ss)] = _molar_concentrations[ss];
    }

    for(unsigned int rxn = 0; rxn < _reactions_system.reaction_set().n_reactions(); rxn++)
    {

        const Antioch::Reaction<CoeffType> & reac = _reactions_system.reaction_set().reaction(rxn);

        VectorStateType dRfwd_dX_s(reac.n_species(), 0.L);
        CoeffType facfwd(1.L);
    
        // pre-fill the participating species partials with the updated rates
        for (unsigned int r=0; r< reac.n_reactants(); r++)
        {
           dRfwd_dX_s[reac.reactant_id(r)] = _rates[rxn];
        }
           // product of concentrations and derivatives term
        for (unsigned int ro=0; ro < reac.n_reactants(); ro++)
        {
           const CoeffType val = 
                Antioch::ant_pow(_mole_concentrations[reac.reactant_id(ro)],
                                 static_cast<int>(reac.reactant_stoichiometric_coefficient(ro)) );

           facfwd *= val;

           const CoeffType dval = 
                 ( static_cast<CoeffType>(reac.reactant_stoichiometric_coefficient(ro))*
                Antioch::ant_pow(_mole_concentrations[reac.reactant_id(ro)],
                                 static_cast<int>(reac.reactant_stoichiometric_coefficient(ro))-1 ) 
                );

           for (unsigned int ri=0; ri<reac.n_reactants(); ri++)
           {
              dRfwd_dX_s[reac.reactant_id(ri)] *= (ri == ro) ? dval : val;
           }
         }

        /// calculate now sources and derivatives

        // reactants contributions
        for (unsigned int r = 0; r < reac.n_reactants(); r++)
          {
            const unsigned int s_id = reac.reactant_id(r);
            const unsigned int r_stoich = reac.reactant_stoichiometric_coefficient(r);
            
            // d/dX_s rate contributions, no need to consider other species than reactants
            for (unsigned int rr=0; rr < reac.n_reactants(); rr++)
              {
                unsigned int s = reac.reactant_id(rr);
                dmole_dX_s[s_id][s] -= (static_cast<CoeffType>(r_stoich) * dRfwd_dX_s[s]);
              }
            // sources, loss term
            mole_sources[reac.reactant_id(r)] -= static_cast<CoeffType>(r_stoich) * facfwd * _rates[rxn];
          }
        
        // product contributions
        for (unsigned int p=0; p < reac.n_products(); p++)
          {
            const unsigned int s_id = reac.product_id(p);
            const unsigned int p_stoich = reac.product_stoichiometric_coefficient(p);
            
            // d/dX_s rate contributions, no need to consider other species than reactants
            for (unsigned int r=0; r < reac.n_reactants(); r++)
              {
                unsigned int s = reac.reactant_id(r);
                dmole_dX_s[s_id][s] += (static_cast<CoeffType>(p_stoich) * dRfwd_dX_s[s]);
              }
            // sources, prod term
            mole_sources[reac.product_id(p)] += static_cast<CoeffType>(p_stoich) * facfwd * _rates[rxn];
          }
    }
  }

  template <typename CoeffType, typename VectorCoeffType>
  inline
  void AtmosphericSteadyState<CoeffType,VectorCoeffType>::solve()
  {

   if(_molar_concentrations[0] < 0.L)first_approximation();

// Newton solver here
// Ax + b = 0
// A is jacobian, b is molar sources

//setup, ionic system syzed
    VectorCoeffType molar_sources;
    std::vector<VectorCoeffType> dmolar_dX_s;

    molar_sources.resize(_ss_species.size());
    dmolar_dX_s.resize(_ss_species.size());
    for(unsigned int s=0; s < _ss_species.size();s++)
    {
      dmolar_dX_s[s].resize(_ss_species.size(),0.L);
    }

// shoot
    CoeffType lim(1.L);
    CoeffType res_mol(0.L);

// physically this precision is ridiculous, which is nice
    CoeffType thresh = std::numeric_limits<CoeffType>::epsilon() * 200.;
//    if(thresh < 1e-10)thresh = 1e-10; 

    unsigned int loop_max(100);
    unsigned int nloop(0);

    while(lim > thresh)
    {

    std::cout << "loop #" << nloop << " ";
    for(unsigned int s = 0; s < _molar_concentrations.size(); s++)
    {
        std::cout << std::scientific << std::setprecision(15)
                  << _molar_concentrations[s] << " ";
    }
    std::cout << std::endl;

      this->compute_sources_and_jacob(molar_sources,dmolar_dX_s);

//TODO what are the different options here?
// neutral hypothesis for the moment   
      this->bring_me_closure(molar_sources,dmolar_dX_s);

      Antioch::set_zero(res_mol);
      for(unsigned int i = 0; i < _ss_species.size(); i++)
      {
        for(unsigned int j = 0; j < _ss_species.size(); j++)
        {
           A(i,j) = dmolar_dX_s[i][j]; //Jacobian
        }
        b(i) = - molar_sources[i]; // - first derivative
        res_mol += Antioch::ant_abs(molar_sources[i]);
      }
      if(res_mol < thresh)break;

      Eigen::PartialPivLU<Eigen::Matrix<CoeffType,Eigen::Dynamic,Eigen::Dynamic> > mypartialPivLu(A);
      Eigen::Matrix<CoeffType,Eigen::Dynamic,1> x(_ss_species.size());
      x = mypartialPivLu.solve(b);


      Antioch::set_zero(lim);
      for(unsigned int s = 0; s < _ss_species.size(); s++)
      {
        _molar_concentrations[s]  += x(s);
        if(_molar_concentrations[s] < 0.)Antioch::set_zero(_molar_concentrations[s]);
        lim += Antioch::ant_abs(x(s)); //absolute increment
      }
      std::cout << std::endl;

      nloop++;
      if(nloop > loop_max)
      {
        std::cerr << "Newton solver failed after " << loop_max << " loops with a residual of " << lim
                  << ", a molar source total of " << res_mol
                  << " and a tolerance of " << thresh << std::endl;
        antioch_error();
      }

    } //solver loop

    std::cout << "solved with densities:\n";
    for(unsigned int s = 0; s < _molar_concentrations.size(); s++)
    {
        std::cout << _molar_concentrations[s] << " ";
    }
    std::cout << std::endl;

  }


  template <typename CoeffType, typename VectorCoeffType>
  template <typename VectorStateType>
  inline
  void AtmosphericSteadyState<CoeffType,VectorCoeffType>::steady_state(VectorStateType & mole_sources) //full system size
  {
    this->solve();
    this->compute_full_sources(mole_sources);

  }

  template <typename CoeffType, typename VectorCoeffType>
  template <typename VectorStateType, typename MatrixStateType>
  inline
  void AtmosphericSteadyState<CoeffType,VectorCoeffType>::steady_state_and_derivs(VectorStateType & mole_sources, MatrixStateType & drate_dn)
  {
    this->solve();
    this->compute_full_sources_and_derivs(mole_sources,drate_dn);
  }

  
  template <typename CoeffType, typename VectorCoeffType>
  inline
  AtmosphericSteadyState<CoeffType,VectorCoeffType>::AtmosphericSteadyState(const std::vector<Antioch::Species> & ss_species, Antioch::KineticsEvaluator<CoeffType> &reactions_system):
      _reactions_system(reactions_system),
      _ss_species(ss_species),
      A(ss_species.size(),ss_species.size()),
      b(ss_species.size())
  {
     _updated_rates.resize(_reactions_system.n_reactions(),0.L);
     _molar_concentrations.resize(_ss_species.size(),-1.L);
     _rates.resize(_reactions_system.n_reactions(),0.L);
     _mole_concentrations.resize(reactions_system.n_species(),-1.L);
     this->build_map();
     return;
  }

  template <typename CoeffType, typename VectorCoeffType>
  inline
  AtmosphericSteadyState<CoeffType,VectorCoeffType>::~AtmosphericSteadyState()
  {
     return;
  }

}

#endif
