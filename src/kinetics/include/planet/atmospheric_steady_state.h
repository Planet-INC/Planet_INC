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

namespace Planet
{
  class AtmosphericSteadyState
  {
      private:
      public:
        AtmosphericSteadyState();
        ~AtmosphericSteadyState();

        //! Newton solver
        template <typename StateType, typename VectorStateType>
        void steady_state(Antioch::KineticsEvaluator<StateType> &reactions_system,
                          const std::vector<Antioch::Species> & ss_species,
                          const Antioch::ChemicalMixture<StateType> &mixture,
                          const StateType &T,
                          VectorStateType & molar_concentrations,
                          VectorStateType & molar_sources) const;

        //! closes the system, ionic system with electric neutrality
        template <typename StateType, typename VectorStateType, typename VectorSolve, typename MatrixSolve>
        void bring_me_closure(MatrixSolve &A,
                              VectorSolve &b,
                              const std::vector<Antioch::Species> & ss_species,
                              const Antioch::ChemicalMixture<StateType> &mixture,
                              const VectorStateType &molar_concentrations) const;

        //! gives a first approximation if needed
        template <typename StateType, typename VectorStateType>
        void first_approximation(const Antioch::ReactionSet<StateType>     &reactions_set, 
                                 const std::vector<Antioch::Species>       & ss_species, 
                                 const Antioch::ChemicalMixture<StateType> & mixture, 
                                 const StateType                           & T, 
                                 VectorStateType                           & molar_concentrations) const;
  };


  template <typename StateType, typename VectorStateType>
  inline
  void AtmosphericSteadyState::first_approximation(const Antioch::ReactionSet<StateType>     &reactions_set, 
                                                   const std::vector<Antioch::Species>       & ss_species, 
                                                   const Antioch::ChemicalMixture<StateType> & mixture, 
                                                   const StateType                           & T, 
                                                   VectorStateType                           & molar_concentrations) const
  {
//all ions to 0
   for(unsigned int s = 0; s < ss_species.size(); s++)
   {
        Antioch::set_zero(molar_concentrations[mixture.species_list_map().at(ss_species[s])]);
   }

//setup
    VectorStateType h_RT_minus_s_R;
    h_RT_minus_s_R.resize(reactions_set.n_species(),0.L); //irreversible

    VectorStateType net_rates;
    VectorStateType kfwd_const;
    VectorStateType kbkwd_const;
    VectorStateType kfwd;
    VectorStateType kbkwd;
    VectorStateType fwd_conc;
    VectorStateType bkwd_conc;
    net_rates.resize(reactions_set.n_reactions());
    kfwd_const.resize(reactions_set.n_reactions());
    kbkwd_const.resize(reactions_set.n_reactions());
    kfwd.resize(reactions_set.n_reactions());
    kbkwd.resize(reactions_set.n_reactions());
    fwd_conc.resize(reactions_set.n_reactions());
    bkwd_conc.resize(reactions_set.n_reactions());

    reactions_set.get_reactive_scheme(T, molar_concentrations,h_RT_minus_s_R,
                                      net_rates, kfwd_const, kbkwd_const,
                                      kfwd, kbkwd, fwd_conc, bkwd_conc);

   const StateType tol = std::numeric_limits<StateType>::epsilon() * 100.L;
//first approx C_i = prod_i / (dloss_dCi) = kfwd_const * fwd_conc / (kfwd_const * conc_{no ss species})
// not a good algorithm...
    StateType sum;
    Antioch::set_zero(sum);
    for(unsigned int s = 0; s < ss_species.size(); s++)
    {
        if(ss_species[s] == Antioch::Species::e)continue;
//prod
        StateType prod(0.L);
        for(unsigned int rxn = 0; rxn < reactions_set.n_reactions(); rxn++)
        {
           const Antioch::Reaction<StateType>& reaction = reactions_set.reaction(rxn);
           for (unsigned int p=0; p<reaction.n_products(); p++)
           {
              if(reaction.product_id(p) == mixture.species_list_map().at(ss_species[s])) //if s is a product
              {
                 prod += kfwd_const[rxn] * fwd_conc[rxn];
                 break; 
              }
           }
        }

// dloss_dCi
        StateType loss(0.L);
        for(unsigned int rxn = 0; rxn < reactions_set.n_reactions(); rxn++)
        {
           const Antioch::Reaction<StateType>& reaction = reactions_set.reaction(rxn);
           StateType conc(1.L);
           bool add(false);
           for (unsigned int r=0; r<reaction.n_reactants(); r++)
           {
              if(reaction.reactant_id(r) == mixture.species_list_map().at(ss_species[s]))//if s is a reactant
              {
                  add = true;
                  continue; 
              }
              conc *= pow( molar_concentrations[reaction.reactant_id(r)],
                              static_cast<int>(reaction.reactant_stoichiometric_coefficient(r)) );
           }
           if(add)loss += conc * kfwd_const[rxn];
        }
        if(loss < tol)loss = 1.L;
        molar_concentrations[mixture.species_list_map().at(ss_species[s])] = prod / loss;// < 1e-4)?1e-4:prod/loss;
        if(molar_concentrations[mixture.species_list_map().at(ss_species[s])] < tol)
            molar_concentrations[mixture.species_list_map().at(ss_species[s])] = tol;
                
        sum += molar_concentrations[mixture.species_list_map().at(ss_species[s])];
    }
    molar_concentrations[mixture.species_list_map().at(Antioch::Species::e)] = sum; //e
    if(sum < tol)
    {
        std::cerr << "Error: A first approximation is global 0" << std::endl;
        antioch_error();
    }
  }


  template <typename StateType, typename VectorStateType, typename VectorSolve, typename MatrixSolve>
  inline
  void AtmosphericSteadyState::bring_me_closure(MatrixSolve &A,
                                                VectorSolve &b,
                                                const std::vector<Antioch::Species> & ss_species,
                                                const Antioch::ChemicalMixture<StateType> &mixture,
                                                const VectorStateType &molar_concentrations) const
  {
        StateType sum;
        Antioch::set_zero(sum);
        unsigned int i_electron = mixture.species_list_map().at(Antioch::Species::e);
        for(unsigned int s = 0; s < ss_species.size(); s++)
        {
           A[i_electron][s] = -1.L;
           sum += molar_concentrations[mixture.species_list_map().at(ss_species[s])];
        }
        A[i_electron][i_electron] = 1.L;
        b[i_electron] = 2.L * molar_concentrations[i_electron] - sum;
  }


  template <typename StateType, typename VectorStateType>
  inline
  void AtmosphericSteadyState::steady_state(Antioch::KineticsEvaluator<StateType> &reactions_system,
                                            const std::vector<Antioch::Species> & ss_species,
                                            const Antioch::ChemicalMixture<StateType> &mixture,
                                            const StateType &T,
                                            VectorStateType & molar_concentrations,
                                            VectorStateType & molar_sources) const
  {

   if(molar_concentrations[mixture.species_list_map().at(ss_species[0])] < 0.L)
        first_approximation(reactions_system.reaction_set(), ss_species, mixture, T, molar_concentrations);

    for(unsigned int s = 0; s < ss_species.size(); s++)
    {
std::cout << mixture.species_inverse_name_map().at(ss_species[s]) << " " 
          << std::setprecision(15) << molar_concentrations[mixture.species_list_map().at(ss_species[s])] << std::endl;
    }
for(unsigned int i = 0; i < 20; i++)std::cout << "*";
std::cout << std::endl;

// Newton solver here
// Ax + b = 0
// A is jacobian, b is what goes to 0 (dc/dt here)
    Eigen::Matrix<StateType,Eigen::Dynamic,Eigen::Dynamic> A(ss_species.size(),ss_species.size());
    Eigen::Matrix<StateType,Eigen::Dynamic,1> b(ss_species.size());

//setup
    VectorStateType h_RT_minus_s_R;
    VectorStateType dh_RT_minus_s_R_dT;
    VectorStateType dmole_dT;
    std::vector<VectorStateType> dmole_dX_s;
    h_RT_minus_s_R.resize(reactions_system.n_species(),0.L); //irreversible
    dh_RT_minus_s_R_dT.resize(reactions_system.n_species());
    molar_sources.resize(reactions_system.n_species());
    dmole_dT.resize(reactions_system.n_species());
    dmole_dX_s.resize(reactions_system.n_species());
    for(unsigned int s=0; s < reactions_system.n_species();s++)
    {
      dmole_dX_s[s].resize(reactions_system.n_species(),0.L);
    }

// shoot
    StateType lim(1.L);
    StateType thresh = std::numeric_limits<StateType>::epsilon() * 500.;
//    if(thresh < 1e-10)thresh = 1e-10; // physically this precision is ridiculous, which is nice
    unsigned int loop_max(50);
    unsigned int nloop(0);
    while(lim > thresh)
    {
      if(nloop > loop_max)antioch_error();
      reactions_system.compute_mole_sources_and_derivs(T, molar_concentrations,
                                                       h_RT_minus_s_R, dh_RT_minus_s_R_dT,
                                                       molar_sources, dmole_dT, dmole_dX_s );

      for(unsigned int s = 0; s < ss_species.size(); s++)
      {
std::cout << mixture.species_inverse_name_map().at(ss_species[s]) << " " 
          << std::setprecision(15) << molar_sources[mixture.species_list_map().at(ss_species[s])] << std::endl;
      }

//TODO what are the different options here?
// neutral hypothesis for the moment   
      this->bring_me_closure(dmole_dX_s,molar_sources,ss_species,mixture,molar_concentrations);

      for(unsigned int i = 0; i < ss_species.size(); i++)
      {
        unsigned int i_ss = mixture.species_list_map().at(ss_species[i]);
        for(unsigned int j = 0; j < ss_species.size(); j++)
        {
           unsigned int j_ss = mixture.species_list_map().at(ss_species[j]);
           A(i,j) = dmole_dX_s[i_ss][j_ss];
        }
        b(i) = - molar_sources[i_ss];
      }

      Eigen::PartialPivLU<Eigen::Matrix<StateType,Eigen::Dynamic,Eigen::Dynamic> > mypartialPivLu(A);
      Eigen::Matrix<StateType,Eigen::Dynamic,1> x(ss_species.size());
      x = mypartialPivLu.solve(b);

      Antioch::set_zero(lim);
      unsigned int i_electron = mixture.species_list_map().at(Antioch::Species::e);
      Antioch::set_zero(molar_concentrations[i_electron]);
      for(unsigned int s = 0; s < ss_species.size(); s++)
      {
        molar_concentrations[mixture.species_list_map().at(ss_species[s])] += x(s);
//std::cout << mixture.species_inverse_name_map().at(ss_species[s]) << " " << std::setprecision(15) << x(s) << std::endl;
        lim += Antioch::ant_abs(x(s));
        if(mixture.species_list_map().at(ss_species[s]) == i_electron)continue;
        molar_concentrations[i_electron] += molar_concentrations[mixture.species_list_map().at(ss_species[s])];
      }
      nloop++;
    }
std::cout << "out after " << nloop << " loops for a limit of " << lim << std::endl;
  }

  
  inline
  AtmosphericSteadyState::AtmosphericSteadyState()
  {
     return;
  }

  inline
  AtmosphericSteadyState::~AtmosphericSteadyState()
  {
     return;
  }

}

#endif
