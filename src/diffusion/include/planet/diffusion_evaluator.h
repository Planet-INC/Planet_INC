//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef _PLANET_DIFFUSION_EVALUATOR_
#define _PLANET_DIFFUSION_EVALUATOR_

//Antioch
#include "antioch/metaprogramming_decl.h"
#include "antioch/cmath_shims.h"
#include "antioch/chemical_mixture.h"

//Planet
#include "planet/binary_diffusion.h"

//C++


namespace Planet
{
  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  class DiffusionEvaluator
  {
     private:
        //!
        DiffusionEvaluator();

        MatrixCoeffType _Dtilde;
        std::vector<std::vector<BinaryDiffusion<CoeffType> > > _diffusion;
        Antioch::ChemicalMixture<CoeffType> &_composition;
        const unsigned int _n_medium;


     public:
        //!
        DiffusionEvaluator(Antioch::ChemicalMixture<CoeffType> &comp);
        //!
        ~DiffusionEvaluator();

        //! update Dtilde
        template<typename StateType, typename VectorStateType, typename MatrixStateType>
        void make_diffusion(const Atmosphere<StateType,VectorStateType,MatrixStateType> &atm);

        //!
        CoeffType diffusion_coefficient(unsigned int s, unsigned int iz) const;

        //!
        VectorCoeffType species_top_to_bottom(unsigned int ineu) const;

        //!
        VectorCoeffType species_bottom_to_top(unsigned int ineu) const;

        //!
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        binary_coefficient(unsigned int i, unsigned int j, const StateType &T, const StateType &P) const
        ANTIOCH_AUTOFUNC(StateType,(_diffusion[i][j].diffusion_model() != DiffusionType::NoData)?
                                             this->binary_coefficient_known(i,j,T,P):
                                                (_composition.M(j) < _composition.M(i))?
                                                        this->binary_coefficient_unknown_ji(i,j,T,P):
                                                        this->binary_coefficient_unknown_ij(i,j,T,P))

        //!
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        binary_coefficient_known(unsigned int i, unsigned int j, const StateType &T, const StateType &P) const
        ANTIOCH_AUTOFUNC(StateType,_diffusion[i][j].binary_coefficient(T,P))

        //!
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        binary_coefficient_unknown_ji(unsigned int i, unsigned int j, const StateType &T, const StateType &P) const
        ANTIOCH_AUTOFUNC(StateType,_diffusion[i][i].binary_coefficient(T,P) * 
                                   Antioch::ant_pow( (_composition.M(j)/_composition.M(i) + StateType(1.L)) / StateType(2.L) )
                                   )
        //!
        template<typename StateType>
        ANTIOCH_AUTO(StateType)
        binary_coefficient_unknown_ij(unsigned int i, unsigned int j, const StateType &T, const StateType &P) const
        ANTIOCH_AUTOFUNC(StateType,_diffusion[i][i].binary_coefficient(T,P) * Antioch::ant_pow(_composition.M(j)/_composition.M(i)))

        //!
        void set_chemical_mixture(Antioch::ChemicalMixture<CoeffType> &comp);

        //!
        void set_binary_coefficient(unsigned int i, unsigned int j, const BinaryDiffusion<CoeffType> &bin_coef);

  };


  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  CoeffType DiffusionEvaluator<CoeffType,MatrixCoeffType>::diffusion_coefficient(unsigned int s, unsigned int iz) const
  {
      return _Dtilde.at(iz).at(s);
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  VectorStateType DiffusionEvaluator<CoeffType,MatrixCoeffType>::species_top_to_bottom(unsigned int ineu) const
  {
     return _Dtilde[ineu];
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  VectorCoeffType DiffusionEvaluator<CoeffType,MatrixCoeffType>::species_bottom_to_top(unsigned int ineu) const
  {
     VectorCoeffType out;
     out.resize(_Dtilde[ineu].size(),0.L);
     for(unsigned int iz = 0; iz < _Dtilde[ineu].size(); iz++)
     {
        out[iz] = _Dtilde[ineu][_Dtilde[ineu].size() - 1 - iz];
     }

     return out;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType, typename MatrixStateType>
  inline
  void DiffusionEvaluator<CoeffType,MatrixCoeffType>::make_diffusion(const Atmosphere<StateType,VectorStateType,MatrixStateType> &atm)
  {
    _Dtilde.resize(_composition.n_species()); //from top to bottom
    std::vector<CoeffType> meanM;

    for(unsigned int s = 0; s < _composition.n_species(); s++)
    {
      meanM.clear();
      meanM.resize(atm.n_altitudes(),0.L);
      for(StateType z = atm.min_alt(); z <= atm.min_alt(); z += atm.step_alt())
      {
        unsigned int iz = atm.altitude(z);
        for(unsigned int i = 0; i < _composition.n_species(); i++)
        {
          if(i == s)continue;
          meanM[iz] += _composition.M(i) * atm.neutral_molar_fraction(i,z);
        }
        meanM[iz] /= StateType(_composition.n_species() - 1);
      }
       _Dtilde[s].resize(atm.n_altitudes(),0.L);
       for(StateType z = atm.min_alt(); z <= atm.min_alt(); z += atm.step_alt())
       {
          unsigned int iz = atm.altitude(z);
          CoeffType n_D;
          Antioch::set_zero(n_D);
          for(unsigned int i = 0; i < _n_medium; i++)
          {
            if(i == s)continue;
            n_D += atm.neutral_molar_density(i,z) / this->binary_coefficient(i,s,atm.temperature(z),atm.pressure(z));
          }
          CoeffType Ds = (atm.total_density(z) - atm.neutral_molar_density(s,z))/n_D;
          _Dtilde[s][iz] = Ds / (CoeffType(1.L) - atm.neutral_molar_density(s,z) / atm.total_density(z) * 
                                  (CoeffType(1.L) - _composition.M(s) / meanM[iz])
                                );
       }
    }
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  DiffusionEvaluator<CoeffType,MatrixCoeffType>::DiffusionEvaluator()
  {
     antioch_error();
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  DiffusionEvaluator<CoeffType,MatrixCoeffType>::~DiffusionEvaluator()
  {
     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  DiffusionEvaluator<CoeffType,MatrixCoeffType>::DiffusionEvaluator(Antioch::ChemicalMixture<CoeffType> &comp):
       _composition(comp),
       _n_medium(2)
  {
     _diffusion.resize(_n_medium); //N2,CH4
     for(unsigned int i = 0; i < _diffusion.size(); i++)
     {
        _diffusion[i].resize(_composition.n_species());
     }
     return; 
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  void DiffusionEvaluator<CoeffType,MatrixCoeffType>::set_binary_coefficient(unsigned int i, unsigned int j, const BinaryDiffusion<CoeffType> &bin_coef)
  {
      _diffusion[i][j] = bin_coef;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  void DiffusionEvaluator<CoeffType,MatrixCoeffType>::set_chemical_mixture(Antioch::ChemicalMixture<CoeffType> &comp)
  {
     _composition = comp;
     return;
  }

}

#endif
