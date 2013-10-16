//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef PLANET_DIFFUSION_EVALUATOR_H
#define PLANET_DIFFUSION_EVALUATOR_H

//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/chemical_mixture.h"

//Planet
#include "planet/molecular_diffusion_evaluator.h"
#include "planet/eddy_diffusion_evaluator.h"

//C++
#include <string>

namespace Planet{

  /*!\class DiffusionEvaluator
 * Stores all kind of diffusions
 *
 */
  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  class DiffusionEvaluator
  {
      private:
       DiffusionEvaluator();

       MolecularDiffusionEvaluator<CoeffType,MatrixCoeffType> _molecular_diff;
       EddyDiffusionEvaluator<CoeffType,VectorCoeffType> _eddy_diff;

      public:
       DiffusionEvaluator(const Antioch::ChemicalMixture<CoeffType> &comp, const CoeffType K0, const VectorCoeffType &ntot);
       ~DiffusionEvaluator();

//eddy coeff
// - make
// - get
// - reset

       //!compute the eddy diffusion
       void make_eddy_diffusion();

       //!returns the eddy diffusion factor K
       const VectorCoeffType eddy_diffusion_coeff() const;

       //!reset the eddy coefficient parameters
       template<typename StateType, typename VectorStateType>
       void reset_eddy_coefficient(const StateType &K0, const VectorStateType &ntot);

//molecular
// - set
// - make
// - get
// - reset
       //!sets a binary molecular diffusion coefficient
       template<typename StateType>
       void set_molecular_diffusion_coeff(unsigned int i, unsigned int j, const BinaryDiffusion<StateType> &bin_coef);

       //!evaluate the molecular diffusion for a given atmosphere
       template <typename StateType, typename VectorStateType, typename MatrixStateType>
       void make_molecular_diffusion(const Atmosphere<StateType,VectorStateType,MatrixStateType> &atm);

       //!returns the molecular diffusion matrix
       const MatrixCoeffType molecular_diffusion() const;

       //!reset the chemical mixture
       template <typename StateType>
       void reset_molecular_diffusion(Antioch::ChemicalMixture<StateType> &comp);


  };

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  DiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::DiffusionEvaluator()
  {
     antioch_error();
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  DiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::~DiffusionEvaluator()
  {
     return;
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  DiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::DiffusionEvaluator(const Antioch::ChemicalMixture<CoeffType> &comp, 
                                                                                      const CoeffType K0, const VectorCoeffType &ntot)
    _molecular_diff(comp),
    _eddy_diff(K0,ntot)
  {
     return;
  }

  
  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  void DiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::make_eddy_diffusion()
  {
      _eddy_diff.make_eddy_diffusion();
      return;
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const VectorCoeffType DiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::eddy_diffusion_coeff()const
  {
      return _eddy_diff.K();
  }
  
  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType>
  inline
  void DiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::reset_eddy_coefficient(const StateType &K0, const VectorStateType &ntot)
  {
     _eddy_diff.reset_K0(K0);
     _eddy_diff.reset_ntot(ntot);
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType>
  inline
  void DiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::set_molecular_diffusion_coeff(unsigned int i, 
                                                                                                      unsigned int j, 
                                                                                                      const BinaryDiffusion<StateType> &bin_coef)
  {
     _molecular_diff.set_binary_coefficient(i,j,bin_coef);
     return;
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType, typename VectorStateType, typename MatrixStateType>
  inline
  void DiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::make_molecular_diffusion(const Atmosphere<StateType,VectorStateType,MatrixStateType> &atm)
  {
     _molecular_diff.make_diffusion(atm);
     return;
  }

  template <typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  const DiffusionEvaluator<CoeffType, VectorCoeffType, MatrixCoeffType>::MatrixCoeffType molecular_diffusion() const
  {
     return _molecular_diff.Dtilde();
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  template<typename StateType>
  inline
  void MolecularDiffusionEvaluator<CoeffType, MatrixCoeffType>::reset_molecular_diffusion(Antioch::ChemicalMixture<StateType> &comp)
  {
     _molecular_diff.reset_chemical_mixture(comp);
     return;
  }
  
}

#endif
