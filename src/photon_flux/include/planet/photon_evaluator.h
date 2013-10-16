//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef PLANET_PHOTON_EVALUATOR_H
#define PLANET_PHOTON_EVALUATOR_H

//Antioch
#include "antioch/species_enum.h"
#include "antioch/particle_flux.h"

//Planet
#include "planet/photon_flux.h"
#include "planet/cross_section.h"

//C++
#include <vector>

namespace Planet
{
  template<typename CoeffType, typename VectorCoeffType>
  class PhotonEvaluator
  {
     private:
        //! no default constructor authorized
        PhotonEvaluator(){antioch_error();return;}

        PhotonFlux<CoeffType,VectorCoeffType> &_hv_flux;
        std::map<Antioch::Species, CrossSection<VectorCoeffType> > _absorbing_species_cs;
        std::vector<Antioch::Species> _absorbing_species;

     public:
        PhotonEvaluator(PhotonFlux<CoeffType,VectorCoeffType> &hv_flux);
        PhotonEvaluator(const PhotonEvaluator<CoeffType,VectorCoeffType> &rhs);
        ~PhotonEvaluator();

        //!\return photon flux
        const PhotonFlux<CoeffType,VectorCoeffType> hv_flux() const;

        //!\return absorbing species
        const std::vector<Antioch::Species> absorbing_species() const;

        //!\return absorbing species cross-section map
        const std::map<Antioch::Species, Antioch::ParticleFlux<VectorCoeffType> > &absorbing_species_cs() const;

        //! cross-section on flux grid
        template<typename VectorStateType>
        void update_cross_section(const VectorStateType &custom_grid);

  };

  template<typename CoeffType, typename VectorCoeffType>
  inline
  PhotonEvaluator(const PhotonEvaluator<CoeffType,VectorCoeffType> &rhs):
   _hv_flux(rhs.hv_flux()),
   _absorbing_species_cs(rhs.absorbing_species_cs()),
   _absorbing_species(rhs.absorbing_species())
  {
     return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  PhotonEvaluator<CoeffType,VectorCoeffType>::PhotonEvaluator(PhotonFlux<CoeffType,VectorCoefftype> &hv_flux):
  _hv_flux(hv_flux)
  {
     _hv_flux.set_evaluator(this);
     return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  const PhotonFlux<CoeffType,VectorCoeffType> PhotonEvaluator<CoeffType,VectorCoeffType>::hv_flux() const
  {
     return _hv_flux;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  const std::vector<Antioch::Species> PhotonEvaluator<CoeffType,VectorCoeffType>::absorbing_species() const
  {
     return _absorbing_species;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  const std::map<Antioch::Species, Antioch::ParticleFlux<VectorCoeffType> > &PhotonEvaluator<CoeffType,VectorCoeffType>::absorbing_species_cs() const
  {
      return _absorbing_species_cs;
  }

  template<typename CoeffType, typename VectorCoeffType>
  template<typename VectorStateType>
  inline
  void PhotonEvaluator<CoeffType,VectorCoeffType>::update_cross_section(const VectorStateType &custom_grid)
  {
     for(unsigned int i = 0; i < _absorbing_species.size(); i++)
     {
        _absorbing_species_cs[_absorbing_species[i]].update_cross_section(custom_grid);
     }

     return;
  }

}

#endif
