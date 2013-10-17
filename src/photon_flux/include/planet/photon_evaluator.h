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
#include <map>

namespace Planet
{
  template<typename CoeffType, typename VectorCoeffType>
  class PhotonEvaluator
  {
     private:
        //! no default constructor authorized
        PhotonEvaluator(){antioch_error();return;}

//parameters & output
        Antioch::ParticleFlux<VectorCoeffType> _phy_at_top;
        std::vector<Antioch::ParticleFlux<VectorCoeffType> > _phy; //alt

//store
        std::map<Antioch::Species, unsigned int> _cross_sections_map;
        std::vector<CrossSection<VectorCoeffType> > _absorbing_species_cs;
        std::vector<Antioch::Species> _absorbing_species;

//dependencies
        PhotonOpacity<CoeffType,VectorCoeffType> &_hv_tau;
        Altitude<CoeffType,VectorCoeffType> &_alt;

        //! cross-section on flux grid
        template<typename VectorStateType>
        void update_cross_section(const VectorStateType &custom_grid);

     public:
        PhotonEvaluator(PhotonOpacity<CoeffType,VectorCoeffType> &hv_tau, Altitude<CoeffType,VectorCoeffType> &alt);
        ~PhotonEvaluator();

        //!sets the photon flux at the top of the atmosphere
        template<typename StateType, typename VectorStateType>
        void set_photon_flux_at_top(const VectorStateType &lambda, const VectorStateType &hv, const StateType &d = 1.L);

        //!\return photon flux
        const std::vector<Antioch::ParticleFlux<VectorCoeffType> > photon_flux() const;

        //!calculate photon flux
        void update_photon_flux();

        //!\return absorbing species
        const std::vector<Antioch::Species> absorbing_species() const;

        //!\return absorbing species cross-section map
        const std::map<Antioch::Species, Antioch::ParticleFlux<VectorCoeffType> > &absorbing_species_cs() const;

        //!adds a photon cross-section
        template<typename VectorStateType>
        void add_cross_section(const VectorStateType &lambda, const VectorStateType &cs, const Antioch::Species &sp);

  };

  template<typename CoeffType, typename VectorCoeffType>
  inline
  PhotonEvaluator<CoeffType,VectorCoeffType>::PhotonEvaluator(PhotonOpacity<CoeffType,VectorCoefftype> &hv_tau, Altitude<CoeffType,VectorCoeffType> &alt):
  _hv_tau(hv_tau),
  _alt(alt)
  {
     return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  template<typename VectorStateType>
  inline
  void PhotonEvaluator<CoeffType,VectorCoeffType>::set_photon_flux_at_top(const VectorStateType &lambda, const VectorStateType &hv, const StateType &d)
  {
     antioch_assert_equal_to(lambda.size(),hv.size());

//phy at top
     _phy_at_top.set_abscissa(lambda);
     VectorCoeffType flux;
     for(unsigned int i = 0; i < hv.size(); i++)
     {
        flux.push_back(hv[i]/(d * d));
     }
     _phy_at_top.set_flux(flux);

//phy at all altitudes
     _phy.clear();
     _phy.resize(_altitude.altitudes().size());
     for(unsigned int iz = 0; iz < _phy.size(); iz++)
     {
        _phy[iz].set_abscissa(lambda);
     }
//cross sections
     this->update_cross_section(lambda)

     return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  template<typename VectorStateType>
  inline
  void PhotonEvaluator<CoeffType,VectorCoeffType>::add_cross_section(const VectorStateType &lambda, const VectorStateType &cs, const Antioch::Species &sp)
  {
     _absorbing_species.push_back(sp);
     _absorbing_species_cs.push_back(CrossSection<VectorCoeffType>(lambda,cs))
     _cross_sections_map[sp] = _absorbing_species.size() - 1;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  const PhotonOpacity<CoeffType,VectorCoeffType> PhotonEvaluator<CoeffType,VectorCoeffType>::photon_flux() const
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
  const std::vector<CrossSection<VectorCoeffType> > &PhotonEvaluator<CoeffType,VectorCoeffType>::absorbing_species_cs() const
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
        _absorbing_species_cs[i].update_cross_section(custom_grid);
     }

     return;
  }

  template<typename CoeffType, typename MatrixCoeffType>
  inline
  void PhotonOpacity<CoeffType,MatrixCoeffType>::update_photon_flux(const CoeffType &a)
  {
     antioch_assert(!_phy_at_top.empty());

     MatrixCoeffType sum_densities;
     sum_densities.resize(_eval->absorbing_species().size(),0.L);

     MatrixCoeffType sigmas;
     sigma.resize(_absorbing_species.size());
     for(unsigned int s = 0; s < sigma.size(); s++)
     {
        sigma[s] =  _absorbing_species_cs[s].cross_section_on_custom_grid();
     }

//
     for(unsigned int s = 0; s < _absorbing_species.size(); s++)
     {
//from top to bottom
       unsigned int iz = _altitude.altitude_map()[_altitude.alt_max()];
       totdens[s][iz]  = _mixture.neutral_molar_fraction()[s][iz] * _mixture.total_density()[iz];
       unsigned int izb = iz;
       for(CoeffType z = _altitude.alt_max() - _altitude.alt_step(); z >= _altitude.alt_min(); z -= _altitude.alt_step())
       {
         iz = _altitude.altitude_map()[z];
         totdens[s][iz] = totdens[s][izb] + _mixture.neutral_molar_fraction()[s][iz] * _mixture.total_density()[iz];
         izb = iz;
       }
     }

//tau
     _hv_tau.update_tau(a, totdens, sigma);

//finally
     for(unsigned int iz = 0; iz < _altitude.altitudes().size(); iz++)
     {
       VectorCoeffType flux;
       flux.resize(_phy_at_top.abscissa().size());
       for(unsigned int ilambda = 0; ilambda < _phy_at_top.abscissa().size(); ilambda++)
       {
         flux[ilambda] = _phy_at_top.flux()[ilambda] * Antioch::ant_exp(- _hv_tau.tau()[iz][ilambda]);
       }
       _phy[nalt].set_flux(flux)
     }

     return; 
  }

  template<typename CoeffType, typename MatrixCoeffType>
  inline
  void PhotonOpacity<CoeffType,MatrixCoeffType>::set_photon_flux_top_atmosphere(const VectorCoeffType &lambda, 
                                                                             const MatrixCoeffType& phyAU, 
                                                                             const CoeffType& dSunTopAtm)
  {
     return;
  }
}

#endif
