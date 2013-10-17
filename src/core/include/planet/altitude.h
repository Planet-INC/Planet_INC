//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef PLANET_ALTITUDE_H
#define PLANET_ALTITUDE_H

//Antioch

//Planet

//C++
#include <map>

namespace Planet
{
  template<typename CoeffType, typename VectorCoeffType>
  class Altitude
  {
     private:
        VectorCoeffType _altitudes;
        std::map<CoeffType,unsigned int> _map_altitudes;
        CoeffType _alt_min;
        CoeffType _alt_max;
        CoeffType _alt_step;

        //! From min, max and step makes the altitudes
        void make_altitudes();

     public:
        Altitude();
        Altitude(const CoeffType &zmin, const CoeffType &zmax, const CoeffType &zstep);
        ~Altitude();

        //! \return minimum altitude
        const CoeffType alt_min() const;

        //! \return maximum altitude
        const CoeffType alt_max() const;

        //! \return step altitude
        const CoeffType alt_step() const;

        //! \return altitude vector
        const VectorCoeffType altitudes() const;

        //! \return altitude map
        const std::map<CoeffType,unsigned int> altitudes_map() const;
  };

  template<typename CoeffType, typename VectorCoeffType>
  inline
  Altitude<CoeffType,VectorCoeffType>::Altitude()
  {
    return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  Altitude<CoeffType,VectorCoeffType>::~Altitude()
  {
    return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  Altitude<CoeffType,VectorCoeffType>::Altitude(const CoeffType &zmin, const CoeffType &zmax, const CoeffType &zstep):
  _alt_min(zmin),
  _alt_max(zmax),
  _alt_step(zstep)
  {
    this->make_altitudes();
    return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  void Altitude<CoeffType,VectorCoeffType>::make_altitudes()
  {
     _altitudes.clear();
     _map_altitudes.clear();
     unsigned int iz(0);
     for(CoeffType z = _alt_min; z <= _alt_max; z += _alt_step)
     {
        _altitudes.push_back(z);
        _map_altitudes[z] = iz;
        iz++;
     }
     return;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  const CoeffType Altitude<CoeffType,VectorCoeffType>::alt_min() const
  {
     return _alt_min;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  const CoeffType Altitude<CoeffType,VectorCoeffType>::alt_step() const
  {
     return _alt_step;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  const CoeffType Altitude<CoeffType,VectorCoeffType>::alt_max() const
  {
     return _alt_max;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  const VectorCoeffType Altitude<CoeffType,VectorCoeffType>::altitudes() const
  {
     return _altitudes;
  }

  template<typename CoeffType, typename VectorCoeffType>
  inline
  const std::map<CoeffType,unsigned int> Altitude<CoeffType,VectorCoeffType>::altitudes_map() const
  {
     return _map_altitudes;
  }

}

#endif
