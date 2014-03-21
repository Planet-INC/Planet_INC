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

#ifndef PLANET_PLANET_INITIAL_GUESS_H
#define PLANET_PLANET_INITIAL_GUESS_H

// libMesh
#include "libmesh/function_base.h"

// Planet
#include "planet/planet_physics_helper.h"

namespace Planet
{
  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  class PlanetInitialGuess : public libMesh::FunctionBase<CoeffType>
  {
  public:

    PlanetInitialGuess(const PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>& physics_helper);
    virtual ~PlanetInitialGuess();
    
    virtual libMesh::AutoPtr<libMesh::FunctionBase<CoeffType> > clone() const;

    virtual CoeffType component( unsigned int i, const libMesh::Point& p, libMesh::Real time = 0.0 );

    virtual void operator()(const libMesh::Point& p, const libMesh::Real time,
                            libMesh::DenseVector<CoeffType>& output);
    
    virtual CoeffType operator()(const libMesh::Point& p, const libMesh::Real time=0.);

  protected:

    const PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>& _physics_helper;

  private:

    PlanetInitialGuess();

  };

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  PlanetInitialGuess<CoeffType,VectorCoeffType,MatrixCoeffType>::PlanetInitialGuess(const PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>& physics_helper)
    : _physics_helper(physics_helper)
  {
    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  PlanetInitialGuess<CoeffType,VectorCoeffType,MatrixCoeffType>::~PlanetInitialGuess()
  {
    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  libMesh::AutoPtr<libMesh::FunctionBase<CoeffType> > PlanetInitialGuess<CoeffType,VectorCoeffType,MatrixCoeffType>::clone() const
  {
    return libMesh::AutoPtr<libMesh::FunctionBase<CoeffType> >( new PlanetInitialGuess<CoeffType,VectorCoeffType,MatrixCoeffType>(*this) );
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  CoeffType PlanetInitialGuess<CoeffType,VectorCoeffType,MatrixCoeffType>::operator()(const libMesh::Point& /*p*/, const libMesh::Real /*time*/)
  {
    libmesh_not_implemented();
    // Dummy
    return 0.0;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  CoeffType PlanetInitialGuess<CoeffType,VectorCoeffType,MatrixCoeffType>::component( unsigned int i, const libMesh::Point& p, libMesh::Real /*time*/ )
  {
    CoeffType value = 0.0;

    CoeffType z = p(0);

    value = _physics_helper.first_guess(i,z);

    return value;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetInitialGuess<CoeffType,VectorCoeffType,MatrixCoeffType>::operator()(const libMesh::Point& p, const libMesh::Real time,
                                                                                 libMesh::DenseVector<CoeffType>& output)
  {
    CoeffType z = p(0);
    
    unsigned int n_species = _physics_helper.composition().neutral_composition().n_species();

    for( unsigned int s = 0; s < n_species; s++ )
      {
        CoeffType value = _physics_helper.first_guess(s,z);
        output(s) = value;
      }

    return;
  }

} // end namespace Planet

#endif // PLANET_PLANET_INITIAL_GUESS_H
