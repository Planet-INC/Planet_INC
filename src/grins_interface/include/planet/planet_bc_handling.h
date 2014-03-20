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

#ifndef PLANET_PLANET_BC_HANDLING_H
#define PLANET_PLANET_BC_HANDLING_H

// GRINS
#include "grins/bc_handling_base.h"


namespace Planet
{
  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  class PlanetBCHandling : public GRINS::BCHandlingBase
  {
  public:

    PlanetBCHandling( const std::string& physics_name, const GetPot& input,
                      const PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>& physics_helper );

    virtual ~PlanetBCHandling();

    virtual int string_to_int( const std::string& bc_type_in ) const;

    virtual void init_bc_data( const libMesh::FEMSystem& system );

    virtual void init_bc_types( const GRINS::BoundaryID bc_id, 
				const std::string& bc_id_string, 
				const int bc_type, 
				const std::string& bc_vars, 
				const std::string& bc_value, 
				const GetPot& input );

    virtual void user_init_dirichlet_bcs( libMesh::FEMSystem* system, 
					  libMesh::DofMap& dof_map,
					  GRINS::BoundaryID bc_id, 
					  GRINS::BCType bc_type ) const;

    virtual void user_apply_neumann_bcs( AssemblyContext& context,
					 const GRINS::CachedValues& cache,
					 const bool request_jacobian,
					 const GRINS::BoundaryID bc_id,
					 const GRINS::BCType bc_type ) const;
    
  protected:
  
    const PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>& _physics_helper;

  private:

    PlanetBCHandling();

  };

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  PlanetBCHandling<CoeffType,VectorCoeffType,MatrixCoeffType>::PlanetBCHandling( const std::string& physics_name, const GetPot& input,
                                                                                 const PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>& physics_helper )
    : GRINS::BCHandlingBase(physics_name,input),
      _physics_helper(physics_helper)
  {
    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  PlanetBCHandling<CoeffType,VectorCoeffType,MatrixCoeffType>::~PlanetBCHandling()
  {
    return;
  }

} // end namespace Planet

#endif // PLANET_PLANET_BC_HANDLING_H
