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

// libMesh
#include "libmesh/fem_system.h"
#include "libmesh/const_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"

// GRINS
#include "grins/bc_handling_base.h"
#include "grins/assembly_context.h"
#include "grins/neumann_func_obj.h"

// Planet
#include "planet/planet_physics_helper.h"
#include "planet/planet_upper_neumann_bc.h"

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

    virtual void user_apply_neumann_bcs( GRINS::AssemblyContext& context,
					 const GRINS::CachedValues& cache,
					 const bool request_jacobian,
					 const GRINS::BoundaryID bc_id,
					 const GRINS::BCType bc_type ) const;
 
    unsigned int n_species() const;

    const std::string& species_name(unsigned int s) const;

  protected:

    void set_species_bc_type( GRINS::BoundaryID bc_id, int bc_type );

    void set_species_bc_values( GRINS::BoundaryID bc_id, const VectorCoeffType& species_values );

    libMesh::Real get_species_bc_value( GRINS::BoundaryID bc_id, unsigned int species ) const;
  
    const PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>& _physics_helper;

    enum PLANET_BC_TYPES{ LOWER_BOUNDARY_DIRICHLET = 0,
                          UPPER_BOUNDARY_NEUMANN };

    std::vector<GRINS::VariableIndex> _species_vars;

    std::vector<std::pair<GRINS::BoundaryID,GRINS::BCType> > _species_bc_map;

    std::map<GRINS::BoundaryID,VectorCoeffType> _species_bc_values;

  private:

    PlanetBCHandling();

  };

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  PlanetBCHandling<CoeffType,VectorCoeffType,MatrixCoeffType>::PlanetBCHandling( const std::string& physics_name, const GetPot& input,
                                                                                 const PlanetPhysicsHelper<CoeffType,VectorCoeffType,MatrixCoeffType>& physics_helper )
    : GRINS::BCHandlingBase(physics_name),
      _physics_helper(physics_helper)
  {
    _species_vars.resize(this->n_species());

    std::string id_str = "Physics/"+_physics_name+"/species_bc_ids";
    std::string bc_str = "Physics/"+_physics_name+"/species_bc_types";
    std::string var_str = "Physics/"+_physics_name+"/species_bc_variables";
    std::string val_str = "Physics/"+_physics_name+"/species_bc_values";
    
    this->read_bc_data( input, id_str, bc_str, var_str, val_str );

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  PlanetBCHandling<CoeffType,VectorCoeffType,MatrixCoeffType>::~PlanetBCHandling()
  {
    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  int PlanetBCHandling<CoeffType,VectorCoeffType,MatrixCoeffType>::string_to_int( const std::string& bc_type_in ) const
  {
    int bc_type_out;

    if( bc_type_in == "lower_boundary_dirichlet" )
      bc_type_out = LOWER_BOUNDARY_DIRICHLET;

    else if( bc_type_in == "upper_boundary_neumann" )
      bc_type_out = UPPER_BOUNDARY_NEUMANN;

    else
      {
	// Call base class to detect any physics-common boundary conditions
	bc_type_out = BCHandlingBase::string_to_int( bc_type_in );
      }

    return bc_type_out;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetBCHandling<CoeffType,VectorCoeffType,MatrixCoeffType>::init_bc_data( const libMesh::FEMSystem& system )
  {
    for( unsigned int s = 0; s < this->n_species(); s++ )
      {
	_species_vars[s] = system.variable_number( "n_"+this->species_name(s) );
      }

    // Now setup upper boundary Neumann bc
    /*! \todo GRINS needs an init_neumann_bc or something to take care of this */
    for( std::map<GRINS::BoundaryID,GRINS::BCType>::const_iterator it = _neumann_bc_map.begin();
         it != _neumann_bc_map.begin(); ++ it)
      {
        GRINS::BCType bc_type = it->second;
        switch(bc_type)
          {
          case(UPPER_BOUNDARY_NEUMANN):
            {
              GRINS::NBCContainer nbc;
              nbc.set_bc_id(it->first);
              for( unsigned int s = 0; s < this->n_species(); s++ )
                {
                  std::tr1::shared_ptr<GRINS::NeumannFuncObj> func( new PlanetUpperNeumannBC<CoeffType,VectorCoeffType,MatrixCoeffType>( _physics_helper, _species_vars, s ) );

                  nbc.add_var_func_pair( _species_vars[s], func );
                }

              this->attach_neumann_bound_func( nbc );
            }
            break;

          default:
            {
              std::cerr << "Error: Invalid Neumann BC type for " << _physics_name
                        << std::endl;
              libmesh_error();
            }
          }
      }

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetBCHandling<CoeffType,VectorCoeffType,MatrixCoeffType>::init_bc_types( const GRINS::BoundaryID bc_id, 
                                                                                   const std::string& bc_id_string, 
                                                                                   const int bc_type, 
                                                                                   const std::string& bc_vars, 
                                                                                   const std::string& bc_value, 
                                                                                   const GetPot& input )
  {
    switch(bc_type)
      {
      case(LOWER_BOUNDARY_DIRICHLET):
	{
          this->set_species_bc_type( bc_id, bc_type );

          VectorCoeffType species_densities(this->n_species());

          _physics_helper.lower_boundary_dirichlet(species_densities);

          this->set_species_bc_values( bc_id, species_densities );
	}
	break;

      case(UPPER_BOUNDARY_NEUMANN):
        {
          this->set_neumann_bc_type( bc_id, bc_type );
        }
        break;

      default:
	{
	  // Call base class to detect any physics-common boundary conditions
	  BCHandlingBase::init_bc_types( bc_id, bc_id_string, bc_type,
                                         bc_vars, bc_value, input );
	}

      }// End switch(bc_type)

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetBCHandling<CoeffType,VectorCoeffType,MatrixCoeffType>::user_init_dirichlet_bcs( libMesh::FEMSystem* /*system*/,
                                                                                             libMesh::DofMap& dof_map,
                                                                                             GRINS::BoundaryID bc_id,
                                                                                             GRINS::BCType bc_type ) const
  {
    switch( bc_type )
      {
      case(LOWER_BOUNDARY_DIRICHLET):
	{
	  std::set<GRINS::BoundaryID> dbc_ids;
	  dbc_ids.insert(bc_id);

	  for( unsigned int s = 0; s < this->n_species(); s++ )
	    {
	      std::vector<GRINS::VariableIndex> dbc_vars(1,_species_vars[s]);

              libMesh::ConstFunction<libMesh::Number> species_func( this->get_species_bc_value(bc_id,s) );

	      libMesh::DirichletBoundary species_dbc( dbc_ids,
						      dbc_vars,
						      &species_func );

	      dof_map.add_dirichlet_boundary( species_dbc );
	    }
	}
	break;

      default:
	{
	  std::cerr << "Error: Invalid Dirichlet BC type for " << _physics_name
		    << std::endl;
	  libmesh_error();
	}
      } //switch( bc_type )

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetBCHandling<CoeffType,VectorCoeffType,MatrixCoeffType>::user_apply_neumann_bcs( GRINS::AssemblyContext& context,
                                                                                            const GRINS::CachedValues& cache,
                                                                                            const bool request_jacobian,
                                                                                            const GRINS::BoundaryID bc_id,
                                                                                            const GRINS::BCType bc_type ) const
  {
    switch( bc_type )
      {
	// General heat flux from user specified function
      case(UPPER_BOUNDARY_NEUMANN):
	{
          for( unsigned int s = 0; s < this->n_species(); s++ )
	    {
              _bound_conds.apply_neumann_normal( context, cache, request_jacobian, _species_vars[s], -1.0,
                                                 this->get_neumann_bound_func( bc_id, _species_vars[s] ) );
            }
	}
	break;
      default:
	{
	  std::cerr << "Error: Invalid Neumann BC type for " << _physics_name
		    << std::endl;
	  libmesh_error();
	}
      } // End switch

    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  unsigned int PlanetBCHandling<CoeffType,VectorCoeffType,MatrixCoeffType>::n_species() const
  {
    return _physics_helper.composition().neutral_composition().n_species();
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  const std::string& PlanetBCHandling<CoeffType,VectorCoeffType,MatrixCoeffType>::species_name(unsigned int s) const
  {
    return _physics_helper.composition().neutral_composition().chemical_species()[s]->species();
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetBCHandling<CoeffType,VectorCoeffType,MatrixCoeffType>::set_species_bc_type( GRINS::BoundaryID bc_id, int bc_type )
  {
    _species_bc_map.push_back( std::make_pair(bc_id,bc_type) );
    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  void PlanetBCHandling<CoeffType,VectorCoeffType,MatrixCoeffType>::set_species_bc_values( GRINS::BoundaryID bc_id, 
                                                                                           const VectorCoeffType& species_values )
  {
    _species_bc_values[bc_id] = species_values;
    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  libMesh::Real PlanetBCHandling<CoeffType,VectorCoeffType,MatrixCoeffType>::get_species_bc_value( GRINS::BoundaryID bc_id, 
                                                                                                   unsigned int species ) const
  {
    return (_species_bc_values.find(bc_id)->second)[species];
  }

} // end namespace Planet

#endif // PLANET_PLANET_BC_HANDLING_H
