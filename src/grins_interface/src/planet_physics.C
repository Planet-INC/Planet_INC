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

// This class
#include "planet/planet_physics.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/cached_values.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"

//Planet
#include "planet/planet_constants.h"

namespace Planet
{
  template <typename CoeffType, typename VectorCoeffType>
  PlanetPhysics<CoeffType,VectorCoeffType>::PlanetPhysics( const GRINS::PhysicsName& physics_name, const GetPot& input )
    : GRINS::Physics(physics_name,input), 
      _n_species( input.vector_variable_size("Physics/Chemistry/species") ),
      _species_FE_family( libMesh::Utility::string_to_enum<libMeshEnums::FEFamily>( input("Physics/Planet/species_FE_family", "LAGRANGE") ) ),
      _species_order( libMesh::Utility::string_to_enum<libMeshEnums::Order>( input("Physics/Planet/species_order", "FIRST") ) )
  {
     _species_var_names.reserve(this->_n_species);
    for( unsigned int i = 0; i < this->_n_species; i++ )
      {
	/*! \todo Make this prefix string an input option */
	std::string var_name = "n_"+std::string(input( "Physics/Chemistry/species", "DIE!", i ));
	_species_var_names.push_back( var_name );
      }

    return;
  }

  template <typename CoeffType, typename VectorCoeffType>
  PlanetPhysics<CoeffType,VectorCoeffType>::~PlanetPhysics()
  {
    return;
  }

  template <typename CoeffType, typename VectorCoeffType>
  void PlanetPhysics<CoeffType,VectorCoeffType>::init_variables( libMesh::FEMSystem* system )
  {
    _species_vars.reserve(this->_n_species);
    for( unsigned int i = 0; i < this->_n_species; i++ )
      {
	_species_vars.push_back( system->add_variable( _species_var_names[i], 
						       this->_species_order, _species_FE_family) );
      }

    return;
  }

  template <typename CoeffType, typename VectorCoeffType>
  void PlanetPhysics<CoeffType,VectorCoeffType>::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    for( unsigned int i = 0; i < this->_n_species; i++ )
      {
	system->time_evolving( _species_vars[i] );
      }

    return;
  }

  template <typename CoeffType, typename VectorCoeffType>
  void PlanetPhysics<CoeffType,VectorCoeffType>::init_context( GRINS::AssemblyContext& context )
  {
    context.get_element_fe(_species_vars[0])->get_JxW();
    context.get_element_fe(_species_vars[0])->get_phi();
    context.get_element_fe(_species_vars[0])->get_dphi();
    context.get_element_fe(_species_vars[0])->get_xyz();

    return;
  }

  template <typename CoeffType, typename VectorCoeffType>
  void PlanetPhysics<CoeffType,VectorCoeffType>::element_time_derivative( bool compute_jacobian,
                                               GRINS::AssemblyContext& context,
                                               GRINS::CachedValues& cache )
  {
    unsigned int n_qpoints = context.get_element_qrule().n_points();

    const unsigned int var = this->_species_vars[0];

    // The number of local degrees of freedom in each variable.
    const unsigned int n_s_dofs = context.get_dof_indices(var).size();
    
    // Element Jacobian * quadrature weights for interior integration.
    const std::vector<libMesh::Real>& JxW =
      context.get_element_fe(var)->get_JxW();
    
    // The pressure shape functions at interior quadrature points.
    const std::vector<std::vector<libMesh::Real> >& s_phi =
      context.get_element_fe(var)->get_phi();

    const std::vector<std::vector<libMesh::RealGradient> >& s_grad_phi =
      context.get_element_fe(var)->get_dphi();
    
    const std::vector<libMesh::Point>& s_qpoint = 
      context.get_element_fe(var)->get_xyz();

    
////here the vectors I need
    std::vector<libMesh::Number> other_altitudes;
    std::vector<std::vector<libMesh::Number> > other_concentrations;
    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        const libMesh::Number r = s_qpoint[qp](0);
        
        libMesh::Real jac = r*r*JxW[qp];

        std::vector<libMesh::Number> molar_concentrations(this->_n_species, 0);
        std::vector<libMesh::Number> dmolar_concentrations_dz(this->_n_species, 0);
        for(unsigned int s=0; s < this->_n_species; s++ )
          {
            molar_concentrations[s] = context.interior_value(this->_species_vars[s],qp);
            dmolar_concentrations_dz[s] = context.interior_gradient(this->_species_vars[s],qp)(0);
          }

        for(unsigned int s=0; s < this->_n_species; s++ )
          {
            const libMesh::Real n_s = context.interior_value(this->_species_vars[s],qp);

            libMesh::DenseSubVector<libMesh::Number> &Fs = 
              context.get_elem_residual(this->_species_vars[s]); // R_{s}

            _helper.compute(molar_concentrations, dmolar_concentrations_dz, // {n}_s, {dn_dz}_s
                            other_altitudes, other_concentrations,  // {n}_s other z, other z
                            r - Constants::Titan::radius<double>() ) ; // z

            libMesh::Real omega = _helper.diffusion_term(s);

            libMesh::Real omega_dot = _helper.chemical_term(s);

            for(unsigned int i=0; i != n_s_dofs; i++)
              {
                Fs(i) += (  omega_dot*s_phi[i][qp] 
                          + 2*omega*n_s*s_phi[i][qp]
                            - omega*n_s*s_grad_phi[i][qp](0) )*jac;

                if( compute_jacobian )
                  {
                    libmesh_not_implemented();
                  }
              }

          }
      }
    
    return;
  }

  template <typename CoeffType, typename VectorCoeffType>
  void PlanetPhysics<CoeffType,VectorCoeffType>::mass_residual( bool compute_jacobian,
                                     GRINS::AssemblyContext& context,
                                     GRINS::CachedValues& cache )
  {
    return;
  }

} // end namespace Planet
