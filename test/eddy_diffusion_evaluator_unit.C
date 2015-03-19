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

//Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/physical_constants.h"
#include "antioch/sigma_bin_converter.h"
#include "antioch/vector_utils.h"

//Planet
#include "planet/eddy_diffusion_evaluator.h"
#include "planet/planet_constants.h"

//C++
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>


template<typename Scalar>
int check_test(Scalar theory, Scalar cal, const std::string &words)
{
  Scalar coeff = (std::numeric_limits<Scalar>::epsilon() > 1e-18L)?10.L:20.L;
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * coeff;

  Scalar dist = (std::abs(theory) < tol)?std::abs(cal - theory):
                                         std::abs((theory-cal)/theory);

  if(dist < tol)return 0;
  std::cout << std::scientific << std::setprecision(20)
            << "failed test: " << words << "\n"
            << "theory: " << theory
            << "\ncalculated: " << cal
            << "\ndifference: " << std::abs((theory-cal)/cal)
            << "\ntolerance: " << tol << std::endl;
  return 1;
}

template<typename Scalar, typename VectorScalar = std::vector<Scalar> >
void read_temperature(VectorScalar &T0, VectorScalar &Tz, const std::string &file)
{
  T0.clear();
  Tz.clear();
  std::string line;
  std::ifstream temp(file);
  getline(temp,line);
  while(!temp.eof())
  {
     Scalar t,tz;
     temp >> t >> tz;
     if(temp.good())
     {
       T0.push_back(t);
       Tz.push_back(tz);
     }
  }
  temp.close();
  return;
}

template<typename Scalar>
Scalar barometry(const Scalar &zmin, const Scalar &z, const Scalar &T, const Scalar &Mm, const Scalar &botdens)
{
   return botdens * Antioch::ant_exp(-(z - zmin)/((Planet::Constants::Titan::radius<Scalar>() + z) * (Planet::Constants::Titan::radius<Scalar>() + zmin) * 1e3 *
                                             Antioch::Constants::Avogadro<Scalar>() * Planet::Constants::Universal::kb<Scalar>() * T / 
                                                        (Planet::Constants::Universal::G<Scalar>() * Planet::Constants::Titan::mass<Scalar>() * Mm))
                              );
}

template <typename Scalar>
int tester(const std::string &input_T)
{

// we require an atmospheric mixture,
// which require a temperature and two
// chemical mixture

//description
  std::vector<std::string> neutrals;
  std::vector<std::string> ions;
  neutrals.push_back("N2");
  neutrals.push_back("CH4");
//ionic system contains neutral system
  ions = neutrals;
  ions.push_back("N2+");
  Scalar MN(14.008e-3L), MC(12.011e-3L), MH(1.008e-3L);
  Scalar MN2 = 2.L*MN , MCH4 = MC + 4.L*MH;
  std::vector<Scalar> Mm;
  Mm.push_back(MN2);
  Mm.push_back(MCH4);

//densities
  std::vector<Scalar> molar_frac;
  molar_frac.push_back(0.96L);
  molar_frac.push_back(0.04L);
  molar_frac.push_back(0.L);
  Scalar dens_tot(1e12L);

//altitudes
  Scalar zmin(600.L),zmax(1400.L),zstep(10.L);

//eddy
  Scalar K0(4.3e6L);//

//neutrals
  Antioch::ChemicalMixture<Scalar> neutral_species(neutrals); 

//ions
  Antioch::ChemicalMixture<Scalar> ionic_species(ions); 

//temperature
  std::vector<Scalar> T0,Tz;
  read_temperature<Scalar>(T0,Tz,input_T);
  Planet::AtmosphericTemperature<Scalar, std::vector<Scalar> > temperature(Tz, T0);

//atmospheric mixture
  Planet::AtmosphericMixture<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > composition(neutral_species, ionic_species, temperature);
  composition.init_composition(molar_frac,dens_tot,zmin,zmax);

//eddy diffusion
  Planet::EddyDiffusionEvaluator<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > eddy_diff(composition,K0);

/************************
 * checks
 ************************/
  molar_frac.pop_back();

  Scalar mean_M;
  Antioch::set_zero(mean_M);
  for(unsigned int s = 0; s < molar_frac.size(); s++)
  {
    mean_M += molar_frac[s] * Mm[s];
  }

  int return_flag(0);


// along lines
  for(Scalar z = zmin; z <= zmax; z += zstep)
  {
     std::stringstream walt;
     walt << z;

     Scalar T     = temperature.neutral_temperature(z);
     Scalar nTot  = barometry(zmin,z,T,mean_M,dens_tot);
     Scalar K     = K0 * Antioch::ant_sqrt(dens_tot/nTot);
     Scalar dK_dT = 0.L;
     Scalar dK_dn = - 0.5L * K / nTot;
     return_flag  = check_test(K, eddy_diff.K(nTot), "eddy diffusion at altitude " + walt.str()) ||
                    check_test(dK_dT, eddy_diff.K_deriv_T(T), "eddy diffusion derived by T at altitude " + walt.str()) ||
                    check_test(dK_dn, eddy_diff.K_deriv_ns(nTot), "eddy diffusion derived by n at altitude " + walt.str()) ||
                    return_flag;
  }

// golden values
  Scalar T    = 250.1L;
  Scalar ntot = 2.75e9L;
  Scalar K_theo     =  8.19977826751209391284227292513405759e7L;
  Scalar dK_theo_dT =  0.L;
  Scalar dK_theo_dn = -1.49086877591128980233495871366073774e-2L;

  return_flag  = check_test(K_theo, eddy_diff.K(ntot), "eddy diffusion") ||
                 check_test(dK_theo_dT, eddy_diff.K_deriv_T(T), "eddy diffusion derived by T") ||
                 check_test(dK_theo_dn, eddy_diff.K_deriv_ns(ntot), "eddy diffusion derived by n") ||
                 return_flag;

  return return_flag;
}

int main(int argc, char** argv)
{
  // Check command line count.
  if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify input file." << std::endl;
      antioch_error();
    }

  return (tester<float>(std::string(argv[1])) ||
          tester<double>(std::string(argv[1])) ||
          tester<long double>(std::string(argv[1])));
}
