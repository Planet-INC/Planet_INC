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
#include "antioch/vector_utils.h"

//Planet
#include "planet/atmospheric_mixture.h"
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
  Scalar coeff(100.);
  if(std::numeric_limits<Scalar>::epsilon() < 1e-17)coeff *= 5.;
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * coeff;
  if(std::abs((theory-cal)/theory) < tol)return 0;
  std::cout << std::scientific  << std::setprecision(20)
            << "failed test: "  << words << "\n"
            << "theory: "       << theory
            << "\ncalculated: " << cal
            << "\ndifference: " << std::abs((theory-cal)/cal)
            << "\ntolerance: "  << tol << std::endl;
  return 1;
}

template<typename Scalar>
Scalar linear_interpolation(const Scalar &z0, const Scalar &z1,
                            const Scalar &T0, const Scalar &T1,
                            const Scalar &z)
{
     Scalar a = (T0 - T1)/(z0 - z1);
     Scalar b = T0 - a * z0;
     return a * z + b;
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
     Scalar t,tz,dt,dtz;
     temp >> t >> tz >> dt >> dtz;
     T0.push_back(t);
     Tz.push_back(tz);
  }
  temp.close();
  return;
}

template <typename Scalar>
Scalar Jeans(const Scalar &m, const Scalar &n, const Scalar &T, const Scalar &z)
{
  Scalar lambda =  Planet::Constants::Universal::G<Scalar>()  * Planet::Constants::Titan::mass<Scalar>() * m / 
                  (Planet::Constants::Universal::kb<Scalar>() * T * (z + Planet::Constants::Titan::radius<Scalar>()) * Scalar(1e3)); 

  return n * Antioch::ant_sqrt(Planet::Constants::Universal::kb<Scalar>() * T / (Scalar(2.) * Planet::Constants::pi<Scalar>() * m))
           * Antioch::ant_exp(-lambda) * (Scalar(1.) + lambda);
}

template <typename Scalar>
int tester(const std::string & input_T)
{
  std::vector<std::string> neutrals;
  std::vector<std::string> ions;
  neutrals.push_back("N2");
  neutrals.push_back("CH4");
//ionic system contains neutral system
  ions = neutrals;
  ions.push_back("N2+");
  Scalar MN(14.008L), MC(12.011L), MH(1.008L);
  Scalar MN2 = 2.L*MN , MCH4 = MC + 4.L*MH;
  std::vector<Scalar> Mm;
  Mm.push_back(MN2);
  Mm.push_back(MCH4);

//densities
  std::vector<Scalar> molar_frac;
  molar_frac.push_back(0.96L);
  molar_frac.push_back(0.04L);
  molar_frac.push_back(0.L);
  Scalar dens_tot(1e12L); //cm-3
  std::vector<Scalar> neutral_molar_concentration;
  neutral_molar_concentration.push_back(9.6e11);
  neutral_molar_concentration.push_back(4e10);

/*******************************
 * first level
 *******************************/
//neutrals
  Antioch::ChemicalMixture<Scalar> neutral_species(neutrals); 
//ions
  Antioch::ChemicalMixture<Scalar> ionic_species(ions); 

/*********************************
 * second level
 *********************************/

//temperature
  std::vector<Scalar> T0,Tz;
  read_temperature<Scalar>(T0,Tz,input_T);
  Planet::AtmosphericTemperature<Scalar, std::vector<Scalar> > temperature(T0, T0, Tz, Tz);

/*********************************
 * third level
 *********************************/

//atmospheric mixture
  Planet::AtmosphericMixture<Scalar,std::vector<Scalar> > composition(neutral_species, ionic_species, temperature);
  composition.init_composition(molar_frac,dens_tot);

/********************
 * checks
 ********************/

  int return_flag(0);

  Scalar zmin = 600.L;
  Scalar zmax = 1400.L;
  Scalar dz   = 10.L;
  for(Scalar z = zmin; z < zmax; z += dz)
  {
    Scalar T = temperature.neutral_temperature(z);

    std::vector<Scalar> scale_heights;
    scale_heights.resize(neutrals.size(),0.L);
    composition.scale_heights(z,scale_heights);

    Scalar M_the(0.L);

    for(unsigned int s = 0; s < neutrals.size(); s++)
    {
       Scalar scale_height = Planet::Constants::Universal::kb<Scalar>() * T / 
                (Planet::Constants::g<Scalar>(Planet::Constants::Titan::radius<Scalar>(), z,Planet::Constants::Titan::mass<Scalar>()) *
                        Mm[s]/Antioch::Constants::Avogadro<Scalar>() * 1e-3L);

       Scalar Jeans_flux = Jeans(Mm[s]/Antioch::Constants::Avogadro<Scalar>() * Scalar(1e-3), //m (kg)
                                 neutral_molar_concentration[s], //n
                                 T,z); //T,alt

       return_flag = return_flag ||
                     check_test(Jeans_flux,composition.Jeans_flux(
                                        composition.neutral_composition().M(s)/Antioch::Constants::Avogadro<Scalar>() * Scalar(1e-3), //g/mol -> kg/mol
                                        neutral_molar_concentration[s], //n, cm-3
                                        T,// T (K)
                                        z), //km -> m
                                "Jeans escape flux of species " + 
                                composition.neutral_composition().species_inverse_name_map().at(composition.neutral_composition().species_list()[s]) + 
                                " at altitude");
                                        
       return_flag = return_flag ||
                     check_test(scale_height, scale_heights[s], "scale height of species at altitude");
       M_the += Mm[s] * molar_frac[s];
    }

    Scalar H_the = Planet::Constants::Universal::kb<Scalar>() * T /
                   (Planet::Constants::g<Scalar>(Planet::Constants::Titan::radius<Scalar>(), z,Planet::Constants::Titan::mass<Scalar>()) *
                    M_the * 1e-3L / Antioch::Constants::Avogadro<Scalar>()); //kb*T/(g(z) * M/Navo)

    Scalar a_the = (Planet::Constants::Titan::radius<Scalar>() + z) / H_the * Scalar(1e3); // to m

    return_flag = return_flag ||
                  check_test(H_the, composition.atmospheric_scale_height(neutral_molar_concentration, z), "atmospheric scale height at altitude");
                  check_test(a_the, composition.a(neutral_molar_concentration, z), "atmospheric a factor at altitude");

  }

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
