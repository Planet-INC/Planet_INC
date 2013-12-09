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

template<typename VectorScalar>
void linear_interpolation(const VectorScalar &temp0, const VectorScalar &alt0,
                          const VectorScalar &alt1, VectorScalar &temp1)
{
  unsigned int j(0);
  typename Antioch::value_type<VectorScalar>::type a;
  typename Antioch::value_type<VectorScalar>::type b;
  temp1.resize(alt1.size());
  for(unsigned int iz = 0; iz < alt1.size(); iz++)
  {
     while(alt0[j] < alt1[iz])
     {
        j++;
        if(!(j < alt0.size()))break;
     }
     if(j == 0)
     {
        Antioch::set_zero(a);
        b = temp0[j];
     }else if(j < alt0.size() - 1)
     {
        a = (temp0[j] - temp0[j-1])/(alt0[j] - alt0[j-1]);
        b = temp0[j] - a * alt0[j];
     }else
     {
        Antioch::set_zero(a);
        b = temp0.back();
     }
     temp1[iz] = a * alt1[iz] + b;
  }
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

//hard sphere radius
  std::vector<Scalar> hard_sphere_radius;
  hard_sphere_radius.push_back(2.0675e-8L * 1e-2L); //N2  in cm -> m
  hard_sphere_radius.push_back(2.3482e-8L * 1e-2L); //CH4 in cm -> m

/*******************************
 * first level
 *******************************/
//altitude
  Planet::Altitude<Scalar,std::vector<Scalar> > altitude(600.,1400.,10.);
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
  std::vector<Scalar> neutral_temperature;
  linear_interpolation(T0,Tz,altitude.altitudes(),neutral_temperature);
  Planet::AtmosphericTemperature<Scalar, std::vector<Scalar> > temperature(neutral_temperature, neutral_temperature, altitude);

/*********************************
 * third level
 *********************************/

//atmospheric mixture
  Planet::AtmosphericMixture<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > composition(neutral_species, ionic_species, altitude, temperature);
  composition.init_composition(molar_frac,dens_tot);
  composition.set_hard_sphere_radius(hard_sphere_radius);

  composition.initialize();

/********************
 * checks
 ********************/

  std::vector<Scalar> exobase;
  Scalar mean_exo(0.L);
  exobase.resize(2);
  std::vector<Scalar> minexobase;
  Scalar mean_minexo(1e303);
  minexobase.resize(2,1e303);

  int return_flag(0);

  for(unsigned int iz = 0; iz < altitude.altitudes().size(); iz++)
  {
    Scalar z = altitude.altitudes()[iz];
    Scalar M_the;
    Antioch::set_zero(M_the);

    for(unsigned int s = 0; s < 2; s++)
    {
       M_the +=  molar_frac[s] * Mm[s];
    }


    Scalar n_tot_the = dens_tot * std::exp(-(z - 600.) / 
                       ( (z + Planet::Constants::Titan::radius<Scalar>())    *
                         (600. + Planet::Constants::Titan::radius<Scalar>()) * 1e3L * //to m
                         ( (Planet::Constants::Universal::kb<Scalar>() * Antioch::Constants::Avogadro<Scalar>()   * temperature.neutral_temperature()[iz]) /
                           (Planet::Constants::Universal::G<Scalar>()  * Planet::Constants::Titan::mass<Scalar>() * M_the * 1e-3L) //to kg/mol
                         )
                       ));

    Scalar free_path_theo(0.L);

    for(unsigned int s = 0; s < 2; s++)
    {
       Scalar scale_height = Planet::Constants::Universal::kb<Scalar>() * temperature.neutral_temperature()[iz] / 
                (Planet::Constants::g<Scalar>(Planet::Constants::Titan::radius<Scalar>(), altitude.altitudes()[iz],Planet::Constants::Titan::mass<Scalar>()) *
                        Mm[s]/Antioch::Constants::Avogadro<Scalar>() * 1e-3L);

       Scalar free_path(1.L);
       Scalar ftmp(0.L);
       for(unsigned int ineu = 0; ineu < 2; ineu++)
       {
          ftmp += 1e6L * n_tot_the * molar_frac[ineu] * Planet::Constants::pi<Scalar>() *
                  Antioch::ant_pow(hard_sphere_radius[ineu] + hard_sphere_radius[s],2) * 
                  Antioch::ant_sqrt(1.L + Mm[s]/Mm[ineu]);
       }
       free_path /= ftmp;
       free_path_theo += free_path * molar_frac[s];
       if(Antioch::ant_abs(free_path - scale_height) < minexobase[s])
       {
          minexobase[s] = Antioch::ant_abs(free_path - scale_height);
          exobase[s] = altitude.altitudes()[iz];
       }

       Scalar Jeans_flux = Jeans(Mm[s]/Antioch::Constants::Avogadro<Scalar>() * Scalar(1e-3), //m (kg)
                                 n_tot_the * molar_frac[s], //n
                                 neutral_temperature[iz],altitude.altitudes()[iz]); //T,alt

//       if(Jeans_flux != 0.)
       {
       return_flag = return_flag ||
                     check_test(Jeans_flux,composition.Jeans_flux(
                                        composition.neutral_composition().M(s)/Antioch::Constants::Avogadro<Scalar>() * Scalar(1e-3), //g/mol -> kg/mol
                                        composition.total_density()[iz] * composition.neutral_molar_fraction()[s][iz], //n, cm-3
                                        temperature.neutral_temperature()[iz],// T (K)
                                        altitude.altitudes()[iz]), //km -> m
                                "Jeans escape flux of species " + 
                                composition.neutral_composition().species_inverse_name_map().at(composition.neutral_composition().species_list()[s]) + 
                                " at altitude");
                                        
       }



       return_flag = return_flag ||
                     check_test(free_path, composition.free_path()[s][iz], "free path of species at altitude") ||
                     check_test(scale_height, composition.scale_height()[s][iz], "scale height of species at altitude");
    }
       

    Scalar H_the = Planet::Constants::Universal::kb<Scalar>() * temperature.neutral_temperature()[iz] /
                   (Planet::Constants::g<Scalar>(Planet::Constants::Titan::radius<Scalar>(), altitude.altitudes()[iz],Planet::Constants::Titan::mass<Scalar>()) *
                    M_the * 1e-3L / Antioch::Constants::Avogadro<Scalar>()); //kb*T/(g(z) * M/Navo)

    if(Antioch::ant_abs(H_the - free_path_theo) < mean_minexo)
    {
        mean_minexo = Antioch::ant_abs(H_the - free_path_theo);
        mean_exo = altitude.altitudes()[iz];
    }

    Scalar a_theo = (Planet::Constants::Titan::radius<Scalar>() + altitude.altitudes()[iz]) / H_the;

    return_flag = return_flag ||
                  check_test(n_tot_the, composition.total_density()[iz], "total density at altitude")                   ||
                  check_test(free_path_theo, composition.atmosphere_free_path()[iz], "atmospheric mean free path at altitude") ||
                  check_test(H_the, composition.atmosphere_scale_height()[iz], "atmospheric scale height at altitude")        ||
                  check_test(a_theo, composition.a_factor()[iz], "a factor at altitude");

  }

  for(unsigned int s = 0; s < 2; s++)
  {
    return_flag = return_flag ||
                  check_test(exobase[s], altitude.altitudes()[composition.exobase()[s]], "exobase of species at altitude");
  }

  return_flag = return_flag ||
                check_test(mean_exo, altitude.altitudes()[composition.atmosphere_exobase()], "exobase of atmosphere at altitude");

  return return_flag;
}

int main(int argc, char** argv)
{
  // Check command line count.
  if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify reaction set XML input file." << std::endl;
      antioch_error();
    }

  return (tester<float>(std::string(argv[1])) ||
          tester<double>(std::string(argv[1])) ||
          tester<long double>(std::string(argv[1])));
}
