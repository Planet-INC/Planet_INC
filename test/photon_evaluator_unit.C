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
#include "planet/photon_evaluator.h"
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
  Scalar coeff = (std::numeric_limits<Scalar>::epsilon() < 1e-12)?1e5:1e3;
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * coeff;
  Scalar criteria = (theory < tol)?std::abs(theory-cal):std::abs((theory-cal)/theory);
  if(criteria < tol)return 0;
  std::cout << std::scientific << std::setprecision(20)
            << "failed test: " << words << "\n"
            << "theory: " << theory
            << "\ncalculated: " << cal
            << "\ndifference: " << std::abs((theory-cal)/cal)
            << "\ntolerance: " << tol << std::endl;
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

template<typename Scalar, typename VectorScalar = std::vector<Scalar> >
void read_crossSection(const std::string &file, VectorScalar &lambda, VectorScalar &sigma)
{
  std::string line;
  std::ifstream sig_f(file);
  getline(sig_f,line);
  while(!sig_f.eof())
  {
     Scalar wv(-1.),sigt;
     sig_f >> wv >> sigt;
     if(!getline(sig_f,line))break;
     if(wv < 0.)break;
     lambda.push_back(wv);//A
     sigma.push_back(sigt);//cm-2/A
  }
  sig_f.close();

  return;
}

template<typename Scalar, typename VectorScalar = std::vector<Scalar> >
void read_hv_flux(VectorScalar &lambda, VectorScalar &phy1AU, const std::string &file)
{
  std::string line;
  std::ifstream flux_1AU(file);
  getline(flux_1AU,line);
  while(!flux_1AU.eof())
  {
     Scalar wv,ir,dirr;
     flux_1AU >> wv >> ir >> dirr;
     if(!lambda.empty() && wv == lambda.back())continue;
     lambda.push_back(wv * 10.L);//nm -> A
     phy1AU.push_back(ir * 1e3L * (wv*1e-9L) / (Antioch::Constants::Planck_constant<Scalar>() * 
                                        Antioch::Constants::light_celerity<Scalar>()));//W/m2/nm -> J/s/cm2/A -> s-1/cm-2/A
  }
  flux_1AU.close();
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

template<typename Scalar, typename VectorScalar>
void calculate_densities(VectorScalar &densities, VectorScalar &sum_dens, 
                        const Scalar &tot_dens, const VectorScalar &molar_frac, const Scalar &Mmean, 
                        const Scalar &zmin,const Scalar &zmax,const Scalar &z, 
                        const Planet::AtmosphericTemperature<Scalar,VectorScalar> &T)
{

   densities.clear();
   densities.resize(molar_frac.size(),0.L);
   Scalar nTot = barometry(zmin,z,T.neutral_temperature(z),Mmean,tot_dens);
   sum_dens.clear();
   sum_dens.resize(molar_frac.size(),0.L);
   Scalar zstep(1.L);
   for(unsigned int s = 0; s < molar_frac.size(); s++)
   {
     densities[s] = molar_frac[s] * nTot;
     for(Scalar ztmp = zmax; ztmp >  z; ztmp -= zstep)
     {
        Scalar nTottmp = barometry(zmin,ztmp,T.neutral_temperature(ztmp),Mmean,tot_dens);
        sum_dens[s] += molar_frac[s] * nTottmp * zstep;
     }
   }

   return;
}

template<typename Scalar>
Scalar a(const Scalar &T, const Scalar &Mmean, const Scalar &z)
{
  return Planet::Constants::g<Scalar>(Planet::Constants::Titan::radius<Scalar>(),z,Planet::Constants::Titan::mass<Scalar>()) * Mmean * 
            (Planet::Constants::Titan::radius<Scalar>() + z) * Scalar(1e3) / (Antioch::Constants::Avogadro<Scalar>() * Planet::Constants::Universal::kb<Scalar>() * T);
}

template<typename Scalar, typename VectorScalar>
void calculate_tau(VectorScalar &opacity, const Planet::Chapman<Scalar> &chapman, 
                   const std::vector<VectorScalar*> &cs, const VectorScalar &lambda_ref,
                   const VectorScalar &sum_dens, const Scalar &a)
{
  opacity.resize(lambda_ref.size());
  Antioch::SigmaBinConverter<VectorScalar> bin_converter;

  for(unsigned int s = 0; s < sum_dens.size(); s++)
  {
      VectorScalar sigma_process;
      bin_converter.y_on_custom_grid(*(cs[2*s]),*(cs[2*s+1]),lambda_ref,sigma_process);
      for(unsigned int il = 0; il < lambda_ref.size(); il++)
      {
        opacity[il] += sigma_process[il] * sum_dens[s];
      }
  }
  for(unsigned int il = 0; il < lambda_ref.size(); il++)
  {
      opacity[il] *= chapman(a) * 1e5; //cm-1.km to no unit
  }

}


template <typename Scalar>
int tester(const std::string &input_T, const std::string &input_hv, 
           const std::string &input_N2, const std::string &input_CH4)
{
//description
  std::vector<std::string> neutrals;
  std::vector<std::string> ions;
  neutrals.push_back("N2");
  neutrals.push_back("CH4");
//ionic system contains neutral system
  ions = neutrals;
  ions.push_back("N2+");
  Scalar MN(14.008L), MC(12.011), MH(1.008L);
  Scalar MN2 = 2.L*MN , MCH4 = MC + 4.L*MH;
  std::vector<Scalar> Mm;
  Mm.push_back(MN2);
  Mm.push_back(MCH4);

//densities
  std::vector<Scalar> molar_frac;
  molar_frac.push_back(0.96L);
  molar_frac.push_back(0.04L);
  molar_frac.push_back(0.L);
  const Scalar dens_tot(1e12L);
  Scalar Mmean(0.L);
  for(unsigned int s = 0; s < Mm.size(); s++)
  {
     Mmean += Mm[s] * molar_frac[s];
  }
  Mmean *= 1e-3; //to kg

//zenith angle
  Scalar chi(120);

//photon flux
  std::vector<Scalar> lambda_hv,phy1AU;
  read_hv_flux<Scalar>(lambda_hv,phy1AU,input_hv);
  Antioch::ParticleFlux<std::vector<Scalar> > phy_at_top;
  phy_at_top.set_abscissa(lambda_hv);
  std::vector<Scalar> flux_at_top(phy1AU.size(),0.);
  for(unsigned int il = 0; il < phy1AU.size(); il++)
  {
     flux_at_top[il] = phy1AU[il] / (Planet::Constants::Saturn::d_Sun<Scalar>() * Planet::Constants::Saturn::d_Sun<Scalar>());
  }
  phy_at_top.set_flux(flux_at_top);

////cross-section
  std::vector<Scalar> lambda_N2,sigma_N2;
  std::vector<Scalar> lambda_CH4, sigma_CH4;
  read_crossSection<Scalar>(input_N2,lambda_N2,sigma_N2);
  read_crossSection<Scalar>(input_CH4,lambda_CH4,sigma_CH4);

//altitudes
  Scalar zmin(600.),zmax(1400.),zstep(10.);

/************************
 * first level
 ************************/

//neutrals
  Antioch::ChemicalMixture<Scalar> neutral_species(neutrals); 

//ions
  Antioch::ChemicalMixture<Scalar> ionic_species(ions); 

//chapman
  Planet::Chapman<Scalar> chapman(chi);

//binary diffusion
//not needed

/************************
 * second level
 ************************/

//temperature
  std::vector<Scalar> T0,Tz;
  read_temperature<Scalar>(T0,Tz,input_T);
  Planet::AtmosphericTemperature<Scalar, std::vector<Scalar> > temperature(T0, T0, Tz, Tz);

//photon opacity
  Planet::PhotonOpacity<Scalar,std::vector<Scalar> > tau(chapman);
  tau.add_cross_section(lambda_N2,  sigma_N2,  Antioch::Species::N2, neutral_species.active_species_name_map().at("N2"));
  tau.add_cross_section(lambda_CH4, sigma_CH4, Antioch::Species::CH4, neutral_species.active_species_name_map().at("CH4"));
  tau.update_cross_section(lambda_hv);

//reaction sets
//not needed

/************************
 * third level
 ************************/

//atmospheric mixture
  Planet::AtmosphericMixture<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > composition(neutral_species, ionic_species, temperature);
  composition.init_composition(molar_frac,dens_tot,zmin,zmax);

//kinetics evaluators
//not needed

/************************
 * fourth level
 ************************/

//photon evaluator
  Planet::PhotonEvaluator<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > photon(phy_at_top,tau,composition);
//here, don't you forget to set the photon flux pointers to the reaction set

//molecular diffusion
//not needed

//eddy diffusion
//not needed

/************************
 * checks
 ************************/

  molar_frac.pop_back();//get the ion outta here

  std::vector<std::vector<Scalar>*> cs;
  cs.push_back(&lambda_N2);
  cs.push_back(&sigma_N2);
  cs.push_back(&lambda_CH4);
  cs.push_back(&sigma_CH4);

  int return_flag(0);

  std::vector<Scalar> flux_at_z(lambda_hv.size(),0.);

  for(Scalar z = zmin; z <= zmax; z += zstep)
  {
    std::vector<Scalar> densities, sum_dens;

    calculate_densities(densities, sum_dens, dens_tot, molar_frac, Mmean, zmin, zmax, z, temperature);

    Scalar T = temperature.neutral_temperature(z);
    Scalar x = a(T,Mmean,z);
    std::vector<Scalar> opacity;
    calculate_tau(opacity,chapman,cs,lambda_hv,sum_dens,x);


    photon.update_photon_flux(densities,sum_dens,z,flux_at_z);
    

    for(unsigned int il = 0; il < lambda_hv.size(); il++)  
    {
        Scalar phy_top  = phy1AU[il] / (Planet::Constants::Saturn::d_Sun<Scalar>() * Planet::Constants::Saturn::d_Sun<Scalar>());
        Scalar phy_theo = phy_top * Antioch::ant_exp(-opacity[il]);
   
        int flag_phy_top = check_test(phy_top, photon.photon_flux_at_top().flux()[il], "phy at top at altitude and wavelength");
        int flag_phy     =  check_test(phy_theo, flux_at_z[il], "phy at altitude and wavelength");
        return_flag = flag_phy_top || flag_phy  || return_flag;
    }
  }

  return return_flag;
}

int main(int argc, char** argv)
{
  // Check command line count.
  if( argc < 5 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify inputs file." << std::endl;
      antioch_error();
    }

  return (tester<float>(std::string(argv[1]),std::string(argv[2]),std::string(argv[3]),std::string(argv[4])) ||
          tester<double>(std::string(argv[1]),std::string(argv[2]),std::string(argv[3]),std::string(argv[4])) ||
          tester<long double>(std::string(argv[1]),std::string(argv[2]),std::string(argv[3]),std::string(argv[4])));
}
