//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

//Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/physical_constants.h"
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
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100.;
  if(std::abs((theory-cal)/theory) < tol)return 0;
  std::cout << "failed test: " << words << "\n"
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
void read_crossSection(const std::string &file, unsigned int nbr, VectorScalar &lambda, VectorScalar &sigma)
{
  std::string line;
  std::ifstream sig_f(file);
  getline(sig_f,line);
  while(!sig_f.eof())
  {
     Scalar wv,sigt,sigbr;
     sig_f >> wv >> sigt;
     for(unsigned int i = 0; i < nbr; i++)sig_f >> sigbr;
     lambda.push_back(wv/10.);//A -> nm
     sigma.push_back(sigt*10.);//cm-2/A -> m-2/nm
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
     lambda.push_back(wv);//nm
     phy1AU.push_back(ir);//W/m2/nm
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

template<typename Scalar, typename VectorScalar, typename MatrixScalar>
void calculate_densities(MatrixScalar &densities, const Scalar &tot_dens, const VectorScalar &molar_frac, 
                        const Scalar &zmin,const Scalar &zmax,const Scalar &zstep, 
                        const VectorScalar &T, const VectorScalar &mm)
{
   unsigned int iz(0);
   Scalar Mm;
   Antioch::set_zero(Mm);
   for(unsigned int s = 0; s < molar_frac.size(); s++)
   {
      Mm += molar_frac[s] * mm[s];
   }
   densities.clear();
   densities.resize(molar_frac.size());
   for(Scalar z = zmin; z <= zmax; z += zstep)
   {
      for(unsigned int s = 0; s < molar_frac.size(); s++)
      {
        densities[s].push_back(molar_frac[s] * barometry(zmin,z,T[iz],Mm,tot_dens));
      }
      iz++;
   }

   return;
}


template <typename Scalar>
int tester()
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
  Scalar dens_tot(1e12L);

//hard sphere radius
  std::vector<Scalar> hard_sphere_radius;
  hard_sphere_radius.push_back(2.0675e-8L * 1e-5L); //N2  in cm -> km
  hard_sphere_radius.push_back(2.3482e-8L * 1e-5L); //CH4 in cm -> km

//zenith angle
  Scalar chi(120);

//photon flux
  std::vector<Scalar> lambda_hv,phy1AU;
  read_hv_flux<Scalar>(lambda_hv,phy1AU,"./input/hv_SSI.dat");

////cross-section
  std::vector<Scalar> lambda_N2,sigma_N2;
  std::vector<Scalar> lambda_CH4, sigma_CH4;
  read_crossSection<Scalar>("./input/N2_hv_cross-sections.dat",3,lambda_N2,sigma_N2);
  read_crossSection<Scalar>("./input/CH4_hv_cross-sections.dat",9,lambda_CH4,sigma_CH4);

//altitudes
  Scalar zmin(600.),zmax(1400.),zstep(10.);

/************************
 * first level
 ************************/

//altitude
  Planet::Altitude<Scalar,std::vector<Scalar> > altitude(zmin,zmax,zstep);

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
  read_temperature<Scalar>(T0,Tz,"input/temperature.dat");
  std::vector<Scalar> neutral_temperature;
  linear_interpolation(T0,Tz,altitude.altitudes(),neutral_temperature);
  Planet::AtmosphericTemperature<Scalar, std::vector<Scalar> > temperature(neutral_temperature, neutral_temperature, altitude);

//photon opacity
  Planet::PhotonOpacity<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > tau(altitude,chapman);

//reaction sets
//not needed

/************************
 * third level
 ************************/

//atmospheric mixture
  Planet::AtmosphericMixture<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > composition(neutral_species, ionic_species, altitude, temperature);
  composition.init_composition(molar_frac,dens_tot);
  composition.set_hard_sphere_radius(hard_sphere_radius);
  composition.initialize();

//kinetics evaluators
//not needed

/************************
 * fourth level
 ************************/

//photon evaluator
  Planet::PhotonEvaluator<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > photon(altitude,tau,composition);
  photon.add_cross_section(lambda_N2,  sigma_N2,  Antioch::Species::N2);
  photon.add_cross_section(lambda_CH4, sigma_CH4, Antioch::Species::CH4);
  photon.set_photon_flux_at_top(lambda_hv, phy1AU, Planet::Constants::Saturn::d_Sun<Scalar>());
  photon.update_photon_flux();

//molecular diffusion
//not needed

//eddy diffusion
//not needed

/************************
 * checks
 ************************/

  std::vector<std::vector<Scalar> > densities;
  std::vector<Scalar> mm;
  mm.push_back(MN2);
  mm.push_back(MCH4);
  calculate_densities(densities,dens_tot,molar_frac,zmin,zmax,zstep,neutral_temperature,mm);

  std::vector<std::vector<Scalar> > opacity;
  calculate_tau(opacity,);


  int return_flag(0);

  for(unsigned int iz = 0; iz < altitude.altitudes().size(); iz++)
  {

  }


  return return_flag;
}

int main()
{

  return (tester<float>()  ||
          tester<double>() ||
          tester<long double>());
}
