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
#include "antioch/physical_constants.h"

//Planet
#include "planet/atmosphere.h"
#include "planet/planet_constants.h"

//C++
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>


template<typename Scalar>
int check_test(Scalar theory, Scalar cal, const std::string &words)
{
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100.;
  Scalar criteria = (theory == 0.L)?std::abs(cal):std::abs((theory-cal)/theory);
  if(criteria < tol)return 0;
  std::cout << "failed test:\n"
            << std::scientific  << std::setprecision(25)
            << "theory: "       << theory
            << "\ncalculated: " << cal
            << "\ndifference: " << criteria
            << "\ntolerance: "  << tol << std::endl;
  std::cout << words << std::endl;
  return 1;
}

template <typename Scalar>
int test()
{
  std::vector<std::string> neutrals;
  std::vector<std::string> ions;
  neutrals.push_back("N2");
  neutrals.push_back("CH4");
  ions.push_back("N2+");
  ions.push_back("e");
  Scalar chi(100.L);
  std::vector<Scalar> lambda_ref,phy1AU,phy_on_top;
  std::ifstream flux_1AU("./input/hv_SSI.dat");
  std::string line;
  getline(flux_1AU,line);
  while(!flux_1AU.eof())
  {
     Scalar wv,ir,dirr;
     flux_1AU >> wv >> ir >> dirr;
     if(!lambda_ref.empty() && wv == lambda_ref.back())continue;
     lambda_ref.push_back(wv);//nm
     phy1AU.push_back(ir);//W/m2/nm
     phy_on_top.push_back(ir/std::pow(Planet::Constants::Saturn::d_Sun<Scalar>(),2));
  }
  flux_1AU.close();

  Planet::Chapman<Scalar> chapman(chi);
  Planet::PhotonFlux<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > photon(chapman);

  photon.set_photon_flux_top_atmosphere(lambda_ref,phy1AU,Planet::Constants::Saturn::d_Sun<Scalar>());

  std::vector<Scalar> molar_frac = {0.96,0.04,0.,0.};
  Scalar dens_tot(1e12);
  Planet::Atmosphere<Scalar, std::vector<Scalar>, std::vector<std::vector<Scalar> > > atm(neutrals,ions,photon);
  atm.make_altitude_grid(600.,1400.,10.);

  std::ifstream temp("input/temperature.dat");
  std::vector<Scalar> T0,Tz;
  getline(temp,line);
  while(!temp.eof())
  {
     Scalar t,tz,dt,dtz;
     temp >> t >> tz >> dt >> dtz;
     T0.push_back(t);
     Tz.push_back(tz);
  }
  temp.close();

  atm.set_temperature(T0,Tz);
  std::vector<Scalar> T = atm.temperature_top_to_bottom();

  std::vector<std::vector<Scalar> > MatrixTotalDensity;
  std::vector<Scalar> tot_dens = atm.total_density_top_to_bottom();

  std::vector<Scalar> lambda_N2,sigma_N2;
  std::vector<std::vector<Scalar> > sigma_rate_N2;
  sigma_rate_N2.resize(3);
  std::ifstream sig_N2("./input/N2_hv_cross-sections.dat");
  std::ifstream sig_CH4("./input/CH4_hv_cross-sections.dat");
  getline(sig_N2,line);
  getline(sig_CH4,line);
  while(!sig_N2.eof())
  {
     Scalar wv,sigt,sig1,sig2,sig3;
     sig_N2 >> wv >> sigt >> sig1 >> sig2 >> sig3;
     lambda_N2.push_back(wv/10.);//A -> nm
     sigma_N2.push_back(sigt*10.);//cm-2/A -> cm-2/nm
     sigma_rate_N2[0].push_back(sig1*10.);
     sigma_rate_N2[1].push_back(sig2*10.);
     sigma_rate_N2[2].push_back(sig3*10.);
  }
  sig_N2.close();
  atm.add_photoabsorption("N2",lambda_N2,sigma_N2);
  std::vector<Scalar> lambda_CH4,sigma_CH4;
  std::vector<std::vector<Scalar> > sigma_rate_CH4;
  sigma_rate_CH4.resize(9);
  while(!sig_CH4.eof())
  {
     Scalar wv,sigt,sig1,sig2,sig3,sig4,sig5,sig6,sig7,sig8,sig9;
     sig_CH4 >> wv >> sigt >> sig1 >> sig2 >> sig3 >> sig4 >> sig5 >> sig6 >> sig7 >> sig8 >> sig9;
     lambda_CH4.push_back(wv/10.);//A -> nm
     sigma_CH4.push_back(sigt*10.);//cm-2/A -> cm-2/nm
     sigma_rate_CH4[0].push_back(sig1*10.);
     sigma_rate_CH4[1].push_back(sig2*10.);
     sigma_rate_CH4[2].push_back(sig3*10.);
     sigma_rate_CH4[3].push_back(sig4*10.);
     sigma_rate_CH4[4].push_back(sig5*10.);
     sigma_rate_CH4[5].push_back(sig6*10.);
     sigma_rate_CH4[6].push_back(sig7*10.);
     sigma_rate_CH4[7].push_back(sig8*10.);
     sigma_rate_CH4[8].push_back(sig9*10.);
  }
  sig_CH4.close();
  atm.add_photoabsorption("CH4",lambda_CH4,sigma_CH4);

  atm.init_composition(molar_frac,dens_tot);

  int return_flag(0);

//Phy at top
  std::vector<Scalar> phy_top;
  phy_top.resize(lambda_ref.size(),0.L);
  for(unsigned int il = 0; il < lambda_ref.size(); il++)
  {
     phy_top[il] = phy1AU[il]/(Planet::Constants::Saturn::d_Sun<Scalar>() * Planet::Constants::Saturn::d_Sun<Scalar>());
     if(check_test(phy_top[il],photon.phy_at_top()[il],"Photon flux at top"))return_flag = 1;
  }

  std::ofstream out("phy_z.dat");
  std::ofstream out_the("phy_z_the.dat");
  out_the << "z N_N N+_N N2+ sCH2_H2 CH3_H CH2_H_H CH4+ CH3+_H CH2+_H2 CH+_H2_H H+_CH3 CH_H2_H" << std::endl;
  std::vector<Scalar> lamb = atm.hv_flux().lambda();
  std::vector<Scalar> sum_over_neutral;
  sum_over_neutral.resize(2,0.L);
  std::vector<Scalar> tau_theo;
  std::vector<Scalar> rate_N2,rate_CH4;
  rate_N2.resize(3);
  rate_CH4.resize(9);
  Scalar MN(14.008L), MC(12.011), MH(1.008L);
  Scalar MN2 = 2.L*MN , MCH4 = MC + 4.L*MH;
  unsigned int ialt(0);
  for(Scalar z = 1400.; z >= 600.; z -= 10.)
  {

    Scalar M_the = molar_frac[0] * MN2 + molar_frac[1] * MCH4;
    Scalar n_tot_the = dens_tot * std::exp(-(z - 600.) / 
                       ( (z + Planet::Constants::Titan::radius<Scalar>())    *
                         (600. + Planet::Constants::Titan::radius<Scalar>()) * 1e3L * //to m
                         ( (Planet::Constants::Universal::kb<Scalar>() * Antioch::Constants::Avogadro<Scalar>() * T[ialt]) /
                           (Planet::Constants::Universal::G<Scalar>() * Planet::Constants::Titan::mass<Scalar>() * M_the * 1e-3L) //to kg/mol
                         )
                       ));
    tau_theo.clear();
    tau_theo.resize(lambda_ref.size(),0.L);
    for(unsigned int i = 0; i < 2; i++)
    {
       sum_over_neutral[i] += molar_frac[i] * n_tot_the;
       for(unsigned int il = 0; il < lambda_ref.size(); il++)
       {
          tau_theo[il] += sum_over_neutral[i] * atm.photon_sigma(i).y_on_custom()[il]; //filtering
       }
    }

    std::vector<Scalar> tau_cal = photon.tau(z,sum_over_neutral);
    for(unsigned int il = 0; il < lambda_ref.size(); il++)
    {
      tau_theo[il] *= chapman.chapman(atm.a(z));
      if(check_test(tau_theo[il],tau_cal[il],"tau at altitude z"))return_flag = 1;
    }

    for(unsigned int ir = 0; ir < 3; ir++)
    {
      rate_N2[ir] = 0.L;
      for(unsigned int il = 0; il < lamb.size(); il++)
      {
        rate_N2[ir] += sigma_rate_N2[ir][il] * phy_on_top[il] * std::exp(-tau_theo[il]);
      }
    }

    for(unsigned int ir = 0; ir < 9; ir++)
    {
      rate_CH4[ir] = 0.L;
      for(unsigned int il = 0; il < lamb.size(); il++)
      {
        rate_CH4[ir] += sigma_rate_CH4[ir][il] * phy_on_top[il] * std::exp(-tau_theo[il]);
      }
    }

    std::vector<Scalar> phy_flux = atm.hv_flux().phy(z);
    for(unsigned int il = 0; il < phy_flux.size(); il++)
    {
       if(check_test(phy_on_top[il] * std::exp(-tau_theo[il]),phy_flux[il],"phy at altitude z and lambda"))return_flag = 1;
       out << lamb[il] << " " << phy_flux[il] << std::endl;
    }

    out_the << z << " ";
    for(unsigned int ir = 0; ir < 3; ir++)
    {
       out_the << rate_N2[ir] << " ";
    }
    for(unsigned int ir = 0; ir < 9; ir++)
    {
       out_the << rate_CH4[ir] << " ";
    }
    out_the << std::endl;
    
    out << std::endl;
    ialt++;
  }

  out.close();
  out_the.close();

  return return_flag;
}

int main()
{

  return (test<float>() ||
          test<double>());// ||
//          test<long double>());
}
