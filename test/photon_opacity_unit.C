//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

//Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/cmath_shims.h"
#include "antioch/sigma_bin_converter.h"
#include "antioch/vector_utils.h"
//Planet
#include "planet/photon_opacity.h"
//C++
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

template<typename Scalar>
int check(const Scalar &test, const Scalar &ref, const Scalar &tol, const std::string &model)
{
  if(Antioch::ant_abs(test - ref)/ref > tol)
  {
     std::cout << std::scientific << std::setprecision(20)
               << "Error in binary diffusion calculations" << std::endl
               << "model is " << model << std::endl
               << "calculated coefficient = " << test << std::endl
               << "solution = " << ref << std::endl
               << "relative error = " << Antioch::ant_abs(test - ref)/ref << std::endl
               << "tolerance = " << tol << std::endl;
     return 1;
  }
  return 0;
}

template<typename Scalar, typename VectorScalar = std::vector<Scalar> >
void read_crossSection(const VectorScalar & lambda, const std::string &file, unsigned int nbr, VectorScalar &sigma)
{
  std::string line;
  std::ifstream sig_f(file);
  getline(sig_f,line);
  VectorScalar sigma_tmp;
  VectorScalar lambda_tmp;
  while(!sig_f.eof())
  {
     Scalar wv,sigt,sigbr;
     sig_f >> wv >> sigt;
     for(unsigned int i = 0; i < nbr; i++)sig_f >> sigbr;
     lambda_tmp.push_back(wv/10.);//A -> nm
     sigma_tmp.push_back(sigt*10.);//cm-2/A -> m-2/nm
  }
  sig_f.close();

  Antioch::SigmaBinConverter<VectorScalar> bin_converter;
  bin_converter.y_on_custom_grid(lambda_tmp,sigma_tmp,lambda,sigma);
  return;
}

template<typename Scalar>
int tester()
{
  Planet::Altitude<Scalar,std::vector<Scalar> > altitude(600,1400,10);
  Scalar chi(120);
  Planet::Chapman<Scalar> chapman(chi);
  Planet::PhotonOpacity<Scalar,std::vector<Scalar>, std::vector<std::vector<Scalar> > > tau(altitude,chapman);
////cross-section
  std::vector<Scalar> lambda;
  for(Scalar l = 0.1; l < 300.; l += 1.)
  {
      lambda.push_back(l);
  }
  std::vector<std::vector<Scalar> > sigma;
  std::vector<Scalar> sigma_tmp;
  read_crossSection<Scalar>(lambda,"./input/N2_hv_cross-sections.dat",3,sigma_tmp);
  sigma.push_back(sigma_tmp);
  read_crossSection<Scalar>(lambda,"./input/CH4_hv_cross-sections.dat",9,sigma_tmp);
  sigma.push_back(sigma_tmp);
  std::vector<Scalar> a;
  a.resize(altitude.altitudes().size(),45.);
//totdens
  std::vector<std::vector<Scalar> > totdens;
  totdens.resize(2);
  totdens[0].resize(altitude.altitudes().size());
  totdens[1].resize(altitude.altitudes().size());
  totdens[0][0] = 0.96L * 1e9 * Scalar(altitude.altitudes().size());
  totdens[1][0] = 0.04L * 1e9 * Scalar(altitude.altitudes().size());
  for(unsigned int iz = 1; iz < altitude.altitudes().size(); iz++)
  {
      totdens[0][iz] = totdens[0][iz - 1] - 0.96L * 1e9;
      totdens[1][iz] = totdens[1][iz - 1] - 0.04L * 1e9;
  }

  tau.update_tau(a,totdens,sigma);

  int return_flag(0);
  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100.;
  
  for(unsigned int iz = 0; iz < altitude.altitudes().size(); iz++)
  {
     for(unsigned int il = 0; il < lambda.size(); il++)
     {
        Scalar tau_exact(0.L);
        for(unsigned int s = 0; s < 2; s++)
        {
           tau_exact += sigma[s][il] * totdens[s][iz];
        }
        tau_exact *= chapman(a[iz]);
        return_flag = return_flag || 
                      check(tau.tau()[iz][il],tau_exact,tol,"tau at altitude and wavelength");
     }
  }
  return return_flag;
}


int main()
{
  return (tester<float>()  ||
          tester<double>() ||
          tester<long double>());
}
