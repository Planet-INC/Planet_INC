//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef _PLANET_THOMAS_SOLVER_
#define _PLANET_THOMAS_SOLVER_

//Antioch
//Planet
//C++

namespace Planet
{
  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  class CrankNicholson
  {
     private:

       //! derive with respect to vec and dx
       void derive(VectorCoeffType &deriv, const VectorCoeffType &vec, const CoeffType &dx);

     public:  
       CrankNicholson();
       ~CrankNicholson();

       void solve(Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType> &atm);

  };

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  void CrankNicholson<CoeffType,VectorCoeffType,MatrixCoeffType>::derive(VectorCoeffType &deriv, const VectorCoeffType &vec, const CoeffType &dx)
  {
     deriv.clear();
     deriv.resize(vec.size(),0.L);
     for(unsigned int ix = 1; ix < vec.size() - 1; ix++)
     {
        deriv[ix] = (vec[ix+1] - vec[i-1]) / (dx + dx);
     }
     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  CrankNicholson<CoeffType,VectorCoeffType,MatrixCoeffType>::CrankNicholson()
  {
     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  CrankNicholson<CoeffType,VectorCoeffType,MatrixCoeffType>::~CrankNicholson()
  {
     return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename MatrixCoeffType>
  inline
  void CrankNicholson<CoeffType,VectorCoeffType,MatrixCoeffType>::solve(Atmosphere<CoeffType,VectorCoeffType,MatrixCoeffType> &atm)
  {

// Matrixes/Vectors of parameters
     MatrixCoeffType species_scale_height;
     MatrixCoeffType neutral_densities;
     MatrixCoeffType diffusion_matrix;
     MatrixCoeffType mean_free_path_matrix;

     VectorCoeffType scale_height = atm.H_bottom_to_top();
     VectorCoeffType K            = atm.K_bottom_to_top();
     VectorCoeffType T            = atm.temperature_bottom_to_top();
     VectorCoeffType nTot         = atm.total_density_bottom_to_top(); 
     VectorCoeffType alpha        = atm.neutral_thermal_coefficient(); 
     VectorCoeffType altitude     = atm.altitudes_bottom_to_top(); 
     VectorCoeffType exobase      = atm.exobase_altitude(); 

     unsigned int Kmaxindex;
     CoeffType Kmax(-1e303);
     for(unsigned int i = 0; i < K.size(); i++)
     {
        if(K[i] < Kmax)
        {
           Kmaxindex = i;
           Kmax = K[i];
        }
     }

     species_scale_height.resize(atm.n_neutral_species());
     neutral_densities.resize(atm.n_neutral_species());
     diffusion_matrix.resize(atm.n_neutral_species());
     for(unsigned int s = 0; s < atm.n_neutral_species(); s++)
     {
       species_scale_height[s]  = atm.H_bottom_to_top(s);
       neutral_densities[s]     = atm.neutral_molar_density_bottom_to_top(s);
       diffusion_matrix[s]      = atm.diffusion_species_bottom_to_top(s);
       mean_free_path_matrix[s] = atm.mean_free_path_bottom_to_top(s);
     }

//derived Vectors
     VectorCoeffType dT_dz;
     VectorCoeffType dK_dz;
     VectorCoeffType dntot_dz;
     VectorCoeffType ddT_ddz;
     VectorCoeffType dscale_height_dz;
//derived Matrixes
     MatrixCoeffType dDs_dz;
     MatrixCoeffType dneu_dz;
     MatrixCoeffType gradH_dz;
     dDs_dz.resize(atm.n_neutral_species());
     dneut_dz.resize(atm.n_neutral_species());
     gradH_dz.resize(atm.n_neutral_species());
                        
     //1-
     //1st order gradient of T, Ds, K, neut, ntot, H and Ha at z(j):
     //2nd order gradient of T at z(j):
     CoeffType dz = atm.step_alt();

     this->derive(dT_dz,T,dz);
     this->derive(dK_dz,K,dz);
     this->derive(dscale_height_dz,scale_height,dz);
     this->derive(dntot_dz,nTot,dz);

     dT_dz[0] = (-CoeffType(3.L) * T[0] + CoeffType(4.L) * T[1] - T[2]) / (dz + dz);
     for(unsigned int i = 1; i < atm.n_altitudes(); i++)
     {
       ddT_ddz[i] = (T[i+1] - CoeffType(2.L) * T[i] + T[i-1]) / (dz * dz);
     }
     for(unsigned int ineu = 0; ineu < atm.n_neutral_species() - 1; ineu++)
     {
        this->derive(dDs_dz[ineu],diffusion_matrix[ineu],dz);
        this->derive(dneut_dz[ineu],neutral_densities[ineu],dz);
        this->derive(gradH_dz[ineu],scale_height[ineu],dz);
     }

//exobase scale height = mean free path
     //N2 exobase used for all the species except H and H2
     CoeffType exobase_alt = exobase[atm.neutral_composition().species_list()[N2]];
     unsigned int iN2 = atm.neutral_composition().species_list()[N2];

     if( ( Antioch::ant_abs(H(iN2,exobase_alt) - atm.mean_free_path(iN2,exobase_alt))/atm.mean_free_path(iN2,exobase_alt)  < 1 ) && 
         (exobase_alt != atm.max_alt() ) && (exobase_alt != atm.min_alt() ) 
       )
     {
        exobase_alt = altitude[exobase_index]; //Exobase altitude in km
     }
     else
     {
        exobase_alt = altitude.back() * Coefftype(9e-1L);   //upper possible limit set for the exobase
        CoeffType tmp(altitude[0] - exobase_alt);
        exobase_index = 0;
        for(unsigned int iz = 1; iz < altitude.size(); iz++)
        {
            if(Antioch::ant_abs(altitude[iz] - exobase_alt) < tmp)
            {
                tmp = Antioch::ant_abs(z - exobase_alt);
                exobase_index = iz;
            }
        }
     }


//matrixes            
     MatrixCoeffType As,Bs,Cs,Csr,a,b,c,d;
     As.resize(atm.n_neutral_species());
     Bs.resize(atm.n_neutral_species());
     Cs.resize(atm.n_neutral_species());
     Csr.resize(atm.n_neutral_species());
     a.resize(atm.n_neutral_species());
     b.resize(atm.n_neutral_species());
     c.resize(atm.n_neutral_species());
     d.resize(atm.n_neutral_species());
     for(unsigned int i = 0; i < atm.n_neutral_species(); i++)
     {
        As[i].resize(atm.n_altitudes(),0.L)
        Bs[i].resize(atm.n_altitudes(),0.L)
        Cs[i].resize(atm.n_altitudes(),0.L)
        Csr[i].resize(atm.n_altitudes(),0.L)
        a[i].resize(atm.n_altitudes(),0.L)
        b[i].resize(atm.n_altitudes(),0.L)
        c[i].resize(atm.n_altitudes(),0.L)
        d[i].resize(atm.n_altitudes(),0.L)
     }

     for(unsigned int nneus = 0; nneus < atm.n_neutrals_species(); nneus++)
     {
        for(unsigned int iz = 1; iz < atm.n_altitudes()-1; iz++)
        {
          CoeffType z = atm.min_alt() + (CoeffType)iz * atm.step_alt();
          As[nneus][iz] = diffusion_matrix[nneus][iz] + K[iz]; //cm2.s-1
          Bs[nneus][iz] = dDs_dz[ineu][iz] + dK_dz[iz] + diffusion_matrix[ineu][iz] * 
                        ( CoeffType(1.L)/scale_height[ineu][iz] + 
                          (CoeffType(1.L) + alpha[nneus] * (CoeffType(1.L) - neutral_densities[nneus][iz]/nTot[iz])) / 
                          T[iz] * dT_dz[iz]
                        ) + K[iz] * (CoeffType(1.L)/scale_height[iz] + CoeffType(1.L) / T[iz] * dT_dz[iz] ); //cm2.s-1.km-1
          Bsr[nneus][iz] = CoeffType(2.L) / (Constants::Titan::radius<CoeffType>() + z ) * 
                          (diffusion_matrix[nneus][iz] + K[iz]); //cm2.s-1.km-1 (additional radial term in spherical coordinates)
          Bs[nneus][iz] += Bsr[nneus][iz];
          Cs[nneus][iz] = dDs_dz[nneus][iz] * 
                          ( CoeffType(1.L)/neutral_scale_height[nneus][iz] + 
                            ( CoeffType(1.L) + alpha[nneus] * (CoeffType(1.L) - neutral_densities[nneus][iz] / nTot[iz])) / 
                              T[iz]  * dT_dz[iz]
                          ) + Ds[nneus][iz] * 
                          ( - gradH_dz[nneus][iz] / Antioch::ant_pow(neutral_scale_height[nneus][iz],2) + alpha[nneus] /
                            nTot[iz] * 
                            (-dneut_dz[nneus][iz] + dnTot_dz[iz](2:lalt) * neutral_densities[nneus][iz] / nTot[iz] ) /
                            T[iz] * dT_dz[iz] + 
                            (CoeffType(1.L) + alpha[nneus] * ( CoeffType(1.L) - neutral_densities[nneus][iz] / nTot[iz])) * 
                            ( (-Antioch::ant_pow(dT_dz[iz],2) / Antioch::ant_pow(T[iz],2) + ddT_ddz[iz] / T[iz] ) )
                          ) + 
                          dK_dz[iz] * ((CoeffType(1.L)/scale_height[iz] + CoeffType(1.L) / T[iz] * dT_dz[nneus]) ) +
                          K[iz] * 
                          ( ( - gradH_dz[iz] / Antioch::ant_pow(scale_height[iz],2) - Antioch::ant_pow(dT_dz[iz],2) /
                                Antioch::ant_pow(T[iz],2) + ddT_ddz[iz] / T[iz]
                            )
                           ); //cm2.s-1.km-2
          Csr[nneus][iz] = CoeffType(2.L) / (Constants::Titan::radius<CoeffType>() + z ) * 
                          ( diffusion_matrix[nneus][iz] * 
                            ( CoeffType(1.L) / neutral_scale_height[nneus][iz] + 
                              (CoeffType(1.L) + alpha[nneus] * (CoeffType(1.L) - neutral_densities[nneus][iz] / nTot[iz])) /
                              T[iz] * dT_dz[nneus][iz]
                            ) + K[iz] * (CoeffType(1.L) / scale_height[iz] + CoeffType(1.L) / T[iz] * dT_dz[iz])
                           );
          Cs[nneus][iz] += Csr[nneus][iz];
          //a*n(j+1)+d*n(j)+b*n(j-1)=c 
          a[nneus][iz] = -As[nneus][iz] / (CoeffType(2.L) * Antioch::ant_pow(dz,2)) - Bs[nneus][iz] / (CoeffType(4.L) * dz); //in s-1.cm2.km-2
          b[nneus][iz] = -As[nneus][iz] / (CoeffType(2.L) * Antioch::ant_pow(dz,2)) + Bs[nneus][iz] / (CoeffType(4.L) * dz); //in s-1.cm2.km-2
          d[nneus][iz] = Antioch::ant_pow(CoeffType(10.),10) / dt + As[nneus][iz] / Antioch::ant_pow(dz,2) - 
                        (Cs[nneus][iz] - Antioch::ant_pow(CoeffType(10.L),10) * Loss[nneus][iz] )/CoeffType(2.L); //in s-1.cm2.km-2  (s-1 -> 10^10 s-1.cm2.km-2)
          c[nneus][iz] = -a[nneus][iz] * neutral_densities[nneus][iz+1] +
                        ( Antioch::ant_pow(CoeffType(10.L),10) / dt * - As[nneus][iz] / Antioch::ant_pow(dz,2) + 
                          (Cs[nneus][iz] - Antioch::ant_pow(10.L,10) * Loss[nneus][iz]) / CoeffType(2.L) 
                        ) * neutral_densities[nneus][iz] - b[nneus][iz] * neutral_densities[nneus][iz-1] + 
                        Antioch::ant_pow(CoeffType(10.L),10) * Prod[nneus][iz]; //in s-1.cm2.km-2.cm-3  (cm-3.s-1 -> 10^10 s-1.cm2.km-2.cm-3)

        } //z loop
     }//neutral loop

     //a.1) Lower Boundary = hydrostatic equation on ntot.
     //Only a slope constraint at the bottom of the profile
     CoeffType ntot1 = nTot[1] * (CoeffType(1.L) + dz / CoeffType(2.L) * (CoeffType(1.L) / scale_height[1] + dT_dz[1] / T[1] ) ) / 
                            (CoeffType(1.L) - dz / CoeffType(2.L) * (CoeffType(1.L) / scale_height[0] + dT_dz[0] / T[0] ));                 

     //For the species with long lifetime, the lower molar fraction boundary condition is kept:
     VectorCoeffType neut1,molfracfloat,Phyfloat;
     neut1.resize(atm.n_neutral_species(),0.L);
     molfracfloat.resize(atm.n_neutral_species(),0.L);
     Phifloat.resize(atm.n_neutral_species(),0.L);
     if(floating(1)>0)
     {
            for(unsigned int lp = 0; lp < floating.size(); lp++)
            {
              Afloat=neut(1:3,floating(lp)).*(Ds(1:3,floating(lp))./H(1:3,floating(lp))+K(1:3)./Ha(1:3));
              Bfloat=1./(Ds(1:3,floating(lp))+K(1:3));
              Cfloat=(neut(1:3,floating(lp))./T(1:3)).*((T(1:3)-T(2:4))/dz).*...
                ((1+alpha(floating(lp)).*(1-neut(1:3,floating(lp))./ntot(1:3))).*Ds(1:3,floating(lp))+K(1:3));
              Phifloat(lp)=(2/dz*(neut(2,floating(lp))-neut(3,floating(lp)))+Bfloat(2)*(Afloat(2)+Cfloat(2))+...
                Bfloat(3)*(Afloat(3)+Cfloat(3)))/(-10^5*(Bfloat(2)+Bfloat(3)*(z(2)/z(3))^2));
              neut0float=neut(2,floating(lp))+1/2*dz*Bfloat(1)*(-10^5*Phifloat(lp)-Afloat(1)-Cfloat(1))+...
                1/2*dz*Bfloat(2)*(-10^5*Phifloat(lp)-Afloat(2)-Cfloat(2));
              molfracfloat(lp)=neut0float/ntot(1); //mean(neut(2:3,floating(lp))./ntot(2:3));
            }
            molfrac0(floating)=molfracfloat;
          }
            
            for lp=1:NbNeutSp
                if molfrac0(i)>0 
                    neut1(lp)=molfrac0(lp)*ntot1;
                else
                    neut1(lp)=neut(2,i)*(1+dz/2*(1/Ha(2)+dTdz(2)/T(2)))/(1-dz/2*(1/Ha(1)+dTdz(1)/T(1)));
                end
            end

 
            //For the species with short lifetime Prod=Loss*neut1
            if Ind_shortlifetime(1)~=0
                for it_slt=1:length(Ind_shortlifetime)
                    neut1(Ind_shortlifetime(it_slt))=Prod(1,Ind_shortlifetime(it_slt))/Loss(1,Ind_shortlifetime(it_slt));
                end
            end         
            c(2,:)=c(2,:)-b(2,:).*neut1;

  }

}

#endif
