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

#ifndef PLANET_ATMOSPHERE_TEMPERATURE_H
#define PLANET_ATMOSPHERE_TEMPERATURE_H

//Antioch
#include "antioch/antioch_asserts.h"

//Planet

//C++


/*before it is merged in Antioch*/
// GSL
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

namespace Planet
{

  /* We copy/paste here the gsl_spliner.h header file of Antioch
     Once it is merged, we need to delete all those lines
  */
  class GSLSpliner;


  /*! \class AtmosphericTemperature
      
      This class contains the models for three temperatures:
        - the neutral species temperature
        - the ionic species temperature
        - the electrons temperature

      We use a spline for the neutral temperature, it is read from
      an input file.
      We use the hypothesis \f$T_\text{ion} = T_\text{neu}\f$.
      We use the simple Titan model:
      \f[
           \begin{split}
             T_\text{e} = 180~\text{K}                             & \text{ if }         z  < 900~\text{km}\\
             T_\text{e} = 180 + \frac{z - 900}{500} * (1150 - 180) & \text{ if } 900 \le z < 1400~\text{km} \\
             T_\text{e} = 1150 + 0.1 * (z - 1400)                  & \text{ else}\\
           \end{split}
      \f]
   */
  template<typename CoeffType, typename VectorCoeffType, typename Spliner = GSLSpliner>
  class AtmosphericTemperature
  {
      private:
        //!no default constructor
        AtmosphericTemperature(){antioch_error();return;}

        Spliner spline;

      public:
        AtmosphericTemperature(const VectorCoeffType &alt_neu, const VectorCoeffType &T_neu);
        ~AtmosphericTemperature();

        //!\return neutral temperature at custom altitude
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
        neutral_temperature(const StateType &z) const
        ANTIOCH_AUTOFUNC(StateType, spline.interpolated_value(z))

        //!\return derivative neutral temperature at custom altitude
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
        dneutral_temperature_dz(const StateType &z) const
        ANTIOCH_AUTOFUNC(StateType, spline.dinterp_dx(z))

        //!\return ionic temperature at custom altitude
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
        ionic_temperature(const StateType &z)   const
        ANTIOCH_AUTOFUNC(StateType, neutral_temperature(z))

        //!\return derivative ionic temperature at custom altitude
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
        dionic_temperature_dz(const StateType &z)   const
        ANTIOCH_AUTOFUNC(StateType, dneutral_temperature_dz(z))

        //!\return electronic temperature at custom altitude
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
        electronic_temperature(const StateType &z)   const
        ANTIOCH_AUTOFUNC(StateType, (z < Antioch::constant_clone(z,900.))?Antioch::constant_clone(z,180.):
                                         (z < Antioch::constant_clone(z,1400.))?Antioch::constant_clone(z,180.) + Antioch::constant_clone(z,2e-3) * (z - Antioch::constant_clone(z,900.)) * Antioch::constant_clone(z,(1150. - 180.)):
                                                                                Antioch::constant_clone(z,1150.) + Antioch::constant_clone(z,0.1) * (z - Antioch::constant_clone(z,1400.))
                        )

        //!\return derivative electronic temperature at custom altitude
        template <typename StateType>
        const ANTIOCH_AUTO(StateType)
        delectronic_temperature_dz(const StateType &z)   const
        ANTIOCH_AUTOFUNC(StateType, (z < Antioch::constant_clone(z,900.))?Antioch::zero_clone(z):
                                         (z < Antioch::constant_clone(z,1400.))?Antioch::constant_clone(z,2e-3):
                                                                                Antioch::constant_clone(z,0.1)
                        )

        //! reset neutral temperature
        void set_neutral_temperature(const VectorCoeffType & alt, const VectorCoeffType &neu);

        //!
        void initialize();

  };

  template<typename CoeffType, typename VectorCoeffType, typename Spliner>
  inline
  AtmosphericTemperature<CoeffType,VectorCoeffType,Spliner>::~AtmosphericTemperature()
  {
    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename Spliner>
  inline
  AtmosphericTemperature<CoeffType,VectorCoeffType,Spliner>::AtmosphericTemperature(const VectorCoeffType &alt_neu, 
                                                                            const VectorCoeffType &T_neu)
        :spline(alt_neu,T_neu)
  {
    return;
  }

  template<typename CoeffType, typename VectorCoeffType, typename Spliner>
  inline
  void AtmosphericTemperature<CoeffType,VectorCoeffType,Spliner>::initialize()
  {
     return; //nothing
  }


  template<typename CoeffType, typename VectorCoeffType, typename Spliner>
  inline
  void AtmosphericTemperature<CoeffType,VectorCoeffType,Spliner>::set_neutral_temperature(const VectorCoeffType & alt_neu, const VectorCoeffType & T_neu)
  {
    spline.spline_delete();
    spline.init(alt_neu,T_neu);
  }



/*
   Here come the gsl_spliner.h header file
*/

  template<bool B>
  struct GSLInterp
  {
      template <typename Scalar>
      Scalar interpolation(const Scalar & x, gsl_spline * spline, gsl_interp_accel * acc)
                {return gsl_spline_eval(spline,x,acc);}

      template <typename Scalar>
      Scalar dinterpolation(const Scalar & x, gsl_spline * spline, gsl_interp_accel * acc)
                {return gsl_spline_eval_deriv(spline,x,acc);}
  };

  template <>
  struct GSLInterp<true>
  {
      template <typename VectorScalar>
      VectorScalar interpolation(const VectorScalar & x, gsl_spline * spline, gsl_interp_accel * acc)
                {
                  VectorScalar out = zero_clone(x);
                  for(unsigned int i =0; i < x.size(); ++i)
                  {
                    out[i] = gsl_spline_eval(spline,x[i],acc);
                  }
                  return out;
                }

      template <typename VectorScalar>
      VectorScalar dinterpolation(const VectorScalar & x, gsl_spline * spline, gsl_interp_accel * acc)
                {
                  VectorScalar out = zero_clone(x);
                  for(unsigned int i =0; i < x.size(); ++i)
                  {
                    out[i] = gsl_spline_eval_deriv(spline,x[i],acc);
                  }
                  return out;
                }
  };

  class GSLSpliner
  {
     public:
       GSLSpliner();
       template <typename VectorCoeffType>
       GSLSpliner(const VectorCoeffType & data_x_point, const VectorCoeffType & data_y_point);
       ~GSLSpliner();

     template <typename VectorCoeffType>
     void spline_init(const VectorCoeffType & data_x_point, const VectorCoeffType & data_y_point);

     void spline_delete();

     template <typename StateType>
     StateType interpolated_value(const StateType & x) const;

     template <typename StateType>
     StateType dinterp_dx(const StateType & x) const;

     private:

       gsl_interp_accel * _acc;
       gsl_spline       * _spline;
  };

  inline
  GSLSpliner::GSLSpliner()
      :_acc(NULL),_spline(NULL)
  {
    _acc =  gsl_interp_accel_alloc();
    return;
  }

  template <typename VectorCoeffType>
  inline
  GSLSpliner::GSLSpliner(const VectorCoeffType & data_x_point, const VectorCoeffType & data_y_point)
        :_acc(NULL),_spline(NULL)
  {
    _acc =  gsl_interp_accel_alloc();
    spline_init(data_x_point, data_y_point);
  }

  GSLSpliner::~GSLSpliner()
  {
    this->spline_delete();
    gsl_interp_accel_free(_acc);
    return;
  }

  inline
  void GSLSpliner::spline_delete()
  {
    gsl_spline_free(_spline);
    return;
  }

  template <typename VectorCoeffType>
  inline
  void GSLSpliner::spline_init(const VectorCoeffType & data_x_point, const VectorCoeffType & data_y_point)
  {
     antioch_assert_equal_to(data_x_point.size(), data_y_point.size());

     _spline = gsl_spline_alloc(gsl_interp_cspline, data_x_point.size());

   // GLS takes only double, raaaaahhhh
     typedef typename Antioch::rebind<VectorCoeffType,double>::type VectorGSLType;
     VectorGSLType gsl_x_point(data_x_point.size(),0);
     VectorGSLType gsl_y_point(data_y_point.size(),0);
     for(unsigned int i = 0; i < data_x_point.size(); i++)
     {
        gsl_x_point[i] = (const double)data_x_point[i];
        gsl_y_point[i] = (const double)data_y_point[i];
     }

     const double * x = &gsl_x_point[0];
     const double * y = &gsl_y_point[0];

     gsl_spline_init(_spline, x, y, data_x_point.size());
  }

  template <typename StateType>
  inline
  StateType GSLSpliner::interpolated_value(const StateType & x) const
  {
     return GSLInterp<Antioch::has_size<StateType>::value>().interpolation(x, _spline, _acc);
  }

  template <typename StateType>
  inline
  StateType GSLSpliner::dinterp_dx(const StateType & x) const
  {
     return GSLInterp<Antioch::has_size<StateType>::value>().dinterpolation(x, _spline, _acc);
  }











}

#endif
