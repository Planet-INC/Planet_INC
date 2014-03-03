//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef PLANET_UNIF_PDF_H
#define PLANET_UNIF_PDF_H

//Antioch
#include "antioch/antioch_asserts.h"

//Planet

//C++

namespace Planet
{
  template <typename CoeffType>
  class UnifPdf: public BasePdf<CoeffType>
  {
      private:
        CoeffType _min;
        CoeffType _max;

      public:
        UnifPdf();
        UnifPdf(const CoeffType &min, const CoeffType &max);
        ~UnifPdf();

        template <typename StateType>
        void set_parameters(const std::vector<StateType> &pars);

        template <typename StateType>
        void set_min(const StateType &min);

        template <typename StateType>
        void set_max(const StateType &max);

        const CoeffType min()   const;

        const CoeffType max()   const;

        const CoeffType value(unsigned int ip = 0) const;
  };

  template <typename CoeffType>
  inline
  UnifPdf<CoeffType>::UnifPdf():
      _min(0.L),
      _max(0.L)
  {
     return;
  }

  template <typename CoeffType>
  inline
  UnifPdf<CoeffType>::UnifPdf(const CoeffType &min, const CoeffType &max):
      _min(min),
      _max(max)
  {
     return;
  }

  template <typename CoeffType>
  inline
  UnifPdf<CoeffType>::~UnifPdf()
  {
     return;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void UnifPdf<CoeffType>::set_min(const StateType &min)
  {
     _min = min;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void UnifPdf<CoeffType>::set_max(const StateType &max)
  {
     _max = max;
  }

  template <typename CoeffType>
  inline
  const CoeffType UnifPdf<CoeffType>::min() const
  {
      return _min;
  }

  template <typename CoeffType>
  inline
  const CoeffType UnifPdf<CoeffType>::max() const
  {
      return _max;
  }

  template <typename CoeffType>
  inline
  const CoeffType UnifPdf<CoeffType>::value(unsigned int ip) const
  {
      return (this->min() + this->max())/2.L;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void UnifPdf<CoeffType>::set_parameters(const std::vector<StateType> &pars)
  {
      antioch_assert_equal_to(pars.size(),2);
      _min = pars[0];
      _max = pars[1];
  }
}

#endif
