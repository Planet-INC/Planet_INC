//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef PLANET_CROSS_SECTION_H
#define PLANET_CROSS_SECTION_H

//Antioch
#include "antioch/sigma_bin_converter.h"

//Planet

//C++

namespace Planet
{
  template<typename VectorCoeffType>
  class CrossSection
  {
      private:
        VectorCoeffType _abscissa;
        VectorCoeffType _cross_section;
        VectorCoeffType _cross_section_on_custom_grid;

        SigmaBinConverter<VectorCoeffType> _converter;

      public:

        CrossSection();
        CrossSection(const CrossSection<VectorCoeffType> &rhs);
        CrossSection(const VectorCoeffType &x, const VectorCoeffType &y);
        ~CrossSection();

        //!sets the abscissa
        template<typename VectorStateType>
        void set_abscissa(const VectorStateType &x);

        //!sets the cross section
        template<typename VectorStateType>
        void set_cross_section(const VectorStateType &cs);

        template<typename VectorStateType>
        void update_cross_section(const VectorStateType &custom_x);

        //!\return the abscissa
        const VectorCoeffType abscissa()      const;

        //!\return the cross-section
        const VectorCoeffType cross_section() const;

        //!\return the cross-section on custom grid
        const VectorCoeffType &cross_section_on_custom_grid() const;
  };

  template<typename VectorCoeffType>
  inline
  CrossSection<VectorCoeffType>::CrossSection(const VectorCoeffType &x, const VectorCoeffType &y):
  _abscissa(x),
  _cross_section(y),
  {
    return;
  }

  template<typename VectorCoeffType>
  inline
  CrossSection<VectorCoeffType>::CrossSection(const CrossSection<VectorCoeffType> &x):
  _abscissa(rhs.abscissa()),
  _cross_section(rhs.cross_section()),
  {
    return;
  }

  template<typename VectorCoeffType>
  template<typename VectorStateType>
  inline
  void CrossSection<VectorCoeffType>::set_abscissa(const VectorStateType &x)
  {
      _abscissa = x;
  }

  template<typename VectorCoeffType>
  template<typename VectorStateType>
  inline
  void CrossSection<VectorCoeffType>::set_cross_section(const VectorStateType &cs)
  {
     _cross_section = cs;
  }

  template<typename VectorCoeffType>
  template<typename VectorStateType>
  inline
  void CrossSection<VectorCoeffType>::update_cross_section(const VectorStateType &custom_x)
  {
     _converter.y_on_custom_grid(_abscissa, _cross_section,
                                 custom_x, _cross_section_on_custom_grid);

     return;
  }

  template<typename VectorCoeffType>
  inline
  const VectorCoeffType CrossSection<VectorCoeffType>::abscissa() const
  {
     return _abscissa;
  }

  template<typename VectorCoeffType>
  inline
  const VectorCoeffType CrossSection<VectorCoeffType>::cross_section() const
  {
     return _cross_section;
  }

  template<typename VectorCoeffType>
  inline
  const VectorCoeffType &CrossSection<VectorCoeffType>::cross_section_on_custom_grid() const
  {
     return _cross_section_on_custom_grid;
  }

}

#endif
