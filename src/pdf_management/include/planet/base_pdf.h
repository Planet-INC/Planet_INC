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

#ifndef PLANET_BASE_PDF_H
#define PLANET_BASE_PDF_H

//Planet
#include "planet/pdf_enum.h"

//C++
#include <iostream>
#include <vector>

namespace Planet
{

//forward declaration
template <typename CoeffType>
class NormPdf;

template <typename CoeffType>
class NorTPdf;

template <typename CoeffType>
class LogNPdf;

template <typename CoeffType>
class LogUPdf;

template <typename CoeffType>
class UnifPdf;

template <typename CoeffType>
class DiriPdf;

template <typename CoeffType>
class DirGPdf;

template <typename CoeffType>
class DirWPdf;

template <typename CoeffType>
class DiUnPdf;

template <typename CoeffType>
class DiUTPdf;

template <typename CoeffType>
class DiOrPdf;


   /*! \class BasePdf
    *
    * Base class for pdf object.
    *
    */
  template <typename CoeffType>
  class BasePdf
  {
      public:

      BasePdf(PDFName::PDFName pdf);

      template <typename StateType>
      BasePdf(const BasePdf<StateType> &rhs);

      virtual ~BasePdf();

      PDFName::PDFName pdf() const;

      const CoeffType value(unsigned int ip = 0) const;

      template <typename StateType>
      void set_parameters(const std::vector<StateType> &pars);

      template <typename StateType>
      void get_parameters(std::vector<StateType> &pars) const;

      void print(std::ostream &out = std::cout) const;

      friend std::ostream& operator<<(std::ostream &out, const BasePdf & bpdf)
      {
        bpdf.print(out);
        return out;
      }

      protected:
      const PDFName::PDFName _my_pdf;

      private:
        BasePdf();

  };

  template <typename CoeffType>
  inline
  BasePdf<CoeffType>::BasePdf(PDFName::PDFName pdf):
       _my_pdf(pdf)
  {
     return;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  BasePdf<CoeffType>::BasePdf(const BasePdf<StateType> &rhs):
    _my_pdf(rhs.pdf())
  {
    return;
  }

  template <typename CoeffType>
  inline
  BasePdf<CoeffType>::~BasePdf()
  {
     return;
  }

  template <typename CoeffType>
  inline
  PDFName::PDFName BasePdf<CoeffType>::pdf() const
  {
     return _my_pdf;
  }
    

  template <typename CoeffType>
  inline
  const CoeffType BasePdf<CoeffType>::value(unsigned int ip) const
  {
       switch(_my_pdf)
       {
        case PDFName::Norm:
          return static_cast<const NormPdf<CoeffType>*>(this)->value(ip);
          break;

        case PDFName::NorT:
          return static_cast<const NorTPdf<CoeffType>* >(this)->value(ip);
          break;

        case PDFName::Unif:
          return static_cast<const UnifPdf<CoeffType>* >(this)->value(ip);
          break;

        case PDFName::LogN:
          return static_cast<const LogNPdf<CoeffType>* >(this)->value(ip);
          break;

        case PDFName::LogU:
          return static_cast<const LogUPdf<CoeffType>* >(this)->value(ip);
          break;

        case PDFName::Diri:
          return static_cast<const DiriPdf<CoeffType>* >(this)->value(ip);
          break;

        case PDFName::DiUn:
          return static_cast<const DiUnPdf<CoeffType>* >(this)->value(ip);
          break;

        case PDFName::DiUT:
          return static_cast<const DiUTPdf<CoeffType>* >(this)->value(ip);
          break;

        case PDFName::DiOr:
          return static_cast<const DiOrPdf<CoeffType>* >(this)->value(ip);
          break;

        case PDFName::DirG:
          return static_cast<const DirGPdf<CoeffType>* >(this)->value(ip);
          break;

        default: //WAT?
          antioch_error();
          break;
       }

    // Dummy
    return 0.;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void BasePdf<CoeffType>::set_parameters(const std::vector<StateType> &pars)
  {
       switch(_my_pdf)
       {
        case PDFName::Norm:
          (static_cast<NormPdf<CoeffType>*>(this))->set_parameters(pars);
          break;

        case PDFName::NorT:
          (static_cast<NorTPdf<CoeffType>* >(this))->set_parameters(pars);
          break;

        case PDFName::Unif:
          (static_cast<UnifPdf<CoeffType>* >(this))->set_parameters(pars);
          break;

        case PDFName::LogN:
          (static_cast<LogNPdf<CoeffType>* >(this))->set_parameters(pars);
          break;

        case PDFName::LogU:
          (static_cast<LogUPdf<CoeffType>* >(this))->set_parameters(pars);
          break;

        case PDFName::Diri:
          (static_cast<DiriPdf<CoeffType>* >(this))->set_parameters(pars);
          break;

        case PDFName::DiUn:
          (static_cast<DiUnPdf<CoeffType>* >(this))->set_parameters(pars);
          break;

        case PDFName::DiUT:
          (static_cast<DiUTPdf<CoeffType>* >(this))->set_parameters(pars);
          break;

        case PDFName::DiOr:
          (static_cast<DiOrPdf<CoeffType>* >(this))->set_parameters(pars);
          break;

        case PDFName::DirG:
          (static_cast<DirGPdf<CoeffType>* >(this))->set_parameters(pars);
          break;

        default: //WAT?
          antioch_error();
          break;
       }
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void BasePdf<CoeffType>::get_parameters(std::vector<StateType> &pars) const
  {
       switch(_my_pdf)
       {
        case PDFName::Norm:
          static_cast<NormPdf<CoeffType>*>(this)->get_parameters(pars);
          break;

        case PDFName::NorT:
          static_cast<NorTPdf<CoeffType>* >(this)->get_parameters(pars);
          break;

        case PDFName::Unif:
          static_cast<UnifPdf<CoeffType>* >(this)->get_parameters(pars);
          break;

        case PDFName::LogN:
          static_cast<LogNPdf<CoeffType>* >(this)->get_parameters(pars);
          break;

        case PDFName::LogU:
          static_cast<LogUPdf<CoeffType>* >(this)->get_parameters(pars);
          break;

        case PDFName::Diri:
          static_cast<DiriPdf<CoeffType>* >(this)->get_parameters(pars);
          break;

        case PDFName::DiUn:
          static_cast<DiUnPdf<CoeffType>* >(this)->get_parameters(pars);
          break;

        case PDFName::DiUT:
          static_cast<DiUTPdf<CoeffType>* >(this)->get_parameters(pars);
          break;

        case PDFName::DiOr:
          static_cast<DiOrPdf<CoeffType>* >(this)->get_parameters(pars);
          break;

        case PDFName::DirG:
          static_cast<DirGPdf<CoeffType>* >(this)->get_parameters(pars);
          break;

        default: //WAT?
          antioch_error();
          break;
       }
  }

  template <typename CoeffType>
  inline
  void BasePdf<CoeffType>::print(std::ostream &out) const
  {
       switch(_my_pdf)
       {
        case PDFName::Norm:
          static_cast<const NormPdf<CoeffType>*>(this)->print(out);
          break;

        case PDFName::NorT:
          static_cast<const NorTPdf<CoeffType>* >(this)->print(out);
          break;

        case PDFName::Unif:
          static_cast<const UnifPdf<CoeffType>* >(this)->print(out);
          break;

        case PDFName::LogN:
          static_cast<const LogNPdf<CoeffType>* >(this)->print(out);
          break;

        case PDFName::LogU:
          static_cast<const LogUPdf<CoeffType>* >(this)->print(out);
          break;

        case PDFName::Diri:
          static_cast<const DiriPdf<CoeffType>* >(this)->print(out);
          break;

        case PDFName::DiUn:
          static_cast<const DiUnPdf<CoeffType>* >(this)->print(out);
          break;

        case PDFName::DiUT:
          static_cast<const DiUTPdf<CoeffType>* >(this)->print(out);
          break;

        case PDFName::DiOr:
          static_cast<const DiOrPdf<CoeffType>* >(this)->print(out);
          break;

        case PDFName::DirG:
          static_cast<const DirGPdf<CoeffType>* >(this)->print(out);
          break;

        default: //WAT?
          antioch_error();
          break;
       }
  }

}

#endif
