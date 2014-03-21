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

#ifndef PLANET_PDF_MANAGEMENT_H
#define PLANET_PDF_MANAGEMENT_H

//Antioch
#include "antioch/antioch_asserts.h" //for error

//Planet
#include "planet/pdf_enum.h"
#include "planet/base_pdf.h"
#include "planet/norm_pdf.h"
#include "planet/nort_pdf.h"
#include "planet/logn_pdf.h"
#include "planet/logu_pdf.h"
#include "planet/unif_pdf.h"
#include "planet/diri_pdf.h"
#include "planet/diun_pdf.h"
#include "planet/dior_pdf.h"
#include "planet/diut_pdf.h"
#include "planet/dirg_pdf.h"

namespace Planet
{
  namespace ManagePDF
  {
     template <typename CoeffType>
     BasePdf<CoeffType> *  create_pdf_pointer(PDFName::PDFName type_pdf);

//

     template <typename CoeffType>
     inline
     BasePdf<CoeffType> *  create_pdf_pointer(PDFName::PDFName type_pdf)
     {

       BasePdf<CoeffType> *  pdf(NULL);
       switch(type_pdf)
       {
        case PDFName::Norm:
          pdf = new NormPdf<CoeffType>();
          break;

        case PDFName::NorT:
          pdf = new NorTPdf<CoeffType>();
          break;

        case PDFName::Unif:
          pdf = new UnifPdf<CoeffType>();
          break;

        case PDFName::LogN:
          pdf = new LogNPdf<CoeffType>();
          break;

        case PDFName::LogU:
          pdf = new LogUPdf<CoeffType>();
          break;

        case PDFName::Diri:
          pdf = new DiriPdf<CoeffType>();
          break;

        case PDFName::DiUn:
          pdf = new DiUnPdf<CoeffType>();
          break;

        case PDFName::DiUT:
          pdf = new DiUTPdf<CoeffType>();
          break;

        case PDFName::DiOr:
          pdf = new DiOrPdf<CoeffType>();
          break;

        case PDFName::DirG:
          pdf = new DirGPdf<CoeffType>();
          break;

        default: //WAT?
          antioch_error();
          break;
       }

       return pdf;
     }

  }//end namespace ManagePDF
}// end namespace Planet

#endif
