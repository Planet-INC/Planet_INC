//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef PLANET_REACTION_PDF_H
#define PLANET_REACTION_PDF_H

//QUESO

//Planet
#include "planet/pdf_enum.h"

//C++
#include <vector>

namespace Planet
{

  template <typename CoeffType>
  class ReactionPDF
  {
       private:
        KineticsBranchingStructure<CoeffType> _branching_ratios;
        std::vector<CoeffType>                _sample;
        PDFName::PDFName                      _pdf;
        std::vector<Antioch::Reaction &>      _reaction;

       public:
         ReactionPDF();
         ~ReactionPDF();

  };


  template <typename CoeffType>
  inline
  ReactionPDF<CoeffType>::ReactionPDF()
  {
     return;
  }

  template <typename CoeffType>
  inline
  ReactionPDF<CoeffType>::~ReactionPDF()
  {
     return;
  }

}

#endif
