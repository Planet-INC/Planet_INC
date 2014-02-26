//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef PLANET_REACTION_PDF_H
#define PLANET_REACTION_PDF_H

//Antioch
#include "antioch/reaction.h"

//Planet
#include "planet/pdf_enum.h"
#include "planet/kinetics_branching_structure.h"

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
