//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef PLANET_KINETICS_STRUCTURE_H
#define PLANET_KINETICS_STRUCTURE_H

//Planet
#include "planet/pdf_enum.h"

//C++
#include <vector>

namespace Planet
{

  class KineticsStructure
  {
       private:
        std::vector<KineticsStructure*> _children;
        KineticsStructure* _mum;
        PDFName::PDFName _pdf;

       public:
         KineticsStructure();
         ~KineticsStructure();

        add_child(KineticsStructure *child);
        set_mother(KineticsStructure *mum);
  };

  inline
  void KineticsStructure::set_mother(KineticsStructure *mum)
  {
      _mum = mum;
  }

  inline
  void KineticsStructure::add_child(KineticsStructure *child)
  {
      _children.push_back(child);
      _children.back()->set_mother(this);
  }

  inline
  KineticsStructure::KineticsStructure():
      _mum(NULL),
      _pdf(PDFName::PDFName::Norm)
  {
     return;
  }

  inline
  KineticsStructure::~KineticsStructure()
  {
      for(unsigned int ic = 0; ic < _children.size(); ic++)
      {
         delete _children[ic];
      }
  }

}

#endif
