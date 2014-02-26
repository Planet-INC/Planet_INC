//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef PLANET_KINETICS_BRANCHING_STRUCTURE_H
#define PLANET_KINETICS_BRANCHING_STRUCTURE_H

//Antioch

//Planet
#include "planet/pdf_enum.h"

//C++
#include <vector>

namespace Planet
{

  template <typename CoeffType>
  class KineticsBranchingStructure
  {
       private:
        std::vector<KineticsBranchingStructure*> _children;
        std::vector<std::vector<CoeffType> >     _sample;
        KineticsBranchingStructure *             _mum;
        unsigned int                             _index_mommy;
        PDFName::PDFName                         _pdf;

       public:
         KineticsBranchingStructure();
         ~KineticsBranchingStructure();

        void add_child( KineticsBranchingStructure *child, unsigned int my_index);
        void set_mother(KineticsBranchingStructure *mum,   unsigned int index_mommy);
        void set_reaction(Antioch::Reaction *reaction);

        KineticStructure * mum();
  };

  template <typename CoeffType>
  inline
  void KineticsBranchingStructure<CoeffType>::set_mother(KineticsBranchingStructure *mum, unsigned int index_mommy)
  {
      _mum = mum;
      _index_mommy = index_mommy;
  }

  template <typename CoeffType>
  inline
  void KineticsBranchingStructure<CoeffType>::set_mother(Antioch::Reaction * reaction)
  {
      KineticsBranchingStructure * kinStruct = this;
      while(kinStruct->mum())kinStruct = kinStruct->mum();
     
      kinStruct->set_reaction(reaction);
  }

  template <typename CoeffType>
  inline
  void KineticsBranchingStructure<CoeffType>::add_child(KineticsBranchingStructure *child, unsigned int my_index)
  {
      _children.push_back(child);
      _children.back()->set_mother(this, my_index);
      _sample.resize(_children.size(),0.L);
  }

  template <typename CoeffType>
  inline
  KineticsBranchingStructure<CoeffType>::KineticsBranchingStructure():
      _mum(NULL),
      _index_mommy(0),
      _pdf(PDFName::PDFName::Norm)
  {
     return;
  }

  template <typename CoeffType>
  inline
  KineticsBranchingStructure<CoeffType>::~KineticsBranchingStructure()
  {
      for(unsigned int ic = 0; ic < _children.size(); ic++)
      {
         delete _children[ic];
      }
  }

  template <typename CoeffType>
  inline
  KineticStructure * KineticStructure<CoeffType>::mum()
  {
      return _mum;
  }

}

#endif
