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

#ifndef PLANET_BRANCHING_RATIO_NODE_H
#define PLANET_BRANCHING_RATIO_NODE_H

//QUESO

//Planet
#include "planet/pdf_management.h"

//C++
#include <string>

namespace Planet
{

  /*! \class BranchingRatioNode
   * 
   * This class is a node in a reaction
   * probabilistic tree. It stores the
   * pdf of the node.
   * 
   */ 
  template <typename CoeffType>
  class BranchingRatioNode
  {
     private:
        PDFName::PDFName     _pdf;
        std::string          _id;
        unsigned int         _n_channels;
        BasePdf<CoeffType> * _pdf_object;

     public:
        BranchingRatioNode();
        template <typename StateType>
        BranchingRatioNode(const BranchingRatioNode<StateType> &rhs);
        ~BranchingRatioNode();

//pdf methods
        void set_pdf(PDFName::PDFName pdf);
        template <typename StateType>
        void set_pdf_parameters(const std::vector<StateType> &pars);

        void set_id(const std::string &id);
        void set_n_channels(unsigned int nchannels);

        template <typename StateType>
        void set_pdf_object(const BasePdf<StateType> * pdf_object);

        PDFName::PDFName pdf()                      const;
        const std::string id()                      const;
        unsigned int n_channels()                   const;
        const BasePdf<CoeffType> * pdf_object_ptr() const;

        const CoeffType value(unsigned int ip) const;

        void print(std::ostream &out = std::cout) const;

        friend std::ostream& operator<<(std::ostream &out, const BranchingRatioNode &brn)
        {
           brn.print(out);
           return out;
        }

        template <typename StateType>
        BranchingRatioNode<CoeffType> &operator=(const BranchingRatioNode<StateType> &rhs);

  };

  template <typename CoeffType>
  inline
  BranchingRatioNode<CoeffType>::BranchingRatioNode():
     _pdf(PDFName::Norm),
     _id("undefined node"),
     _n_channels(0),
     _pdf_object(NULL)
  {
     return;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  BranchingRatioNode<CoeffType>::BranchingRatioNode(const BranchingRatioNode<StateType> &rhs)
  {
     *this = rhs;
  }

  
  template <typename CoeffType>
  template <typename StateType>
  inline
  BranchingRatioNode<CoeffType> &BranchingRatioNode<CoeffType>::operator=(const BranchingRatioNode<StateType> &rhs)
  {
     if(this != &rhs)
     {
       _pdf = rhs.pdf();
       _id  = rhs.id();
       _n_channels = rhs.n_channels();
       if(_pdf_object)delete _pdf_object;
       _pdf_object = NULL;
       _pdf_object = ManagePDF::create_pdf_pointer<CoeffType>(_pdf);
       std::vector<CoeffType> pdf_pars;
       rhs.pdf_object_ptr()->get_parameters(pdf_pars);
       _pdf_object->set_parameters(pdf_pars);
     }

     return *this;
  }

  template <typename CoeffType>
  inline
  BranchingRatioNode<CoeffType>::~BranchingRatioNode()
  {
     if(_pdf_object)delete _pdf_object;
     _pdf_object = NULL;
     return;
  }

  template <typename CoeffType>
  inline
  void BranchingRatioNode<CoeffType>::set_pdf(PDFName::PDFName pdf)
  {
     if(_pdf_object)delete _pdf_object;

     _pdf = pdf;
     _pdf_object = ManagePDF::create_pdf_pointer<CoeffType>(_pdf);
  }

  template <typename CoeffType>
  inline
  void BranchingRatioNode<CoeffType>::set_id(const std::string &id)
  {
     _id = id;
  }

  template <typename CoeffType>
  inline
  void BranchingRatioNode<CoeffType>::set_n_channels(unsigned int nchannels)
  {
     _n_channels = nchannels;
  }

  template <typename CoeffType>
  inline
  PDFName::PDFName BranchingRatioNode<CoeffType>::pdf() const
  {
     return _pdf;
  }

  template <typename CoeffType>
  inline
  const std::string BranchingRatioNode<CoeffType>::id() const
  {
     return _id;
  }

  template <typename CoeffType>
  inline
  unsigned int BranchingRatioNode<CoeffType>::n_channels() const
  {
     return _n_channels;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void BranchingRatioNode<CoeffType>::set_pdf_object(const BasePdf<StateType> * pdf_object)
  {
      _pdf_object = pdf_object;
  }

  template <typename CoeffType>
  inline
  const CoeffType BranchingRatioNode<CoeffType>::value(unsigned int ip) const
  {
     antioch_assert(_pdf_object);
     return _pdf_object->value(ip);
  }

  template <typename CoeffType>
  inline
  const BasePdf<CoeffType> * BranchingRatioNode<CoeffType>::pdf_object_ptr() const
  {
     return _pdf_object;
  }

  template <typename CoeffType>
  inline
  void BranchingRatioNode<CoeffType>::print(std::ostream &out) const
  {
     out << "Node " << _id << std::endl;
     (_pdf_object)?out << "Distribution " << *_pdf_object:
                   out << "No distribution associated";
     out << std::endl;
     out << _n_channels << " channels" << std::endl;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void BranchingRatioNode<CoeffType>::set_pdf_parameters(const std::vector<StateType> &pars)
  {
      antioch_assert(_pdf_object);
      _pdf_object->set_parameters(pars);
  }
}

#endif
