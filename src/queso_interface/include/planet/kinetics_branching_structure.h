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

#ifndef PLANET_KINETICS_BRANCHING_STRUCTURE_H
#define PLANET_KINETICS_BRANCHING_STRUCTURE_H

//Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/kinetics_parsing.h"

//Planet
#include "planet/pdf_management.h"
#include "planet/branching_ratio_node.h"

//C++
#include <vector>
#include <string>
#include <map>

namespace Planet
{

  /*! \class KineticsBranchingStructure
   * 
   * This class stores a whole reaction, and
   * provides the reaction in the partial rate
   * constants form.
   * This is where the pdf informations are kept, the
   * structure is a probabilistic tree, with nodes and
   * pathes for each channels.
   * 
   * 
   */
  template <typename CoeffType>
  class KineticsBranchingStructure
  {
       private:
        std::vector<PDFName::PDFName>                 _pdf_k;
        std::vector<BasePdf<CoeffType>*>              _pdf_k_object;

        std::vector<BranchingRatioNode<CoeffType> >   _nodes;
        std::vector<std::vector<std::string> >        _track_br_nodes;    // nodes of channel
        std::vector<std::vector<unsigned int> >       _track_br_in_nodes; // parameters ib channel in nodes
        std::map<std::string,PDFName::PDFName>        _pdf_map;
        std::map<std::string,unsigned int>            _nodes_map;

       public:
        KineticsBranchingStructure();
        ~KineticsBranchingStructure();

        //! pdf map
        const std::map<std::string,PDFName::PDFName> pdf_map() const;

        //! pdf map
        const std::map<std::string,unsigned int> nodes_map()    const;

        //! adds already defined node
        template <typename StateType>
        void add_node(const BranchingRatioNode<StateType> &node);

        //! adds new empty node, needs a name
        void create_node(const std::string &id);

        //! modifiable reference to node
        BranchingRatioNode<CoeffType> & node(const std::string &id);

        //!sets the pdf for k
        void set_k_pdf(const std::vector<PDFName::PDFName> & pdf);

        //!pdf object for k
        std::vector<BasePdf<CoeffType> *> pdf_k_object();

        //!set the number of channels
        void set_n_channels(unsigned int nbr);

        //! add a node and channel in node to a br path
        void add_to_br_path(unsigned int ibr, const std::string &name_node, unsigned int node_channel);

        Antioch::KineticsType<CoeffType>* create_rate_constant(unsigned int ibr, Antioch::KineticsModel::KineticsModel kineticsModel);

  };

  template <typename CoeffType>
  inline
  KineticsBranchingStructure<CoeffType>::KineticsBranchingStructure()
  {
     _pdf_map["Norm"] = PDFName::Norm;
     _pdf_map["NorT"] = PDFName::NorT;
     _pdf_map["Unif"] = PDFName::Unif;
     _pdf_map["LogN"] = PDFName::LogN;
     _pdf_map["LogU"] = PDFName::LogU;
     _pdf_map["Diri"] = PDFName::Diri;
     _pdf_map["DiUn"] = PDFName::DiUn;
     _pdf_map["DiUT"] = PDFName::DiUT;
     _pdf_map["DirG"] = PDFName::DirG;
     _pdf_map["DiOr"] = PDFName::DiOr;
     return;
  }

  template <typename CoeffType>
  inline
  KineticsBranchingStructure<CoeffType>::~KineticsBranchingStructure()
  {
    for(unsigned int ip = 0; ip < _pdf_k_object.size(); ip++)
    {
      ManagePDF::delete_pdf_pointer(_pdf_k_object[ip],_pdf_k[ip]);
    }
    return;
  }

  template <typename CoeffType>
  inline
  const std::map<std::string,PDFName::PDFName> KineticsBranchingStructure<CoeffType>::pdf_map() const
  {
      return _pdf_map;
  }

  template <typename CoeffType>
  template <typename StateType>
  inline
  void KineticsBranchingStructure<CoeffType>::add_node(const BranchingRatioNode<StateType> &node)
  {
      _nodes.push_back(node);
      _nodes_map[node.id()] = _nodes.size() - 1;
  }

      
  template <typename CoeffType>
  inline
  void KineticsBranchingStructure<CoeffType>::create_node(const std::string &id)
  {
    BranchingRatioNode<CoeffType> node;
    node.set_id(id);
    this->add_node(node);
  }

  template <typename CoeffType>
  inline
  BranchingRatioNode<CoeffType> & KineticsBranchingStructure<CoeffType>::node(const std::string &id)
  {
      return _nodes[_nodes_map.at(id)];
  }

  template <typename CoeffType>
  inline
  void KineticsBranchingStructure<CoeffType>::set_k_pdf(const std::vector<PDFName::PDFName> &pdf)
  {
      if(!_pdf_k_object.empty())
      {
         for(unsigned int i = 0; i < _pdf_k_object.size(); i++)
         {
           ManagePDF::delete_pdf_pointer(_pdf_k_object[i],_pdf_k[i]);
         }
      }
      _pdf_k.resize(pdf.size());
      _pdf_k_object.resize(pdf.size(),NULL);
      for(unsigned int ip = 0; ip < _pdf_k_object.size(); ip++)
      {
           _pdf_k[ip]        = pdf[ip];
           _pdf_k_object[ip] = ManagePDF::create_pdf_pointer<CoeffType>(_pdf_k[ip]);
      }
  }

  template <typename CoeffType>
  inline
  std::vector<BasePdf<CoeffType> *> KineticsBranchingStructure<CoeffType>::pdf_k_object()
  {
     return _pdf_k_object;
  }

  template <typename CoeffType>
  inline
  void KineticsBranchingStructure<CoeffType>::set_n_channels(unsigned int nbr)
  {
     _track_br_nodes.resize(nbr);
     _track_br_in_nodes.resize(nbr);
  }

  template <typename CoeffType>
  inline
  void KineticsBranchingStructure<CoeffType>::add_to_br_path(unsigned int ibr, const std::string &name_node, unsigned int node_channel)
  {
     antioch_assert_less(ibr,_track_br_nodes.size());
     _track_br_nodes[ibr].push_back(name_node);
     _track_br_in_nodes[ibr].push_back(node_channel);
  }

  template <typename CoeffType>
  inline
  Antioch::KineticsType<CoeffType>* KineticsBranchingStructure<CoeffType>::create_rate_constant(unsigned int ibr, 
                                                                                             Antioch::KineticsModel::KineticsModel kineticsModel)
  {
     antioch_assert(!_pdf_k_object.empty());

     std::vector<CoeffType> data;
     data.resize(_pdf_k_object.size());
     for(unsigned int ik = 0; ik < _pdf_k_object.size(); ik++)
     {
       data[ik] = _pdf_k_object[ik]->value();// A
     }

     // Cf, branching ratio calculation
     for(unsigned int inode = 0; inode < _track_br_nodes[ibr].size(); inode++)
     {
std::cout << "node path " << inode << std::endl;
std::cout << "node name is " << _track_br_nodes[ibr][inode] << std::endl;
std::cout << "node index " << _nodes_map.at(_track_br_nodes[ibr][inode]) << std::endl;
std::cout << "parameter #" <<  _track_br_in_nodes[ibr][inode] << std::endl;
       data[0] *= _nodes[_nodes_map.at(_track_br_nodes[ibr][inode])].value(_track_br_in_nodes[ibr][inode]);
     }
     
std::cout << "huh?" << std::endl;
     return Antioch::build_rate<CoeffType>(data,kineticsModel);
  }

  template <typename CoeffType>
  inline
  const std::map<std::string,unsigned int> KineticsBranchingStructure<CoeffType>::nodes_map()    const
  {
     return _nodes_map;
  }

}

#endif
