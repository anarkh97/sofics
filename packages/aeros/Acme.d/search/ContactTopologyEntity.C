// This file is part of a modified version of ACME: Algorithms for Contact in
// a Multiphysics Environment, derived from version 2.7f
//
// Copyright (C) 2007 Sandia Corporation
// Copyright (C) 2011 Stanford University
//
// ACME is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ACME is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with ACME.  If not, see <http://www.gnu.org/licenses/>.


#ifndef ContactTopologyEntity_C_
#define ContactTopologyEntity_C_

#include "ContactTopologyEntity.h"
#include "ContactParOStream.h"
#include "ContactNode.h"
#include "ContactShellNode.h"
#include "ContactEdge.h"
#include "ContactFace.h"
#include "ContactElement.h"
#include "contact_assert.h"
#include <cstddef>
#include <cstring>

template<typename DataType>
ContactTopologyEntity<DataType>::ContactTopologyEntity(DataType *data_array_, 
                                             const ContactType base_type_)
  : temp_tag(0),temp_tag1(0),in_proximity(0),
    host_array_index(-1),proc_array_index(-1),enf_array_index(-1),
    primary_proc_array_index(-1),owner_proc_array_index(-1),block_id(-1),
    base_type(base_type_),data_array(data_array_), 
    global_id(0,-1),shared(false),owner(-1),secondary_owner(-1)
{ }

template<typename DataType>
ContactTopologyEntity<DataType>::ContactTopologyEntity( int Block_ID, 
			                      int host_index_in_block,
                                              DataType *data_array_, 
                                              const ContactType base_type_)
  : temp_tag(0),temp_tag1(0),in_proximity(0),
    host_array_index(host_index_in_block), proc_array_index(-1),
    enf_array_index(-1),primary_proc_array_index(-1), 
    owner_proc_array_index(-1),block_id(Block_ID),
    base_type(base_type_),
    data_array(data_array_),global_id(0,-1),
    shared(false), owner(-1),secondary_owner(-1)
{ }

template<typename DataType>
ContactTopologyEntity<DataType>::~ContactTopologyEntity()
{}

template<typename DataType>
void 
ContactTopologyEntity<DataType>::Display(ContactParOStream& postream)
{
  postream<<"ContactEntity: "<<global_id<<"\n";
  postream<<"               entity type:             ";
  switch (Base_Type()) {
  case CT_NODE:
    postream<<"CT_NODE\n";
    break;
  case CT_SHELL_NODE:
    postream<<"CT_SHELL_NODE\n";
    break;
  case CT_EDGE:
    postream<<"CT_EDGE\n";
    break;
  case CT_FACE:
    {
      postream<<"CT_FACE\n";
      postream<<"               face type:               ";
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(this);
      switch (face->FaceType()) {
      case ContactSearch::QUADFACEL4:
        postream<<"QUADFACEL4\n";
        break;
      case ContactSearch::QUADFACEQ8:
        postream<<"QUADFACEQ8\n";
        break;
      case ContactSearch::QUADFACEQ9:
        postream<<"QUADFACEQ9\n";
        break;
      case ContactSearch::TRIFACEL3:
        postream<<"TRIFACEL3\n";
        break;
      case ContactSearch::TRIFACEQ6:
        postream<<"TRIFACEQ6\n";
        break;
      case ContactSearch::SHELLQUADFACEL4:
        postream<<"SHELLQUADFACEL4\n";
        break;
      case ContactSearch::SHELLTRIFACEL3:
        postream<<"SHELLTRIFACEL3\n";
        break;
      case ContactSearch::LINEFACEL2:
        postream<<"LINEFACEL2\n";
        break;
      case ContactSearch::LINEFACEQ3:
        postream<<"LINEFACEQ3\n";
        break;
      default:
        postream<<"UNKNOWN\n";
        break;
      }
    }
    break;
  case CT_ELEM:
    postream<<"CT_ELEM\n";
    break;
  case CT_ELEMENT: 
    {
      postream<<"CT_ELEMENT\n";
      postream<<"               element type:            ";
      ContactElement* element = static_cast<ContactElement*>(this);
      switch (element->ElementType()) {
      case ContactSearch::CARTESIANHEXELEMENTL8:
        postream<<"CARTESIANHEXELEMENTL8\n";
        break;
      case ContactSearch::HEXELEMENTL8:
        postream<<"HEXELEMENTL8\n";
        break;
      default:
        postream<<"UNKNOWN\n";
        break;
      }
    }
    break;
  case CT_ANALYTIC_SURFACE:
    postream<<"CT_ANALYTIC_SURFACE\n";
    break;
  case CT_UNKNOWN:
    postream<<"CT_UNKNOWN\n";
    break;
  default:
    POSTCONDITION(0);
    break;
  }
  switch (Base_Type()) {
  case CT_NODE:
  case CT_SHELL_NODE:
    {
      postream<<"               node type:               ";
      ContactNode<Real>* node = static_cast<ContactNode<Real>*>(this);
      switch (node->NodeType()) {
      case ContactSearch::NODE:
        postream<<"NODE\n";
        break;
      case ContactSearch::POINT:
        postream<<"POINT\n";
        break;
      default:
        postream<<"UNKNOWN\n";
        break;
      }
      postream<<"               physical type:           ";
      switch (node->Physical_Type()) {
      case ContactNode<Real>::CONTINUUM_NODE:
        postream<<"CONTINUUM_NODE\n";
        break;
      case ContactNode<Real>::MIXED_NODE:
        postream<<"MIXED_NODE\n";
        break;
      case ContactNode<Real>::SHELL_NODE:
        postream<<"SHELL_NODE\n";
        break;
      default:
        postream<<"UNKNOWN\n";
        break;
      }
    }
    break;
  default:
    break;
  }
  postream<<"               block_id:                "<<block_id<<"\n";
  postream<<"               host_block_array_index:  "<<host_array_index<<"\n";
  postream<<"               host_global_array_index: "<<host_global_array_index<<"\n";
  postream<<"               proc_array_index:        "<<proc_array_index<<"\n";
  postream<<"               host_gid[0]:             "<<global_id.HiInt()<<"\n";
  postream<<"               host_gid[1]:             "<<global_id.LoInt()<<"\n";
  postream<<"               ownership:               ";
  switch (ownership) {
  case OWNED:
    postream<<"OWNED\n";
    break;
  default:  
    postream<<"NOT_OWNED\n";
    break;
  }
  postream<<"               shared:                  ";
  if (shared) {
    postream<<"SHARED\n";
  } else { 
    postream<<"NOT_SHARED\n";
  }
  postream<<"               primary_owner:           "<<owner<<"\n";
  postream<<"               secondary_owner:         "<<secondary_owner<<"\n";
  postream<<"               owner_array_index:       "<<owner_proc_array_index<<"\n";
  switch (Base_Type()) {
  case CT_NODE:
  case CT_SHELL_NODE: {
    ContactNode<Real>* node = static_cast<ContactNode<Real>*>(this);
    postream<<"               exodus_id:               "<<node->Exodus_ID()<<"\n";
    if (Base_Type()==CT_SHELL_NODE) {
      ContactShellNode* snode = static_cast<ContactShellNode*>(this);
      postream<<"               shell_node_base_id:      "<<snode->Shell_Node_Base_ID()<<"\n";
    }
    postream<<"               entity key:              "
            <<node->Entity_Key()<<"\n";
    postream<<"               number NFI:              "
            <<node->Number_NodeFace_Interactions()<<"\n";
    postream<<"               number NSI:              "
            <<node->Number_NodeSurface_Interactions()<<"\n";
    } break;
  case CT_EDGE: {
    postream<<"               node connectivity:      ";
    ContactEdge<Real>* edge = static_cast<ContactEdge<Real>*>(this);
    for (int i=0; i<edge->Nodes_Per_Edge(); ++i) {
      postream<<" "<<edge->Node(i)->Global_ID();
    }
    postream<<"\n";
    postream<<"               face connectivity:      ";
    if (shared) {
       postream<<" "<<edge->Face(0)->Global_ID();
    } else {
      for (int i=0; i<edge->Number_Face_Connections(); ++i) {
        postream<<" "<<edge->Face(i)->Global_ID();
      }
    }
    postream<<"\n";
    } break;
  case CT_FACE: {
    postream<<"               node connectivity:      ";
    ContactFace<Real>* face = static_cast<ContactFace<Real>*>(this);
    for (int i=0; i<face->Nodes_Per_Face(); ++i) {
      postream<<" "<<face->Node(i)->Global_ID();
    }
    postream<<"\n";
    postream<<"               face connectivity:      ";
    if (face->NumberOfNeighbors()>0) {
      for (int i=0; i<face->Edges_Per_Face(); ++i) {
        if (face->NeighborInfo()[i].owner>=0) {
          ContactHostGlobalID GID( face->NeighborInfo()[i].host_gid[0], 
                                   face->NeighborInfo()[i].host_gid[1] );
          postream<<" "<<i<<":"<<GID;
        }
      }
    }
    postream<<"\n";
    } break;
  case CT_ELEMENT: {
    postream<<"               connectivity:           ";
    ContactElement* element = static_cast<ContactElement*>(this);
    for (int i=0; i<element->Nodes_Per_Element(); ++i) {
      postream<<" "<<element->Node(i)->Global_ID();
    }
    postream<<"\n";
    } break;
  default:
    POSTCONDITION(0);
    break;
  }
}

#endif  // #define ContactTopologyEntity_C_
