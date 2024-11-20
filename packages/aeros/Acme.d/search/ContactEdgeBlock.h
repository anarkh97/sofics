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


#ifndef ContactEdgeBlock_h_
#define ContactEdgeBlock_h_

#include "ContactEdge.h"
#include "ContactBlockEntityList.h"

class ContactEdgeBlock {

 public:

  ContactEdgeBlock( ContactSearch::ContactEdge_Type, int, int, int, int&,
		    ContactTopology* );
  ContactEdgeBlock( ContactSearch::ContactEdge_Type, int,
		    ContactTopology* );
  ~ContactEdgeBlock();
  
  void Delete_Edge_List( );
  
  void Delete_Edges( );

  inline ContactSearch::ContactEdge_Type Type() { return type;};
  inline int Entity_Key() { return entity_key; };
  inline int ID() { return id; };
  inline int Number_of_Edges() { return number_of_edges; };
  ContactBlockEntityList* EdgeList() { return edge_list; };

  void Insert_Edge( ContactEdge<Real>* );
  void Delete_Edge( ContactEdge<Real>* );
#ifndef CONTACT_NO_MPI
  void Insert_Edge( char* );
#endif

 private:
  
  ContactEdgeBlock(ContactEdgeBlock&);
  ContactEdgeBlock& operator=(ContactEdgeBlock&);
  
  ContactTopology* topology;

  int number_of_edges;
  int num_edges_added;
  ContactSearch::ContactEdge_Type type;
  int entity_key;
  int id;

  ContactBlockEntityList* edge_list;

};

#endif //_ContactEdgeBlock_h_
