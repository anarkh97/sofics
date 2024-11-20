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


#ifndef ContactLineEdgeL2_h_
#define ContactLineEdgeL2_h_

#include "ContactEdge.h"

class ContactFixedSizeAllocator;
template<typename DataType> class ContactNode;

template<typename DataType>
class ContactLineEdgeL2 : public ContactEdge<DataType> {

 public:
  ContactLineEdgeL2( int block_index = -1, int index_in_block = -1 );
  static ContactLineEdgeL2<DataType>* new_ContactLineEdgeL2( ContactFixedSizeAllocator&,
						   int block_index=-1, 
						   int index_in_block = -1 );
  ~ContactLineEdgeL2();

 protected:
 private:
  ContactNode<DataType>* nodes[2];
};

#endif // ContactLineEdgeL2_h_
