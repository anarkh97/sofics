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


#ifndef ContactLineEdgeQ3_h_
#define ContactLineEdgeQ3_h_

#include "ContactEdge.h"

class ContactFixedSizeAllocator;
template<typename DataType> class ContactNode;

class ContactLineEdgeQ3 : public ContactEdge<Real> {

 public:
  ContactLineEdgeQ3( int block_index=-1, int host_index_in_block=-1 );
  static ContactLineEdgeQ3* new_ContactLineEdgeQ3( ContactFixedSizeAllocator&,
						   int block_index=-1, 
						   int host_index_in_block = -1 );
  ~ContactLineEdgeQ3();

 protected:
 private:
  ContactNode<Real>* nodes[3];
};

#endif // ContactLineEdgeQ3_h_
