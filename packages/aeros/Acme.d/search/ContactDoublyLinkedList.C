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


#include "ContactEntity.h"
#include "ContactEdge.h"
#include "ContactNode.h"
#include "ContactDoublyLinkedList.h"
#include <new>

using namespace std;

ContactTopologyDLL::ContactTopologyDLL() {
  current = data.begin();
  num_entities = 0;
}

ContactTopologyDLL::~ContactTopologyDLL() {
}

int ContactTopologyDLL::Append( ContactTopologyEntity<Real>* entity ) {
  data.push_back(entity);
  num_entities++;
  return data.size() - 1;
}

void ContactTopologyDLL::Clear() {
  data.clear();
  current = data.begin();
  num_entities = 0;
}

void ContactTopologyDLL::Display(ContactParOStream& postream)
{
  IteratorStart();
  while( ContactTopologyEntity<Real>* entity=IteratorForward() ){
    entity->Display(postream);
  }
}

