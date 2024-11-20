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


#include "ContactPolygon.h"
#include "ContactFixedSizeAllocator.h"
#include "allocators.h"
#include <cstddef>
#include <cstring>
#include <new>

ContactPolyVert::ContactPolyVert(ContactPolyVertType Type, 
                                 int ID, Real xi, Real eta)
{
  type    = Type;
  type_id = -1;
  id      = ID;
  shared  = NULL;
  next    = NULL;
  x       = xi;
  y       = eta;
}

void ContactPolyVert_SizeAllocator( ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactPolyVert),
                100,  // block size
                0);   // initial block size
  alloc.Set_Name( "ContactPolyVert allocator" );
}

ContactPolyVert* ContactPolyVert::new_ContactPolyVert(
                                               ContactFixedSizeAllocator& alloc,
                                               ContactPolyVertType Type, 
                                               int ID, Real xi, Real eta)
{
  return new (alloc.New_Frag())ContactPolyVert(Type, ID, xi, eta);
}

ContactPolyVert::~ContactPolyVert() {}

void ContactPolyVert::DetermineType(int face_type, Real tol)
{
  switch (type) {
  case PERIMETER:
    switch (face_type) {
    case ContactSearch::TRIFACEL3:
      if (x<(-1.0+tol) && y<(-1.0+tol)) {
        type    = PVERTEX;
        type_id = 0;
      } else if (x>(1.0-tol) && y<(-1.0+tol)) {
        type    = PVERTEX;
        type_id = 1;
      } else if (x>(1.0-tol) && y>(1.0-tol)) {
        type    = PVERTEX;
        type_id = 2;
      } else if (y<(-1.0+tol)) {
        type    = PEDGE;
        type_id = 0;
      } else if (x<(-1.0+tol)) {
        type    = PEDGE;
        type_id = 2;
      } else if (-tol<x+y && x+y<tol) {
        type    = PEDGE;
        type_id = 1;
      }
      break;
    case ContactSearch::QUADFACEL4:
      if (x<(-1.0+tol) && y<(-1.0+tol)) {
        type    = PVERTEX;
        type_id = 0;
      } else if (x>(1.0-tol) && y<(-1.0+tol)) {
        type    = PVERTEX;
        type_id = 1;
      } else if (x>(1.0-tol) && y>(1.0-tol)) {
        type    = PVERTEX;
        type_id = 2;
      } else if (x<(-1.0+tol) && y>(1.0-tol)) {
        type    = PVERTEX;
        type_id = 3;
      } else if (y<(-1.0+tol)) {
        type    = PEDGE;
        type_id = 0;
      } else if (x>(1.0-tol)) {
        type    = PEDGE;
        type_id = 1;
      } else if (y>(1.0-tol)) {
        type    = PEDGE;
        type_id = 2;
      } else if (x<(-1.0+tol)) {
        type    = PEDGE;
        type_id = 3;
      }
      break;
    default:
      POSTCONDITION(0);
      break;
    }
    break;
  case SLAVE:
    switch (face_type) {
    case ContactSearch::TRIFACEL3:
      if (x<(-1.0+tol) && y<(-1.0+tol)) {
        type_id = 0;
      } else if (x>(1.0-tol) && y<(-1.0+tol)) {
        type_id = 1;
      } else if (x<(-1.0+tol) && y>(1.0-tol)) {
        type_id = 2;
      }
      break;
    case ContactSearch::QUADFACEL4:
      if (x<(-1.0+tol) && y<(-1.0+tol)) {
        type_id = 0;
      } else if (x>(1.0-tol) && y<(-1.0+tol)) {
        type_id = 1;
      } else if (x>(1.0-tol) && y>(1.0-tol)) {
        type_id = 2;
      } else if (x<(-1.0+tol) && y>(1.0-tol)) {
        type_id = 3;
      }
      break;
    default:
      POSTCONDITION(0);
      break;
    }
    break;
  default:
    break;
  }
}

//==============================================================================

ContactPolyEdge::ContactPolyEdge(int ID)
{
  id    = ID;
  vert1 = NULL;
  vert2 = NULL;
  poly1 = NULL;
  poly2 = NULL;
}

void ContactPolyEdge_SizeAllocator( ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactPolyEdge),
                100,  // block size
                0);   // initial block size
  alloc.Set_Name( "ContactPolyEdge allocator" );
}

ContactPolyEdge* ContactPolyEdge::new_ContactPolyEdge(
                                               ContactFixedSizeAllocator& alloc,
                                               int ID)
{
  return new (alloc.New_Frag())ContactPolyEdge(ID);
}

ContactPolyEdge::~ContactPolyEdge() { }

//==============================================================================

ContactPoly::ContactPoly(int ID, int Nverts)
{
  // polys are closed, so nverts == nedges
  id     = ID;
  nverts = Nverts;
  nedges = Nverts;
  verts  = new ContactPolyVert* [Nverts];
  edges  = new ContactPolyEdge* [Nverts];
}

void ContactPoly_SizeAllocator( ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactPoly),
                100,  // block size
                0);   // initial block size
  alloc.Set_Name( "ContactPoly allocator" );
}

ContactPoly* ContactPoly::new_ContactPoly(ContactFixedSizeAllocator& alloc,
                                          int ID, int Nverts)
{
  return new (alloc.New_Frag())ContactPoly(ID, Nverts);
}

ContactPoly::~ContactPoly() 
{
  delete [] verts;
  delete [] edges;
}

void ContactPoly::EquivalanceVertices()
{
  int i0  = 0;
  int cnt = 1;
  for (int i=1; i<nverts; ++i) {
    if (verts[i]!=verts[i0]) {
      verts[++i0] = verts[i];
      cnt++;
    }
  }
  nverts = cnt;
  nedges = cnt;
}
