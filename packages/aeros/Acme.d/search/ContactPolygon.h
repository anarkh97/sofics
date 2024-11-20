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


#ifndef ContactPolygon_h_
#define ContactPolygon_h_

#include "Contact_Defines.h"
#include "ContactEntity.h"
#include "ContactSearch.h"

class ContactPolyVert;
class ContactPolyEdge;
class ContactPoly;
class ContactFixedSizeAllocator;

class ContactPolyVert {

  public:
    enum ContactPolyVertType{ UNKNOWN, INTERIOR, PERIMETER, 
                              PEDGE, PVERTEX, SLAVE};
    ContactPolyVert(ContactPolyVertType, int, Real, Real);
    static ContactPolyVert* new_ContactPolyVert( ContactFixedSizeAllocator&,
				                 ContactPolyVertType, 
                                                 int, Real, Real);
    ~ContactPolyVert();
    
    void   DetermineType(int, Real);
    inline void Type(ContactPolyVertType t) { type=t; };
    inline ContactPolyVertType Type() { return type; };
    inline void TypeID(int t) { type_id=t; };
    inline int  TypeID() { return type_id; };
    inline int  ID() { return id; };
    inline void ID(int ID) { id=ID; };
    inline Real X() { return x; };
    inline Real Y() { return y; };
    inline void Shared(ContactPolyVert* v) { shared=v; };
    inline ContactPolyVert* Shared() { return shared; };
    inline void Next(ContactPolyVert* v) { next=v; };
    inline ContactPolyVert* Next() { return next; };
  
  private:
  
    ContactPolyVert(ContactParOStream&);
    ContactPolyVert& operator=(ContactParOStream&);
    
    ContactPolyVertType type;
    int  type_id;
    int  id;
    Real x;
    Real y;
    ContactPolyVert *shared;
    ContactPolyVert *next;
};

class ContactPolyEdge {

  public:
    ContactPolyEdge(int);
    static ContactPolyEdge* new_ContactPolyEdge( ContactFixedSizeAllocator&,
                                                 int);
    ~ContactPolyEdge();
  
    ContactPolyEdge(ContactPolyEdge&);
    ContactPolyEdge& operator=(ContactPolyEdge&);
    
    inline int  ID() { return id; };
    inline void ID(int ID) { id=ID; };
    inline void Vertex1(ContactPolyVert* vert) { vert1=vert; };
    inline void Vertex2(ContactPolyVert* vert) { vert2=vert; };
    inline void Poly1(ContactPoly* poly) { poly1=poly; };
    inline void Poly2(ContactPoly* poly) { poly2=poly; };
    inline ContactPolyVert* Vertex1() { return vert1; };
    inline ContactPolyVert* Vertex2() { return vert2; };
    inline ContactPoly*     Poly1() { return poly1; };
    inline ContactPoly*     Poly2() { return poly2; };
    
  private:
    int id;
    ContactPolyVert* vert1;
    ContactPolyVert* vert2;
    ContactPoly*     poly1;
    ContactPoly*     poly2;
};

class ContactPoly {

  public:
    ContactPoly(int, int);
    static ContactPoly* new_ContactPoly( ContactFixedSizeAllocator&, int, int);
    ~ContactPoly();
  
    ContactPoly(ContactPoly&);
    ContactPoly& operator=(ContactPoly&);
    
    inline int  ID() { return id; };
    inline void ID(int ID) { id=ID; };
    inline int  NumVerts() { return nverts; };
    inline void NumVerts(int n) { nverts=n; };
    inline void Vertex(int index, ContactPolyVert* vertex) 
           { verts[index] = vertex; };
    inline ContactPolyVert* Vertex(int index) { return verts[index]; };
           
    inline int  NumEdges() { return nedges; };
    inline void NumEdges(int n) { nedges=n; };
    inline void Edge(int index, ContactPolyEdge* edge) 
           { edges[index] = edge; };
    inline ContactPolyEdge* Edge(int index) { return edges[index]; };
    void   EquivalanceVertices();
    
  private:
    int id;
    int nverts;
    int nedges;
    ContactPolyVert** verts;
    ContactPolyEdge** edges;
};

#endif
