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


#include <algorithm>

#include "allocators.h"
#include "ContactErrors.h"
#include "ContactFaceBlock.h"
#include "ContactFaceFaceInteraction.h"
#include "ContactFaceCoverageInteraction.h"
#include "ContactSearch.h"
#include "ContactTopology.h"
#include "ContactParOStream.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactPolygon.h"
#include "ContactTriFaceL3.h"
#include "CString.h"

#ifndef CONTACT_NO_MPI
#include "mpi.h"
#include "ContactZoltan.h"
#endif

#include <iostream>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <new>

struct ContactFaceFaceGraphNode {
  int       id;
  ContactFaceFaceGraphNode *next;
  ContactFaceFaceGraphNode(int v, ContactFaceFaceGraphNode* n) {id=v, next=n;};
};

struct ContactFaceFaceGraphHead {
  int       count;
  int       visit;
  ContactFaceFaceGraphNode *head;
};

typedef int Tuple[3];

extern "C" int compare_tuples(const void*, const void*);
void AddPolyVertex(ContactPolyVert*, ContactPolyVert**, int, Real);
void SortContactPolyVertices(ContactPolyVert**, int, int);
int  CompareVerts(ContactPolyVert*, ContactPolyVert*, int);
void MergeColinearEdges(ContactPolyVert**, int&, int&);
static ContactFaceFaceGraphNode* new_ContactFaceFaceGraphNode(
                                        ContactFixedSizeAllocator&,
                                        int v, ContactFaceFaceGraphNode* n);

void 
ContactSearch::Process_Face_Coverage()
{
  Real t3_node_positions[3][3];
  t3_node_positions[0][0] = -1.0;
  t3_node_positions[0][1] = -1.0;
  t3_node_positions[0][2] =  0.0;
  t3_node_positions[1][0] =  1.0;
  t3_node_positions[1][1] = -1.0;
  t3_node_positions[1][2] =  0.0;
  t3_node_positions[2][0] = -1.0;
  t3_node_positions[2][1] =  1.0;
  t3_node_positions[2][2] =  0.0;
  
  int number_of_faces = primary_topology->Number_of_Faces();
  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(primary_topology->FaceList()->EntityList());
  for (int i=0; i<number_of_faces; ++i) {
    ContactFace<Real>* face = Faces[i];
    ContactInteractionDLL<Real>* interactions = face->Get_FaceFace_Interactions();
    if(interactions != NULL) {
      interactions->IteratorStart();
      while (ContactInteractionEntity<Real>* interaction=interactions->IteratorForward()){
	int j, k;
	Real tol = 1.0e-6;
	ContactFaceFaceInteraction<Real>* cffi = 
	  static_cast<ContactFaceFaceInteraction<Real>*> (interaction);
	//ContactFace<Real>* face = cffi->SlaveFace();
	POSTCONDITION(face==cffi->SlaveFace());
#if CONTACT_DEBUG_PRINT_LEVEL>=8
	postream << "\n=======================================================\n";
	postream << "Slave face is " << face->Global_ID();
	switch( face->FaceType() ){
	case ContactSearch::QUADFACEL4 :
	  postream << "  Type is QUADFACEL4\n";
	  break;
	case ContactSearch::TRIFACEL3 :
	  postream << "  Type is TRIFACEL3\n";
	  break;
	default:
	  postream << "  Type is UNKNOWN\n";
	  break;
	}
#endif      

	// Determine the total number of polygons and vertices
	int nverts = face->Nodes_Per_Face();
	interactions->IteratorStart();
	while ((interaction=interactions->IteratorForward())){
	  cffi = static_cast<ContactFaceFaceInteraction<Real>*> (interaction);
	  nverts += cffi->NumEdges();
	}
	// add one to the number of face/face polys for the slave face
	int npolys = interactions->NumEntities()+1;
      
	// Convert the interactions to planar geometry
	int poly_edge_count=0;
	int vindex=0, pindex=0;
	ContactPolyVert** vlist = new ContactPolyVert* [nverts];
	ContactPoly**     plist = new ContactPoly* [npolys];
      
	// Add the nodes for the slave face
	switch( face->FaceType() ){
	case ContactSearch::QUADFACEL4 :
	  vlist[vindex] = ContactPolyVert::new_ContactPolyVert(
							       allocators[ALLOC_ContactPolyVert],
							       ContactPolyVert::SLAVE, vindex,
							       -1.0, -1.0); vindex++;
	  vlist[vindex] = ContactPolyVert::new_ContactPolyVert(
							       allocators[ALLOC_ContactPolyVert],
							       ContactPolyVert::SLAVE, vindex,  
							       1.0, -1.0); vindex++;
	  vlist[vindex] = ContactPolyVert::new_ContactPolyVert(
							       allocators[ALLOC_ContactPolyVert],
							       ContactPolyVert::SLAVE, vindex,  
							       1.0,  1.0); vindex++;
	  vlist[vindex] = ContactPolyVert::new_ContactPolyVert(
							       allocators[ALLOC_ContactPolyVert],
							       ContactPolyVert::SLAVE, vindex, 
							       -1.0,  1.0); vindex++;
	  break;
	case ContactSearch::TRIFACEL3 :
	  vlist[vindex] = ContactPolyVert::new_ContactPolyVert(
							       allocators[ALLOC_ContactPolyVert],
							       ContactPolyVert::SLAVE, vindex,
							       -1.0, -1.0); vindex++;
	  vlist[vindex] = ContactPolyVert::new_ContactPolyVert(
							       allocators[ALLOC_ContactPolyVert],
							       ContactPolyVert::SLAVE, vindex,  
							       1.0, -1.0); vindex++;
	  vlist[vindex] = ContactPolyVert::new_ContactPolyVert(
							       allocators[ALLOC_ContactPolyVert],
							       ContactPolyVert::SLAVE, vindex,  
							       -1.0,  1.0); vindex++;
	  break;
	default:
	  POSTCONDITION(false);
	  break;
	}
      
	// Add the nodes from the interactions
	k = 0;
	interactions->IteratorStart();
	while ((interaction=interactions->IteratorForward())){
	  cffi = static_cast<ContactFaceFaceInteraction<Real>*> (interaction);
	  ContactFaceFaceVertex<Real>* vertices = cffi->Get_Vertices();
	  plist[pindex] = ContactPoly::new_ContactPoly(allocators[ALLOC_ContactPoly],
						       pindex, cffi->NumEdges());
	  for (j=0; j<cffi->NumEdges(); j++, vindex++) {
	    int j0 = j-1;
	    if (j0<0) j0 = cffi->NumEdges()-1;
	    Real x = vertices[j].slave_x;
	    Real y = vertices[j].slave_y;
	    int  edge_node = 0;
	    //if (vertices[j].slave_edge_id>0 || 
	    //    vertices[j0].slave_edge_id>0) edge_node = 1;
	    if( face->FaceType() == ContactSearch::TRIFACEL3){
	      ContactTriFaceL3<Real>* t3face = static_cast<ContactTriFaceL3<Real>*>(face);
	      Real local_coords[3];
	      Real global_coords[3];
	      local_coords[0] = x;
	      local_coords[1] = y;
	      local_coords[2] = 1.0-x-y;
	      t3face->Compute_Global_Coords(t3_node_positions, local_coords, global_coords);
	      x = global_coords[0];
	      y = global_coords[1];
	      if (-tol<x+y && x+y<tol) edge_node = 1;
	    }
	    if (x<(-1.0+tol) || x>(1.0-tol) || 
		y<(-1.0+tol) || y>(1.0-tol)) edge_node = 1;
	    if (edge_node) {
	      vlist[vindex] = ContactPolyVert::new_ContactPolyVert(
								   allocators[ALLOC_ContactPolyVert],
								   ContactPolyVert::PERIMETER, 
								   vindex, x, y);
	    } else {
	      vlist[vindex] = ContactPolyVert::new_ContactPolyVert(
								   allocators[ALLOC_ContactPolyVert],
								   ContactPolyVert::INTERIOR, 
								   vindex, x, y);
	    }
	    plist[pindex]->Vertex(j, vlist[vindex]);
	  }
	  poly_edge_count += plist[pindex]->NumEdges();
	  pindex++;
	}
	int pslave = npolys-1;
#if CONTACT_DEBUG_PRINT_LEVEL>=8
	postream << "\nOriginal Vertex List --\n";
	for (j=0; j<nverts; ++j) {
	  ContactPolyVert* vertex = vlist[j];
	  postream << "  " << j << ":  (" 
		   << vertex->X() << ", " 
		   << vertex->Y() << ")  type = ";
	  switch (vertex->Type()) {
	  case ContactPolyVert::INTERIOR:
	    postream << "INTERIOR";
	    break;
	  case ContactPolyVert::PERIMETER:
	    postream << "PERIMETER";
	    break;
	  case ContactPolyVert::PEDGE:
	    postream << "PEDGE (" << vertex->TypeID() << ")";
	    break;
	  case ContactPolyVert::PVERTEX:
	    postream << "PVERTEX (" << vertex->TypeID() << ")";
	    break;
	  case ContactPolyVert::SLAVE:
	    postream << "SLAVE (" << vertex->TypeID() << ")";
	    break;
	  default:
	    postream << "UNKNOWN";
	    break;
	  }
	  postream << "\n";
	}
#endif
      
	// Equivalence all the nodes.  First, create a uniform spatial 
	// subdivision table.  Then for each node, determine it's 
	// position in the table and then test it against all the nodes
	// in it's cell and the neighboring cells.
	int nsize = 25;
	ContactPolyVert** table = new ContactPolyVert* [nsize*nsize];
	for (j=0; j<nsize*nsize; ++j) {
	  table[j] = NULL;
	}
	Real tol2 = tol*tol;
	for (j=0; j<nverts; ++j) {
	  ContactPolyVert *vert = vlist[j];
	  AddPolyVertex(vert, table, nsize, tol2);
	}
	delete [] table;
	for (j=0; j<nverts; ++j) {
	  vlist[j]->DetermineType(face->FaceType(), tol);
	}
#if CONTACT_DEBUG_PRINT_LEVEL>=8
	postream << "\nEquivalanced Vertex List --\n";
	for (j=0; j<nverts; ++j) {
	  ContactPolyVert* vertex = vlist[j];
	  postream << "  " << j << ":  (" 
		   << vertex->X() << ", " 
		   << vertex->Y() << ")  type = ";
	  switch (vertex->Type()) {
	  case ContactPolyVert::INTERIOR:
	    postream << "INTERIOR";
	    break;
	  case ContactPolyVert::PERIMETER:
	    postream << "PERIMETER";
	    break;
	  case ContactPolyVert::PEDGE:
	    postream << "PEDGE (" << vertex->TypeID() << ")";
	    break;
	  case ContactPolyVert::PVERTEX:
	    postream << "PVERTEX (" << vertex->TypeID() << ")";
	    break;
	  case ContactPolyVert::SLAVE:
	    postream << "SLAVE (" << vertex->TypeID() << ")";
	    break;
	  default:
	    postream << "UNKNOWN";
	    break;
	  }
	  postream << "  ==> " << vertex->Shared()->ID() << "\n";
	}
#endif
      
	// adjust the polygons to point to the equivalenced vertices
	for (j=0; j<pslave; ++j) {
	  ContactPoly* poly = plist[j];
	  for (k=0; k<poly->NumVerts(); ++k) {
	    ContactPolyVert* vertex = poly->Vertex(k)->Shared();
	    poly->Vertex(k, vertex);
	  }
	  poly->EquivalanceVertices();
	}
      
	SortContactPolyVertices(vlist, nverts, face->FaceType());
	int num_active_verts    = 0;
	int num_perimeter_verts = 0;
	for (j=0; j<nverts; ++j) {
	  vlist[j]->ID(j);
	  if (vlist[j]->Type()!=ContactPolyVert::UNKNOWN) {
	    num_active_verts++;
	    if (vlist[j]->Type()!=ContactPolyVert::INTERIOR) {
	      num_perimeter_verts++;
	    }
	  }
	}
      
	// add the polygon for the slave face using all the perimeter vertices
	plist[pslave] = ContactPoly::new_ContactPoly(
						     allocators[ALLOC_ContactPoly],
						     pslave, num_perimeter_verts);
	for (j=0; j<num_perimeter_verts; ++j) {
	  plist[pslave]->Vertex(j, vlist[j]);
	}
	poly_edge_count += plist[pslave]->NumEdges();
      
	// generate edge tuples, process all the polygon edges in
	// a clockwise direction except for the slave face polygon
	// which is processed in a counter-clockwise direction
	int ntuples = 0;
	Tuple* edge_tuples = new Tuple[poly_edge_count];
	for (j=0; j<pslave; ++j) {
	  ContactPoly* poly = plist[j];
	  int v2 = poly->NumEdges()-1;
	  int v1 = (v2+1)%poly->NumEdges();
	  for (k=0; k<poly->NumEdges(); 
	       k++, v2--, v1=(v2+1)%poly->NumEdges()) {
	    edge_tuples[ntuples][0] = std::min(poly->Vertex(v1)->ID(), 
					       poly->Vertex(v2)->ID());
	    edge_tuples[ntuples][1] = std::max(poly->Vertex(v1)->ID(), 
					       poly->Vertex(v2)->ID());
	    edge_tuples[ntuples][2] = ntuples;
	    ntuples++;
	  }
	}
	ContactPoly* poly = plist[pslave];
	for (k=0; k<poly->NumEdges(); ++k) {
	  int v1 = k;
	  int v2 = (k+1)%poly->NumVerts();
	  edge_tuples[ntuples][0] = std::min(poly->Vertex(v1)->ID(), 
					     poly->Vertex(v2)->ID());
	  edge_tuples[ntuples][1] = std::max(poly->Vertex(v1)->ID(), 
					     poly->Vertex(v2)->ID());
	  edge_tuples[ntuples][2] = ntuples;
	  ntuples++;
	}

	// sort the edge tuples
	std::qsort((void*)edge_tuples, ntuples, sizeof(Tuple), compare_tuples);

	// set up a transfer function from edge tuple number to final edge number
	int  current  = 0;
	int  nedges   = 0;
	int* transfer = new int [ntuples];
	for (j=0; j<ntuples; ++j) {
	  if ((edge_tuples[j][0] != edge_tuples[current][0]) ||
	      (edge_tuples[j][1] != edge_tuples[current][1])) {
	    current = j;
	    nedges++;
	  }
	  transfer[edge_tuples[j][2]] = nedges;
	}
	nedges++;
	delete [] edge_tuples;
	ContactPolyEdge** elist = new ContactPolyEdge* [nedges];
	for (j=0; j<nedges; ++j) {
	  elist[j] = ContactPolyEdge::new_ContactPolyEdge(allocators[ALLOC_ContactPolyEdge],j);
	}
      
	// set up the edges, process all the polygon edges in
	// a clockwise direction except for the slave face polygon
	// which is processed in a counter-clockwise direction
	ntuples = 0;

	for (j=0; j<pslave; ++j) {
	  poly   = plist[j];
	  int v2 = poly->NumEdges()-1;
	  int v1 = (v2+1)%poly->NumEdges();
	  for (k=0; k<poly->NumEdges(); 
	       k++, v2--, v1=(v2+1)%poly->NumEdges()) {
	    ContactPolyEdge* edge = elist[transfer[ntuples++]];
	    poly->Edge(k,edge);
	    if (edge->Poly1() == NULL) {
	      edge->Vertex1(poly->Vertex(v1));
	      edge->Vertex2(poly->Vertex(v2));
	      edge->Poly1(poly);
	    } else if (edge->Poly2() == NULL) {
	      edge->Poly2(poly);
	    }
	  }
	  interactions->IteratorForward();;
	}
	poly = plist[pslave];
	for (k=0; k<poly->NumEdges(); ++k) {
	  int v1 = k;
	  int v2 = (k+1)%poly->NumEdges();
	  ContactPolyEdge* edge = elist[transfer[ntuples++]];
	  poly->Edge(k,edge);
	  if (edge->Poly1() == NULL) {
	    edge->Vertex1(poly->Vertex(v1));
	    edge->Vertex2(poly->Vertex(v2));
	    edge->Poly1(poly);
	  } else if (edge->Poly2() == NULL) {
	    edge->Poly2(poly);
	  }
	}
	delete [] transfer;
          
#if CONTACT_DEBUG_PRINT_LEVEL>=8
	postream << "\nActive Vertex List --\n";
	postream << "Number of active vertices    = " << num_active_verts << "\n";
	postream << "Number of perimeter vertices = " << num_perimeter_verts << "\n";
	for (j=0; j<num_active_verts; ++j) {
	  ContactPolyVert* vertex = vlist[j];
	  postream << "  " << j << ":  (" 
		   << vertex->X() << ", " 
		   << vertex->Y() << ")  type = ";
	  switch (vertex->Type()) {
	  case ContactPolyVert::INTERIOR:
	    postream << "INTERIOR\n";
	    break;
	  case ContactPolyVert::PERIMETER:
	    postream << "PERIMETER\n";
	    break;
	  case ContactPolyVert::PEDGE:
	    postream << "PEDGE (" << vertex->TypeID() << ")\n";
	    break;
	  case ContactPolyVert::PVERTEX:
	    postream << "PVERTEX (" << vertex->TypeID() << ")\n";
	    break;
	  case ContactPolyVert::SLAVE:
	    postream << "SLAVE (" << vertex->TypeID() << ")\n";
	    break;
	  default:
	    postream << "UNKNOWN\n";
	    break;
	  }
	}
      
	postream << "\nActive Edge List --\n";
	for (j=0; j<nedges; ++j) {
	  ContactPolyEdge* edge = elist[j];
	  postream << "  " << j << ":  Vertices = ("  
		   << edge->Vertex1()->ID() << ", " 
		   << edge->Vertex2()->ID() << ")   Polys = (";
	  if (edge->Poly1() != NULL) {
	    postream<<edge->Poly1()->ID()<<", ";
	  } else {
	    postream<<"NULL"<<", ";
	  }
	  if (edge->Poly2() != NULL) {
	    postream<<edge->Poly2()->ID()<<")\n";
	  } else {
	    postream<<"NULL)\n";
	  }
	}
      
	postream << "\nActive Poly List --\n";
	for (j=0; j<npolys; ++j) {
	  ContactPoly* poly = plist[j];
	  postream << "  " << j << ":  Verts = ";
	  for (k=0; k<poly->NumVerts(); ++k) {
	    postream << "  " << poly->Vertex(k)->ID();
	  }
	  postream << "\n      Edges = ";
	  for (k=0; k<poly->NumVerts(); ++k) {
	    postream << "  " << poly->Edge(k)->ID();
	  }
	  postream << "\n";
	}
	postream.flush();
#endif

	// Construct a directed graph if all the edges of 
	// all the polygons that are not multiply connected. 
	struct ContactFaceFaceGraphHead *graph = 
	  new struct ContactFaceFaceGraphHead [num_active_verts];
	for (j=0; j<num_active_verts; ++j) {
	  graph[j].count = 0;
	  graph[j].visit = 0;
	  graph[j].head  = NULL;
	}
	int max_count = 0;
	for (j=0; j<nedges; ++j) {
	  ContactPolyEdge* edge = elist[j];
	  if (edge->Poly1()!=NULL && edge->Poly2()==NULL) {
	    ContactPolyVert* vert0 = edge->Vertex1();
	    ContactPolyVert* vert1 = edge->Vertex2();
	    int id0 = vert0->ID();
	    int id1 = vert1->ID();
	    graph[id0].count++;
	    max_count = std::max(max_count,graph[id0].count);
	    if (graph[id0].head) {
	      graph[id0].head = new_ContactFaceFaceGraphNode(allocators[ALLOC_ContactFaceFaceGraphNode],id1,graph[id0].head);
	    } else {
	      graph[id0].head = new_ContactFaceFaceGraphNode(allocators[ALLOC_ContactFaceFaceGraphNode],id1,NULL);
	    }
	  }
	}
#if CONTACT_DEBUG_PRINT_LEVEL>=8
	postream << "\nVertex Graph List --\n";
	for (j=0; j<num_active_verts; ++j) {
	  postream << j << ": N=" << graph[j].count << "  =>";
	  ContactFaceFaceGraphNode* link = graph[j].head;
	  while (link) {
	    postream << "  " << link->id;
	    link = link->next;
	  }
	  postream << "\n";
	}
	postream.flush();
#endif

	int count  = 0;
	int index0 = -1;
	int index1 = -1;
	// Traverse the graph to build the partially exposed face polygon(s)
	if (max_count>1) {
	  for (j=0; j<num_active_verts; ++j) {
	    if (graph[j].count==max_count) break;
	  }
	  do {
	    for (k=0; k<graph[j].count; ++k) {
	      int n;
	      ContactFaceCoverageInteraction* cfci =
		ContactFaceCoverageInteraction::new_ContactFaceCoverageInteraction(
										   allocators[ALLOC_ContactFaceCoverageInteraction], face );
	      face->Store_FaceCoverage_Interaction( cfci );
	      count++;
	      index0 = j;
	      index1 = j;
	      vlist[index1]->Next(NULL);
	      ContactFaceFaceGraphNode* link = graph[j].head;
	      for (n=0; n<graph[j].visit; ++n) {
		link = link->next;
	      }
	      graph[j].visit++;
	      int id = link->id;
	      while (id!=j) {
		vlist[id]->Next(NULL);
		vlist[index1]->Next(vlist[id]);
		index1 = id;
		count++;
		link = graph[id].head;
		for (n=0; n<graph[id].visit; ++n) {
		  link = link->next;
		}
		graph[id].visit++;
		id = link->id;
	      }
	      MergeColinearEdges(vlist, index0, count);
	      ContactPolyVert* llnode = vlist[index0];
	      while (llnode) {
		Real x = llnode->X();
		Real y = llnode->Y();
		if (face->FaceType()==TRIFACEL3) {
		  Real local_coords[3];
		  Real global_coords[3];
		  global_coords[0] = x;
		  global_coords[1] = y;
		  global_coords[2] = 0.0;
		  ContactTriFaceL3<Real>* t3face = static_cast<ContactTriFaceL3<Real>*>(face);
		  t3face->Compute_Local_Coords(t3_node_positions, global_coords, local_coords);
		  x = local_coords[0];
		  y = local_coords[1];
		}
		cfci->AddVertex(x, y);
		llnode = llnode->Next();
	      }
	      count  = 0;
	      index0 = -1;
	      index1 = -1;
	    }
	    for (j=0; j<num_active_verts; ++j) {
	      if (graph[j].count-graph[j].visit>1) break;
	    }
	    if (j==num_active_verts) {
	      for (j=0; j<num_active_verts; ++j) {
		if (graph[j].count-graph[j].visit==1) break;
	      }
	    }
	  } while (j!=num_active_verts);
	}
	if (max_count>=1) {
	  for (j=0; j<num_active_verts; ++j) {
	    while (graph[j].count-graph[j].visit>0) {
	      ContactFaceCoverageInteraction* cfci =
		ContactFaceCoverageInteraction::new_ContactFaceCoverageInteraction(
										   allocators[ALLOC_ContactFaceCoverageInteraction], face );
	      face->Store_FaceCoverage_Interaction( cfci );
	      count++;
	      index0 = j;
	      index1 = j;
	      vlist[index1]->Next(NULL);
	      ContactFaceFaceGraphNode* link = graph[j].head;
	      for (k=0; k<graph[j].visit; ++k) {
		link = link->next;
	      }
	      graph[j].visit++;
	      int id = link->id;
	      while (id!=j) {
		vlist[id]->Next(NULL);
		vlist[index1]->Next(vlist[id]);
		index1 = id;
		count++;
		link = graph[id].head;
		for (k=0; k<graph[id].visit; ++k) {
		  link = link->next;
		}
		graph[id].visit++;
		id = link->id;
	      }
#if CONTACT_DEBUG_PRINT_LEVEL>=8
	      postream << "\nInitial vertex list --\n";
	      postream << "  ";
	      ContactPolyVert* llnode1 = vlist[index0];
	      while (llnode1) {
		postream << "  " << llnode1->ID();
		llnode1 = llnode1->Next();
	      }
	      postream << "\n";
	      postream.flush();
#endif
	      MergeColinearEdges(vlist, index0, count);
#if CONTACT_DEBUG_PRINT_LEVEL>=8
	      postream << "\nMerged vertex list --\n";
	      postream << "  ";
	      llnode1 = vlist[index0];
	      while (llnode1) {
		postream << "  " << llnode1->ID();
		llnode1 = llnode1->Next();
	      }
	      postream << "\n";
	      postream.flush();
#endif
	      ContactPolyVert* llnode = vlist[index0];
	      while (llnode) {
		Real x = llnode->X();
		Real y = llnode->Y();
		if (face->FaceType()==TRIFACEL3) {
		  Real local_coords[3];
		  Real global_coords[3];
		  global_coords[0] = x;
		  global_coords[1] = y;
		  global_coords[2] = 0.0;
		  ContactTriFaceL3<Real>* t3face = static_cast<ContactTriFaceL3<Real>*>(face);
		  t3face->Compute_Local_Coords(t3_node_positions, global_coords, local_coords);
		  x = local_coords[0];
		  y = local_coords[1];
		}
		cfci->AddVertex(x, y);
		llnode = llnode->Next();
	      }
#if CONTACT_DEBUG_PRINT_LEVEL>=8
	      postream << "\ncfci num vertices = " << cfci->NumVertices() << "\n";
	      postream.flush();
#endif
	      count  = 0;
	      index0 = -1;
	      index1 = -1;
	    }
	  }
	} else {
	  // Fully covered face, create a coverage interaction with no polys
	  ContactFaceCoverageInteraction* cfci =
	    ContactFaceCoverageInteraction::new_ContactFaceCoverageInteraction(
									       allocators[ALLOC_ContactFaceCoverageInteraction], face );
	  face->Store_FaceCoverage_Interaction( cfci );
	}
      
	// Cleanup
	for (j=0; j<num_active_verts; ++j) {
	  ContactFaceFaceGraphNode* link = graph[j].head;
	  while (link) {
	    ContactFaceFaceGraphNode* next = link->next;
	    //delete link;
	    allocators[ALLOC_ContactFaceFaceGraphNode].Delete_Frag( link );
	    link = next;
	  }
	}
	delete [] graph;
      
	for (j=0; j<nedges; ++j) {
	  elist[j]->~ContactPolyEdge();
	  allocators[ALLOC_ContactPolyEdge].Delete_Frag(elist[j]);
	}
	delete [] elist;
      
	for (j=0; j<nverts; ++j) {
	  vlist[j]->~ContactPolyVert();
	  allocators[ALLOC_ContactPolyVert].Delete_Frag(vlist[j]);
	}
	delete [] vlist;
      
	for (j=0; j<npolys; ++j) {
	  plist[j]->~ContactPoly();
	  allocators[ALLOC_ContactPoly].Delete_Frag(plist[j]);
	}
	delete [] plist;
      }
    }
  }
}

void AddPolyVertex(ContactPolyVert* vert, 
                   ContactPolyVert** table, 
                   int nsize, Real tol2)
{
  int ix = (int)std::floor(0.5*(vert->X()+1.1)*(nsize-1))%nsize;
  int iy = (int)std::floor(0.5*(vert->Y()+1.1)*(nsize-1))%nsize;
  int x0 = std::max(0,ix-1);
  int x1 = std::min(nsize-1,ix+1);
  int y0 = std::max(0,iy-1);
  int y1 = std::min(nsize-1,iy+1);
  Real min_d = BIGNUM;
  ContactPolyVert* min_ptr = NULL;
  for (int x=x0; x<=x1; x++) {
    for (int y=y0; y<=y1; y++) {
      int index = y*nsize+x;
      ContactPolyVert *ptr;
      for (ptr=table[index]; ptr!=NULL; ptr=ptr->Next()) {
        Real dx = ptr->X() - vert->X();
        Real dy = ptr->Y() - vert->Y();
        Real d  = dx*dx + dy*dy;
        if (d <= min_d) {
          min_d   = d;
          min_ptr = ptr;
        }
      }
    }
  }
  if ((min_ptr && min_d <= tol2)) {
    // match found, set this vertex to the matching vertex
    if (min_ptr->Type()==ContactPolyVert::SLAVE) {
      min_ptr->Type(vert->Type());
      min_ptr->TypeID(vert->TypeID());
    }
    vert->Shared(min_ptr);
    vert->Type(ContactPolyVert::UNKNOWN);
  } else {
    // no match, add this vertex to the table
    int index = iy*nsize+ix;
    vert->Next(table[index]);
    table[index] = vert;
    vert->Shared(vert);
  }
}

extern "C" int compare_tuples(const void *t1, const void *t2)
{
  int    status;
  Tuple* tuple1;
  Tuple* tuple2;
  
  tuple1 = (Tuple*)t1;
  tuple2 = (Tuple*)t2;
  status = (((*tuple1)[0] < (*tuple2)[0]) ? -1 :
	    (((*tuple1)[0] > (*tuple2)[0]) ?  1 :
	     (((*tuple1)[1] < (*tuple2)[1]) ? -1 :
	      (((*tuple1)[1] > (*tuple2)[1]) ?  1 : 0))));
  return status;
}

void SortContactPolyVertices(ContactPolyVert** vlist, int nverts, int face_type)
{
  if (nverts>1) {
    int i, j, k, n;
    ContactPolyVert* vert;
    k = (nverts>>1)+1;
    n = nverts;
    for (;;) {
      if (k>1) {
        vert = vlist[--k-1];
      } else {
        vert       = vlist[n-1];
        vlist[n-1] = vlist[0];
        if (--n == 1) {
          vlist[0] = vert;
          break;
        }
      }
      i = k;
      j = k<<1;
      while (j<=n) {
        if ((j<n) && (CompareVerts(vlist[j-1],vlist[j],face_type)<0)) ++j;
        if (CompareVerts(vert,vlist[j-1],face_type)<0) {
          vlist[i-1] = vlist[j-1];
          j += (i=j);
        } else {
          j = n+1;
        }
      }
      vlist[i-1] = vert;
    }
  }
}

int CompareVerts(ContactPolyVert* v1, ContactPolyVert* v2, int face_type)
{
  // 1st list perimeter nodes in ccw direction
  // 2nd list interior nodes
  // 3rd list all others in any order
  int type1 = v1->Type();
  int type2 = v2->Type();
  if (type1!=ContactPolyVert::UNKNOWN && 
      type2==ContactPolyVert::UNKNOWN) return -1;
  if (type1==ContactPolyVert::UNKNOWN && 
      type2==ContactPolyVert::UNKNOWN) return  0;
  if (type1==ContactPolyVert::UNKNOWN && 
      type2!=ContactPolyVert::UNKNOWN) return  1;

  if (type1!=ContactPolyVert::INTERIOR && 
      type2==ContactPolyVert::INTERIOR) return -1;
  if (type1==ContactPolyVert::INTERIOR && 
      type2!=ContactPolyVert::INTERIOR) return  1;
  if (type1==ContactPolyVert::INTERIOR && 
      type2==ContactPolyVert::INTERIOR) {
    Real x1 = v1->X();
    Real y1 = v1->Y();
    Real x2 = v2->X();
    Real y2 = v2->Y();
    if (x1<x2)  return -1;
    if (x1>x2)  return  1;
    if (x1==x2) {
      if (y1<y2)  return -1;
      if (y1>y2)  return  1;
      return 0;
    }
  }
  
  if ((v1->Type()==ContactPolyVert::PVERTEX ||
       v1->Type()==ContactPolyVert::SLAVE) &&
      (v2->Type()==ContactPolyVert::PVERTEX ||
       v2->Type()==ContactPolyVert::SLAVE)) {
    int node1 = v1->TypeID();
    int node2 = v2->TypeID();
    if (node1==node2) return  0;
    if (node1<node2)  return -1;
    if (node1>node2)  return  1;
  }
  if (v1->Type()==ContactPolyVert::PEDGE &&
      v2->Type()==ContactPolyVert::PEDGE) {
    int edge1 = v1->TypeID();
    int edge2 = v2->TypeID();
    if (edge1<edge2)  return -1;
    if (edge1>edge2)  return  1;
    if (edge1==edge2) {
      Real x1 = v1->X();
      Real y1 = v1->Y();
      Real x2 = v2->X();
      Real y2 = v2->Y();
      switch (face_type) {
      case ContactSearch::TRIFACEL3:
        switch (edge1) {
        case 0:
          if (x1<x2)  return -1;
          if (x1>x2)  return  1;
          if (x1==x2) return  0;
          break;
        case 1:
          if (x1>x2)  return -1;
          if (x1<x2)  return  1;
          if (x1==x2) return  0;
          break;
        case 2:
          if (y1>y2)  return -1;
          if (y1<y2)  return  1;
          if (y1==y2) return  0;
          break;
        }
        break;
      case ContactSearch::QUADFACEL4:
        switch (edge1) {
        case 0:
          if (x1<x2)  return -1;
          if (x1>x2)  return  1;
          if (x1==x2) return  0;
          break;
        case 1:
          if (y1<y2)  return -1;
          if (y1>y2)  return  1;
          if (y1==y2) return  0;
          break;
        case 2:
          if (x1>x2)  return -1;
          if (x1<x2)  return  1;
          if (x1==x2) return  0;
          break;
        case 3:
          if (y1>y2)  return -1;
          if (y1<y2)  return  1;
          if (y1==y2) return  0;
          break;
        }
        break;
      }
    }
  }
  if (v1->Type()==ContactPolyVert::PEDGE) {
    int edge1 = v1->TypeID();
    int node2 = v2->TypeID();
    switch (face_type) {
    case ContactSearch::TRIFACEL3:
      if (edge1==node2) return  1;
      if (edge1<node2)  return -1;
      if (edge1>node2)  return  1;
      break;
    case ContactSearch::QUADFACEL4:
      if (edge1==node2) return  1;
      if (edge1<node2)  return -1;
      if (edge1>node2)  return  1;
      break;
    }
  }
  if (v2->Type()==ContactPolyVert::PEDGE) {
    int node1 = v1->TypeID();
    int edge2 = v2->TypeID();
    switch (face_type) {
    case ContactSearch::TRIFACEL3:
      if (node1==edge2) return -1;
      if (node1<edge2)  return -1;
      if (node1>edge2)  return  1;
      break;
    case ContactSearch::QUADFACEL4:
      if (node1==edge2) return -1;
      if (node1<edge2)  return -1;
      if (node1>edge2)  return  1;
      break;
    }
  }
  
  // Should never get here but must have a return value
  POSTCONDITION( 0 );
  return -1024;
}

void MergeColinearEdges(ContactPolyVert** vlist, int& index0, int& count)
{
  if (count<4) return;

  int ndelete = 0;
  ContactPolyVert* llnode=vlist[index0];
  while (llnode) {
    ContactPolyVert* vert1 = llnode;
    ContactPolyVert* vert2 = vert1->Next();
    if (!vert2) {
      vert2 = vlist[index0];
    }
    ContactPolyVert* vert3 = vert2->Next();
    if (!vert3) {
      vert3 = vlist[index0];
    }
    Real dx1  = vert1->X() - vert2->X();
    Real dy1  = vert1->Y() - vert2->Y();
    Real len1 = std::sqrt(dx1*dx1+dy1*dy1);
    dx1 /= len1;
    dy1 /= len1;
    Real dx2  = vert3->X() - vert2->X();
    Real dy2  = vert3->Y() - vert2->Y();
    Real len2 = std::sqrt(dx2*dx2+dy2*dy2);
    dx2 /= len2;
    dy2 /= len2;
    Real dot_product = dx1*dx2+dy1*dy2;
    if (dot_product<(-1.0+1.e-6)) {
      // delete vertex2
      if (vert2==vlist[index0]) {
        index0 = vert3->ID();
      } else if (vert3==vlist[index0]) {
        vert1->Next(NULL);
      } else {
        vert1->Next(vert3);
      }
      vert2 = vert3;
      //ContactPolyVert* vert3 = vert2->Next();
      vert3 = vert2->Next();
      if (!vert3) {
        vert3 = vlist[index0];
      }
      ndelete++;
    } else {
      llnode = llnode->Next();
    }
  }
  count -= ndelete;
}

ContactFaceFaceGraphNode* new_ContactFaceFaceGraphNode(
                                        ContactFixedSizeAllocator& alloc,
                                        int v, ContactFaceFaceGraphNode* n)
{
  return new (alloc.New_Frag())ContactFaceFaceGraphNode(v, n);
}

void ContactFaceFaceGraphNode_SizeAllocator( ContactFixedSizeAllocator& alloc)
{
  alloc.Resize( sizeof(ContactFaceFaceGraphNode),
                100,  // block size
                0);   // initial block size
  alloc.Set_Name( "ContactFaceFaceGraphNode allocator" );
}
