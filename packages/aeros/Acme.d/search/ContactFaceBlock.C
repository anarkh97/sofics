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


#include "ContactFaceBlock.h"
#include "ContactQuadFaceL4.h"
#include "ContactQuadFaceQ8.h"
#include "ContactQuadFaceQ9.h"
#include "ContactTriFaceL3.h"
#include "ContactTriFaceQ6.h"
#include "ContactShellQuadFaceL4.h"
#include "ContactShellTriFaceL3.h"
#include "ContactLineFaceL2.h"
#include "ContactLineFaceQ3.h"
#include "ContactSearch.h"
#include "ContactTopology.h"
#include "ContactFixedSizeAllocator.h"

#include <iostream>

ContactFaceBlock::ContactFaceBlock( ContactSearch::ContactFace_Type Type,
				    int block_index,
				    int Number_of_Faces,
				    int& Next_ID,
                                    int* host_ids,
                                    ContactTopology* Topology )
  : topology(Topology)
{
  PRECONDITION( Number_of_Faces >= 0 );
  PRECONDITION( block_index >= 0 );
  
  ContactSearch* search = topology->Search();

  int j;

  master          = false;
  slave           = false;
  face_list       = new ContactBlockEntityList(Topology->Search()->Get_Allocators(), Number_of_Faces);
  type            = Type;
  id              = block_index;
  entity_key      = topology->Number_of_Element_Blocks()+block_index;
  number_of_faces = Number_of_Faces;
  num_faces_added = Number_of_Faces;

  int hostID = 0;

  edge_type = ContactSearch::Face_Edge_Type(type);

  switch( type ){
  case ContactSearch::LINEFACEL2:
    for( j=0 ; j<number_of_faces ; ++j ) {
      ContactFace<Real>* face = ContactLineFaceL2::new_ContactLineFaceL2(
           search->Get_Allocators(),
           block_index,hostID++,entity_key);
      face->Global_ID(host_ids[2*j],host_ids[2*j+1]);
      face_list->Insert(face);
      ++Next_ID;
    }
    break;
  case ContactSearch::LINEFACEQ3:
    for( j=0 ; j<number_of_faces ; ++j ) {
      ContactFace<Real>* face = ContactLineFaceQ3::new_ContactLineFaceQ3(
           search->Get_Allocators(),
           block_index,hostID++,entity_key);
      face->Global_ID(host_ids[2*j],host_ids[2*j+1]);
      face_list->Insert(face);
      ++Next_ID;
    }
    break;
  case ContactSearch::QUADFACEL4:
    for( j=0 ; j<number_of_faces ; ++j ) {
      ContactFace<Real>* face = ContactQuadFaceL4<Real>::new_ContactQuadFaceL4(
           search->Get_Allocators(),
           block_index,hostID++,entity_key );
      face->Global_ID(host_ids[2*j],host_ids[2*j+1]);
      face_list->Insert(face);
      ++Next_ID;
    }
    break;
  case ContactSearch::QUADFACEQ8:
    for( j=0 ; j<number_of_faces ; ++j ) {
      ContactFace<Real>* face = ContactQuadFaceQ8::new_ContactQuadFaceQ8(
           search->Get_Allocators(),
           block_index,hostID++,entity_key );
      face->Global_ID(host_ids[2*j],host_ids[2*j+1]);
      face_list->Insert(face);
      ++Next_ID;
    }
    break;
  case ContactSearch::QUADFACEQ9:
    for( j=0 ; j<number_of_faces ; ++j ) {
      ContactFace<Real>* face = ContactQuadFaceQ9::new_ContactQuadFaceQ9(
           search->Get_Allocators(),
           block_index,hostID++,entity_key );
      face->Global_ID(host_ids[2*j],host_ids[2*j+1]);
      face_list->Insert(face);
      ++Next_ID;
    }
    break;
  case ContactSearch::TRIFACEL3:
    for( j=0 ; j<number_of_faces ; ++j ) {
      ContactFace<Real>* face = ContactTriFaceL3<Real>::new_ContactTriFaceL3(
           search->Get_Allocators(),
           block_index,hostID++,entity_key );
      face->Global_ID(host_ids[2*j],host_ids[2*j+1]);
      face_list->Insert(face);
      ++Next_ID;
    }
    break;
  case ContactSearch::TRIFACEQ6:
    for( j=0 ; j<number_of_faces ; ++j ) {
      ContactFace<Real>* face = ContactTriFaceQ6::new_ContactTriFaceQ6(
           search->Get_Allocators(),
           block_index,hostID++,entity_key );
      face->Global_ID(host_ids[2*j],host_ids[2*j+1]);
      face_list->Insert(face);
      ++Next_ID;
    }
    break;
  case ContactSearch::SHELLQUADFACEL4:
    for( j=0 ; j<number_of_faces ; ++j ) {
      ContactFace<Real>* face = ContactShellQuadFaceL4<Real>::new_ContactShellQuadFaceL4(
           search->Get_Allocators(),
           block_index,hostID++,entity_key );
      face->Global_ID(host_ids[2*j],host_ids[2*j+1]);
      face_list->Insert(face);
      ++Next_ID;
    }
    break;
  case ContactSearch::SHELLTRIFACEL3:
    for( j=0 ; j<number_of_faces ; ++j ) {
      ContactFace<Real>* face = ContactShellTriFaceL3<Real>::new_ContactShellTriFaceL3(
           search->Get_Allocators(),
           block_index,hostID++,entity_key );
      face->Global_ID(host_ids[2*j],host_ids[2*j+1]);
      face_list->Insert(face);
      ++Next_ID;
    }
    break;
  default:
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Unknown face type in ContactFaceBlock" << std::endl;
#endif
    POSTCONDITION( 0 );
  }  
}

ContactFaceBlock::ContactFaceBlock( ContactSearch::ContactFace_Type Type,
				    int block_index, ContactTopology* Topology )
  : topology(Topology)
{
  master          = false;
  slave           = false;
  face_list       = new ContactBlockEntityList(Topology->Search()->Get_Allocators());
  type            = Type;
  id              = block_index;
  entity_key      = block_index;
  number_of_faces = 0;
  num_faces_added = 0;
  edge_type = ContactSearch::Face_Edge_Type(type);
}

ContactFaceBlock::~ContactFaceBlock()
{
  Delete_Faces();
  delete face_list;
}

void
ContactFaceBlock::Add_Faces( int Number_of_Faces, int* host_ids )
{
  PRECONDITION( Number_of_Faces >= 0 );
  
  ContactSearch* search = topology->Search();

  int j;

  number_of_faces += Number_of_Faces;

  edge_type = ContactSearch::Face_Edge_Type(type);

  switch( type ){
  case ContactSearch::LINEFACEL2:
    for( j=0 ; j<Number_of_Faces ; ++j ) {
      ContactFace<Real>* face = ContactLineFaceL2::new_ContactLineFaceL2(
           search->Get_Allocators(),id,-1,entity_key);
      face->Global_ID(host_ids[2*j],host_ids[2*j+1]);
      face_list->Insert(face);
    }
    break;
  case ContactSearch::LINEFACEQ3:
    for( j=0 ; j<Number_of_Faces ; ++j ) {
      ContactFace<Real>* face = ContactLineFaceQ3::new_ContactLineFaceQ3(
           search->Get_Allocators(),id,-1,entity_key);
      face->Global_ID(host_ids[2*j],host_ids[2*j+1]);
      face_list->Insert(face);
    }
    break;
  case ContactSearch::QUADFACEL4:
    for( j=0 ; j<Number_of_Faces ; ++j ) {
      ContactFace<Real>* face = ContactQuadFaceL4<Real>::new_ContactQuadFaceL4(
           search->Get_Allocators(),id,-1,entity_key );
      POSTCONDITION(face);
      face->Global_ID(host_ids[2*j],host_ids[2*j+1]);
      face_list->Insert(face);
    }
    break;
  case ContactSearch::QUADFACEQ8:
    for( j=0 ; j<Number_of_Faces ; ++j ) {
      ContactFace<Real>* face = ContactQuadFaceQ8::new_ContactQuadFaceQ8(
           search->Get_Allocators(),id,-1,entity_key );
      face->Global_ID(host_ids[2*j],host_ids[2*j+1]);
      face_list->Insert(face);
    }
    break;
  case ContactSearch::QUADFACEQ9:
    for( j=0 ; j<Number_of_Faces ; ++j ) {
      ContactFace<Real>* face = ContactQuadFaceQ9::new_ContactQuadFaceQ9(
           search->Get_Allocators(),id,-1,entity_key );
      face->Global_ID(host_ids[2*j],host_ids[2*j+1]);
      face_list->Insert(face);
    }
    break;
  case ContactSearch::TRIFACEL3:
    for( j=0 ; j<Number_of_Faces ; ++j ) {
      ContactFace<Real>* face = ContactTriFaceL3<Real>::new_ContactTriFaceL3(
           search->Get_Allocators(),id,-1,entity_key );
      face->Global_ID(host_ids[2*j],host_ids[2*j+1]);
      face_list->Insert(face);
    }
    break;
  case ContactSearch::TRIFACEQ6:
    for( j=0 ; j<Number_of_Faces ; ++j ) {
      ContactFace<Real>* face = ContactTriFaceQ6::new_ContactTriFaceQ6(
           search->Get_Allocators(),id,-1,entity_key );
      face->Global_ID(host_ids[2*j],host_ids[2*j+1]);
      face_list->Insert(face);
    }
    break;
  case ContactSearch::SHELLQUADFACEL4:
    for( j=0 ; j<Number_of_Faces ; ++j ) {
      ContactFace<Real>* face = ContactShellQuadFaceL4<Real>::new_ContactShellQuadFaceL4(
           search->Get_Allocators(),id,-1,entity_key );
      face->Global_ID(host_ids[2*j],host_ids[2*j+1]);
      face_list->Insert(face);
    }
    break;
  case ContactSearch::SHELLTRIFACEL3:
    for( j=0 ; j<Number_of_Faces ; ++j ) {
      ContactFace<Real>* face = ContactShellTriFaceL3<Real>::new_ContactShellTriFaceL3(
           search->Get_Allocators(),id,-1,entity_key );
      face->Global_ID(host_ids[2*j],host_ids[2*j+1]);
      face_list->Insert(face);
    }
    break;
  default:
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Unknown face type in ContactFaceBlock" << std::endl;
#endif
    POSTCONDITION( 0 );
  }
}

void ContactFaceBlock::Delete_Face_List()
{
  Delete_Faces();
  face_list->CleanUp();
  number_of_faces = 0;
  num_faces_added = 0;
}

void ContactFaceBlock::Delete_Faces()
{
  if( face_list ){
    ContactFixedSizeAllocator* alloc=ContactSearch::Get_ContactFaceAlloc(type, 
         topology->Search()->Get_Allocators());
    POSTCONDITION(alloc!=NULL);
    face_list->IteratorStart();
    while (ContactTopologyEntity<Real>* entity = face_list->IteratorForward()) {
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
      face->~ContactFace<Real>();
      alloc->Delete_Frag(face);
    }
  }
}

void ContactFaceBlock::Insert_Face( ContactFace<Real>* face )
{
  PRECONDITION( face );
  face_list->Insert(face);
  number_of_faces = face_list->NumEntities();
}

#ifndef CONTACT_NO_MPI
void ContactFaceBlock::Insert_Face( char* buffer )
{
  PRECONDITION( buffer );
  face_list->Insert(buffer);
  number_of_faces = face_list->NumEntities();
}
void ContactFaceBlock::Insert_Face_ForSecondary( char* buffer )
{
  PRECONDITION( buffer );
  face_list->Insert_ForSecondary(buffer);
  number_of_faces = face_list->NumEntities();
}
#endif

void ContactFaceBlock::Delete_Face( ContactFace<Real>* face )
{
  PRECONDITION( face );
  face_list->Delete(face);
  number_of_faces = face_list->NumEntities();
}

void ContactFaceBlock::ComputeBoundingBox(int nconfigs,
                                          VariableHandle POSITION1,
                                          VariableHandle POSITION2,
                                          MPI_Comm& SearchComm)
{
  if( face_list ){
    local_bounding_box.Reset();
    face_list->IteratorStart();
    while (ContactTopologyEntity<Real>* entity = face_list->IteratorForward()) {
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
      ContactBoundingBox object_box;
      for (int k=0; k<face->Nodes_Per_Face(); ++k) {
        ContactNode<Real>* node = face->Node(k);
        Real* position = node->Variable(POSITION1);
        object_box.add_point(position);
        if (nconfigs>1) {
          position = node->Variable(POSITION2);
          object_box.add_point(position);
        }
      }
      local_bounding_box.add_box(object_box);
    }
  }
  global_bounding_box.set_box(local_bounding_box);
  contact_global_boundingbox(global_bounding_box, SearchComm);
}

