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


#ifndef CONTACT_NO_MPI
#include "ContactEnfZoltan.h"
#include "zoltan.h"
#include "ContactEnforcement.h"
#include "ContactZoltan.h"
#include "ContactTopologyEntityHash.h"
#include "ContactNode.h"
#include "ContactEdge.h"
#include "ContactQuadFaceL4.h"
#include "ContactQuadFaceQ8.h"
#include "ContactTriFaceL3.h"
#include "ContactTriFaceQ6.h"
#include "ContactShellNode.h"
#include "ContactShellQuadFaceL4.h"
#include "ContactShellTriFaceL3.h"
#include "ContactLineFaceL2.h"
#include "ContactLineEdgeL2.h"
#include "ContactLineEdgeQ3.h"
#include "ContactHexElementL8.h"
#include "ContactAsymComm.h"
#include "ContactSymComm.h"
#include "Contact_Communication.h"
#include <iostream>
#include <cstring>


ContactEnfZoltan::ContactEnfZoltan( MPI_Comm communicator, 
				    Zoltan_Struct* Zoltan_PTR,
				    int& error )
    : num_import(0), import_global_ids(NULL), import_local_ids(NULL), import_procs(NULL), 
      num_export(0), export_global_ids(NULL), export_local_ids(NULL), export_procs(NULL)
{
  PRECONDITION( Zoltan_PTR );
  Zoltan_Ptr_ = Zoltan_PTR;
  error = ZOLTAN_OK;
}

ContactEnfZoltan::~ContactEnfZoltan(){}

void ContactEnfZoltan::Set_MultiNodeCallBacks( ContactEnforcement* enforcement )
{
  Zoltan_Set_Obj_Size_Multi_Fn( Zoltan_Ptr_,      
				ContactEnfMigrateNodeSizes,
				reinterpret_cast<void*>(enforcement) );
  
  Zoltan_Set_Pack_Obj_Multi_Fn( Zoltan_Ptr_,
				ContactEnfMigratePackNodes,
				reinterpret_cast<void*>(enforcement) );

  Zoltan_Set_Unpack_Obj_Multi_Fn( Zoltan_Ptr_,
				  ContactEnfMigrateUnpackNodes,
				  reinterpret_cast<void*>(enforcement) );
}


void ContactEnfZoltan::Set_MultiSharedNodeCallBacks( ContactEnforcement* 
						     enforcement )
{
   Zoltan_Set_Obj_Size_Multi_Fn( Zoltan_Ptr_,      
				ContactEnfMigrateSharedNodeSizes,
				reinterpret_cast<void*>(enforcement) );

  Zoltan_Set_Pack_Obj_Multi_Fn( Zoltan_Ptr_,
				ContactEnfMigratePackSharedNodes,
				reinterpret_cast<void*>(enforcement) );
  
  Zoltan_Set_Unpack_Obj_Multi_Fn( Zoltan_Ptr_,
				  ContactEnfMigrateUnpackSharedNodes,
				  reinterpret_cast<void*>(enforcement) );
}

void ContactEnfZoltan::Set_MultiFaceCallBacks( ContactEnforcement* enforcement )
{
  Zoltan_Set_Obj_Size_Multi_Fn( Zoltan_Ptr_,      
				ContactEnfMigrateFaceSizes,
				reinterpret_cast<void*>(enforcement) );

  Zoltan_Set_Pack_Obj_Multi_Fn( Zoltan_Ptr_,
				ContactEnfMigratePackFaces,
				reinterpret_cast<void*>(enforcement) );

  Zoltan_Set_Unpack_Obj_Multi_Fn( Zoltan_Ptr_,
				  ContactEnfMigrateUnpackFaces,
				  reinterpret_cast<void*>(enforcement) );
}


void 
ContactEnfZoltan::Set_MultiElementCallBacks( ContactEnforcement* enforcement )
{
  Zoltan_Set_Obj_Size_Multi_Fn( Zoltan_Ptr_,      
				ContactEnfMigrateElementSizes,
				reinterpret_cast<void*>(enforcement) );

  Zoltan_Set_Pack_Obj_Multi_Fn( Zoltan_Ptr_,
				ContactEnfMigratePackElements,
				reinterpret_cast<void*>(enforcement) );

  Zoltan_Set_Unpack_Obj_Multi_Fn( Zoltan_Ptr_,
				  ContactEnfMigrateUnpackElements,
				  reinterpret_cast<void*>(enforcement) );
}


int ContactEnfZoltan::Compute_Destinations ( int            Num_import,
					     ZOLTAN_ID_PTR  Import_global_ids,
					     ZOLTAN_ID_PTR  Import_local_ids,
					     int*           Import_procs,
					     int*           Num_export,
					     ZOLTAN_ID_PTR* Export_global_ids,
					     ZOLTAN_ID_PTR* Export_local_ids,
					     int**          Export_procs )
{
  return Zoltan_Compute_Destinations( Zoltan_Ptr_, 
				      Num_import, Import_global_ids, 
				      Import_local_ids, Import_procs, 
				      Num_export, Export_global_ids, 
				      Export_local_ids, Export_procs );
}


int ContactEnfZoltan::Help_Migrate ( int           Num_import,
                                     ZOLTAN_ID_PTR Import_global_ids,
                                     ZOLTAN_ID_PTR Import_local_ids,
                                     int*          Import_procs,
                                     int           Num_export,
                                     ZOLTAN_ID_PTR Export_global_ids,
                                     ZOLTAN_ID_PTR Export_local_ids,
                                     int*          Export_procs )
{
  return Zoltan_Help_Migrate( Zoltan_Ptr_,
			      Num_import, Import_global_ids,
			      Import_local_ids, Import_procs,
			      Num_export, Export_global_ids,
			      Export_local_ids, Export_procs );
}                                                                               


/******************************************************************************
 ******************************************************************************

                   N O D E    B A S E D   F U N C T I O N S

 ******************************************************************************
 ******************************************************************************/

void ContactEnfMigrateNodeSizes(void *data,  
				int num_gid_entries, 
                                int num_lid_entries, 
				int num_ids,
				ZOLTAN_ID_PTR gids, 
                                ZOLTAN_ID_PTR lids,
				int* sizes, 
                                int *ierr)
{
  ContactEnforcement* enforcement = (ContactEnforcement *)data;
  ContactSearch*      search      = enforcement->Search();
  ContactTopology*    topology    = search->Get_Primary_Topology(); 
  ZOLTAN_ID_PTR lid = lids;
  for( int i=0 ; i<num_ids ; ++i ){
    int entity_index  = ContactZoltanLID::Index(lid);
    ContactNode<Real>* node = static_cast<ContactNode<Real>*>
                        (topology->NodeList()->Find(entity_index));
    POSTCONDITION(node);
    sizes[i] = node->Size();
    lid += ZOLTAN_LID_SIZE;
  }
  *ierr = ZOLTAN_OK;
}

void ContactEnfMigratePackNodes(void *data,  
				int num_gid_entries, int num_lid_entries,
				int num_ids,
				ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
				int* proc, int* sizes, int* indices,
				char *Buf, int *ierr) 
{
  ContactEnforcement* enforcement = (ContactEnforcement*) data;
  ContactSearch*      search      = enforcement->Search();
  ContactTopology*    topology    = search->Get_Primary_Topology();
  ZOLTAN_ID_PTR lid = lids;
  for( int i=0 ; i<num_ids ; ++i ){
    int entity_index  = ContactZoltanLID::Index(lid);
    ContactNode<Real>* node = static_cast<ContactNode<Real>*>
                        (topology->NodeList()->Find(entity_index));
    POSTCONDITION(node);

    char* buf = Buf + indices[i];
    node->Pack(buf);
    lid += ZOLTAN_LID_SIZE;
  }
  *ierr = ZOLTAN_OK;
}

void ContactEnfMigrateUnpackNodes(void *data, int num_gid_entries,
				  int num_ids, ZOLTAN_ID_PTR gids, 
				  int* sizes, int* indices,
				  char *Buf, int *ierr)
{
  ContactEnforcement* enforcement = (ContactEnforcement *)data;
  ContactSearch*      search      = enforcement->Search();
  ContactTopology*    topology    = search->Get_Primary_Topology();
  ContactNode<Real>**       node_list   = enforcement->Phantom_Nodes();
  ContactFixedSizeAllocator* allocs = search->Get_Allocators();
  
  int num_primary_nodes = topology->Number_of_Nodes();
  int& i = enforcement->Number_Imported_Phantom_Nodes();
  for( int j=0 ; j<num_ids ; ++j ){
    char* buf = Buf + indices[j];
    int entity_type = *((int*)buf);
    switch( entity_type ) {
    case CT_NODE:
      {
	node_list[i] = ContactNode<Real>::new_ContactNode(allocs,
		         ContactSearch::NODE );
	break;
      }
    case CT_SHELL_NODE:
      {
	node_list[i] = ContactShellNode::new_ContactShellNode(allocs,
		         ContactSearch::NODE );
	break;
      }
    default:
      {
        POSTCONDITION(false);
        break;
      }
    }
    node_list[i]->ProcArrayIndex(num_primary_nodes+i);
    node_list[i]->Unpack(buf); 
    ++i;
  }
  *ierr = ZOLTAN_OK;
}


/******************************************************************************
 ******************************************************************************

                   F U N C T I O N S   F O R   S H A R E D    N O D E S

 *****************************************************************************
 ******************************************************************************/

void ContactEnfMigrateSharedNodeSizes(void *data, 
				      int num_gid_entries, 
                                      int num_lid_entries, 
				      int num_ids,
				      ZOLTAN_ID_PTR gid, 
                                      ZOLTAN_ID_PTR lid,
				      int* sizes, int *ierr)
{
  ContactEnforcement* enforcement = (ContactEnforcement *)data;
  int size = enforcement->Max_Num_Shared_Procs()*sizeof(int);
  for( int i=0 ; i<num_ids ; ++i )
    sizes[i] = size;
  *ierr = ZOLTAN_OK;
}


void ContactEnfMigratePackSharedNodes(void *data, 
				      int num_gid_entries, int num_lid_entries,
				      int num_ids,
				      ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, 
				      int* procs, int* sizes, int* indices,
				      char *Buf, int *ierr)
{
  ContactEnforcement* enforcement = (ContactEnforcement*) data;
  ContactSearch*      search      = enforcement->Search();
  ContactTopology*    topology    = search->Get_Primary_Topology();
  int max_num_shared_procs = enforcement->Max_Num_Shared_Procs();

  ZOLTAN_ID_PTR lid = lids;
  for( int i=0 ; i<num_ids ; ++i ){
    int entity_index  = ContactZoltanLID::Index(lid);
    ContactNode<Real>* node = static_cast<ContactNode<Real>*>
                        (topology->NodeList()->Find(entity_index));
    POSTCONDITION(node);
    int index = node->temp_tag;
    POSTCONDITION( 0<= index);
    int* i_buffer = reinterpret_cast<int*>(Buf+indices[i]);
    for (int j = 0; j<max_num_shared_procs; ++j) {
      i_buffer[j] = enforcement->Shared_Node_Export_Info()[index*max_num_shared_procs+j];  
    }
    lid += num_lid_entries;
  }
  *ierr = ZOLTAN_OK;
}

void ContactEnfMigrateUnpackSharedNodes(void *data, 
					int num_gid_entries, int num_ids,
					ZOLTAN_ID_PTR gids, 
					int* sizes, int* indices,
					char *Buf, int *ierr)
{
  ContactEnforcement* enforcement = (ContactEnforcement *)data;
  int max_num_shared_procs = enforcement->Max_Num_Shared_Procs();

 //
  //  Build a sorted list of objects to add
  //
  int data_to_add_size = num_ids;
  GIDSortClass *data_to_add = new GIDSortClass[data_to_add_size];
  GIDSortClass *temp_array  = new GIDSortClass[data_to_add_size];
  ZOLTAN_ID_PTR gid = gids;
  for (int i = 0; i < data_to_add_size; ++i) {
    ContactHostGlobalID global_id(gid);
    data_to_add[i].key = global_id;
    data_to_add[i].old_index = -1;
    gid += num_gid_entries;
  }
  mergeSort(data_to_add, temp_array, data_to_add_size);
  delete [] temp_array;

  //
  //  Build a sorted list of the existing array
  //
  int existing_array_size = enforcement->Get_Count_Shared_Node_Import();
  ContactHostGlobalID *shared_node_import_gid = enforcement->Get_Shared_Node_Import_Gid();
  GIDSortClass *existing_array = new GIDSortClass[existing_array_size];
                temp_array     = new GIDSortClass[existing_array_size];
  for (int i = 0; i < existing_array_size; ++i) {
    existing_array[i].key = shared_node_import_gid[i];
    existing_array[i].old_index = i;
  }
  mergeSort(existing_array, temp_array, existing_array_size);
  delete [] temp_array;



  gid = gids;
  for( int i=0 ; i<num_ids ; ++i ){
    ContactHostGlobalID global_id(gid);

    int index = enforcement->Insert_Shared_Node_Import_GID(global_id,
                                                           data_to_add, 
                                                           data_to_add_size, 
                                                           existing_array, 
                                                           existing_array_size);
    int* i_buffer = reinterpret_cast<int*>(Buf+indices[i]);
  
    for (int j = 0; j < max_num_shared_procs ; ++j) {
      int proc = i_buffer[j];
      // negative number is used to denote the end of the list
      if ( proc < 0 ) break;
      enforcement->Shared_Node_Import_Info()[index*max_num_shared_procs+j]=proc;
    }
    gid += num_gid_entries;
  }
  
  delete [] data_to_add;
  delete [] existing_array;

  *ierr = ZOLTAN_OK;
}

/******************************************************************************
 ******************************************************************************

                   F A C E    B A S E D   F U N C T I O N S

 ******************************************************************************
 ******************************************************************************/

void ContactEnfMigrateFaceSizes(void *data,  
				int num_gid_entries, 
                                int num_lid_entries, 
				int num_ids,
				ZOLTAN_ID_PTR gids, 
                                ZOLTAN_ID_PTR lids,
				int* sizes, 
                                int *ierr)
{
  ContactEnforcement* enforcement = (ContactEnforcement *)data;
  ContactSearch*      search      = enforcement->Search();
  ContactTopology*    topology    = search->Get_Primary_Topology(); 
  
  for( int i=0 ; i<num_ids ; ++i ){
    ZOLTAN_ID_PTR lid = &(lids[i*num_lid_entries]);
    int entity_index  = ContactZoltanLID::Index(lid);
    ContactFace<Real>* face = static_cast<ContactFace<Real>*>
                        (topology->FaceList()->Find(entity_index));
    POSTCONDITION(face);
    sizes[i] = face->Size();
  }
  *ierr = ZOLTAN_OK;
}

void ContactEnfMigratePackFaces(void *data, 
				int num_gid_entries, int num_lid_entries, 
				int num_ids,
				ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
				int* proc, int* sizes, int* indices,
				char *Buf, int *ierr)
{
  ContactEnforcement* enforcement = (ContactEnforcement*) data;
  ContactSearch*      search      = enforcement->Search();
  ContactTopology*    topology    = search->Get_Primary_Topology(); 
  
  for( int i=0 ; i<num_ids ; ++i ){
    ZOLTAN_ID_PTR lid = &(lids[i*num_lid_entries]);
    int entity_index  = ContactZoltanLID::Index(lid);
    ContactFace<Real>* face = static_cast<ContactFace<Real>*>
                        (topology->FaceList()->Find(entity_index));
    POSTCONDITION(face);
    char* buf = &(Buf[indices[i]]);
    face->Pack(buf);
  }
  *ierr = ZOLTAN_OK;
}


void ContactEnfMigrateUnpackFaces(void *data, int num_gid_entries, 
				  int num_ids, ZOLTAN_ID_PTR gid, 
				  int* sizes, int* indices,
				  char *Buf, int *ierr)
{
  ContactEnforcement* enforcement   = (ContactEnforcement *)data;
  ContactSearch*      search        = enforcement->Search();
  ContactTopology*    topology      = search->Get_Primary_Topology();
  ContactFace<Real>**       face_list     = enforcement->Phantom_Faces();
  ContactFixedSizeAllocator* allocs = search->Get_Allocators();
  
  int num_primary_faces = topology->Number_of_Faces();
  int& i = enforcement->Number_Imported_Phantom_Faces();
  for( int j=0 ; j<num_ids ; ++j ){
    char* buf = &(Buf[indices[j]]);
    PRECONDITION( (reinterpret_cast<int*>(buf))[0] == CT_FACE );
    ContactSearch::ContactFace_Type face_type = (ContactSearch::ContactFace_Type)(reinterpret_cast<int*>(buf))[1];

    face_list[i] = ContactSearch::New_ContactFace(face_type, allocs);
    face_list[i]->ProcArrayIndex(num_primary_faces+i); 
    face_list[i]->Unpack(buf); 
    i++;
  }
  *ierr = ZOLTAN_OK;
}

/******************************************************************************
 ******************************************************************************

                 E L E M E N T    B A S E D   F U N C T I O N S

 ******************************************************************************
 ******************************************************************************/

void ContactEnfMigrateElementSizes(void *data,  
				   int num_gid_entries, 
                                   int num_lid_entries, 
				   int num_ids,
				   ZOLTAN_ID_PTR gids, 
                                   ZOLTAN_ID_PTR lids,
				   int* sizes, 
                                   int *ierr)
{
  ContactEnforcement* enforcement    = (ContactEnforcement *)data;
  ContactSearch*      search         = enforcement->Search();
  ContactTopology*    topology       = search->Get_Primary_Topology(); 

  for( int i=0 ; i<num_ids ; ++i ){
    ZOLTAN_ID_PTR lid = &(lids[i*num_lid_entries]);
    int entity_index  = ContactZoltanLID::Index(lid);
    ContactElement* element = static_cast<ContactElement*>
                        (topology->ElemList()->Find(entity_index));
    POSTCONDITION(element);
    sizes[i] = element->Size();
  }
  *ierr = ZOLTAN_OK;
}

void ContactEnfMigratePackElements(void *data, 
				   int num_gid_entries, int num_lid_entries, 
				   int num_ids,
				   ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
				   int* proc, int* sizes, int* indices,
				   char *Buf, int *ierr)
{
  ContactEnforcement* enforcement    = (ContactEnforcement*) data;
  ContactSearch*      search         = enforcement->Search();
  ContactTopology*    topology       = search->Get_Primary_Topology(); 

  for( int i=0 ; i<num_ids ; ++i ){
    ZOLTAN_ID_PTR lid = &(lids[i*num_lid_entries]);
    int entity_index  = ContactZoltanLID::Index(lid);
    ContactElement* element = static_cast<ContactElement*>
                        (topology->ElemList()->Find(entity_index));
    POSTCONDITION(element);
    char* buf = &(Buf[indices[i]]);
    element->Pack(buf);
  }
  *ierr = ZOLTAN_OK;
}

void ContactEnfMigrateUnpackElements(void *data, int num_gid_entries,
				     int num_ids, ZOLTAN_ID_PTR gid,
				     int* sizes, int* indices,
				     char *Buf, int *ierr)
{
  ContactEnforcement* enforcement   = (ContactEnforcement *)data;
  ContactSearch*      search        = enforcement->Search();
  ContactTopology*    topology      = search->Get_Primary_Topology();
  ContactElement**    element_list  = enforcement->Phantom_Elements();
  ContactFixedSizeAllocator* allocs = search->Get_Allocators();
  
  int num_primary_elements = topology->Number_of_Elements();
  int& i = enforcement->Number_Imported_Phantom_Elements();

  for( int j=0 ; j<num_ids ; ++j ){
    char* buf = &(Buf[indices[j]]);
    PRECONDITION((reinterpret_cast<int*>(buf))[0]==CT_ELEMENT);
    int element_type = (reinterpret_cast<int*>(buf))[1];
    switch( element_type ){
    case ContactSearch::CARTESIANHEXELEMENTL8 :
      element_list[i] = ContactCartesianHexElementL8::
	new_ContactCartesianHexElementL8(allocs);
      break;
    case ContactSearch::HEXELEMENTL8 :
      element_list[i] = ContactHexElementL8::new_ContactHexElementL8(allocs);
      break;
    }
    element_list[i]->ProcArrayIndex(num_primary_elements+i); 
    element_list[i]->Unpack(buf); 
    i++;
  }
  *ierr = ZOLTAN_OK;
}

#endif
