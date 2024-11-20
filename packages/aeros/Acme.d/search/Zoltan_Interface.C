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

#include "zoltan.h"
#include "Zoltan_Interface.h"
#include "ContactSearch.h"
#include "ContactTopology.h"
#include "ContactZoltan.h"
#include "ContactZoltanID.h"
#include "ContactHostGlobalID.h"
#include "ContactTopologyEntityHash.h"
#include "ContactNode.h"
#include "ContactEdge.h"
#include "ContactElementBlock.h"
#include "ContactFaceBlock.h"
#include "ContactEdgeBlock.h"
#include "ContactNodeBlock.h"
#include "ContactHexElementL8.h"
#include "ContactQuadFaceL4.h"
#include "ContactQuadFaceQ8.h"
#include "ContactQuadFaceQ9.h"
#include "ContactTriFaceL3.h"
#include "ContactTriFaceQ6.h"
#include "ContactShellNode.h"
#include "ContactShellQuadFaceL4.h"
#include "ContactShellTriFaceL3.h"
#include "ContactLineFaceL2.h"
#include "ContactLineEdgeL2.h"
#include "ContactLineEdgeQ3.h"
#include "ContactNodeFaceInteraction.h"
#include "ContactNodeSurfaceInteraction.h"
#include "ContactFaceFaceInteraction.h"
#include "ContactFaceCoverageInteraction.h"
#include "ContactElementElementInteraction.h"
#include "ContactZoltanComm.h"
#include "Contact_Communication.h"
#include "ContactFixedSizeAllocator.h"
#include <iostream>
#include <cstdio>
#include <cstring>
#include <map>

using namespace std;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//      The following three routines are for the RCB decomposition.
//
//-------------------------------------------------------------------------

int ContactQueryNumObjects(void *data, int *ierr)
{
  ContactSearch*   search   = (ContactSearch *)data;
  ContactTopology* topology = search->Get_Primary_Topology(); 
  int number_of_nodes = topology->Number_of_Nodes();
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList()); 
   
#if CONTACT_DEBUG_PRINT_LEVEL>=3
  ContactParOStream& postream = search->ParOStream();
  postream<<"    Step "<<search->StepNumber()<<", RCB -----------------------------------\n";
#endif
  int NumberOfLocalObjects = 0;
  for (int i=0; i<number_of_nodes; ++i) {
    if (Nodes[i]->IsGhostRCB()) {
      ++NumberOfLocalObjects;
      //postream<<"      tagging node "<<Nodes[i]->Exodus_ID()<<"\n";
    }
  }
  *ierr = ZOLTAN_OK;
#if CONTACT_DEBUG_PRINT_LEVEL>=3
  postream<<"      Tagged "<<NumberOfLocalObjects<<" nodes for RCB\n";
#endif
  return NumberOfLocalObjects;
}

void ContactQueryObjectList(void *data, 
                            int num_gid_entries, int num_lid_entries, 
                            ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
                            int wdim, float *objwgts, int *ierr)
{
  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  PRECONDITION(num_lid_entries==ZOLTAN_LID_SIZE);
  ContactSearch*   search   = (ContactSearch *)data;
  ContactTopology* topology = search->Get_Primary_Topology(); 
  int number_of_nodes = topology->Number_of_Nodes();
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());  
  
  int k=0;
  float *wts = objwgts;
  ZOLTAN_ID_PTR lid = lids;
  ZOLTAN_ID_PTR gid = gids;
  if (wdim!=0) {
    for (int i=0; i<number_of_nodes; ++i) {
      if (Nodes[i]->IsGhostRCB()) {
        Nodes[i]->ZoltanLID(CT_NODE, lid);
        Nodes[i]->ZoltanGID(CT_NODE, gid);
        for (int j=0; j<wdim; ++j) {
          *wts++ = 1.0;
        }
        lid += num_lid_entries;
        gid += num_gid_entries;
        ++k;
      }
    }
  } else {
    for (int i=0; i<number_of_nodes; ++i) {
      if (Nodes[i]->IsGhostRCB()) {
        Nodes[i]->ZoltanLID(CT_NODE, lid);
        Nodes[i]->ZoltanGID(CT_NODE, gid);
        lid += num_lid_entries;
        gid += num_gid_entries;
        ++k;
      }
    }
  }
  *ierr = ZOLTAN_OK;
}     

int ContactQueryNumGeomObjects(void *data, int *ierr)
{
  *ierr = ZOLTAN_OK;
  ContactSearch* search = (ContactSearch *)data;
  return search->Get_Primary_Topology()->Dimensionality();
}

void ContactQueryGeomMultiValues(void *data, 
                                 int num_gid_entries, 
                                 int num_lid_entries,
                                 int num_obj, 
                                 ZOLTAN_ID_PTR gids, 
                                 ZOLTAN_ID_PTR lids,
                                 int num_dim, 
                                 double *geom_vec, 
                                 int *ierr)
{
  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  PRECONDITION(num_lid_entries==ZOLTAN_LID_SIZE);
  ContactSearch*   search   = (ContactSearch *)data;
  ContactTopology* topology = search->Get_Primary_Topology();
  ContactZoltan*   zoltan   = search->Get_Zoltan();   
  VariableHandle   POSITION = zoltan->Position();
  ContactNode<Real>**    Nodes    = 
    reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());
  
  PRECONDITION(num_dim==topology->Dimensionality());

  double *zp = geom_vec;
  ZOLTAN_ID_PTR lid = lids;
  for (int i=0; i<num_obj; ++i) {
    int entity_index = ContactZoltanLID::Index(lid);
    PRECONDITION(entity_index>=0);
    PRECONDITION(entity_index<topology->Number_of_Nodes());
    POSTCONDITION(Nodes[entity_index]);
    Real* position = Nodes[entity_index]->Variable(POSITION);
    for (int j=0; j<num_dim; ++j) {
      *zp++ = *position++;
    }
    lid += num_lid_entries;
  }
  *ierr = ZOLTAN_OK;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//      The following three routines are for migrating topological
//      entities from the primary topology to the secondary topology.
//
//-------------------------------------------------------------------------

void ContactMigrateEntityExportSizes(void *data,  
			             int num_gid_entries, 
                                     int num_lid_entries, 
			             int num_ids,
			             ZOLTAN_ID_PTR gids, 
                                     ZOLTAN_ID_PTR lids,
			             int* sizes,
			             int *ierr )
{
  ContactSearch*   search = (ContactSearch *)data;

  #if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
  Real bytes_nodes = search->Bytes_For_Nodes();
  Real bytes_faces = search->Bytes_For_Faces();
  Real bytes_elems = search->Bytes_For_Elems();
  #endif

  #ifdef CONTACT_TIMINGS
  search->Timer()->Start_Timer( search->Owner_help_migrate_size_time() );
  #endif
  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  PRECONDITION(num_lid_entries==ZOLTAN_LID_SIZE);
  ContactTopology* topology = search->Get_Primary_Topology();
  ContactNode<Real>** nodes = 
    reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());
  ContactFace<Real>** faces = 
    reinterpret_cast<ContactFace<Real>**>(topology->FaceList()->EntityList());
  ContactElement** elements = 
    reinterpret_cast<ContactElement**>(topology->ElemList()->EntityList());

  int edgebased_physical_face  = (search->PhysicalFaceAlgorithm()) == ContactSearch::PF_EDGE_BASED;
  int enable_off_face_tracking = search->OffFaceTrackingStatus();

  ZOLTAN_ID_PTR lid = lids;
  for( int i=0 ; i<num_ids ; ++i ){
    int entity_type   = ContactZoltanLID::Type(lid);
    int entity_index  = ContactZoltanLID::Index(lid);
    int& size         = sizes[i];
    switch (entity_type) {
    case CT_NODE:
    case CT_SHELL_NODE:
      {
	ContactNode<Real>* node = nodes[entity_index];
        POSTCONDITION(node);
        int state = (node->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE))?1:-2;
        #ifdef CONTACT_OLD_XFER
	size = node->Size(state);
        #else
	size = node->Size_ForSecondary(state);
        #endif
        #if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
        bytes_nodes += size;
        #endif
	break;
      }
    case CT_EDGE:
      POSTCONDITION(false);
      break;
    case CT_FACE:
      {
	ContactFace<Real>* face = faces[entity_index];
        POSTCONDITION(face);
        #ifdef CONTACT_OLD_XFER
	size = face->Size(enable_off_face_tracking);
        #else
	size = face->Size_ForSecondary(enable_off_face_tracking,
                                       edgebased_physical_face);
        #endif
        #if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
        bytes_faces += size;
        #endif
	break;
      }
    case CT_ELEMENT:
      {
	ContactElement* element = elements[entity_index];
        POSTCONDITION(element);
        #ifdef CONTACT_OLD_XFER
	size = element->Size();
        #else
	size = element->Size_ForSecondary();
        #endif
        #if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
        bytes_elems += size;
        #endif
	break;
      }
    default:
      POSTCONDITION(false);
      break;
    }
    lid += num_lid_entries;
  }
  *ierr = ZOLTAN_OK;
  #ifdef CONTACT_TIMINGS
  search->Timer()->Stop_Timer( search->Owner_help_migrate_size_time() );
  #endif
  #if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
  search->Bytes_For_Nodes(bytes_nodes);
  search->Bytes_For_Faces(bytes_faces);
  search->Bytes_For_Elems(bytes_elems);
  #endif
}

void ContactMigrateEntityExportPack(void *data, 
			            int num_gid_entries, 
                                    int num_lid_entries,
			            int num_ids,
			            ZOLTAN_ID_PTR gids, 
                                    ZOLTAN_ID_PTR lids,
			            int* procs, int* sizes, 
                                    int* indices, char *Buf, 
                                    int *ierr)
{
  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  PRECONDITION(num_lid_entries==ZOLTAN_LID_SIZE);
  ContactSearch*   search   = (ContactSearch *)data;
  ContactTopology* topology = search->Get_Primary_Topology();
  ContactNode<Real>** nodes = 
    reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());
  ContactFace<Real>** faces = 
    reinterpret_cast<ContactFace<Real>**>(topology->FaceList()->EntityList());
  ContactElement** elements = 
    reinterpret_cast<ContactElement**>(topology->ElemList()->EntityList());
  #ifdef CONTACT_TIMINGS
  search->Timer()->Start_Timer( search->Owner_help_migrate_pack_time() );
  #endif

  int edgebased_physical_face  = (search->PhysicalFaceAlgorithm()) == ContactSearch::PF_EDGE_BASED;
  int enable_off_face_tracking = search->OffFaceTrackingStatus();
  
  //ContactParOStream& postream = search->ParOStream();
  
  ZOLTAN_ID_PTR lid = lids;
  for( int i=0 ; i<num_ids ; ++i ){
    int entity_type  = ContactZoltanLID::Type(lid);
    int entity_index = ContactZoltanLID::Index(lid);
    char* buf        = Buf+indices[i];
    switch (entity_type) {
    case CT_NODE:
    case CT_SHELL_NODE:
      {
        ContactNode<Real>* node = nodes[entity_index];
	POSTCONDITION(node);
        int state = (node->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE))?1:-2;
        #ifdef CONTACT_OLD_XFER
        node->Pack(buf, state);
        #else
        node->Pack_ForSecondary(buf, state);
        #endif
      }
      break;
    case CT_EDGE:
      POSTCONDITION(false);
      break;
    case CT_FACE:
      {
        ContactFace<Real>* face = faces[entity_index];
	POSTCONDITION(face);
	#ifdef CONTACT_OLD_XFER
	face->Pack(buf,enable_off_face_tracking);
	#else
	face->Pack_ForSecondary(buf,enable_off_face_tracking,
                                edgebased_physical_face);
	#endif
      }
      break;
    case CT_ELEMENT:
      {
        ContactElement* element = elements[entity_index];
	POSTCONDITION(element);
	#ifdef CONTACT_OLD_XFER
	element->Pack(buf);
	#else
	element->Pack_ForSecondary(buf);
	#endif
      }
      break;
    default:
      POSTCONDITION(false);
      break;
    }
    lid += num_lid_entries;
  }
  *ierr = ZOLTAN_OK;
  #ifdef CONTACT_TIMINGS
  search->Timer()->Stop_Timer( search->Owner_help_migrate_pack_time() );
  #endif

}

void ContactMigrateEntityUnpack(void *data, 
                                int num_gid_entries,
			        int num_ids, 
                                ZOLTAN_ID_PTR gids, 
			        int* sizes,  
                                int* indices,
			        char *Buf, 
                                int *ierr)
{
  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  ContactSearch*   search   = (ContactSearch *)data;
  ContactTopology* topology = search->Get_Secondary_Topology(); 
  #ifdef CONTACT_TIMINGS
  search->Timer()->Start_Timer( search->Owner_help_migrate_unpack_time() );
  #endif

  //ContactParOStream& postream = search->ParOStream();
  for( int i=0 ; i<num_ids ; ++i ){
    char* buf       = Buf + indices[i];
    int*  ibuf      = reinterpret_cast<int*>(buf);
    int entity_type = ibuf[0];
    int block       = ibuf[ContactTopologyEntity<Real>::BLOCK_ID];
    REMEMBER(ContactHostGlobalID GID(&gids[i*num_gid_entries]));
    switch (entity_type) {
    case CT_NODE:
    case CT_SHELL_NODE:
      PRECONDITION( block>=0 && block<topology->Number_of_Node_Blocks() );
      #ifdef CONTACT_OLD_XFER
      topology->Node_Block(block)->Insert_Node(buf);
      #else
      topology->Node_Block(block)->Insert_Node_ForSecondary(buf);
      #endif
      break;
    case CT_EDGE:
      POSTCONDITION(false);
      break;
    case CT_FACE: 
      PRECONDITION( block>=0 && block<topology->Number_of_Face_Blocks() );
      #ifdef CONTACT_OLD_XFER
      topology->Face_Block(block)->Insert_Face(buf);
      #else
      topology->Face_Block(block)->Insert_Face_ForSecondary(buf);
      #endif
      break;
    case CT_ELEMENT:
      PRECONDITION( block>=0 && block<topology->Number_of_Element_Blocks() );
      #ifdef CONTACT_OLD_XFER
      topology->Element_Block(block)->Insert_Element(buf);
      #else
      topology->Element_Block(block)->Insert_Element_ForSecondary(buf);
      #endif
      break;
    default:
      POSTCONDITION(false);
      break;
    }
  }
  *ierr = ZOLTAN_OK;
  #ifdef CONTACT_TIMINGS
  search->Timer()->Stop_Timer( search->Owner_help_migrate_unpack_time() );
  #endif
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//      The following three routines are for migrating the interactions
//      from the secondary topology back to the primary topology.
//
//-------------------------------------------------------------------------

void ContactMigrateInteractionSizes(void *data, 
				    int num_gid_entries, 
                                    int num_lid_entries,
				    int num_ids,
				    ZOLTAN_ID_PTR gids, 
                                    ZOLTAN_ID_PTR lids,
				    int* sizes, 
                                    int *ierr)
{
  ContactSearch* search = (ContactSearch *)data;

#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
  Real bytes_nodes = search->Bytes_For_Nodes();
  Real bytes_faces = search->Bytes_For_Faces();
  Real bytes_elems = search->Bytes_For_Elems();
#endif

#ifdef CONTACT_TIMINGS
  search->Timer()->Start_Timer( search->Interaction_help_migrate_size_time() );
#endif
  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  PRECONDITION(num_lid_entries==ZOLTAN_LID_SIZE);
  ContactTopology* topology = search->Get_Secondary_Topology();
  ContactNode<Real>** nodes = 
    reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());
  ContactFace<Real>** faces = 
    reinterpret_cast<ContactFace<Real>**>(topology->FaceList()->EntityList());
  ContactElement** elements = 
    reinterpret_cast<ContactElement**>(topology->ElemList()->EntityList());
  
  int state0 = 0;
  ZOLTAN_ID_PTR lid = lids;
  for( int i=0 ; i<num_ids ; ++i ){
    int entity_type  = ContactZoltanLID::Type(lid);
    int entity_index = ContactZoltanLID::Index(lid);
    switch (entity_type) {
    case CT_NODE:
    case CT_SHELL_NODE:
      {
        ContactNode<Real>* node = nodes[entity_index];
	POSTCONDITION(node);
	sizes[i] = node->Size_Interactions(state0);
#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
        bytes_nodes += sizes[i];
#endif
      }
      break;
    case CT_EDGE:
      POSTCONDITION(false);
      break;
    case CT_FACE:
      {
        ContactFace<Real>* face = faces[entity_index];
	POSTCONDITION(face);
	sizes[i] = face->Size_Interactions(state0);
#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
        bytes_faces += sizes[i];
#endif
      }
      break;
    case CT_ELEMENT:
      {
        ContactElement* element = elements[entity_index];
	POSTCONDITION(element);
	sizes[i] = element->Size_Interactions(state0);
#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
        bytes_elems += sizes[i];
#endif
      }
      break;
    default:
      POSTCONDITION(false);
      break;
    }
    lid += num_lid_entries;
  }
  *ierr = ZOLTAN_OK;
#ifdef CONTACT_TIMINGS
  search->Timer()->Stop_Timer( search->Interaction_help_migrate_size_time() );
#endif
#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
  search->Bytes_For_Nodes(bytes_nodes);
  search->Bytes_For_Faces(bytes_faces);
  search->Bytes_For_Elems(bytes_elems);
#endif
}

void ContactMigratePackInteractions(void *data, 
				    int num_gid_entries, int num_lid_entries,
				    int num_ids,
				    ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
				    int* proc, int* sizes, int* indices,
				    char *Buf, int *ierr)
{
  ContactSearch* search = (ContactSearch *)data;
#ifdef CONTACT_TIMINGS
  search->Timer()->Start_Timer( search->Interaction_help_migrate_pack_time() );
#endif
  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  PRECONDITION(num_lid_entries==ZOLTAN_LID_SIZE);
  ContactTopology* topology = search->Get_Secondary_Topology();
  ContactNode<Real>** nodes = 
    reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());
  ContactFace<Real>** faces = 
    reinterpret_cast<ContactFace<Real>**>(topology->FaceList()->EntityList());
  ContactElement** elements = 
    reinterpret_cast<ContactElement**>(topology->ElemList()->EntityList());
  
  int state0 = 0;
  ZOLTAN_ID_PTR lid = lids;
  for( int i=0 ; i<num_ids ; ++i ){
    char* buf        = &(Buf[indices[i]]);
    int entity_type  = ContactZoltanLID::Type(lid);
    int entity_index = ContactZoltanLID::Index(lid);
    switch (entity_type) {
    case CT_NODE:
    case CT_SHELL_NODE:
      {
        ContactNode<Real>* node = nodes[entity_index];
	POSTCONDITION(node);
	node->Pack_Interactions( buf, state0 );
      }
      break;
    case CT_FACE:
      {
        ContactFace<Real>* face = faces[entity_index];
	POSTCONDITION(face);
	face->Pack_Interactions( buf, state0 );
      }
      break;
    case CT_ELEMENT:
      {
        ContactElement* element = elements[entity_index];
	POSTCONDITION(element);
	element->Pack_Interactions( buf, state0 );
      }
      break;
    default:
      POSTCONDITION(false);
      break;
    }
    lid += num_lid_entries;
  }
  *ierr = ZOLTAN_OK;
#ifdef CONTACT_TIMINGS
  search->Timer()->Stop_Timer( search->Interaction_help_migrate_pack_time() );
#endif
}

void ContactMigrateUnpackInteractions(void *data, int num_gid_entries,
				      int num_ids, ZOLTAN_ID_PTR gids, 
				      int* sizes, int* indices,
				      char *Buf, int *ierr)
{
  ContactSearch* search = (ContactSearch *)data;
#ifdef CONTACT_TIMINGS
  search->Timer()->Start_Timer( search->Interaction_help_migrate_unpack_time() );
#endif
  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  ContactTopology* dest_topology = search->Get_Primary_Topology();

  ZOLTAN_ID_PTR gid = gids;
  for( int i=0 ; i<num_ids ; ++i ){
    ContactInteractionEntity<Real>* interaction;
    char* buf       = &(Buf[indices[i]]);
    int entity_type = ContactZoltanGID::Type(gid);
    ContactHostGlobalID GID(gid);
    switch (entity_type) {
    case CT_NODE:
    case CT_SHELL_NODE:
      {
        int j;
        ContactNode<Real>* node = static_cast<ContactNode<Real>*>
                            (dest_topology->NodeList()->Find(GID));
	POSTCONDITION( node );
	node->Unpack_Interactions(buf, 0);
	if (node->Number_Interactions()>0) {
          ContactNodeEntityInteraction** interactions = 
            node->Get_NodeEntity_Interactions();
          for (j=0; j<node->Number_NodeEntity_Interactions(); ++j) {
            switch (interactions[j]->Get_Type()) {
            case ContactNodeEntityInteraction::NODE_FACE_INTERACTION: {
              ContactNodeFaceInteraction*  cnfi = 
                static_cast<ContactNodeFaceInteraction*>(interactions[j]);
              cnfi->Connect_Face( dest_topology );
              } break;
            case ContactNodeEntityInteraction::NODE_SURFACE_INTERACTION: {
              ContactNodeSurfaceInteraction* cnsi = 
                static_cast<ContactNodeSurfaceInteraction*>(interactions[j]);
              cnsi->Connect_Surface( dest_topology );
              } break;
            }
          }
        }
      }
      break;
    case CT_FACE:
      {
        ContactFace<Real>* face = static_cast<ContactFace<Real>*>
                            (dest_topology->FaceList()->Find(GID));
	POSTCONDITION( face );
	face->Unpack_Interactions(buf);
        ContactInteractionDLL<Real>* interactions = face->Get_FaceFace_Interactions();
        if(interactions == NULL) continue;
        interactions->IteratorStart();
        while ((interaction=interactions->IteratorForward())){
          ContactFaceFaceInteraction<Real>* cffi = 
            static_cast<ContactFaceFaceInteraction<Real>*>(interaction);
          cffi->Connect_MasterFace( dest_topology );
        }
      }
      break;
    case CT_ELEMENT:
      {
        ContactElement* element = static_cast<ContactElement*>
                            (dest_topology->ElemList()->Find(GID));
	POSTCONDITION( element );
	element->Unpack_Interactions(buf);
        ContactInteractionDLL<Real>* interactions = element->Get_ElementElement_Interactions();
        interactions->IteratorStart();
        while ((interaction=interactions->IteratorForward())){
          ContactElementElementInteraction* ceei = 
            static_cast<ContactElementElementInteraction*>(interaction);
          ceei->Connect_MasterElement( dest_topology );
        }
      }
      break;
    default:
      std::cerr<<"Entity Type = "<<entity_type<<std::endl;
      POSTCONDITION(false);
      break;
    }
    gid += num_gid_entries;
  }
  *ierr = ZOLTAN_OK;
#ifdef CONTACT_TIMINGS
  search->Timer()->Stop_Timer( search->Interaction_help_migrate_unpack_time() );
#endif
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//      The following three routines are for querying other processors
//      to update references to entities on that processor from other
//      processors.
//
//-------------------------------------------------------------------------

int ContactHostidQuerySize(void *data,  
                           int num_gid_entries, 
                           int num_lid_entries, 
                           ZOLTAN_ID_PTR gid, 
                           ZOLTAN_ID_PTR lid,
                           int *ierr)
{
  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  PRECONDITION(num_lid_entries==ZOLTAN_LID_SIZE);
  int size = sizeof(ContactInteractionEntity<Real>::entity_data);
  *ierr = ZOLTAN_OK;
  return size;
}

void ContactHostidQueryPack(void *data, 
                            int num_gid_entries, int num_lid_entries, 
                            ZOLTAN_ID_PTR gid, ZOLTAN_ID_PTR lid,
		            int proc, int size, char *buf, int *ierr)
{
  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  PRECONDITION(num_lid_entries==ZOLTAN_LID_SIZE);
  ContactSearch*   search   = (ContactSearch *)data;
  ContactTopology* topology = search->Get_Primary_Topology();
  
  int* i_buffer = reinterpret_cast<int*>(buf);
  
  int entity_type = ContactZoltanGID::Type(gid);
  ContactHostGlobalID global_id( ContactZoltanGID::Hi(gid), 
                                 ContactZoltanGID::Lo(gid) );
  
  switch (entity_type) {
  case CT_NODE:
    {
      ContactNode<Real>* node = static_cast<ContactNode<Real>*>
                           (topology->NodeList()->Find(global_id));
      *i_buffer++ = CT_NODE;
      if (node) {
        *i_buffer++ = node->Global_ID().HiInt();
        *i_buffer++ = node->Global_ID().LoInt();
        *i_buffer++ = node->Owner(); 
        *i_buffer++ = node->BlockID(); 
        *i_buffer++ = node->HostGlobalArrayIndex(); 
        //*i_buffer++ = node->ProcArrayIndex();   
        *i_buffer++ = node->OwnerProcArrayIndex();
      } else {
        *i_buffer++ = -1;
        *i_buffer++ = -1;
        *i_buffer++ = -1;
        *i_buffer++ = -1;
        *i_buffer++ = -1;
        //*i_buffer++ = -1;
        *i_buffer++ = -1;
      }
    }
    break;
  case CT_FACE:
    {
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>
                           (topology->FaceList()->Find(global_id));
      *i_buffer++ = CT_FACE;
      if (face) {
        *i_buffer++ = face->Global_ID().HiInt();
        *i_buffer++ = face->Global_ID().LoInt();
        *i_buffer++ = face->Owner(); 
        *i_buffer++ = face->BlockID(); 
        *i_buffer++ = face->HostGlobalArrayIndex(); 
        //*i_buffer++ = face->ProcArrayIndex();   
        *i_buffer++ = face->OwnerProcArrayIndex();
      } else {
        *i_buffer++ = -1;
        *i_buffer++ = -1;
        *i_buffer++ = -1;
        *i_buffer++ = -1;
        *i_buffer++ = -1;
        //*i_buffer++ = -1;
        *i_buffer++ = -1;
      }
    }
    break;
  case CT_ELEMENT:
    {
      ContactElement* element = static_cast<ContactElement*>
                                (topology->ElemList()->Find(global_id));
      *i_buffer++ = CT_ELEMENT;
      if (element) {
        *i_buffer++ = element->Global_ID().HiInt();
        *i_buffer++ = element->Global_ID().LoInt();
        *i_buffer++ = element->Owner(); 
        *i_buffer++ = element->BlockID(); 
        *i_buffer++ = element->HostGlobalArrayIndex(); 
        //*i_buffer++ = element->ProcArrayIndex();   
        *i_buffer++ = element->OwnerProcArrayIndex();
      } else {
        *i_buffer++ = -1;
        *i_buffer++ = -1;
        *i_buffer++ = -1;
        *i_buffer++ = -1;
        *i_buffer++ = -1;
        //*i_buffer++ = -1;
        *i_buffer++ = -1;
      }
    }
    break;
  default:
    POSTCONDITION(false);
    break;
  }
  POSTCONDITION((char*)i_buffer<=buf+size);
  *ierr = ZOLTAN_OK;
}

void ContactHostidQueryUnpack(void *data, int num_gid_entries,
                              ZOLTAN_ID_PTR gid, int size,  
                              char *buf, int *ierr)
{
  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  int* i_buffer = reinterpret_cast<int*>(buf);
  if (i_buffer[1] != -1 &&   //entity_data->HostGlobalID[0]
      i_buffer[2] != -1 ) {  //entity_data->HostGlobalID[1]
    ContactSearch*   search   = (ContactSearch *)data;
    ContactTopology* topology = search->Get_Primary_Topology();

    std::vector<ContactInteractionEntity<Real>::entity_data *> *query_linklist = topology->QueryLinkList();
    ContactInteractionEntity<Real>::entity_data* entity_data = 
        new ContactInteractionEntity<Real>::entity_data[1];
    entity_data->type                      = *i_buffer++;
    entity_data->host_gid[0]               = *i_buffer++;
    entity_data->host_gid[1]               = *i_buffer++;
    entity_data->owner                     = *i_buffer++;
    entity_data->block_id                  = *i_buffer++;
    entity_data->index_in_host_array       = *i_buffer++;
    //entity_data->index_in_proc_array       = *i_buffer++;
    entity_data->index_in_owner_proc_array = *i_buffer++;
    query_linklist->push_back( entity_data );
    POSTCONDITION((char*)i_buffer<=buf+size);
  }
  *ierr = ZOLTAN_OK;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//      The following three routines are for migrating topological
//      entities for incremental dynamic load balancing.
//
//-------------------------------------------------------------------------

void ContactDynamicLoadBalanceSize(void *data,  
                                   int num_gid_entries, 
                                   int num_lid_entries,
			           int num_ids, 
                                   ZOLTAN_ID_PTR gids, 
                                   ZOLTAN_ID_PTR lids,
			           int* sizes,
                                   int *ierr)
{
  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  PRECONDITION(num_lid_entries==ZOLTAN_LID_SIZE);
  ContactSearch*   search   = (ContactSearch *)data;
  ContactTopology* topology = search->Get_Primary_Topology();

  int state1 = 1;
  ZOLTAN_ID_PTR gid = gids;
  for( int i=0 ; i<num_ids ; ++i ){
    int entity_type = ContactZoltanGID::Type(gid);
    int& size       = sizes[i];
    ContactHostGlobalID GID(gid);
    switch (entity_type) {
    case CT_NODE:
    case CT_SHELL_NODE:
      {
	ContactNode<Real>* node = NULL;
        for(int j=0 ; j<topology->Number_of_Node_Blocks() ; ++j ){
          node = static_cast<ContactNode<Real>*>
                 (topology->Node_Block(j)->NodeList()->Find(GID));
          if ( node ) break;
        }
        POSTCONDITION(node);
	size = node->Size(state1);
	break;
      }
    case CT_FACE:
      {
        ContactFace<Real>* face = NULL;
        for(int j=0 ; j<topology->Number_of_Face_Blocks() ; ++j ){
          face = static_cast<ContactFace<Real>*>
                 (topology->Face_Block(j)->FaceList()->Find(GID));
          if ( face ) break;
        }
        POSTCONDITION(face);
	size = face->Size(1);
	break;
      }
    case CT_ELEMENT:
      {
        ContactElement* element = NULL;
        for(int j=0 ; j<topology->Number_of_Element_Blocks() ; ++j ){
          element = static_cast<ContactElement*>
                    (topology->Element_Block(j)->ElemList()->Find(GID));
          if ( element ) break;
        }
        POSTCONDITION(element);
	size = element->Size();
	break;
      }
    default:
      POSTCONDITION(false);
      break;
    }
    gid += num_gid_entries;
  }
  *ierr = ZOLTAN_OK;
}

void ContactDynamicLoadBalancePack(void *data, 
                                   int num_gid_entries, int num_lid_entries,
			           int num_ids, 
                                   ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
			           int* procs, int* sizes, int* indices,
			           char *Buf, int *ierr)
{
  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  PRECONDITION(num_lid_entries==ZOLTAN_LID_SIZE);
  ContactSearch*   search   = (ContactSearch *)data;
  ContactTopology* topology = search->Get_Primary_Topology();
  
  ZOLTAN_ID_PTR gid = gids;
  for( int i=0 ; i<num_ids ; ++i ){
    int entity_type = ContactZoltanGID::Type(gid);
    char* buf       = Buf+indices[i];
    ContactHostGlobalID GID(gid);
    switch (entity_type) {
    case CT_NODE:
    case CT_SHELL_NODE:
      {
	ContactNode<Real>* node = NULL;
        for(int j=0 ; j<topology->Number_of_Node_Blocks() ; ++j ){
          node = static_cast<ContactNode<Real>*>
                 (topology->Node_Block(j)->NodeList()->Find(GID));
          if ( node ) break;
        }
	POSTCONDITION(node);
        int state = (node->Secondary_Owner()==procs[i])?1:-2;
        node->Pack(buf, state);
      }
      break;
    case CT_FACE:
      {
        ContactFace<Real>* face = NULL;
        for(int j=0 ; j<topology->Number_of_Face_Blocks() ; ++j ){
          face = static_cast<ContactFace<Real>*>
                 (topology->Face_Block(j)->FaceList()->Find(GID));
          if ( face ) break;
        }
	POSTCONDITION(face);
	face->Pack(buf,1);
      }
      break;
    case CT_ELEMENT:
      {
        ContactElement* element = NULL;
        for(int j=0 ; j<topology->Number_of_Element_Blocks() ; ++j ){
          element = static_cast<ContactElement*>
                    (topology->Element_Block(j)->ElemList()->Find(GID));
          if ( element ) break;
        }
	POSTCONDITION(element);
	element->Pack(buf);
      }
      break;
    default:
      POSTCONDITION(false);
      break;
    }
    gid += num_gid_entries;
  }
  *ierr = ZOLTAN_OK;
}

void ContactDynamicLoadBalanceUnpack(void *data, int num_gid_entries,
				     int num_ids, ZOLTAN_ID_PTR gids, 
				     int* sizes,  int* indices,
				     char *Buf, int *ierr)
{
  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  ContactSearch*   search   = (ContactSearch *)data;
  ContactTopology* topology = search->Get_Primary_Topology(); 

  ZOLTAN_ID_PTR gid = gids;
  for( int i=0 ; i<num_ids ; ++i ){
    int entity_type = ContactZoltanGID::Type(gid);
    char* buf       = Buf + indices[i];
    int*  ibuf      = reinterpret_cast<int*>(buf);
    int block       = ibuf[ContactTopologyEntity<Real>::BLOCK_ID];
    switch (entity_type) {
    case CT_NODE:
    case CT_SHELL_NODE:
      PRECONDITION( ibuf[0] == CT_NODE || 
                    ibuf[0] == CT_SHELL_NODE);
      PRECONDITION( block>=0 && block<topology->Number_of_Node_Blocks() );
      topology->Node_Block(block)->Insert_Node(buf);
    break;
    case CT_FACE:
      PRECONDITION( ibuf[0] == CT_FACE );
      PRECONDITION( block>=0 && block<topology->Number_of_Face_Blocks() );
      topology->Face_Block(block)->Insert_Face(buf);
    break;
    case CT_ELEMENT:
      PRECONDITION( ibuf[0] == CT_ELEMENT );
      PRECONDITION( block>=0 && block<topology->Number_of_Element_Blocks() );
      topology->Element_Block(block)->Insert_Element(buf);
    break;
    default:
      POSTCONDITION(false);
      break;
    }
    gid += num_gid_entries;
  }
  *ierr = ZOLTAN_OK;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//      The following three routines are for migrating topological
//      entities for ghosting in the primary topology.
//
//-------------------------------------------------------------------------

void ContactGhostingExportSizes(void *data,  
			        int num_gid_entries, 
                                int num_lid_entries, 
			        int num_ids,
			        ZOLTAN_ID_PTR gids, 
                                ZOLTAN_ID_PTR lids,
			        int* sizes,
			        int *ierr )
{
  ContactSearch* search = (ContactSearch *)data;

#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
  //ContactParOStream& postream = search->ParOStream();
  Real bytes_nodes = search->Bytes_For_Nodes();
  Real bytes_faces = search->Bytes_For_Faces();
  Real bytes_elems = search->Bytes_For_Elems();
#endif

  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  PRECONDITION(num_lid_entries==ZOLTAN_LID_SIZE);
  ContactTopology* topology = search->Get_Primary_Topology();
  ContactNode<Real>**    nodes    = reinterpret_cast<ContactNode<Real>**>   (topology->NodeList()->EntityList());
  ContactFace<Real>**    faces    = reinterpret_cast<ContactFace<Real>**>   (topology->FaceList()->EntityList());
  ContactElement** elements = reinterpret_cast<ContactElement**>(topology->ElemList()->EntityList());

  int edgebased_physical_face  = (search->PhysicalFaceAlgorithm()) == ContactSearch::PF_EDGE_BASED;
  int enable_off_face_tracking = search->OffFaceTrackingStatus();

  int state1 = 1;
  ZOLTAN_ID_PTR lid = lids;
  for( int i=0 ; i<num_ids ; ++i ){
    int entity_type  = ContactZoltanLID::Type(lid);
    int entity_index = ContactZoltanLID::Index(lid);
    int& size        = sizes[i];
    switch (entity_type) {
    case CT_NODE:
    case CT_SHELL_NODE:
      {
	ContactNode<Real>* node = nodes[entity_index];
        POSTCONDITION(node);
#ifdef CONTACT_OLD_XFER
	size = node->Size(state1);
#else
	size = node->Size_ForSecondary(state1);
#endif
#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
        bytes_nodes += size;
#endif
      }
      break;
    case CT_EDGE:
      POSTCONDITION(false);
      break;
    case CT_FACE:
      {
	ContactFace<Real>* face = faces[entity_index];
        POSTCONDITION(face);
#ifdef CONTACT_OLD_XFER
	size = face->Size(enable_off_face_tracking);
#else
	size = face->Size_ForSecondary(enable_off_face_tracking,
                                       edgebased_physical_face);
#endif
#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
        bytes_faces += size;
#endif
      }
      break;
    case CT_ELEMENT:
      {
	ContactElement* element = elements[entity_index];
        POSTCONDITION(element);
#ifdef CONTACT_OLD_XFER
	size = element->Size();
#else
	size = element->Size_ForSecondary();
#endif
#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
        bytes_elems += size;
#endif
      }
      break;
    default:
      POSTCONDITION(false);
      break;
    }
    lid += num_lid_entries;
  }
  *ierr = ZOLTAN_OK;
#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
  search->Bytes_For_Nodes(bytes_nodes);
  search->Bytes_For_Faces(bytes_faces);
  search->Bytes_For_Elems(bytes_elems);
#endif
}

void ContactGhostingExportPack(void *data, 
			       int num_gid_entries, 
                               int num_lid_entries,
			       int num_ids,
			       ZOLTAN_ID_PTR gids, 
                               ZOLTAN_ID_PTR lids,
			       int* procs, int* sizes, 
                               int* indices, char *Buf, 
                               int *ierr)
{
  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  PRECONDITION(num_lid_entries==ZOLTAN_LID_SIZE);
  ContactSearch*     search   = (ContactSearch *)data;
  ContactTopology*   topology = search->Get_Primary_Topology();
  ContactNode<Real>** nodes =       reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());
  ContactFace<Real>** faces =       reinterpret_cast<ContactFace<Real>**>(topology->FaceList()->EntityList());
  ContactElement** elements = reinterpret_cast<ContactElement**>(topology->ElemList()->EntityList());

  int edgebased_physical_face  = (search->PhysicalFaceAlgorithm()) == ContactSearch::PF_EDGE_BASED;
  int enable_off_face_tracking = search->OffFaceTrackingStatus();
  
  ZOLTAN_ID_PTR lid = lids;
  for( int i=0 ; i<num_ids ; ++i ){
    int entity_type  = ContactZoltanLID::Type(lid);
    int entity_index = ContactZoltanLID::Index(lid);
    char* buf         = Buf+indices[i];
    switch (entity_type) {
    case CT_NODE:
    case CT_SHELL_NODE:
      {
        ContactNode<Real>* node = nodes[entity_index];
	POSTCONDITION(node);
#ifdef CONTACT_OLD_XFER
        node->Pack(buf, -2);
#else
        node->Pack_ForSecondary(buf, -2);
#endif
      }
      break;
    case CT_EDGE:
      POSTCONDITION(false);
      break;
    case CT_FACE:
      {
        ContactFace<Real>* face = faces[entity_index];
	POSTCONDITION(face);
#ifdef CONTACT_OLD_XFER
	face->Pack(buf,enable_off_face_tracking);
#else
	face->Pack_ForSecondary(buf,enable_off_face_tracking,
                                edgebased_physical_face);
#endif
      }
      break;
    case CT_ELEMENT:
      {
        ContactElement* element = elements[entity_index];
	POSTCONDITION(element);
#ifdef CONTACT_OLD_XFER
	element->Pack(buf);
#else
	element->Pack_ForSecondary(buf);
#endif
      }
      break;
    default:
      POSTCONDITION(false);
      break;
    }
    lid += num_lid_entries;
  }
  *ierr = ZOLTAN_OK;
}

void ContactGhostingImportSizes(void *data,  
			        int num_gid_entries, 
                                int num_lid_entries, 
			        int num_ids,
			        ZOLTAN_ID_PTR gids, 
                                ZOLTAN_ID_PTR lids,
			        int* sizes,
			        int *ierr )
{
  ContactSearch*   search   = (ContactSearch *)data;

#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
  ContactParOStream& postream = search->ParOStream();
  Real bytes_nodes = 0.0;
  Real bytes_faces = 0.0;
  Real bytes_elems = 0.0;
#endif

  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  PRECONDITION(num_lid_entries==ZOLTAN_LID_SIZE);
  ContactTopology* topology = search->Get_Primary_Topology();

  int edgebased_physical_face  = (search->PhysicalFaceAlgorithm()) == ContactSearch::PF_EDGE_BASED;
  int enable_off_face_tracking = search->OffFaceTrackingStatus();

  int state1 = 1;
  ZOLTAN_ID_PTR gid = gids;
  for( int i=0 ; i<num_ids ; ++i ){
    int& size       = sizes[i];
    int entity_type = ContactZoltanGID::Type(gid);
    ContactHostGlobalID GID(gid);
    switch (entity_type) {
    case CT_NODE:
      {
        ContactNode<Real>* node = static_cast<ContactNode<Real>*>
                            (topology->NodeList()->Find(GID));
        POSTCONDITION(node);
#ifdef CONTACT_OLD_XFER
	size = node->Size(state1);
#else
	size = node->Size_ForSecondary(state1);
#endif
#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
        bytes_nodes += size;
#endif
	break;
      }
    case CT_SHELL_NODE:
      {
        ContactShellNode* node = static_cast<ContactShellNode*>
                            (topology->NodeList()->Find(GID));
        POSTCONDITION(node);
#ifdef CONTACT_OLD_XFER
	size = node->Size(state1);
#else
	size = node->Size_ForSecondary(state1);
#endif
#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
        bytes_nodes += size;
#endif
	break;
      }
    case CT_EDGE:
      POSTCONDITION(false);
      break;
    case CT_FACE:
      {
        ContactFace<Real>* face = static_cast<ContactFace<Real>*>
                            (topology->FaceList()->Find(GID));
        POSTCONDITION(face);
#ifdef CONTACT_OLD_XFER
	size = face->Size(enable_off_face_tracking);
#else
	size = face->Size_ForSecondary(enable_off_face_tracking,
                                       edgebased_physical_face);
#endif
#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
        bytes_faces += size;
#endif
	break;
      }
    case CT_ELEMENT:
      {
        ContactElement* element = static_cast<ContactElement*>
                            (topology->ElemList()->Find(GID));
        POSTCONDITION(element);
#ifdef CONTACT_OLD_XFER
	size = element->Size();
#else
	size = element->Size_ForSecondary();
#endif
#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
        bytes_elems += size;
#endif
	break;
      }
    default:
      POSTCONDITION(false);
      break;
    }
    gid += num_gid_entries;
  }
  *ierr = ZOLTAN_OK;
#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
  bytes_nodes /= 1024;
  bytes_faces /= 1024;
  bytes_elems /= 1024;
  postream<<"    Data Sent by Zoltan when ghosting in the primary topology...\n";
  postream<<"      nodes:  "<<bytes_nodes<<" Kb\n";
  postream<<"      faces:  "<<bytes_faces<<" Kb\n";
  postream<<"      elems:  "<<bytes_elems<<" Kb\n";
#endif
}

void ContactGhostingImportPack(void *data, 
			       int num_gid_entries, 
                               int num_lid_entries,
			       int num_ids,
			       ZOLTAN_ID_PTR gids, 
                               ZOLTAN_ID_PTR lids,
			       int* procs, int* sizes, 
                               int* indices, char *Buf, 
                               int *ierr)
{
  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  PRECONDITION(num_lid_entries==ZOLTAN_LID_SIZE);
  ContactSearch*     search   = (ContactSearch *)data;
  ContactTopology*   topology = search->Get_Primary_Topology();

  int edgebased_physical_face  = (search->PhysicalFaceAlgorithm()) == ContactSearch::PF_EDGE_BASED;
  int enable_off_face_tracking = search->OffFaceTrackingStatus();
  
  ZOLTAN_ID_PTR gid = gids;
  for( int i=0 ; i<num_ids ; ++i ){
    char* buf       = Buf+indices[i];
    int entity_type = ContactZoltanGID::Type(gid);
    ContactHostGlobalID GID(gid);
    switch (entity_type) {
    case CT_NODE:
      {
        ContactNode<Real>* node = static_cast<ContactNode<Real>*>
                            (topology->NodeList()->Find(GID));
	POSTCONDITION(node);
#ifdef CONTACT_OLD_XFER
        node->Pack(buf, -2);
#else
        node->Pack_ForSecondary(buf, -2);
#endif
      }
      break;
    case CT_SHELL_NODE:
      {
        ContactShellNode* node = static_cast<ContactShellNode*>
                            (topology->NodeList()->Find(GID));
	POSTCONDITION(node);
#ifdef CONTACT_OLD_XFER
        node->Pack(buf, -2);
#else
        node->Pack_ForSecondary(buf, -2);
#endif
      }
      break;
    case CT_EDGE:
      POSTCONDITION(false);
      break;
    case CT_FACE:
      {
        ContactFace<Real>* face = static_cast<ContactFace<Real>*>
                            (topology->FaceList()->Find(GID));
	POSTCONDITION(face);
#ifdef CONTACT_OLD_XFER
	face->Pack(buf,enable_off_face_tracking);
#else
	face->Pack_ForSecondary(buf,enable_off_face_tracking,
                                edgebased_physical_face);
#endif
      }
      break;
    case CT_ELEMENT:
      {
        ContactElement* element = static_cast<ContactElement*>
                            (topology->ElemList()->Find(GID));
	POSTCONDITION(element);
#ifdef CONTACT_OLD_XFER
	element->Pack(buf);
#else
	element->Pack_ForSecondary(buf);
#endif
      }
      break;
    default:
      POSTCONDITION(false);
      break;
    }
    gid += num_gid_entries;
  }
  *ierr = ZOLTAN_OK;
}

void ContactGhostingUnpack(void *data, 
                           int num_gid_entries,
			   int num_ids, 
                           ZOLTAN_ID_PTR gids, 
			   int* sizes,  int* indices,
			   char *Buf, int *ierr)
{
  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  ContactSearch*   search   = (ContactSearch *)data;
  ContactTopology* topology = search->Get_Primary_Topology();
  

  for( int i=0 ; i<num_ids ; ++i ){
    char* buf       = Buf + indices[i];
    int*  ibuf      = reinterpret_cast<int*>(buf);
    int block       = ibuf[ContactTopologyEntity<Real>::BLOCK_ID];
    ContactHostGlobalID GID(ibuf[ContactTopologyEntity<Real>::GID_HI],
                            ibuf[ContactTopologyEntity<Real>::GID_LO]);
    switch (ibuf[0]) {
    case CT_NODE:
    case CT_SHELL_NODE:
      PRECONDITION( block>=0 && block<topology->Number_of_Node_Blocks() );
      if (topology->Node_Block(block)->NodeList()->Find(GID)==NULL)
#ifdef CONTACT_OLD_XFER
        topology->Ghosted_Node_Block(block)->Insert_Node(buf);
#else
        topology->Ghosted_Node_Block(block)->Insert_Node_ForSecondary(buf);
#endif
      break;
    case CT_EDGE:
      POSTCONDITION(false);
      break;
    case CT_FACE:
      
      PRECONDITION( block>=0 && block<topology->Number_of_Face_Blocks() );
      if (topology->Face_Block(block)->FaceList()->Find(GID)==NULL)
#ifdef CONTACT_OLD_XFER
        topology->Ghosted_Face_Block(block)->Insert_Face(buf);
#else
        topology->Ghosted_Face_Block(block)->Insert_Face_ForSecondary(buf);
#endif
      break;
    case CT_ELEMENT:
      PRECONDITION( block>=0 && block<topology->Number_of_Element_Blocks() );
      if (topology->Element_Block(block)->ElemList()->Find(GID)==NULL)
#ifdef CONTACT_OLD_XFER
        topology->Ghosted_Element_Block(block)->Insert_Element(buf);
#else
        topology->Ghosted_Element_Block(block)->Insert_Element_ForSecondary(buf);
#endif
      break;
    default:
      POSTCONDITION(false);
      break;
    }
  }
  
  *ierr = ZOLTAN_OK;
}

void ContactUpdateGhostingExportSizes(void *data,  
			        int num_gid_entries, 
                                int num_lid_entries, 
			        int num_ids,
			        ZOLTAN_ID_PTR gids, 
                                ZOLTAN_ID_PTR lids,
			        int* sizes,
			        int *ierr )
{
  ContactSearch* search = (ContactSearch *)data;

#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
  //ContactParOStream& postream = search->ParOStream();
  Real bytes_nodes = search->Bytes_For_Nodes();
  Real bytes_faces = search->Bytes_For_Faces();
  Real bytes_elems = search->Bytes_For_Elems();
#endif

  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  PRECONDITION(num_lid_entries==ZOLTAN_LID_SIZE);
  ContactTopology* topology = search->Get_Primary_Topology();
  ContactNode<Real>**    nodes    = reinterpret_cast<ContactNode<Real>**>   (topology->NodeList()->EntityList());
  ContactFace<Real>**    faces    = reinterpret_cast<ContactFace<Real>**>   (topology->FaceList()->EntityList());
  ContactElement** elements = reinterpret_cast<ContactElement**>(topology->ElemList()->EntityList());

  ZOLTAN_ID_PTR lid = lids;
  for( int i=0 ; i<num_ids ; ++i ){
    int entity_type  = ContactZoltanLID::Type(lid);
    int entity_index = ContactZoltanLID::Index(lid);
    int& size        = sizes[i];
    switch (entity_type) {
    case CT_NODE:
    case CT_SHELL_NODE:
      {
	ContactNode<Real>* node = nodes[entity_index];
        POSTCONDITION(node);
	size = node->Size_ForDataUpdate();
#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
        bytes_nodes += size;
#endif
      }
      break;
    case CT_EDGE:
      POSTCONDITION(false);
      break;
    case CT_FACE:
      {
	ContactFace<Real>* face = faces[entity_index];
        POSTCONDITION(face);
	size = face->Size_ForDataUpdate();
#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
        bytes_faces += size;
#endif
      }
      break;
    case CT_ELEMENT:
      {
	ContactElement* element = elements[entity_index];
        POSTCONDITION(element);
	size = element->Size_ForDataUpdate();
#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
        bytes_elems += size;
#endif
      }
      break;
    default:
      POSTCONDITION(false);
      break;
    }
    lid += num_lid_entries;
  }
  *ierr = ZOLTAN_OK;
#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
  search->Bytes_For_Nodes(bytes_nodes);
  search->Bytes_For_Faces(bytes_faces);
  search->Bytes_For_Elems(bytes_elems);
#endif
}

void ContactUpdateGhostingExportPack(void *data, 
			       int num_gid_entries, 
                               int num_lid_entries,
			       int num_ids,
			       ZOLTAN_ID_PTR gids, 
                               ZOLTAN_ID_PTR lids,
			       int* procs, int* sizes, 
                               int* indices, char *Buf, 
                               int *ierr)
{
  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  PRECONDITION(num_lid_entries==ZOLTAN_LID_SIZE);
  ContactSearch*   search   = (ContactSearch *)data;
  ContactTopology* topology = search->Get_Primary_Topology();
  ContactNode<Real>**    nodes    = reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());
  ContactFace<Real>**    faces    = reinterpret_cast<ContactFace<Real>**>(topology->FaceList()->EntityList());
  ContactElement** elements = reinterpret_cast<ContactElement**>(topology->ElemList()->EntityList());
  
  ZOLTAN_ID_PTR lid = lids;
  for( int i=0 ; i<num_ids ; ++i ){
    int entity_type  = ContactZoltanLID::Type(lid);
    int entity_index = ContactZoltanLID::Index(lid);
    char* buf         = Buf+indices[i];
    switch (entity_type) {
    case CT_NODE:
    case CT_SHELL_NODE:
      {
        ContactNode<Real>* node = nodes[entity_index];
	POSTCONDITION(node);
        node->Pack_ForDataUpdate(buf);
      }
      break;
    case CT_EDGE:
      POSTCONDITION(false);
      break;
    case CT_FACE:
      {
        ContactFace<Real>* face = faces[entity_index];
	POSTCONDITION(face);
	face->Pack_ForDataUpdate(buf);
      }
      break;
    case CT_ELEMENT:
      {
        ContactElement* element = elements[entity_index];
	POSTCONDITION(element);
	element->Pack_ForDataUpdate(buf);
      }
      break;
    default:
      POSTCONDITION(false);
      break;
    }
    lid += num_lid_entries;
  }
  *ierr = ZOLTAN_OK;
}

void ContactUpdateGhostingImportSizes(void *data,  
			        int num_gid_entries, 
                                int num_lid_entries, 
			        int num_ids,
			        ZOLTAN_ID_PTR gids, 
                                ZOLTAN_ID_PTR lids,
			        int* sizes,
			        int *ierr )
{
  ContactSearch* search = (ContactSearch *)data;

#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
  //ContactParOStream& postream = search->ParOStream();
  Real bytes_nodes = search->Bytes_For_Nodes();
  Real bytes_faces = search->Bytes_For_Faces();
  Real bytes_elems = search->Bytes_For_Elems();
#endif

  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  PRECONDITION(num_lid_entries==ZOLTAN_LID_SIZE);
  ContactTopology* topology = search->Get_Primary_Topology();

  ZOLTAN_ID_PTR gid = gids;
  for( int i=0 ; i<num_ids ; ++i ){
    int entity_type  = ContactZoltanGID::Type(gid);
    int& size        = sizes[i];
    ContactHostGlobalID GID(gid);
    switch (entity_type) {
    case CT_NODE:
    case CT_SHELL_NODE:
      {
        ContactNode<Real>* node = static_cast<ContactNode<Real>*>
                            (topology->NodeList()->Find(GID));
        POSTCONDITION(node);
	size = node->Size_ForDataUpdate();
#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
        bytes_nodes += size;
#endif
      }
      break;
    case CT_EDGE:
      POSTCONDITION(false);
      break;
    case CT_FACE:
      {
        ContactFace<Real>* face = static_cast<ContactFace<Real>*>
                            (topology->FaceList()->Find(GID));
        POSTCONDITION(face);
	size = face->Size_ForDataUpdate();
#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
        bytes_faces += size;
#endif
      }
      break;
    case CT_ELEMENT:
      {
        ContactElement* element = static_cast<ContactElement*>
                            (topology->ElemList()->Find(GID));
        POSTCONDITION(element);
	size = element->Size_ForDataUpdate();
#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
        bytes_elems += size;
#endif
      }
      break;
    default:
      POSTCONDITION(false);
      break;
    }
    gid += num_gid_entries;
  }
  *ierr = ZOLTAN_OK;
#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
  search->Bytes_For_Nodes(bytes_nodes);
  search->Bytes_For_Faces(bytes_faces);
  search->Bytes_For_Elems(bytes_elems);
#endif
}

void ContactUpdateGhostingImportPack(void *data, 
			       int num_gid_entries, 
                               int num_lid_entries,
			       int num_ids,
			       ZOLTAN_ID_PTR gids, 
                               ZOLTAN_ID_PTR lids,
			       int* procs, int* sizes, 
                               int* indices, char *Buf, 
                               int *ierr)
{
  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  PRECONDITION(num_lid_entries==ZOLTAN_LID_SIZE);
  ContactSearch*   search   = (ContactSearch *)data;
  ContactTopology* topology = search->Get_Primary_Topology();
  
  ZOLTAN_ID_PTR gid = gids;
  for( int i=0 ; i<num_ids ; ++i ){
    int entity_type = ContactZoltanGID::Type(gid);
    char* buf       = Buf+indices[i];
    ContactHostGlobalID GID(gid);
    switch (entity_type) {
    case CT_NODE:
    case CT_SHELL_NODE:
      {
        ContactNode<Real>* node = static_cast<ContactNode<Real>*>
                            (topology->NodeList()->Find(GID));
	POSTCONDITION(node);
        node->Pack_ForDataUpdate(buf);
      }
      break;
    case CT_EDGE:
      POSTCONDITION(false);
      break;
    case CT_FACE:
      {
        ContactFace<Real>* face = static_cast<ContactFace<Real>*>
                            (topology->FaceList()->Find(GID));
	POSTCONDITION(face);
	face->Pack_ForDataUpdate(buf);
      }
      break;
    case CT_ELEMENT:
      {
        ContactElement* element = static_cast<ContactElement*>
                            (topology->ElemList()->Find(GID));
	POSTCONDITION(element);
	element->Pack_ForDataUpdate(buf);
      }
      break;
    default:
      POSTCONDITION(false);
      break;
    }
    gid += num_gid_entries;
  }
  *ierr = ZOLTAN_OK;
}

void ContactUpdateGhostingUnpack(void *data, 
                           int num_gid_entries,
			   int num_ids, 
                           ZOLTAN_ID_PTR gids, 
			   int* sizes,  int* indices,
			   char *Buf, int *ierr)
{
  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  ContactSearch*   search   = (ContactSearch *)data;
  ContactTopology* topology = search->Get_Primary_Topology();

  ZOLTAN_ID_PTR gid = gids;
  for( int i=0 ; i<num_ids ; ++i ){
    char* buf       = Buf + indices[i];
    int*  ibuf      = reinterpret_cast<int*>(buf);
    int block       = ibuf[0];
    int entity_type = ContactZoltanGID::Type(gid);
    ContactHostGlobalID GID(gid);
    switch (entity_type) {
    case CT_NODE:
    case CT_SHELL_NODE:
      PRECONDITION( block>=0 && block<topology->Number_of_Node_Blocks() );
      {
        ContactNode<Real>* node = reinterpret_cast<ContactNode<Real>*>(topology->Ghosted_Node_Block(block)->NodeList()->Find(GID));
        POSTCONDITION(node);
        node->Unpack_ForDataUpdate(buf);
      }
      break;
    case CT_EDGE:
      POSTCONDITION(false);
      break;
    case CT_FACE:
      PRECONDITION( block>=0 && block<topology->Number_of_Face_Blocks() );
      {
        ContactFace<Real>* face = reinterpret_cast<ContactFace<Real>*>(topology->Ghosted_Face_Block(block)->FaceList()->Find(GID));
        POSTCONDITION(face);
        face->Unpack_ForDataUpdate(buf);
      }
      break;
    case CT_ELEMENT:
      PRECONDITION( block>=0 && block<topology->Number_of_Element_Blocks() );
      {
        ContactElement* element = reinterpret_cast<ContactElement*>(topology->Ghosted_Element_Block(block)->ElemList()->Find(GID));
        POSTCONDITION(element);
        element->Unpack_ForDataUpdate(buf);
      }
      break;
    default:
      POSTCONDITION(false);
      break;
    }
    gid += num_gid_entries;
  }
  *ierr = ZOLTAN_OK;
}

void ContactUpdateTiedImportSizes(void *data,  
                                    int num_gid_entries, 
                                    int num_lid_entries, 
                                    int num_ids,
                                    ZOLTAN_ID_PTR gids, 
                                    ZOLTAN_ID_PTR lids,
                                    int* sizes,
                                    int *ierr)
{
  ContactSearch* search = (ContactSearch *)data;

  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  PRECONDITION(num_lid_entries==ZOLTAN_LID_SIZE);
  ContactTopology* topology = search->Get_Primary_Topology();

  ZOLTAN_ID_PTR gid = gids;
  for( int i=0 ; i<num_ids ; ++i ){
    int entity_type  = ContactZoltanGID::Type(gid);
    ContactHostGlobalID GID(gid);
    switch (entity_type) {
    case CT_FACE:
      {
        ContactFace<Real>* face = static_cast<ContactFace<Real>*>
                            (topology->FaceList()->Find(GID));
        POSTCONDITION(face);
        face->SetContextBit(ContactTopologyEntity<Real>::TIED);
        // face_type, face normal & nodal coordinates
        sizes[i] = sizeof(Real)*(1+3*(1+face->Nodes_Per_Face()));
      }
      break;
    default:
      POSTCONDITION(false);
      break;
    }
    gid += num_gid_entries;
  }
  *ierr = ZOLTAN_OK;
}

void ContactUpdateTiedImportPack(void *data, 
                                   int num_gid_entries, 
                                   int num_lid_entries,
                                   int num_ids,
                                   ZOLTAN_ID_PTR gids, 
                                   ZOLTAN_ID_PTR lids,
                                   int* procs, int* sizes, 
                                   int* indices, char *Buf, 
                                   int *ierr)
{
  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  PRECONDITION(num_lid_entries==ZOLTAN_LID_SIZE);
  ContactSearch*   search   = (ContactSearch *)data;
  ContactTopology* topology = search->Get_Primary_Topology();
  VariableHandle FACE_NORMAL =
    topology->Variable_Handle( ContactTopology::Face_Normal );
  VariableHandle POSITION = topology->Variable_Handle( ContactTopology::Current_Position );
  
  if (search->NumConfigs()>1) {
    POSITION = topology->Variable_Handle( ContactTopology::Predicted_Position );
  }
    
  ZOLTAN_ID_PTR gid = gids;
  for( int i=0 ; i<num_ids ; ++i ){
    int entity_type = ContactZoltanGID::Type(gid);
    char* buf       = Buf+indices[i];
    ContactHostGlobalID GID(gid);
    switch (entity_type) {
    case CT_FACE:
      {
        ContactFace<Real>* face = static_cast<ContactFace<Real>*>
                            (topology->FaceList()->Find(GID));
        POSTCONDITION(face);
        int j = 0;
        Real* r_buf = reinterpret_cast<Real*>( buf );
        Real* face_normal = face->Variable(FACE_NORMAL);
        r_buf[j++] = face->FaceType();
        r_buf[j++] = face_normal[0];
        r_buf[j++] = face_normal[1];
        r_buf[j++] = face_normal[2];
        for(int k=0 ; k<face->Nodes_Per_Face(); ++k ){
          Real* coords = face->Node(k)->Variable(POSITION);
          r_buf[j++] = coords[0];
          r_buf[j++] = coords[1];
          r_buf[j++] = coords[2];
        }
      }
      break;
    default:
      POSTCONDITION(false);
      break;
    }
    gid += num_gid_entries;
  }
  *ierr = ZOLTAN_OK;
}

void ContactUpdateTiedUnpack(void *data, 
                             int num_gid_entries,
			     int num_ids, 
                             ZOLTAN_ID_PTR gids, 
			     int* sizes,  int* indices,
			     char *Buf, int *ierr)
{
  PRECONDITION(num_gid_entries==ZOLTAN_GID_SIZE);
  ContactSearch*   search   = (ContactSearch *)data;
  ContactTopology* topology = search->Get_Primary_Topology();
  ContactFixedSizeAllocator* allocators = search->Get_Allocators();
  VariableHandle POSITION = topology->Variable_Handle( ContactTopology::Current_Position );
  
  if (search->NumConfigs()>1) {
    POSITION = topology->Variable_Handle( ContactTopology::Predicted_Position );
  }

  //need to create a map to return the appropriate 
  //location in Buf for a given face gid
  std::map<ContactHostGlobalID, int> lookup;
  for( int i=0 ; i<num_ids ; ++i ){
    ContactHostGlobalID GID(&gids[i*num_gid_entries]);
    lookup[GID] = i;
  }
  
  MPI_Comm SearchComm = search->Get_Comm();
  int my_proc         = contact_processor_number( SearchComm );
  int number_of_nodes = topology->Number_of_Nodes();
  ContactNode<Real>** Nodes = reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());
  for (int i=0; i<number_of_nodes; ++i) {
    ContactNode<Real>* node = Nodes[i];
    if(node->Ownership() == ContactTopologyEntity<Real>::OWNED){
      int num_tied_at_node = 0;
      ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions(ContactSearch::STATE1);
      for (int j=0; j<node->Number_NodeEntity_Interactions(ContactSearch::STATE1); ++j) {
        ContactNodeEntityInteraction* cnei = interactions[j];
        if( cnei->Is_Tied() || cnei->Is_InfSlip() ){
          if (cnei->Get_Type()==ContactNodeEntityInteraction::NODE_FACE_INTERACTION) {
            ContactNodeFaceInteraction *cnfi = static_cast<ContactNodeFaceInteraction*>(cnei);
              if (cnfi->FaceEntityData()->owner != my_proc) {
              ContactHostGlobalID Face_GID( cnfi->FaceEntityData()->host_gid[0], 
                                            cnfi->FaceEntityData()->host_gid[1] );
              PRECONDITION(lookup.find(Face_GID) != lookup.end());
              int index    = lookup.find(Face_GID)->second;
              Real* buffer = (Real*)(Buf+indices[index]);
              Real  gap0   = cnei->Get_Gap_Initial();
              Real* norm   = cnei->Get_Physical_Face_Normal();
              ContactNodeEntityInteraction *nei = cnei->New_Instance(allocators);
              nei->Set_Source(ContactNodeEntityInteraction::RETRIEVED_TIED);
              nei->Update_Tied_Interaction( POSITION, buffer, gap0, norm );
              node->Store_NodeEntity_Interaction( num_tied_at_node++, nei );
            }
          }
        }//end is tied if
      }//end loop over all interactions on node
    }//end ownership if
  }//end loop over all nodes
  *ierr = ZOLTAN_OK;
}

#endif



