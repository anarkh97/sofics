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

#include "ContactCommList.h"
#include "ContactTopology.h"
#include "ContactVector.h"
#include "ContactParOStream.h"

using namespace std;

//*******************************************************************************************************
//
//  Default constructor, just initialize everything to NULL or zero
//
//*******************************************************************************************************

ContactCommList::ContactCommList() :
  num_comm_partners(0),
  num_entities(0),
  comm_proc_ids(NULL),
  num_to_proc(NULL),
  offset(NULL),
  entity_list(NULL),
  entity_index_list(NULL)
{}

//*******************************************************************************************************
//
//  Allocate list using stored sizes
//
//*******************************************************************************************************
void ContactCommList::Allocate() {
  if(num_comm_partners) {
    comm_proc_ids = new int[num_comm_partners];
    num_to_proc   = new int[num_comm_partners];
    offset        = new int[num_comm_partners];
    entity_list   = new ContactTopologyEntity<Real>*[num_entities];
  }
}

//*******************************************************************************************************
//
//  Print the current communication list to the terminal
//
//*******************************************************************************************************
void ContactCommList::Print(const char* line_prefix) {
  cout<<line_prefix<<"num_com_partners: "<<num_comm_partners<<endl;
  cout<<line_prefix<<"num_entities: "<<num_entities<<endl;
  for(int i = 0; i < num_comm_partners; ++i) {
    cout<<line_prefix<<"Processor List "<<i<<endl;
    cout<<line_prefix<<"  comm_proc_id: "<<comm_proc_ids[i]<<endl;
    cout<<line_prefix<<"  num_to_proc: "<<num_to_proc[i]<<endl;
    cout<<line_prefix<<"  offset: "<<offset[i]<<endl;
    cout<<line_prefix<<"  Entities"<<endl;
    for(int j = 0; j < num_to_proc[i]; ++j) {
      int index = entity_index_list[offset[i] + j];
      ContactTopologyEntity<Real> *entity = entity_list[offset[i] + j];
      cout<<line_prefix<<"    entity_id: "<<entity->Global_ID()<<", index = "<<index<<endl;
    }
  }
}

void ContactCommList::Print(const char* line_prefix, ContactParOStream& os) {
  os<<line_prefix<<"num_com_partners: "<<num_comm_partners<<"\n";
  os<<line_prefix<<"num_entities: "<<num_entities<<"\n";
  for(int i = 0; i < num_comm_partners; ++i) {
    os<<line_prefix<<"Processor List "<<i<<"\n";
    os<<line_prefix<<"  comm_proc_id: "<<comm_proc_ids[i]<<"\n";
    os<<line_prefix<<"  num_to_proc: "<<num_to_proc[i]<<"\n";
    os<<line_prefix<<"  offset: "<<offset[i]<<"\n";
    os<<line_prefix<<"  Entities"<<"\n";
    for(int j = 0; j < num_to_proc[i]; ++j) {
      int index = entity_index_list[offset[i] + j];
      ContactTopologyEntity<Real> *entity = entity_list[offset[i] + j];
      os<<line_prefix<<"    entity_id: "<<entity->Global_ID()<<", index = "<<index<<"\n";
    }
  }
}

//*******************************************************************************************************
//
//  Calculate the offsets for the lists
//
//*******************************************************************************************************
void ContactCommList::CalculateOffset() {
  if( num_comm_partners ){
    offset[0] = 0;
    for( int i=1 ; i<num_comm_partners ; ++i )
      offset[i] = offset[i-1] + num_to_proc[i-1];
  }
}



//*******************************************************************************************************
//
//  Sort all the entries in a communication list by entity global id.  The sorting is nessecary for
//  merge operations and also should yeild better cache coherency during use.
//
//*******************************************************************************************************
void ContactCommList::Sort() {
  for(int ilist = 0 ; ilist < num_comm_partners; ++ilist) {
    ContactTopology::SortEntityList(num_to_proc[ilist], entity_list + offset[ilist]);
  } 
}


//*******************************************************************************************************
//
//  Initialize the current list as the result of the addition of two input lists
//
//*******************************************************************************************************
void ContactCommList::CreateFromListAddition(const ContactCommList &com1, 
                                        const ContactCommList &com2) {
  //
  //  Create temporay data structures for use in the communication merge
  //  NKC, really should replace all this with STL vectors, darn lack of 
  //  templates!!!!!!!
  //
  int num_procs;
  MPI_Comm_size( MPI_COMM_WORLD, &num_procs );

  ACME::Int_Vector com1_list_size(num_procs,0);
  ACME::Int_Vector com2_list_size(num_procs,0);
  ACME::Int_Vector com1_offset(num_procs,-1);
  ACME::Int_Vector com2_offset(num_procs,-1);
  ACME::Int_Vector new_com_list_size(num_procs);

  //
  //  Determine the maximum communiation data array sizes by adding the per processor sizes
  //  of the two input communication objects
  //
  for(int iproc = 0; iproc < com1.num_comm_partners; ++iproc) {
    int proc_num = com1.comm_proc_ids[iproc];
    com1_list_size[proc_num] = com1.num_to_proc[iproc];
    com1_offset[proc_num] = com1.offset[iproc];
  }
  for(int iproc = 0; iproc < com2.num_comm_partners; ++iproc) {
    int proc_num = com2.comm_proc_ids[iproc];
    com2_list_size[proc_num] = com2.num_to_proc[iproc];
    com2_offset[proc_num] = com2.offset[iproc];
  }
  //
  //  Allocate the temporary communication object arrays
  //  entity_list_temp is a matrix of ContactTopologyEntity<Real> pointers.
  //  The matrix has num_procs rows and num_to_send/receive length for each row.
  //
  ContactTopologyEntity<Real>*** entity_list_temp = new ContactTopologyEntity<Real>**[num_procs];
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    new_com_list_size[iproc] = com1_list_size[iproc] + com2_list_size[iproc];
    entity_list_temp[iproc] = new ContactTopologyEntity<Real>*[new_com_list_size[iproc]];
  }
  //
  //  Loop over the lists from the input comm specs to created a combined sorted
  //  list.  Delete any duplicate entries and keep track of the total number of
  //  objects in the lists.
  //
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    int com1_index = 0;
    int com2_index = 0;
    int new_com_index = 0;
    while(true) {
      ContactTopologyEntity<Real>* com1_entity = NULL;
      ContactTopologyEntity<Real>* com2_entity = NULL;
      if(com1_index < com1_list_size[iproc]) com1_entity = com1.entity_list[com1_offset[iproc] + com1_index];
      if(com2_index < com2_list_size[iproc]) com2_entity = com2.entity_list[com2_offset[iproc] + com2_index];
      if(com1_entity && com2_entity) {
        if(com1_entity->Global_ID() < com2_entity->Global_ID()) {
          entity_list_temp[iproc][new_com_index++] = com1_entity;
          com1_index++;
        } else if (com1_entity->Global_ID() > com2_entity->Global_ID()) {
          entity_list_temp[iproc][new_com_index++] = com2_entity;
          com2_index++;
        } else {
          entity_list_temp[iproc][new_com_index++] = com1_entity;
          com1_index++;
          com2_index++; 
          new_com_list_size[iproc]--;
        }
      } else if(com1_entity) {
        entity_list_temp[iproc][new_com_index++] = com1_entity;
        com1_index++;
      } else if(com2_entity) {
        entity_list_temp[iproc][new_com_index++] = com2_entity;
        com2_index++;
      } else {
        break;
      }
    }
  }
  //
  //  Create the actual AsymComm data, convert the processor lists into the correct format and allocate everything to
  //  the exact required size.
  //
  //  Step 1: Detemine the actual number of communication partners and total number of entities
  //
  num_comm_partners = 0;
  num_entities      = 0;
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(new_com_list_size[iproc] > 0) {
      ++num_comm_partners;
      num_entities += new_com_list_size[iproc];
    } 
  }
  //
  //  Allocate all communication specification data
  //
  comm_proc_ids = new int[num_comm_partners];
  num_to_proc   = new int[num_comm_partners];
  offset        = new int[num_comm_partners];
  entity_list   = new ContactTopologyEntity<Real>*[num_entities];
  //
  //  Fill the communication data
  //
  int comm_partner = 0;
  int total = 0;

  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(new_com_list_size[iproc] > 0) {
      comm_proc_ids[comm_partner] = iproc;
      num_to_proc[comm_partner] = new_com_list_size[iproc];
      offset[comm_partner] = total;
      for(int ientity = 0; ientity < num_to_proc[comm_partner]; ++ientity) {
        entity_list[ientity+total] = entity_list_temp[iproc][ientity];
      }
      total += num_to_proc[comm_partner];
      ++comm_partner;
    }
  }
  //
  //  Delete all the temporary arrays
  //
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    delete [] entity_list_temp[iproc];
  }
  delete [] entity_list_temp;
}


//*******************************************************************************************************
//
//  Initialize the current list from passed in Zoltan data
//
//*******************************************************************************************************
void ContactCommList::CreateFromIDs(const int num,
                               const LB_ID_PTR gids,
                               const LB_ID_PTR lids,
			       const int* procs,
                               ContactTopologyEntityHash& entity_hash) {

  // Allocate some scratch
  int*  iscratch = new int[num];
  //
  // Build the data information
  //
  num_entities = num;
  num_comm_partners = 0;
  int comm_proc;
  // Build a unique list of processors to communicate with
  for(int i=0 ; i<num_entities ; ++i ){
    comm_proc = procs[i];
    bool Proc_Exists = false;
    for(int j=0 ; j<num_comm_partners ; ++j ){
      if( iscratch[j] == comm_proc ) Proc_Exists = true;
    }
    if( !Proc_Exists ) iscratch[num_comm_partners++] = comm_proc;
  }
  // Now build the actual import data arrays
  if( num_comm_partners ){
    comm_proc_ids = new int[num_comm_partners];
    num_to_proc = new int[num_comm_partners];
    offset = new int[num_comm_partners];
    entity_list = new ContactTopologyEntity<Real>*[num_entities];
    std::memcpy( comm_proc_ids, iscratch, 
	    num_comm_partners*sizeof(int) );
    int entity_index = 0;
    for(int i=0 ; i<num_comm_partners ; ++i ){
      comm_proc = comm_proc_ids[i];
      offset[i] = entity_index;
      num_to_proc[i] = 0;
      for(int j=0 ; j<num_entities ; ++j ){
	if( procs[j] == comm_proc ){
	  ContactHostGlobalID GID( &gids[j*ZOLTAN_GID_SIZE] );
	  entity_list[entity_index++] = (ContactTopologyEntity<Real>*) entity_hash.find( GID );
	  ++num_to_proc[i];
	}
      }
    }
    POSTCONDITION( entity_index == num_entities );
  } else {
    comm_proc_ids = NULL;
    num_to_proc = NULL;
    offset = NULL;
    entity_list = NULL;
  }
    
  if( iscratch ) delete [] iscratch;
}


//*******************************************************************************************************
//
//  Initialize the current list from passed in Zoltan data
//
//*******************************************************************************************************
void ContactCommList::CreateFromIDs(const int num,
                               const LB_ID_PTR gids,
                               const LB_ID_PTR lids,
			       const int* procs,
                               ContactTopologyEntityHash& entity_hash,
                               ContactTopologyEntityList& lookup_list) {

  // Allocate some scratch
  int*  iscratch = new int[num];
  //
  // Build the data information
  //
  num_entities = num;
  num_comm_partners = 0;
  int comm_proc;
  // Build a unique list of processors to communicate with
  for(int i=0 ; i<num_entities ; ++i ){
    comm_proc = procs[i];
    bool Proc_Exists = false;
    for(int j=0 ; j<num_comm_partners ; ++j ){
      if( iscratch[j] == comm_proc ) Proc_Exists = true;
    }
    if( !Proc_Exists ) iscratch[num_comm_partners++] = comm_proc;
  }
  // Now build the actual import data arrays
  if( num_comm_partners ){
    comm_proc_ids = new int[num_comm_partners];
    num_to_proc = new int[num_comm_partners];
    offset = new int[num_comm_partners];
    entity_list = new ContactTopologyEntity<Real>*[num_entities];
    std::memcpy( comm_proc_ids, iscratch, 
	    num_comm_partners*sizeof(int) );
    int entity_index = 0;
    for(int i=0 ; i<num_comm_partners ; ++i ){
      comm_proc = comm_proc_ids[i];
      offset[i] = entity_index;
      num_to_proc[i] = 0;
      for(int j=0 ; j<num_entities ; ++j ){
	if( procs[j] == comm_proc ){
	  ContactHostGlobalID GID( &gids[j*ZOLTAN_GID_SIZE] );
          ContactTopologyEntity<Real> *entity = lookup_list.Find(GID);
          if(entity == NULL) entity = (ContactTopologyEntity<Real>*)entity_hash.find(GID );
	  entity_list[entity_index++] = entity;
	  ++num_to_proc[i];
	}
      }
    }
    POSTCONDITION( entity_index == num_entities );
  } else {
    comm_proc_ids = NULL;
    num_to_proc = NULL;
    offset = NULL;
    entity_list = NULL;
  }
    
  if( iscratch ) delete [] iscratch;
}


//*******************************************************************************************************
//
//  Initialize the current list from passed in Zoltan data
//
//*******************************************************************************************************
void ContactCommList::CreateFromIDs(const int num,
                               const LB_ID_PTR gids,
                               const LB_ID_PTR lids,
			       const int* procs,
                               ContactTopologyEntityList* lookup_list, 
                               ContactType base_type) {
  num_entities = 0;
  for(int i=0 ; i<num ; ++i ){
    int entity_type = ContactZoltanLID::Type(&lids[i*ZOLTAN_LID_SIZE]);
    if (entity_type==base_type) ++num_entities;
    if (entity_type==CT_NODE && base_type==CT_SHELL_NODE) ++num_entities;
  }
  
  
  // Allocate some scratch
  int*  iscratch = new int[num_entities];
  //
  // Build the data information
  //            
  num_comm_partners = 0;
  int comm_proc;
  // Build a unique list of processors to communicate with
  for(int i=0 ; i<num ; ++i ){
    int entity_type = ContactZoltanGID::Type(&gids[i*ZOLTAN_GID_SIZE]);
    if (entity_type==base_type ||
        (entity_type==CT_NODE && base_type==CT_SHELL_NODE) ) {
      comm_proc = procs[i];
      bool Proc_Exists = false;
      for(int j=0 ; j<num_comm_partners ; ++j ){
        if( iscratch[j] == comm_proc ) Proc_Exists = true;
      }
      if( !Proc_Exists ) iscratch[num_comm_partners++] = comm_proc;
    }
  }
  
  
  // Now build the actual import data arrays
  if( num_comm_partners ){
    comm_proc_ids = new int[num_comm_partners];
    num_to_proc = new int[num_comm_partners];
    offset = new int[num_comm_partners];
    entity_list = new ContactTopologyEntity<Real>*[num_entities];
    std::memcpy( comm_proc_ids, iscratch, 
	         num_comm_partners*sizeof(int) );
    int entity_index = 0;
    for(int i=0 ; i<num_comm_partners ; ++i ){
      comm_proc = comm_proc_ids[i];
      offset[i] = entity_index;
      num_to_proc[i] = 0;
      for(int j=0 ; j<num ; ++j ){
        int entity_type = ContactZoltanGID::Type(&gids[j*ZOLTAN_GID_SIZE]);
        if (entity_type==base_type ||
           (entity_type==CT_NODE && base_type==CT_SHELL_NODE) ) {
	  if( procs[j] == comm_proc ){
	    ContactHostGlobalID GID( &gids[j*ZOLTAN_GID_SIZE] );
            ContactTopologyEntity<Real>* entity = lookup_list->Find(GID);
            POSTCONDITION(entity);
	    entity_list[entity_index++] = entity;
	    ++num_to_proc[i];
	  }
	}
      }
    }
    POSTCONDITION( entity_index == num_entities );
  } else {
    comm_proc_ids = NULL;
    num_to_proc = NULL;
    offset = NULL;
    entity_list = NULL;
  }
    
  if( iscratch ) delete [] iscratch;
}

void ContactCommList::CreateFromIDs(const int num,
                                    const LB_ID_PTR gids,
                                    const LB_ID_PTR lids,
			            const int* procs,
                                    ContactTopologyEntityList& lookup_list) {
  num_entities = num;
  
  // Allocate some scratch
  int*  iscratch = new int[num_entities];
  //
  // Build the data information
  //            
  num_comm_partners = 0;
  int comm_proc;
  // Build a unique list of processors to communicate with
  for(int i=0 ; i<num ; ++i ){
    comm_proc = procs[i];
    bool Proc_Exists = false;
    for(int j=0 ; j<num_comm_partners ; ++j ){
      if( iscratch[j] == comm_proc ) Proc_Exists = true;
    }
    if( !Proc_Exists ) iscratch[num_comm_partners++] = comm_proc;
  }
  
  // Now build the actual import data arrays
  if( num_comm_partners ){
    comm_proc_ids = new int[num_comm_partners];
    num_to_proc   = new int[num_comm_partners];
    offset        = new int[num_comm_partners];
    entity_list   = new ContactTopologyEntity<Real>*[num_entities];
    std::memcpy( comm_proc_ids, iscratch, 
	         num_comm_partners*sizeof(int) );
    int entity_index = 0;
    for(int i=0 ; i<num_comm_partners ; ++i ){
      comm_proc      = comm_proc_ids[i];
      offset[i]      = entity_index;
      num_to_proc[i] = 0;
      for(int j=0 ; j<num_entities ; ++j ){
	if( procs[j] == comm_proc ){
	  ContactHostGlobalID GID( &gids[j*ZOLTAN_GID_SIZE] );
          ContactTopologyEntity<Real> *entity = lookup_list.Find(GID);
	  entity_list[entity_index++] = entity;
	  ++num_to_proc[i];
          POSTCONDITION(entity != NULL);
	}
      }
    }
    POSTCONDITION( entity_index == num_entities );
  } else {
    comm_proc_ids = NULL;
    num_to_proc = NULL;
    offset = NULL;
    entity_list = NULL;
  }
    
  if( iscratch ) delete [] iscratch;
}


//*******************************************************************************************************
//
//  Create the list from basic list initialization data
//
//*******************************************************************************************************

ContactCommList::ContactCommList(const int  Num_Comm_Partners,
                                 const int* Num_to_Proc,
                                 const int* Comm_Proc_IDs,
                                 ContactTopologyEntity<Real>** Entity_List) :
  entity_index_list(NULL) {

  int i;
  num_comm_partners = Num_Comm_Partners;
  if( num_comm_partners ){
    comm_proc_ids = new int[num_comm_partners];
    std::memcpy( comm_proc_ids, Comm_Proc_IDs,
	    num_comm_partners*sizeof(int) );
    num_to_proc = new int[num_comm_partners];
    std::memcpy( num_to_proc, Num_to_Proc,
	    num_comm_partners*sizeof(int) );
    offset = new int[num_comm_partners];
    offset[0] = 0;
    for( i=1 ; i<num_comm_partners ; ++i )
      offset[i] = offset[i-1] + num_to_proc[i-1];
  } else {
    comm_proc_ids = NULL;
    num_to_proc = NULL;
    offset = NULL;
  }
  num_entities = 0;
  for( i=0 ; i<num_comm_partners ; ++i ) 
    num_entities += Num_to_Proc[i];
  if( num_entities ){
    entity_list = new ContactTopologyEntity<Real>*[num_entities];
    std::memcpy( entity_list, Entity_List,
	    num_entities*sizeof(ContactTopologyEntity<Real>*) );
  } else
    entity_list = NULL;
}


//*******************************************************************************************************
//
//  Allocate the list from basic list initialization data
//
//*******************************************************************************************************

ContactCommList::ContactCommList(const int  Num_Comm_Partners,
                                 const int* Num_to_Proc,
                                 const int* Comm_Proc_IDs) :
  entity_index_list(NULL) 
{
  num_comm_partners = Num_Comm_Partners;
  if( num_comm_partners ){
    comm_proc_ids = new int[num_comm_partners];
    std::memcpy( comm_proc_ids, Comm_Proc_IDs,
	    num_comm_partners*sizeof(int) );
    num_to_proc = new int[num_comm_partners];
    std::memcpy( num_to_proc, Num_to_Proc,
	    num_comm_partners*sizeof(int) );
    offset = NULL;
  } else {
    comm_proc_ids = NULL;
    num_to_proc = NULL;
    offset = NULL;
  }
  num_entities = 0;
  entity_list = NULL;
}

//*******************************************************************************************************
//
//  The communication list contains a list of communication entities.  However, it is much
//  faster to access data via straight indexes.  This routine calculates those indexes for
//  enforcement
// 
//*******************************************************************************************************

void ContactCommList::Set_Index_From_EnfArrayIndex() {
  if(entity_index_list) delete [] entity_index_list;
  entity_index_list = new int[num_entities];
  for(int i = 0; i < num_entities; ++i) {
    entity_index_list[i] = entity_list[i]->EnfArrayIndex();
  }
}


//*******************************************************************************************************
//
//  Destroy any alloacted data
//
//*******************************************************************************************************

ContactCommList::~ContactCommList(){
  delete [] comm_proc_ids;
  delete [] num_to_proc;
  delete [] offset;
  delete [] entity_list;
  if(entity_index_list) delete [] entity_index_list;
}

#endif
