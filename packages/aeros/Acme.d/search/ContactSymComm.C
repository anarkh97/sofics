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


#include "ContactSymComm.h"
#include "ContactTopologyEntity.h"
#include "ContactTopologyEntityHash.h"
#include "ContactTopologyEntityList.h"
#include "ContactNode.h"
#include "ContactTopology.h"
#include "Contact_Communication.h"
#include <cstring>

#ifndef CONTACT_NO_MPI

ContactSymComm::ContactSymComm( int Num_Comm_Partners,
				const int* Comm_Proc_Ids,
				const int* Num_to_Proc,
				ContactTopologyEntity<Real>** Entity_List )
{
  PRECONDITION( Num_Comm_Partners >= 0 );
  PRECONDITION( Comm_Proc_Ids || !Num_Comm_Partners );
  PRECONDITION( Num_to_Proc   || !Num_Comm_Partners );
  PRECONDITION( Entity_List   || !Num_Comm_Partners );

  num_comm_partners = Num_Comm_Partners;
  if( Num_Comm_Partners ){
    allocated_procs = num_comm_partners;
    comm_proc_ids = new int[num_comm_partners];
    std::memcpy( comm_proc_ids, Comm_Proc_Ids, sizeof(int)*num_comm_partners );
    num_to_proc = new int[num_comm_partners];
    std::memcpy( num_to_proc, Num_to_Proc, sizeof(int)*num_comm_partners );
    offset = new int[num_comm_partners];
    num_entities = 0;
    for( int i=0 ; i<num_comm_partners ; ++i ){
      offset[i] = num_entities;
      num_entities += num_to_proc[i];
    }
    allocated_entities = num_entities;
    entity_list = new ContactTopologyEntity<Real>*[num_entities];
    std::memcpy( entity_list, Entity_List, sizeof(ContactTopologyEntity<Real>*)*num_entities );
  } else {
    num_entities = 0;
    comm_proc_ids = NULL;
    num_to_proc = NULL;
    entity_list = NULL;
    offset = NULL;
  }
  // We have assumed in swapadd that the processor ids we communicate to are
  // in ascending order.
  //KHB: ADD CODING TO CHECK
}


ContactSymComm::ContactSymComm()
{

  allocated_procs = 0;
  allocated_entities = 0;
  num_comm_partners = 0;
  num_entities = 0;
  comm_proc_ids = NULL;
  num_to_proc = NULL;
  offset = NULL;
  entity_list = NULL;
}



ContactSymComm::~ContactSymComm()
{
  delete [] comm_proc_ids;
  delete [] num_to_proc;
  delete [] entity_list;
  delete [] offset;
}

void ContactSymComm::Build_Comm_Plan( ContactSymComm* GhostSymComm,
				      int max_shared_procs,
				      int size_shared_node_info,
				      int* shared_node_info,
				      ContactHostGlobalID* shared_node_gid,
				      ContactTopologyEntityHash* node_hash,
				      MPI_Comm& comm )
{
  // Problem: We have the original GhostSymComm which specifies all the
  //          communication in the primary decomposition.  We have now 
  //          created "phantom" copies of the nodes (the owning processor
  //          created the phantoms).  We have constructed another data
  //          structure which contains a full list of processors for every
  //          node which had a phantom created.  This has been communicated
  //          to all processors that have a ghost or phantom copy (a processor
  //          can't have both a phantom and a ghost copy of a node).  This
  //          data structure has the following
  //
  //             max_shared_procs  The maximum number of processors any node
  //                               is shared between.  This is a globally
  //                               consistent number.
  //
  //             size_shared_nodes
  //                               The number of nodes on this processor
  //                               (including phantoms) that have phantoms
  //                               on some processor.
  //                               
  //             shared_node_gid   A one-dimensional array (of size
  //                               size_shared_nodes) that has the global id
  //                               for each node in shared_node_info.
  //
  //             shared_node_info  A two-dimensional array (of size
  //                               [size_shared_node_info,max_shared_procs]
  //                               which has a list of processors that share
  //                               this node (this includes this processor).
  //                               A negative value indicates the end of the 
  //                               list for a given node (since it may be 
  //                               shared by less than the maximum number).
  //
  // Need:    We need to construct a new comm plan (in this object) that has
  //          both the ghost and phantom copies of the nodes represented so
  //          we can do swapadds.
  //
  // Notes:   The size_shared_node_info array is modified in this routine
  //          and temp_tag is used for temporary storage.
  //

  int i,j,k,l;

  int my_proc = contact_processor_number( comm );
  int num_procs = contact_number_of_processors( comm );

  // Remove my processor ID from node_sharing_info and compact the list.
  for( i=0 ; i<size_shared_node_info ; ++i ){
    REMEMBER( bool MY_PROC_FOUND = false );
    for( j=0 ; j<max_shared_procs ; ++j ){
      if( shared_node_info[i*max_shared_procs+j] == my_proc ){
	REMEMBER( MY_PROC_FOUND = true );
        for( k=j ; k<max_shared_procs-1 ; ++k ){
	  shared_node_info[i*max_shared_procs+k] = 
	    shared_node_info[i*max_shared_procs+k+1];
	}
	shared_node_info[i*max_shared_procs+max_shared_procs-1] = -12345;
	break;
      }
    }
    POSTCONDITION( MY_PROC_FOUND );
  }

  //
  //  Set temp_tag = -1234 for all nodes that are in the GhostSymComm
  //
  for( i=0 ; i<GhostSymComm->Num_Comm_Partners() ; ++i ){
    int n_to_proc = GhostSymComm->Num_to_Proc( i );
    ContactTopologyEntity<Real>** e_list = GhostSymComm->Entity_List( i );
    for( j=0 ; j<n_to_proc ; ++j ){
      e_list[j]->temp_tag = -1234;
    }
  }

  // Set temp_tag for each node in node_sharing_info to the index in 
  // node_sharing_info
  for( i=0 ; i<size_shared_node_info ; ++i ){
    ContactNode<Real>* node = static_cast<ContactNode<Real> *>
      (node_hash->find( shared_node_gid[i] ));
    POSTCONDITION( node );
    node->temp_tag = i;
  }

  // Now run the GhostSymComm and std::remove any node/proc pair that exists in the 
  // GhostSymComm.
  for( i=0 ; i<GhostSymComm->Num_Comm_Partners() ; ++i ){
    int n_to_proc = GhostSymComm->Num_to_Proc( i );
    ContactTopologyEntity<Real>** e_list = GhostSymComm->Entity_List( i );
    int comm_proc = GhostSymComm->Comm_Proc_ID( i );
    for( j=0 ; j<n_to_proc ; ++j ){
      if( e_list[j]->temp_tag != -1234 ){
	int index = e_list[j]->temp_tag;
	PRECONDITION( index >= 0   &&   index < size_shared_node_info );
	// Remove the proc in this comm plan from node_sharing_info
	REMEMBER( bool Found = false );
	for( k=0 ; k<max_shared_procs ; ++k ){
	  if( comm_proc == shared_node_info[index*max_shared_procs+k] ){
	    REMEMBER( Found = true );
	    for( l=k ; l<max_shared_procs-1 ; ++l ){
	      shared_node_info[index*max_shared_procs+l] = 
		shared_node_info[index*max_shared_procs+l+1];
	    }
	    shared_node_info[index*max_shared_procs+max_shared_procs-1] = 
	      -12345;
	    break;
	  }
	}
	POSTCONDITION( Found );
      }
    }
  }

  // Now we have removed all duplicate data (i.e., a node/proc pair is only 
  // listed once in the union of node_sharing_info & GhostSymComm.  Determine 
  // the size ofthe arrays.
  int* num_to_each_proc = new int[num_procs];
  std::memset( num_to_each_proc, 0, num_procs*sizeof(int) );
  num_entities = GhostSymComm->Size();
  for( i=0 ; i<size_shared_node_info ; ++i ){
    // There MUST be at least one proc not in the GhostSymComm for each node
    PRECONDITION( shared_node_info[i*max_shared_procs] >= 0 );
    for( j=0 ; j<max_shared_procs ; ++j ){
      int proc_id = shared_node_info[i*max_shared_procs+j];
      if( proc_id >= 0 ){
	++num_entities;
	// If we haven't seen this proc yet, add it
	if( num_to_each_proc[proc_id] == 0 ) ++num_comm_partners;
	++num_to_each_proc[proc_id];
      } else
	break;
    }
  }
  for( i=0 ; i<GhostSymComm->Num_Comm_Partners() ; ++i ){
    int proc_id = GhostSymComm->Comm_Proc_ID( i );
    POSTCONDITION( proc_id >= 0);
    // If we haven't seen this proc yet, add one to the number of comm partners
    if( num_to_each_proc[proc_id] == 0 ) ++num_comm_partners;
    num_to_each_proc[proc_id] += GhostSymComm->Num_to_Proc( i );
  }
  
  if( !num_comm_partners ) {
    delete [] num_to_each_proc;
    return;
  }   

  // Allocate the arrays
  comm_proc_ids = new int[num_comm_partners];
  // num_to_each_proc is already allocated and large enough so use it
  num_to_proc = num_to_each_proc;
  offset = new int[num_comm_partners];
  entity_list = new ContactTopologyEntity<Real>*[num_entities];
  
  // Fill in the processor and the offset arrays 
  int proc_index = 0;
  int count_nodes = 0;
  for( i=0 ; i<num_procs ; ++i ){
    if( num_to_proc[i] > 0 ){
      PRECONDITION( i != my_proc );
      comm_proc_ids[proc_index] = i;
      num_to_proc[proc_index] = num_to_each_proc[i];
      offset[proc_index] = count_nodes;
      count_nodes += num_to_each_proc[i];
      ++proc_index;
    }
  }
  POSTCONDITION( proc_index == num_comm_partners );
  POSTCONDITION( count_nodes == num_entities );

  // Build the entity list.  Note: don't worry about sorting them at this stage.
  // We will sort them as the final step with the Sort() member function.
  //
  // Just copy the GhostSymComm data for starters
  std::memset( entity_list, 0, num_entities*sizeof(ContactTopologyEntity<Real>*) );
  for( i=0 ; i<GhostSymComm->Num_Comm_Partners() ; ++i ){
    int proc_id = GhostSymComm->Comm_Proc_ID( i );
    // find this processor in the list
    REMEMBER( bool PROC_FOUND = false );
    for( j=0 ; j<num_comm_partners ; ++j ){
      if( comm_proc_ids[j] == proc_id ){
	//	PRECONDITION((offset[j+1]-offset[j])>=GhostSymComm->Num_to_Proc( i ) );
	std::memcpy( entity_list+offset[j], GhostSymComm->Entity_List( i ),
		GhostSymComm->Num_to_Proc( i )*sizeof(ContactTopologyEntity<Real>*) );
	REMEMBER( PROC_FOUND = true );
	break;
      }
    }
    POSTCONDITION( PROC_FOUND );
  }
  //
  // Now add the phantoms
  //
  for( i=0 ; i<size_shared_node_info ; ++i ){
    // There MUST be at least one proc not in the GhostSymComm for each node
    PRECONDITION( shared_node_info[i*max_shared_procs] >= 0 );
    for( j=0 ; j<max_shared_procs ; ++j ){
      int proc_id = shared_node_info[i*max_shared_procs+j];
      if( proc_id < 0 ) break;
      REMEMBER( bool PROC_FOUND = false );
      for( k=0 ; k<num_comm_partners ; ++k ){
	if( comm_proc_ids[k] == proc_id ){
	  REMEMBER( PROC_FOUND = true );
	  // find the location to insert this node 
	  // (NULL as used as unfilled locations)
	  REMEMBER( bool INSERTED = false );
	  ContactNode<Real>* node = static_cast<ContactNode<Real> *>
	    (node_hash->find( shared_node_gid[i]));
	  POSTCONDITION( node );
	  int begin = offset[k];
	  int end;
	  if (k+1!=num_comm_partners) end = offset[k+1];
	  else end = num_entities;
	  for( l=begin ; l<end ; ++l ){
	    // The node better not be in the list already
	    PRECONDITION( entity_list[l] != node );
	    if( entity_list[l] == NULL ){
	      entity_list[l] = node;
	      REMEMBER( INSERTED = true );
	      break;
	    }
	  }
	  POSTCONDITION( INSERTED );
	}
      }
      POSTCONDITION( PROC_FOUND );
    }
  }

  // Now we added everything so all that remains is to sort the list
  Sort();
}

void ContactSymComm::Build_Comm_Plan( ContactSymComm* GhostSymComm,
				      int max_shared_procs,
				      int size_shared_node_info,
				      int* shared_node_info,
				      ContactHostGlobalID* shared_node_gid,
				      ContactTopologyEntityList* node_hash,
				      MPI_Comm& comm )
{
  // Problem: We have the original GhostSymComm which specifies all the
  //          communication in the primary decomposition.  We have now 
  //          created "phantom" copies of the nodes (the owning processor
  //          created the phantoms).  We have constructed another data
  //          structure which contains a full list of processors for every
  //          node which had a phantom created.  This has been communicated
  //          to all processors that have a ghost or phantom copy (a processor
  //          can't have both a phantom and a ghost copy of a node).  This
  //          data structure has the following
  //
  //             max_shared_procs  The maximum number of processors any node
  //                               is shared between.  This is a globally
  //                               consistent number.
  //
  //             size_shared_nodes
  //                               The number of nodes on this processor
  //                               (including phantoms) that have phantoms
  //                               on some processor.
  //                               
  //             shared_node_gid   A one-dimensional array (of size
  //                               size_shared_nodes) that has the global id
  //                               for each node in shared_node_info.
  //
  //             shared_node_info  A two-dimensional array (of size
  //                               [size_shared_node_info,max_shared_procs]
  //                               which has a list of processors that share
  //                               this node (this includes this processor).
  //                               A negative value indicates the end of the 
  //                               list for a given node (since it may be 
  //                               shared by less than the maximum number).
  //
  // Need:    We need to construct a new comm plan (in this object) that has
  //          both the ghost and phantom copies of the nodes represented so
  //          we can do swapadds.
  //
  // Notes:   The size_shared_node_info array is modified in this routine
  //          and temp_tag is used for temporary storage.
  //

  int i,j,k,l;

  int my_proc = contact_processor_number( comm );
  int num_procs = contact_number_of_processors( comm );

  // Remove my processor ID from node_sharing_info and compact the list.
  for( i=0 ; i<size_shared_node_info ; ++i ){
    REMEMBER( bool MY_PROC_FOUND = false );
    for( j=0 ; j<max_shared_procs ; ++j ){
      if( shared_node_info[i*max_shared_procs+j] == my_proc ){
	REMEMBER( MY_PROC_FOUND = true );
        for( k=j ; k<max_shared_procs-1 ; ++k ){
	  shared_node_info[i*max_shared_procs+k] = 
	    shared_node_info[i*max_shared_procs+k+1];
	}
	shared_node_info[i*max_shared_procs+max_shared_procs-1] = -12345;
	break;
      }
    }
    POSTCONDITION( MY_PROC_FOUND );
  }

  // Set temp_tag = -1234 for all nodes that are in the GhostSymComm
  for( i=0 ; i<GhostSymComm->Num_Comm_Partners() ; ++i ){
    int n_to_proc = GhostSymComm->Num_to_Proc( i );
    ContactTopologyEntity<Real>** e_list = GhostSymComm->Entity_List( i );
    for( j=0 ; j<n_to_proc ; ++j ){
      e_list[j]->temp_tag = -1234;
    }
  }

  // Set temp_tag for each node in node_sharing_info to the index in 
  // node_sharing_info
  for( i=0 ; i<size_shared_node_info ; ++i ){
    ContactNode<Real>* node = static_cast<ContactNode<Real> *>
      (node_hash->Find( shared_node_gid[i] ) );
    POSTCONDITION( node );
    node->temp_tag = i;
  }

  // Now run the GhostSymComm and std::remove any node/proc pair that exists in the 
  // GhostSymComm.
  for( i=0 ; i<GhostSymComm->Num_Comm_Partners() ; ++i ){
    int n_to_proc = GhostSymComm->Num_to_Proc( i );
    ContactTopologyEntity<Real>** e_list = GhostSymComm->Entity_List( i );
    int comm_proc = GhostSymComm->Comm_Proc_ID( i );
    for( j=0 ; j<n_to_proc ; ++j ){
      if( e_list[j]->temp_tag != -1234 ){
	int index = e_list[j]->temp_tag;
	PRECONDITION( index >= 0   &&   index < size_shared_node_info );
	// Remove the proc in this comm plan from node_sharing_info
	REMEMBER( bool Found = false );
	for( k=0 ; k<max_shared_procs ; ++k ){
	  if( comm_proc == shared_node_info[index*max_shared_procs+k] ){
	    REMEMBER( Found = true );
	    for( l=k ; l<max_shared_procs-1 ; ++l ){
	      shared_node_info[index*max_shared_procs+l] = 
		shared_node_info[index*max_shared_procs+l+1];
	    }
	    shared_node_info[index*max_shared_procs+max_shared_procs-1] = 
	      -12345;
	    break;
	  }
	}
	POSTCONDITION( Found );
      }
    }
  }

  // Now we have removed all duplicate data (i.e., a node/proc pair is only 
  // listed once in the union of node_sharing_info & GhostSymComm.  Determine 
  // the size ofthe arrays.
  int* num_to_each_proc = new int[num_procs];
  std::memset( num_to_each_proc, 0, num_procs*sizeof(int) );
  num_entities = GhostSymComm->Size();
  for( i=0 ; i<size_shared_node_info ; ++i ){
    // There MUST be at least one proc not in the GhostSymComm for each node
    PRECONDITION( shared_node_info[i*max_shared_procs] >= 0 );
    for( j=0 ; j<max_shared_procs ; ++j ){
      int proc_id = shared_node_info[i*max_shared_procs+j];
      if( proc_id >= 0 ){
	++num_entities;
	// If we haven't seen this proc yet, add it
	if( num_to_each_proc[proc_id] == 0 ) ++num_comm_partners;
	++num_to_each_proc[proc_id];
      } else
	break;
    }
  }
  for( i=0 ; i<GhostSymComm->Num_Comm_Partners() ; ++i ){
    int proc_id = GhostSymComm->Comm_Proc_ID( i );
    POSTCONDITION( proc_id >= 0);
    // If we haven't seen this proc yet, add one to the number of comm partners
    if( num_to_each_proc[proc_id] == 0 ) ++num_comm_partners;
    num_to_each_proc[proc_id] += GhostSymComm->Num_to_Proc( i );
  }
  
  if( !num_comm_partners ) {
    delete [] num_to_each_proc;
    return;
  }   

  // Allocate the arrays
  comm_proc_ids = new int[num_comm_partners];
  // num_to_each_proc is already allocated and large enough so use it
  num_to_proc = num_to_each_proc;
  offset = new int[num_comm_partners];
  entity_list = new ContactTopologyEntity<Real>*[num_entities];
  
  // Fill in the processor and the offset arrays 
  int proc_index = 0;
  int count_nodes = 0;
  for( i=0 ; i<num_procs ; ++i ){
    if( num_to_proc[i] > 0 ){
      PRECONDITION( i != my_proc );
      comm_proc_ids[proc_index] = i;
      num_to_proc[proc_index] = num_to_each_proc[i];
      offset[proc_index] = count_nodes;
      count_nodes += num_to_each_proc[i];
      ++proc_index;
    }
  }
  POSTCONDITION( proc_index == num_comm_partners );
  POSTCONDITION( count_nodes == num_entities );

  // Build the entity list.  Note: don't worry about sorting them at this stage.
  // We will sort them as the final step with the Sort() member function.
  //
  // Just copy the GhostSymComm data for starters
  std::memset( entity_list, 0, num_entities*sizeof(ContactTopologyEntity<Real>*) );
  for( i=0 ; i<GhostSymComm->Num_Comm_Partners() ; ++i ){
    int proc_id = GhostSymComm->Comm_Proc_ID( i );
    // find this processor in the list
    REMEMBER( bool PROC_FOUND = false );
    for( j=0 ; j<num_comm_partners ; ++j ){
      if( comm_proc_ids[j] == proc_id ){
	//	PRECONDITION((offset[j+1]-offset[j])>=GhostSymComm->Num_to_Proc( i ) );
	std::memcpy( entity_list+offset[j], GhostSymComm->Entity_List( i ),
		GhostSymComm->Num_to_Proc( i )*sizeof(ContactTopologyEntity<Real>*) );
	REMEMBER( PROC_FOUND = true );
	break;
      }
    }
    POSTCONDITION( PROC_FOUND );
  }
  //
  // Now add the phantoms
  for( i=0 ; i<size_shared_node_info ; ++i ){
    // There MUST be at least one proc not in the GhostSymComm for each node
    PRECONDITION( shared_node_info[i*max_shared_procs] >= 0 );
    for( j=0 ; j<max_shared_procs ; ++j ){
      int proc_id = shared_node_info[i*max_shared_procs+j];
      if( proc_id < 0 ) break;
      REMEMBER( bool PROC_FOUND = false );
      for( k=0 ; k<num_comm_partners ; ++k ){
	if( comm_proc_ids[k] == proc_id ){
	  REMEMBER( PROC_FOUND = true );
	  // find the location to insert this node 
	  // (NULL as used as unfilled locations)
	  REMEMBER( bool INSERTED = false );
	  ContactNode<Real>* node = static_cast<ContactNode<Real> *>
	    (node_hash->Find( shared_node_gid[i] ));
	  POSTCONDITION( node );
	  int begin = offset[k];
	  int end;
	  if (k+1!=num_comm_partners) end = offset[k+1];
	  else end = num_entities;
	  for( l=begin ; l<end ; ++l ){
	    // The node better not be in the list already
	    PRECONDITION( entity_list[l] != node );
	    if( entity_list[l] == NULL ){
	      entity_list[l] = node;
	      REMEMBER( INSERTED = true );
	      break;
	    }
	  }
	  POSTCONDITION( INSERTED );
	}
      }
      POSTCONDITION( PROC_FOUND );
    }
  }

  // Now we added everything so all that remains is to sort the list
  Sort();
}


void ContactSymComm::Build_Subset_Comm_Plan_Using_Temp_Tag( ContactSymComm& SymComm)
{
  // This function builds a subset comm plan using temp_tag as a true (1) or
  // false (0) flag if the entity should be added.

  // We can allocate all the data to the size of the original comm plan since we
  // know the subset plan will either be smaller or equal in size (it can never
  // be larger).

  allocated_procs = SymComm.Num_Comm_Partners();
  if( allocated_procs ){
    comm_proc_ids = new int[allocated_procs];
    std::memset(comm_proc_ids, 0, allocated_procs*sizeof(int));
    num_to_proc = new int[allocated_procs];
    std::memset(num_to_proc, 0, allocated_procs*sizeof(int));
    offset = new int[allocated_procs];
    std::memset(offset, 0, allocated_procs*sizeof(int));
  }
  allocated_entities = SymComm.Size();
  if( allocated_entities ){
    entity_list = new ContactTopologyEntity<Real>*[allocated_entities];
  } 

  for( int i=0 ; i<allocated_procs ; ++i ){
    int n_to_proc = SymComm.Num_to_Proc( i );
    ContactTopologyEntity<Real>** full_entity_list_to_proc = SymComm.Entity_List(i);
    bool this_proc_counted = false;
    for( int j=0 ; j<n_to_proc ; ++j ){
      if( full_entity_list_to_proc[j]->temp_tag == 1 ){
	if( !this_proc_counted ){
	  offset[num_comm_partners] = num_entities;
	  comm_proc_ids[num_comm_partners++] = SymComm.Comm_Proc_ID(i);
	  this_proc_counted = true;
	}
	num_to_proc[num_comm_partners-1] += 1;
	entity_list[num_entities++] = full_entity_list_to_proc[j];
      }
    }
  }
}


void ContactSymComm::Sort()
{
  // Sort the nodes to each processor by Global_ID.  This guarantees the
  // same ordering on all processors.
  for( int i=0 ; i<num_comm_partners ; ++i ){
    int n_to_proc = num_to_proc[i];
    ContactTopologyEntity<Real>** list = entity_list+offset[i];
    ContactTopology::SortEntityList( n_to_proc, list );
  }
}



void ContactSymComm::Check(MPI_Comm Comm) {

  int i;

  // This coding checks the SymComm to make sure that it is reasonable
  
  // ** Local Checks **
  //  Make sure comm procs are in ascending order and are in the std::right range.
  //  The assumption of ascending order is necessary for the current
  //  implementation of swapadd.

  REMEMBER(int num_procs = contact_number_of_processors(Comm));
#if CONTACT_DEBUG_PRINT_LEVEL>=1
  int my_proc   = contact_processor_number(Comm);
#endif

  for ( i = 1; i < num_comm_partners; ++i) {
    PRECONDITION(comm_proc_ids[i] >= 0 && comm_proc_ids[i] < num_procs);
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    if ( comm_proc_ids[i] <= comm_proc_ids[i-1] )
      std::cerr << "  WARNING: on processor " << my_proc 
	   << " comm proc id " << i << " = " << comm_proc_ids[i]
	   << " is less than previous comm proc: " << comm_proc_ids[i-1]
	   << std::endl;
#endif
      POSTCONDITION( comm_proc_ids[i] > comm_proc_ids[i-1] );
  }

  // Now check that individual elements in each list are in ascending order
  for ( i = 0; i < num_comm_partners; ++i) {
    for (int j = 1; j < num_to_proc[i]; ++j) {
#if CONTACT_DEBUG_PRINT_LEVEL>=1 || defined(CONTACT_DEBUG)
      ContactTopologyEntity<Real> * entity_1 = entity_list[offset[i]+j];
      ContactTopologyEntity<Real> * entity_2 = entity_list[offset[i]+j-1];
#endif
#if CONTACT_DEBUG_PRINT_LEVEL>=1
      if ( entity_1->Global_ID() < entity_2->Global_ID() ) 
	std::cerr << "  WARNING: on processor " << my_proc 
	     << " entity " << j << " in group " << i 
	     << " is less than previous entity. " << std::endl;
#endif
      POSTCONDITION( entity_1->Global_ID() > entity_2->Global_ID() );
    }
  }

  // ** Parallel Checks **
  //  nothing std::right now

}


#endif
