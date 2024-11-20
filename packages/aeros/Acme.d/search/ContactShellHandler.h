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


#ifndef ContactShellHandler_h_
#define ContactShellHandler_h_

class ContactErrors;
class ContactCommBuffer;
template<typename DataType> class ContactNode;
class ContactSymComm;
class ContactCommBuffer;
class ContactTopology;
class ContactBlockEntityList;
class ContactHostGlobalID;
class ContactFixedSizeAllocator;

#include "ContactSearch.h"
#include "contact_assert.h"
#include "ContactEntity.h"
#include "ContactParOStream.h"

#ifndef CONTACT_NO_MPI
#include "mpi.h"
class ContactSymComm;
class ContactAsymComm;
#endif

// ASG Temporary to do print_map
#include <iostream>
#include <fstream>

class ContactShellHandler {

 public:

  enum ShellID_Resolution{SEQUENTIAL_IDS, USE_OLD_HANDLER};

  ContactShellHandler( ContactErrors* Errors,
		       int Dimensionality, 
		       int number_node_blocks,
      const ContactSearch::ContactNode_Type* node_block_types, 
		       const int* host_number_nodes_per_block,
		       int*& number_nodes_per_block,
		       const int* host_exodus_node_ids,
		       int*& exodus_node_ids,
		       const int* host_global_node_ids,
		       int*& global_node_ids,
		       ContactType** node_entity_types,
		       const double* coords,
		       int number_face_blocks, 
      const ContactSearch::ContactFace_Type* face_block_types,
		       const int* number_faces_per_block, 
		       const int* global_face_ids,
		       const int* host_face_connectivity,
		       int*& face_connectivity,
		       const double* face_lofting_factors,
		       int number_comm_partners, 
		       const int* comm_proc_ids,
		       const int* host_number_nodes_to_partner, 
		       int*& number_nodes_to_partner, 
		       const int* host_comm_nodes,
		       int*& comm_nodes,
		       MPI_Comm& communicator,
		       ContactSearch* search,
		       ContactTopology* Topology,
		       ContactCommBuffer* CommBuffer,
		       ShellID_Resolution res_method,
		       ContactShellHandler * old_handler,
		       ContactSearch::ContactErrorCode& error_code );
  ~ContactShellHandler();

  void Set_NodeBlk_Position( int, ContactBlockEntityList*, 
                             VariableHandle, const Real* );
  void Set_NodeBlk_KinConstr( int, ContactBlockEntityList*, 
                              VariableHandle, VariableHandle, 
			      const int*, const Real* );
  void Set_NodeBlk_RemainingGap( int, ContactBlockEntityList*, 
                                 VariableHandle, const Real* );
  void Set_NodeBlk_GhostingGap( int, ContactBlockEntityList*, 
                                VariableHandle, const Real* );
 
  void Complete_Initialization();

  void Loft_Nodes( int num_configs, VariableHandle, VariableHandle,
		   VariableHandle, VariableHandle, VariableHandle );

  int Number_Host_Code_Nodes() { return number_host_code_nodes; };

  void GID_of_Host_Code_Node(int host_node, int & first, int & second) { 
    first = orig_global_node_ids[2*host_node];
    second = orig_global_node_ids[2*host_node+1];
    return; 
  };

  inline int Max_Num_Acme_Nodes_for_Host_Node() {
    return max_num_acme_nodes_for_host_code_node;
  };

  inline int Num_Acme_Nodes_for_Host_Node( int node ) {
    PRECONDITION (node < number_host_code_nodes && node >= 0);
    return num_acme_nodes_for_host_code_node[node];
  };
  
  inline int Acme_Node_for_Host_Node( int node, int num ) {
    PRECONDITION (node < number_host_code_nodes && node >= 0);
    PRECONDITION (num < Num_Acme_Nodes_for_Host_Node(node) && num >= 0);
    return host_to_acme_node_map[offset_for_node[node] + num];
  };
  
  void
  Acme_NodeGID_for_Host_Node( int node, int num,
                              int& hi, int& lo );
  
  inline int* Acme_to_Host_Node_Map() { return acme_to_host_node_map; };

  void NumberNodes();

  ContactParOStream postream;
  //-- use this to convert postream into std::ofstream; also change .C
  // std::ofstream postream;
  
 private:
  
  ContactShellHandler(ContactShellHandler&);
  ContactShellHandler& operator=(ContactShellHandler&);

  void Initialize_Node_Maps();
  void Create_Node( int host_id, int acme_id );
  void Renumber_Nodes( int number_node_blocks,
		       const int* host_number_nodes_per_block,
		       int* number_nodes_per_block,
		       const int& number_face_blocks,
		       const ContactSearch::ContactFace_Type* face_block_types,
		       const int* number_faces_per_block,
		       int* face_connectivity );
  void Complete_Node_Maps();

  void compute_inv_conn(const int & num_host_code_nodes,
			const int & number_face_blocks,
			const	int * face_connectivity,
			const	int * conn_offsets,
	          const ContactSearch::ContactFace_Type* face_block_types,
			const int * number_faces_per_block,
			const double * face_lofting_factors,
			const int& number_comm_partners,
			const int* comm_proc_ids,
			const int* number_nodes_to_partner,
			const int* comm_nodes,
			int *& inv_conn,
			int *& offset_inv_conn,
			int *& inv_conn_count,
			int *& shell_face_blk_list,
			int *& shell_face_id_list,
			double *& shell_loft_factors,
			int & num_shell_faces,
			int & num_shell_nodes,
			bool parallel_consistent,
			ContactSearch * search );

  bool process_face(const int shell_face,
		    bool * face_process_flag,
		    bool * face_node_process_flag,
		    const ContactSearch::ContactFace_Type* face_block_types,
		    int * face_conn_changeable,
		    int * face_conn_reference,
		    const int * conn_offsets,
		    const int * inv_conn,
		    const int * offset_inv_conn,
		    const int * inv_conn_count,
		    const int * shell_face_blk_list,
		    const int * shell_face_id_list,
		    const double * shell_loft_factors,
		    int * num_nds_per_blk,
		    const double * coords,
		    ContactSearch * search );

  int process_edge(const int shell_face,
		   const int anchor_node,
		   const int shell_node,
		   const int node1, 
		   const int node2,
		   const int * inv_conn,
		   const int * offset_inv_conn,
		   const int * inv_conn_count,
		   const int * conn_offsets,
		   const int * shell_face_blk_list,
		   const int * shell_face_id_list,
		   const double * shell_loft_factors,
		   int * face_conn_changeable,
		   int * face_conn_reference,
		   const ContactSearch::ContactFace_Type* face_block_types,
		   bool * face_node_process_flag,
		   const double * coords,
		   ContactSearch * search );
  
  bool find_opposite(const int shell_face,
		     const ContactSearch::ContactFace_Type* face_block_types,
		     const int * conn_offsets,
		     const int * face_connectivity,
		     const int * offset_inv_conn,
		     const int * inv_conn_count,
		     const int * inv_conn,
		     const int * shell_face_blk_list,
		     const int * shell_face_id_list,
		     ContactSearch * search);

  int find_best_face(int node1,
		     int node2,
		     int face,
		     const ContactSearch::ContactFace_Type* face_block_types,
		     const int * conn_offsets,
		     const int * face_connectivity,
		     const int * shell_face_blk_list,
		     const int * shell_face_id_list,
		     const double * shell_loft_factors,
		     int num_adj_faces,
		     const int * adj_faces,
		     const double * coords,
		     ContactSearch * search);
  
#ifndef CONTACT_NO_MPI
  ContactSearch::ContactErrorCode ghost_shell_faces(
			 const int& number_comm_partners,
			 const int* comm_proc_ids,
			 const int* number_nodes_to_partner,
			 const int* comm_nodes,
			 const int& num_orig_shell_faces,
			 const int& num_orig_face_blks,
			 const int* inv_conn,
			 const int* inv_conn_count,
			 const int* offset_inv_conn,
			 const int* face_connectivity,
			 const int* conn_offsets,
     const ContactSearch::ContactFace_Type* face_block_types,
			 const int* number_orig_faces_per_blk,
			 const int* global_face_ids,
			 const double* coords,
			 const int * shell_face_blk_list,
			 const int * shell_face_id_list,
			 const double * shell_loft_factors,
			 int  & full_number_face_blocks,
			 int *& full_face_connectivity,
			 int *& full_conn_offsets,
           ContactSearch::ContactFace_Type*& full_face_block_types,
			 int *& full_number_faces_per_blk,
			 int *& full_global_face_ids,
			 int & full_num_host_code_nodes,
			 double *& full_coords,
			 double *& full_shell_loft_factors,
			 ContactSearch* search);

  ContactSearch::ContactErrorCode
  parallel_node_resolution(const int* full_number_faces_per_block,
			   const int* full_updated_inv_conn,
			   const int* full_updated_offset_inv_conn,
			   const int* full_updated_inv_conn_count,
			   const int* full_updated_shell_face_blk,
			   const int* full_updated_shell_face_id,
			   const int* full_face_global_ids,
			   const int& number_comm_partners,
			   const int* comm_proc_ids,
			   const int* host_number_nodes_to_partner,
			   const int* comm_nodes,
			   int* number_nodes_to_partner,
			   int* new_comm_nodes,
			   int* exodus_node_ids,
			   int* global_node_ids,
			   ShellID_Resolution res_method,
			   ContactShellHandler * old_handler,
			   ContactSearch * search);

#endif // CONTACT_NO_MPI

  enum NodeStatus {OWNED_NOT_SHARED=0, OWNED_SHARED, NOT_OWNED_SHARED};

  void number_nodes_sequentially(const int num_owned_host_code_nodes,
				 const int num_owned_newly_created_nodes,
				 const NodeStatus * node_status,
				 int* exodus_node_ids,
				 int* global_node_ids);
  
  void number_nodes_using_handler(const int num_owned_host_code_nodes,
				  const int num_owned_newly_created_nodes,
				  const NodeStatus * node_status,
				  const int* inv_conn,
				  const int* inv_conn_count,
				  const int* offset_inv_conn,
				  const int* full_shell_face_blk,
				  const int* full_shell_face_id,
				  const int* number_faces_per_block,
				  const int* global_face_ids,
				  ContactShellHandler * old_handler,
				  int* exodus_node_ids,
				  int* global_node_ids);
  
  void Mark_Shell_Nodes();

  int m_number_faces;
  int m_full_number_faces;
  int number_acme_nodes;
  int* host_nodes_per_block;
  int* block_offset_for_nodes;
  int* offset_for_node;
  int* num_acme_nodes_for_host_code_node;
  int  max_num_acme_nodes_for_host_code_node;
  int* new_acme_host_node_map;
  int* host_to_acme_node_map;

  int dimensionality;
  int number_host_code_nodes;
  int* acme_to_host_node_map;
  int* num_created_shell_nodes;
  int next_new_id;
  int* opposite_face;
  int* shell_face_id;
  int* shell_face_blk;
  double* shell_loft_factors;
  bool* is_a_tab_node;

  int* orig_node_block_types;
  int* orig_number_nodes_per_block;
  int* orig_exodus_node_ids;
  int* orig_global_node_ids;
  int* orig_face_connectivity;
  int* orig_comm_proc_ids;
  int* orig_number_nodes_to_partner;
  int  orig_number_comm_nodes;
  int* orig_comm_nodes;

  int* new_number_nodes_per_block;
  int* new_exodus_node_ids;
  int* new_global_node_ids;
  int* new_face_connectivity;
  int* new_number_nodes_to_partner;
  int* new_comm_nodes;

  int number_lofted_nodes;
  ContactShellNode** lofted_nodes;
  bool lofting_computed;

  // shell node communication data structures
  int my_proc_id;
#ifndef CONTACT_NO_MPI
  ContactSymComm* comm_plan;
  ContactAsymComm* comm_plan_to_owner;
  ContactAsymComm* comm_plan_to_ghost;
  ContactAsymComm* import_loft_face_data_comm_plan;
  int* number_lofting_ghosts;
  int* offset_to_buffer_offsets;
  int* buffer_offsets;
  void Build_Shell_Comm_Plans();
  void Build_Loft_Nodes_Data_Structures();
#endif

  ContactTopology * topology;
  ContactCommBuffer * comm_buffer;
  MPI_Comm SearchComm;

  // Data and structures for holding created node in integer form until we
  // are ready to allocate the true arrays (which can't be done until we know
  // the total size).
  struct created_node {
    int host_id;
    int acme_id;
    created_node* next_node;
  };
  ContactFixedSizeAllocator* allocator;
  created_node* first_node;
  created_node* last_node;

};

#endif // #ifndef ContactShellHandler_h_
