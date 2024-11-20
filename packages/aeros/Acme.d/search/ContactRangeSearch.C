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

  
#include <ContactRangeSearch.h>
#include <ContactBoundingBox.h>
#include <ContactBoundingBoxHierarchy.h>
#include <ContactBoundingBoxHierarchy_Int.h>
#include <Contact_Communication.h>
#include <Contact_Defines.h>
#include <ContactVector.h>
#include <ContactTopology.h>
#include <ContactNodeFaceInteraction.h>
#include <ContactSearch.h>
#include <ContactSearchData.h>
#include "search_methods.h"
#include "contact_sorting.h"

using namespace std;

ContactNodeFaceRangeSearch::ContactNodeFaceRangeSearch(SearchType my_type, 
                                                       ContactTopology *my_topo,
                                                       const int num_active_nodes_,
                                                       int *node_search_status,
                                                       ContactSearch *my_search,
                                                       VariableHandle &NODE_COORD_START_,
                                                       VariableHandle &NODE_COORD_END_) :
  node_hierarchy(NULL),
  node_block_BB_hierarchy(NULL),
  node_hierarchy_ptrs(NULL),
  node_key_array(NULL),
  scratch(NULL),
  index(NULL),
  rank(NULL),
  rank2(NULL),
  position(NULL),
  node_map(NULL),  
  search_topology(my_topo),
  search(my_search),
  search_type(my_type),
  num_active_nodes(num_active_nodes_),
  NODE_COORD_START(NODE_COORD_START_), 
  NODE_COORD_END(NODE_COORD_END_) {
  //
  //  Save useful pointer and variable references
  //
  Nodes = reinterpret_cast<ContactNode<Real>**>(search_topology->NodeList()->EntityList());
  NODE_NORMAL        = search_topology->Variable_Handle( ContactTopology::Node_Normal        );
  FACE_NORMAL        = search_topology->Variable_Handle( ContactTopology::Face_Normal        );
  REMAINING_GAP      = search_topology->Variable_Handle( ContactTopology::Remaining_Gap      );
  //
  //  Branch on search type, set up the required data for the 
  //  search based on search type
  //
  int number_of_nodes = search_topology->Number_of_Nodes();
  
  int num_configs = 1;
  if (NODE_COORD_START!=NODE_COORD_END) num_configs = 2;

  if(search_type == BASIC_SEARCH) {
    //
    //  Build a bounding box hierarchy containing all active nodes
    //
    int  auto_tol   = search->AutoTol();
    Real box_scale  = search->BoxInflation();
    ObjectBoundingBox *node_boxes = new ObjectBoundingBox[num_active_nodes];
    int cur_node_box = 0;
    for (int inode=0; inode<number_of_nodes; ++inode) { 
      if( node_search_status[inode] ) {
        ContactNode<Real> *node = Nodes[inode];
        node_boxes[cur_node_box].set_object_number(inode);
        node->ComputeBoundingBoxForSearch(num_configs,
                                          REMAINING_GAP,
                                          NODE_COORD_START,
                                          NODE_COORD_END,
                                          auto_tol,
                                          box_scale,
                                          node_boxes[cur_node_box]);
	++cur_node_box;
      }
    }

    int hierarchy_size = 2 * num_active_nodes - 1;
    if(hierarchy_size > 0) {
      node_hierarchy = new ObjectBoundingBoxHierarchy[hierarchy_size];
      ObjectBoundingBoxHierarchy::create_hierarchy(node_hierarchy, node_boxes, num_active_nodes);
    }
    delete [] node_boxes;
  } else if(search_type == KEY_AWARE_SEARCH) {
    //
    //  Build a bounding box hierarchy containing all active nodes
    //
    int  auto_tol  = search->AutoTol();
    Real box_scale = search->BoxInflation();
    ObjectBoundingBox *node_boxes = new ObjectBoundingBox[num_active_nodes];
    int cur_node_box = 0;
    node_entity_keys.resize(number_of_nodes, -2);
    for (int inode=0; inode<number_of_nodes; ++inode) { 
      if( node_search_status[inode] ) {
	ContactNode<Real> *node = Nodes[inode];
	node_entity_keys[inode] = node->Get_Owning_Entity();
	node_boxes[cur_node_box].set_object_number(inode);
        node->ComputeBoundingBoxForSearch(num_configs,
                                          REMAINING_GAP,
                                          NODE_COORD_START,
                                          NODE_COORD_END,
                                          auto_tol,
                                          box_scale,
                                          node_boxes[cur_node_box]);
	++cur_node_box;
      }
    }
    int hierarchy_size = 2 * num_active_nodes - 1;
    if(hierarchy_size > 0) {
      node_hierarchy = new ObjectBoundingBoxHierarchy[hierarchy_size];
      ObjectBoundingBoxHierarchy::create_hierarchy(node_hierarchy, node_boxes, num_active_nodes);
    }
    delete [] node_boxes;

    //
    //  Calculate owned node keys for each node of the tree
    //
    node_key_array = new int[hierarchy_size]; 
    node_hierarchy->calc_keys_recurse(node_key_array, node_entity_keys);
  } else if(search_type == ENTITY_SEARCH) {
    //
    //  Build one bounding box hierarchy for each contact entity block.  A node will either have its
    //  own entity key or will be contained in the entity of a face the node is attached to.  
    //
    //  If a node has faces that belong to multiple entities, add that node to an additional special node tree
    //  (at location[0]).
    //  It is assumed that the case of nodes belonging to multiple entities is rare, thus if these nodes
    //  do occur they are processed seperately from the other nodes to avoid slowing down the far more common,
    //  one entity per node case.
    //
    //  First build entity sorted lists of contact nodes
    //
    int num_search_entities = search->Search_Data()->Num_Search_Entities();

    vector< vector<ContactNode<Real>*> > node_block_entities(num_search_entities + 1);
    vector< vector<int> >          node_global_index(num_search_entities + 1);

    for (int inode=0; inode<number_of_nodes; ++inode) { 
      if( node_search_status[inode] ) {
	ContactNode<Real> *node = Nodes[inode];
	int node_key = node->Get_Owning_Entity();
	if(node_key >= 0) {
	  node_block_entities[node_key+1].push_back(node);
	  node_global_index[node_key+1].push_back(inode);
	} else {
	  node_block_entities[0].push_back(node);
	  node_global_index[0].push_back(inode);
	}
      }
    }
    //
    //  Determine the maximum number of nodes for any block list, and predict the
    //  total size of all bounding box hierarchy
    //
    int max_node_list_size = 0;
    int total_hierarchy_size = 0;
    for(int iblock = 0; iblock < num_search_entities+1; ++iblock) {
      const int num_nodes = node_block_entities[iblock].size();
      if(num_nodes > 0) {
	total_hierarchy_size += num_nodes*2 - 1;
	if(num_nodes > max_node_list_size) max_node_list_size = num_nodes;
      }
    } 
    //
    //  Allocate the bounding box hierarchy lists
    //
    node_hierarchy                = new ObjectBoundingBoxHierarchy[total_hierarchy_size];
    node_hierarchy_ptrs           = new ObjectBoundingBoxHierarchy*[num_search_entities+1];
    ObjectBoundingBox *node_boxes = new ObjectBoundingBox[max_node_list_size];
    ACME::ObjectBoundingBox_Vector   node_block_boxes;
    //
    //  Build a bounding box hierarchy for each node list
    //
    int  auto_tol  = search->AutoTol();
    Real box_scale = search->BoxInflation();
    node_hierarchy_ptrs[0] = node_hierarchy;
    for(int iblock = 0; iblock < num_search_entities+1; ++iblock) {
      vector<ContactNode<Real>*> &nodes = node_block_entities[iblock];
      const int num_nodes = nodes.size();
      //
      //  Step 1, compute the bounding box for each node in the list
      //
      for(int inode = 0; inode < num_nodes; ++inode) {
	ContactNode<Real> *node = nodes[inode];
        node->ComputeBoundingBoxForSearch(num_configs,
                                          REMAINING_GAP,
                                          NODE_COORD_START,
                                          NODE_COORD_END,
                                          auto_tol,
                                          box_scale,
                                          node_boxes[inode]);
	node_boxes[inode].set_object_number(node_global_index[iblock][inode]);
      }
      //
      //  Step 2, create the nodal bounding box hierarchy for the current entity
      //
      if(num_nodes > 0) { 
	ObjectBoundingBoxHierarchy::create_hierarchy(node_hierarchy_ptrs[iblock], node_boxes, num_nodes);
	if(iblock+1 <= num_search_entities) node_hierarchy_ptrs[iblock+1] = node_hierarchy_ptrs[iblock] + num_nodes * 2 - 1;
	//
	//  Add the bounding box for the current entity to the list of entity bounding boxes
	//
	ObjectBoundingBox node_block_box(*(node_hierarchy_ptrs[iblock]), iblock);
	node_block_boxes.push_back(node_block_box);
      } else {
	if(iblock+1 <= num_search_entities) node_hierarchy_ptrs[iblock+1] = node_hierarchy_ptrs[iblock];
      }
      if(num_nodes == 0) node_hierarchy_ptrs[iblock] = NULL;
    }
    delete [] node_boxes;
    //
    //  Build a bounding box hierarchy of the nodal entity bounding boxes
    //
    if(node_block_boxes.size() > 0) {
      node_block_BB_hierarchy = new ObjectBoundingBoxHierarchy[node_block_boxes.size() * 2 - 1];
      ObjectBoundingBoxHierarchy::create_hierarchy(node_block_BB_hierarchy, node_block_boxes.get_buffer(), node_block_boxes.size());
    }
  }
}

ContactNodeFaceRangeSearch::ContactNodeFaceRangeSearch(SearchType my_type, 
                                                       ContactTopology *my_topo,
                                                       const int num_active_nodes_,
                                                       int *node_search_status,
                                                       ContactSearch *my_search,
                                                       VariableHandle &NODE_COORD_START_,
                                                       VariableHandle &NODE_COORD_END_,
                                                       Real capture_tol) :
  node_hierarchy(NULL),
  node_block_BB_hierarchy(NULL),
  node_hierarchy_ptrs(NULL),
  node_key_array(NULL),
  scratch(NULL),
  index(NULL),
  rank(NULL),
  rank2(NULL),
  position(NULL),
  node_map(NULL),  
  search_topology(my_topo),
  search(my_search),
  search_type(my_type),
  num_active_nodes(num_active_nodes_),
  NODE_COORD_START(NODE_COORD_START_), 
  NODE_COORD_END(NODE_COORD_END_) {
  //
  //  Save useful pointer and variable references
  //
  Nodes = reinterpret_cast<ContactNode<Real>**>(search_topology->NodeList()->EntityList());
  NODE_NORMAL        = search_topology->Variable_Handle( ContactTopology::Node_Normal        );
  FACE_NORMAL        = search_topology->Variable_Handle( ContactTopology::Face_Normal        );
  REMAINING_GAP      = search_topology->Variable_Handle( ContactTopology::Remaining_Gap      );
  //
  //  Branch on search type, set up the required data for the 
  //  search based on search type
  //
  int number_of_nodes = search_topology->Number_of_Nodes();
  
  int num_configs = 1;
  if (NODE_COORD_START!=NODE_COORD_END) num_configs = 2;

  if(search_type == BASIC_SEARCH) {
    //
    //  Build a bounding box hierarchy containing all active nodes
    //
    int  auto_tol   = search->AutoTol();
    Real box_scale  = search->BoxInflation();
    ObjectBoundingBox *node_boxes = new ObjectBoundingBox[num_active_nodes];
    int cur_node_box = 0;
    for (int inode=0; inode<number_of_nodes; ++inode) { 
      if( node_search_status[inode] ) {
        ContactNode<Real> *node = Nodes[inode];
        node_boxes[cur_node_box].set_object_number(inode);
        node->ComputeBoundingBoxForSearch(num_configs,
                                          REMAINING_GAP,
                                          NODE_COORD_START,
                                          NODE_COORD_END,
                                          auto_tol,
                                          box_scale,
                                          capture_tol,
                                          node_boxes[cur_node_box]);
	++cur_node_box;
      }
    }

    int hierarchy_size = 2 * num_active_nodes - 1;
    if(hierarchy_size > 0) {
      node_hierarchy = new ObjectBoundingBoxHierarchy[hierarchy_size];
      ObjectBoundingBoxHierarchy::create_hierarchy(node_hierarchy, node_boxes, num_active_nodes);
    }
    delete [] node_boxes;
  } else if(search_type == KEY_AWARE_SEARCH) {
    //
    //  Build a bounding box hierarchy containing all active nodes
    //
    int  auto_tol  = search->AutoTol();
    Real box_scale = search->BoxInflation();
    ObjectBoundingBox *node_boxes = new ObjectBoundingBox[num_active_nodes];
    int cur_node_box = 0;
    node_entity_keys.resize(number_of_nodes, -2);
    for (int inode=0; inode<number_of_nodes; ++inode) { 
      if( node_search_status[inode] ) {
	ContactNode<Real> *node = Nodes[inode];
	node_entity_keys[inode] = node->Get_Owning_Entity();
	node_boxes[cur_node_box].set_object_number(inode);
        node->ComputeBoundingBoxForSearch(num_configs,
                                          REMAINING_GAP,
                                          NODE_COORD_START,
                                          NODE_COORD_END,
                                          auto_tol,
                                          box_scale,
                                          capture_tol,
                                          node_boxes[cur_node_box]);
	++cur_node_box;
      }
    }
    int hierarchy_size = 2 * num_active_nodes - 1;
    if(hierarchy_size > 0) {
      node_hierarchy = new ObjectBoundingBoxHierarchy[hierarchy_size];
      ObjectBoundingBoxHierarchy::create_hierarchy(node_hierarchy, node_boxes, num_active_nodes);
    }
    delete [] node_boxes;

    //
    //  Calculate owned node keys for each node of the tree
    //
    node_key_array = new int[hierarchy_size]; 
    node_hierarchy->calc_keys_recurse(node_key_array, node_entity_keys);
  } else if(search_type == ENTITY_SEARCH) {
    //
    //  Build one bounding box hierarchy for each contact entity block.  A node will either have its
    //  own entity key or will be contained in the entity of a face the node is attached to.  
    //
    //  If a node has faces that belong to multiple entities, add that node to an additional special node tree
    //  (at location[0]).
    //  It is assumed that the case of nodes belonging to multiple entities is rare, thus if these nodes
    //  do occur they are processed seperately from the other nodes to avoid slowing down the far more common,
    //  one entity per node case.
    //
    //  First build entity sorted lists of contact nodes
    //
    int num_search_entities = search->Search_Data()->Num_Search_Entities();

    vector< vector<ContactNode<Real>*> > node_block_entities(num_search_entities + 1);
    vector< vector<int> >          node_global_index(num_search_entities + 1);

    for (int inode=0; inode<number_of_nodes; ++inode) { 
      if( node_search_status[inode] ) {
	ContactNode<Real> *node = Nodes[inode];
	int node_key = node->Get_Owning_Entity();
	if(node_key >= 0) {
	  node_block_entities[node_key+1].push_back(node);
	  node_global_index[node_key+1].push_back(inode);
	} else {
	  node_block_entities[0].push_back(node);
	  node_global_index[0].push_back(inode);
	}
      }
    }
    //
    //  Determine the maximum number of nodes for any block list, and predict the
    //  total size of all bounding box hierarchy
    //
    int max_node_list_size = 0;
    int total_hierarchy_size = 0;
    for(int iblock = 0; iblock < num_search_entities+1; ++iblock) {
      const int num_nodes = node_block_entities[iblock].size();
      if(num_nodes > 0) {
	total_hierarchy_size += num_nodes*2 - 1;
	if(num_nodes > max_node_list_size) max_node_list_size = num_nodes;
      }
    } 
    //
    //  Allocate the bounding box hierarchy lists
    //
    node_hierarchy                = new ObjectBoundingBoxHierarchy[total_hierarchy_size];
    node_hierarchy_ptrs           = new ObjectBoundingBoxHierarchy*[num_search_entities+1];
    ObjectBoundingBox *node_boxes = new ObjectBoundingBox[max_node_list_size];
    ACME::ObjectBoundingBox_Vector   node_block_boxes;
    //
    //  Build a bounding box hierarchy for each node list
    //
    int  auto_tol  = search->AutoTol();
    Real box_scale = search->BoxInflation();
    node_hierarchy_ptrs[0] = node_hierarchy;
    for(int iblock = 0; iblock < num_search_entities+1; ++iblock) {
      vector<ContactNode<Real>*> &nodes = node_block_entities[iblock];
      const int num_nodes = nodes.size();
      //
      //  Step 1, compute the bounding box for each node in the list
      //
      for(int inode = 0; inode < num_nodes; ++inode) {
	ContactNode<Real> *node = nodes[inode];
        node->ComputeBoundingBoxForSearch(num_configs,
                                          REMAINING_GAP,
                                          NODE_COORD_START,
                                          NODE_COORD_END,
                                          auto_tol,
                                          box_scale,
                                          capture_tol,
                                          node_boxes[inode]);
	node_boxes[inode].set_object_number(node_global_index[iblock][inode]);
      }
      //
      //  Step 2, create the nodal bounding box hierarchy for the current entity
      //
      if(num_nodes > 0) { 
	ObjectBoundingBoxHierarchy::create_hierarchy(node_hierarchy_ptrs[iblock], node_boxes, num_nodes);
	if(iblock+1 <= num_search_entities) node_hierarchy_ptrs[iblock+1] = node_hierarchy_ptrs[iblock] + num_nodes * 2 - 1;
	//
	//  Add the bounding box for the current entity to the list of entity bounding boxes
	//
	ObjectBoundingBox node_block_box(*(node_hierarchy_ptrs[iblock]), iblock);
	node_block_boxes.push_back(node_block_box);
      } else {
	if(iblock+1 <= num_search_entities) node_hierarchy_ptrs[iblock+1] = node_hierarchy_ptrs[iblock];
      }
      if(num_nodes == 0) node_hierarchy_ptrs[iblock] = NULL;
    }
    delete [] node_boxes;
    //
    //  Build a bounding box hierarchy of the nodal entity bounding boxes
    //
    if(node_block_boxes.size() > 0) {
      node_block_BB_hierarchy = new ObjectBoundingBoxHierarchy[node_block_boxes.size() * 2 - 1];
      ObjectBoundingBoxHierarchy::create_hierarchy(node_block_BB_hierarchy, node_block_boxes.get_buffer(), node_block_boxes.size());
    }
  }
}

bool ContactNodeFaceRangeSearch::valid_search() {

  if(search_type == BASIC_SEARCH) {
    if(node_hierarchy == NULL) return false;
  } else if(search_type == KEY_AWARE_SEARCH) {
    if(node_hierarchy == NULL) return false;
  } else if(search_type == ENTITY_SEARCH) {
    if(node_block_BB_hierarchy == NULL) return false;;
  }
  return true;
}

void ContactNodeFaceRangeSearch::search_for_overlap(ContactBoundingBox &face_search_box,
                                                    ACME::Int_Vector &node_keys,
                                                    ACME::ContactNode_Vector &node_list,
                                                    std::vector<bool> &valid_inter) {

  node_list.empty();
  node_keys.empty();
  if(search_type == BASIC_SEARCH) {
    PRECONDITION(node_hierarchy);
    ObjectBoundingBoxHierarchy::search_for_overlap(node_hierarchy,
                                                   face_search_box,
                                                   Nodes,
                                                   node_list);
  } else if(search_type == KEY_AWARE_SEARCH) {
    PRECONDITION(node_hierarchy);
    ObjectBoundingBoxHierarchy::search_for_overlap_key_aware(node_hierarchy,
			      			      	     face_search_box,
							     Nodes,
							     node_keys,
							     node_list,
							     valid_inter,
                                                             node_key_array);
  } else if(search_type == ENTITY_SEARCH) {
    PRECONDITION(node_block_BB_hierarchy);
    ObjectBoundingBoxHierarchy::search_for_overlap_hybrid_tree(node_block_BB_hierarchy,
                                                               node_hierarchy_ptrs,
							       face_search_box,
							       Nodes,
							       node_keys,
							       node_list,
							       valid_inter);
  }
}

void ContactNodeFaceRangeSearch::process_analytic_surfaces(ContactNodeFaceInteraction::InteractionSource Process_Method,
                                                           VariableHandle POSITION, VariableHandle POSITION_2) {
  if(search_type == BASIC_SEARCH) {
    search->Process_Analytic_Surfaces(node_hierarchy, Nodes, Process_Method, POSITION, POSITION_2);
  } else if(search_type == KEY_AWARE_SEARCH) {
    search->Process_Analytic_Surfaces(node_hierarchy, Nodes, Process_Method, POSITION, POSITION_2);
  } else if(search_type == ENTITY_SEARCH) {
    search->Process_Analytic_Surfaces_New(node_hierarchy_ptrs, Nodes, Process_Method, POSITION, POSITION_2);
  }
}


ContactNodeFaceRangeSearch::~ContactNodeFaceRangeSearch() {
  if(node_hierarchy) delete [] node_hierarchy;
  if(node_key_array) delete [] node_key_array;
  if(node_block_BB_hierarchy) delete [] node_block_BB_hierarchy;
  if(node_hierarchy_ptrs) delete [] node_hierarchy_ptrs;
  if(scratch) delete [] scratch;
  if(index) delete [] index;
  if(rank) delete [] rank;
  if(rank2) delete [] rank2;
  if(position) delete [] position;
  if(node_map) delete [] node_map;   
}



namespace ACME {

  //
  //  Determine what type of range search to perform.
  //
  //  NKC Note, need better hueristics to pick search, probally best to just occasinally 
  //  determine search effort and switch to a lower effort method when available.  Best search
  //  type will vary by problem set up, processor, and analysis time.  However, all should
  //  yield identical awnswers so can switch dynamically.
  //
  //    0) Basic Search:  Do a all to all search and later throw out anything that the user
  //    said shouldn't interact.  Use this search when the interaction matrix is full.
  //
  //    1) Key Aware Search:  For this search each object in the node hierarchy is given a
  //    key as well as a bounding box.  Use of keys allow throwing out interactions based
  //    on what the user speficies as potential interations.  This is particularlly useful
  //    for cases where there is little self contact or where there are several pure
  //    master-slave interations.  This search is arbitrarly used when the interaction matrix 
  //    is not totally full.
  //
  //    2) Enity Search:  Construct a object tree for every entity.  Search only trees 
  //    specified as able to interact against each other.  This type of search can
  //    be very benifical in cases without self contact.  This search is arbitrarly used if
  //    the interaction matrix contains no more than 10% self contact
  //
  //    3) Old Search: Use ACME's old list range based search method.  The list range search
  //    is cheaper to setup than the tree based searches.  However, it searching is significantly
  //    slower than any of the tree based searches.  Additionaly the old Search has a O(N^1.5) complexity
  //    compared to the tree searches O(N Log N) thus it may become extremly expensive for large data sets.
  //
  ContactNodeFaceRangeSearch::SearchType determine_range_search_type(ContactSearchData *search_data) {
    const int num_search_entities = search_data->Num_Search_Entities();
    int num_self_interactions = 0;
    int matrix_size = 0;
    int active_interactions = 0;
    for(int irow = 0; irow < num_search_entities; ++irow) {
      for(int icol = 0; icol < num_search_entities; ++icol) {
        if((int)search_data->Get_Search_Data(ContactSearch::INTERACTION_TYPE, irow, icol) != ContactSearch::NO_INTERACTION) {
          active_interactions++;
          if(irow == icol) num_self_interactions++;
	}
        matrix_size++;
      }
    }
    if(num_self_interactions < 0.10 *(Real) num_search_entities) {
      return ContactNodeFaceRangeSearch::ENTITY_SEARCH;
    } else if(active_interactions < 0.25 * matrix_size) {
      return ContactNodeFaceRangeSearch::KEY_AWARE_SEARCH;
    } else {
      return ContactNodeFaceRangeSearch::BASIC_SEARCH;
    }
  }

  //-------------------------------------------------------------------------------------------------------------------------
  //  Overlapping sphere ghost takes as input a number of spheres owned by the current processor.
  //  It returns the spheres on other processors that need to be ghosted.
  //  NOTE, for this routine to be efficent the objects must be in an approximetly RCB like decomposition.  Meaning, all 
  //  the objects for a given processor should be contained in a compact axis alligned box.  If this condition is not met
  //  much more ghosting than nessecary may occur.
  //
  int Overlapping_Sphere_Ghost(const int num_sphere, 
                               const Real* sphere_data,
		               MPI_Comm &mpi_communicator,
			       Int_Vector &ghost_indexes,
                               Int_Vector &ghost_procs) {
    //
    //  Clear the return data arrays
    //
    ghost_indexes.clear();
    ghost_procs.clear();
#ifdef CONTACT_NO_MPI
    return 0;
#else
    //
    //  Determine the total number of processors involved in the communication and the current processor number
    //
    int num_procs(1);
    int current_proc(0);
    MPI_Comm_size(mpi_communicator, &num_procs);
    MPI_Comm_rank(mpi_communicator, &current_proc);
    if(num_procs == 1) return 0;
    //
    //  Compute the processor local bounding box.  Store the box in a unique entry in a global processor bounding box array.
    //
    ObjectBoundingBox current_processor_box;
    for(int isphere = 0; isphere < num_sphere; ++isphere) {
      const Real *sphere_coords = sphere_data + isphere * 4;
      const Real *sphere_rad    = sphere_data + isphere * 4 + 3;
      current_processor_box.add_sphere(sphere_coords, *sphere_rad);

    }
    ObjectBoundingBox *proc_box_array = new ObjectBoundingBox[num_procs];
    proc_box_array[current_proc] = current_processor_box;

    //
    //  Do a global communication to communicate all processor bounding boxes to all processors in the group
    //
    ObjectBoundingBox::global_box_combine(proc_box_array, num_procs, mpi_communicator);
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      proc_box_array[iproc].set_object_number(iproc);
    }
    //
    //  Create a hierarchy of processor bounding boxes.  This hierarchy will be used to search for overlaps between processors and
    //  objects.
    //
    int hierarchy_size = 2 * num_procs - 1;
    ObjectBoundingBoxHierarchy *proc_box_hierarchy = NULL;
    if(hierarchy_size > 0) {
      proc_box_hierarchy = new ObjectBoundingBoxHierarchy[hierarchy_size];
      ObjectBoundingBoxHierarchy::create_hierarchy(proc_box_hierarchy, proc_box_array, num_procs);
    }    
    delete [] proc_box_array;
    if(proc_box_hierarchy == NULL) return 0;
    //
    //  Determine potential communication partners.  If this processor bounding box overlaps another processor's bounding box
    //  Then that processor is a potential communication partner
    //
    int *proc_list = new int[num_procs];
    Int_Vector comm_partners(num_procs, 0);
    int num_proc_overlaps = 0;
    ObjectBoundingBoxHierarchy::search_for_overlap_loop(proc_box_hierarchy,current_processor_box, proc_list, num_proc_overlaps);
    for(int i = 0; i < num_proc_overlaps; ++i) {
      if(proc_list[i] == current_proc) continue;
      comm_partners[proc_list[i]] = 1;
    }
    //
    //  Loop over all objects and compute potential objects to ghosts.  Add potential objects to the sending lists
    //
    Real_Vector *send_spheres = new Real_Vector[num_procs];
    Real_Vector *recv_spheres = new Real_Vector[num_procs];
    Int_Vector *send_indexes = new Int_Vector[num_procs];
    for(int isphere = 0; isphere < num_sphere; ++isphere) {
      ContactBoundingBox object_box;
      const Real *sphere_coords = sphere_data + isphere * 4;
      const Real *sphere_rad    = sphere_data + isphere * 4 + 3; 
      object_box.add_sphere(sphere_coords, *sphere_rad);
      int num_overlaps = 0;
      ObjectBoundingBoxHierarchy::search_for_overlap_loop(proc_box_hierarchy,object_box, proc_list, num_overlaps);
      for(int i = 0; i < num_overlaps; ++i) { 
        int overlapping_proc = proc_list[i];
        if(overlapping_proc == current_proc) continue;
        send_indexes[overlapping_proc].push_back(isphere);
        send_spheres[overlapping_proc].push_back(sphere_coords[0]);
        send_spheres[overlapping_proc].push_back(sphere_coords[1]);
        send_spheres[overlapping_proc].push_back(sphere_coords[2]);
        send_spheres[overlapping_proc].push_back(*sphere_rad);



      }
    }
    //
    //  Exchange sphere data
    //
    Parallel_Data_Exchange(send_spheres, recv_spheres, comm_partners, mpi_communicator);
    //
    //  On each processor, search each sphere in the send list against each sphere in the receive list.  
    //  The spheres in the send list that overlap spheres in the receive list must be ghosted.
    //
    Int_Vector interaction_list;
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(comm_partners[iproc]) {
        Int_Vector first_interaction;
        Int_Vector last_interaction;
        int num_send_spheres = send_spheres[iproc].size()/4;
        int num_recv_spheres = recv_spheres[iproc].size()/4;
        Overlapping_Sphere_Search(num_send_spheres,
                                  send_spheres[iproc].get_buffer(),
                                  num_recv_spheres,
                                  recv_spheres[iproc].get_buffer(),
                                  interaction_list,
                                  first_interaction,
                                  last_interaction);

	for(int isphere = 0; isphere < num_send_spheres; ++isphere) {
          if(first_interaction[isphere] < last_interaction[isphere]) {
            ghost_indexes.push_back(send_indexes[iproc][isphere]);
            ghost_procs.push_back(iproc);
	  }
	}
      }
    }
    delete [] send_indexes;
    delete [] proc_list;
    delete [] recv_spheres;
    delete [] send_spheres;
    delete [] proc_box_hierarchy;

    return 0;
#endif

  }
  //-----------------------------------------------------------------------------------------------------------------
  //
  //  Overlapping sphere search takes as input a number of spheres defined by center points and radii.  
  //  It finds all interactions such that spheres in group1 interact with spheres in group2.  This routine 
  //  assumes exactly three dimensional spheres.
  //
  int Overlapping_Sphere_Search(const int num_sphere1, 
                                const Real* sphere_data1,
                                const int num_sphere2, 
                                const Real* sphere_data2,
			        Int_Vector& interaction_list,
			        Int_Vector& first_interaction,
                                Int_Vector& last_interaction) {
    //
    //  Set the known return vector sizes
    //
    first_interaction.resize(num_sphere1, -1);
    last_interaction.resize(num_sphere1, -1);
    interaction_list.empty();
    //
    //  Build bounding boxes around all the spheres in group 2.  Use the Integerized bounding boxes so at the same
    //  time calculate a box that encompases the entire domain.
    //
    DomainBox domain_box;
    ContactBoundingBox *sphere_boxes = new ContactBoundingBox[num_sphere2];
    for(int isphere = 0; isphere < num_sphere2; ++isphere) {
      const Real *sphere_coords = sphere_data2 + isphere*4;
      const Real *sphere_rad    = sphere_data2 + isphere*4 + 3;
      sphere_boxes[isphere].add_sphere(sphere_coords, *sphere_rad);
      domain_box.add_box(sphere_boxes[isphere]);
    }
    //
    //  Transform the Real bounding boxes into integer bounding boxes
    //
    domain_box.calculate_constants();
    ObjectBoundingBox_Int *sphere_boxes_int = new ObjectBoundingBox_Int[num_sphere2];
    for(int isphere = 0; isphere < num_sphere2; ++isphere) {
      sphere_boxes_int[isphere].set_box(sphere_boxes[isphere], domain_box);
      sphere_boxes_int[isphere].set_object_number(isphere);
    }
    delete [] sphere_boxes;
    //
    //  Compute the bounding box hierarchy.  This is a hierachical tree representation off all the
    //  object bounding boxes.
    //
    int hierarchy_size = 2 * num_sphere2 - 1;
    ObjectBoundingBoxHierarchy_Int *sphere_hierarchy = NULL;
    if(hierarchy_size > 0) {
      sphere_hierarchy = new ObjectBoundingBoxHierarchy_Int[hierarchy_size];
      ObjectBoundingBoxHierarchy_Int::create_hierarchy(sphere_hierarchy, sphere_boxes_int, num_sphere2);    
    }    
    delete [] sphere_boxes_int;
    if(sphere_hierarchy == NULL) return 0;
    //
    //  Create an array to store interactions returned by the recursive search routines.  There are at maximum
    //  N interactions per object when searching N objects
    //
    int *overlap_list = new int[num_sphere2];
    //
    //  Loop over all spheres in group1 and search them against those objects in group2
    //
    for(int isphere = 0; isphere < num_sphere1; ++isphere) {
      const Real *sphere_coords = sphere_data1 + isphere*4;
      const Real *sphere_rad    = sphere_data1 + isphere*4 + 3;
      ContactBoundingBox search_sphere;
      search_sphere.add_sphere(sphere_coords, *sphere_rad);
      ContactBoundingBox_Int search_sphere_int;
      search_sphere_int.set_box(search_sphere, domain_box);
      //
      //  Extract the index numbers of all spheres that interact with sphere isphere
      //
      int list_size = 0;
      ObjectBoundingBoxHierarchy_Int::search_for_overlap_loop(sphere_hierarchy,search_sphere_int, 
                                                   overlap_list,
                                                   list_size);
      //
      //  Record the value of the first interaction for the current sphere
      //
      first_interaction[isphere] = interaction_list.size();
      //
      //  Add the current interactions to the interaction list, do the acutall sphere overlap calculation
      //  here.  (The overlap_loop call calcualtes overlapping bounding boxes, not true spheres.)
      //
      const Real *idata = sphere_data1 + isphere * 4;
      for(int ilist = 0; ilist < list_size; ++ilist) {
        int jsphere         = overlap_list[ilist];
        const Real *jdata   = sphere_data2 + jsphere * 4;
        if(Sphere_Sphere_Overlap(idata, jdata)){
          interaction_list.push_back(jsphere);
	}
      }
      //
      //  Record the value of the last interaction for the current sphere.
      //
      last_interaction[isphere] = interaction_list.size();
    }

    delete [] overlap_list;
    delete [] sphere_hierarchy;

    return 0; 
  }


  //-------------------------------------------------------------------------------------------------------------------------
  //  Sphere_Point_Ghost takes as input a number of spheres and points owned by the current processor.
  //  It returns the points that need to be ghosted to other processors.
  //  NOTE, for this routine to be efficent the objects must be in an approximetly RCB like decomposition.  Meaning, all 
  //  the objects for a given processor should be contained in a compact axis alligned box.  If this condition is not met
  //  much more ghosting than nessecary may occur.
  //
  int Sphere_Point_Ghost(const int num_sphere, 
                         const Real* sphere_data,
                         const int num_point,
                         const Real* point_data,
	                 MPI_Comm &mpi_communicator,
			 Int_Vector &ghost_indexes,
                         Int_Vector &ghost_procs) {
    //
    //  Clear the return data arrays
    //
    ghost_indexes.clear();
    ghost_procs.clear();
#ifdef CONTACT_NO_MPI
    return 0;
#else
    //
    //  Determine the total number of processors involved in the communication and the current processor number
    //
    int num_procs(1);
    int current_proc(0);
    MPI_Comm_size(mpi_communicator, &num_procs);
    MPI_Comm_rank(mpi_communicator, &current_proc);
    if(num_procs == 1) return 0;
    //
    //  Compute the processor local bounding boxes for spheres and for points.  
    //  Store the boxes in unique entries in a global processor bounding box array.
    //
    ObjectBoundingBox sphere_processor_box;
    ObjectBoundingBox point_processor_box;
    for(int isphere = 0; isphere < num_sphere; ++isphere) {
      const Real *sphere_coords = sphere_data + isphere * 4;
      const Real *sphere_rad    = sphere_data + isphere * 4 + 3;
      sphere_processor_box.add_sphere(sphere_coords, *sphere_rad);
    }
    for(int ipoint = 0; ipoint < num_point; ++ipoint) {
      const Real *point_coords = point_data + ipoint * 3;
      point_processor_box.add_point(point_coords);
    }

    ObjectBoundingBox *sphere_proc_box_array = new ObjectBoundingBox[num_procs];
    ObjectBoundingBox *point_proc_box_array  = new ObjectBoundingBox[num_procs];

    sphere_proc_box_array[current_proc] = sphere_processor_box;
    point_proc_box_array[current_proc]  = point_processor_box;
    //
    //  Do a global communication to communicate all processor bounding boxes to all processors in the group
    //
    ObjectBoundingBox::global_box_combine(sphere_proc_box_array, num_procs, mpi_communicator);
    ObjectBoundingBox::global_box_combine(point_proc_box_array,  num_procs, mpi_communicator);

    for(int iproc = 0; iproc < num_procs; ++iproc) {
      sphere_proc_box_array[iproc].set_object_number(iproc);
      point_proc_box_array[iproc].set_object_number(iproc);
    }
    //
    //  Create a hierarchies of sphere and point processor bounding boxes.  
    //  This hierarchy will be used to search for overlaps between processors and
    //  objects.
    //
    int hierarchy_size = 2 * num_procs - 1;
    ObjectBoundingBoxHierarchy *sphere_box_hierarchy = NULL;
    ObjectBoundingBoxHierarchy *point_box_hierarchy = NULL;

    if(hierarchy_size > 0) {
      sphere_box_hierarchy = new ObjectBoundingBoxHierarchy[hierarchy_size];
      point_box_hierarchy = new ObjectBoundingBoxHierarchy[hierarchy_size];
      ObjectBoundingBoxHierarchy::create_hierarchy(sphere_box_hierarchy, sphere_proc_box_array, num_procs);
      ObjectBoundingBoxHierarchy::create_hierarchy(point_box_hierarchy, point_proc_box_array, num_procs);
    }    
    delete [] sphere_proc_box_array;
    delete [] point_proc_box_array;
    //
    //  Determine potential communication partners.  If this processores point box overlaps another processor's sphere box
    //  Then that processor is a potential communication partner
    //
    int *proc_list = new int[num_procs];
    Int_Vector comm_partners(num_procs, 0);
    int num_overlaps = 0;
    ObjectBoundingBoxHierarchy::search_for_overlap_loop(sphere_box_hierarchy,point_processor_box, proc_list, num_overlaps);
    for(int i = 0; i < num_overlaps; ++i) {
      if(proc_list[i] == current_proc) continue;
      comm_partners[proc_list[i]] = 1;
    }
    //
    //  Loop over all objects and compute potential objects to ghosts.  Add potential objects to the sending lists
    //
    Real_Vector *send_points = new Real_Vector[num_procs];
    Int_Vector *send_point_indexes = new Int_Vector[num_procs];
    for(int ipoint = 0; ipoint < num_point; ++ipoint) {
      ContactBoundingBox object_box;
      const Real *point_coords = point_data + ipoint * 3;
      object_box.add_point(point_coords);
      num_overlaps = 0;
      ObjectBoundingBoxHierarchy::search_for_overlap_loop(sphere_box_hierarchy,object_box, proc_list, num_overlaps);
      for(int i = 0; i < num_overlaps; ++i) { 
        int overlapping_proc = proc_list[i];
        if(overlapping_proc == current_proc) continue;
        send_points[overlapping_proc].push_back(point_coords[0]);
        send_points[overlapping_proc].push_back(point_coords[1]);
        send_points[overlapping_proc].push_back(point_coords[2]);
      }
    }

    Real_Vector *send_spheres = new Real_Vector[num_procs];
    Real_Vector *recv_spheres = new Real_Vector[num_procs];
    for(int isphere = 0; isphere < num_sphere; ++isphere) {
      ContactBoundingBox object_box;
      const Real *sphere_coords = sphere_data + isphere * 4;
      const Real *sphere_rad    = sphere_data + isphere * 4 + 3;
      object_box.add_sphere(sphere_coords, *sphere_rad);
      num_overlaps = 0;
      ObjectBoundingBoxHierarchy::search_for_overlap_loop(point_box_hierarchy,object_box, proc_list, num_overlaps);
      for(int i = 0; i < num_overlaps; ++i) { 
        int overlapping_proc = proc_list[i];
        if(overlapping_proc == current_proc) continue;
        send_spheres[overlapping_proc].push_back(sphere_coords[0]);
        send_spheres[overlapping_proc].push_back(sphere_coords[1]);
        send_spheres[overlapping_proc].push_back(sphere_coords[2]);
        send_spheres[overlapping_proc].push_back(*sphere_rad);
      }
    }
    //
    //  Exchange sphere data
    //
    Parallel_Data_Exchange(send_spheres, recv_spheres, comm_partners, mpi_communicator);
    //
    //  On each processor, search each point in the send list against each sphere in the receive list.  
    //  The points in the send list that overlap spheres in the receive list must be ghosted.
    //
    Int_Vector interaction_list;
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(comm_partners[iproc]) {
        Int_Vector first_interaction;
        Int_Vector last_interaction;
        int num_send_points  = send_points [iproc].size()/3;
        int num_recv_spheres = recv_spheres[iproc].size()/4;
        Sphere_Point_Search(num_recv_spheres,
                            recv_spheres[iproc].get_buffer(),
                            num_send_points,
                            send_points[iproc].get_buffer(),
                            interaction_list,
                            first_interaction,
                            last_interaction);

        Int_Vector send_point_marker(num_send_points, 0);
	for(int isphere = 0; isphere < num_recv_spheres; ++isphere) {
          for(int ipoint = first_interaction[isphere]; ipoint < last_interaction[isphere]; ++ipoint) {
            int point_num = interaction_list[ipoint];
            send_point_marker[point_num] = 1;
	  }
	}
        for(int ipoint = 0; ipoint < num_send_points; ++ipoint) {
          if(send_point_marker[ipoint] == 1) {
            ghost_indexes.push_back(send_point_indexes[iproc][ipoint]);
            ghost_procs.push_back(iproc);
	  }
	}
      }
    }
    delete [] send_point_indexes;
    delete [] proc_list;
    delete [] send_points;
    delete [] recv_spheres;
    delete [] send_spheres;
    delete [] sphere_box_hierarchy;
    delete [] point_box_hierarchy;

    return 0;
#endif
  }

  //-----------------------------------------------------------------------------------------------------------------
  //
  //
  //  Sphere_Point_Search takes as input a number of spheres defined by center points and radii.  
  //  And a series of points.  The search finds intersections between the points and the spheres
  //
  int Sphere_Point_Search(const int num_sphere, 
                          const Real* sphere_data,
                          const int num_point, 
                          const Real* point_data,
			  Int_Vector& interaction_list,
		          Int_Vector& first_interaction,
                          Int_Vector& last_interaction) {
    //
    //  Set the known return vector sizes
    //
    first_interaction.resize(num_sphere, -1);
    last_interaction.resize(num_sphere, -1);
    interaction_list.empty();
    //
    //  Build bounding boxes around all the points.  Use the Integerized bounding boxes so at the same
    //  time calculate a box that encompases the entire domain.
    //
    DomainBox domain_box;
    ContactBoundingBox *point_boxes = new ContactBoundingBox[num_point];
    for(int ipoint = 0; ipoint < num_point; ++ipoint) {
      const Real *point_coords = point_data + ipoint*3;
      point_boxes[ipoint].add_point(point_coords);
      domain_box.add_box(point_boxes[ipoint]);
    }
    //
    //  Transform the Real bounding boxes into integer bounding boxes
    //
    domain_box.calculate_constants();
    ObjectBoundingBox_Int *point_boxes_int = new ObjectBoundingBox_Int[num_point];
    for(int ipoint = 0; ipoint < num_point; ++ipoint) {
      point_boxes_int[ipoint].set_box(point_boxes[ipoint], domain_box);
      point_boxes_int[ipoint].set_object_number(ipoint);
    }
    delete [] point_boxes;
    //
    //  Compute the bounding box hierarchy.  This is a hierachical tree representation off all the
    //  object bounding boxes.
    //
    int hierarchy_size = 2 * num_point - 1;
    ObjectBoundingBoxHierarchy_Int *point_hierarchy = NULL;
    if(hierarchy_size > 0) {
      point_hierarchy = new ObjectBoundingBoxHierarchy_Int[hierarchy_size];
      ObjectBoundingBoxHierarchy_Int::create_hierarchy(point_hierarchy, point_boxes_int, num_point);    
    }    
    delete [] point_boxes_int;
    if(point_hierarchy == NULL) return 0;
    //
    //  Create an array to store interactions returned by the recursive search routines.  There are at maximum
    //  N interactions per object when searching N objects
    //
    int *overlap_list = new int[num_point];
    //
    //  Loop over all spheres in group1 and search them against those objects in group2
    //
    for(int isphere = 0; isphere < num_sphere; ++isphere) {
      const Real *sphere_coords = sphere_data + isphere*4;
      const Real *sphere_rad    = sphere_data + isphere*4 + 3;
      ContactBoundingBox search_sphere;
      search_sphere.add_sphere(sphere_coords, *sphere_rad);
      ContactBoundingBox_Int search_sphere_int;
      search_sphere_int.set_box(search_sphere, domain_box);
      //
      //  Extract the index numbers of all spheres that interact with sphere isphere
      //
      int list_size = 0;
      ObjectBoundingBoxHierarchy_Int::search_for_overlap_loop(point_hierarchy,search_sphere_int, 
                                                  overlap_list,
                                                  list_size);
      //
      //  Record the value of the first interaction for the current sphere
      //
      first_interaction[isphere] = interaction_list.size();
      //
      //  Add the current interactions to the interaction list, do the acutall sphere overlap calculation
      //  here.  (The overlap_loop call calcualtes overlapping bounding boxes, not true spheres.)
      //
      const Real *isphere_data = sphere_data + isphere * 4;
      for(int ilist = 0; ilist < list_size; ++ilist) {
        int jpoint         = overlap_list[ilist];
        const Real *jpoint_data   = point_data + jpoint * 3;
        if(Sphere_Point_Overlap(isphere_data, jpoint_data)){
          interaction_list.push_back(jpoint);
	}
      }
      //
      //  Record the value of the last interaction for the current sphere.
      //
      last_interaction[isphere] = interaction_list.size();
    }

    delete [] overlap_list;
    delete [] point_hierarchy;

    return 0; 
  }

  //--------------------------------------------------------------------------------------------------


  //-------------------------------------------------------------------------------------------------------------------------
  //  BoxA_BoxB_Ghost takes as input a number of box set A and box set B owned by the current processor.
  //  It returns the boxes in set B that need to be ghosted to other processors.
  //  NOTE, for this routine to be efficent the objects must be in an approximetly RCB like decomposition.  Meaning, all 
  //  the objects for a given processor should be contained in a compact axis alligned box.  If this condition is not met
  //  much more ghosting than nessecary may occur.
  //
  int BoxA_BoxB_Ghost(const int num_boxA, 
                      const Real* boxA_data,
                      const int num_boxB,
                      const Real* boxB_data,
	              MPI_Comm &mpi_communicator,
		      Int_Vector &ghost_indexes,
                      Int_Vector &ghost_procs) {
    //
    //  Clear the return data arrays
    //
    ghost_indexes.clear();
    ghost_procs.clear();
#ifdef CONTACT_NO_MPI
    return 0;
#else
    //
    //  Determine the total number of processors involved in the communication and the current processor number
    //
    int num_procs(1);
    int current_proc(0);
    MPI_Comm_size(mpi_communicator, &num_procs);
    MPI_Comm_rank(mpi_communicator, &current_proc);
    if(num_procs == 1) return 0;
    //
    //  Compute the processor local bounding boxes for the box sets
    //  Store the boxes in unique entries in a global processor bounding box array.
    //
    ObjectBoundingBox boxA_processor_box;
    ObjectBoundingBox boxB_processor_box;
    for(int iboxA = 0; iboxA < num_boxA; ++iboxA) {
      const Real *corner1 = boxA_data + iboxA * 6 + 0;
      const Real *corner2 = boxA_data + iboxA * 6 + 3;
      boxA_processor_box.add_point(corner1);
      boxA_processor_box.add_point(corner2);
    }
    for(int iboxB = 0; iboxB < num_boxB; ++iboxB) {
      const Real *corner1 = boxB_data + iboxB * 6 + 0;
      const Real *corner2 = boxB_data + iboxB * 6 + 3;
      boxB_processor_box.add_point(corner1);
      boxB_processor_box.add_point(corner2);
    }

    ObjectBoundingBox *boxA_proc_box_array = new ObjectBoundingBox[num_procs];
    ObjectBoundingBox *boxB_proc_box_array = new ObjectBoundingBox[num_procs];

    boxA_proc_box_array[current_proc] = boxA_processor_box;
    boxB_proc_box_array[current_proc] = boxB_processor_box;
    //
    //  Do a global communication to communicate all processor bounding boxes to all processors in the group
    //
    ObjectBoundingBox::global_box_combine(boxA_proc_box_array, num_procs, mpi_communicator);
    ObjectBoundingBox::global_box_combine(boxB_proc_box_array, num_procs, mpi_communicator);

    for(int iproc = 0; iproc < num_procs; ++iproc) {
      boxA_proc_box_array[iproc].set_object_number(iproc);
      boxB_proc_box_array[iproc].set_object_number(iproc);
    }
    //
    //  Create a hierarchies of boxA and boxB processor bounding boxes.  
    //  This hierarchy will be used to search for overlaps between processors and
    //  objects.
    //
    int hierarchy_size = 2 * num_procs - 1;
    ObjectBoundingBoxHierarchy *boxA_box_hierarchy = NULL;
    ObjectBoundingBoxHierarchy *boxB_box_hierarchy = NULL;

    if(hierarchy_size > 0) {
      boxA_box_hierarchy = new ObjectBoundingBoxHierarchy[hierarchy_size];
      boxB_box_hierarchy = new ObjectBoundingBoxHierarchy[hierarchy_size];
      ObjectBoundingBoxHierarchy::create_hierarchy(boxA_box_hierarchy, boxA_proc_box_array, num_procs);
      ObjectBoundingBoxHierarchy::create_hierarchy(boxB_box_hierarchy, boxB_proc_box_array, num_procs);
    }    
    delete [] boxA_proc_box_array;
    delete [] boxB_proc_box_array;
    //
    //  Determine potential communication partners.  If this processores boxB box overlaps another processor's boxA box
    //  Then that processor is a potential communication partner
    //
    //  NKC, may want to break out send and recive for a little more optimization, but tricky....
    //
    int *proc_list = new int[num_procs];
    Int_Vector comm_partners(num_procs, 0);
    int num_overlaps = 0;

    ObjectBoundingBoxHierarchy::search_for_overlap_loop(boxA_box_hierarchy, boxB_processor_box, proc_list, num_overlaps);
    for(int i = 0; i < num_overlaps; ++i) {
      if(proc_list[i] == current_proc) continue;
      comm_partners[proc_list[i]] = 1;
    }

    num_overlaps = 0;
    ObjectBoundingBoxHierarchy::search_for_overlap_loop(boxB_box_hierarchy, boxA_processor_box, proc_list, num_overlaps);
    for(int i = 0; i < num_overlaps; ++i) {
      if(proc_list[i] == current_proc) continue;
      comm_partners[proc_list[i]] = 1;
    }
    //
    //  Loop over all objects and compute potential objects to ghosts.  Add potential objects to the sending lists
    //
    Real_Vector *send_boxBs = new Real_Vector[num_procs];
    Int_Vector  *send_boxB_indexes = new Int_Vector[num_procs];
    for(int iboxB = 0; iboxB < num_boxB; ++iboxB) {
      ContactBoundingBox object_box;
      const Real *boxB_corner1 = boxB_data + iboxB * 6 + 0;
      const Real *boxB_corner2 = boxB_data + iboxB * 6 + 3;
      object_box.add_point(boxB_corner1);
      object_box.add_point(boxB_corner2);
      num_overlaps = 0;
      ObjectBoundingBoxHierarchy::search_for_overlap_loop(boxA_box_hierarchy,object_box, proc_list, num_overlaps);
      for(int i = 0; i < num_overlaps; ++i) { 
        int overlapping_proc = proc_list[i];
        if(overlapping_proc == current_proc) continue;
        send_boxBs[overlapping_proc].push_back(boxB_corner1[0]);
        send_boxBs[overlapping_proc].push_back(boxB_corner1[1]);
        send_boxBs[overlapping_proc].push_back(boxB_corner1[2]);
        send_boxBs[overlapping_proc].push_back(boxB_corner2[0]);
        send_boxBs[overlapping_proc].push_back(boxB_corner2[1]);
        send_boxBs[overlapping_proc].push_back(boxB_corner2[2]);
        send_boxB_indexes[overlapping_proc].push_back(iboxB);
      }
    }

    Real_Vector *send_boxAs = new Real_Vector[num_procs];
    Real_Vector *recv_boxAs = new Real_Vector[num_procs];
    for(int iboxA = 0; iboxA < num_boxA; ++iboxA) {
      ContactBoundingBox object_box;
      const Real *boxA_corner1 = boxA_data + iboxA * 6 + 0;
      const Real *boxA_corner2 = boxA_data + iboxA * 6 + 3;
      object_box.add_point(boxA_corner1);
      object_box.add_point(boxA_corner2);
      num_overlaps = 0;
      ObjectBoundingBoxHierarchy::search_for_overlap_loop(boxB_box_hierarchy,object_box, proc_list, num_overlaps);
      for(int i = 0; i < num_overlaps; ++i) { 
        int overlapping_proc = proc_list[i];
        if(overlapping_proc == current_proc) continue;
        send_boxAs[overlapping_proc].push_back(boxA_corner1[0]);
        send_boxAs[overlapping_proc].push_back(boxA_corner1[1]);
        send_boxAs[overlapping_proc].push_back(boxA_corner1[2]);
        send_boxAs[overlapping_proc].push_back(boxA_corner2[0]);
        send_boxAs[overlapping_proc].push_back(boxA_corner2[1]);
        send_boxAs[overlapping_proc].push_back(boxA_corner2[2]);
      }
    }
    //
    //  Exchange boxA data
    //
    Parallel_Data_Exchange(send_boxAs, recv_boxAs, comm_partners, mpi_communicator);
    //
    //  On each processor, search each boxB in the send list against each boxA in the receive list.  
    //  The boxBs in the send list that overlap boxAs in the receive list must be ghosted.
    //
    Int_Vector interaction_list;
    for(int iproc = 0; iproc < num_procs; ++iproc) {
      if(comm_partners[iproc]) {
        Int_Vector first_interaction;
        Int_Vector last_interaction;
        int num_send_boxBs = send_boxBs[iproc].size()/6;
        int num_recv_boxAs = recv_boxAs[iproc].size()/6;
        BoxA_BoxB_Search(num_recv_boxAs,
                         recv_boxAs[iproc].get_buffer(),
                         num_send_boxBs,
                         send_boxBs[iproc].get_buffer(),
                         interaction_list,
                         first_interaction,
                         last_interaction);

        Int_Vector send_boxBs_marker(num_send_boxBs, 0);
	for(int iboxA = 0; iboxA < num_recv_boxAs; ++iboxA) {
          for(int iboxB = first_interaction[iboxA]; iboxB < last_interaction[iboxA]; ++iboxB) {
            int boxB_num = interaction_list[iboxB];
            send_boxBs_marker[boxB_num] = 1;
	  }
	}
        for(int iboxB = 0; iboxB < num_send_boxBs; ++iboxB) {
          if(send_boxBs_marker[iboxB] == 1) {
            ghost_indexes.push_back(send_boxB_indexes[iproc][iboxB]);
            ghost_procs.push_back(iproc);
	  }
	}
      }
    }
    delete [] send_boxB_indexes;
    delete [] proc_list;
    delete [] send_boxBs;
    delete [] recv_boxAs;
    delete [] send_boxAs;
    delete [] boxA_box_hierarchy;
    delete [] boxB_box_hierarchy;

    return 0;
#endif
  }

  //-----------------------------------------------------------------------------------------------------------------
  //
  //
  //  BoxA_BoxB_Search takes as input a number of boxAs and boxBs.  
  //  The search finds intersections between the boxes in set B and the boxes in set A.
  //
  int BoxA_BoxB_Search(const int num_boxA, 
                       const Real* boxA_data,
                       const int num_boxB, 
                       const Real* boxB_data,
	               Int_Vector& interaction_list,
		       Int_Vector& first_interaction,
                       Int_Vector& last_interaction) {
    //
    //  Set the known return vector sizes
    //
    first_interaction.resize(num_boxA, -1);
    last_interaction.resize(num_boxA, -1);
    interaction_list.empty();
    //
    //  Build bounding boxes around all the boxBs.  Use the Integerized bounding boxes so at the same
    //  time calculate a box that encompases the entire domain.
    //
    DomainBox domain_box;
    ContactBoundingBox *boxB_boxes = new ContactBoundingBox[num_boxB];
    for(int iboxB = 0; iboxB < num_boxB; ++iboxB) {
      const Real *boxB_corner1 = boxB_data + iboxB*6 + 0;
      const Real *boxB_corner2 = boxB_data + iboxB*6 + 3;
      boxB_boxes[iboxB].add_point(boxB_corner1);
      boxB_boxes[iboxB].add_point(boxB_corner2);
      domain_box.add_box(boxB_boxes[iboxB]);
    }
    //
    //  Transform the Real bounding boxes into integer bounding boxes
    //
    domain_box.calculate_constants();
    ObjectBoundingBox_Int *boxB_boxes_int = new ObjectBoundingBox_Int[num_boxB];
    for(int iboxB = 0; iboxB < num_boxB; ++iboxB) {
      boxB_boxes_int[iboxB].set_box(boxB_boxes[iboxB], domain_box);
      boxB_boxes_int[iboxB].set_object_number(iboxB);
    }
    delete [] boxB_boxes;
    //
    //  Compute the bounding box hierarchy.  This is a hierachical tree representation off all the
    //  object bounding boxes.
    //
    int hierarchy_size = 2 * num_boxB - 1;
    ObjectBoundingBoxHierarchy_Int *boxB_hierarchy = NULL;
    if(hierarchy_size > 0) {
      boxB_hierarchy = new ObjectBoundingBoxHierarchy_Int[hierarchy_size];
      ObjectBoundingBoxHierarchy_Int::create_hierarchy(boxB_hierarchy, boxB_boxes_int, num_boxB);    
    }    
    delete [] boxB_boxes_int;
    if(boxB_hierarchy == NULL) return 0;
    //
    //  Create an array to store interactions returned by the recursive search routines.  There are at maximum
    //  N interactions per object when searching N objects
    //
    int *overlap_list = new int[num_boxB];
    //
    //  Loop over all boxAs in group1 and search them against those objects in group2
    //
    for(int iboxA = 0; iboxA < num_boxA; ++iboxA) {
      const Real *boxA_corner1 = boxA_data + iboxA*6 + 0;
      const Real *boxA_corner2 = boxA_data + iboxA*6 + 3;
      ContactBoundingBox search_boxA;
      search_boxA.add_point(boxA_corner1);
      search_boxA.add_point(boxA_corner2);
      ContactBoundingBox_Int search_boxA_int;
      search_boxA_int.set_box(search_boxA, domain_box);
      //
      //  Extract the index numbers of all boxAs that interact with boxA iboxA
      //
      int list_size = 0;
      ObjectBoundingBoxHierarchy_Int::search_for_overlap_loop(boxB_hierarchy,search_boxA_int, 
                                                  overlap_list,
                                                  list_size);
      //
      //  Record the value of the first interaction for the current boxA
      //
      first_interaction[iboxA] = interaction_list.size();
      //
      //  Add the current interactions to the interaction list, do the acutall boxA overlap calculation
      //  here.  (The overlap_loop call calcualtes overlapping bounding boxes, not true boxAs.)
      //
      for(int ilist = 0; ilist < list_size; ++ilist) {
        int jboxB         = overlap_list[ilist];
        interaction_list.push_back(jboxB);
      }
      //
      //  Record the value of the last interaction for the current boxA.
      //
      last_interaction[iboxA] = interaction_list.size();
    }

    delete [] overlap_list;
    delete [] boxB_hierarchy;

    return 0; 
  }

  //--------------------------------------------------------------------------------------------------

  //
  //  Helper routine, determine if two spheres overlap.
  //
  bool Sphere_Sphere_Overlap(const Real* isphere_data,
                             const Real* jsphere_data) {
    const Real *icoords = isphere_data;
    const Real *irad    = isphere_data + 3;

    const Real *jcoords = jsphere_data;
    const Real *jrad    = jsphere_data + 3;

    Real dist_squared = (icoords[0] - jcoords[0])*(icoords[0] - jcoords[0]) +
                        (icoords[1] - jcoords[1])*(icoords[1] - jcoords[1]) +
	                (icoords[2] - jcoords[2])*(icoords[2] - jcoords[2]);
    Real rad_squared  = (*irad + *jrad) * (*irad + *jrad);

    if(dist_squared < rad_squared) {
      return true;
    } else {
      return false;
    } 
  }

  bool Sphere_Point_Overlap(const Real* sphere_data,
                            const Real* point_data) {
    const Real *icoords = sphere_data;
    const Real *irad    = sphere_data + 3;

    const Real *jcoords = point_data;

    Real dist_squared = (icoords[0] - jcoords[0])*(icoords[0] - jcoords[0]) +
                        (icoords[1] - jcoords[1])*(icoords[1] - jcoords[1]) +
	                (icoords[2] - jcoords[2])*(icoords[2] - jcoords[2]);
    Real rad_squared  = (*irad) * (*irad);

    if(dist_squared < rad_squared) {
      return true;
    } else {
      return false;
    } 
  }
} // end namespace ACME
