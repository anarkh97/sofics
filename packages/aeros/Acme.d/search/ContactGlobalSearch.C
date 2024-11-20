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


#include "allocators.h"
#include "contact_assert.h"
#include "ContactErrors.h"
#include "ContactFaceBlock.h"
#include "ContactNodeNodeInteraction.h"
#include "ContactNodeFaceInteraction.h"
#include "ContactNodeSurfaceInteraction.h"
#include "ContactFaceFaceInteraction.h"
#include "ContactElementElementInteraction.h"
#include "ContactElement.h"
#include "ContactHexElementL8.h"
#include "ContactWedgeElementL6.h"
#include "ContactShellNode.h"
#include "ContactSearch.h"
#include "ContactTopology.h"
#include "contact_sorting.h"
#include "CString.h"
#include "ContactSearchData.h"
#include "search_methods.h"
#include "ContactParOStream.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactSequentialAllocator.h"
#include "ContactBoundingBox.h"
#include "ContactBoundingBoxHierarchy.h"
#include "ContactBoundingBoxHierarchy_Int.h"
#include "ContactRangeSearch.h"
#include "contact_tolerances.h"
#if (MAX_FFI_DERIVATIVES > 0)
#include "ContactLineEdgeL2.h"
#include "ContactQuadFaceL4.h"
#include "ContactTriFaceL3.h"
#include "ContactShellQuadFaceL4.h"
#include "ContactShellTriFaceL3.h"
#include "ContactHexElementL8.h"
#include "ContactWedgeElementL6.h"
#endif

Real intersection_volume(Real thex1[24][4][3],Real thex2[24][4][3]);


#ifndef CONTACT_NO_MPI
#include "mpi.h"
#include "ContactZoltan.h"
#endif

#include <iostream>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>

using namespace std;

extern "C" {
  double intersect_3d(double *t0, double *t1, double *t2, double *t3,
                      double *xmesh, double *ymesh, double *zmesh, 
                      double *areas);
}

void
ContactSearch::GlobalSearch(SearchType search_type, 
                            Real gap_tol, int num_configs,
                            VariableHandle CURRENT_POSITION, 
                            VariableHandle PREDICTED_POSITION,
                            VariableHandle AUGMENTED_POSITION )
{
  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    timer.Start_Timer( global_search_time );
  #endif
  
  #if CONTACT_DEBUG_PRINT_LEVEL>=2
    postream<<"  Performing global search\n";
    postream<<"    global track = "<<global_tracking_interval<<"\n";
    postream<<"    track step   = "<<tracking_step<<"\n";
    #else
      #ifdef CONTACT_HEARTBEAT 
        if (contact_processor_number(SearchComm)==0) {
          std::cout<<"  Performing global search\n";
          std::cout<<"    global track = "<<global_tracking_interval<<"\n";
          std::cout<<"    track step   = "<<tracking_step<<"\n";
        }
      #endif
  #endif

  // Retrieve any tied interactions
#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    timer.Start_Timer( search_retr_tied_time );
#endif
  switch (num_configs) {
  case 1:
    Retrieve_Tied_Interactions_From_Primary(CURRENT_POSITION);
    break;
  case 2:
    Retrieve_Tied_Interactions_From_Primary(PREDICTED_POSITION);
    break;
  case 3:
    Retrieve_Tied_Interactions_From_Primary(PREDICTED_POSITION);
    break;
  default:
    POSTCONDITION(0);
    break;
  }
#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    timer.Stop_Timer( search_retr_tied_time );
#endif

  bool did_search = false;
  if (search_data->SearchForInteractions()) {

    did_search = true;

    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    timer.Start_Timer( search_topology_create_time );
    #endif
    
    switch (num_configs) {
    case 1:
      Create_Search_Topology(CURRENT_POSITION);
      break;
    case 2:
      Create_Search_Topology(PREDICTED_POSITION);
      break;
    case 3:
      Create_Search_Topology(AUGMENTED_POSITION);
      break;
    }
    if (error_code != NO_ERROR) return;

    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    timer.Start_Timer( on_processor_search_time );
    #endif

    // Get the scratch memory
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      timer.Start_Timer( search_scratch_time );
    #endif
    search_topology->Get_Scratch();
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      timer.Stop_Timer( search_scratch_time );
    #endif

    // Build the physical face list for each node
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      timer.Start_Timer( search_phys_face_time );
    #endif
    switch (physical_face_algorithm) {
    case PF_NONE:
      Build_Physical_Face_List_None(ContactSearch::SEARCH, 
                                    ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE, false);
      break;
    case PF_FACE_BASED:
      Build_Physical_Face_List_FaceWalk(ContactSearch::SEARCH, 
                                        ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE, false);
      break;
    case PF_EDGE_BASED:
      Build_Physical_Face_List_EdgeWalk(ContactSearch::SEARCH, 
                                        ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE, false);
      break;
    default:
      POSTCONDITION(0);
      break;
    }
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      timer.Stop_Timer( search_phys_face_time );
    #endif

    // Retrieve any glued interactions
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      timer.Start_Timer( search_retr_tied_time );
    #endif
    switch (num_configs) {
    case 1:
        Retrieve_Glued_Interactions(CURRENT_POSITION, ContactSearch::SEARCH);
        break;
    case 2:
        Retrieve_Glued_Interactions(PREDICTED_POSITION, ContactSearch::SEARCH);
        break;
    case 3:
        Retrieve_Glued_Interactions(PREDICTED_POSITION, ContactSearch::SEARCH);
        break;
    default:
        POSTCONDITION(0);
    break;
    }
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      timer.Stop_Timer( search_retr_tied_time );
    #endif

    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      timer.Start_Timer( search_serial_setup_time );
    #endif

    if (do_node_node_search) {
      Global_NodeNodeSearch(search_type, num_configs,
                            CURRENT_POSITION, 
                            PREDICTED_POSITION,
                            AUGMENTED_POSITION);
    }

    if (do_node_face_search) {
      Global_NodeFaceSearch(search_type, 
                            gap_tol, num_configs,
                            CURRENT_POSITION, 
                            PREDICTED_POSITION,
                            AUGMENTED_POSITION);
    }

    if (do_face_face_search || do_coverage_search) {
      Global_FaceFaceSearch(search_type, num_configs,
                            CURRENT_POSITION, 
                            PREDICTED_POSITION,
                            AUGMENTED_POSITION);
    }

    if (do_elem_elem_search) {
      Global_ElemElemSearch(search_type, num_configs,
                            CURRENT_POSITION, 
                            PREDICTED_POSITION,
                            AUGMENTED_POSITION);
    }

    #if CONTACT_DEBUG_PRINT_LEVEL>=5
      if (do_node_face_search) {
        if( contact_processor_number(SearchComm) == 0 ) {
          postream << "\n\n\n***** >>Node/Entity interactions before Interaction Definition******* \n\n";
        }
        search_topology->Display_NodeEntity_Interactions(postream);
      }
    #endif

    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      timer.Start_Timer( search_id_time );
    #endif
    Interaction_Definition( num_configs,  ContactSearch::SEARCH, ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE );
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      timer.Stop_Timer( search_id_time );
    #endif

    #if CONTACT_DEBUG_PRINT_LEVEL>=4
    if (do_node_face_search) {
      if( contact_processor_number(SearchComm) == 0 ) {
        postream << "\n\n\n***** Node/Entity interactions after Interaction Definition******* \n\n";
      }
      search_topology->Display_NodeEntity_Interactions(postream);
    }
    #endif

    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      timer.Start_Timer( search_release_scratch_time );
    #endif
    if (did_search) search_topology->Release_Scratch();
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      timer.Stop_Timer( search_release_scratch_time );
    #endif

    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      timer.Stop_Timer( on_processor_search_time );
    #endif

    if (did_search) Define_Primary_Interactions();

    if( contact_number_of_processors(SearchComm)>1 && secondary_topology ){
      #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
        timer.Start_Timer( cleanup_secondary_time );
      #endif
      if (did_search) search_topology->CleanUp();
      #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
        timer.Stop_Timer( cleanup_secondary_time );
      #endif
    }

  } // end of if (have_only_tied_interactions)
    
  #if CONTACT_DEBUG_PRINT_LEVEL>=2
    if (do_node_face_search) {
      primary_topology->Display_NodeEntity_Interactions_Summary(
          postream, ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE, (char*)"    ", 0);
    }
    #else
      #ifdef CONTACT_HEARTBEAT 
      if (do_node_face_search) {
        primary_topology->Display0_NodeEntity_Interactions_Summary(
            ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE, (char*)"    ", 0);
      }
      #endif
  #endif
    
  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    timer.Stop_Timer( global_search_time );
  #endif
}

void
ContactSearch::NewGlobalSearch(SearchType search_type, 
                               Real gap_tol, int num_configs,
                               VariableHandle CURRENT_POSITION, 
                               VariableHandle PREDICTED_POSITION,
                               VariableHandle AUGMENTED_POSITION )
{
  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    timer.Start_Timer( global_search_time );
  #endif
  
  #if CONTACT_DEBUG_PRINT_LEVEL>=2
    postream<<"  Performing new global search\n";
    postream<<"    global track = "<<global_tracking_interval<<"\n";
    postream<<"    track step   = "<<tracking_step<<"\n";
  #else
    #ifdef CONTACT_HEARTBEAT 
      if (contact_processor_number(SearchComm)==0) {
        std::cout<<"  Performing new global search\n";
        std::cout<<"    global track = "<<global_tracking_interval<<"\n";
        std::cout<<"    track step   = "<<tracking_step<<"\n";
      }
    #endif
  #endif

  if (no_ghosting==0) {
  if (tracking_step==0) {
    VariableHandle POSITION = -1;
    switch (num_configs) {
    case 1:
      POSITION = CURRENT_POSITION;
      break;
    case 2:
      POSITION = PREDICTED_POSITION;
      break;
    case 3:
      POSITION = AUGMENTED_POSITION;
      break;
    }
    #if !defined(CONTACT_NO_MPI)
    if( contact_number_of_processors(SearchComm)>1 ){
      #ifdef CONTACT_TIMINGS
        timer.Start_Timer( cleanup_secondary_time );
      #endif
      primary_topology->DeleteGhosting();
      #ifdef CONTACT_TIMINGS
        timer.Stop_Timer( cleanup_secondary_time );
      #endif
      if(do_node_face_search &&
         !do_face_face_search &&
         !do_node_node_search &&
         !do_coverage_search &&
         !do_elem_elem_search) {
        primary_topology->DoCaptureGhosting_New_NodeFace(POSITION, reasonable_gap);
      } else {  
        primary_topology->DoGhosting(POSITION, reasonable_gap);
      }
    }  
    #endif
  } else {
    #if !defined(CONTACT_NO_MPI)
    if( contact_number_of_processors(SearchComm)>1 ){
      primary_topology->UpdateGhosting();
    }
    #endif
  }}
  search_topology = primary_topology; 

  // Retrieve any tied interactions
#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    timer.Start_Timer( search_retr_tied_time );
#endif
  switch (num_configs) {
  case 1:
    Retrieve_Tied_Interactions_From_Primary(CURRENT_POSITION);
    break;
  case 2:
    Retrieve_Tied_Interactions_From_Primary(PREDICTED_POSITION);
    break;
  case 3:
    Retrieve_Tied_Interactions_From_Primary(PREDICTED_POSITION);
    break;
  default:
    POSTCONDITION(0);
    break;
  }
#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    timer.Stop_Timer( search_retr_tied_time );
#endif

  bool did_search = false;
  if (search_data->SearchForInteractions()) {

    did_search = true;

    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      timer.Start_Timer( on_processor_search_time );
    #endif

    ///if (tracking_step==0) {
      // Get the scratch memory
      #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
        timer.Start_Timer( search_scratch_time );
      #endif
      search_topology->Get_Scratch();
      #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
        timer.Stop_Timer( search_scratch_time );
      #endif
    //}

    // Build the physical face list for each node
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      timer.Start_Timer( search_phys_face_time );
    #endif 
    if (do_node_face_search) {
      switch (physical_face_algorithm) {
      case PF_NONE:
        Build_Physical_Face_List_None(ContactSearch::SEARCH, 
                                      ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE, false);
        break;
      case PF_FACE_BASED:
        if (tracking_type==GLOBAL_TRACKING && tracking_step>0) {
          Build_Physical_Face_List_FaceWalk(ContactSearch::SEARCH, 
                                            ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE, true);
        } else {
          Build_Physical_Face_List_FaceWalk(ContactSearch::SEARCH, 
                                            ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE, false);
        }
        break;
      case PF_EDGE_BASED:
        Build_Physical_Face_List_EdgeWalk(ContactSearch::SEARCH, 
                                          ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE, false);
        break;
      default:
        POSTCONDITION(0);
        break;
      }
    }
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      timer.Stop_Timer( search_phys_face_time );
    #endif

    // Retrieve any tied interactions
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      timer.Start_Timer( search_retr_tied_time );
    #endif
    switch (num_configs) {
    case 1:
      Retrieve_Glued_Interactions(CURRENT_POSITION, ContactSearch::SEARCH);
      break;
    case 2:
      Retrieve_Glued_Interactions(PREDICTED_POSITION, ContactSearch::SEARCH);
      break;
    case 3:
      Retrieve_Glued_Interactions(PREDICTED_POSITION, ContactSearch::SEARCH);
      break;
    default:
      POSTCONDITION(0);
      break;
    }
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      timer.Stop_Timer( search_retr_tied_time );
    #endif

    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      timer.Start_Timer( search_serial_setup_time );
    #endif

    if (do_node_node_search) {
      Global_NodeNodeSearch(search_type, num_configs,
                            CURRENT_POSITION, 
                            PREDICTED_POSITION,
                            AUGMENTED_POSITION);
    }

    if (do_node_face_search) {
      Global_NodeFaceSearch(search_type, 
                            gap_tol, num_configs,
                            CURRENT_POSITION, 
                            PREDICTED_POSITION,
                            AUGMENTED_POSITION);
    }

    if (do_face_face_search || do_coverage_search) {
      Global_FaceFaceSearch(search_type, num_configs,
                            CURRENT_POSITION, 
                            PREDICTED_POSITION,
                            AUGMENTED_POSITION);
    }

    if (do_elem_elem_search) {
      Global_ElemElemSearch(search_type, num_configs,
                            CURRENT_POSITION, 
                            PREDICTED_POSITION,
                            AUGMENTED_POSITION);
    }

    #if CONTACT_DEBUG_PRINT_LEVEL>=5
      if (do_node_face_search) {
        if( contact_processor_number(SearchComm) == 0 ) {
          postream << "\n\n\n***** >>Node/Entity interactions before Interaction Definition******* \n\n";
        }
        search_topology->Display_NodeEntity_Interactions(postream);
      }
    #endif

    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      timer.Start_Timer( search_id_time );
    #endif
    Interaction_Definition( num_configs,  ContactSearch::SEARCH, ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE );
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      timer.Stop_Timer( search_id_time );
    #endif

    #if CONTACT_DEBUG_PRINT_LEVEL>=4
      if (do_node_face_search) {
        if( contact_processor_number(SearchComm) == 0 ) {
          postream << "\n\n\n***** Node/Entity interactions after Interaction Definition******* \n\n";
        }
        search_topology->Display_NodeEntity_Interactions(postream);
      }
    #endif

    //if (global_tracking_interval>1 && tracking_step>=global_tracking_interval-1) {
      #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
        timer.Start_Timer( search_release_scratch_time );
      #endif
      if (did_search) search_topology->Release_Scratch();
      #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
        timer.Stop_Timer( search_release_scratch_time );
      #endif
    //}

    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      timer.Stop_Timer( on_processor_search_time );
    #endif

    //Define_Primary_Interactions();

    #if !defined(CONTACT_NO_MPI)
    if (no_ghosting ==0) {
    if (contact_number_of_processors(SearchComm)>1 && tracking_step==0) {
      #ifdef CONTACT_TIMINGS
        timer.Start_Timer( cleanup_secondary_time );
      #endif
      if (did_search) primary_topology->UpdateGhostingSetupForNoSecondary();
      #ifdef CONTACT_TIMINGS
        timer.Stop_Timer( cleanup_secondary_time );
      #endif
    }}
    #endif

  } // end of if (have_only_tied_interactions)
    
  #if CONTACT_DEBUG_PRINT_LEVEL>=2
    if (do_node_face_search) {
      primary_topology->Display_NodeEntity_Interactions_Summary(
          postream, ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE, (char*)"    ", 0);
    }
  #else
    #ifdef CONTACT_HEARTBEAT 
      if (do_node_face_search) {
      primary_topology->Display0_NodeEntity_Interactions_Summary(
          ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE, (char*)"    ", 0);
      }
    #endif
  #endif
   
  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    timer.Stop_Timer( global_search_time );
  #endif
}

void
ContactSearch::Global_NodeNodeSearch(SearchType search_type, int num_configs,
                                     VariableHandle CURRENT_POSITION, 
                                     VariableHandle PREDICTED_POSITION,
                                     VariableHandle AUGMENTED_POSITION )
{
#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Start_Timer( search_serial_setup_time );
#endif
  VariableHandle POSITION;
  switch (num_configs) {
  case 1:
    POSITION = CURRENT_POSITION;
    break;
  case 2:
    POSITION = PREDICTED_POSITION;
    break;
  case 3:
    POSITION = AUGMENTED_POSITION;
    break;
  }
  VariableHandle NODE_RADIUS =
    search_topology->Variable_Handle( ContactTopology::Node_Radius );
  
  int number_of_nodes = search_topology->Number_of_Nodes();
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(search_topology->NodeList()->EntityList());
  int* list    = new int[number_of_nodes];
  Real* dist   = new Real[number_of_nodes];
  Real node_rad_tol = search_data->Max_Search_Tolerance();
  //
  //  Build a bounding box hierarchy containing all nodes.
  //
  ObjectBoundingBox *node_boxes = new ObjectBoundingBox[number_of_nodes];
  for(int inode = 0; inode < number_of_nodes; ++inode) {
    Real* n_position = Nodes[inode]->Variable(POSITION);
    Real node_rad = (*(Nodes[inode]->Variable(NODE_RADIUS))) + node_rad_tol;
    node_boxes[inode].set_object_number(inode);
    node_boxes[inode].add_sphere(n_position, node_rad);
  }

  int hierarchy_size = 2 * number_of_nodes - 1;
  ObjectBoundingBoxHierarchy *node_hierarchy = NULL;
  if(hierarchy_size > 0) {
    node_hierarchy = new ObjectBoundingBoxHierarchy[hierarchy_size];
    ObjectBoundingBoxHierarchy::create_hierarchy(node_hierarchy, node_boxes, number_of_nodes);
  }
  delete [] node_boxes;

#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Stop_Timer( search_serial_setup_time );
  timer.Start_Timer( search_main_time );
#endif
  //
  //  Do a cache optimized loop of the nodes in the same order as the hierarchical searching arrays
  //
  ObjectBoundingBoxHierarchy *last_ptr = &(node_hierarchy[hierarchy_size-1]);
  for(int inode = 0; inode < hierarchy_size; ++inode) {
    const ObjectBoundingBoxHierarchy &cur_leaf = node_hierarchy[inode];
    const int &object_num = cur_leaf.get_object_number();
    if(object_num < 0) continue;
    ContactNode<Real>* node0 = Nodes[object_num];
    int node0_key = node0->Entity_Key();
    if (node0_key<0) continue;
    int list_size = 0;

    node_hierarchy->search_for_overlap_recurse_sym(cur_leaf, list, list_size, last_ptr);

    list_size = Cull_Node_List(list_size, list, node0, dist);
    for(int j=0 ; j<list_size ; ++j ){
      ContactNode<Real>* node1 = Nodes[list[j]];
      int node1_key = node1->Entity_Key();

      if(node1_key >= 0) {
        if(node1->Ownership() == ContactTopologyEntity<Real>::OWNED) {
          ContactNodeNodeInteraction* cnni0 =
            ContactNodeNodeInteraction::new_ContactNodeNodeInteraction(
                             allocators[ALLOC_ContactNodeNodeInteraction],
                             node1, node0, node1_key, dist[j] );
          node1->Add_NodeNode_Interaction( cnni0 );
	}
      }

      //
      //  NKC NOTE, generally probably only want an assymetric portion of interactions.
      //  Since A interacts with B, B interacts with A can halve the storage and significantly 
      //  reduce computation time
      //
      if(node0->Ownership() == ContactTopologyEntity<Real>::OWNED) {
        ContactNodeNodeInteraction* cnni1 =
          ContactNodeNodeInteraction::new_ContactNodeNodeInteraction(
                           allocators[ALLOC_ContactNodeNodeInteraction],
                           node0, node1, node0_key, dist[j] );
        node0->Add_NodeNode_Interaction( cnni1 );
      }
    }
  }

  if(node_hierarchy != NULL) delete [] node_hierarchy;


#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Stop_Timer( search_main_time );
#endif

#if CONTACT_DEBUG_PRINT_LEVEL>=6
  postream << "\n\n\n***** Node-Node interactions before Interaction Definition *******\n\n";
  search_topology->Display_NodeNode_Interactions(postream);
#endif

#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Start_Timer( search_id_time );
#endif
  Interaction_Definition( num_configs, ContactSearch::SEARCH, ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE );
#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Stop_Timer( search_id_time );
#endif
  
#if CONTACT_DEBUG_PRINT_LEVEL>=6
  postream << "\n\n\n***** Node-Node Interactions  *****\n\n";
  search_topology->Display_NodeNode_Interactions(postream);
#endif

  delete [] list;
  delete [] dist;
}

void
ContactSearch::Global_NodeFaceSearch(SearchType search_type, 
                                     Real gap_tol, int num_configs,
                                     VariableHandle CURRENT_POSITION, 
                                     VariableHandle PREDICTED_POSITION,
                                     VariableHandle AUGMENTED_POSITION )
{
  #ifdef CONTACT_HEARTBEAT 
    if (contact_processor_number(SearchComm)==0) {
      cout<<"    Global_NodeFaceSearch\n";
    }
  #endif
  VariableHandle NODE_NORMAL = 
    primary_topology->Variable_Handle( ContactTopology::Node_Normal );
  VariableHandle FACE_NORMAL = 
  primary_topology->Variable_Handle( ContactTopology::Face_Normal );

  int number_of_nodes = search_topology->Number_of_Nodes();
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(search_topology->NodeList()->EntityList());
  int number_of_faces = search_topology->Number_of_Faces();
  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(search_topology->FaceList()->EntityList());

  //
  //  Allocate data required by the search
  //
  Real user_search_tol     = search_data->Max_Search_Tolerance();
  Real user_tang_tol       = search_data->Max_Search_Tangential_Tolerance();
  int  num_search_entities = search_data->Num_Search_Entities();

  vector<int>  physical_faces(number_of_nodes);
  vector<Real> norm_search_tol(num_search_entities);
  vector<Real> tang_search_tol(num_search_entities);
  vector<bool> valid_inter(num_search_entities+1);

  int*  node_search_status = new int[number_of_nodes];
  int*  face_search_status = new int[number_of_faces];
  
  if (tracking_type==GLOBAL_TRACKING && tracking_step>0) {
    build_nodeface_search_status_proximity(node_search_status, face_search_status);
  } else {
    build_nodeface_search_status(node_search_status, face_search_status);
  }
  int num_analytic_surfaces = primary_topology->Number_of_Analytic_Surfaces();
  #if CONTACT_DEBUG_PRINT_LEVEL>=2
    postream<<"      number of active nodes = "<<num_active_nodes<<"\n";
    postream<<"      number of active faces = "<<num_active_faces<<"\n";
  #else
  #ifdef CONTACT_HEARTBEAT 
    int number_nodes1 = contact_global_sum(num_active_nodes, SearchComm );
    int number_faces1 = contact_global_sum(num_active_faces, SearchComm );
    if (contact_processor_number(SearchComm)==0) {
      cout<<"      number of active nodes = "<<number_nodes1<<"\n";
      cout<<"      number of active faces = "<<number_faces1<<"\n";
    }
  #endif
  #endif
  
  ContactTopologyEntityList* face_list = search_topology->FaceList();
  const int num_face_blocks = search_topology->Number_of_Face_Blocks();
  
  VariableHandle NODE_COORD_START;
  VariableHandle NODE_COORD_END;
  switch (num_configs) {
  case 1:
    NODE_COORD_START = CURRENT_POSITION;
    NODE_COORD_END   = CURRENT_POSITION;
    break;
  case 2:
    NODE_COORD_START = CURRENT_POSITION;
    NODE_COORD_END   = PREDICTED_POSITION;
    break;
  case 3:
    NODE_COORD_START = CURRENT_POSITION;
    NODE_COORD_END   = AUGMENTED_POSITION;
    break;
  }
  
  if (tracking_type==GLOBAL_TRACKING) {
    int num_faces_in_proximity = 0;
    int num_nodes_in_proximity = 0;
    if (tracking_step==0) {
      //
      // do expanded 'capture' search here
      //
      #if CONTACT_DEBUG_PRINT_LEVEL>=2
        postream<<"      Initializing capture search\n";
      #else
      #ifdef CONTACT_HEARTBEAT 
        if (contact_processor_number(SearchComm)==0) {
          cout<<"      Initializing capture search\n";
        }
      #endif
      #endif
      //
      //  Determine the type of range search to use
      //
      ContactNodeFaceRangeSearch::SearchType range_search_type = ACME::determine_range_search_type(search_data);
      //
      //  Construct the range search helper object and 
      //  do any initialization that the object requires
      //
      Real capture_tol = capture_motion;
      ContactNodeFaceRangeSearch capture_search(range_search_type, 
   						search_topology, 
   						num_active_nodes, 
   						node_search_status,
   						this,
   						NODE_COORD_START,
   						NODE_COORD_END,
   						capture_tol);
      #if CONTACT_DEBUG_PRINT_LEVEL>=2 || defined(CONTACT_HEARTBEAT)
      int max_proximity_nodes = 0;
      #endif  
      if(capture_search.valid_search()) {
        #if CONTACT_DEBUG_PRINT_LEVEL>=2
          postream<<"      Performing capture search, tol = "<<capture_tol<<"\n";
        #else
        #ifdef CONTACT_HEARTBEAT 
          if (contact_processor_number(SearchComm)==0) {
            cout<<"      Performing capture search, tol = "<<capture_tol<<"\n";
          }
        #endif
        #endif
   	//
   	//  Loop over all the face blocks in the mesh
   	//
   	ACME::Int_Vector	 node_keys(num_active_nodes);
   	ACME::ContactNode_Vector node_list(num_active_nodes);
   	int global_face_index = 0;
   	Real box_expand = user_search_tol+capture_tol;
   	for(int iface_block = 0; iface_block < num_face_blocks; ++iface_block) {
   	  ContactFace<Real>** BlockFaces = reinterpret_cast<ContactFace<Real>**>(face_list->BlockEntityList(iface_block));
   	  int num_faces = face_list->BlockNumEntities(iface_block);
   	  //
   	  //  Determine the search tolerances and valid interactions for each face, node block pair
   	  //
   	  valid_inter[0] = true;
   	  for(int inode_block = 0; inode_block < num_search_entities; ++inode_block) {
   	    Real norm_tol = search_data->Get_Search_Data(ContactSearch::SEARCH_NORMAL_TOLERANCE, inode_block, iface_block );
   	    Real tang_tol = search_data->Get_Search_Data(ContactSearch::SEARCH_TANGENTIAL_TOLERANCE, inode_block, iface_block );
   	    norm_search_tol[inode_block] = norm_tol;
   	    tang_search_tol[inode_block] = tang_tol;
   	    if((int)search_data->Get_Search_Data(ContactSearch::INTERACTION_TYPE, inode_block, iface_block ) == NO_INTERACTION) {
   	      valid_inter[inode_block + 1] = false;
   	    } else {
   	      valid_inter[inode_block + 1] = true;
   	    }
   	  }
   	  //
   	  //  Loop over all of the faces contained in this face block, compute interactions on each
   	  //
   	  for(int iface = 0; iface < num_faces; ++iface, ++global_face_index) {
   	    //
   	    //  If a face is not active for the search at the current time step, 
   	    //  do not perform any searching, move on to the next face
   	    //
   	    if(!face_search_status[global_face_index]) continue;
   	    ContactBoundingBox face_current_box;
   	    ContactBoundingBox face_predicted_box;
   	    ContactBoundingBox face_search_box;
   	    ContactFace<Real> *face = BlockFaces[iface];
   	    face->ComputeBoundingBoxForSearch(num_configs,
   					      NODE_COORD_START,
   					      NODE_COORD_END, 
   					      auto_tol,
   					      box_inflation,
   					      box_expand,
   					      face_current_box,
   					      face_predicted_box,
   					      face_search_box);
   	    //
   	    //  Determine all potential face/node overlaps using a hybrid range search.  Perform a tree
   	    //  search on the entity bounding boxes follwed immediatly by a tree search on the 
   	    //  relavent node entity trees
   	    //
   	    capture_search.search_for_overlap(face_search_box,
   					      node_keys,
   					      node_list,
   					      valid_inter);
   	    if (node_list.size()==0) continue;
   	    //
   	    //  Perform additional list culling
   	    //
   	    int list_size = Cull_Node_List_New1(node_list,
   						face,
   						NODE_NORMAL, 
   						FACE_NORMAL, 
   						num_configs, 
   						face_current_box,
   						face_predicted_box,
   						valid_inter,
   						NODE_COORD_START,
   						NODE_COORD_END);
   	    if (list_size>0) {
              #if CONTACT_DEBUG_PRINT_LEVEL>=2 || defined(CONTACT_HEARTBEAT)
              max_proximity_nodes = std::max(max_proximity_nodes,list_size);
              #endif
              face->in_proximity = 1;
              for( int k=0 ; k<list_size ; ++k ){
   		node_list[k]->in_proximity = 1;
   	      }
   	    }
   	  }
   	}
      }
      for (int n=0; n<number_of_nodes; ++n) {
   	if (Nodes[n]->in_proximity) ++num_nodes_in_proximity;
      }
      for (int n=0; n<number_of_faces; ++n) {
   	if (Faces[n]->in_proximity) ++num_faces_in_proximity;
      }
      build_nodeface_search_status_proximity(node_search_status, face_search_status);
      #if CONTACT_DEBUG_PRINT_LEVEL>=2
      postream<<"        number of nodes in proximity = "<<num_nodes_in_proximity<<"\n";
      postream<<"        number of faces in proximity = "<<num_faces_in_proximity<<"\n";
      postream<<"        max nodes in proximity/face  = "<<max_proximity_nodes<<"\n";
      postream<<"        avg nodes in proximity/face  = "<<num_nodes_in_proximity/num_faces_in_proximity<<"\n";
      #else
        #ifdef CONTACT_HEARTBEAT 
        int global_max_proximity_nodes = contact_global_maximum(max_proximity_nodes,SearchComm);
        if (contact_processor_number(SearchComm)==0) {
          cout<<"        number of nodes in proximity = "<<num_nodes_in_proximity<<"\n";
          cout<<"        number of faces in proximity = "<<num_faces_in_proximity<<"\n";
          cout<<"        max nodes in proximity/face  = "<<global_max_proximity_nodes<<"\n";
          cout<<"        avg nodes in proximity/face  = "<<num_nodes_in_proximity/num_faces_in_proximity<<"\n";
        }
        #endif
      #endif
    }
  }
  if (num_active_nodes && (num_active_faces || num_analytic_surfaces)) {
    #if CONTACT_DEBUG_PRINT_LEVEL>=2
      if (tracking_type==GLOBAL_TRACKING) {
        postream<<"      Performing proximity search\n";
        postream<<"        number of active nodes = "<<num_active_nodes<<"\n";
        postream<<"        number of active faces = "<<num_active_faces<<"\n";
      } else {
        postream<<"      Performing range search check\n";
      }
    #else
    #ifdef CONTACT_HEARTBEAT 
      if (contact_processor_number(SearchComm)==0) {
        if (tracking_type==GLOBAL_TRACKING) {
          cout<<"      Performing proximity search\n";
          cout<<"        number of active nodes = "<<num_active_nodes<<"\n";
          cout<<"        number of active faces = "<<num_active_faces<<"\n";
        } else {
          cout<<"      Performing range search check\n";
        }
      }
    #endif
    #endif
    //
    //  Determine the type of range search to use
    //
    ContactNodeFaceRangeSearch::SearchType range_search_type = ACME::determine_range_search_type(search_data);
    //
    //  Construct the range search helper object and 
    //  do any initialization that the object requires
    //
    ContactNodeFaceRangeSearch range_search(range_search_type, 
  					    search_topology, 
  					    num_active_nodes, 
  					    node_search_status,
  					    this,
  					    NODE_COORD_START,
  					    NODE_COORD_END);   

#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    timer.Stop_Timer( search_serial_setup_time );
    timer.Start_Timer( search_main_time );
#endif

    if(range_search.valid_search()) {
      #if CONTACT_DEBUG_PRINT_LEVEL>=2
        postream<<"      Performing range search\n";
      #else
      #ifdef CONTACT_HEARTBEAT 
        if (contact_processor_number(SearchComm)==0) {
          cout<<"      Performing range search\n";
        }
      #endif
      #endif
      //
      //  Loop over all the face blocks in the mesh
      //
      ACME::Int_Vector         node_keys(num_active_nodes);
      ACME::ContactNode_Vector node_list(num_active_nodes);
      if (num_active_faces) {
        #if CONTACT_DEBUG_PRINT_LEVEL>=2
        postream<<"      Performing range search for faces\n";
        #else
        #ifdef CONTACT_HEARTBEAT 
        if (contact_processor_number(SearchComm)==0) {
          cout<<"      Performing range search for faces\n";
        }
        #endif
        #endif
        int global_face_index = 0;
        for(int iface_block = 0; iface_block < num_face_blocks; ++iface_block) {
          ContactFace<Real>** BlockFaces = reinterpret_cast<ContactFace<Real>**>(face_list->BlockEntityList(iface_block));
          int num_faces = face_list->BlockNumEntities(iface_block);
          //
          //  Determine the search tolerances and valid interactions for each face, node block pair
          //
          valid_inter[0] = true;
          for(int inode_block = 0; inode_block < num_search_entities; ++inode_block) {
            Real norm_tol = search_data->Get_Search_Data(ContactSearch::SEARCH_NORMAL_TOLERANCE, inode_block, iface_block );
            Real tang_tol = search_data->Get_Search_Data(ContactSearch::SEARCH_TANGENTIAL_TOLERANCE, inode_block, iface_block );
            norm_search_tol[inode_block] = norm_tol;
            tang_search_tol[inode_block] = tang_tol;
            if((int)search_data->Get_Search_Data(ContactSearch::INTERACTION_TYPE, inode_block, iface_block ) == NO_INTERACTION) {
              valid_inter[inode_block + 1] = false;
            } else {
              valid_inter[inode_block + 1] = true;
            }
          }
          //
          //  Loop over all of the faces contained in this face block, compute interactions on each
          //
          for(int iface = 0; iface < num_faces; ++iface, ++global_face_index) {
            //
            //  If a face is not active for the search at the current time step, do not perform any searching,
            //  move on to the next face
            //
            if(!face_search_status[global_face_index]) continue;
            ContactBoundingBox face_current_box;
            ContactBoundingBox face_predicted_box;
            ContactBoundingBox face_search_box;
            ContactFace<Real> *face = BlockFaces[iface];
            face->ComputeBoundingBoxForSearch(num_configs,
                                              NODE_COORD_START,
                                              NODE_COORD_END, 
                                              auto_tol,
                                              box_inflation,
                                              user_search_tol,
                                              face_current_box,
                                              face_predicted_box,
                                              face_search_box);
            //
            //  Determine all potential face/node overlaps using a hybrid range search.  Perform a tree
            //  search on the entity bounding boxes follwed immediatly by a tree search on the 
            //  relavent node entity trees
            //
            range_search.search_for_overlap(face_search_box,
                                            node_keys,
                                            node_list,
                                            valid_inter);
            if (node_list.size()==0) continue;
            
            //
            //  Perform additional list culling
            //
            int list_size = Cull_Node_List_New(node_list,
                                               &(physical_faces[0]), 
                                               face,
                                               NODE_NORMAL, 
                                               FACE_NORMAL, 
                                               num_configs, 
                                               face_current_box,
                                               face_predicted_box,
                                               node_keys,
                                               norm_search_tol,
                                               tang_search_tol,
                                               valid_inter,
                                               NODE_COORD_START,
                                               NODE_COORD_END);
            
            if (list_size==0) continue;
            
            Real REL_TANG_TOL_VAR = REL_TANG_TOL;
            switch (search_type) {
            case STATIC1CONFIG:
            case STATIC1CONFIG_TIED:
              process_node_face_interactions(list_size, 
                                             node_list, 
                                             face, 
                                             CURRENT_POSITION,
                                             Nodes,
                                             NODE_NORMAL,
                                             REL_TANG_TOL_VAR,
                                             user_tang_tol,
                                             gap_tol,
                                             &(physical_faces[0]),
                                             node_keys,
                                             ContactNodeEntityInteraction::CLOSEST_POINT_PROJECTION_1 );  
              break;    
            case STATIC2CONFIG:
              process_node_face_interactions(list_size, 
                                             node_list, 
                                             face, 
                                             PREDICTED_POSITION,
                                             Nodes,
                                             NODE_NORMAL,
                                             REL_TANG_TOL_VAR,
                                             user_tang_tol,
                                             gap_tol,
                                             &(physical_faces[0]),
                                             node_keys,
                                             ContactNodeEntityInteraction::CLOSEST_POINT_PROJECTION_2 );  
              break;    
            case DYNAMIC2CONFIG:
              ComputeInteractions(list_size, 
                                  node_list,
                                  face,
                                  CURRENT_POSITION,
                                  AUGMENTED_POSITION,
                                  PREDICTED_POSITION,
                                  NODE_NORMAL,
                                  user_search_tol,
                                  user_tang_tol,
                                  REL_TANG_TOL_VAR,
                                  gap_tol,
                                  &(physical_faces[0]),
                                  node_keys); 
              break;
            default:
              POSTCONDITION(0);   
              break;    
            }
          }
        }
      }
      if (num_analytic_surfaces) {
        #if CONTACT_DEBUG_PRINT_LEVEL>=2
        postream<<"      Performing range search for analytic surfaces\n";
        #else
        #ifdef CONTACT_HEARTBEAT 
        if (contact_processor_number(SearchComm)==0) {
          cout<<"      Performing range search for analytic surfaces\n";
        }
        #endif
        #endif
        switch (search_type) {
        case STATIC1CONFIG:
        case STATIC1CONFIG_TIED:
          range_search.process_analytic_surfaces(ContactNodeEntityInteraction::CLOSEST_POINT_PROJECTION_1, CURRENT_POSITION);
          break;
        case STATIC2CONFIG:
          range_search.process_analytic_surfaces(ContactNodeEntityInteraction::CLOSEST_POINT_PROJECTION_2, PREDICTED_POSITION);
          break;
        case DYNAMIC2CONFIG:
          range_search.process_analytic_surfaces(ContactNodeFaceInteraction::MOVING_INTERSECTION, PREDICTED_POSITION, CURRENT_POSITION);
          break;
        default:
          POSTCONDITION(0);   
          break;   
        }   
      }
    }
  }
  delete [] node_search_status;
  delete [] face_search_status;

#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Stop_Timer( search_main_time );
#endif
}

void
ContactSearch::Global_FaceFaceSearch(SearchType search_type, int num_configs,
                                     VariableHandle CURRENT_POSITION, 
                                     VariableHandle PREDICTED_POSITION,
                                     VariableHandle AUGMENTED_POSITION )
{
#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Start_Timer( search_serial_setup_time );
#endif
  VariableHandle POSITION;
  switch (num_configs) {
  case 1:
    POSITION = CURRENT_POSITION;
    break;
  case 2:
    POSITION = PREDICTED_POSITION;
    break;
  case 3:
    POSITION = AUGMENTED_POSITION;
    break;
  }
  VariableHandle NODE_NORMAL = 
    primary_topology->Variable_Handle( ContactTopology::Node_Normal );
  VariableHandle FACE_NORMAL = 
    primary_topology->Variable_Handle( ContactTopology::Face_Normal );
  VariableHandle CENTROID =
    primary_topology->Variable_Handle( ContactTopology::Centroid );
  VariableHandle CHARACTERISTIC_LENGTH =
    primary_topology->Variable_Handle( ContactTopology::Characteristic_Length );
  
  int number_of_nodes = search_topology->Number_of_Nodes();
  int number_of_faces = search_topology->Number_of_Faces();
  int number_of_edges = search_topology->Number_of_Edges();
  
  //
  //  Determine the user search tolerance
  //
  Real user_search_tol = std::max(search_data->Max_Search_Normal_Tolerance(),
                             search_data->Max_Search_Tangential_Tolerance() );
  //
  //  Build a bounding box hierarchy containing all faces
  //
  ContactFace<Real>** Faces = reinterpret_cast<ContactFace<Real>**>(search_topology->FaceList()->EntityList());
  ObjectBoundingBox *face_boxes = new ObjectBoundingBox[number_of_faces];
  ObjectBoundingBox *face_boxes_tmp = new ObjectBoundingBox[number_of_faces];

#if (MAX_FFI_DERIVATIVES > 0)
  ContactFace<ActiveScalar>::Initialize_Lookup_Arrays();
  ContactEdge<ActiveScalar>::Initialize_Lookup_Arrays();
#endif

  for(int iface = 0; iface < number_of_faces; ++iface) {
    ContactFace<Real>* face = Faces[iface];
    face_boxes[iface].set_object_number(iface);
    for(int inode = 0; inode < face->Nodes_Per_Face(); ++inode) {
      //
      //  Expand the face bounding box to encompass all of the faces nodes
      //
      Real *node_position = face->Node(inode)->Variable(POSITION);
      face_boxes[iface].add_point(node_position);
    }
    //
    //  Expand the face bounding box by half of the tolerance value.  Since a symmetric tree on tree search is performed
    //  Using half the search tolerance on all objects is equivalent to using the whole search tolerance on just the search
    //  object
    //
    face_boxes[iface].add_tolerance(user_search_tol/2.0);
    face_boxes_tmp[iface] = face_boxes[iface];
  }
  int hierarchy_size = 2 * number_of_faces - 1;
  ObjectBoundingBoxHierarchy *face_hierarchy = NULL;

  if(hierarchy_size > 0) {
    face_hierarchy = new ObjectBoundingBoxHierarchy[hierarchy_size];
    ObjectBoundingBoxHierarchy::create_hierarchy(face_hierarchy, face_boxes_tmp, number_of_faces);
  }    

  // Process any face-face interactions
  int* list = new int[number_of_faces];
  
  // Loop over all the faces
  int old_block_id = -1;
  int new_block_id = -1;
  ContactElem<Real>* master_element = NULL;
#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Stop_Timer( search_serial_setup_time );
#endif

#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Start_Timer( search_main_time );
#endif
  //
  //  Do a cache optimized loop of the nodes in the same order as the hierarchical searching arrays
  //
  //
  //  NKC NOTE, though slower, switched this over to a direct loop to ensure 
  //  that the regression tests get identical awnswers in serial and parallel
  //
  //for(int iface = 0; iface < hierarchy_size; ++iface) {
  //const ObjectBoundingBoxHierarchy &cur_leaf = face_hierarchy[iface];
  //const int &object_num = cur_leaf.get_object_number();
  //if(object_num < 0) continue;
  //ContactFace<Real>* master_face = Faces[object_num];

  for(int iface=0 ; iface<number_of_faces ; ++iface ){
    ContactFace<Real>* master_face = Faces[iface];

    

    new_block_id = master_face->BlockID();
    if (new_block_id!=old_block_id) {
      // Delete the old master element if it exists
      if (master_element!=NULL) {
   	master_element->DeleteTopology(allocators);
   	master_element->~ContactElem<Real>();
   	switch (search_topology->Face_Block(old_block_id)->Type()) {
   	case TRIFACEL3:
        case SHELLTRIFACEL3:
   	  allocators[ALLOC_ContactWedgeElemL6].Delete_Frag(master_element);
   	  break;
   	case QUADFACEL4:
        case SHELLQUADFACEL4:
   	  allocators[ALLOC_ContactHexElemL8].Delete_Frag(master_element);
   	  break;
   	default:
   	  PRECONDITION(0);
   	  break;  
   	}
        master_element = NULL;
      }
      // Contstruct the master element from the master face type
      switch (search_topology->Face_Block(new_block_id)->Type()) {
      case TRIFACEL3:
      case SHELLTRIFACEL3:
   	master_element = ContactWedgeElemL6<Real>::new_ContactWedgeElemL6( 
   				      allocators[ALLOC_ContactWedgeElemL6],
   				      ContactSearch::WEDGEELEML6, 0);
   	break;
      case QUADFACEL4:
      case SHELLQUADFACEL4:
   	master_element = ContactHexElemL8<Real>::new_ContactHexElemL8( 
   				      allocators[ALLOC_ContactHexElemL8],
   				      ContactSearch::HEXELEML8, 0);
   	break;
      default:
   	PRECONDITION(0);
   	break;  
      }
      master_element->BuildTopology(number_of_nodes,
   				    number_of_edges,
   				    number_of_faces,
   				    allocators);
      old_block_id = new_block_id;  
    }
    //
    // Determine all geometrically possible interactions between this face and all other faces in the mesh
    //
    int list_size = 0;
    if(face_hierarchy != NULL) ObjectBoundingBoxHierarchy::search_for_overlap_loop(face_hierarchy,face_boxes[iface], list, list_size);
    if( list_size ){
      if (no_parallel_consistency==INACTIVE) Sort(list_size, list);
      for (int ilist=0; ilist<list_size; ++ilist) {
        ContactFace<Real>* slave_face = Faces[list[ilist]];
	//
	//  Cull out the found interactions based on various criteria
	//
        //  Test 1: Cull the face if no interaction is requested between these enetities
	//
    	int  master_face_key = master_face->Entity_Key();
    	int  slave_face_key  = slave_face->Entity_Key();
        int Interaction_Type = (int)
             search_data->Get_Search_Data( ContactSearch::INTERACTION_TYPE,
                                           slave_face_key, master_face_key );
        if( Interaction_Type == NO_INTERACTION ) continue;
	//
	//  Test 2: Cull out the faces that are not owned by this processor
	//
        if( slave_face->Ownership() != ContactTopologyEntity<Real>::OWNED ) continue;
        //
	//  Test 3: Cull out the face if it is the same as the master face
	//
        if (slave_face->ProcArrayIndex() == master_face->ProcArrayIndex()) continue;
        //
	//  Test 4: Cull out the face if is not "opposed" to the master face
	//
        Real* face1_normal = slave_face->Variable(FACE_NORMAL);
        Real* face2_normal = master_face->Variable(FACE_NORMAL);
        Real dot_product = 0;
        for( int idim=0 ; idim<dimensionality ; ++idim ) dot_product += face1_normal[idim]*face2_normal[idim];
        if( dot_product >= 0 ) continue;
	//
	//  Add the face-face interaction
	//

        Real normal_search_tol = search_data->Get_Search_Data( ContactSearch::SEARCH_NORMAL_TOLERANCE,
    				                               slave_face_key, master_face_key );
  	master_element->UpdateTopology(master_face, POSITION, 
    				       FACE_NORMAL, NODE_NORMAL,
    				       normal_search_tol);

        ContactFaceFaceInteraction<Real> *cffi = Face_Face_Search(slave_face, master_face, master_element, POSITION, allocators);
        if(cffi) {
          if(compute_partials == ACTIVE) {
            switch(computed_partials_order) {
              default : {
                if(slave_face->FaceType() == TRIFACEL3 && master_face->FaceType() == TRIFACEL3) { 
                  // delete cffi and recompute interaction and 1st order (hand-coded) partial derivatives
                  cffi->~ContactFaceFaceInteraction<Real>();
                  allocators[ALLOC_ContactFaceFaceInteraction].Delete_Frag(cffi);
                  cffi = Partial_Face_Face_Search(slave_face, master_face, master_element, POSITION, normal_search_tol, allocators);
                }
                else {
#if (MAX_FFI_DERIVATIVES > 0)
                  // compute 1st order partial derivatives using automatic differentiation
                  // convert the slave face to ActiveScalar
                  ContactFace<ActiveScalar> *active_slave_face;
                  switch (slave_face->FaceType()) {
                  case TRIFACEL3: 
                    active_slave_face = ContactTriFaceL3<ActiveScalar>::new_ContactTriFaceL3(active_allocators);
                    break;
                  case QUADFACEL4:
                    active_slave_face = ContactQuadFaceL4<ActiveScalar>::new_ContactQuadFaceL4(active_allocators);
                    break;
                  case SHELLTRIFACEL3:
                    active_slave_face = ContactShellTriFaceL3<ActiveScalar>::new_ContactShellTriFaceL3(active_allocators);
                    break;
                  case SHELLQUADFACEL4:
                    active_slave_face = ContactShellQuadFaceL4<ActiveScalar>::new_ContactShellQuadFaceL4(active_allocators);
                    break;
                  default:
                    PRECONDITION(0);
                    break;
                  }
                  for( int i=0; i<active_slave_face->Nodes_Per_Face(); ++i ) {
                    ContactNode<ActiveScalar>* node = ContactNode<ActiveScalar>::new_ContactNode(active_allocators, slave_face->Node(i)->NodeType());
                    Real* pos = slave_face->Node(i)->Variable(POSITION);
                    ActiveScalar* active_pos = node->Variable(POSITION);
                    for( int j=0 ; j<dimensionality ; ++j ) {
                      active_pos[j] = InitActiveScalar(MAX_FFI_DERIVATIVES, 3*i+j, pos[j]);
                    }
                    active_slave_face->ConnectNode(i, node);
                    node->Connect_Face(active_slave_face); // ?
                  }
                  for( int i=0; i<active_slave_face->Edges_Per_Face(); ++i ) {
                    ContactEdge<ActiveScalar>* edge = ContactLineEdgeL2<ActiveScalar>::new_ContactLineEdgeL2(
                                  active_allocators[ContactSearch::ALLOC_ContactLineEdgeL2],
                                  ContactSearch::LINEEDGEL2);
          
                    ContactNode<Real> *node[3];
                    slave_face->Get_Edge_Nodes(i, node);
                    ContactNode<ActiveScalar> *active_node[3];
                    active_slave_face->Get_Edge_Nodes(i, active_node);
                    if(slave_face->Edge(i)->Node(0) == node[0]) {
                      for (int ii=0; ii<edge->Nodes_Per_Edge(); ++ii)
                        edge->ConnectNode(ii,active_node[ii]);
                    }
                    else { // some edges can be reversed
                      for (int ii=0; ii<edge->Nodes_Per_Edge(); ++ii)
                        edge->ConnectNode(edge->Nodes_Per_Edge()-ii-1,active_node[ii]);
                    }
                    active_slave_face->ConnectEdge(i, edge);
                    edge->ConnectFace(active_slave_face);
                  }
                  active_slave_face->Compute_Centroid(POSITION, CENTROID);
                  active_slave_face->Compute_Normal(POSITION, FACE_NORMAL);
                  active_slave_face->Compute_CharacteristicLength(POSITION, CHARACTERISTIC_LENGTH);
          
                  // convert the master face
                  ContactFace<ActiveScalar> *active_master_face;
                  switch (master_face->FaceType()) {
                  case TRIFACEL3:
                    active_master_face = ContactTriFaceL3<ActiveScalar>::new_ContactTriFaceL3(active_allocators);
                    break;
                  case QUADFACEL4:
                    active_master_face = ContactQuadFaceL4<ActiveScalar>::new_ContactQuadFaceL4(active_allocators);
                    break;
                  case SHELLTRIFACEL3:
                    active_master_face = ContactShellTriFaceL3<ActiveScalar>::new_ContactShellTriFaceL3(active_allocators);
                    break;
                  case SHELLQUADFACEL4:
                    active_master_face = ContactShellQuadFaceL4<ActiveScalar>::new_ContactShellQuadFaceL4(active_allocators);
                    break;
                  default:
                    PRECONDITION(0);
                    break;
                  }
                  for( int i=0; i<active_master_face->Nodes_Per_Face(); ++i ) {
                    ContactNode<ActiveScalar>* node = ContactNode<ActiveScalar>::new_ContactNode(active_allocators, master_face->Node(i)->NodeType());
                    Real* pos = master_face->Node(i)->Variable(POSITION);
                    ActiveScalar* active_pos = node->Variable(POSITION);
                    for( int j=0 ; j<dimensionality ; ++j ) {
                      active_pos[j] = InitActiveScalar(MAX_FFI_DERIVATIVES, MAX_FFI_DERIVATIVES/2+3*i+j, pos[j]);
                    }
                    active_master_face->ConnectNode(i, node);
                    node->Connect_Face(active_master_face);
                  }
                  for( int i=0; i<active_master_face->Edges_Per_Face(); ++i ) {
                    ContactEdge<ActiveScalar>* edge = ContactLineEdgeL2<ActiveScalar>::new_ContactLineEdgeL2(
                                  active_allocators[ContactSearch::ALLOC_ContactLineEdgeL2],
                                  ContactSearch::LINEEDGEL2);
          
                    ContactNode<Real> *node[3];
                    master_face->Get_Edge_Nodes(i, node);
                    ContactNode<ActiveScalar> *active_node[3];
                    active_master_face->Get_Edge_Nodes(i, active_node);
                    if(master_face->Edge(i)->Node(0) == node[0]) {
                      for (int ii=0; ii<edge->Nodes_Per_Edge(); ++ii)
                        edge->ConnectNode(ii,active_node[ii]);
                    }
                    else { // some edges can be reversed
                      for (int ii=0; ii<edge->Nodes_Per_Edge(); ++ii)
                        edge->ConnectNode(edge->Nodes_Per_Edge()-ii-1,active_node[ii]);
                    }
                    active_master_face->ConnectEdge(i, edge);
                    edge->ConnectFace(active_master_face);
                  }
                  active_master_face->Compute_Centroid(POSITION, CENTROID);
                  active_master_face->Compute_Normal(POSITION, FACE_NORMAL);
                  active_master_face->Compute_CharacteristicLength(POSITION, CHARACTERISTIC_LENGTH);
          
                  // convert the master element
                  ContactElem<ActiveScalar> *active_master_element;
                  switch (master_element->Elem_Type()) {
                  case HEXELEML8:
                    active_master_element = ContactHexElemL8<ActiveScalar>::new_ContactHexElemL8(
                                                active_allocators[ALLOC_ContactHexElemL8],
                                                ContactSearch::HEXELEML8, 0);
                    break;
                  case WEDGEELEML6:
                    active_master_element = ContactWedgeElemL6<ActiveScalar>::new_ContactWedgeElemL6(
                                                active_allocators[ALLOC_ContactWedgeElemL6],
                                                ContactSearch::WEDGEELEML6, 0);
                    break;
                  default:
                    PRECONDITION(0);
                    break;
                  }
                  active_master_element->BuildTopology(number_of_nodes,
                                                       number_of_edges,
                                                       number_of_faces,
                                                       active_allocators);
                  active_master_element->UpdateTopology(active_master_face, POSITION,
                                                        FACE_NORMAL, NODE_NORMAL,
                                                        normal_search_tol);
          
                  // compute the interaction and derivatives
                  ContactFaceFaceInteraction<ActiveScalar> *active_cffi = Face_Face_Search(active_slave_face, active_master_face, active_master_element,
                                                                                           POSITION, active_allocators);
                  // copy the derivatives from active_cffi to cffi
                  if(active_cffi) {
                    cffi->Set_Derivatives(active_cffi); 

                    // delete active_cffi
                    active_cffi->~ContactFaceFaceInteraction<ActiveScalar>();
                    active_allocators[ALLOC_ContactFaceFaceInteraction].Delete_Frag(active_cffi);
                  }
                  else {
                    cffi->~ContactFaceFaceInteraction<Real>();
                    allocators[ALLOC_ContactFaceFaceInteraction].Delete_Frag(cffi);
                    cffi = 0;
                  }

                  // delete active_slave_face
                  for( int i=0; i<active_slave_face->Edges_Per_Face(); ++i ) { 
                    ContactEdge<ActiveScalar>* edge = active_slave_face->Edge(i);
                    edge->~ContactEdge<ActiveScalar>();
                    active_allocators[ALLOC_ContactLineEdgeL2].Delete_Frag(edge);
                  }
                  for( int i=0; i<active_slave_face->Nodes_Per_Face(); ++i ) { 
                    ContactNode<ActiveScalar>* node = active_slave_face->Node(i);
                    node->~ContactNode<ActiveScalar>();
                    active_allocators[ALLOC_ContactNode].Delete_Frag(node);
                  }
                  active_slave_face->~ContactFace<ActiveScalar>();
                  switch (slave_face->FaceType()) {
                  case TRIFACEL3:
                    active_allocators[ALLOC_ContactTriFaceL3].Delete_Frag(active_slave_face);
                    break;
                  case QUADFACEL4:
                    active_allocators[ALLOC_ContactQuadFaceL4].Delete_Frag(active_slave_face);
                    break;
                  case SHELLTRIFACEL3:
                    active_allocators[ALLOC_ContactShellTriFaceL3].Delete_Frag(active_slave_face);
                    break;
                  case SHELLQUADFACEL4:
                    active_allocators[ALLOC_ContactShellQuadFaceL4].Delete_Frag(active_slave_face);
                    break;
                  default:
                    PRECONDITION(0);
                    break;
                  }
          
                  // delete active_master_face
                  for( int i=0; i<active_master_face->Edges_Per_Face(); ++i ) {
                    ContactEdge<ActiveScalar>* edge = active_master_face->Edge(i);
                    edge->~ContactEdge<ActiveScalar>();
                    active_allocators[ALLOC_ContactLineEdgeL2].Delete_Frag(edge);
                  }
                  for( int i=0; i<active_master_face->Nodes_Per_Face(); ++i ) {
                    ContactNode<ActiveScalar>* node = active_master_face->Node(i);
                    node->~ContactNode<ActiveScalar>();
                    active_allocators[ALLOC_ContactNode].Delete_Frag(node);
                  }
                  active_master_face->~ContactFace<ActiveScalar>();
                  switch (master_face->FaceType()) {
                  case TRIFACEL3:
                    active_allocators[ALLOC_ContactTriFaceL3].Delete_Frag(active_master_face);
                    break;
                  case QUADFACEL4:
                    active_allocators[ALLOC_ContactQuadFaceL4].Delete_Frag(active_master_face);
                    break;
                  case SHELLTRIFACEL3:
                    active_allocators[ALLOC_ContactShellTriFaceL3].Delete_Frag(active_master_face);
                    break;
                  case SHELLQUADFACEL4:
                    active_allocators[ALLOC_ContactShellQuadFaceL4].Delete_Frag(active_master_face);
                    break;
                  default:
                    PRECONDITION(0);
                    break;
                  }
          
                  // delete active_master_element
                  active_master_element->DeleteTopology(active_allocators);
                  active_master_element->~ContactElem<ActiveScalar>();
                  switch (master_element->Elem_Type()) {
                  case HEXELEML8:
                    active_allocators[ALLOC_ContactHexElemL8].Delete_Frag(active_master_element);
                    break;
                  case WEDGEELEML6:
                    active_allocators[ALLOC_ContactWedgeElemL6].Delete_Frag(active_master_element);
                    break;
                  default:
                    PRECONDITION(0);
                    break;
                  }
#endif
                }
              } break;
              case 2 :
                // delete cffi and recompute interaction with 1st and 2nd order (hand-coded) partial derivatives
                cffi->~ContactFaceFaceInteraction<Real>();
                allocators[ALLOC_ContactFaceFaceInteraction].Delete_Frag(cffi);
                cffi = Second_Partial_Face_Face_Search(slave_face, master_face, master_element, POSITION, normal_search_tol, allocators);
                break;
            }
          }
          if(cffi) slave_face->Store_FaceFace_Interaction(cffi);
        }
      }
    }
  }
  delete [] face_boxes;
  delete [] face_boxes_tmp;
  if(face_hierarchy != NULL) delete [] face_hierarchy;
  delete [] list;

#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Stop_Timer( search_main_time );
#endif

#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Start_Timer( search_release_scratch_time );
#endif
  // Delete the old master element if it exists
  if (master_element!=NULL) {
    master_element->DeleteTopology(allocators);
    master_element->~ContactElem<Real>();
    switch (search_topology->Face_Block(old_block_id)->Type()) {
    case TRIFACEL3:
    case SHELLTRIFACEL3:
      allocators[ALLOC_ContactWedgeElemL6].Delete_Frag(master_element);
      break;
    case QUADFACEL4:
    case SHELLQUADFACEL4:
      allocators[ALLOC_ContactHexElemL8].Delete_Frag(master_element);
      break;
    default:
      PRECONDITION(0);
      break;  
    }
    master_element = NULL;
  }
#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Stop_Timer( search_release_scratch_time );
#endif
  
#if CONTACT_DEBUG_PRINT_LEVEL>=6
  postream << "\n\n\n***** Face-Face Interactions *****\n";
  search_topology->Display_FaceFace_Interactions(postream);
#endif
}

void
ContactSearch::Global_ElemElemSearch(SearchType search_type, int num_configs,
                                     VariableHandle CURRENT_POSITION, 
                                     VariableHandle PREDICTED_POSITION,
                                     VariableHandle AUGMENTED_POSITION )
{
#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Start_Timer( search_serial_setup_time );
#endif

  VariableHandle POSITION;
  switch (num_configs) {
  case 1:
    POSITION = CURRENT_POSITION;
    break;
  case 2:
    POSITION = PREDICTED_POSITION;
    break;
  case 3:
    POSITION = AUGMENTED_POSITION;
    break;
  default:
    POSTCONDITION(0);   
    break;   
  }
  
  int number_of_elements = search_topology->Number_of_Elements();
  ContactElement** Elements = 
    reinterpret_cast<ContactElement**>(search_topology->ElemList()->EntityList());
  
  // Process any element-element interactions
  int* list   = new int[number_of_elements];
  
#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Stop_Timer( search_serial_setup_time );
#endif

  ObjectBoundingBox *element_boxes = new ObjectBoundingBox[number_of_elements];
  int num_master_elements = 0;
  ContactTopologyEntityList* elem_list = search_topology->ElemList();
  for (int i=0; i<search_topology->Number_of_Element_Blocks(); ++i) {
    if (search_data->Is_ElemBlock_ElemElemMaster(i)) {
      int nelems = elem_list->BlockNumEntities(i);
      ContactElement** Elems = 
        reinterpret_cast<ContactElement**>(elem_list->BlockEntityList(i));
      for (int j=0; j<nelems; ++j){
        ContactElement* element = Elems[j];
        element_boxes[num_master_elements].set_object_number(element->ProcArrayIndex());
        for(int inode = 0; inode < element->Nodes_Per_Element(); ++inode) {
          Real *node_position = element->Node(inode)->Variable(POSITION);
          element_boxes[num_master_elements].add_point(node_position);
        }
#if CONTACT_DEBUG_PRINT_LEVEL>=10
        postream<<"Master element "<<num_master_elements<<", GID "<<element->Global_ID()<<"\n";
        postream<<"  min("<<element_boxes[num_master_elements].get_x_min()<<", "
                          <<element_boxes[num_master_elements].get_y_min()<<", "
                          <<element_boxes[num_master_elements].get_z_min()<<")\n";
        postream<<"  max("<<element_boxes[num_master_elements].get_x_max()<<", "
                          <<element_boxes[num_master_elements].get_y_max()<<", "
                          <<element_boxes[num_master_elements].get_z_max()<<")\n";
#endif
        element->temp_tag = num_master_elements++;
      }
    }
  }
#if CONTACT_DEBUG_PRINT_LEVEL>=3
  postream << "number of master elements = "<<num_master_elements<<"\n";
#endif

  int hierarchy_size = 2 * number_of_elements - 1;
  ObjectBoundingBoxHierarchy *element_hierarchy = NULL;
  if(hierarchy_size > 0) {
    element_hierarchy = new ObjectBoundingBoxHierarchy[hierarchy_size];
    ObjectBoundingBoxHierarchy::create_hierarchy(element_hierarchy, 
                                                 element_boxes, 
                                                 number_of_elements);
  }
  delete [] element_boxes;
  
#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Stop_Timer( search_serial_setup_time );
#endif

#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Start_Timer( search_main_time );
#endif
  // use
  Real thex[24][4][3];
  Real shex[24][4][3];
  int ntets = 0;
  if( no_warped_volume == ContactSearch::INACTIVE )
    ntets = 24;
  else
    ntets = 5;
  int requested_ntets = ntets;
  // Loop over all the elements
  for (int i=0; i<search_topology->Number_of_Element_Blocks(); ++i) {
    if (!search_data->Is_ElemBlock_ElemElemSlave(i)) continue;
    int nelems = elem_list->BlockNumEntities(i);
    ContactElement** Elems = 
      reinterpret_cast<ContactElement**>(elem_list->BlockEntityList(i));
    for (int j=0; j<nelems; ++j){
      ContactBoundingBox object_box;
      ContactElement* slave_element = Elems[j];
      if(slave_element->Ownership() != ContactTopologyEntity<Real>::OWNED) continue;
      int  slave_element_key = slave_element->Entity_Key();
      //
      // Determine all geometrically possible interactions between 
      // this element and all other elements in the mesh
      //
      for(int inode = 0; inode < slave_element->Nodes_Per_Element(); ++inode) {
        Real *node_position = slave_element->Node(inode)->Variable(POSITION);
        object_box.add_point(node_position);
      }
      int list_size = 0;
      ObjectBoundingBoxHierarchy::search_for_overlap_loop(element_hierarchy,object_box, list, list_size);
      if( list_size ){
        if (no_parallel_consistency==INACTIVE) Sort(list_size, list);
        for (int k=0; k<list_size; ++k) {
          ContactElement* master_element = Elements[list[k]];

          //
          //  Cull out the found interactions based on various criteria
          //
          //  Test 1: Cull the element if no interaction is requested between these entities
          //

          int  master_element_key = master_element->Entity_Key();
          int Interaction_Type = (int)
               search_data->Get_Search_Data( ContactSearch::INTERACTION_TYPE,
                                             slave_element_key, master_element_key );
          if( Interaction_Type == NO_INTERACTION ) continue;
          //
          //  Test 2: Cull out the element if it is the same as the master element
          //
          if (slave_element->ProcArrayIndex() == master_element->ProcArrayIndex()) continue;
          //
          //  Add the element-element interaction
          //
	  
          if (slave_element->ElementType() ==HEXELEMENTL8 &&
              master_element->ElementType()==HEXELEMENTL8) {
	    int tets24 = 24;
            master_element->TetDice(tets24,thex, POSITION);
            slave_element->TetDice(tets24,shex, POSITION);

	    Real volume = intersection_volume(thex,shex);
	    if (volume > 0.0) {
	      ContactElementElementInteraction* ceei =
                ContactElementElementInteraction::new_ContactElementElementInteraction(
		allocators[ALLOC_ContactElementElementInteraction],
	        slave_element, master_element, volume );
	      slave_element->Store_ElementElement_Interaction(ceei);
	    }
          }
	  else{
	    
	    ContactElement* cartesian_hex = NULL;
	    if (slave_element->ElementType() ==HEXELEMENTL8 &&
		master_element->ElementType()==CARTESIANHEXELEMENTL8) {
	      cartesian_hex = master_element;
	      slave_element->TetDice(ntets, thex, POSITION);
	    }
	    if (slave_element->ElementType() ==CARTESIANHEXELEMENTL8 &&
		master_element->ElementType()==HEXELEMENTL8) {
	      cartesian_hex = slave_element;
	      master_element->TetDice(ntets, thex, POSITION);
	    }
	    if (slave_element->ElementType() ==CARTESIANHEXELEMENTL8 &&
		master_element->ElementType()==CARTESIANHEXELEMENTL8) {
	      cartesian_hex = slave_element;
	      ntets = 5; // cartesian elements can't be warped so set to 5 will 
	      // later set back to requested_ntets (see below)
	      master_element->TetDice(ntets, thex, POSITION);
	    }
	   	    
	    int n;
	    Real mesh[3][2];
	    ContactNode<Real>* node = cartesian_hex->Node(0);
	    Real* xyz  = node->Variable(POSITION);
	    mesh[0][0] = mesh[0][1] = xyz[0];
	    mesh[1][0] = mesh[1][1] = xyz[1];
	    mesh[2][0] = mesh[2][1] = xyz[2];
	    for (n=1; n<cartesian_hex->Nodes_Per_Element(); ++n) {
	      node = cartesian_hex->Node(n);
	      xyz  = node->Variable(POSITION);
	      mesh[0][0] = std::min(mesh[0][0], xyz[0]);
	      mesh[0][1] = std::max(mesh[0][1], xyz[0]);
	      mesh[1][0] = std::min(mesh[1][0], xyz[1]);
	      mesh[1][1] = std::max(mesh[1][1], xyz[1]);
	      mesh[2][0] = std::min(mesh[2][0], xyz[2]);
	      mesh[2][1] = std::max(mesh[2][1], xyz[2]);
	    }
            
	    Real volume = 0.0;
	    for (n=0; n<ntets; ++n) {
	      Real areas[1];
	      Real vol = intersect_3d(&thex[n][0][0], &thex[n][1][0], 
				      &thex[n][2][0], &thex[n][3][0],
				      &mesh[0][0], &mesh[1][0], &mesh[2][0],
				      areas);
	      volume += vol;
	    }
	    ntets = requested_ntets; // if changed because both are cartesian then switch back
	    if (volume>0.0) {
	      ContactElementElementInteraction* ceei =
                ContactElementElementInteraction::new_ContactElementElementInteraction(
		allocators[ALLOC_ContactElementElementInteraction],
		slave_element, master_element, volume );
	      slave_element->Store_ElementElement_Interaction(ceei);
	    }
	  }
        }
      }
    }
  }
  if(element_hierarchy != NULL) delete [] element_hierarchy;
  delete [] list;

#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Stop_Timer( search_main_time );
#endif

#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Start_Timer( search_release_scratch_time );
#endif

#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Stop_Timer( search_release_scratch_time );
#endif
  
#if CONTACT_DEBUG_PRINT_LEVEL>=6
  postream << "\n\n\n***** Element-Element Interactions *****\n";
  search_topology->Display_ElementElement_Interactions(postream);
#endif
}

void ContactSearch::build_nodeface_search_status(int *node_search_status,
                                                 int *face_search_status) {
  int i, ii;
  
  ContactTopologyEntityList* node_list = search_topology->NodeList();
  ContactTopologyEntityList* face_list = search_topology->FaceList();
  
  int nnodes = node_list->BlockNumEntities(0);
  ContactNode<Real>** Nodes  = 
    reinterpret_cast<ContactNode<Real>**>(node_list->BlockEntityList(0));
  for (i=0; i<nnodes; ++i){
    ContactNode<Real>* node = Nodes[i];
    node->temp_tag = 0;
    //bool test2 = node->Ownership()==ContactTopologyEntity<Real>::OWNED;
    //postream<<"    1> checking node "<<node->Exodus_ID()<<"\n";
    //postream<<"         Owner         "<<node->Secondary_Owner()<<"\n";
    //postream<<"         Ownership     "<<test2<<"\n";
    //postream<<"         StatusFlag    "<<node->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)<<"\n";
    //postream<<"         Physical_Type "<<node->Physical_Type()<<"\n";
    if (!node->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)) continue;
    if (node->Ownership()     != ContactTopologyEntity<Real>::OWNED)           continue;
    if (node->Physical_Type() == ContactNode<Real>::SHELL_TAB_NODE)            continue;
    if (node->Physical_Type() == ContactNode<Real>::MIXED_TAB_NODE)            continue;
    node->temp_tag = 1;
    //postream<<"    1> tagging node "<<node->Exodus_ID()<<"\n";
  } 
  for (i=1; i<search_topology->Number_of_Node_Blocks(); ++i) {
    if (!search_data->Is_NodeBlock_NodeFaceSlave(i)) continue;
    nnodes = node_list->BlockNumEntities(i);
    Nodes  = reinterpret_cast<ContactNode<Real>**>(node_list->BlockEntityList(i));
    for (int j=0; j<nnodes; ++j){
      ContactNode<Real>* node = Nodes[j];
      node->temp_tag = 0;
      if (!node->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)) continue;
      if (node->Ownership()     != ContactTopologyEntity<Real>::OWNED)           continue;
      if (node->Physical_Type() == ContactNode<Real>::SHELL_TAB_NODE)            continue;
      if (node->Physical_Type() == ContactNode<Real>::MIXED_TAB_NODE)            continue;
      node->temp_tag = 2;
      //postream<<"    2> tagging node "<<node->Exodus_ID()<<"\n";
    } 
  }
     
  for (i=0; i<search_topology->Number_of_Face_Blocks(); ++i) {
    if (!search_data->Is_FaceBlock_NodeFaceSlave(i)) continue;
    int nfaces = face_list->BlockNumEntities(i);
    ContactFace<Real>** Faces  = 
      reinterpret_cast<ContactFace<Real>**>(face_list->BlockEntityList(i));
    for (int j=0; j<nfaces; ++j){
      ContactFace<Real>* face = Faces[j];
      //postream<<"    3a> tagging face "<<face->Global_ID()<<"\n";
      int num_nodes = face->Nodes_Per_Face();
      for (int k=0; k<num_nodes; ++k) {
        ContactNode<Real>* node = face->Node(k);
        if (node->temp_tag==1) {
        node->temp_tag=2;
        //postream<<"    3b> tagging node "<<node->Exodus_ID()<<"\n";
        }
      }
    }
  }
  
  num_active_nodes = 0;
  for (ii=0, i=0; i<search_topology->Number_of_Node_Blocks(); ++i) {
    nnodes = node_list->BlockNumEntities(i);
    Nodes  = reinterpret_cast<ContactNode<Real>**>(node_list->BlockEntityList(i));
    for (int j=0; j<nnodes; ++j){
      if (Nodes[j]->temp_tag==2) {
        node_search_status[ii++] = 1;
        //postream<<"    4> tagging node "<<Nodes[j]->Exodus_ID()<<"\n";
        ++num_active_nodes;
      } else {
        node_search_status[ii++] = 0;
      }
    } 
  }
  
  num_active_faces = 0;
  for (ii=0, i=0; i<search_topology->Number_of_Face_Blocks(); ++i) {
    int nfaces = face_list->BlockNumEntities(i);
    if (search_data->Is_FaceBlock_NodeFaceMaster(i)) {
      for (int j=0; j<nfaces; ++j){
        face_search_status[ii++] = 1;
        ++num_active_faces;
      }
    } else {
      for (int j=0; j<nfaces; ++j){
        face_search_status[ii++] = 0;
      }
    }
  }
}

void ContactSearch::build_nodeface_search_status_proximity(int *node_search_status,
                                                           int *face_search_status) {
  int i, ii;
  
  ContactTopologyEntityList* node_list = search_topology->NodeList();
  ContactTopologyEntityList* face_list = search_topology->FaceList();
  
  int nnodes = node_list->BlockNumEntities(0);
  ContactNode<Real>** Nodes  = 
    reinterpret_cast<ContactNode<Real>**>(node_list->BlockEntityList(0));
  for (i=0; i<nnodes; ++i){
    ContactNode<Real>* node = Nodes[i];
    node->temp_tag = 0;
    if (!node->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)) continue;
    if (node->Ownership()     != ContactTopologyEntity<Real>::OWNED)           continue;
    if (node->Physical_Type() == ContactNode<Real>::SHELL_TAB_NODE)            continue;
    if (node->Physical_Type() == ContactNode<Real>::MIXED_TAB_NODE)            continue;
    node->temp_tag = 1;
  } 
  for (i=1; i<search_topology->Number_of_Node_Blocks(); ++i) {
    if (!search_data->Is_NodeBlock_NodeFaceSlave(i)) continue;
    nnodes = node_list->BlockNumEntities(i);
    Nodes  = reinterpret_cast<ContactNode<Real>**>(node_list->BlockEntityList(i));
    for (int j=0; j<nnodes; ++j){
      ContactNode<Real>* node = Nodes[j];
      node->temp_tag = 0;
      if (!node->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)) continue;
      if (node->Ownership()     != ContactTopologyEntity<Real>::OWNED)           continue;
      if (node->Physical_Type() == ContactNode<Real>::SHELL_TAB_NODE)            continue;
      if (node->Physical_Type() == ContactNode<Real>::MIXED_TAB_NODE)            continue;
      node->temp_tag = 2;
    } 
  }
     
  for (i=0; i<search_topology->Number_of_Face_Blocks(); ++i) {
    if (!search_data->Is_FaceBlock_NodeFaceSlave(i)) continue;
    int nfaces = face_list->BlockNumEntities(i);
    ContactFace<Real>** Faces  = 
      reinterpret_cast<ContactFace<Real>**>(face_list->BlockEntityList(i));
    for (int j=0; j<nfaces; ++j){
      ContactFace<Real>* face = Faces[j];
      int num_nodes = face->Nodes_Per_Face();
      for (int k=0; k<num_nodes; ++k) {
        ContactNode<Real>* node = face->Node(k);
        if (node->temp_tag==1) node->temp_tag=2;
      }
    }
  }
  
  num_active_nodes = 0;
  for (ii=0, i=0; i<search_topology->Number_of_Node_Blocks(); ++i) {
    nnodes = node_list->BlockNumEntities(i);
    Nodes  = reinterpret_cast<ContactNode<Real>**>(node_list->BlockEntityList(i));
    for (int j=0; j<nnodes; ++j){
      if (Nodes[j]->temp_tag==2 && Nodes[j]->in_proximity) {
        node_search_status[ii++] = 1;
        ++num_active_nodes;
      } else {
        node_search_status[ii++] = 0;
      }
    } 
  }
    
  num_active_faces = 0;
  for (ii=0, i=0; i<search_topology->Number_of_Face_Blocks(); ++i) {
    int nfaces = face_list->BlockNumEntities(i);
    ContactFace<Real>** Faces = reinterpret_cast<ContactFace<Real>**>(face_list->BlockEntityList(i));
    if (search_data->Is_FaceBlock_NodeFaceMaster(i)) {
      for (int j=0; j<nfaces; ++j ){
         if (Faces[j]->in_proximity) {
           face_search_status[ii++] = 1;
            ++num_active_faces;
         } else {
           face_search_status[ii++] = 0;
         }
      }
    } else {
      for (int j=0; j<nfaces; ++j){
        face_search_status[ii++] = 0;
      }
    }
  }
}

//
//  Take a list of nodes and use it to compute a list of interactions
//
void ContactSearch::ComputeInteractions(const int list_size, 
                                        ACME::ContactNode_Vector &node_list,
                                        ContactFace<Real> *face,
                                        const VariableHandle &CURRENT_POSITION,
                                        const VariableHandle &AUGMENTED_POSITION,
                                        const VariableHandle &PREDICTED_POSITION,
                                        const VariableHandle &NODE_NORMAL,
                                        Real &user_search_tol,   // not used
                                        Real &user_tang_tol,     // CTOL2
                                        Real &relative_tang_tol, // CTOL1
                                        Real &gap_tol,
                                        int *physical_faces,
                                        const ACME::Int_Vector &node_keys) 
{
  int last_facet_type = -1;
  for(int k=0 ; k<list_size ; ++k ){
    ContactNode<Real> *node = node_list[k];

    int npairs = 1;
    int nfacets;


    ContactNodeFaceInteraction::InteractionSource Process_Method = 
      (ContactNodeFaceInteraction::InteractionSource)Dynamic_Process_Method( node, face );

    if(Process_Method==ContactNodeFaceInteraction::MOVING_INTERSECTION) {
      if(last_facet_type != 0) {
	last_facet_type = 0;
	face->FacetDecomposition(nfacets,
				 ms_coordinates_c, ms_normals_c, CURRENT_POSITION,
				 ms_coordinates_a, ms_normals_a, AUGMENTED_POSITION,
				 ms_coordinates_p, ms_normals_p, PREDICTED_POSITION );
      }
    } else {
      if(last_facet_type != 1){
	last_facet_type = 1;
	face->FacetDecomposition(nfacets,
				 ms_coordinates_a, ms_normals_a, AUGMENTED_POSITION,
				 ms_coordinates_p, ms_normals_p, PREDICTED_POSITION);
      }
    }

#ifdef CONTACT_DEBUG_NODE
    bool PRINT_THIS_NODE = primary_topology->Is_a_Debug_Node( node );
    if( PRINT_THIS_NODE ){
      postream << "  Processing Potential Interaction for Node ("
	       << node->Exodus_ID() << ") with Face "
	       << face->Global_ID() << "\n";
    }
#endif

    Real* sn_coordinates_c = node->Variable(CURRENT_POSITION);
    Real* sn_coordinates_a = node->Variable(AUGMENTED_POSITION);
    Real* sn_coordinates_p = node->Variable(PREDICTED_POSITION);
    //    Real* node_normal = node->Variable(NODE_NORMAL);
    for (int nn=0; nn<nfacets; ++nn) {
      pushback_dir_flag[nn] = 2; // normal pushback
    }
    for (int nn=0; nn<nfacets*data_size; ++nn) {
      ctrcl_facets[nn] = 0.0;
    }
    if(Process_Method==ContactNodeFaceInteraction::MOVING_INTERSECTION){
      if (dimensionality==3) {
	if (auto_tol) {
          for (int nn=0; nn<nfacets; ++nn) {
	    FORTRAN(cnodetri_movsrch_aug_noauto)( npairs,
						  sn_coordinates_c,
						  &ms_coordinates_c[nn*9], 
						  &ms_normals_c[nn*3],
						  sn_coordinates_p,
						  &ms_coordinates_p[nn*9], 
						  &ms_normals_p[nn*3],
						  sn_coordinates_a,
						  &ms_coordinates_a[nn*9], 
						  &ms_normals_a[nn*3],
						  &ctrcl_facets[nn*data_size],
						  gap_tol,
						  user_tang_tol );
	  }
	} else {
          for (int nn=0; nn<nfacets; ++nn) {
	    FORTRAN(cnodetri_movsrch_aug_noauto)( npairs,
						  sn_coordinates_c,
						  &ms_coordinates_c[nn*9], 
						  &ms_normals_c[nn*3],
						  sn_coordinates_p,
						  &ms_coordinates_p[nn*9], 
						  &ms_normals_p[nn*3],
						  sn_coordinates_a,
						  &ms_coordinates_a[nn*9], 
						  &ms_normals_a[nn*3],
						  &ctrcl_facets[nn*data_size],
						  gap_tol,
						  user_tang_tol );
	  }
	}
      }
      face->FacetDynamicRestriction(nfacets, ctrcl_facets, ctrcl);
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE && ctrcl[0] > 0 ){
	postream << " 	    Global Coords: (" 
		 << ctrcl[1] << "," << ctrcl[2] << "," 
		 << ctrcl[3] << ")\n";
	postream << " 	    Gap:	   " 
		 << ctrcl[4] << "\n";
	postream << " 	    Pushback Dir:  (" 
		 << ctrcl[5] << "," << ctrcl[6] << "," 
		 << ctrcl[7] << ")\n";
	postream << " 	    Normal Dir:    (" 
		 << ctrcl[8] << "," << ctrcl[9] << "," 
		 << ctrcl[10] << ")\n";
	postream << " 	    Location: "
		 << (int) ctrcl[11] << "\n";
	postream << " 	    Time - In/Out: "
		 << ctrcl[12] << "\n";
	postream << "  Interaction Being Sent to Process Interaction\n";
      }
#endif
    }
    // process interaction(s) flagged for CPP
    bool CPProj_Search = false;
    if( Process_Method == ContactNodeFaceInteraction::CLOSEST_POINT_PROJECTION_2 ) CPProj_Search = true;
    if( Process_Method == ContactNodeFaceInteraction::MOVING_INTERSECTION ){
      if( ctrcl[0] == -1 ){
	CPProj_Search = true;
	Process_Method = ContactNodeFaceInteraction::
	  CLOSEST_POINT_PROJECTION_2_FROM_MOVING;
      }
    }
    if ( CPProj_Search ) {
      if (auto_tol) {
      for (int nn=0; nn<nfacets; ++nn) {
	FORTRAN(cnodetriangle_cpproj_aug_noauto)(npairs,
					  sn_coordinates_p,
					  &ms_coordinates_p[nn*9], 
					  &ms_normals_p[nn*3],
					  sn_coordinates_a,
					  &ms_coordinates_a[nn*9], 
					  &ms_normals_a[nn*3],
					  &pushback_dir_flag[nn],
					  ctrl,
					  &ctrcl_facets[nn*data_size],
                                          relative_tang_tol,
                                          user_tang_tol,
                                          gap_tol);
#ifdef CONTACT_DEBUG_NODE
	if( PRINT_THIS_NODE ){
	  if( ctrcl_facets[nn*data_size] == 0 )
	    postream << "    Rejected with subtriangle " << nn << "\n";
	  else
	    postream << "    Accepted with subtriangle " << nn << "\n";
	}
#endif
      }
      } else {
      for (int nn=0; nn<nfacets; ++nn) {
	FORTRAN(cnodetriangle_cpproj_aug_noauto)(npairs,
					  sn_coordinates_p,
					  &ms_coordinates_p[nn*9], 
					  &ms_normals_p[nn*3],
					  sn_coordinates_a,
					  &ms_coordinates_a[nn*9], 
					  &ms_normals_a[nn*3],
					  &pushback_dir_flag[nn],
					  ctrl,
					  &ctrcl_facets[nn*data_size],
                                          relative_tang_tol,
                                          user_tang_tol,
                                          gap_tol);
      }
      }
      face->FacetStaticRestriction(nfacets, ms_coordinates_p, 
				   ms_normals_p, ctrcl_facets, ctrcl);
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE && ctrcl[0] > 0 ){
	postream << " 	    Global Coords: (" 
		 << ctrcl[1] << "," << ctrcl[2] << "," 
		 << ctrcl[3] << ")\n";
	postream << " 	    Gap:	   " 
		 << ctrcl[4] << "\n";
	postream << " 	    Pushback Dir:  (" 
		 << ctrcl[5] << "," << ctrcl[6] << "," 
		 << ctrcl[7] << ")\n";
	postream << " 	    Normal Dir:    (" 
		 << ctrcl[8] << "," << ctrcl[9] << "," 
		 << ctrcl[10] << ")\n";
	postream << " 	    Location: "
		 << (int) ctrcl[11] << "\n";
	postream << " 	    Time - In/Out: "
		 << ctrcl[12] << "\n";
	postream << "  Interaction Being Sent to Process Interaction\n";
      }
#endif
    }
    
    if( ctrcl[0] == 1 ){
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE )
	postream << "  Interaction Being Sent to Process Interaction\n";
#endif
      // The interaction is valid, find the physical face if there is one.
      // if there is not a physical face, set pf_norm to the node normal
      Real* pf_norm=NULL;
      switch( physical_faces[k] ){
      case 0:
	pf_norm = PHYSICAL_FACE_NORMAL_1.Get_Scratch(node->ProcArrayIndex());
	break;
      case 1:
	pf_norm = PHYSICAL_FACE_NORMAL_2.Get_Scratch(node->ProcArrayIndex());
	break;
      case 2:
	pf_norm = PHYSICAL_FACE_NORMAL_3.Get_Scratch(node->ProcArrayIndex());
	break;
      default:
        pf_norm = node->Variable(NODE_NORMAL);
	break;
      }
      bool is_tied=(search_data->Get_Search_Data(INTERACTION_TYPE, 
						 node_keys[k],
						 face->Entity_Key())
		    == TIED_INTERACTION ) ? true : false;
      bool is_infSlip=(search_data->Get_Search_Data(INTERACTION_TYPE,
						    node_keys[k],
						    face->Entity_Key())
		       == INFINITESIMAL_SLIP_INTERACTION ) ? true : false;
      bool is_tracked = (tracking_type==LOCAL_TRACKING);
      bool is_glued=(search_data->Get_Search_Data(INTERACTION_TYPE,
						    node_keys[k],
						    face->Entity_Key())
		       == GLUED_INTERACTION ) ? true : false;
      ContactNodeFaceInteraction* cnfi = 
	ContactNodeFaceInteraction::new_ContactNodeFaceInteraction( 
								   allocators[ALLOC_ContactNodeFaceInteraction],
								   Process_Method,node,face,ctrcl,node_keys[k],
								   pf_norm, is_tied, is_infSlip, is_glued, is_tracked,
                                                                   enable_off_face_tracking, PREDICTED_POSITION );
      // If no interaction exists for this node, store this one,
      // otherwise, either compete the interactions or accumulate 
      // this one
      int action = 0;
      Process_Interaction( PREDICTED_POSITION, cnfi, action );
    }
  }
}
