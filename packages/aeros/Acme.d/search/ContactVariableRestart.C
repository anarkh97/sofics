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


#include "ContactAnalyticPlane.h"
#include "ContactAnalyticSphere.h"
#include "ContactAnalyticCylinderInside.h"
#include "ContactAnalyticCylinderOutside.h"
#include "ContactCommBuffer.h"
#include "ContactErrors.h"
#include "ContactNodeBlock.h"
#include "ContactFaceBlock.h"
#include "ContactEdgeBlock.h"
#include "ContactNodeFaceInteraction.h"
#include "ContactNodeSurfaceInteraction.h"
#include "ContactFaceFaceInteraction.h"
#include "ContactFaceCoverageInteraction.h"
#include "ContactElementElementInteraction.h"
#include "ContactSearch.h"
#include "ContactSearchData.h"
#include "ContactTable.h"
#include "ContactTopology.h"
#include "ContactShellHandler.h"
#include "ContactShellNode.h"
#include "ContactSymComm.h"
#include "ContactFixedSizeAllocator.h"
#include <cstring>

#ifndef CONTACT_NO_MPI
#include "mpi.h"
#include "ContactZoltan.h"
#include "ContactZoltanComm.h"
#endif

using namespace std;

int ContactSearch::Number_General_Restart_Variables()
{
  return 25;
}

int ContactSearch::Number_Nodal_Restart_Variables()
{
  int num_vars = 0;
  // count variables for each interaction
  num_vars += 8;  // Face owner, block_id, host gid (2 ints),
                  // is tied flag, is inf slip flag, is_glued, is tracked flag;
  num_vars +=   ContactNodeFaceInteraction::NUMBER_SCALAR_VARS;
  num_vars += 3*ContactNodeFaceInteraction::NUMBER_VECTOR_VARS;
  // multiply # vars per interaction by max number of interactions
  num_vars *= Max_Interactions_Per_Node();
  // add non-interaction data
  num_vars += 6;  // remaining gap, node ghost gap vector variables
  // multiply all of the above by max number of shell nodes created
  // from a host code node
  if (primary_topology->Have_Shells()) {
    num_vars += 3; // previous lofting
    num_vars *=
      primary_topology->Shell_Handler()->Max_Num_Acme_Nodes_for_Host_Node();
  }

  return num_vars;
}

int ContactSearch::Number_Edge_Restart_Variables()
{
  return 0;
}

int ContactSearch::Number_Face_Restart_Variables()
{
  ContactTopologyEntity<Real>* entity;
  
  int num_ffi_l = primary_topology->Number_FaceFace_Interactions();
  int num_fci_l = primary_topology->Number_FaceCoverage_Interactions();

  int sum_input[2], sum_output[2];
  sum_input[0] = num_ffi_l;
  sum_input[1] = num_fci_l;
  contact_global_sum(sum_input, sum_output, 2, SearchComm);
  int num_ffi = sum_output[0];
  int num_fci = sum_output[1];

  int num_elem_vars = 0;
  if (num_ffi>0 || num_fci>0) {
    num_elem_vars = 2; // num_ffi and num_fci
    if (num_ffi>0) {
      int max_ffi_l = 0;
      int max_ffi_verts_l = 0;
      primary_topology->FaceList()->IteratorStart();
      while ((entity=primary_topology->FaceList()->IteratorForward())) {
        ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
        max_ffi_l = std::max(max_ffi_l,face->Number_FaceFace_Interactions());
        ContactInteractionDLL<Real>* interactions = face->Get_FaceFace_Interactions();
        if(interactions == NULL) continue;
        interactions->IteratorStart();
        while (ContactInteractionEntity<Real>* interaction=interactions->IteratorForward()){
          ContactFaceFaceInteraction<Real>* cffi = 
                 static_cast<ContactFaceFaceInteraction<Real>*> (interaction);
          max_ffi_verts_l = std::max(max_ffi_verts_l,cffi->NumEdges()+1);
        }
      }

      sum_input[0] = max_ffi_l;
      sum_input[1] = max_ffi_verts_l;
      contact_global_maximum(sum_input, sum_output, 2, SearchComm);
      int max_ffi       = sum_output[0];
      int max_ffi_verts = sum_output[1];
      num_elem_vars    += max_ffi*(3+6*max_ffi_verts);
    }
    if (num_fci>0) {
      int max_fci_l = 0;
      int max_fci_verts_l = 0;
      primary_topology->FaceList()->IteratorStart();
      while ((entity=primary_topology->FaceList()->IteratorForward())) {
        ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
        max_fci_l = std::max(max_fci_l,face->Number_FaceCoverage_Interactions());
        ContactInteractionDLL<Real>* interactions = face->Get_FaceCoverage_Interactions();
        if(interactions != NULL) {
          interactions->IteratorStart();
          while (ContactInteractionEntity<Real>* interaction=interactions->IteratorForward()){
            ContactFaceCoverageInteraction* cfci = 
                   static_cast<ContactFaceCoverageInteraction*> (interaction);
            max_fci_verts_l = std::max(max_fci_verts_l,cfci->NumVertices());
          }
	}
      }

      sum_input[0] = max_fci_l;
      sum_input[1] = max_fci_verts_l;
      contact_global_maximum(sum_input, sum_output, 2, SearchComm);
      int max_fci       = sum_output[0];
      int max_fci_verts = sum_output[1];
      num_elem_vars    += max_fci*(1+2*max_fci_verts);
    }
  }
  return num_elem_vars;
}

int ContactSearch::Number_Element_Restart_Variables()
{
  int num_elem_vars = 0;
  
  int num_eei_l = primary_topology->Number_ElementElement_Interactions();
  int num_eei   = contact_global_sum(num_eei_l, SearchComm );
  if (num_eei>0) {
    num_elem_vars = 1; // num_eei
    int max_eei_l = 0;
    ContactTopologyEntity<Real>* entity;
    primary_topology->ElemList()->IteratorStart();
    while ((entity=primary_topology->ElemList()->IteratorForward())) {
      ContactElement* element = static_cast<ContactElement*>(entity);
      max_eei_l = std::max(max_eei_l,
                      element->Number_ElementElement_Interactions());
    }
    int max_eei    = contact_global_maximum(max_eei_l, SearchComm );;
    num_elem_vars += max_eei*5;
  }
  return num_elem_vars;
}

ContactSearch::ContactErrorCode 
ContactSearch::Extract_General_Restart_Variable( Real* data )
{
  data[ 0] = multiple_interaction_status;
  data[ 1] = normal_smoothing_status;
  data[ 2] = compute_node_areas;
  data[ 3] = partition_gap_status;
  data[ 4] = old_dynamic_search;
  data[ 5] = tracking_type;
  data[ 6] = enable_off_face_tracking;
  data[ 7] = no_secondary;
  data[ 8] = search_cull;
  data[ 9] = no_warped_volume;
  data[10] = no_parallel_consistency;
  data[11] = orig_sharp_smooth_angle;
  data[12] = normal_smoothing_distance;
  data[13] = smoothing_resolution;
  data[14] = global_search_cull;
  data[15] = auto_tol;
  data[16] = aggressive_tolerances;
  data[17] = skip_physical_faces;
  data[18] = physical_face_algorithm;
  data[19] = box_inflation;
  data[20] = gap_inflation;
  data[21] = step_number;
  data[22] = tracking_step;
  data[23] = global_tracking_interval;
  data[24] = initialized_tied;
  return NO_ERROR;
}

ContactSearch::ContactErrorCode 
ContactSearch::Implant_General_Restart_Variable( Real* data )
{
  multiple_interaction_status = (Search_Option_Status) data[0];
  normal_smoothing_status     = (Search_Option_Status) data[1];
  compute_node_areas          = (Search_Option_Status) data[2];
  partition_gap_status        = (Search_Option_Status) data[3];
  old_dynamic_search          = (Search_Option_Status) data[4];
  tracking_type               = (Track_Type) data[5];
  enable_off_face_tracking    = (int) data[6];
  no_secondary                = (Search_Option_Status) data[7];
  search_cull                 = (Search_Option_Status) data[8];
  no_warped_volume            = (Search_Option_Status) data[9];
  no_parallel_consistency     = (Search_Option_Status) data[10];
  orig_sharp_smooth_angle          = data[11];
  sharp_smooth_curvature      = ComputeCurvatureFromAngle(orig_sharp_smooth_angle);
  normal_smoothing_distance   = data[12];
  smoothing_resolution        = (Smoothing_Resolution) data[13];
  global_search_cull          = (Search_Cull) data[14];
  auto_tol                    = (Search_Option_Status) data[15];
  aggressive_tolerances       = (Search_Option_Status) data[16];
  skip_physical_faces         = (Search_Option_Status) data[17];
  physical_face_algorithm     = (PF_Algorithm) data[18];
  box_inflation               = data[19];
  gap_inflation               = data[20];
  step_number                 = (int) data[21];
  tracking_step               = (int) data[22];
  global_tracking_interval    = (int) data[23];
  initialized_tied            = (bool) data[24];



  primary_topology->Set_Search_Options( multiple_interaction_status,
					normal_smoothing_status,
					sharp_smooth_curvature,
                                        no_parallel_consistency ,
                                        auto_tol);
  secondary_topology->Set_Search_Options( multiple_interaction_status,
					  normal_smoothing_status,
					  sharp_smooth_curvature,
                                          no_parallel_consistency ,
                                        auto_tol);
  restart = true;
  initialized = false;
  if (tracking_type==NO_TRACKING) {
    enable_tracking = INACTIVE;
  } else {
    enable_tracking = ACTIVE;
    if (tracking_type==GLOBAL_TRACKING) tracking_step=0;
  }
  return NO_ERROR;
}

ContactSearch::ContactErrorCode 
ContactSearch::Extract_Nodal_Restart_Variable( int n , Real* data, int *node_id )
{
  if( n>Number_Nodal_Restart_Variables() || n<=0  ){
    errors->Add_Error_Message("Variable Index is unreasonable");
    return INVALID_ID;
  }

  bool have_shells = primary_topology->Have_Shells();
  int nnod;
  if ( have_shells ) 
    nnod = primary_topology->Shell_Handler()->Number_Host_Code_Nodes();
  else 
    nnod = primary_topology->Number_of_Nodes();
  std::memset( data, 0, nnod*sizeof(Real) );

  VariableHandle REMAINING_GAP =
    primary_topology->Variable_Handle( ContactTopology::Remaining_Gap );
  VariableHandle NODE_GHOST_GAP =
    primary_topology->Variable_Handle( ContactTopology::Node_Ghost_Gap );

  int max_interactions = Max_Interactions_Per_Node();
  int vars_per_interaction = 8 +
    ContactNodeFaceInteraction::NUMBER_SCALAR_VARS +
    3*ContactNodeFaceInteraction::NUMBER_VECTOR_VARS;
  int vars_non_interaction = 6;
  int vars_add_shell = 3;
  int vars_per_shell_node = vars_per_interaction * max_interactions +
    vars_non_interaction + vars_add_shell;

  // here is where we figure out what we are writing for this requested
  // variable. For a pure hex mesh, the ordering is:
  //     info for interaction 1
  //     info for interaction 2
  //     info for interaction 3
  //     extra node info (non-interaction data)
  // For a mesh with shells, the ordering is:
  //     info for interaction 1 for first shell node from host node
  //     info for interaction 2 for first shell node from host node
  //     info for interaction 3 for first shell node from host node
  //     extra node info for first shell node from host node
  //     lofting for first shell node from host node
  //       ...
  //     info for interaction 1 for max shell node from host node
  //     info for interaction 2 for max shell node from host node
  //     info for interaction 3 for max shell node from host node
  //     extra node info for max shell node from host node
  //     lofting for max shell node from host node

  int host_to_acme_node = 0;
  int interaction_num = 0;
  int var_num = n - 1; // convert to C numbering
  bool output_interaction_data = false;
  bool output_lofting = false;
  bool output_non_interaction_data = false;

  // find which acme node (if we have shells)
  while (var_num >= vars_per_shell_node) {
    PRECONDITION(have_shells);
    host_to_acme_node ++;
    var_num -= vars_per_shell_node;
  }

  // find what variable to output
  if ( var_num < vars_per_interaction * max_interactions) {
    // this is interaction data
    // find which interaction
    output_interaction_data = true;
    while ( var_num >= vars_per_interaction) {
      interaction_num ++;
      var_num -= vars_per_interaction;
    }
  } else {
    var_num -= vars_per_interaction * max_interactions;
    if ( var_num < vars_non_interaction ) {
      // this is node data (non-interaction data)
      output_non_interaction_data = true;
    } else {
      // this is shell-specific data (lofting)
      PRECONDITION(have_shells);
      output_lofting = true;
      var_num -= vars_non_interaction;
      POSTCONDITION(var_num < vars_add_shell);
    }
  }

  //
  // loop over nodes to fill data array for this piece of data
  //
  for( int i=0 ; i<nnod ; ++i){
    ContactNode<Real>* node=NULL;

    // Compute ids of nodes to pass to host code and get pointer to node
    if (have_shells) {

      // get shell node base id from first related node
      int hi, lo;
      primary_topology->Shell_Handler()->
        Acme_NodeGID_for_Host_Node(i,0,hi,lo);
      ContactHostGlobalID first_gid(hi,lo);
      node = static_cast<ContactNode<Real>*>
        (primary_topology->NodeList()->Find(first_gid));
      ContactShellNode * shell_node =
        dynamic_cast<ContactShellNode*>(node);
      if (NULL != shell_node) {
        node_id[i]=shell_node->Shell_Node_Base_ID();
      } else {
        node_id[i]=node->Exodus_ID();
      }
      node = NULL;

      // We are looping over the max number of shell nodes related to
      // this host code node. If we don't have one for this host code
      // node, return 0.0 for the data.
      if ( primary_topology->Shell_Handler()->Num_Acme_Nodes_for_Host_Node(i)
           <=  host_to_acme_node) {
        data[i] = 0.0;
        continue;
      }

      // find the ACME node related to this host code node. If it is a hex,
      // there is only one; if it is a shell, there are at least two.
      primary_topology->Shell_Handler()->
        Acme_NodeGID_for_Host_Node(i,host_to_acme_node,hi,lo);
      ContactHostGlobalID gid(hi,lo);
      node = static_cast<ContactNode<Real>*>
        (primary_topology->NodeList()->Find(gid));

    } else {

      // getting the node is easy if you don't have shells
      node = static_cast<ContactNode<Real>*>
               (primary_topology->NodeList()->Find(i));
      node_id[i]=node->Exodus_ID();
    }



    POSTCONDITION(node);

    if ( output_lofting ) {
      // output lofting if that is what we need to do
      PRECONDITION (have_shells);
      PRECONDITION (var_num < 3);
      ContactShellNode * shell_node = NULL;
      shell_node = dynamic_cast<ContactShellNode*>(node);
      if ( shell_node ) {
        data[i] = shell_node->Previous_Lofting()[var_num];
      } else {
        data[i] = 0.0;
      }
    } else if ( output_non_interaction_data ) {
      // output non-interaction data
      PRECONDITION (var_num < vars_non_interaction);
      if (var_num == 0) {
        Real * remain_gap = node->Variable( REMAINING_GAP);
        data[i] = remain_gap[0];
      } else if (var_num == 1) {
        Real * remain_gap = node->Variable( REMAINING_GAP);
        data[i] = remain_gap[1];
      } else if (var_num == 2) {
        Real * remain_gap = node->Variable( REMAINING_GAP);
        data[i] = remain_gap[2];
      } else if (var_num == 3) {
        Real * ghost_gap = node->Variable( NODE_GHOST_GAP);
        data[i] = ghost_gap[0];
      } else if (var_num == 4) {
        Real * ghost_gap = node->Variable( NODE_GHOST_GAP);
        data[i] = ghost_gap[1];
      } else if (var_num == 5) {
        Real * ghost_gap = node->Variable( NODE_GHOST_GAP);
        data[i] = ghost_gap[2];
      }
    } else if ( output_interaction_data ) {
      // Get the correct interaction
      if (interaction_num<node->Number_NodeEntity_Interactions() &&
          node->Get_NodeEntity_Interaction( interaction_num )->Get_Type()==ContactNodeEntityInteraction::NODE_FACE_INTERACTION) {
        ContactNodeFaceInteraction* cnfi =
          static_cast<ContactNodeFaceInteraction*>(node->Get_NodeEntity_Interaction( interaction_num ));
        if( var_num==0 ){ // output the on-processor face id
          // Increment the on processor ID by one so 0 can be used to denote
          // no interaction.  We will have to std::remove the +1 upon reading the
          //data.
          data[i] = cnfi->FaceEntityData()->owner + 1;
        } else if( var_num==1 ){
          data[i] = cnfi->FaceEntityData()->block_id;
        } else if( var_num==2 ){
          data[i] = cnfi->FaceEntityData()->host_gid[0];
        } else if( var_num==3 ){
          data[i] = cnfi->FaceEntityData()->host_gid[1];
        } else if( var_num==4 ){ // output the is_tied flag
          if (cnfi->Is_Tied()) data[i] = 1;
          else data[i] = 0;
        } else if( var_num==5 ){ // output the is_infSlip flag
          if (cnfi->Is_InfSlip()) data[i] = 1;
          else data[i] = 0;
        } else if( var_num==6 ){ // output the is_glued flag
          if (cnfi->Is_Glued()) data[i] = 1;
          else data[i] = 0;
        } else if( var_num==7 ){ // output the is_tracked flag
          if (cnfi->Is_Tracked()) data[i] = 1;
          else data[i] = 0;
        } else { // output all variables stored on the interaction
          VariableHandle VH = (VariableHandle) var_num-8;
          data[i] = *(cnfi->Variable( VH ));
        }
      } // if on cnfi
    } // if on type of output variable
  } // end loop on nodes
  return NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactSearch::Implant_Nodal_Restart_Variable( int n, Real* data )
{
  if( n>Number_Nodal_Restart_Variables()  || n<=0 ){
    errors->Add_Error_Message("Variable Index is unreasonable");
    return INVALID_ID;
  }

  bool have_shells = primary_topology->Have_Shells();
  int nnod;
  if ( have_shells ) 
    nnod = primary_topology->Shell_Handler()->Number_Host_Code_Nodes();
  else  nnod = primary_topology->Number_of_Nodes();

  VariableHandle REMAINING_GAP =
    primary_topology->Variable_Handle( ContactTopology::Remaining_Gap );
  VariableHandle NODE_GHOST_GAP =
    primary_topology->Variable_Handle( ContactTopology::Node_Ghost_Gap );

  int max_interactions = Max_Interactions_Per_Node();
  int vars_per_interaction = 8 +
    ContactNodeFaceInteraction::NUMBER_SCALAR_VARS +
    3*ContactNodeFaceInteraction::NUMBER_VECTOR_VARS;
  int vars_non_interaction = 6;
  int vars_add_shell = 3;
  int vars_per_shell_node = vars_per_interaction * max_interactions +
    vars_non_interaction + vars_add_shell;

  // here is where we figure out what we are reading for this requested
  // variable. For a pure hex mesh, the ordering is:
  //     info for interaction 1
  //     info for interaction 2
  //     info for interaction 3
  //     extra node info (non-interaction data)
  // For a mesh with shells, the ordering is:
  //     info for interaction 1 for first shell node from host node
  //     info for interaction 2 for first shell node from host node
  //     info for interaction 3 for first shell node from host node
  //     extra node info for first shell node from host node
  //     lofting for first shell node from host node
  //       ...
  //     info for interaction 1 for max shell node from host node
  //     info for interaction 2 for max shell node from host node
  //     info for interaction 3 for max shell node from host node
  //     extra node info for max shell node from host node
  //     lofting for max shell node from host node

  int host_to_acme_node = 0;
  int interaction_num = 0;
  int var_num = n - 1; // convert to C numbering
  bool input_interaction_data = false;
  bool input_lofting = false;
  bool input_non_interaction_data = false;

  // find which acme node (if we have shells)
  while (var_num >= vars_per_shell_node) {
    PRECONDITION(have_shells);
    host_to_acme_node ++;
    var_num -= vars_per_shell_node;
  }

  // find what variable to output
  if ( var_num < vars_per_interaction * max_interactions) {
    // this is interaction data
    // find which interaction
    input_interaction_data = true;
    while ( var_num >= vars_per_interaction) {
      interaction_num ++;
      var_num -= vars_per_interaction;
    }
  } else {
    var_num -= vars_per_interaction * max_interactions;
    if ( var_num < vars_non_interaction ) {
      // this is node data (non-interaction data)
      input_non_interaction_data = true;
    } else {
      // this is shell-specific data (lofting)
      PRECONDITION(have_shells);
      input_lofting = true;
      var_num -= vars_non_interaction;
      POSTCONDITION(var_num < vars_add_shell);
    }
  }

  //
  // loop over nodes to fill data array for this piece of data
  //
  for( int i=0 ; i<nnod ; ++i){
    ContactNode<Real>* node=NULL;
    int host_node_index = -1;

    if (have_shells) {
      if ( primary_topology->Shell_Handler()->Num_Acme_Nodes_for_Host_Node(i)
           <=  host_to_acme_node) {
        data[i] = 0.0;
        continue;
      }
      int hi, lo;
      primary_topology->Shell_Handler()->
        Acme_NodeGID_for_Host_Node(i,host_to_acme_node,hi,lo);
      ContactHostGlobalID gid(hi,lo);
      node = static_cast<ContactNode<Real>*>(primary_topology->NodeList()->Find(gid));

      // also get base node -- it is one that points to proper location
      // in host code indexing
      ContactShellNode * shell_node = dynamic_cast<ContactShellNode*>(node);
      if (shell_node) {
        host_node_index =
          primary_topology->Shell_Handler()->Acme_to_Host_Node_Map()[node->HostGlobalArrayIndex()];
      } else {
        host_node_index = node->HostGlobalArrayIndex();
      }
    } else {
      node = static_cast<ContactNode<Real>*>(primary_topology->NodeList()->Find(i));
      host_node_index = node->HostGlobalArrayIndex();
    }
    POSTCONDITION(node);



    if ( input_lofting ) {
      // input lofting for shell node if that is what we need to do
      PRECONDITION (have_shells);
      PRECONDITION (var_num < 3);
      ContactShellNode * shell_node = NULL;
      shell_node = dynamic_cast<ContactShellNode*>(node);
      if ( shell_node ) {
        shell_node->Previous_Lofting()[var_num] = data[host_node_index];
      }
    } else if ( input_non_interaction_data ) {
      // input non-interaction data
      PRECONDITION (var_num < vars_non_interaction);
      if (var_num == 0) {
        Real * remain_gap = node->Variable( REMAINING_GAP);
        remain_gap[0] = data[i];
      } else if (var_num == 1) {
        Real * remain_gap = node->Variable( REMAINING_GAP);
        remain_gap[1] = data[i];
      } else if (var_num == 2) {
        Real * remain_gap = node->Variable( REMAINING_GAP);
        remain_gap[2] = data[i];
      } else if (var_num == 3) {
        Real * ghost_gap = node->Variable( NODE_GHOST_GAP);
        ghost_gap[0] = data[i];
      } else if (var_num == 4) {
        Real * ghost_gap = node->Variable( NODE_GHOST_GAP);
        ghost_gap[1] = data[i];
      } else if (var_num == 5) {
        Real * ghost_gap = node->Variable( NODE_GHOST_GAP);
        ghost_gap[2] = data[i];
      }
    } else if ( input_interaction_data ) {
      // input data for an interaction
      ContactNodeFaceInteraction* cnfi;
      if( var_num == 0 ){
        if( data[host_node_index] != 0 ){
          cnfi = ContactNodeFaceInteraction::new_ContactNodeFaceInteraction(Get_Allocator( ALLOC_ContactNodeFaceInteraction ),node );

          // Subtract 1 to account for the +1 in the Extract function.
          cnfi->FaceEntityData()->owner = ((int) data[host_node_index]-1);
          PRECONDITION(!node->Get_NodeEntity_Interaction( interaction_num ));
          node->Store_NodeEntity_Interaction( interaction_num, cnfi );
        }
      } else if( var_num == 1 ){
        if (interaction_num<node->Number_NodeEntity_Interactions() &&
            node->Get_NodeEntity_Interaction( interaction_num )->Get_Type()==ContactNodeEntityInteraction::NODE_FACE_INTERACTION) {
          cnfi = static_cast<ContactNodeFaceInteraction*>(node->Get_NodeEntity_Interaction( interaction_num ));
          PRECONDITION( cnfi != NULL || data[host_node_index] == 0.0 );
          cnfi->FaceEntityData()->block_id = ((int) data[host_node_index]);
        }
      } else if( var_num == 2 ){
        if (interaction_num<node->Number_NodeEntity_Interactions() &&
            node->Get_NodeEntity_Interaction( interaction_num )->Get_Type()==ContactNodeEntityInteraction::NODE_FACE_INTERACTION) {
          cnfi = static_cast<ContactNodeFaceInteraction*>(node->Get_NodeEntity_Interaction( interaction_num ));
          PRECONDITION( cnfi != NULL || data[host_node_index] == 0.0 );
          cnfi->FaceEntityData()->host_gid[0] = ((int) data[host_node_index]);
        }
      } else if( var_num == 3 ){
        if (interaction_num<node->Number_NodeEntity_Interactions() &&
            node->Get_NodeEntity_Interaction( interaction_num )->Get_Type()==ContactNodeEntityInteraction::NODE_FACE_INTERACTION) {
          cnfi = static_cast<ContactNodeFaceInteraction*>(node->Get_NodeEntity_Interaction( interaction_num ));
          PRECONDITION( cnfi != NULL || data[host_node_index] == 0.0 );
          cnfi->FaceEntityData()->host_gid[1] = ((int) data[host_node_index]);
          cnfi->Connect_Face( *primary_topology->FaceList() );
        }
      } else if( var_num == 4 ){
        if (interaction_num<node->Number_NodeEntity_Interactions() &&
            node->Get_NodeEntity_Interaction( interaction_num )->Get_Type()==ContactNodeEntityInteraction::NODE_FACE_INTERACTION) {
          cnfi = static_cast<ContactNodeFaceInteraction*>(node->Get_NodeEntity_Interaction( interaction_num ));
          PRECONDITION( cnfi != NULL || data[host_node_index] == 0.0 );
          if ( data[host_node_index] == 1 ) cnfi->Is_Tied(true);
          else cnfi->Is_Tied(false);
        }
      } else if( var_num == 5 ){
        if (interaction_num<node->Number_NodeEntity_Interactions() &&
            node->Get_NodeEntity_Interaction( interaction_num )->Get_Type()==ContactNodeEntityInteraction::NODE_FACE_INTERACTION) {
          cnfi = static_cast<ContactNodeFaceInteraction*>(node->Get_NodeEntity_Interaction( interaction_num ));
          PRECONDITION( cnfi != NULL || data[host_node_index] == 0.0 );
          if ( data[host_node_index] == 1 ) cnfi->Is_InfSlip(true);
          else cnfi->Is_InfSlip(false);
        }
      } else if( var_num == 6 ){
        if (interaction_num<node->Number_NodeEntity_Interactions() &&
            node->Get_NodeEntity_Interaction( interaction_num )->Get_Type()==ContactNodeEntityInteraction::NODE_FACE_INTERACTION) {
          cnfi = static_cast<ContactNodeFaceInteraction*>(node->Get_NodeEntity_Interaction( interaction_num ));
          PRECONDITION( cnfi != NULL || data[host_node_index] == 0.0 );
          if ( data[host_node_index] == 1 ) cnfi->Is_Glued(true);
          else cnfi->Is_Glued(false);
        }
      } else if( var_num == 7 ){
        if (interaction_num<node->Number_NodeEntity_Interactions() &&
            node->Get_NodeEntity_Interaction( interaction_num )->Get_Type()==ContactNodeEntityInteraction::NODE_FACE_INTERACTION) {
          cnfi = static_cast<ContactNodeFaceInteraction*>(node->Get_NodeEntity_Interaction( interaction_num ));
          PRECONDITION( cnfi != NULL || data[host_node_index] == 0.0 );
          if ( data[host_node_index] == 1 ) cnfi->Is_Tracked(true);
          else cnfi->Is_Tracked(false);
        }
      } else if( var_num <  vars_per_interaction ){
        if (interaction_num<node->Number_NodeEntity_Interactions() &&
            node->Get_NodeEntity_Interaction( interaction_num )->Get_Type()==ContactNodeEntityInteraction::NODE_FACE_INTERACTION) {
          cnfi = static_cast<ContactNodeFaceInteraction*>(node->Get_NodeEntity_Interaction( interaction_num ));
          PRECONDITION( cnfi != NULL || data[host_node_index] == 0.0 );
          VariableHandle VH = (VariableHandle) var_num-8;
          *cnfi->Variable( VH ) = data[host_node_index];
        }
      } // end switch on var_num
    } // end if on type of input variable
  } // end loop on nodes
  restart = true;
  return NO_ERROR;
}

ContactSearch::ContactErrorCode 
ContactSearch::Extract_Edge_Restart_Variable( int n , Real* data )
{
  if( n>Number_Edge_Restart_Variables()  || n<=0 ){
    errors->Add_Error_Message("Variable Index is unreasonable");
    return INVALID_ID;
  }
  return NO_ERROR;
}

ContactSearch::ContactErrorCode 
ContactSearch::Implant_Edge_Restart_Variable( int n , Real* data )
{
  if( n>Number_Edge_Restart_Variables()  || n<=0 ){
    errors->Add_Error_Message("Variable Index is unreasonable");
    return INVALID_ID;
  }
  restart = true;
  return NO_ERROR;
}

ContactSearch::ContactErrorCode 
ContactSearch::Extract_Face_Restart_Variable( int n , Real* data )
{
  ContactTopologyEntity<Real>* entity;
  ContactInteractionEntity<Real>* interaction;
  if( n>Number_Face_Restart_Variables()  || n<=0 ){
    errors->Add_Error_Message("Variable Index is unreasonable");
    return INVALID_ID;
  }
  // Initialize the data
  int number_of_faces = primary_topology->Number_of_Faces();
  std::memset( data, 0, number_of_faces*sizeof(Real) );

  int i;
  if (n==1) {
    i = 0;
    primary_topology->FaceList()->IteratorStart();
    while ((entity=primary_topology->FaceList()->IteratorForward())) {
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
      data[i++] = face->Number_FaceFace_Interactions();
    }
  } else if (n==2) {
    i = 0;
    primary_topology->FaceList()->IteratorStart();
    while ((entity=primary_topology->FaceList()->IteratorForward())) {
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
      data[i++] = face->Number_FaceCoverage_Interactions();
    }
  } else {
    int N=n-1-2;
    int max_ffi_l = 0;
    int max_ffi_verts_l = 0;
    primary_topology->FaceList()->IteratorStart();
    while ((entity=primary_topology->FaceList()->IteratorForward())) {
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
      max_ffi_l = std::max(max_ffi_l,face->Number_FaceFace_Interactions());
      ContactInteractionDLL<Real>* interactions = face->Get_FaceFace_Interactions();
      if(interactions == NULL) continue;
      interactions->IteratorStart();
      while ((interaction=interactions->IteratorForward())){
        ContactFaceFaceInteraction<Real>* cffi = 
          static_cast<ContactFaceFaceInteraction<Real>*>(interaction);
        max_ffi_verts_l = std::max(max_ffi_verts_l,cffi->NumEdges()+1);
      }
    }
    int sum_input[2], sum_output[2];
    sum_input[0] = max_ffi_l;
    sum_input[1] = max_ffi_verts_l;
    contact_global_maximum(sum_input, sum_output, 2, SearchComm);
    int max_ffi       = sum_output[0];
    int max_ffi_verts = sum_output[1];
    
    int max_fci_verts_l = 0;
    primary_topology->FaceList()->IteratorStart();
    while ((entity=primary_topology->FaceList()->IteratorForward())) {
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
      ContactInteractionDLL<Real>* interactions = face->Get_FaceCoverage_Interactions();
      if(interactions != NULL) {
        interactions->IteratorStart();
        while ((interaction=interactions->IteratorForward())){
          ContactFaceCoverageInteraction* cfci = 
            static_cast<ContactFaceCoverageInteraction*>(interaction);
          max_fci_verts_l = std::max(max_fci_verts_l,cfci->NumVertices());
        }
      }
    }
    int max_fci_verts = contact_global_maximum(max_fci_verts_l, SearchComm );
    
    int ffi_num = -1;
    int fci_num = -1;
    
    if (max_ffi>0 && N<max_ffi*(3+6*max_ffi_verts)) {
      int offset = N;
      ffi_num    = offset/(3+6*max_ffi_verts);
    } else {
      int offset = N-max_ffi*(3+6*max_ffi_verts);
      fci_num    = offset/(1+2*max_fci_verts);
    }
  
    if (ffi_num>=0) {
      i = 0;
      ContactFaceFaceInteraction<Real>* cffi;
      primary_topology->FaceList()->IteratorStart();
      while ((entity=primary_topology->FaceList()->IteratorForward())) {
        ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
        int cnt=0;
        ContactInteractionDLL<Real>* interactions = face->Get_FaceFace_Interactions();
        if(interactions != NULL) {
          interactions->IteratorStart();
          while ((interaction=interactions->IteratorForward())){
            if (cnt == ffi_num) break;
            cnt++;
          }
          if (interaction) {
            cffi = static_cast<ContactFaceFaceInteraction<Real>*> (interaction);
            int offset = N-ffi_num*(3+6*max_ffi_verts);
            if (offset==0) {
              data[i] = cffi->MasterFaceEntityData()->host_gid[0];
            } else if (offset==1) {
              data[i] = cffi->MasterFaceEntityData()->host_gid[1];
            } else if (offset==2) {
              data[i] = cffi->NumEdges()+1;
            } else {
              offset -= 3;
              int ffi_vertex = offset/6;
              int ffi_data   = offset-6*ffi_vertex;
              ContactFaceFaceVertex<Real>* vertices = cffi->Get_Vertices();
              switch (ffi_data) {
              case 0:
                data[i] = vertices[ffi_vertex].slave_x;
                break;
              case 1:
                data[i] = vertices[ffi_vertex].slave_y;
                break;
              case 2:
                data[i] = vertices[ffi_vertex].master_x;
                break;
              case 3:
                data[i] = vertices[ffi_vertex].master_y;
                break;
              case 4:
                data[i] = vertices[ffi_vertex].slave_edge_id;
                break;
              case 5:
                data[i] = vertices[ffi_vertex].master_edge_flag;
                break;
              }
            }
          }
	}
        i++;
      }
    }
    if (fci_num>=0) {
      i = 0;
      ContactFaceCoverageInteraction* cfci;
      primary_topology->FaceList()->IteratorStart();
      while ((entity=primary_topology->FaceList()->IteratorForward())) {
        ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
        int cnt=0;
        ContactInteractionDLL<Real>* interactions = face->Get_FaceCoverage_Interactions();
        if(interactions != NULL) {
          interactions->IteratorStart();
          while ((interaction=interactions->IteratorForward())){
            if (cnt == fci_num) break;
            cnt++;
          }
          if (interaction) {
            cfci = static_cast<ContactFaceCoverageInteraction*> (interaction);
            int offset = N-max_ffi*(3+6*max_ffi_verts)-fci_num*(1+2*max_fci_verts);
            if (offset==0) {
              data[i] = cfci->NumVertices();
            } else {
              offset--;
              int fci_vertex = offset/2;
              int fci_data   = offset-2*fci_vertex;
              ContactFaceCoverageVertex* ll_node;
              int kk=0;
              for( ll_node=cfci->Head(); ll_node; 
                   ll_node=ll_node->next, kk++ ){
                if (kk==fci_vertex) break;
              }
              switch (fci_data) {
              case 0:
                data[i] = ll_node->slave_x;
                break;
              case 1:
                data[i] = ll_node->slave_y;
                break;
              }
            }
          }
	}
        i++;
      }
    }
  }
  return NO_ERROR;
}

ContactSearch::ContactErrorCode 
ContactSearch::Implant_Face_Restart_Variable( int n , Real* data )
{
  ContactTopologyEntity<Real>* entity;
  if( n>Number_Face_Restart_Variables()  || n<=0 ){
    errors->Add_Error_Message("Variable Index is unreasonable");
    return INVALID_ID;
  }

  ContactFaceFaceInteraction<Real>* cffi;
  ContactFaceCoverageInteraction* cfci;
  int i;
  if (n==1) {
    i = 0;
    primary_topology->FaceList()->IteratorStart();
    while ((entity=primary_topology->FaceList()->IteratorForward())) {
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
      int num_faceface_interactions = (int) data[i++];
      for( int j=0 ; j<num_faceface_interactions ; ++j){
        cffi = ContactFaceFaceInteraction<Real>::new_ContactFaceFaceInteraction(
                allocators[ALLOC_ContactFaceFaceInteraction]);
        face->Store_FaceFace_Interaction( cffi );
      }
    }
  } else if (n==2) {
    i = 0;
    primary_topology->FaceList()->IteratorStart();
    while ((entity=primary_topology->FaceList()->IteratorForward())) {
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
      int num_facecoverage_interactions = (int) data[i++];
      for( int j=0 ; j<num_facecoverage_interactions ; ++j){
        cfci = ContactFaceCoverageInteraction::new_ContactFaceCoverageInteraction(
                allocators[ALLOC_ContactFaceCoverageInteraction]);
        face->Store_FaceCoverage_Interaction( cfci );
      }
    }
  } else {
    int N=n-1-2;
    int max_ffi_l = 0;
    int max_ffi_verts_l = 0;
    primary_topology->FaceList()->IteratorStart();
    while ((entity=primary_topology->FaceList()->IteratorForward())) {
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
      max_ffi_l = std::max(max_ffi_l, face->Number_FaceFace_Interactions());
      ContactInteractionDLL<Real>* interactions = face->Get_FaceFace_Interactions();
      if(interactions == NULL) continue;
      interactions->IteratorStart();
      while (ContactInteractionEntity<Real>* interaction=interactions->IteratorForward()){
        cffi = static_cast<ContactFaceFaceInteraction<Real>*> (interaction);
        max_ffi_verts_l = std::max(max_ffi_verts_l,cffi->NumEdges()+1);
      }
    }
    int sum_input[2], sum_output[2];
    sum_input[0] = max_ffi_l;
    sum_input[1] = max_ffi_verts_l;
    contact_global_maximum(sum_input, sum_output, 2, SearchComm);
    int max_ffi       = sum_output[0];
    int max_ffi_verts = sum_output[1];

    int max_fci_verts_l = 0;
    primary_topology->FaceList()->IteratorStart();
    while ((entity=primary_topology->FaceList()->IteratorForward())) {
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
      ContactInteractionDLL<Real>* interactions = face->Get_FaceCoverage_Interactions();
      if(interactions != NULL) {
        interactions->IteratorStart();
        while (ContactInteractionEntity<Real>* interaction=interactions->IteratorForward()){
          cfci = static_cast<ContactFaceCoverageInteraction*> (interaction);
          max_fci_verts_l = std::max(max_fci_verts_l,cfci->NumVertices());
        }
      }
    }
    int max_fci_verts = contact_global_maximum(max_fci_verts_l, SearchComm );
    
    int ffi_num  = -1;
    int fci_num  = -1;
    
    if (max_ffi>0 && N<max_ffi*(3+6*max_ffi_verts)) {
      int offset = N;
      ffi_num    = offset/(3+6*max_ffi_verts);
    } else {
      int offset = N-max_ffi*(3+6*max_ffi_verts);
      fci_num    = offset/(1+2*max_fci_verts);
    }
  
    if (ffi_num>=0) {
      i = 0;
      primary_topology->FaceList()->IteratorStart();
      while ((entity=primary_topology->FaceList()->IteratorForward())) {
        ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
        int cnt=0;
        ContactInteractionEntity<Real>* interaction = NULL;
        ContactInteractionDLL<Real>* interactions = face->Get_FaceFace_Interactions();
        if(interactions != NULL) {
          interactions->IteratorStart();
          while ((interaction=interactions->IteratorForward())){
            if (cnt == ffi_num) break;
            cnt++;
          }
          if (interaction) {
            cffi = static_cast<ContactFaceFaceInteraction<Real>*> (interaction);
            int offset = N-ffi_num*(3+6*max_ffi_verts);
            if (offset==0) {
              cffi->MasterFaceEntityData()->host_gid[0] = (int) data[i];
            } else if (offset==1) {
              cffi->MasterFaceEntityData()->host_gid[1] = (int) data[i];
            } else if (offset==2) {
              cffi->NumEdges((int) data[i]-1);
            } else {
              offset -= 3;
              int ffi_vertex = offset/6;
              int ffi_data   = offset-6*ffi_vertex;
              ContactFaceFaceVertex<Real>* vertices = cffi->Get_Vertices();
              switch (ffi_data) {
              case 0:
                vertices[ffi_vertex].slave_x = data[i];
                break;
              case 1:
                vertices[ffi_vertex].slave_y = data[i];
                break;
              case 2:
                vertices[ffi_vertex].master_x = data[i];
                break;
              case 3:
                vertices[ffi_vertex].master_y = data[i];
                break;
              case 4:
                vertices[ffi_vertex].slave_edge_id = (int)data[i];
                break;
              case 5:
                vertices[ffi_vertex].master_edge_flag = (int)data[i];
                break;
              }
            }
          }
	}
        i++;
      }
    }
    if (fci_num>=0) {
      i = 0;
      primary_topology->FaceList()->IteratorStart();
      while ((entity=primary_topology->FaceList()->IteratorForward())) {
        ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
        int cnt = 0;
        ContactInteractionEntity<Real>* interaction = NULL;
        ContactInteractionDLL<Real>* interactions = face->Get_FaceCoverage_Interactions();
        if(interactions != NULL) {
          interactions->IteratorStart();
          while ((interaction=interactions->IteratorForward())){
            if (!interaction) break;
            cnt++;
          }
          if (interaction) {
            cfci = static_cast<ContactFaceCoverageInteraction*> (interaction);
            int offset = N-max_ffi*(3+6*max_ffi_verts)-fci_num*(1+2*max_fci_verts);
            if (offset>0) {
              offset--;
              int fci_vertex = offset/2;
              int fci_data   = offset-2*fci_vertex;
              while (fci_vertex>=cfci->NumVertices()) {
                cfci->AddVertex(0.0, 0.0);
              }
              ContactFaceCoverageVertex* ll_node=cfci->Head();
              for (int kk=0; kk<fci_vertex; kk++, ll_node=ll_node->next) {
                if (!interaction) break;
              }
              if (ll_node) {
                switch (fci_data) {
                case 0:
                  ll_node->slave_x = data[i];
                  break;
                case 1:
                  ll_node->slave_y = data[i];
                  break;
                }
              }
            }
          }
	}
        i++;
      }
    }
  }
  
  restart = true;
  return NO_ERROR;
}

ContactSearch::ContactErrorCode 
ContactSearch::Extract_Element_Restart_Variable( int n , Real* data )
{
  ContactTopologyEntity<Real>* entity;
  ContactInteractionEntity<Real>* interaction;
  if( n>Number_Element_Restart_Variables() || n<=0 ){
    errors->Add_Error_Message("Variable Index is unreasonable");
    return INVALID_ID;
  }
  
  // Initialize the data
  int number_of_elements = primary_topology->Number_of_Elements();
  std::memset( data, 0, number_of_elements*sizeof(Real) );
  
  int i;
  if (n==1) {
    i = 0;
    primary_topology->ElemList()->IteratorStart();
    while ((entity=primary_topology->ElemList()->IteratorForward())) {
      ContactElement* element = static_cast<ContactElement*>(entity);
      data[i++] = element->Number_ElementElement_Interactions();
    }
  } else {
    i = 0;
    int eei_num = (n-2)/5;
    int offset  = (n-2)-eei_num*5;
    ContactElementElementInteraction* ceei;
    primary_topology->ElemList()->IteratorStart();
    while ((entity=primary_topology->ElemList()->IteratorForward())) {
      ContactElement* element = static_cast<ContactElement*>(entity);
      if (element->Number_ElementElement_Interactions()>eei_num) {
        int cnt=0;
        ContactInteractionDLL<Real>* interactions = element->Get_ElementElement_Interactions();
        interactions->IteratorStart();
        while ((interaction=interactions->IteratorForward())){
          if (cnt == eei_num) {
            ceei = static_cast<ContactElementElementInteraction*> (interaction);
            switch (offset) {
            case 0:
              data[i] = ceei->SlaveElementEntityData()->host_gid[0];
              break;
            case 1:
              data[i] = ceei->SlaveElementEntityData()->host_gid[1];
              break;
            case 2:
              data[i] = ceei->MasterElementEntityData()->host_gid[0];
              break;
            case 3:
              data[i] = ceei->MasterElementEntityData()->host_gid[1];
              break;
            case 4:
              data[i] = ceei->Scalar_Var(ContactElementElementInteraction::VOLUME);
              break;
            }
            break;
          }
          cnt++;
        }
      }
      i++;
    }
  }
  return NO_ERROR;
}

ContactSearch::ContactErrorCode 
ContactSearch::Implant_Element_Restart_Variable( int n , Real* data )
{
  ContactTopologyEntity<Real>* entity;
  ContactInteractionEntity<Real>* interaction;
  
  if( n>Number_Element_Restart_Variables()  || n<=0 ){
    errors->Add_Error_Message("Variable Index is unreasonable");
    return INVALID_ID;
  }
  
  ContactElementElementInteraction* ceei;
  int i;
  if (n==1) {
    i = 0;
    primary_topology->ElemList()->IteratorStart();
    while ((entity=primary_topology->ElemList()->IteratorForward())) {
      ContactElement* element = static_cast<ContactElement*>(entity);
      int num_elementelement_interactions = (int) data[i];
      for( int j=0 ; j<num_elementelement_interactions ; ++j){
        ceei = ContactElementElementInteraction::new_ContactElementElementInteraction(
                allocators[ALLOC_ContactElementElementInteraction]);
        element->Store_ElementElement_Interaction( ceei );
      }
      i++;
    }
  } else {
    i = 0;
    int eei_num = (n-2)/5;
    int offset  = (n-2)-eei_num*5;
    primary_topology->ElemList()->IteratorStart();
    while ((entity=primary_topology->ElemList()->IteratorForward())) {
      ContactElement* element = static_cast<ContactElement*>(entity);
      if (element->Number_ElementElement_Interactions()>eei_num) {
        int cnt=0;
        ContactInteractionDLL<Real>* interactions = element->Get_ElementElement_Interactions();
        interactions->IteratorStart();
        while ((interaction=interactions->IteratorForward())){
          if (cnt == eei_num) {
            ceei = static_cast<ContactElementElementInteraction*> (interaction);
            switch (offset) {
            case 0:
              ceei->SlaveElementEntityData()->host_gid[0] = (int) data[i];
              break;
            case 1:
              ceei->SlaveElementEntityData()->host_gid[1] = (int) data[i];
              break;
            case 2:
              ceei->MasterElementEntityData()->host_gid[0] = (int) data[i];
              break;
            case 3:
              ceei->MasterElementEntityData()->host_gid[1] = (int) data[i];
              break;
            default:
	      VariableHandle VH = (VariableHandle) offset-4;
	      *ceei->Variable( VH ) = data[i];
              break;
            }
            break;
          }
          cnt++;
        }
      }
      i++;
    }
  }
  
  restart = true;
  return NO_ERROR;
}

ContactSearch::ContactErrorCode ContactSearch::Complete_Restart() {
  int total_num_procs = contact_number_of_processors( SearchComm );

  //
  //  Need to recalculate the index_in_host_array, index_in_proc_array, and index_in_owner_proc_array values
  //
  //  Each processor knows this data about faces it owns, need to make this data parallel consistent
  //  Start by building a list of global face ids that this processor needs to know about
  //
  //  Create an vector of faces to request.  Each face is identified with a two integer id number
  //
  vector< vector<int> >needed_global_face_ids(total_num_procs);
  vector< vector<int> >face_ids_to_send(total_num_procs);

  primary_topology->NodeList()->IteratorStart();
  while (ContactTopologyEntity<Real>* entity=primary_topology->NodeList()->IteratorForward()) {  
    ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);


    ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions();
    for (int i=0; i<node->Number_NodeEntity_Interactions(); ++i) {
      if (interactions[i]->Get_Type()!=ContactNodeEntityInteraction::NODE_FACE_INTERACTION) continue;
      ContactNodeFaceInteraction* cnfi = static_cast<ContactNodeFaceInteraction*>(interactions[i]);
      ContactInteractionEntity<Real>::entity_data *face_entity_data = cnfi->FaceEntityData();
      int global_id_0 = face_entity_data->host_gid[0];
      int global_id_1 = face_entity_data->host_gid[1];
      int owning_proc = face_entity_data->owner;
      needed_global_face_ids[owning_proc].push_back(global_id_0);
      needed_global_face_ids[owning_proc].push_back(global_id_1);
    }
  }
  //
  //  Exchange data arrays so that processors that own a face know which processer needs information about 
  //  that face
  //
#ifndef CONTACT_NO_MPI
  ACME::Parallel_Data_Exchange(needed_global_face_ids, face_ids_to_send, SearchComm);
#endif
  //
  //  Fill the face data export array for the locally owned faces.
  //
  vector< vector<int> >face_data_export(total_num_procs);
  vector< vector<int> >face_data_import(total_num_procs);
  for(int iproc = 0; iproc < total_num_procs; ++iproc) {
    for(int iface = 0; iface < face_ids_to_send[iproc].size(); iface += 2) {
      int global_id_0 = face_ids_to_send[iproc][iface + 0];
      int global_id_1 = face_ids_to_send[iproc][iface + 1];

      //
      //  Look up the current face
      //
      ContactHostGlobalID gid(global_id_0, global_id_1);
      ContactTopologyEntity<Real> *face = primary_topology->FaceList()->Find(gid);
      face_data_export[iproc].push_back(face->HostArrayIndex());
      face_data_export[iproc].push_back(face->ProcArrayIndex());
      face_data_export[iproc].push_back(face->OwnerProcArrayIndex());
    }
  }
  //
  //  Copy the export data arrays to the import data arrays
  //
#ifndef CONTACT_NO_MPI
  ACME::Parallel_Data_Exchange(face_data_export, face_data_import, SearchComm);
#endif
  //
  //  Use the imported data to finalize the face interaction indexes
  //
  vector<int> import_list_index(total_num_procs, 0);
  primary_topology->NodeList()->IteratorStart();
  while (ContactTopologyEntity<Real>* entity=primary_topology->NodeList()->IteratorForward()) {  
    ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
    ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions();
    for (int i=0; i<node->Number_NodeEntity_Interactions(); ++i) {
      if (interactions[i]->Get_Type()!=ContactNodeEntityInteraction::NODE_FACE_INTERACTION) continue;
      ContactNodeFaceInteraction* cnfi = static_cast<ContactNodeFaceInteraction*>(interactions[i]);
      ContactInteractionEntity<Real>::entity_data *face_entity_data = cnfi->FaceEntityData();
      int owning_proc = face_entity_data->owner;
      face_entity_data->index_in_host_array       = face_data_import[owning_proc][import_list_index[owning_proc]+0];
      //face_entity_data->index_in_proc_array       = face_data_import[owning_proc][import_list_index[owning_proc]+1];
      face_entity_data->index_in_owner_proc_array = face_data_import[owning_proc][import_list_index[owning_proc]+2];
      import_list_index[owning_proc] += 3;
    }
  }
  //
  // We need to connect the face* for the interactions whose faces are 
  // on processor.  Otherwise we will dereference NULL pointers in the
  // enforcement.
  //
  primary_topology->NodeList()->IteratorStart();
  while (ContactTopologyEntity<Real>* entity=primary_topology->NodeList()->IteratorForward()) {
    ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
    ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions();
    for (int i=0; i<node->Number_NodeEntity_Interactions(); ++i) {
      if (interactions[i]->Get_Type()==ContactNodeEntityInteraction::NODE_FACE_INTERACTION) {
        ContactNodeFaceInteraction* cnfi = static_cast<ContactNodeFaceInteraction*>(interactions[i]);
        cnfi->Connect_Face( *primary_topology->FaceList() );
      }
    }
  }
  primary_topology->FaceList()->IteratorStart();
  while (ContactTopologyEntity<Real>* entity=primary_topology->FaceList()->IteratorForward()) {
    ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
    ContactInteractionDLL<Real>* interactions = face->Get_FaceFace_Interactions();
    if(interactions == NULL) continue;
    interactions->IteratorStart();
    while (ContactInteractionEntity<Real>* interaction=interactions->IteratorForward()){
      ContactFaceFaceInteraction<Real>* cffi = 
        static_cast<ContactFaceFaceInteraction<Real>*>(interaction);
      cffi->Connect_SlaveFace( *primary_topology->FaceList() );
      cffi->Connect_MasterFace( *primary_topology->FaceList() );
    }
  }
  primary_topology->FaceList()->IteratorStart();
  while (ContactTopologyEntity<Real>* entity=primary_topology->FaceList()->IteratorForward()) {
    ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
    ContactInteractionDLL<Real>* interactions = face->Get_FaceCoverage_Interactions();
    if(interactions != NULL) {
      interactions->IteratorStart();
      while (ContactInteractionEntity<Real>* interaction=interactions->IteratorForward()){
        ContactFaceCoverageInteraction* cfci = 
          static_cast<ContactFaceCoverageInteraction*>(interaction);
        cfci->Connect_SlaveFace( *primary_topology->FaceList() );
      }
    }
  }
  primary_topology->ElemList()->IteratorStart();
  while (ContactTopologyEntity<Real>* entity=primary_topology->ElemList()->IteratorForward()) {
    ContactElement* element = static_cast<ContactElement*>(entity);
    ContactInteractionDLL<Real>* interactions = element->Get_ElementElement_Interactions();
    interactions->IteratorStart();
    while (ContactInteractionEntity<Real>* interaction=interactions->IteratorForward()){
      ContactElementElementInteraction* ceei = 
             static_cast<ContactElementElementInteraction*> (interaction);
      ceei->Connect_SlaveElement( *primary_topology->ElemList() );
      ceei->Connect_MasterElement( *primary_topology->ElemList() );
    }
  }
  // ASG Note 6/8/06:
  //
  //    The following note and code has been in ACME for quite a while.
  //    I think I put it in long ago due to a problem from calore or
  //    adagio in their contact implementation when doing restart.
  //    However, this same call breaks shell restart because it
  //    changes the lofting, especially the previous lofting. I am
  //    leaving it here as a place holder in case adagio or calore
  //    report a problem with restart.

//   // Recompute all the surface geometry in case the code does an
//   // enforcement before a search.  We will use the CURRENT_POSITION
//   // since its required for all searches.
//
//   VariableHandle CURRENT_POSITION =
//     primary_topology->Variable_Handle( ContactTopology::Current_Position );
//   primary_topology->Compute_Surface_Geometry( CURRENT_POSITION ,1);

  restart = true;
  return NO_ERROR;
}
