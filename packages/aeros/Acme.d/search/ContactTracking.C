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


#include "ContactAnalyticCylinderInside.h"
#include "ContactAnalyticCylinderOutside.h"
#include "ContactAnalyticPlane.h"
#include "ContactAnalyticSphere.h"
#include "ContactAnalyticSurface.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactNodeFaceInteraction.h"
#include "ContactNodeSurfaceInteraction.h"
#include "ContactParOStream.h"
#include "ContactSearch.h"
#include "ContactSearchData.h"
#include "ContactTopology.h"
#include "search_methods.h"

#include <cmath>

ContactSearch::ContactErrorCode 
ContactSearch::Dynamic_Tracking_Search_2_Configuration()
{
  errors->Add_Error_Message(
	"Dyanmic_Tracking_Search_2_Configuration not yet implemented");
  return UNIMPLEMENTED_FUNCTION;

#if 0
  int i,j,k;

#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  Real start_serial_search_time = MPI_Wtime();
#endif

  error_code = NO_ERROR;

  primary_topology->Update_State();

  // Get the variable handles for everything I need
  // The variables are the same in either decomposition so use primary
  VariableHandle CURRENT_POSITION = 
    primary_topology->Variable_Handle( ContactTopology::Current_Position );
  VariableHandle PREDICTED_POSITION = 
    primary_topology->Variable_Handle( ContactTopology::Predicted_Position );
  VariableHandle NODE_NORMAL = 
    primary_topology->Variable_Handle( ContactTopology::Node_Normal );
  VariableHandle FACE_NORMAL = 
    primary_topology->Variable_Handle( ContactTopology::Face_Normal );

  primary_topology->Set_Search_Options( multiple_interaction_status,
					normal_smoothing_status,
					sharp_smooth_curvature );
  Compute_Surface_Geometry(PREDICTED_POSITION, false);

  Create_Search_Topology(PREDICTED_POSITION);

  // Get the scratch memory
  search_topology->Register_Physical_Face_Scratch_Variables();
  search_topology->Get_Scratch();

  // Build the physical face list for each node
  Build_Physical_Face_List();
  
#ifdef CONTACT_DEBUG_NODE
  if( debug_node_set ){
    if( primary_topology->Debug_Node_Global_ID().Owning_Processor() == 
	contact_processor_number( SearchComm ) ){
      postream << "\nDebug Node Information for Node " 
	       << primary_topology->Debug_Node()->HostCodeGID() << "\n";
    }
  }
#endif
  
  // Retrieve any tied interactions
  Retrieve_Tied_Interactions();

  ContactInteractions* Interactions = search_topology->Interactions();
  
  for( j=0 ; j<Max_Interactions_Per_Node() ; ++j){
    for( k=0 ; k<search_topology->Number_of_Node_Blocks() ; ++k){
      ContactNodeBlock* node_block = search_topology->Node_Block(k);
      int num_nodes_in_block = node_block->Number_of_Nodes();
      ContactNode<Real>** nodes = node_block->Nodes();
      for( i=0 ; i<num_nodes_in_block ; ++i){
        ContactNodeFaceInteraction* cnfi = 
          nodes[i]->Get_NodeFace_Interaction( j, 1 );
        if( cnfi ){
	  Update_Interaction( cnfi->Node(), cnfi->Face(),
			      PREDICTED_POSITION,
			      NODE_NORMAL,
			      FACE_NORMAL );
	}
      }
    }
  }
#if CONTACT_DEBUG_PRINT_LEVEL>=5
  if( contact_processor_number(SearchComm) == 0 )
    std::cout << "\n\n\n***** Interactions before Interaction Definition******* \n" 
	 << std::endl;
  search_topology->Interactions()->Display_NodeFace_Interactions(postream);
  postream.flush();
#endif

  Interaction_Definition();

#if CONTACT_DEBUG_PRINT_LEVEL>=4
  if( contact_processor_number(SearchComm) == 0 )
    std::cout << "\n\n\n***** Interactions after Interaction Definition******* \n" 
	 << std::endl;
  search_topology->Interactions()->Display_NodeFace_Interactions(postream);
  postream.flush();
#endif

  // Process any analytic Surfaces
 
#if CONTACT_DEBUG_PRINT_LEVEL>=4
  if( contact_processor_number(SearchComm) == 0 )
    std::cout << "\n\n\nNode Surface Interactions" << std::endl;
  search_topology->Interactions()->Display_NodeSurface_Interactions(postream);
  postream.flush();
#endif

  search_topology->Release_Scratch();

#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  serial_search_time += MPI_Wtime() - start_serial_search_time;
#endif

  Define_Primary_Interactions();
  
  if( secondary_topology ){
    secondary_topology->CleanUp();
  }

#ifdef CONTACT_DEBUG_NODE
  primary_topology->Display_Debug_Node( postream );
  postream.flush();
#endif

  return (ContactErrorCode) contact_global_error_check( error_code, communicator );
#endif
}



void ContactSearch::Update_Interaction( ContactNode<Real>* node, ContactFace<Real>* face,
					VariableHandle POSITION, 
					VariableHandle NODE_NORMAL,
					VariableHandle FACE_NORMAL )
{
  int i;
  
  // We follow the following steps to find the new interactions
  // 1) Process with the current face
  // 2) If its inside the face with no close edges, we're done.
  // 3) If it has close edges (which will is true if we are outside),
  //    then
  //    a) If it has one close edge, then go across the edge and process
  //       with that edge (no recursion for Arne's special case yet)
  //    b) If it has two close edges, go to the node shared by these edges
  //       and process with all faces connected to this node (considering
  //       the multiple interaction status and sharp-non_sharp angle)


  //
  // ***********************
  // *    S T E P   1      *
  // ***********************
  //
  Process_Tracking_Face( face, node, POSITION, NODE_NORMAL, FACE_NORMAL );
  
  int number_close_edges=0,edge_0,edge_1;
  Real local_coords[3];
  face->Get_Close_Edges( local_coords, number_close_edges, edge_0, edge_1 );
  ContactFace<Real>* neighbor;
  switch( number_close_edges ){
  case 0:
    // We are strictly inside this face and should only process this face
    return;
  case 1:
    // Grab the neighbor face and process
    neighbor = face->Edge( edge_0 )->Neighbor_Face( face );
    Process_Tracking_Face( neighbor, node, POSITION, NODE_NORMAL, FACE_NORMAL );
    break;
  case 2:
    // Grab the node at the shared edge
    ContactEdge<Real>* edge0 = face->Edge( edge_0 );
    ContactNode<Real>* node0_e0 = edge0->Node(0);
    ContactNode<Real>* node1_e0 = edge0->Node(edge0->Nodes_Per_Edge()-1);
    ContactEdge<Real>* edge1 = face->Edge( edge_1 );
    ContactNode<Real>* node0_e1 = edge1->Node(0);
    ContactNode<Real>* common_node;
    if( node0_e1 == node0_e0 || node0_e1 == node1_e0 )
      common_node = node0_e1;
    else
      common_node = edge1->Node(edge1->Nodes_Per_Edge()-1);
    for( i=0 ; i<common_node->Number_Face_Connections() ; ++i){
      neighbor = common_node->GetFace( i );
      if( neighbor != face ){
	Process_Tracking_Face( neighbor, node, POSITION, NODE_NORMAL, 
			       FACE_NORMAL );
      }
    }
    break;
  }
}



void ContactSearch::Process_Tracking_Face( ContactFace<Real>* face, 
					   ContactNode<Real>* node,
					   VariableHandle POSITION,
					   VariableHandle NODE_NORMAL,
					   VariableHandle FACE_NORMAL )
{
#if 0
  int i,ii,nn;
 // KHB: For efficiency, this should be done once and passed in
  Compute_Max_Relative_Node_Motion( max_node_motion, search_topology );
  Real mag_node_motion = 0;
  for( i=0 ; i<dimensionality ; ++i) 
    mag_node_motion += max_node_motion[i]*max_node_motion[i];
  mag_node_motion = std::sqrt( mag_node_motion );
  Real user_search_tol = std::max(search_data->Max_Search_Normal_Tolerance(),
			     search_data->Max_Search_Tangential_Tolerance() );
  Real ctoln2 = std::max(mag_node_motion,user_search_tol);
  // The computed areas in cnodetriange_cpproj for a point exactly on a
  // corner (or possibly on an edge) may yield a very small negative number.
  // CTOLT1 is used in an if-test to avoid bad consequences from this
  // possibility.  A better approach to handle this is desirable. (RMS)
  Real CTOLT1 = 1.e-12;
  Real CTOLT2 = search_data->Max_Search_Tangential_Tolerance();

  // KHB: This memory needs to be handled a better way.
  Real ms_normals[3*24];
  Real ctrl[24];
  Real ctrcl[12];
  Real ctrcl_facets[384];
  Real coordinates[288];
  int pushback_dir[32];
  int nfacets;
  face->FacetDecomposition(nfacets,coordinates,ms_normals,POSITION);
  Real* coords_node = node->Variable(POSITION);
  Real* node_normal = node->Variable(NODE_NORMAL);
  for( nn=0 ; nn<nfacets ; nn++ ){
    pushback_dir[nn] = 2; // normal pushback
  }
  for( ii=0 ; ii<nfacets ; ii++ ){
    ctrcl_facets[ii*12] = 0.0;
  }
  for( nn=0 ; nn<nfacets ; nn++ ){
    if (dimensionality==3) {
      int npairs = 1;
      FORTRAN(cnodetriangle_cpproj)(npairs, coords_node, 
				    &coordinates[nn*9],
				    &ms_normals[nn*3],
				    &pushback_dir[nn],
				    ctrl,
				    &ctrcl_facets[nn*12],
				    CTOLT1,CTOLT2,ctoln2);
    }
  }
  face->FacetStaticRestriction(nfacets, coordinates, ms_normals, 
			       ctrcl_facets, ctrcl);
  if( ctrcl[0] == 1 ){
    ContactNodeFaceInteraction* cnfi = 
      ContactNodeFaceInteraction::new_ContactNodeFaceInteraction( 
		allocators[ALLOC_ContactNodeFaceInteraction],
		ContactNodeFaceInteraction::CLOSEST_POINT_PROJECTION_2, 
		node, face, ctrcl, node->Entity_Key(0),
		face->Variable(FACE_NORMAL) );
    Process_Interaction( cnfi );
  }
#endif
}
