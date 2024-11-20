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
#include "ContactUtilities.h"
#include "ContactAnalyticCylinderInside.h"
#include "ContactAnalyticCylinderOutside.h"
#include "ContactAnalyticPlane.h"
#include "ContactAnalyticSphere.h"
#include "ContactAnalyticSurface.h"
#include "ContactCommBuffer.h"
#include "ContactEdgeBlock.h"
#include "ContactEnforcement.h"
#include "ContactErrors.h"
#include "ContactFaceBlock.h"
#include "ContactNodeFaceInteraction.h"
#include "ContactNodeSurfaceInteraction.h"
#include "ContactFaceFaceInteraction.h"
#include "ContactElementElementInteraction.h"
#include "ContactElement.h"
#include "ContactLineEdgeL2.h"
#include "ContactLineEdgeQ3.h"
#include "ContactLineFaceL2.h"
#include "ContactLineFaceQ3.h"
#include "ContactShellNode.h"
#include "ContactTriFaceL3.h"
#include "ContactQuadFaceL4.h"
#include "ContactQuadFaceQ8.h"
#include "ContactQuadFaceQ9.h"
#include "ContactTriFaceQ6.h"
#include "ContactShellQuadFaceL4.h"
#include "ContactShellTriFaceL3.h"
#include "ContactHexElementL8.h"
#include "ContactWedgeElementL6.h"
#include "ContactSearch.h"
#include "ContactTopology.h"
#include "contact_sorting.h"
#include "CString.h"
#include "ContactSearchData.h"
#include "Contact_Communication.h"
#include "ContactSymComm.h"
#include "search_methods.h"
#include "ContactParOStream.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactSequentialAllocator.h"
#include "ContactTable.h"
#include "ContactBoundingBox.h"
#include "ContactBoundingBoxHierarchy.h"
#include "ContactBoundingBoxHierarchy_Int.h"
#include "contact_tolerances.h"

#ifndef CONTACT_NO_MPI
#include "mpi.h"
#include "ContactZoltan.h"
#include "Zoltan_Interface.h"
#include "ContactZoltanCommUtils.h"
#endif

#include <algorithm>
#include <iostream>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

using namespace std;

using acme::Cross;
using acme::Dot;
using acme::Magnitude;

#define NO_CONTACT_HSFC

char* ACME_Version() { return (char*) "v2.7f"; };
char* ACME_VersionDate() { return (char*) "08/10/07"; };

int ACME_MPI_Compatibility(int host_compile)
{
  if (host_compile != MPI_COMPILE) {
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    std::cerr << "Host and library MPI compilations are incompatible" << std::endl;
#endif
    return 1;
  }
  return 0;
}





ContactSearch::ContactSearch( int  Dimensionality,
                              int  Number_of_States,
                              int  Number_of_Analytic_Surfaces,
                              int  Number_of_Node_Blocks,
                              const ContactNode_Type* Node_Block_Types,
                              const int* Number_Nodes_in_Blocks,
			      const int* Node_Exodus_IDs,
                              const int* Node_Global_IDs,
			      const double* coords,
                              int  Number_of_Face_Blocks,
                              const ContactFace_Type* Face_Block_Types,
                              const int* Number_Faces_in_Blocks,
                              const int* Face_Global_IDs,
                              const int* Face_Connectivity,
			      const Real* Face_Lofting_Factors,
                              int  Number_of_Element_Blocks,
                              const ContactElement_Type* Element_Block_Types,
                              const int* Number_Elements_in_Blocks,
                              const int* Element_Global_IDs,
                              const int* Element_Connectivity,
                              int  Number_of_Nodal_Comm_Partners,
                              const int* Node_Comm_Proc_IDs,
                              const int* Number_Nodes_To_Partner,
                              const int* Communication_Nodes,
                              MPI_Comm& mpi_communicator,
                              ContactErrorCode& error )
  : initialized(false), 
    initialized_tied(false), 
    initialized_context(false), 
    initializing_tied(false),
    step_number(0),
    postream( mpi_communicator ),
    timer( mpi_communicator ),
#ifndef CONTACT_NO_MPI
    zoltan(NULL),
#endif
    num_tracked_nodes(0),
    global_tracking_interval(0),
    tracking_step(0),
    restart(false)
{
  capture_motion    = 0.0;
  delete_ghosting   = 0;
  ws_size           = 1;
  max_facets        = 32;
  data_size         = ContactFace<Real>::LENGTH;
  ctrl              = new Real[24*ws_size];
  ctrcl             = new Real[data_size*ws_size];
  ctrcl_facets      = new Real[data_size*max_facets*ws_size];
  pushback_dir_flag = new int[max_facets];
  ms_coordinates_c  = new Real[9*max_facets*ws_size];
  ms_coordinates_a  = new Real[9*max_facets*ws_size];
  ms_coordinates_p  = new Real[9*max_facets*ws_size];
  ms_normals_c      = new Real[3*max_facets*ws_size];
  ms_normals_a      = new Real[3*max_facets*ws_size];
  ms_normals_p      = new Real[3*max_facets*ws_size];


  // compute initial internal tolerances
  box_inflation = BOX_INFLATION_FACTOR;
  gap_inflation = GAP_INFLATION_FACTOR+1.0;

  //
  //  Initialize static face data arrays for fast lookup of face type
  //  information from the base class
  //
  ContactFace<Real>::Initialize_Lookup_Arrays();
  ContactEdge<Real>::Initialize_Lookup_Arrays();

  int i;
  error              = NO_ERROR;
  errors             = new ContactErrors();
  SearchComm         = mpi_communicator;
  comm_buffer        = new ContactCommBuffer();
  primary_topology   = NULL;
  secondary_topology = NULL;
  search_data        = NULL;
  tables             = NULL;
  num_tables         = 0;
  
  number_registered_enforcements = 0;

  enforcement = NULL;

  if( !(Dimensionality==2 || Dimensionality==3) ){
    errors->Add_Error_Message("Invalid Dimesion");
    error = INVALID_DATA;
    return;
  }
  if( Number_of_Node_Blocks < 1 ){
    errors->Add_Error_Message("Number of Node Blocks must be positive");
    error = INVALID_DATA;
    return;
  }
  if( Number_of_Face_Blocks<0 ){
    errors->Add_Error_Message("# of Face Blocks must be positive");
    error = INVALID_DATA;
    return;
  }
  if( Number_of_Element_Blocks < 0 ){
    errors->Add_Error_Message("# of Element Blocks must be positive");
    error = INVALID_DATA;
    return;
  }
  if( Number_of_Analytic_Surfaces < 0 ) {
    errors->Add_Error_Message("# of Analytic Surfaces must be positive");
    error = INVALID_DATA;
    return;
  }
  if( Number_of_Face_Blocks > 0 && Number_of_Element_Blocks > 0 ) {
    errors->Add_Error_Message("Cannot have both Element and Face blocks");
    error = INVALID_DATA;
    return;
  }
  
  #ifndef CONTACT_NO_MPI
  for (i=0; i<Number_of_Node_Blocks; ++i) {
    int type     = Node_Block_Types[i];
    int max_type = contact_global_maximum(type, SearchComm);
    int min_type = contact_global_minimum(type, SearchComm);
    if (max_type!=min_type) {
      std::sprintf(message,"Node type for node block %d is not consistent across all processors",i);
      errors->Add_Error_Message(message);
      error = INVALID_DATA;
    }
  }
  for (i=0; i<Number_of_Face_Blocks; ++i) {
    int type     = Face_Block_Types[i];
    int max_type = contact_global_maximum(type, SearchComm);
    int min_type = contact_global_minimum(type, SearchComm);
    if (max_type!=min_type) {
      std::sprintf(message,"Face type for face block %d is not consistent across all processors",i);
      errors->Add_Error_Message(message);
      error = INVALID_DATA;
    }
  }
  for (i=0; i<Number_of_Element_Blocks; ++i) {
    int type     = Element_Block_Types[i];
    int max_type = contact_global_maximum(type, SearchComm);
    int min_type = contact_global_minimum(type, SearchComm);
    if (max_type!=min_type) {
      std::sprintf(message,"Elem type for element block %d is not consistent across all processors",i);
      errors->Add_Error_Message(message);
      error = INVALID_DATA;
    }
  }
  if (error != NO_ERROR) return;
  #endif

#ifdef CONTACT_NO_MPI
  // in serial debug mode, do preconditions to affirm that if there 
  // are a non-zero number of node, face, or element blocks, then the 
  // corresponding arrays of information about the node/face blocks 
  // are valid.
  PRECONDITION( Number_of_Node_Blocks == 0 || Node_Block_Types );
  PRECONDITION( Number_of_Node_Blocks == 0 || Node_Exodus_IDs );
  PRECONDITION( Number_of_Node_Blocks == 0 || Node_Global_IDs );
  PRECONDITION( Number_of_Node_Blocks == 0 || Number_Nodes_in_Blocks );
  PRECONDITION( Number_of_Face_Blocks == 0 || Face_Block_Types );
  PRECONDITION( Number_of_Face_Blocks == 0 || Face_Global_IDs );
  PRECONDITION( Number_of_Face_Blocks == 0 || Number_Faces_in_Blocks );
  PRECONDITION( Number_of_Face_Blocks == 0 || Face_Connectivity );
  PRECONDITION( Number_of_Element_Blocks == 0 || Element_Connectivity );
#endif

  Set_Up_Allocators();
  
  dimensionality = Dimensionality;
  num_states = Number_of_States;

  primary_topology = new ContactTopology( errors, 
                                          dimensionality,
					  Number_of_Analytic_Surfaces, 
					  Number_of_Node_Blocks,
					  Node_Block_Types,
					  Number_Nodes_in_Blocks,
					  Node_Exodus_IDs,
					  Node_Global_IDs,
					  coords,
					  Number_of_Face_Blocks, 
					  Face_Block_Types,
					  Number_Faces_in_Blocks,
					  Face_Global_IDs, 
					  Face_Connectivity,
					  Face_Lofting_Factors,
					  Number_of_Element_Blocks, 
					  Element_Block_Types,
					  Number_Elements_in_Blocks,
					  Element_Global_IDs, 
					  Element_Connectivity,
					  Number_of_Nodal_Comm_Partners, 
					  Node_Comm_Proc_IDs, 
					  Number_Nodes_To_Partner,
					  Communication_Nodes,
					  comm_buffer,
					  SearchComm,
                                          this,
					  error );

  primary_topology->topology_type = ContactTopology::PRIMARY;

  // If the topology is bad, simply return
  error = (ContactErrorCode) contact_global_error_check( error, SearchComm );
  if( error ) return;

  // Store the coordinates into the topology
  int c_offset = 0;
  for( i=0 ; i<Number_of_Node_Blocks ; ++i ){
    primary_topology->Set_NodeBlk_Positions( i+1, coords+c_offset );
    c_offset += dimensionality*Number_Nodes_in_Blocks[i];
  }

  search_data = new ContactSearchData(primary_topology);

  // Set all the default search options
  multiple_interaction_status = INACTIVE;
  normal_smoothing_status     = INACTIVE;
  compute_node_areas          = INACTIVE;
  partition_gap_status        = INACTIVE;
  new_tied_enforcement        = INACTIVE;
  old_dynamic_search          = INACTIVE;
  enable_tracking             = INACTIVE;
  enable_off_face_tracking    = INACTIVE;
  no_secondary                = INACTIVE;
  no_ghosting                 = INACTIVE;
  keep_ghosting               = INACTIVE;
  search_cull                 = INACTIVE;
  no_warped_volume            = INACTIVE;
  no_parallel_consistency     = INACTIVE;
  smoothing_resolution        = USE_NODE_NORMAL;
  orig_sharp_smooth_angle     = 60.0;
  sharp_smooth_curvature      = ComputeCurvatureFromAngle(orig_sharp_smooth_angle);
  normal_smoothing_distance   = 0.0;
  global_search_cull          = SLAVE_CULL;
  tracking_type               = NO_TRACKING;
  auto_tol                    = INACTIVE;
  aggressive_tolerances       = INACTIVE;
  skip_physical_faces         = INACTIVE;
  edge_physical_faces         = INACTIVE;
  physical_face_algorithm     = PF_FACE_BASED;
  shell_simple_lofting        = INACTIVE;
  compute_partials            = INACTIVE;
  computed_partials_order     = 0;

  //===================================================================
  // Create a topology object with the current entity block structures
  //===================================================================
  int num_edge_blocks = primary_topology->Number_of_Edge_Blocks();
  ContactEdge_Type* edge_block_types = new ContactEdge_Type[num_edge_blocks];
  for( i=0 ; i<num_edge_blocks ; ++i )
    edge_block_types[i] = primary_topology->Edge_Block(i)->Type();

  secondary_topology = new ContactTopology( errors,
                                            dimensionality,
                                            Number_of_Analytic_Surfaces, 
                                            Number_of_Node_Blocks,
                                            Node_Block_Types,
                                            num_edge_blocks,
                                            edge_block_types,
                                            Number_of_Face_Blocks, 
                                            Face_Block_Types,
                                            Number_of_Element_Blocks, 
                                            Element_Block_Types,
                                            SearchComm,
                                            this );

  secondary_topology->topology_type = ContactTopology::SECONDARY;

  delete [] edge_block_types;

#ifndef CONTACT_NO_MPI
  if (contact_number_of_processors(SearchComm)>1) {
    create_zoltan_object(error);
  }
#endif

  // Initizalize all the "registered" data to zero
  // This only needs to be done once here (not each time a new object is
  // created in the secondary decomposition) so its not done as part
  // of the constructor.
  int number_of_nodes = primary_topology->Number_of_Nodes();
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(primary_topology->NodeList()->EntityList());
  for (i=0; i<number_of_nodes; ++i) {
    Nodes[i]->Initialize_Memory();
  }
  int number_of_edges = primary_topology->Number_of_Edges();
  ContactEdge<Real>** Edges = 
    reinterpret_cast<ContactEdge<Real>**>(primary_topology->EdgeList()->EntityList());
  for (i=0; i<number_of_edges; ++i) {
    Edges[i]->Initialize_Memory();
  }
  int number_of_faces = primary_topology->Number_of_Faces();
  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(primary_topology->FaceList()->EntityList());
  for (i=0; i<number_of_faces; ++i) {
    Faces[i]->Initialize_Memory();
  }
  int number_of_elements = primary_topology->Number_of_Elements();
  ContactElement** Elements = 
    reinterpret_cast<ContactElement**>(primary_topology->ElemList()->EntityList());
  for (i=0; i<number_of_elements; ++i) {
    Elements[i]->Initialize_Memory();
  }

#ifdef CONTACT_TIMINGS
  Register_Timers();
#endif

#if CONTACT_DEBUG_PRINT_LEVEL>=2
  postream << "Search Object is Constructed\n";
  postream << "  Number of Nodes    = " 
	   << primary_topology->Number_of_Nodes() << "\n";
  postream << "  Number of Edges    = " 
	   << primary_topology->Number_of_Edges() << "\n";
  postream << "  Number of Faces    = " 
	   << primary_topology->Number_of_Faces() << "\n";
  postream << "  Number of Elements = " 
	   << primary_topology->Number_of_Elements() << "\n";
  postream << "  Number of Surfaces = " 
	   << primary_topology->Number_of_Analytic_Surfaces() << "\n";
  postream.flush();
#endif


  reasonable_gap        = 0.0;
  max_node_motion[0]	= 0.0;
  max_node_motion[1]	= 0.0;
  max_node_motion[2]	= 0.0;
  max_node_displacement = 0.0;
  max_remaining_gaps[0] = 0.0;
  max_remaining_gaps[1] = 0.0;
  max_remaining_gaps[2] = 0.0;
  max_remaining_gap_mag = 0.0;
}

ContactSearch::~ContactSearch()
{
#ifdef CONTACT_TIMINGS
  if( contact_processor_number( SearchComm ) == 0 ) {
    char* offon[2]  = { "off", "on" };
    char* sm_res[2] = { "NODE_BASED", "EDGE_BASED" };
    char* pf_alg[3] = { "NONE", "FACE_BASED", "EDGE_BASED" };
    char* cull[2]   = { "NONE", "SLAVE" };
    char* track[3]  = { "NONE", "LOCAL", "GLOBAL" };
    std::cout << "ACME Options/Data ----------------------------------------\n";
    
    std::cout << "  No Secondary              = " << offon[no_secondary] << "\n";
    std::cout << "  Search Cull               = " << cull[global_search_cull] << "\n";
    std::cout << "  Keep Ghosting             = " << offon[keep_ghosting] << "\n";
    std::cout << "  No Parallel Consistency   = " << offon[no_parallel_consistency] << "\n";
    
    if (do_node_face_search) {
    std::cout << "  Tracking Type             = " << track[tracking_type] << "\n";
    if (tracking_type==GLOBAL_TRACKING) {
    std::cout << "    Tracking Interval       = " << global_tracking_interval << "\n";
    }
    if (tracking_type==LOCAL_TRACKING) {
    std::cout << "    Tracking Interval       = " << global_tracking_interval << "\n";
    std::cout << "    Off Face Tracking       = " << offon[enable_off_face_tracking] << "\n";
    }
    std::cout << "  Old Dynamic Search        = " << offon[old_dynamic_search] << "\n";
    std::cout << "  Physical Face Algorithm   = " << pf_alg[physical_face_algorithm]     << "\n";
    std::cout << "  Multiple interactions     = " << offon[multiple_interaction_status] << "\n";
    if (multiple_interaction_status){
    std::cout << "    Sharp/Smooth Angle  = " << orig_sharp_smooth_angle << "\n";
    }
    std::cout << "  Normal smoothing          = " << offon[normal_smoothing_status] << "\n";
    if (normal_smoothing_status){
    std::cout << "  Sharp/Smooth Angle    = " << orig_sharp_smooth_angle << "\n";
    std::cout << "  Normal Smoothing Distance = " << normal_smoothing_distance << "\n";
    std::cout << "  Smoothing Resolution      = " << sm_res[smoothing_resolution] << "\n";
    }
    std::cout << "  Partition Gap Status      = " << offon[partition_gap_status] << "\n";
    std::cout << "  Compute Node Areas        = " << offon[compute_node_areas] << "\n";
    std::cout << "  Auto Tolerance            = " << offon[auto_tol] << "\n";
    if (auto_tol) {
    std::cout << "    aggressive_tolerances   = " << offon[aggressive_tolerances]   << "\n";
    }
    }
    if (do_elem_elem_search) {
      std::cout << "  No warped volume          = " << offon[no_warped_volume] << "\n";
    }
    if (primary_topology->Have_Shells())
       std::cout << "  Simple shell lofting      = " << offon[compute_node_areas] << "\n";
  }
  timer.Output_Timings();
  contact_global_sync(SearchComm);
#endif
  postream.flush();  // just to make sure
  if( primary_topology   ) delete primary_topology;
  if( secondary_topology ) delete secondary_topology;
  if( search_data ) delete search_data;
  if( errors ) delete errors;
  if( comm_buffer ) delete comm_buffer;
#ifndef CONTACT_NO_MPI
  if (contact_number_of_processors(SearchComm)>1) {
    if( zoltan ) delete zoltan;
  }
#endif
  
  if (allocators) {
    for (int i=0; i<ALLOC_NUM_ALLOCATED_ENTITIES; ++i) {
      allocators[i].Purge();
    }
    delete [] allocators;
  }
  delete scratch_allocator;
#if (MAX_FFI_DERIVATIVES > 0)
  if (active_allocators) {
    for (int i=0; i<ALLOC_NUM_ALLOCATED_ENTITIES; ++i) {
      active_allocators[i].Purge();
    }
    delete [] active_allocators;
  }
#endif
  if( num_tables ){
    for( int i=0 ; i<num_tables ; ++i ) delete tables[i];
    delete [] tables;
  }
  //
  //  Remove any references to the search object from the enforcement object
  //
  if( enforcement ) {
    for( int i=0 ; i<number_registered_enforcements ; ++i ) {
       enforcement[i]->Remove_Search_Reference(this);
    }
    delete [] enforcement;
  }

  primary_topology   = NULL;
  secondary_topology = NULL;
  search_data        = NULL;
  errors             = NULL;
  comm_buffer        = NULL;
#ifndef CONTACT_NO_MPI
  zoltan             = NULL;
#endif
  allocators         = NULL;
#if (MAX_FFI_DERIVATIVES > 0)
  active_allocators  = NULL;
#endif
  tables             = NULL;
  enforcement        = NULL;
  
  if(ctrl) delete [] ctrl;
  if(ctrcl) delete [] ctrcl;
  if(ctrcl_facets) delete [] ctrcl_facets;
  if(pushback_dir_flag) delete [] pushback_dir_flag;
  if(ms_coordinates_c) delete [] ms_coordinates_c;
  if(ms_coordinates_a) delete [] ms_coordinates_a;
  if(ms_coordinates_p) delete [] ms_coordinates_p;
  if(ms_normals_c) delete [] ms_normals_c;
  if(ms_normals_a) delete [] ms_normals_a;
  if(ms_normals_p) delete [] ms_normals_p;
}

void ContactSearch::Reset_Search() {
  step_number   = 0;
  tracking_step = 0;
}

void
ContactSearch::UpdateSearch( const int* num_node_deaths_per_block, 
		             const int* node_deaths_global_ids,
                             const int* num_face_deaths_per_block, 
	                     const int* face_deaths_global_ids,
                             const int* num_element_deaths_per_block, 
	                     const int* element_deaths_global_ids,
                             
                             const int* num_node_births_per_block, 
	                     const int* node_births_exodus_ids,
	                     const int* node_births_global_ids,
	                     const int* number_face_births_per_block, 
	                     const int* face_births_global_ids,
	                     const int* face_births_connectivity,
	                     const int* number_element_births_per_block, 
	                     const int* element_births_global_ids,
	                     const int* element_births_connectivity,
                             
	                     const int  num_node_exports,
                             const int* node_export_id_list,
                             const int* node_export_pid,
	                     const int  num_face_exports,
                             const int* face_export_id_list,
                             const int* face_export_pid,
	                     const int  num_element_exports,
                             const int* element_export_id_list,
                             const int* element_export_pid,
                             
                             const int  num_nodes,
                             const int* node_host_ids,
                             const int  num_faces,
                             const int* face_host_ids,
                             const int  num_elements,
                             const int* element_host_ids,
                             
	                     const int  num_comm_partners,
	                     const int* comm_proc_id,
	                     const int* number_nodes_to_partner,
	                     const int* comm_node,
                             
	                     ContactErrorCode& error )
{
  error = NO_ERROR;

  tables = NULL;
  num_tables = 0;

  if( primary_topology->Have_Shells() ){
    errors->Add_Error_Message(
      "Element Birth/Death, Adaptivity & DLB not yet supported with shells" );
    error = UNIMPLEMENTED_FUNCTION;
    return;
  }
  
  primary_topology->UpdateTopology( errors, 
  
                                    num_node_deaths_per_block, 
                                    node_deaths_global_ids,
                                    num_face_deaths_per_block, 
		                    face_deaths_global_ids,
                                    num_element_deaths_per_block, 
		                    element_deaths_global_ids,
                                    
                                    num_node_births_per_block, 
		                    node_births_global_ids, 
		                    node_births_exodus_ids,
		                    number_face_births_per_block, 
		                    face_births_global_ids,
		                    face_births_connectivity,
		                    number_element_births_per_block, 
		                    element_births_global_ids,
		                    element_births_connectivity,
                                    
		                    num_node_exports,
                                    node_export_id_list,
                                    node_export_pid,
		                    num_face_exports,
                                    face_export_id_list,
                                    face_export_pid,
		                    num_element_exports,
                                    element_export_id_list,
                                    element_export_pid,
                             
                                    num_nodes,
                                    node_host_ids,
                                    num_faces,
                                    face_host_ids,
                                    num_elements,
                                    element_host_ids,
                                    
		                    num_comm_partners,
		                    comm_proc_id,
		                    number_nodes_to_partner,
		                    comm_node,
	
                                    postream,
				    error );

  // If the topology is bad, simply return
  error = (ContactErrorCode) contact_global_error_check( error, SearchComm );
  if( error ) {
    return;
  }
  
  // Initizalize all the "registered" data to zero
  // This only needs to be done once here (not each time a new object is
  // created in the secondary decomposition) so its not done as part
  // of the constructor.
  int i;
  int number_of_nodes = primary_topology->Number_of_Nodes();
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(primary_topology->NodeList()->EntityList());
  for (i=0; i<number_of_nodes; ++i) {
    Nodes[i]->Initialize_Memory();
  }
  int number_of_edges = primary_topology->Number_of_Edges();
  ContactEdge<Real>** Edges = 
    reinterpret_cast<ContactEdge<Real>**>(primary_topology->EdgeList()->EntityList());
  for (i=0; i<number_of_edges; ++i) {
    Edges[i]->Initialize_Memory();
  }
  int number_of_faces = primary_topology->Number_of_Faces();
  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(primary_topology->FaceList()->EntityList());
  for (i=0; i<number_of_faces; ++i) {
    Faces[i]->Initialize_Memory();
  }
  int number_of_elements = primary_topology->Number_of_Elements();
  ContactElement** Elements = 
    reinterpret_cast<ContactElement**>(primary_topology->ElemList()->EntityList());
  for (i=0; i<number_of_elements; ++i) {
    Elements[i]->Initialize_Memory();
  }

#ifdef CONTACT_TIMINGS
  Register_Timers();
#endif

  int  num_new_nodes = primary_topology->Number_of_Nodes();
  int* old_to_new_node_map = primary_topology->UpdateNodeMap();
  for( i=0 ; i<number_registered_enforcements ; ++i ) {
     enforcement[i]->Update_For_Topology_Change( num_new_nodes,
                                                 old_to_new_node_map );
  }
  
#if CONTACT_DEBUG_PRINT_LEVEL>=2
  postream << "Search Object is Updated\n";
  postream << "  Number of Nodes    = " 
	   << primary_topology->Number_of_Nodes() << "\n";
  postream << "  Number of Edges    = " 
	   << primary_topology->Number_of_Edges() << "\n";
  postream << "  Number of Faces    = " 
	   << primary_topology->Number_of_Faces() << "\n";
  postream << "  Number of Elements = " 
	   << primary_topology->Number_of_Elements() << "\n";
  postream.flush();
#endif

}

#ifndef CONTACT_NO_MPI
void ContactSearch::create_zoltan_object( ContactErrorCode& error )
{
  //========================
  // Create a zoltan object 
  //========================
  int zoltan_error;
  zoltan = new ContactZoltan(SearchComm, zoltan_error);
  if( zoltan_error == ZOLTAN_FATAL ){
    error = ZOLTAN_ERROR;
    errors->Add_Error_Message("Error Creating Zoltan Object");
    return;
  }
#ifdef CONTACT_HSFC
  zoltan->Set_Method( (char*)"HSFC" );
#else
  zoltan->Set_Method( (char*)"RCB" );
  zoltan->Set_Param( (char*)"RCB_REUSE",       (char*)"1" );
#endif
  // Check_Geom causes Zoltan to do unnecessary checking, so turn it off
  // For us this is
  //   1) check dot weights > 0 (we don't specify them)
  //   2) Make sure #dots_in = #dots_out
  //   3) Check that the dots are loadbalanced within a tolerance
  //   4) Check all points are within the RCB boxes
  zoltan->Set_Param( (char*)"CHECK_GEOM",      (char*)"0" );
  zoltan->Set_Param( (char*)"KEEP_CUTS",       (char*)"1" );
  zoltan->Set_Param( (char*)"DEBUG_LEVEL",     (char*)"0" );
#if CONTACT_DEBUG_PRINT_LEVEL>=10
  zoltan->Set_Param( (char*)"RETURN_LISTS",    (char*)"ALL");
#else
  zoltan->Set_Param( (char*)"RETURN_LISTS",    (char*)"NONE");
#endif
  char num_lid_entries[2];
  char num_gid_entries[2];
  std::sprintf(num_lid_entries,"%d",ZOLTAN_LID_SIZE);
  std::sprintf(num_gid_entries,"%d",ZOLTAN_GID_SIZE);
  zoltan->Set_Param( (char*)"NUM_LID_ENTRIES", num_lid_entries );
  zoltan->Set_Param( (char*)"NUM_GID_ENTRIES", num_gid_entries );
#if defined(__PUMAGON__) || defined(CONTACT_ZOLTAN_TFLOPS_SPECIAL)
  zoltan->Set_Param( (char*)"TFLOPS_SPECIAL",  (char*)"1" );
#endif
}
#endif

void
ContactSearch::Set_Zoltan_Param(char *zoltan_parameter, char *zoltan_value)
{
#ifndef CONTACT_NO_MPI
  if ( zoltan ) {
    zoltan->Set_Param(zoltan_parameter, zoltan_value);
  }
#endif
}

void
ContactSearch::Set_Zoltan_Method(char *zoltan_method)
{
#ifndef CONTACT_NO_MPI
  if ( zoltan ) {
    zoltan->Set_Method(zoltan_method);
  }
#endif
}

void
ContactSearch::Set_Up_Allocators()
{
  allocators = new ContactFixedSizeAllocator[ALLOC_NUM_ALLOCATED_ENTITIES];
  ContactNode_SizeAllocator<Real>(allocators[ALLOC_ContactNode]);
  ContactLineEdgeL2_SizeAllocator<Real>(allocators[ALLOC_ContactLineEdgeL2]);
  ContactLineEdgeQ3_SizeAllocator(allocators[ALLOC_ContactLineEdgeQ3]);
  ContactLineFaceL2_SizeAllocator(allocators[ALLOC_ContactLineFaceL2]);
  ContactLineFaceQ3_SizeAllocator(allocators[ALLOC_ContactLineFaceQ3]);
  ContactQuadFaceL4_SizeAllocator<Real>(allocators[ALLOC_ContactQuadFaceL4]);
  ContactQuadFaceQ8_SizeAllocator(allocators[ALLOC_ContactQuadFaceQ8]);
  ContactQuadFaceQ9_SizeAllocator(allocators[ALLOC_ContactQuadFaceQ9]);
  ContactTriFaceL3_SizeAllocator<Real>(allocators[ALLOC_ContactTriFaceL3]);
  ContactTriFaceQ6_SizeAllocator(allocators[ALLOC_ContactTriFaceQ6]);
  ContactHexElemL8_SizeAllocator<Real>(allocators[ALLOC_ContactHexElemL8]);
  ContactWedgeElemL6_SizeAllocator<Real>(allocators[ALLOC_ContactWedgeElemL6]);
  ContactNodeNodeInteraction_SizeAllocator(
	 allocators[ALLOC_ContactNodeNodeInteraction]);
  ContactNodeFaceInteraction_SizeAllocator(
	 allocators[ALLOC_ContactNodeFaceInteraction]);
  ContactNodeSurfaceInteraction_SizeAllocator(
	 allocators[ALLOC_ContactNodeSurfaceInteraction]);
  ContactFaceFaceInteraction_SizeAllocator<Real>(
	 allocators[ALLOC_ContactFaceFaceInteraction]);
  ContactFaceCoverageInteraction_SizeAllocator(
	 allocators[ALLOC_ContactFaceCoverageInteraction]);
  ContactElementElementInteraction_SizeAllocator(
	 allocators[ALLOC_ContactElementElementInteraction]);
  ContactPolyVert_SizeAllocator(allocators[ALLOC_ContactPolyVert]);
  ContactPolyEdge_SizeAllocator(allocators[ALLOC_ContactPolyEdge]);
  ContactPoly_SizeAllocator(allocators[ALLOC_ContactPoly]);
  ContactFaceFaceGraphNode_SizeAllocator(allocators[ALLOC_ContactFaceFaceGraphNode]);
  ContactCartesianHexElementL8_SizeAllocator(allocators[ALLOC_ContactCartesianHexElementL8]);
  ContactHexElementL8_SizeAllocator(allocators[ALLOC_ContactHexElementL8]);
  ContactShellQuadFaceL4_SizeAllocator<Real>(allocators[ALLOC_ContactShellQuadFaceL4]);
  ContactShellTriFaceL3_SizeAllocator<Real>(allocators[ALLOC_ContactShellTriFaceL3]);
  ContactShellNode_SizeAllocator(allocators[ALLOC_ContactShellNode]);

  // Right now set the block size to 10000
  scratch_allocator = new ContactSequentialAllocator( 100000,1 ) ;

#if (MAX_FFI_DERIVATIVES > 0)
  active_allocators = new ContactFixedSizeAllocator[ALLOC_NUM_ALLOCATED_ENTITIES];
  ContactNode_SizeAllocator<ActiveScalar>(active_allocators[ALLOC_ContactNode]);
  ContactLineEdgeL2_SizeAllocator<ActiveScalar>(active_allocators[ALLOC_ContactLineEdgeL2]);
  ContactLineEdgeQ3_SizeAllocator(active_allocators[ALLOC_ContactLineEdgeQ3]);
  ContactLineFaceL2_SizeAllocator(active_allocators[ALLOC_ContactLineFaceL2]);
  ContactLineFaceQ3_SizeAllocator(active_allocators[ALLOC_ContactLineFaceQ3]);
  ContactQuadFaceL4_SizeAllocator<ActiveScalar>(active_allocators[ALLOC_ContactQuadFaceL4]);
  ContactQuadFaceQ8_SizeAllocator(active_allocators[ALLOC_ContactQuadFaceQ8]);
  ContactQuadFaceQ9_SizeAllocator(active_allocators[ALLOC_ContactQuadFaceQ9]);
  ContactTriFaceL3_SizeAllocator<ActiveScalar>(active_allocators[ALLOC_ContactTriFaceL3]);
  ContactTriFaceQ6_SizeAllocator(active_allocators[ALLOC_ContactTriFaceQ6]);
  ContactHexElemL8_SizeAllocator<ActiveScalar>(active_allocators[ALLOC_ContactHexElemL8]);
  ContactWedgeElemL6_SizeAllocator<ActiveScalar>(active_allocators[ALLOC_ContactWedgeElemL6]);
  ContactNodeNodeInteraction_SizeAllocator(
         active_allocators[ALLOC_ContactNodeNodeInteraction]);
  ContactNodeFaceInteraction_SizeAllocator(
         active_allocators[ALLOC_ContactNodeFaceInteraction]);
  ContactNodeSurfaceInteraction_SizeAllocator(
         active_allocators[ALLOC_ContactNodeSurfaceInteraction]);
  ContactFaceFaceInteraction_SizeAllocator<ActiveScalar>(
         active_allocators[ALLOC_ContactFaceFaceInteraction]);
  ContactFaceCoverageInteraction_SizeAllocator(
         active_allocators[ALLOC_ContactFaceCoverageInteraction]);
  ContactElementElementInteraction_SizeAllocator(
         active_allocators[ALLOC_ContactElementElementInteraction]);
  ContactPolyVert_SizeAllocator(active_allocators[ALLOC_ContactPolyVert]);
  ContactPolyEdge_SizeAllocator(active_allocators[ALLOC_ContactPolyEdge]);
  ContactPoly_SizeAllocator(active_allocators[ALLOC_ContactPoly]);
  ContactFaceFaceGraphNode_SizeAllocator(active_allocators[ALLOC_ContactFaceFaceGraphNode]);
  ContactCartesianHexElementL8_SizeAllocator(active_allocators[ALLOC_ContactCartesianHexElementL8]);
  ContactHexElementL8_SizeAllocator(active_allocators[ALLOC_ContactHexElementL8]);
  ContactShellQuadFaceL4_SizeAllocator<ActiveScalar>(active_allocators[ALLOC_ContactShellQuadFaceL4]);
  ContactShellTriFaceL3_SizeAllocator<ActiveScalar>(active_allocators[ALLOC_ContactShellTriFaceL3]);
  ContactShellNode_SizeAllocator(active_allocators[ALLOC_ContactShellNode]);
#endif
}

int ContactSearch::Number_of_Errors()
{
  return errors->Number_of_Errors();
}

const char* ContactSearch::Error_Message( int i )
{
  // The interface uses a FORTRAN numbering
  return errors->Error_Message( i-1 );
}

ContactSearch::ContactErrorCode
ContactSearch::Check_Search_Data_Size( int size_data_per_pair,
				       int number_of_entity_keys )
{
  bool Error = false;
  if( size_data_per_pair != NSIZSD ){
    Error = true;
    std::sprintf(message,"Size of Data per Pair in Search_Data does not match");
    errors->Add_Error_Message(message);
    std::sprintf(message,"    Size of Data per Pair should be %d",NSIZSD);
    errors->Add_Error_Message(message);
    std::sprintf(message,"    Host code thinks the size is %d",size_data_per_pair);
    errors->Add_Error_Message(message);
  }
  if( number_of_entity_keys != search_data->Num_Search_Entities() ){
    Error = true;
    std::sprintf(message,"Number of Entity Keys does not match");
    errors->Add_Error_Message(message);
    std::sprintf(message,"   Search Library thinks the value is %d",search_data->Num_Search_Entities());
    errors->Add_Error_Message(message);
    std::sprintf(message,"   Host Code thinks the value is %d",
	    number_of_entity_keys);
    errors->Add_Error_Message(message);
  }
  if( Error ) return INVALID_DATA;
  return NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactSearch::Add_Debug_Node( int Exodus_ID )
{
#ifdef CONTACT_DEBUG_NODE
  ContactErrorCode ec = primary_topology->Add_Debug_Node( Exodus_ID );
  return ec;
#else
  std::sprintf(message,"ACME was compiled without -DCONTACT_DEBUG_NODE.");
  errors->Add_Error_Message(message);
  std::sprintf(message,"  ContactSearch::Debug_Node( int Exodus_ID ) is not available.");
  errors->Add_Error_Message(message);
  return UNIMPLEMENTED_FUNCTION;
#endif
}


ContactSearch::ContactErrorCode
ContactSearch::Set_Search_Option( Search_Option option, 
				  Search_Option_Status status,
				  Real* data )
{
  switch( option ){
  case MULTIPLE_INTERACTIONS:
    multiple_interaction_status = status;
    if( status == ACTIVE ){
      PRECONDITION( data );

      orig_sharp_smooth_angle = *data;
      if( orig_sharp_smooth_angle < 0.0 ){
	std::sprintf(message,"Sharp-Smooth angle is negative.");
	errors->Add_Error_Message(message);
	return INVALID_DATA;
      }
      sharp_smooth_curvature = ComputeCurvatureFromAngle(orig_sharp_smooth_angle);
    }
    break;
  case NORMAL_SMOOTHING:
    normal_smoothing_status = status;
    if( status == ACTIVE ){
      PRECONDITION( data );
      orig_sharp_smooth_angle = data[0];
      if( orig_sharp_smooth_angle < 0.0 ){
	std::sprintf(message,"Sharp-Smooth angle is negative.");
	errors->Add_Error_Message(message);
	return INVALID_DATA;
      }
      sharp_smooth_curvature = ComputeCurvatureFromAngle(orig_sharp_smooth_angle);
      normal_smoothing_distance = data[1];
      smoothing_resolution = (ContactSearch::Smoothing_Resolution ) data[2];
      int error = 0;
      if( normal_smoothing_distance > 1.0 || normal_smoothing_distance <= 0.0 ){
	error = 1;
	std::sprintf(message,"Normal Smoothing Distance is out of bounds.");
	errors->Add_Error_Message(message);
      }
      if( error == 1 ) return INVALID_DATA;
    }
    break;
    
  case COMPUTE_NODE_AREAS:
    compute_node_areas = status;
    break;
  case PARTITION_GAP:
    partition_gap_status = status;
    break;
  case NEW_TIED_ENFORCEMENT:
    new_tied_enforcement = status;
    break;
  case OLD_DYNAMIC_SEARCH:
    old_dynamic_search = status;
    break;
  case ENABLE_TRACKING:
    enable_tracking = status;
    if( status == ACTIVE ){
      PRECONDITION( data );
      tracking_type            = (Track_Type) data[0];
      global_tracking_interval = (int) data[1];
      if (tracking_type==LOCAL_TRACKING) {
        enable_off_face_tracking = (int) data[2];
      } else {
        enable_off_face_tracking = 0;
        no_secondary             = ACTIVE;
      }
    } else {
      tracking_type            = NO_TRACKING;
      global_tracking_interval = 0;
      enable_off_face_tracking = 0;
    }
    break;
  case NO_SECONDARY:
    no_secondary = status;
    break;
  case NO_GHOSTING:
    no_ghosting = status;
    break;
  case GLOBAL_SEARCH_CULL:
    search_cull = status;
    if( status == ACTIVE ){
      PRECONDITION( data );
      global_search_cull = (Search_Cull) data[0];
    } else {
      global_search_cull = NO_CULL;
    }
    break;
  case NO_WARPED_VOLUME:
    no_warped_volume = status;
    break;
  case NO_PARALLEL_CONSISTENCY:
    no_parallel_consistency = status;
    if( status == ACTIVE ){
      primary_topology->RebuildTopologyEntityList(no_parallel_consistency);
    }
    break;
  case AUTO_TOL:
    auto_tol = status;
    box_inflation = BOX_INFLATION_FACTOR;
    gap_inflation = GAP_INFLATION_FACTOR+1.0;
// this doesn't work with the fortran Search_Interface
//    if (status == ACTIVE && data) {
//       box_inflation = data[0]*BOX_INFLATION_FACTOR;
//       gap_inflation = data[1]*GAP_INFLATION_FACTOR+1.0;
//    }
    break;
  case AGGRESSIVE_TOLERANCES:
    aggressive_tolerances = status;
    if (status == ACTIVE) {
       box_inflation = 0.1*BOX_INFLATION_FACTOR;
       gap_inflation = 0.1*GAP_INFLATION_FACTOR+1.0;
    }
    else if (data) {
       box_inflation = data[0]*BOX_INFLATION_FACTOR;
       gap_inflation = data[1]*GAP_INFLATION_FACTOR+1.0;
    }
    break;
  case SKIP_PHYSICAL_FACES:
    skip_physical_faces = status;
    if( status == INACTIVE ){
      physical_face_algorithm = (PF_Algorithm) data[0];
    } else {
      physical_face_algorithm = PF_NONE;
    }
    break;
  case EDGE_PHYSICAL_FACES:
    edge_physical_faces = status;
    if( status == ACTIVE ){
      physical_face_algorithm = PF_EDGE_BASED;
    } else {
      physical_face_algorithm = PF_FACE_BASED;
    }
    break;
  case SHELL_SIMPLE_LOFTING:
    shell_simple_lofting = status;
    break;
  case COMPUTE_PARTIALS:
    compute_partials = status;
    if( status == ACTIVE ){
      computed_partials_order = (int) data[0];
    } else {
      computed_partials_order = 0;
    }
  default:
    POSTCONDITION(0);
    break;
  }
  return NO_ERROR;
}

void
ContactSearch::Get_Search_Option( Search_Option option, 
                                  Search_Option_Status& status,
				  Real* data )
{
  status = INACTIVE;
  
  switch( option ){
  case MULTIPLE_INTERACTIONS:
    status = multiple_interaction_status;
    if( status == ACTIVE && data!=NULL){
      PRECONDITION( data );
      data[0] = orig_sharp_smooth_angle;
    }
    break;
  case NORMAL_SMOOTHING:
    status = normal_smoothing_status;
    if( status == ACTIVE && data!=NULL ){
      PRECONDITION( data );
      data[0] = orig_sharp_smooth_angle;
      data[1] = normal_smoothing_distance;
      data[2] = smoothing_resolution;
    }
    break;
  case COMPUTE_NODE_AREAS:
    status = compute_node_areas;
    break;
  case PARTITION_GAP:
    status = partition_gap_status;
    break;
  case NEW_TIED_ENFORCEMENT:
    status = new_tied_enforcement;
    break;
  case OLD_DYNAMIC_SEARCH:
    status = old_dynamic_search;
    break;
  case ENABLE_TRACKING:
    status = enable_tracking;
    if( status == ACTIVE && data!=NULL ){
      PRECONDITION( data );
      data[0] = tracking_type;
      data[1] = global_tracking_interval;
      data[2] = enable_off_face_tracking;
    }
    break;
  case NO_SECONDARY:
    status = no_secondary;
    break;
  case NO_GHOSTING:
    status = no_ghosting;
    break;
  case GLOBAL_SEARCH_CULL:
    status = search_cull;
    if( status == ACTIVE && data!=NULL ){
      PRECONDITION( data );
      data[0] = global_search_cull;
    }
    break;
  case NO_WARPED_VOLUME:
    status = no_warped_volume;
    break;
  case NO_PARALLEL_CONSISTENCY:
    status = no_parallel_consistency;
    break;
  case AUTO_TOL:
    status = auto_tol;
    if (status == ACTIVE && data!=NULL) {
      PRECONDITION( data );
      data[0] = box_inflation;
      data[1] = gap_inflation;
    }
    break;
  case AGGRESSIVE_TOLERANCES:
    status = aggressive_tolerances;
    if (status == ACTIVE && data!=NULL) {
      PRECONDITION( data );
      data[0] = box_inflation;
      data[1] = gap_inflation;
    }
    break;
  case SKIP_PHYSICAL_FACES:
    status = skip_physical_faces;
    if( status == INACTIVE && data!=NULL ){
      PRECONDITION( data );
      data[0] = physical_face_algorithm;
    }
    break;
  case EDGE_PHYSICAL_FACES:
    status = edge_physical_faces;
    if( status == INACTIVE && data!=NULL ){
      PRECONDITION( data );
      data[0] = physical_face_algorithm;
    }
    break;
  case SHELL_SIMPLE_LOFTING:
    status = shell_simple_lofting;
    break;
  case COMPUTE_PARTIALS:
    status = compute_partials;
    if( status == ACTIVE && data!=NULL ){
      PRECONDITION( data );
      data[0] = computed_partials_order;
    }
    break;
  default:
    POSTCONDITION(0);
    status = INACTIVE;
    break;
  }
}

void ContactSearch::Register_Enforcement( ContactEnforcement* enf )
{
  if( number_registered_enforcements ){
    ContactEnforcement** temp = enforcement;
    enforcement = new ContactEnforcement*[number_registered_enforcements+1];
    std::memcpy( enforcement, temp, 
	    number_registered_enforcements*sizeof(ContactEnforcement*) );
    delete [] temp;
  } else
    enforcement = new ContactEnforcement*[1];
  enforcement[number_registered_enforcements++] = enf;
}

void ContactSearch::Delete_Enforcement( ContactEnforcement* enf )
{
  PRECONDITION (NULL != enf);
  PRECONDITION (number_registered_enforcements > 0);

  if ( number_registered_enforcements > 1) {
    ContactEnforcement** temp = enforcement;
    enforcement = new ContactEnforcement*[number_registered_enforcements-1];
    int count = 0;
    for (int i =0; i < number_registered_enforcements; ++i) {
      if ( temp[i] == enf ) continue;
      enforcement[count] = temp[i];
      ++count;
    }
    number_registered_enforcements --;
    POSTCONDITION ( count == number_registered_enforcements);
    delete [] temp;
  } else {
    number_registered_enforcements = 0;
  }
}

ContactSearch::ContactErrorCode
ContactSearch::Exodus_Output( int Exodus_ID, Real Time )
{
  ContactEnforcement* enf = NULL;
  if( number_registered_enforcements ) enf = enforcement[0];
  return primary_topology->Exodus_Output(Exodus_ID, 
					 Time, 
					 errors, 
					 search_data,
					 normal_smoothing_status,
					 multiple_interaction_status,
					 orig_sharp_smooth_angle, 
					 normal_smoothing_distance, 
					 smoothing_resolution,
					 compute_node_areas,
					 enf );
}

ContactSearch::ContactErrorCode
ContactSearch::Set_Node_Block_Remaining_Gap( int node_blk_id,
					     const Real* gap)
{
  ContactErrorCode error = 
    primary_topology->Set_NodeBlk_RemainingGap( node_blk_id, gap);
  return error;
}

ContactSearch::ContactErrorCode
ContactSearch::Set_Node_Block_Ghosting_Gap( int node_blk_id,
					    const Real* gap)
{
  ContactErrorCode error = 
    primary_topology->Set_NodeBlk_GhostingGap( node_blk_id, gap);
  return error;
}

ContactSearch::ContactErrorCode
ContactSearch::Set_Node_Block_Kinematic_Constraints( int node_blk_id,
					     const int* constraints_per_node,
					     const Real* constraint_vector)
{
  ContactErrorCode error = 
    primary_topology->Set_NodeBlk_KinConstr( node_blk_id, 
                                             constraints_per_node, 
                                             constraint_vector);
  return error;
}



ContactSearch::ContactErrorCode 
ContactSearch::Set_Node_Block_Configuration(ContactNode_Configuration config,
					    int node_blk_id, 
					    const Real* positions)
{
  ContactErrorCode error = NO_ERROR;
  switch (config) {
  case (CURRENT_CONFIG):
    error = primary_topology->Set_NodeBlk_Positions( node_blk_id, positions );
    break;
  case (PREDICTED_CONFIG):
    error = primary_topology->Set_NodeBlk_Positions_2( node_blk_id, positions );
    break;
  }
  return error;
}

ContactSearch::ContactErrorCode 
ContactSearch::Set_Node_Block_Attributes(Node_Block_Attribute attr,
					 int node_blk_id, 
					 const Real* attributes)
{
  ContactErrorCode error;
  error = primary_topology->Set_NodeBlk_Attributes( attr, node_blk_id, 
						    attributes );
  return error;
}

ContactSearch::ContactErrorCode
ContactSearch::Set_Face_Block_Attributes( Face_Block_Attribute attr,
					  int face_block_id,
					  const Real* attributes )
{
  ContactErrorCode error;
  error = primary_topology->Set_FaceBlk_Attributes( attr, face_block_id,
						    attributes );
  return error;
}


ContactSearch::ContactErrorCode 
ContactSearch::Add_Analytic_Surface(AnalyticSurface_Type surface_type,
                                    const Real* data)
{
  int surface_key = primary_topology->Number_of_Added_Analytic_Surfaces() +
                    primary_topology->Number_of_Node_Blocks() - 1 +
                    primary_topology->Number_of_Face_Blocks();
  int surface_id  = primary_topology->Number_of_Added_Analytic_Surfaces();
  ContactAnalyticSurface* surface;
  switch( surface_type ){
  case( PLANE ) :
    surface = new ContactAnalyticPlane( surface_id, surface_key, data );
    primary_topology->Add_Analytic_Surface( surface );
    break;
  case( SPHERE ) :
    surface = new ContactAnalyticSphere( surface_id, surface_key, data );
    primary_topology->Add_Analytic_Surface( surface );    
    break;
  case( CYLINDER_INSIDE ) :
    surface = new ContactAnalyticCylinderInside( surface_id, surface_key, data );
    primary_topology->Add_Analytic_Surface( surface );
    break;
  case( CYLINDER_OUTSIDE ) :
    surface = new ContactAnalyticCylinderOutside( surface_id, surface_key, data );
    primary_topology->Add_Analytic_Surface( surface );
    break;
  default:
    return UNKNOWN_TYPE;
  }
  return surface->Check_for_Errors( errors );
}

ContactSearch::ContactErrorCode
ContactSearch::Set_Analytic_Surface_Configuration( int id, const Real* data )
{
  return primary_topology->Set_Analytic_Surface_Configuration( id, data );
}

void ContactSearch::Set_Search_Data( const Real* data )
{
  search_data->Set_Search_Data( data );
}

void ContactSearch::Set_Search_Data( Search_Data_Index srch_index,
                                     int node_key, int face_key, Real value)
{
  search_data->Set_Search_Data(srch_index, node_key, face_key, value);
}

ContactSearch::ContactErrorCode ContactSearch::Delete_All_Interactions(){ 
  primary_topology->Delete_All_Interactions(); 
  step_number         = 0;
  tracking_step       = 0;
  initialized         = false;
  initialized_tied    = false;
  initialized_context = false;
  initializing_tied   = false;
  return NO_ERROR; 
};


ContactSearch::ContactErrorCode 
ContactSearch::Add_Table( int ID, int Num_Points, 
			  Real* abscissa, Real* ordinate )
{
  ContactErrorCode ec = NO_ERROR;

  // Since we don't expect to have tables in most problems and would have
  // a small number in any case, I'll just reallocate the array each time
  // one is added
  if( ID <= 0 ){
    ec = INVALID_ID;
    errors->Add_Error_Message("Table ID must be positive");
  }
  if( Num_Points <= 0 ){
    ec = INVALID_DATA;
    errors->Add_Error_Message("Number of points for a table must be positive");
  }
  PRECONDITION( abscissa && ordinate );
  if( ec == NO_ERROR ){
    ContactTable** temp = tables;
    tables = new ContactTable*[num_tables+1];
    if( num_tables > 0 ){
      std::memcpy( tables, temp, num_tables*sizeof(ContactTable*) );
      delete [] temp;
    }
    tables[num_tables++] = new ContactTable(ID,Num_Points,abscissa,ordinate );
  }
  return ec;
}

const ContactTable* ContactSearch::Get_Table( int ID )
{
  ContactTable* table = NULL;
  for( int i=0 ; i<num_tables ; ++i ){
    if( tables[i]->ID() == ID ){
      table = tables[i];
      break;
    }
  }
  return table;
}

ContactSearch::ContactErrorCode
ContactSearch::Interaction_Definition( int num_configs, ContactSearch::Topology use_topology,
                                       ContactTopologyEntity<Real>::SearchContext status )
{
  ContactTopology *topology(0);
  switch(use_topology){
    case ContactSearch::PRIMARY:
      topology = primary_topology;
      break;
    case ContactSearch::SECONDARY :
      topology = secondary_topology;
      break;
    case ContactSearch::SEARCH :
      topology = search_topology;
      break;
    default:
      PRECONDITION(0);
      break;
  }

  // Loop over all nodes and check if they have interactions and 
  // have the correct status.  If so, process them.

  // Get the variable handles I need
  VariableHandle CURRENT_POSITION =
    topology->Variable_Handle( ContactTopology::Current_Position );
  VariableHandle PREDICTED_POSITION = 
    topology->Variable_Handle( ContactTopology::Predicted_Position );
  VariableHandle AUGMENTED_POSITION = 
    topology->Variable_Handle( ContactTopology::Augmented_Position );
  VariableHandle NODE_NORMAL =
    topology->Variable_Handle( ContactTopology::Node_Normal );
  VariableHandle FACE_NORMAL =
    topology->Variable_Handle( ContactTopology::Face_Normal );
  VariableHandle CURVATURE = 
    topology->Variable_Handle( ContactTopology::Curvature );
  VariableHandle NUM_KIN_CONSTR =
    topology->Variable_Handle( ContactTopology::Num_Kin_Constr );
  VariableHandle KIN_CONSTR_VECTOR =
    topology->Variable_Handle( ContactTopology::Kin_Constr_Vector );
  VariableHandle REMAINING_GAP = 
    topology->Variable_Handle( ContactTopology::Remaining_Gap );

  Real tolmin = time_to_contact_tolerance_min;
  Real tolmax = time_to_contact_tolerance_max;
  Real time_to_contact;

  bool need_to_compact_nodeface_interactions = false;

  int i,j,kk;
  PRECONDITION( Max_Interactions_Per_Node() == 3 );
  ContactNodeEntityInteraction* nei_list[3];
  int number_of_nodes = topology->Number_of_Nodes();
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());
  for (i=0; i<number_of_nodes; ++i) {
      ContactNode<Real>* node = Nodes[i];
#ifdef CONTACT_DEBUG_NODE
      bool PRINT_THIS_NODE = primary_topology->Is_a_Debug_Node( node );
#endif
      //if this node is not locally owned then skip it.
      if (node->Ownership() != ContactTopologyEntity<Real>::OWNED) continue;
      
      //if this node doesn't have any node/entity interactions then skip it.
      if(node->Number_NodeEntity_Interactions() == 0) continue;
      
      int num_interactions = 0;
      ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions();
      if (initializing_tied) {
        for (j=0; j<node->Number_NodeEntity_Interactions(); ++j) {
          nei_list[num_interactions++] = interactions[j];
        }
      } else {
        for (j=0; j<node->Number_NodeEntity_Interactions(); ++j) {
          ContactNodeEntityInteraction* cnei = interactions[j];
          if( !(cnei->Is_Tied() || cnei->Is_InfSlip()) ){
            nei_list[num_interactions++] = cnei;
          }
        }
      }
      for( j=num_interactions; j<Max_Interactions_Per_Node(); ++j ){
        nei_list[j] = NULL;
      }

      if (num_interactions==0) continue;

      for( j=0 ; j<num_interactions ; ++j ){
        ContactNodeEntityInteraction* cnei =  nei_list[j];

        Real gap = cnei->Get_Gap_Cur();
#ifdef CONTACT_DEBUG_NODE
	if( PRINT_THIS_NODE ){
          if(cnei->Get_Type() ==ContactNodeEntityInteraction::NODE_FACE_INTERACTION ) {
            ContactNodeFaceInteraction*  cnfi = 
              static_cast<ContactNodeFaceInteraction*>(cnei);
            postream << "Processing Interaction_Definition for Node ("
                     << cnfi->Node()->Exodus_ID() << ") with Face "
                     << cnfi->Face()->Global_ID() << "\n";
	  } else if(cnei->Get_Type() ==ContactNodeEntityInteraction::NODE_SURFACE_INTERACTION){
            ContactNodeSurfaceInteraction* cnsi = 
             static_cast<ContactNodeSurfaceInteraction*>(cnei);
            postream << "Processing Interaction_Definition for Node ("
                     << cnsi->Node()->Exodus_ID() << ") with Analytic surface "
                     << cnsi->SurfaceID() << "\n";
          }
          postream << "      gap    = " << gap << "\n";
          Real* pbdir = cnei->Get_Pushback_Dir();
          postream << "      pbdir  = " << pbdir[0] << " " << pbdir[1]
                   << " " << pbdir[2] << "\n";
          Real* coords = cnei->Vector_Var(ContactNodeFaceInteraction::COORDINATES);
          postream << "      coords = " << coords[0] << " " 
                   << coords[1] << " " << coords[2] << "\n";
	}
#endif
        bool valid_interaction = true;
        int entity_key = cnei->Get_Entity_Key();
	int node_key   = cnei->Get_Node_Key();
        Real gap_tol = search_data->Get_Search_Data(SEARCH_NORMAL_TOLERANCE, node_key, entity_key);
        if (auto_tol) {
          Real node_reasonable_gap = max_node_displacement;
          Real *rem_gap = node->Variable(REMAINING_GAP);
          node_reasonable_gap += Magnitude(rem_gap);
          node_reasonable_gap *= gap_inflation;
          gap_tol = std::max(node_reasonable_gap,gap_tol);
        }
	ContactNodeEntityInteraction::InteractionSource source = cnei->Get_Source();  
        switch( source ){
	case ContactNodeEntityInteraction::UNKNOWN_SOURCE:
	  POSTCONDITION( 0 );
	  break;
        case ContactNodeEntityInteraction::CLOSEST_POINT_PROJECTION_1:
        case ContactNodeEntityInteraction::CLOSEST_POINT_PROJECTION_2:
          if( gap > gap_tol ){
	    valid_interaction = false;
#ifdef CONTACT_DEBUG_NODE
	    if( PRINT_THIS_NODE )
	      postream << "  Throwing OUT Interaction: " << gap << " > "
		       << gap_tol << "\n";
#endif
	  }
          break;
	case ContactNodeEntityInteraction::CLOSEST_POINT_PROJECTION_2_FROM_MOVING:
	  if( gap > gap_tol || gap < -gap_tol ){
	    valid_interaction = false;
#ifdef CONTACT_DEBUG_NODE
	    if( PRINT_THIS_NODE ){
	      if( gap > gap_tol )
		postream << "  Throwing OUT Interaction (CPP from Moving): " 
			 << gap << " > " << gap_tol << "\n";
	      else 
		postream << "  Throwing OUT Interaction (CPP from Moving): " 
			 << gap << " < "<< -gap_tol << "\n";
	    }
#endif
	  }
	  break;
        case ContactNodeEntityInteraction::MOVING_INTERSECTION:
          {
            time_to_contact = cnei->Get_Time_To_Contact();
            int tracking    = global_tracking_interval;
            Real time_tolmin = tolmin;
            Real time_tolmax = tolmax;
            if ((tracking_type==LOCAL_TRACKING) && tracking>1) {
              time_tolmin -= (tracking-1)*0.5;
              time_tolmax += (tracking-1)*0.5;
            }
            if( gap > gap_tol || 
                time_to_contact < time_tolmin || 
                time_to_contact > time_tolmax ){
	      valid_interaction = false;
#ifdef CONTACT_DEBUG_NODE
	      if( PRINT_THIS_NODE ){
	        postream << "  Throwing OUT Interaction: \n";
	        postream << "        Gap: " << gap << " > Toler: " 
	                 << gap_tol << "\n";
	        postream << "        Time: " << time_to_contact << " < Tolmin:"
	                 << time_tolmin << "\n";
	        postream << "        Time: " << time_to_contact << " > Tolmax:"
	                 << time_tolmax << "\n";   
	      }
#endif
	    }
	  }
          break;
        default:
          break;
        }

        if( valid_interaction ) {
#ifdef CONTACT_DEBUG_NODE
	  if( PRINT_THIS_NODE ) postream << "  Accepting Interaction\n";
#endif
          //There is no corellation to this for surfaces, 
          //so doing this only for the node face interactions
          if(cnei->Get_Type()==ContactNodeEntityInteraction::NODE_FACE_INTERACTION){
            ContactNodeFaceInteraction* cnfi = static_cast<ContactNodeFaceInteraction *>(cnei);
            ContactFace<Real>* face  = cnfi->Face();
            Real* Coordinates  = cnfi->Get_Coordinates();
            Real* Contactpoint = cnfi->Get_Contact_Point();
	    VariableHandle COORD_HANDLE;
            switch (source) {
            case ContactNodeEntityInteraction::CLOSEST_POINT_PROJECTION_1 :
	      COORD_HANDLE = CURRENT_POSITION;
              face->Compute_Local_Coordinates(0.0,CURRENT_POSITION, 
                    CURRENT_POSITION,FACE_NORMAL, Contactpoint, Coordinates );
              break;
            case ContactNodeEntityInteraction::CLOSEST_POINT_PROJECTION_2 :
	    case ContactNodeEntityInteraction::
	      CLOSEST_POINT_PROJECTION_2_FROM_MOVING :
	      COORD_HANDLE = PREDICTED_POSITION;
	      face->Compute_Local_Coordinates(1.0,CURRENT_POSITION, 
		    PREDICTED_POSITION,FACE_NORMAL, Contactpoint, Coordinates );
              break;
            case ContactNodeEntityInteraction::MOVING_INTERSECTION :
	      COORD_HANDLE = PREDICTED_POSITION;
	      face->Compute_Local_Coordinates(1.0,CURRENT_POSITION, 
                    PREDICTED_POSITION,FACE_NORMAL, Contactpoint, Coordinates );
              break;
            default:
              break;
            }
	    //
	    // BUG_KHP: Compare this code with the same code in
	    // BUG_KHP: ContactNodeFaceInteraction.
	    //

	    // Recompute the gap and the pushback direction based on the
	    // actual face using the shape functions.  For warped faces,
	    // this can be off by enough to cause problems.
	    //
	    // (KHB 5/8/03) This is causing problems for Adagio in the
	    // case where the node is "on a warped surface".  See the
	    // bug #10505.
	    Real* node_pp;
	    Real& gap_cur = cnfi->Get_Gap_Cur();
	    if( std::fabs(gap_cur) > 0.0  && 
		!(source  == ContactNodeFaceInteraction::RETRIEVED_TIED  ||
                  source  == ContactNodeFaceInteraction::RETRIEVED_GLUED )){
	      Real* n_dir = cnfi->Get_Normal();
	      cnfi->Face()->Compute_Normal(COORD_HANDLE, n_dir, 
					   cnfi->Get_Coordinates());
	      Real shape_functions[MAX_NODES_PER_FACE];
	      cnfi->Face()->Evaluate_Shape_Functions( 
						     cnfi->Get_Coordinates(),
						     shape_functions );
	      Contactpoint[0] = Contactpoint[1] = Contactpoint[2] = 0.0;
	      for( int ii=0 ; ii<cnfi->Face()->Nodes_Per_Face() ; ++ii ){
		Real* pp = cnfi->Face()->Node(ii)->Variable(COORD_HANDLE);
		Contactpoint[0] += shape_functions[ii]*pp[0];
		Contactpoint[1] += shape_functions[ii]*pp[1];
		Contactpoint[2] += shape_functions[ii]*pp[2];
	      }
	      Real pbdir_new[3];
	      node_pp       = cnfi->Node()->Variable(COORD_HANDLE);
	      pbdir_new[0]  = Contactpoint[0] - node_pp[0];
	      pbdir_new[1]  = Contactpoint[1] - node_pp[1];
	      pbdir_new[2]  = Contactpoint[2] - node_pp[2];

	      Real true_gap = Magnitude( pbdir_new);

	      Real* pb_dir  = cnfi->Get_Pushback_Dir();

	      Real pb_dot_n = Dot(pbdir_new,n_dir);

	      if( true_gap ){
		pbdir_new[0] /= true_gap;
		pbdir_new[1] /= true_gap;
		pbdir_new[2] /= true_gap;
		if( pb_dot_n > 0.0 ){
		  true_gap *= -1.0;
		} else {
		  pbdir_new[0] *= -1.0;
		  pbdir_new[1] *= -1.0;
		  pbdir_new[2] *= -1.0;
		}	      
	      } else {
		pbdir_new[0] = pb_dir[0];
		pbdir_new[1] = pb_dir[1];
		pbdir_new[2] = pb_dir[2];
	      }

#ifdef CONTACT_DEBUG_NODE
	      if( PRINT_THIS_NODE ){
		postream << "  Contactpoint = " << Contactpoint[0] << " "
			 << Contactpoint[1] << " " << Contactpoint[2] << "\n";
		postream << "  node_pp       = " << node_pp[0] << " " 
			 << node_pp[1] << " " << node_pp[2] << "\n";
		postream << "  Resetting Gap from " << gap_cur << " to " 
			 << true_gap << "\n";
		postream << "  Resetting PBDir from (" << pb_dir[0] << ","
			 << pb_dir[1] << "," << pb_dir[2] << ") to ("
			 << pbdir_new[0] << "," << pbdir_new[1] << ","
			 << pbdir_new[2] << ")\n";
	      }
#endif	      
	      gap_cur = true_gap;
	      pb_dir[0] = pbdir_new[0];
	      pb_dir[1] = pbdir_new[1];
	      pb_dir[2] = pbdir_new[2];
            }
	    if( !(cnfi->Is_Tied() || cnfi->Is_InfSlip() || cnfi->Is_Glued()) ){
              if (!((face->FaceType()==ContactSearch::LINEFACEL2) ||
                    (face->FaceType()==ContactSearch::LINEFACEQ3))) {
                //MWG: need to fix this, not yet implemented for line faces (2d)

	        if( num_configs == 3 )
	          cnfi->Modify_for_Curvature( AUGMENTED_POSITION,
	                                      FACE_NORMAL,
	                                      CURVATURE,
                                              sharp_smooth_curvature,
                                              multiple_interaction_status );
	        else
	          cnfi->Modify_for_Curvature( CURRENT_POSITION,
	                                      FACE_NORMAL,
	                                      CURVATURE,
                                              sharp_smooth_curvature,
                                              multiple_interaction_status );

	        if( normal_smoothing_status == ACTIVE )
	          cnfi->Modify_for_Normal_Smoothing( CURRENT_POSITION,
	                                             NODE_NORMAL,
	                                             FACE_NORMAL,
	                                             CURVATURE,
	                                             smoothing_resolution,
	                                             normal_smoothing_distance,
                                                     sharp_smooth_curvature);
              }
	    }
          }

	  // Partition the gap between g_cur*dt and gap_old
	  if( num_configs != 1 && partition_gap_status == ACTIVE &&
              !(cnei->Is_Tied() || cnei->Is_InfSlip() || cnei->Is_Glued()) ){
	    Partition_Between_GapCur_and_GapOld( cnei );
	  } else {
            cnei->Get_Gap_Old() = 0;
	  }

          bool Would_Violate_Constraints;
	  cnei->Modify_for_Kinematic_Constraints( NUM_KIN_CONSTR,
						  KIN_CONSTR_VECTOR,
						  Would_Violate_Constraints);
	  if( Would_Violate_Constraints ){
#ifdef CONTACT_DEBUG_NODE
	    if( PRINT_THIS_NODE ) postream << 
     "  Throwing out Interaction: Would Violate Kinematics Constraints\n";
#endif
	    need_to_compact_nodeface_interactions = true;
	    cnei->Node()->Delete_NodeEntity_Interaction( cnei ); 
	    nei_list[j] = NULL;
	    continue;
	  }

          if( compute_node_areas == ACTIVE && status==ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE ){
            if (!(source == ContactNodeFaceInteraction::RETRIEVED_TIED || source == ContactNodeFaceInteraction::RETRIEVED_GLUED) ){
              Real& slave_node_area = cnei->Get_Node_Area();
              slave_node_area = 0.0;
              int num_physical_faces = (int) *(NUMBER_PHYSICAL_FACES.Get_Scratch(i));
              if( num_physical_faces > 0 ){
                int cur_physical_face = 0;
                Real* pf_normal = cnei->Get_Physical_Face_Normal();
                Real* pf_normals[3];
                if( num_physical_faces == 1){
                  pf_normals[0] = PHYSICAL_FACE_NORMAL_1.Get_Scratch(i);
                  pf_normal[0]  = pf_normals[0][0];
                  pf_normal[1]  = pf_normals[0][1];
                  pf_normal[2]  = pf_normals[0][2];
                }
                else{
                  pf_normals[0] = PHYSICAL_FACE_NORMAL_1.Get_Scratch(i);
                  pf_normals[1] = PHYSICAL_FACE_NORMAL_2.Get_Scratch(i);
                  pf_normals[2] = PHYSICAL_FACE_NORMAL_3.Get_Scratch(i);
                  for(int k = 0; k < num_physical_faces; ++k){
                    if(pf_normals[k][0] == pf_normal[0] &&
                       pf_normals[k][1] == pf_normal[1] &&
                       pf_normals[k][2] == pf_normal[2] ){
                      cur_physical_face = k;
                    }
                  }
                }
                int num_connected_faces = node->Number_Face_Connections();

	        for( int iface=0 ; iface<num_connected_faces ; ++iface ){
                  if(node->GetFacePFIndex(iface) != cur_physical_face) continue;
		  ContactFace<Real>* Face = node->GetFace(iface);
		  Real slave_node_areas[MAX_NODES_PER_FACE];
		  PRECONDITION( Face->Nodes_Per_Face() <= MAX_NODES_PER_FACE );
		  VariableHandle POSITION = CURRENT_POSITION;
		  if( num_configs == 3 ) POSITION = AUGMENTED_POSITION;
		  Face->Compute_Node_Areas( POSITION, FACE_NORMAL,
					    &(slave_node_areas[0]) );
		  int inode;
		  Real* ms_face_normal = cnei->Get_Normal();
		  Real* sl_face_normal = Face->Variable(FACE_NORMAL);
		  Real projection_factor = 
		    ms_face_normal[0]*sl_face_normal[0] +
		    ms_face_normal[1]*sl_face_normal[1] +
		    ms_face_normal[2]*sl_face_normal[2] ;
		  projection_factor *= -1.0;
		  for( inode=0 ; inode<Face->Nodes_Per_Face() ; ++inode ){
		    if( Face->Node(inode) == node ){
		      slave_node_area += 
		        projection_factor*slave_node_areas[inode];
		      break;
		    }
		  }
		  POSTCONDITION( inode<Face->Nodes_Per_Face() );
	        }
              }
            }
          }
        } else {
          // Remove the interaction
	  need_to_compact_nodeface_interactions = true;
          cnei->Node()->Delete_NodeEntity_Interaction( cnei );
	  nei_list[j] = NULL;
        }

      }//end for num interactions
      if( num_interactions > 1 ) {
      // Now check for the case of co-linear constraints.  
      // Compact the list
	int num_int = 0;
	bool colinear[4];
	bool any_colinear = false;
	for( j=0 ; j<num_interactions ; ++j )
	  if( nei_list[j] ) nei_list[num_int++] = nei_list[j];
	for( j=0 ; j<num_int-1 ; ++j ){
	  if( nei_list[j] ){
	    colinear[j] = false;
	    PRECONDITION( nei_list[j] );
	    Real* pb_dir = nei_list[j]->Get_Pushback_Dir();
	    for( kk=j+1 ; kk<num_int ; ++kk ){
	      if( nei_list[kk] ){
		Real* pbd = nei_list[kk]->Get_Pushback_Dir();
		Real dot = Dot(pb_dir, pbd);
		if( dot >= 0.99 ){
		  any_colinear = true;
		  colinear[j] = true;
		  colinear[kk] = true;
		} else {
		  colinear[kk] = false;
		}
	      }
	    }
	    if( any_colinear ){
	      // Modify the interactions that are co-linear to be along
	      // the surface normal 
	      for( kk=j ; kk<num_int ; ++kk ){
		if( colinear[kk] ){
		  PRECONDITION( nei_list[kk] );
		  pb_dir = nei_list[kk]->Get_Pushback_Dir();
		  Real& gap_cur = nei_list[kk]->Get_Gap_Cur();
		  Real& gap_old = nei_list[kk]->Get_Gap_Old();
		  Real total_gap = gap_cur + gap_old;
#ifdef CONTACT_DEBUG_NODE
		  if( PRINT_THIS_NODE ){
		    postream << "Modifying interaction " << kk 
			     << " because its colinear with other interactions" 
			     << "\n";
		    postream << "    Original Data:" << "\n";
		    postream << "       GapCur = " << gap_cur << "\n";
		    postream << "       GapOld = " << gap_old << "\n";
		    postream << "       PBDIR  = " << pb_dir[0] << " " 
			     << pb_dir[1] << " " << pb_dir[2] << "\n";
		  }
#endif
		  gap_old = 0.0;
                  //I'm not 100% certain that this next line is right, but I think it is
                  //BMP 8/31/2005
                  Real* f_normal = nei_list[kk]->Get_Normal();
		  Real dp = f_normal[0]*pb_dir[0] + f_normal[1]*pb_dir[1] +
  		            f_normal[2]*pb_dir[2];
		  gap_cur = dp*total_gap;
		  pb_dir[0] = f_normal[0];
		  pb_dir[1] = f_normal[1];
		  pb_dir[2] = f_normal[2];
		  if( num_configs > 1 && partition_gap_status == ACTIVE ) 	     
		    Partition_Between_GapCur_and_GapOld( nei_list[kk] );
#ifdef CONTACT_DEBUG_NODE
		  if( PRINT_THIS_NODE ){
		    postream << "    Original Data:" << "\n";
		    postream << "       GapCur = " 
			     <<  nei_list[kk]->Get_Gap_Cur()
			     << "\n";
		    postream << "       GapOld = " 
			     <<  nei_list[kk]->Get_Gap_Old()
			     << "\n";
		    postream << "       PBDIR = " << pb_dir[0] << " " 
			     << pb_dir[1] << " " << pb_dir[2] << "\n";
		  }
#endif
		}
	      }
	    }
	  }
	}
      }
  }
#ifdef CONTACT_DEBUG_NODE
  //postream.flush();
#endif
  return (NO_ERROR);
}

void 
ContactSearch::Process_Interaction( VariableHandle NODE_COORDS,
				    ContactNodeEntityInteraction* cnei,
				    int& action )
{
  // action = -1, interaction not added
  //        >= 0, interaction replaced/added return interaction number
  action = -1;

  ContactNode<Real>* node = cnei->Node();
#ifdef CONTACT_DEBUG_NODE
  bool PRINT_THIS_NODE = primary_topology->Is_a_Debug_Node( node );
  if (PRINT_THIS_NODE)
    postream<<"    Process_Interaction() for node "<<node->Global_ID()<<"\n";
#endif
  ContactNodeEntityInteraction* cnei_stored=NULL;
  switch( multiple_interaction_status ){
  case( INACTIVE ):{
    action = 0;
    cnei_stored = node->Get_NodeEntity_Interaction( 0 );
    if( !cnei_stored ){
      // no previous interactions so simply store this one
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE ) 
	postream << "    Storing Interaction (none stored)\n";
#endif
      node->Store_NodeEntity_Interaction( 0, cnei);
      return;
    } else {
      // compete them
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE )
	postream << "    Competing Against Stored Face \n"; 
#endif
      if( Compete_Interactions( NODE_COORDS, cnei_stored, cnei ) == cnei ){
	node->Store_NodeEntity_Interaction( 0, cnei);
#ifdef CONTACT_DEBUG_NODE
	if( PRINT_THIS_NODE ) postream << "      Chose candidate face\n";
#endif
      } else {
	cnei->Delete_Frag(allocators);
#ifdef CONTACT_DEBUG_NODE
	if( PRINT_THIS_NODE ) postream << "      Chose stored face\n";
#endif
      }
      return;
    }
  }
  case( ACTIVE ):{
    for( int i=0 ; i<Max_Interactions_Per_Node() ; ++i ){
      cnei_stored = node->Get_NodeEntity_Interaction( i );
      if( !cnei_stored ){
	// No interaction is stored in location i and we aren't connected to
	// previous interactions so simply store this one
#ifdef CONTACT_DEBUG_NODE
	if( PRINT_THIS_NODE ) 
	  postream << "    Storing Interaction (none stored)\n";
#endif
	node->Store_NodeEntity_Interaction( i, cnei);
        action = i;
	return;
      } else {
	// If interaction i is contacting the same physical face, these
	// two interactions should be competed, otherwise they should be
	// accumulated.
	Real* stored_pf_normal    = cnei_stored->Get_Physical_Face_Normal();
	Real* candidate_pf_normal = cnei->Get_Physical_Face_Normal();
	Real dot = Dot(stored_pf_normal, candidate_pf_normal);
	if( dot >= 0.99 ){
          action = i;
#ifdef CONTACT_DEBUG_NODE
	  if( PRINT_THIS_NODE )
	    postream << "    Competing Against Stored Face \n"; 
#endif
	  if( Compete_Interactions( NODE_COORDS, cnei_stored, cnei ) == cnei ){
	    node->Store_NodeEntity_Interaction( i, cnei);
#ifdef CONTACT_DEBUG_NODE
	    if( PRINT_THIS_NODE ) postream << "    Chose candidate face\n";
#endif
	  } else {
  	    cnei->Delete_Frag(allocators);
#ifdef CONTACT_DEBUG_NODE
	    if( PRINT_THIS_NODE ) postream << "    Chose stored face\n";
#endif
	  }
	  return;
	}
      }
    }
#if CONTACT_DEBUG_PRINT_LEVEL>=1
    int id = 0;
    if (node->Is_a_Shell_Node()) {
      id = static_cast<ContactShellNode*> (node)->Shell_Node_Base_ID();
    } else {
      id = node->Exodus_ID();
    }
    std::cerr << "P"<<contact_processor_number( SearchComm )
              << ": Exceeding the maximum number of interactions at node global ID: "
              << node->Global_ID()
              << " host ID: "
              << id
              << "\n";
    ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions();
    for( int k=0 ; k<Max_Interactions_Per_Node() ; ++k ){
      if (interactions[k]->Get_Type()==ContactNodeEntityInteraction::NODE_FACE_INTERACTION) {
        ContactNodeFaceInteraction* cnfi = static_cast<ContactNodeFaceInteraction *>(interactions[k]);
        std::cerr << "   Node entity key "<<cnfi->Get_Node_Key()<<" has interaction with Face "
                  <<cnfi->Face()->Global_ID()<<", entity key "<<cnfi->Face()->Entity_Key()<<"\n";
        Real* pfn = cnfi->Vector_Var(ContactNodeFaceInteraction::PHYSICAL_FACE_NORMAL);
        std::cerr << "     PF Normal = ("<<pfn[0]<<", "<<pfn[1]<<", "<<pfn[2]<<")\n";
      }
    }
    if (cnei->Get_Type()==ContactNodeEntityInteraction::NODE_FACE_INTERACTION) {
      ContactNodeFaceInteraction* cnfi_new = static_cast<ContactNodeFaceInteraction *>(cnei);
      std::cerr << "   Node entity key "<<cnfi_new->Get_Node_Key()<<" trying to add Face "
                <<cnfi_new->Face()->Global_ID()<<", entity key "<<cnfi_new->Face()->Entity_Key()<<"\n";
      Real* pfn_new = cnfi_new->Vector_Var(ContactNodeFaceInteraction::PHYSICAL_FACE_NORMAL);
      std::cerr << "     PF Normal = ("<<pfn_new[0]<<", "<<pfn_new[1]<<", "<<pfn_new[2]<<")\n";
      std::cerr << "\nOld interactions:\n";
      interactions = node->Get_NodeEntity_Interactions(1);
      for( int k=0 ; k<Max_Interactions_Per_Node() ; ++k ){
        if (interactions[k]->Get_Type()==ContactNodeEntityInteraction::NODE_FACE_INTERACTION) {
          ContactNodeFaceInteraction* cnfi = static_cast<ContactNodeFaceInteraction *>(interactions[k]);
          std::cerr << "   Node entity key "<<cnfi->Get_Node_Key()<<" has interaction with Face "
                    <<cnfi->Face()->Global_ID()<<", entity key "<<cnfi->Face()->Entity_Key()<<"\n";
          Real* pfn = cnfi->Vector_Var(ContactNodeFaceInteraction::PHYSICAL_FACE_NORMAL);
          std::cerr << "     PF Normal = ("<<pfn[0]<<", "<<pfn[1]<<", "<<pfn[2]<<")\n";
        }
      }
    }
              
#endif
    POSTCONDITION( 0 );
  }
  }
}

void ContactSearch::Display()
{
#if CONTACT_DEBUG_PRINT_LEVEL>=6
  if( primary_topology ) primary_topology->Display(postream);
  postream.flush();
#endif
}

ContactFixedSizeAllocator& ContactSearch::Get_Allocator( AllocatorType t )
{
  PRECONDITION(t < ALLOC_NUM_ALLOCATED_ENTITIES);
  return allocators[t];
}

void ContactSearch::Size_NodeNode_Interactions(
        int& num_interactions,  // number of interactions being returned
        int& size               // size of data array per interaction
)
{
  primary_topology->Size_NodeNode_Interactions( num_interactions, size );
}

void ContactSearch::Get_NodeNode_Interactions( int* slave_node_block_ids,
                                               int* slave_node_indexes_in_block,
                                               int* master_node_block_ids,  
                                               int* master_node_indexes_in_block,
                                               int* master_node_proc,
                                               Real* interaction_data
                                               )
{
  primary_topology->Get_NodeNode_Interactions( slave_node_block_ids, 
					       slave_node_indexes_in_block,
					       master_node_block_ids, 
					       master_node_indexes_in_block,
					       master_node_proc, 
					       interaction_data );
}

void ContactSearch::Size_NodeFace_Interactions(
        int& num_interactions,  // number of interactions being returned
        int& size               // size of data array per interaction
)
{
  primary_topology->Size_NodeFace_Interactions( num_interactions, size );
}

void ContactSearch::Get_NodeFace_Interactions( int* node_block_ids,
                                               int* node_indexes_in_block,
					       int* node_keys,
                                               int* face_block_ids,  
                                               int* face_indexes_in_block,
                                               int* face_proc,
                                               Real* interaction_data
                                               )
{
  primary_topology->Get_NodeFace_Interactions( node_block_ids, 
					       node_indexes_in_block,
					       node_keys,
					       face_block_ids, 
					       face_indexes_in_block,
					       face_proc, 
					       interaction_data );
}

void ContactSearch::Size_NodeSurface_Interactions(
        int& num_interactions,  // number of interactions being returned
        int& size               // size of data array per interaction
)
{
  primary_topology->Size_NodeSurface_Interactions(num_interactions, size);
}

void ContactSearch::Get_NodeSurface_Interactions( int* node_block_ids,
                                                  int* node_indexes_in_block,
                                                  int* surface_ids,  
                                                  Real* interaction_data )
{
  primary_topology->Get_NodeSurface_Interactions( node_block_ids, 
						  node_indexes_in_block,
						  surface_ids, 
						  interaction_data );
}

void ContactSearch::Size_FaceFace_Interactions(
        int& num_interactions,  // number of interactions being returned
        int& size               // size of whole data array
)
{
  primary_topology->Size_FaceFace_Interactions( num_interactions, size );
}

void ContactSearch::Get_FaceFace_Interactions( int* slave_face_block_ids,
                                               int* slave_face_indexes_in_block,
                                               int* slave_face_proc,
                                               int* master_face_block_ids,  
                                               int* master_face_indexes_in_block,
                                               int* master_face_proc,
                                               int* interaction_index,
                                               Real* interaction_data )
{
  primary_topology->Get_FaceFace_Interactions( slave_face_block_ids, 
					       slave_face_indexes_in_block,
					       slave_face_proc,
					       master_face_block_ids, 
					       master_face_indexes_in_block,
					       master_face_proc, 
                                               interaction_index,
					       interaction_data );
}

void ContactSearch::Size_FaceCoverage_Interactions(
        int& num_interactions,  // number of interactions being returned
        int& size               // size of whole data array
)
{
  primary_topology->Size_FaceCoverage_Interactions( num_interactions, 
                                                    size );
}

void ContactSearch::Get_FaceCoverage_Interactions( int* face_block_ids,
                                                   int* face_indexes_in_block,
                                                   int* interaction_index,
                                                   Real* interaction_data )
{
  primary_topology->Get_FaceCoverage_Interactions( face_block_ids, 
					           face_indexes_in_block,
                                                   interaction_index,
					           interaction_data );
}

void ContactSearch::Size_ElementElement_Interactions(
        int& num_interactions,  // number of interactions being returned
        int& size               // size of whole data array
)
{
  primary_topology->Size_ElementElement_Interactions( num_interactions, 
                                                      size );
}

void ContactSearch::Get_ElementElement_Interactions( int* slave_element_block_ids,
                                                     int* slave_element_indexes_in_block,
                                                     int* master_element_block_ids,  
                                                     int* master_element_indexes_in_block,
                                                     int* master_element_proc,
                                                     Real* interaction_data )
{
  primary_topology->Get_ElementElement_Interactions( slave_element_block_ids, 
					             slave_element_indexes_in_block,
					             master_element_block_ids, 
					             master_element_indexes_in_block,
					             master_element_proc, 
					             interaction_data );
}

int ContactSearch::Cull_Node_List_New( ACME::ContactNode_Vector &node_list, 
				       int* physical_faces,
				       ContactFace<Real>* face,
				       VariableHandle NODE_NORMAL,
				       VariableHandle FACE_NORMAL,
				       int number_of_configurations, 
                                       const ContactBoundingBox &current_face_box,
                                       const ContactBoundingBox &predicted_face_box,
                                       ACME::Int_Vector &node_keys,
                                       const vector<Real> &norm_search_tol,
                                       const vector<Real> &tang_search_tol,
                                       const vector<bool> &valid_inter,
                                       VariableHandle NODE_COORD_START,
                                       VariableHandle NODE_COORD_END)
{
  // This function culls down a node list for a given pair based on the 
  // following criteria
  //   1) Is the node owned by this processor
  //  2a) Is the node connected to the face
  //  2b) Is the node a shell node and is its "mating" node connected to this
  //      face  
  //  3a) Is the average node normal "opposed" to the face
  //       We want to do the stregnth of contact check this way because of
  //       the following possibility
  //                                      
  //           Physical Face Normal 1       
  //          / \          -              
  //           |          /| Node normal  
  //           |         /
  //      --------------
  //                   |
  //                   |
  //                   |----->  Physical Face Normal 2
  //                   |
  //                   |
  //
  //       Paired Face
  //          /
  //         |                           In this case, the node is not
  //        /                            opposed to the node and shouldn't
  //        |   -\                       be considered.  However, it is 
  //       /      -\                     opposed to physical face 1.
  //       |        ->  Face Normal
  //
  //  3b) Is the node a "tab" node
  //  3c) Is the face opposed to the 'most opposed' physical face.
  //   4) Does the user want this node interacting with the face
  //   5) Using a Search Box specific to the face-node pair
  //      (Note: We used the largest search tolerance to get the pairs)
  //   6) Is the node tied to this surface already
  //
  int i,j,ii;

  int nodes_per_face = face->Nodes_Per_Face();
  PRECONDITION( nodes_per_face <= MAX_NODES_PER_FACE );
  ContactNode<Real>** face_nodes = face->Nodes();
  Real* face_normal = face->Variable(FACE_NORMAL);

  // Get the variable handles I need
  VariableHandle REMAINING_GAP = search_topology->Variable_Handle( ContactTopology::Remaining_Gap );

  bool edge_normal_defined = false;
  Real edge_normals[MAX_NODES_PER_FACE][3];

  int pairs = 0;
  int list_size = node_list.size();
  for( i=0 ; i<list_size ; ++i ){
    ContactNode<Real> *node = node_list[i];

#ifdef CONTACT_DEBUG_NODE
    bool PRINT_THIS_NODE = primary_topology->Is_a_Debug_Node( node );
    if( PRINT_THIS_NODE ){
      VariableHandle CURRENT_POSITION =   
	search_topology->Variable_Handle( ContactTopology::Current_Position );
      VariableHandle PREDICTED_POSITION = 
	search_topology->Variable_Handle( ContactTopology::Predicted_Position);
      VariableHandle AUGMENTED_POSITION = 
	search_topology->Variable_Handle( ContactTopology::Augmented_Position);
      postream << "Debug Node ("
               << node->Exodus_ID() << ") is paired with face " 
	       << face->Global_ID() << "\n";
      Real* fnv = face->Variable(FACE_NORMAL);
      Real* nnv = node->Variable(NODE_NORMAL);
      postream << "      Face Normal = " << fnv[0] << " " 
	       << fnv[1] << " " << fnv[2] << "\n";
      for( int idbgn=0 ; idbgn<face->Nodes_Per_Face() ; idbgn++ ){
	ContactNode<Real>* nd = face->Node(idbgn); 
	int id;
	if (nd->Is_a_Shell_Node()) {
	  id = static_cast<ContactShellNode*> (nd)->Shell_Node_Base_ID();
        } else {
	  id = nd->Exodus_ID();
        }
	Real* c = nd->Variable(CURRENT_POSITION);
	postream << "        Node (" 
                 << nd->Exodus_ID() << ") -- host id ("<< id <<")\n";
	postream << "            Current:   "
		 << c[0] << " " << c[1] << " " << c[2] << "\n";
	if( number_of_configurations >= 2) { 
	  Real * pc = nd->Variable( PREDICTED_POSITION );
	  postream << "            Predicted: "
		   << pc[0] << " " << pc[1] << " " << pc[2] << "\n";
	}
	if( number_of_configurations == 3) { 
	  Real * ac = nd->Variable( AUGMENTED_POSITION );
	  postream << "            Augmented: "
		   << ac[0] << " " << ac[1] << " " << ac[2] << "\n";
	}
      }
      postream << "      Slave Node Normal = " << nnv[0] << " " 
	       << nnv[1] << " " << nnv[2] << "\n";
      Real* c = node->Variable(CURRENT_POSITION); 
      postream << "      Slave Node -- host id (" 
               << node->Exodus_ID() << ")\n";
      postream << "            Current:   "
               << c[0] << " " << c[1] << " " << c[2] << "\n";
      if( number_of_configurations >= 2) { 
        Real * pc = node->Variable( PREDICTED_POSITION );
        postream << "            Predicted: "
                 << pc[0] << " " << pc[1] << " " << pc[2] << "\n";
      }
      if( number_of_configurations == 3) { 
        Real * ac = node->Variable( AUGMENTED_POSITION );
        postream << "            Augmented: "
                 << ac[0] << " " << ac[1] << " " << ac[2] << "\n";
      }
    }
#endif

    //
    // Test 1: Cull out the nodes that make up the face
    //
    for( j=0 ; j<nodes_per_face ; ++j ){
      if( node == face_nodes[j] ){
        break;
      }
    }
    if(j<nodes_per_face) {
      continue;
    }

    if( node->Is_a_Shell_Node() && ( face->FaceType() == SHELLQUADFACEL4 || face->FaceType() == SHELLTRIFACEL3 ) ){
      //
      // Test 2b: Is this a shell node whose "mating" node is connected to
      //          this face?
      int shell_node_base_id =  (static_cast<ContactShellNode*>(node))->Shell_Node_Base_ID();
      for( j=0 ; j<nodes_per_face ; ++j ){
        int face_shell_node_base_id = (static_cast<ContactShellNode*>(face_nodes[j]))->Shell_Node_Base_ID();
	if (face_shell_node_base_id == shell_node_base_id) {
	  break;
	}
      }
      if(j<nodes_per_face) continue;      
    }
    //
    // Test 2: Is the node "opposed" to the face if it has a face normal
    //
    if( node->NodeType() != ContactSearch::POINT ){
      Real* node_normal = node->Variable(NODE_NORMAL);

      Real dot_product = Dot(face_normal, node_normal);

      if( dot_product>=0.0 ) {
#ifdef CONTACT_DEBUG_NODE
        if( PRINT_THIS_NODE ){
          postream<<"  Node culled, node normal does not oppose face normal\n";
        }
#endif
        continue;
      }
    }
    //
    // Test 3: Is the node "opposed" to the physical face 
    // Compute the physical face and node data key
    //
    int physical_face = -1;
    int node_key;
    if(node_keys.size() <= i) {
      node_key = node->Entity_Key();
      if(!valid_inter[node_key+1]) continue;
    } else {
      node_key = node_keys[i];
    }
    int num_physical_faces = (int) *(NUMBER_PHYSICAL_FACES.Get_Scratch(node->ProcArrayIndex()));
    if( num_physical_faces>0 ){
      Real most_opposed = 2.0;
      Real* pf_normals[3];
      pf_normals[0] = PHYSICAL_FACE_NORMAL_1.Get_Scratch(node->ProcArrayIndex());
      pf_normals[1] = PHYSICAL_FACE_NORMAL_2.Get_Scratch(node->ProcArrayIndex());
      pf_normals[2] = PHYSICAL_FACE_NORMAL_3.Get_Scratch(node->ProcArrayIndex());
      Find_Physical_Face(num_physical_faces, pf_normals, face_normal, physical_face, most_opposed);
      if( most_opposed< 0.0) {
        if( node_key == -1 ){
	  //
          // Look for the first valid face entity key in this physical face
          // Test 4: Does the user want this node interacting with this face
	  //
          bool found_valid = false;
          int num_connected_faces = node->Number_Face_Connections();
          for( ii=0 ; ii<num_connected_faces ; ++ii ){
            if(node->GetFacePFIndex(ii) != physical_face) continue;
            node_key = node->GetFace(ii)->Entity_Key();
            if( valid_inter[node_key + 1] ){
	      found_valid = true;
	      break;
	    }
	  }
  	  if( !found_valid ){

#ifdef CONTACT_DEBUG_NODE
            if( PRINT_THIS_NODE ){
              postream<<"  Node culled, face has multiple physical faces and no valid physical face found\n";
            }
#endif

	    continue;
	  }
	}
      } else {
#ifdef CONTACT_DEBUG_NODE
        if( PRINT_THIS_NODE ){
          postream<<"  Node culled, no valid physical face found\n";
        }
#endif
        continue;
      }
    }
    //
    //  Test 4:  Compute a temporal bounding hull for the face and determine if the node lies inside of
    //           that hull.  The bounding hull will ultimately take into account:
    //              1) Non-axis alligned edges and faces
    //              2) Rigid body motion of the face and node
    //              3) True node-face pair normal and tangential tolerances.
    //
    if (!edge_normal_defined) {
      //
      //  Define edge normal directions for the face.  This only needs to be done once for the face even if testing against multiple nodes.
      //  Only do this calculation if required, but if it is required save the results for the entire loop.
      //
      for(int iedge = 0; iedge < nodes_per_face; ++iedge) {
        Real *start_coord1 = face_nodes[iedge]->Variable(NODE_COORD_START);
        int edge_node2 = (iedge + 1)%nodes_per_face;
        Real *start_coord2 = face_nodes[edge_node2]->Variable(NODE_COORD_START);
        Real edge_vec[3];
        edge_vec[0] = start_coord2[0] - start_coord1[0];        
        edge_vec[1] = start_coord2[1] - start_coord1[1];        
        edge_vec[2] = start_coord2[2] - start_coord1[2];        
	//
	//  Compute the edge normal as the cross product of the edge vector and the face normal.
	//  Note, here we compute a non-unit vector.  The factor that the normal normal vector
	//  is not unit length will be taken into account later.  This is done to avoid sqrt
	//  calls.
	//
        Real *edge_normal = edge_normals[iedge];
	Cross(edge_vec, face_normal, edge_normal);
      }
      edge_normal_defined = true; 
    }
    //
    //  Consider a moving coordinate system that's location is defined at the contact node.  In this coordinate
    //  system the contact node will appear fixed, and will always sit at the coordinate system origin.  In this
    //  coordinate system the face will have two positions, a start time postion and an end time position.  
    //  Consider that the moving 2D face is swept into a 3D prisim over time.  If the node lies anywhere either
    //  inside of that prisim, or within tolerance of the edge of the prisim, then the node should be considered
    //  for contact.
    //
    //  The prisim is made of two planes with normal equal to the face normal, and num_edges planes that point
    //  out from the prism sides.  This following code chunk will compute how far away a the node is from each of
    //  these planes.  These distances are given by max_norm_plus, max_norm_minus, and edge_max_dist.  Need to
    //  set all values to very negative to get the max.
    //
    Real max_norm_plus  = -HUGE_VAL;
    Real max_norm_minus = -HUGE_VAL;
    Real edge_max_dist_times_edge_length[MAX_NODES_PER_FACE];
    for(int iedge = 0; iedge < nodes_per_face; ++iedge) {
      edge_max_dist_times_edge_length[iedge] = -HUGE_VAL;
    }
    //
    //  The prism has 2 times number of face nodes nodes.  One set of nodes in the start configuration and a 
    //  second set of nodes in the end configuration.  Loop over each of these nodes and determine what prisim
    //  plane they define.  Compute the most restritive set of prisim bounding planes that completely contain
    //  the prisim.  Note, prism planes are defined by the previosly calculated plane normals and the furthest
    //  nodal point along that normal.
    //
    Real *node_start = node->Variable(NODE_COORD_START);
    Real *node_end   = node->Variable(NODE_COORD_END);
    for(int inode = 0; inode < nodes_per_face; ++inode) {
      //
      //  Define prisim nodes in the moving coordinate system where the moving node is the origin.
      //
      Real *start_coord = face_nodes[inode]->Variable(NODE_COORD_START);
      Real *end_coord   = face_nodes[inode]->Variable(NODE_COORD_END);

      Real adjusted_start[3] = {start_coord[0] - node_start[0],
                                start_coord[1] - node_start[1],
                                start_coord[2] - node_start[2]};

      Real adjusted_end[3]   = {end_coord[0]   - node_end[0],
                                end_coord[1]   - node_end[1],
                                end_coord[2]   - node_end[2]};
      //
      //  Define two planes in the +/- face normal direction to encompass the prism, norm_plus and
      //  norm_minus.  Find the maximal values in norm_plus and norm_minus
      //    
      Real norm_plus_start = Dot(adjusted_start, face_normal);
      Real norm_plus_end =   Dot(adjusted_end, face_normal);

      max_norm_plus  = max(max_norm_plus , norm_plus_start);
      max_norm_plus  = max(max_norm_plus , norm_plus_end);
      max_norm_minus = max(max_norm_minus, -norm_plus_start);
      max_norm_minus = max(max_norm_minus, -norm_plus_end  );
      //
      //  Define nodes_per_face edge planes of the prism.  Find the maximal nodal values in
      //  each plane. edge_normal is not a unit vector, so really here we are calculating the
      //  maximum distance from the plane times the edge length.  This will be accounted for later,
      //  again this is all to avoid expensive sqrt calls.
      //
      for(int iedge = 0; iedge < nodes_per_face; ++iedge) {
        Real *edge_normal = edge_normals[iedge];
        edge_max_dist_times_edge_length[iedge] = max(Dot(adjusted_start,edge_normal), edge_max_dist_times_edge_length[iedge]);
        edge_max_dist_times_edge_length[iedge] = max(Dot(adjusted_end,  edge_normal), edge_max_dist_times_edge_length[iedge]);
      }
    }
    //
    //  Determine the tolerance at which a point outside the bounding face hull is considered close enough to the
    //  face to be included in the interaction definitions. 
    //  NKC, note, not sure I really need to add gap to tangential tolerance, but it seems much safer to do so.
    //
    Real tolerance_norm_dist;
    Real tolerance_tang_dist;
    if (auto_tol) {
      Real max_box_dimension = current_face_box.max_dimension();
      tolerance_norm_dist = std::max(box_inflation*max_box_dimension, norm_search_tol[node_key]);
      tolerance_tang_dist = std::max(box_inflation*max_box_dimension, tang_search_tol[node_key]);
    } else {
      tolerance_norm_dist = norm_search_tol[node_key];
      tolerance_tang_dist = tang_search_tol[node_key];
    } 
    Real *gap_vector = node->Variable(REMAINING_GAP);
    Real gap = Magnitude(gap_vector);
    tolerance_norm_dist += gap;
    tolerance_tang_dist += gap;
    //
    //  Determine if the contact point (always at (0,0,0)) lies outside of any of the defined planes.  If it is
    //  outside of any plane, throw out the interaction immediatly, no need to check other planes.
    //
    bool inside = false;
    while(true) {
      //
      //  Normal tolerances simple, if the plane is futher than tolerance throw out the interactions
      //
      if(-max_norm_plus > tolerance_norm_dist) break;
      if(-max_norm_minus > tolerance_norm_dist) break;
      int iedge;
      for(iedge = 0; iedge < nodes_per_face; ++iedge) {
	//
	//  Edge tolerances fancy.  Need to take into account the fact that edge_normal are not unit vectors.
	//  The true calulation that we want is:
	//
	//  [1] if(-edge_max_dist_times_edge_length[iedge] / edge_length[iedge] > tolerance_tang_dist) break;
	//
	//  This is transformed as follows:
	//    [2] tolerance_tang_dist is always > 0.
	//    [3] so if edge_max_dist_times_edge_length[iedge] > 0 then [1] will always be false.  
	//        In this case we can skip the rest of the calculation
	//    [4] Given [3], then -(edge_max_dist_times_edge_length[iedge]/edge_length[iedge]) >= 0
	//    [5] If A > 0 && B > 0  then   (A > B) iff (A^2 > B^2)
	//    [6] Using [2], [4], and [5], we can square both sides of inequality in [1].
	//    [7] edge_length*edge_length = (sqrt(x*x + y*y + z*z))^2 = x*x + y*y + z*z
	//    [8] Plug in everything above to compare the squared inequality.  Thus, the sqrts are avoid, close
	//        to doubling the ultimate calculation speed.
	//
	//  Further optimization notes:
	//    Often nodes will lie unambigously inside of the prisim (edge_max_dist_times_edge_length[iedge]>0).  
	//    It is only for those nodes outside the prism that the below calculation operates.  Thus I do not
	//    precompute anything (such as the ||edge_normal||^2 or tolerance^2.
	//
        Real *edge_normal = edge_normals[iedge];
        if(edge_max_dist_times_edge_length[iedge]>0) {
	} else {
          if(((edge_max_dist_times_edge_length[iedge]*edge_max_dist_times_edge_length[iedge])/
             (edge_normal[0]*edge_normal[0] + edge_normal[1] * edge_normal[1] + edge_normal[2] * edge_normal[2])) > 
	     tolerance_tang_dist*tolerance_tang_dist) break;
	}
      }
      if(iedge != nodes_per_face) break;
      inside = true;
      break;
    }
    if(!inside) {
      //
      //  If the node is not inside the prisim, cull the interaction
      //
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE ){
        postream<<"  Node culled, node is outside non-axis aligned face bounding hull\n";
      }
#endif

      continue;
    }
    //
    // Test 5: Was the node previously "tied" to a face on this surface
    //         Look for interactions in the "old" state.
    //
    if( multiple_interaction_status == ACTIVE ){
        ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions(1);
        for (j=0; j<node->Number_NodeEntity_Interactions(1); ++j) {
          ContactNodeEntityInteraction* cnei = interactions[j];
          if( cnei->Is_Tied() || cnei->Is_InfSlip() || cnei->Is_Glued()){
            Real* pfs = cnei->Get_Physical_Face_Normal();
            Real* pfc = NULL;
            switch( physical_face ){
            case 0:
              pfc = PHYSICAL_FACE_NORMAL_1.Get_Scratch(node->ProcArrayIndex());
              break;
            case 1:
              pfc = PHYSICAL_FACE_NORMAL_2.Get_Scratch(node->ProcArrayIndex());
              break;
            case 2:
              pfc = PHYSICAL_FACE_NORMAL_3.Get_Scratch(node->ProcArrayIndex());
              break;
            default:
              pfc = node->Variable(NODE_NORMAL);
              break;
            }
            Real dot = pfs[0]*pfc[0] + pfs[1]*pfc[1] + pfs[2]*pfc[2];
#ifdef CONTACT_DEBUG_NODE
            if( PRINT_THIS_NODE ){
              postream << "    Tied physical face normal   : "
                     << pfs[0] << ", "
                     << pfs[1] << ", "
                     << pfs[2] << "\n"
                     << "    Current physical face normal: "
                     << pfc[0] << ", "
                     << pfc[1] << ", "
                     << pfc[2] << "\n"
                     << "    Dot product of pfs and pfc = " << dot << "\n"
                     << "    Physical face: " << physical_face << "\n";
            }
#endif
            if( dot > 0.98 ){
	      break;
	    }
          }
        }
        if(j<node->Number_NodeEntity_Interactions(1)) continue;
    } else {
      // For single constraints, keep the stored tied over anything new
      if(node->Number_NodeEntity_Interactions(1)) {
        ContactNodeEntityInteraction* cnei = node->Get_NodeEntity_Interactions(1)[0];
        if(cnei){
          if( cnei->Is_Tied() || cnei->Is_InfSlip() || cnei->Is_Glued()){
	    continue;
          }
	}
      }
    }
#ifdef CONTACT_DEBUG_NODE
    if( PRINT_THIS_NODE ){
      postream<<"  Node was not culled\n";
    }
#endif


    node_list[pairs] = node_list[i];
    node_keys[pairs] = node_key;
    physical_faces[pairs] = physical_face;
    ++pairs;
  }
  return pairs;
}

int ContactSearch::Cull_Node_List( int list_size, int* list, 
				   ContactNode<Real>* mnode, Real* distance)
{
  // This function culls down a node list for a given pair based on the 
  // following criteria
  //   1) Is the node different than the slave node
  //   2) Is the node owned by this processor (moved up to proximity search)
  //   3) Does the user want this node interacting with the node (moved up to
  //      proximity search)
  //   4) Check sphere/sphere intersection
  //
  int i,k;

  // Get the variable handles I need
  VariableHandle NODE_RADIUS =
    search_topology->Variable_Handle( ContactTopology::Node_Radius );
  VariableHandle CURRENT_POSITION =
    search_topology->Variable_Handle( ContactTopology::Current_Position );

  int pairs = 0;
  for( i=0 ; i<list_size ; ++i ){
    ContactNode<Real>* snode = static_cast<ContactNode<Real>*>
                         (search_topology->NodeList()->Find(list[i]));
    bool valid = true;
    
#ifdef CONTACT_DEBUG_NODE
    bool PRINT_THIS_NODE = primary_topology->Is_a_Debug_Node( snode );
    postream<<"  processing sub-node "<<snode->Global_ID()<<"\n";
#endif
    
    int mnode_key = mnode->Entity_Key();
    int snode_key = snode->Entity_Key();
    
    // Test 3b: Does the user want this slave node interacting with this master node
    if (mnode_key == -1 || snode_key == -1) {
      valid = false;
    } else {
      int Interaction_Type = (int)
        search_data->Get_Search_Data( ContactSearch::INTERACTION_TYPE, 
     				      snode_key, mnode_key );
      if( Interaction_Type == NO_INTERACTION ) {
#ifdef CONTACT_DEBUG_NODE
        if( PRINT_THIS_NODE ) postream << "    Culled: test 3b, "<<snode_key<<" vs "<<mnode_key<<"\n";
#endif
        valid = false;
      }
    }
    if (!valid) continue;
    
    // Test 4: check sphere/sphere intersection
    Real current_dist  = 0.0;
    Real radius_0      = *mnode->Variable(NODE_RADIUS);
    Real radius_1      = *snode->Variable(NODE_RADIUS);
    Real critical_dist = (radius_0+radius_1)*(radius_0+radius_1);
    
    Real* Position0 = mnode->Variable(CURRENT_POSITION);
    Real* Position1 = snode->Variable(CURRENT_POSITION);
    for( k=0 ; k<dimensionality ; ++k ){
      current_dist += (Position1[k] - Position0[k]) * (Position1[k] - Position0[k]);
    }
    if (current_dist >= critical_dist) {
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE ) postream << "    Culled: test 4, dist = "<<std::sqrt(current_dist)<<"\n";
#endif
      valid = false;
    }
    if (!valid) continue;
#ifdef CONTACT_DEBUG_NODE
    if( PRINT_THIS_NODE ) postream << " Keep: dist = "<<std::sqrt(current_dist)<<"\n";
#endif
    distance[pairs] = std::sqrt(current_dist);

    list[pairs] = list[i];
    ++pairs;
  }
  return pairs;
}

int ContactSearch::Cull_Node_List_New1( ACME::ContactNode_Vector &node_list, 
				        ContactFace<Real>* face,
				        VariableHandle NODE_NORMAL,
				        VariableHandle FACE_NORMAL,
				        int number_of_configurations, 
                                        const ContactBoundingBox &current_face_box,
                                        const ContactBoundingBox &predicted_face_box,
                                        const vector<bool> &valid_inter,
                                        VariableHandle NODE_COORD_START,
                                        VariableHandle NODE_COORD_END)
{
  // This function culls down a node list for a given pair based on the 
  // following criteria
  //   1) Is the node connected to the face
  //   2) Is the node a shell node and is its 
  //      "mating" node connected to this face  
  //   3) Does the user want this node interacting with the face
  //
  int i,j;

  int nodes_per_face = face->Nodes_Per_Face();
  PRECONDITION( nodes_per_face <= MAX_NODES_PER_FACE );
  // Get the variable handles I need

  ContactNode<Real>** face_nodes = face->Nodes();

  int pairs = 0;
  int list_size = node_list.size();
  for( i=0 ; i<list_size ; ++i ){
    ContactNode<Real> *node = node_list[i];

#ifdef CONTACT_DEBUG_NODE
    bool PRINT_THIS_NODE = primary_topology->Is_a_Debug_Node( node );
    if( PRINT_THIS_NODE ){
      VariableHandle CURRENT_POSITION =   
	search_topology->Variable_Handle( ContactTopology::Current_Position );
      VariableHandle PREDICTED_POSITION = 
	search_topology->Variable_Handle( ContactTopology::Predicted_Position);
      VariableHandle AUGMENTED_POSITION = 
	search_topology->Variable_Handle( ContactTopology::Augmented_Position);
      postream << "Debug Node ("
               << node->Exodus_ID() << ") is paired with face " 
	       << face->Global_ID() << "\n";
      Real* fnv = face->Variable(FACE_NORMAL);
      Real* nnv = node->Variable(NODE_NORMAL);
      postream << "      Face Normal = " << fnv[0] << " " 
	       << fnv[1] << " " << fnv[2] << "\n";
      for( int idbgn=0 ; idbgn<face->Nodes_Per_Face() ; idbgn++ ){
	ContactNode<Real>* nd = face->Node(idbgn);
	int id;
	if (nd->Is_a_Shell_Node()) {
	  id = static_cast<ContactShellNode*> (nd)->Shell_Node_Base_ID();
        } else {
	  id = nd->Exodus_ID();
        }
	Real* c = nd->Variable(CURRENT_POSITION);
	postream << "        Node (" 
                 << nd->Exodus_ID() << ") -- host id ("<< id <<")\n";
	postream << "            Current:   "
		 << c[0] << " " << c[1] << " " << c[2] << "\n";
	if( number_of_configurations >= 2) { 
	  Real * pc = nd->Variable( PREDICTED_POSITION );
	  postream << "            Predicted: "
		   << pc[0] << " " << pc[1] << " " << pc[2] << "\n";
	}
	if( number_of_configurations == 3) { 
	  Real * ac = nd->Variable( AUGMENTED_POSITION );
	  postream << "            Augmented: "
		   << ac[0] << " " << ac[1] << " " << ac[2] << "\n";
	}
      }
      postream << "      Slave Node Normal = " << nnv[0] << " " 
	       << nnv[1] << " " << nnv[2] << "\n";
      Real* c = node->Variable(CURRENT_POSITION);
      postream << "      Slave Node -- host id (" 
               << node->Exodus_ID() << ")\n";
      postream << "            Current:   "
               << c[0] << " " << c[1] << " " << c[2] << "\n";
      if( number_of_configurations >= 2) { 
        Real * pc = node->Variable( PREDICTED_POSITION );
        postream << "            Predicted: "
                 << pc[0] << " " << pc[1] << " " << pc[2] << "\n";
      }
      if( number_of_configurations == 3) { 
        Real * ac = node->Variable( AUGMENTED_POSITION );
        postream << "            Augmented: "
                 << ac[0] << " " << ac[1] << " " << ac[2] << "\n";
      }
    }
#endif

    //
    // Test 1: Cull out the nodes that make up the face
    //
    for( j=0 ; j<nodes_per_face ; ++j ){
      if( node == face_nodes[j] ){
#ifdef CONTACT_DEBUG_NODE
	if( PRINT_THIS_NODE ) postream << "    Culled: Part of Face\n";
#endif
        break;
      }
    }
    if(j<nodes_per_face) continue;

    //
    // Test 2: Is the node a shell node and is its 
    //         "mating" node connected to this face
    //
    if( node->Is_a_Shell_Node() && ( face->FaceType() == SHELLQUADFACEL4 || face->FaceType() == SHELLTRIFACEL3 ) ){
      //
      // Test 2b: Is this a shell node whose "mating" node is connected to
      //          this face?
      int shell_node_base_id =  (static_cast<ContactShellNode*>(node))->Shell_Node_Base_ID();
      for( j=0 ; j<nodes_per_face ; ++j ){
        int face_shell_node_base_id = (static_cast<ContactShellNode*>(face_nodes[j]))->Shell_Node_Base_ID();
	if (face_shell_node_base_id == shell_node_base_id) {
#ifdef CONTACT_DEBUG_NODE
	  if( PRINT_THIS_NODE ) postream << "    Culled: Opposing Shell Face\n";
#endif
	  break;
	}
      }
      if(j<nodes_per_face) continue;      
    }
    
    //
    // Test 3: Does the user want this node interacting with the face
    //
    /*
    bool found_valid = false;
    int num_physical_faces = (int) *(NUMBER_PHYSICAL_FACES.Get_Scratch(node->ProcArrayIndex()));
    for (int kk=0; kk<num_physical_faces; ++kk) {
      int start_face = (int) (PHYSICAL_FACE_LIST_INDEX.Get_Scratch(node->ProcArrayIndex()))[kk];
      int end_face;
      if( kk < num_physical_faces - 1){
        end_face = (int)(PHYSICAL_FACE_LIST_INDEX.Get_Scratch(node->ProcArrayIndex()))[kk+1];
      } else {
        end_face = node->Number_Face_Connections();
      }
      for( ii=start_face ; ii<end_face ; ++ii ){
        int node_key = node->Face(ii)->Entity_Key();
        if( valid_inter[node_key + 1] ){
	  found_valid = true;
          break;
        }
      }
      if( found_valid ) break;
    }
    if( found_valid ){
      node_list[pairs] = node_list[i];
      ++pairs;
    }
    */
    node_list[pairs] = node_list[i];
    ++pairs;
  }
  return pairs;
}

void ContactSearch::Compute_Nodes_in_Box(int& ne,
                                         Real* const xmin,
                                         Real* const xmax,
                                         Real* position,
                                         int* index,
                                         int& number_of_nodes,
                                         Real* scratch,
                                         int* rank2,
                                         int& list_size,
                                         int* list) {
  int ilo[MAX_DIMENSIONALITY], iup[MAX_DIMENSIONALITY];
  for(int j = 0; j < dimensionality; ++j) {
    int offset = number_of_nodes*j;
    FORTRAN(contact_get_bound)( ne,&position[offset],&index[offset],
                                number_of_nodes,&xmin[j],&xmax[j],
                                number_of_nodes,&ilo[j],&iup[j],
                                number_of_nodes,scratch );
  }
  list_size = 0;
  FORTRAN(contact_make_list)( number_of_nodes,index,rank2,iup,ilo,ne,
                              list,list_size,ne,dimensionality );
}

void
ContactSearch::Set_Initial_Gap()
{
#if CONTACT_DEBUG_PRINT_LEVEL>=10
  postream << " Set_Initial_Gap() \n";   
#endif
  ContactTopology *topology = primary_topology;

  int number_of_nodes = topology->Number_of_Nodes();
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());
  for (int i=0; i<number_of_nodes; ++i) {
    ContactNode<Real>* node = Nodes[i];
    if(node->Ownership() == ContactTopologyEntity<Real>::OWNED){
      ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions();
      for (int j=0; j<node->Number_NodeEntity_Interactions(); ++j) {
        ContactNodeEntityInteraction* cnei = interactions[j];
        cnei->Set_Initial_Gap();
#if CONTACT_DEBUG_PRINT_LEVEL>=10
        Real gap0 = cnei->Get_Gap_Initial();
        postream << "  node: (" << i << ", " << j << ") initial gap : " << gap0 << "\n";
#endif
      }//end loop over all interactions on node
    }//end status flag if
  }//end loop over all nodes
}

void 
ContactSearch::Retrieve_Glued_Interactions(VariableHandle POSITION, 
                                           ContactSearch::Topology use_topology,
                                           int no_physical_faces)
{
  // If there are no glued interactions requested by the user, then we could 
  // simply return but this is hard to tell with enforcement models so just 
  // check.  What really needs to be done is to add a check of whether or not 
  // the user has specified glued interactions in the search data or if the 
  // enforcement has any models that change state to/from a glued interaction.

  ContactTopology *topology(0);
  switch(use_topology){
    case(ContactSearch::PRIMARY):
      topology = primary_topology;
      break;
    case(ContactSearch::SECONDARY):
      topology = secondary_topology;
      break;
    case(ContactSearch::SEARCH):
      topology = search_topology;
      break;
  }

  VariableHandle FACE_NORMAL = 
    topology->Variable_Handle( ContactTopology::Face_Normal );

  // There is the possibility some of the interactions are glued so loop through
  // all of them and "copy" ones from the old state to the new state
  #if CONTACT_DEBUG_PRINT_LEVEL>=2 || defined(CONTACT_HEARTBEAT)
  int cnt = 0;
  #endif
  int number_of_nodes = topology->Number_of_Nodes();
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());
  if (no_physical_faces) {
    for (int i=0; i<number_of_nodes; ++i) {
      ContactNode<Real>* node = Nodes[i];
      if(node->Ownership() == ContactTopologyEntity<Real>::OWNED){
        int num_glued_at_node = 0;
        ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions(1);
        for (int j=0; j<node->Number_NodeEntity_Interactions(1); ++j) {
          ContactNodeEntityInteraction* cnei = interactions[j];
          if( cnei->Is_Glued() ){
            ContactNodeEntityInteraction *nei = cnei->New_Instance(allocators);
            Real  gap0    = cnei->Get_Gap_Initial();
            Real* pf_norm = cnei->Get_Physical_Face_Normal();
            nei->Set_Source(ContactNodeEntityInteraction::RETRIEVED_GLUED);
            nei->Update_Tied_Interaction( POSITION, FACE_NORMAL,
                                          gap0, pf_norm );
            node->Store_NodeEntity_Interaction( num_glued_at_node++, nei );
            #if CONTACT_DEBUG_PRINT_LEVEL>=2 || defined(CONTACT_HEARTBEAT)
            ++cnt;
            #endif
          }//end is tied if
        }//end loop over all interactions on node
      }//end ownership if
    }//end loop over all nodes
  } else {
    for (int i=0; i<number_of_nodes; ++i) {
      ContactNode<Real>* node = Nodes[i];
      if(node->Ownership() == ContactTopologyEntity<Real>::OWNED){
        int num_glued_at_node = 0;
        ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions(1);
        for (int j=0; j<node->Number_NodeEntity_Interactions(1); ++j) {
          ContactNodeEntityInteraction* cnei = interactions[j];
          if( cnei->Is_Glued() ){
            ContactNodeEntityInteraction *nei = cnei->New_Instance(allocators);
            Real gap0 = nei->Get_Gap_Initial();
            nei->Set_Source(ContactNodeEntityInteraction::RETRIEVED_GLUED);
            nei->Update_Tied_Interaction( POSITION, FACE_NORMAL,
                                          gap0,
                                          NUMBER_PHYSICAL_FACES,
                                          PHYSICAL_FACE_NORMAL_1,
                                          PHYSICAL_FACE_NORMAL_2,
                                          PHYSICAL_FACE_NORMAL_3 );
            node->Store_NodeEntity_Interaction( num_glued_at_node++, nei );
            #if CONTACT_DEBUG_PRINT_LEVEL>=2 || defined(CONTACT_HEARTBEAT)
            ++cnt;
            #endif
          }//end is tied if
        }//end loop over all interactions on node
      }//end ownership if
    }//end loop over all nodes
  }
#if CONTACT_DEBUG_PRINT_LEVEL>=2
  postream<<"    Retrieved "<<cnt<<" glued interactions\n";
  #else
  #ifdef CONTACT_HEARTBEAT
  int knt = contact_global_sum(cnt, SearchComm );
  if (contact_processor_number(SearchComm)==0) {
    std::cout<<"    Retrieved "<<knt<<" glued interactions\n";
  }
  #endif
#endif
}

#if 0
//disabling this as it is no longer used, but leaving it here incase we need to resurect it later
//bmp version 2.7  3/1/2007
void 
ContactSearch::Retrieve_Tied_Interactions(VariableHandle POSITION, 
                                          ContactSearch::Topology use_topology, 
                                          ContactTopologyEntity<Real>::SearchContext status)
{
  // If there are no tied interactions requested by the user, then we could 
  // simply return but this is hard to tell with enforcement models so just check

  ContactTopology *topology(0);
  switch(use_topology){
    case(ContactSearch::PRIMARY):
      topology = primary_topology;
      break;
    case(ContactSearch::SECONDARY):
      topology = secondary_topology;
      break;
    case(ContactSearch::SEARCH):
      topology = search_topology;
      break;
  }

  VariableHandle FACE_NORMAL = 
    topology->Variable_Handle( ContactTopology::Face_Normal );

  // There is the possibility some of the interactions are tied so loop through
  // all of them and "copy" ones from the old state to the new state
#if CONTACT_DEBUG_PRINT_LEVEL>=2
  int cnt = 0;
#endif
  int number_of_nodes = topology->Number_of_Nodes();
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());
  for (int i=0; i<number_of_nodes; ++i) {
    ContactNode<Real>* node = Nodes[i];
    //if(node->CheckContext(status) &&
    //   node->Ownership() == ContactTopologyEntity<Real>::OWNED){
    if(node->Ownership() == ContactTopologyEntity<Real>::OWNED){
      int num_tied_at_node = 0;
      ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions(1);
      for (int j=0; j<node->Number_NodeEntity_Interactions(1); ++j) {
        ContactNodeEntityInteraction* cnei = interactions[j];
        if( cnei->Is_Tied() || cnei->Is_InfSlip() ){
          Real gap0 = cnei->Get_Gap_Initial();
#if CONTACT_DEBUG_PRINT_LEVEL>=10
          postream << "  node: (" << i << ", " << j << ") updating initial gap : " << gap0 << "\n";
#endif
          ContactNodeEntityInteraction *nei = cnei->New_Instance(allocators);
          nei->Set_Source(ContactNodeEntityInteraction::RETRIEVED_TIED);
          nei->Update_Tied_Interaction( POSITION, FACE_NORMAL, gap0,
                                        NUMBER_PHYSICAL_FACES,
                                        PHYSICAL_FACE_NORMAL_1,
                                        PHYSICAL_FACE_NORMAL_2,
                                        PHYSICAL_FACE_NORMAL_3 );
          node->Store_NodeEntity_Interaction( num_tied_at_node++, nei );
#if CONTACT_DEBUG_PRINT_LEVEL>=2
          ++cnt;
#endif
        }//end is tied if
      }//end loop over all interactions on node
    }//end status flag if
  }//end loop over all nodes
#if CONTACT_DEBUG_PRINT_LEVEL>=2
  postream<<"    Retrieved "<<cnt<<" tied interactions\n";
#endif
}

#endif

void 
ContactSearch::Retrieve_Tied_Interactions_From_Primary(VariableHandle POSITION)
{
  // If there are no tied interactions requested by the user, then simply return

  if (!(search_data->Have_Tied_Interactions())) return;

  //ghost all off processor faces and tie them to the interaction.
  ContactTopology *topology(primary_topology);

  VariableHandle FACE_NORMAL = 
    topology->Variable_Handle( ContactTopology::Face_Normal );

  // There is the possibility some of the interactions are tied so loop through
  // all of them and "copy" ones from the old state to the new state  
  #if CONTACT_DEBUG_PRINT_LEVEL>=2 || defined(CONTACT_HEARTBEAT)
  int cnt = 0;
  #endif
  int my_proc         = contact_processor_number( SearchComm );
  int number_of_nodes = topology->Number_of_Nodes();
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());
  for (int i=0; i<number_of_nodes; ++i) {
    ContactNode<Real>* node = Nodes[i];
    if(node->Ownership() == ContactTopologyEntity<Real>::OWNED){
      int num_tied_at_node = 0;
      ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions(STATE1);
      for (int j=0; j<node->Number_NodeEntity_Interactions(STATE1); ++j) {
        ContactNodeEntityInteraction* cnei = interactions[j];
        if( cnei->Is_Tied() || cnei->Is_InfSlip() ){
          #if CONTACT_DEBUG_PRINT_LEVEL>=2 || defined(CONTACT_HEARTBEAT)
          ++cnt;
          #endif
          node->SetContextBit(ContactTopologyEntity<Real>::TIED);
          node->ClearContextBit(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE);
          node->ClearContextBit(ContactTopologyEntity<Real>::TRACK_SEARCH_SLAVE);
          if (cnei->Get_Type()==ContactNodeEntityInteraction::NODE_FACE_INTERACTION) {
            ContactNodeFaceInteraction *cnfi = static_cast<ContactNodeFaceInteraction*>(cnei);
            if (cnfi->FaceEntityData()->owner != my_proc) continue;
            POSTCONDITION(cnfi->Face());
          }
          Real  gap0 = cnei->Get_Gap_Initial();
          Real* norm = cnei->Get_Physical_Face_Normal();
          ContactNodeEntityInteraction *nei = cnei->New_Instance(allocators);
          nei->Set_Source(ContactNodeEntityInteraction::RETRIEVED_TIED);
          nei->Update_Tied_Interaction( POSITION, FACE_NORMAL,gap0, norm );
          node->Store_NodeEntity_Interaction( num_tied_at_node++, nei );
        }//end is tied if
      }//end loop over all interactions on node
    }//end ownership if
  }//end loop over all nodes
#ifndef CONTACT_NO_MPI
  topology->UpdateTiedFaces();
#endif
  #if CONTACT_DEBUG_PRINT_LEVEL>=2
  postream<<"    Retrieved "<<cnt<<" tied interactions\n";
  #else
  #ifdef CONTACT_HEARTBEAT
  int knt = contact_global_sum(cnt, SearchComm );
  if (contact_processor_number(SearchComm)==0) {
    std::cout<<"    Retrieved "<<knt<<" tied interactions\n";
  }
  #endif
  #endif
}

void 
ContactSearch::Retrieve_Tracked_Interactions(ContactSearch::SearchType search_type, Real gap_tol)
{
  #if CONTACT_DEBUG_PRINT_LEVEL>=2 || defined(CONTACT_HEARTBEAT)
  int cnt = 0;
  #endif
  ContactTopology *topology = primary_topology;
  
  VariableHandle CURRENT_POSITION   = topology->Variable_Handle( ContactTopology::Current_Position );
  VariableHandle PREDICTED_POSITION = topology->Variable_Handle( ContactTopology::Predicted_Position );
  VariableHandle POSITION           = -1;
  switch (search_type) {
  case STATIC1CONFIG:
  case STATIC1CONFIG_TIED:
    POSITION = CURRENT_POSITION;
    break;
  case STATIC2CONFIG:
    POSITION = PREDICTED_POSITION;
    break;
  case DYNAMIC2CONFIG:
    // MWG:  Never use augmented position??
    POSITION = PREDICTED_POSITION;
    break;
  default:
    POSTCONDITION(0);
    break;
  }
  POSTCONDITION(POSITION>=0);
  
  int number_of_nodes = topology->Number_of_Nodes();
  ContactNode<Real>** Nodes = 
      reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());
  ContactTopologyEntityList* face_list = topology->FaceList();
  // There is the possibility some of the interactions are tied and/or tracked so loop 
  // through all of them and copy all the tied and valid tracked interactions from the
  // old state to the new state
  int my_proc = contact_processor_number( SearchComm );
  for (int i=0; i<number_of_nodes; ++i) {
    ContactNode<Real>* node = Nodes[i];
    #ifdef CONTACT_DEBUG_NODE
    bool PRINT_THIS_NODE = primary_topology->Is_a_Debug_Node( node );
    #endif
    // Only process nodes that are owned by the local processor 
    // and are tagged to be processed by the tracked search.
    if(node->CheckContext(ContactTopologyEntity<Real>::TRACK_SEARCH_SLAVE) &&
       node->Ownership() == ContactTopologyEntity<Real>::OWNED){
      int interaction_num     = 0;
      int num_tracked_at_node = 0;
      int num_interactions    = node->Number_NodeEntity_Interactions(ContactSearch::STATE1);
      ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions(ContactSearch::STATE1);
      #ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE && num_interactions>0){
        postream << "    processing "<<num_interactions<<" interactions for node " << node->Global_ID() << "\n";
      }
      #endif
      for (int j=0; j<num_interactions; ++j) {
        int old_interaction_num = j;
        ContactNodeEntityInteraction* cnei = interactions[j];
        if( cnei->Is_Tied() || cnei->Is_InfSlip() || cnei->Is_Glued()) continue;
        if (cnei->Get_Type()==ContactNodeEntityInteraction::NODE_FACE_INTERACTION) {
          ContactNodeFaceInteraction *cnfi = static_cast<ContactNodeFaceInteraction*>(cnei);
          int status0 = 0;
          int status1 = 0;
          ContactFace<Real>* face = cnfi->Face();
          int nfaces = cnfi->NumSharedFaces();
          bool check_shared_edge_faces = false;
          #ifdef CONTACT_DEBUG_NODE
          if( PRINT_THIS_NODE ){
            postream << "      processing interaction " << j << " (sliding) with face "<<face->Global_ID()<<", nfaces = "<<nfaces<<"\n";
          }
          #endif
          status0 = CheckForNodeFaceInteraction(node, face, search_type,
                                                cnfi,
                                                old_interaction_num,
                                                interaction_num,
                                                gap_tol);
          #ifdef CONTACT_DEBUG_NODE
          if( PRINT_THIS_NODE ){
            postream << "      return status = "<<status0<<"\n";
          }
          #endif
          if (status0) {
            // interaction detected
            #ifdef CONTACT_DEBUG_NODE
            if( PRINT_THIS_NODE ){
              postream << "      interaction found, number = "<<interaction_num <<"\n";
            }
            #endif
            if (nfaces>0) {
              ContactNodeFaceInteraction* nfi = 
                 static_cast<ContactNodeFaceInteraction*>
                    (node->Get_NodeEntity_Interaction(interaction_num));
              POSTCONDITION(nfi);
              if (nfi->Location()==2) {
                check_shared_edge_faces = true;
                #ifdef CONTACT_DEBUG_NODE
                if( PRINT_THIS_NODE ){
                  postream << "      interaction was outside so checking "<<nfaces<<" neighboring faces\n";
                }
                #endif
              }
            }
          } else {
            // interaction not detected, find nearest neighbor to check based on old interaction
            if (nfaces>0) {
              check_shared_edge_faces = true;
              #ifdef CONTACT_DEBUG_NODE
              if( PRINT_THIS_NODE ){
                postream << "      no interaction was found so checking "
                         <<nfaces<<" neighboring faces\n";
              }
              #endif
            }
          }
          if (check_shared_edge_faces) {
            for (int k=0; k<nfaces; ++k) {
              ContactTopologyEntity<Real>::connection_data *face_info = cnfi->SharedFaceData(k);
              ContactHostGlobalID GID( face_info->host_gid[0], face_info->host_gid[1] );
              ContactFace<Real>* NeighborFace = NULL;
              if (face_info->owner == my_proc) {
                NeighborFace = static_cast<ContactFace<Real>*>(face_list->Find(GID));
                POSTCONDITION(NeighborFace);
              } else {
                int block    = face_info->owner_proc_array_index;
                NeighborFace = static_cast<ContactFace<Real>*>
                   (topology->Ghosted_Face_Block(block)->FaceList()->Find( GID ));
                POSTCONDITION(NeighborFace);
              }
              POSTCONDITION(NeighborFace);
              #ifdef CONTACT_DEBUG_NODE
              if( PRINT_THIS_NODE ){
                postream << "      tracking check on face "<<NeighborFace->Global_ID()<<"\n";
              }
              #endif
              status1 = CheckForNodeFaceInteraction(node, NeighborFace, search_type,
                                                    cnfi,
                                                    old_interaction_num,
                                                    interaction_num,
                                                    gap_tol);
              #ifdef CONTACT_DEBUG_NODE
              if( PRINT_THIS_NODE ){
                postream << "      return status = "<<status1<<"\n";
              }
              #endif
            }
          }
          if (status0 || status1) {
            // interaction detected
            ++num_tracked_at_node;
            ++interaction_num;
            #if CONTACT_DEBUG_PRINT_LEVEL>=2 || defined(CONTACT_HEARTBEAT)
            ++cnt;
            #endif
          } else {
            #ifdef CONTACT_DEBUG_NODE
            if( PRINT_THIS_NODE ){
              postream<<"      >> dropped a tracked interaction\n";
            }
            #endif
          }
        }//end cnfi != NULL if
      }//end loop over all interactions on node
    }//end status flag if
  }//end loop over all nodes 
  #if CONTACT_DEBUG_PRINT_LEVEL>=2
  postream<<"    Retrieved "<<cnt<<" tracked interactions\n";
  #else
  #ifdef CONTACT_HEARTBEAT
  int knt = contact_global_sum(cnt, SearchComm );
  if (contact_processor_number(SearchComm)==0) {
    std::cout<<"    Retrieved "<<knt<<" tracked interactions\n";
  }
  #endif
  #endif
}

int
ContactSearch::CheckForNodeFaceInteraction(ContactNode<Real>* node, ContactFace<Real>* face,
                                           ContactSearch::SearchType search_type,
                                           ContactNodeFaceInteraction* old_cnfi, 
                                           int old_interaction_num,
                                           int new_interaction_num,
                                           Real gap_tol)
{
  // Only called when doing tracked search.  Note, since physical faces are not calculated
  // in the track search, the checks here are not as involved as the checks in the the
  // global search.
  // return status = 0, no interaction at location new_interaction_num
  //               = 1, interaction at location new_interaction_num
  int status = 0;
  ContactTopology *topology = primary_topology;
  
  Real user_search_tol  = search_data->Max_Search_Tolerance();
  Real user_tang_tol    = search_data->Max_Search_Tangential_Tolerance();
  Real REL_TANG_TOL_VAR = REL_TANG_TOL;
  VariableHandle FACE_NORMAL = 
    topology->Variable_Handle( ContactTopology::Face_Normal );
  VariableHandle NODE_NORMAL = 
    topology->Variable_Handle( ContactTopology::Node_Normal );
  
  VariableHandle CURRENT_POSITION   = topology->Variable_Handle( ContactTopology::Current_Position );
  VariableHandle PREDICTED_POSITION = topology->Variable_Handle( ContactTopology::Predicted_Position );
  VariableHandle AUGMENTED_POSITION = topology->Variable_Handle( ContactTopology::Augmented_Position );
  VariableHandle POSITION;
  switch (search_type) {
  case STATIC1CONFIG:
  case STATIC1CONFIG_TIED:
    POSITION = CURRENT_POSITION;
    break;
  case STATIC2CONFIG:
    POSITION = PREDICTED_POSITION;
    break;
  case DYNAMIC2CONFIG:
    POSITION = PREDICTED_POSITION;
    break;
  default:
    POSTCONDITION(0);
    break;
  }
  int node_key = old_cnfi->Get_Node_Key();
  int Interaction_Type = (int)
    search_data->Get_Search_Data( ContactSearch::INTERACTION_TYPE, 
                                  node_key, face->Entity_Key() );
#ifdef CONTACT_DEBUG_NODE
  bool PRINT_THIS_NODE = primary_topology->Is_a_Debug_Node( node );
  if( PRINT_THIS_NODE ){
    postream << "      Processing Potential Interaction for Node ("
             << node->Exodus_ID() << ") with Face "
             << face->Global_ID() << "\n";
    postream << "        node_key = "<<node_key<<"\n";
    postream << "        face_key = "<<face->Entity_Key()<<"\n";
    postream << "        type     = "<<Interaction_Type<<"\n";
  }
#endif
  if( Interaction_Type == NO_INTERACTION ) return status;
  
  if( node->NodeType() != ContactSearch::POINT ){
    Real* node_normal = node->Variable(NODE_NORMAL);
    Real* face_normal = face->Variable(FACE_NORMAL);
    Real dot_product = 0;
    for( int k=0 ; k<dimensionality ; ++k )
      dot_product += face_normal[k]*node_normal[k];
    if( (node->Physical_Type()!=ContactNode<Real>::SHELL_NODE && dot_product>=0.0 ) ||
        (node->Physical_Type()==ContactNode<Real>::SHELL_NODE && dot_product> 0.0 ) ){
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE ){
        if (node->Physical_Type()==ContactNode<Real>::SHELL_NODE && dot_product>0.0) {
          postream << "        node_ptype == SHELL_NODE and dot_product = "<<dot_product<<"\n";
        }
        if (node->Physical_Type()!=ContactNode<Real>::SHELL_NODE && dot_product>=0.0) {
          postream << "        node_ptype != SHELL_NODE and dot_product = "<<dot_product<<"\n";
        }
      }
#endif
      return status;
    }
  }

  int nn;
  int nfacets = 0;
  ContactNodeFaceInteraction::InteractionSource Process_Method =
    ContactNodeFaceInteraction::UNKNOWN_SOURCE;
  switch (search_type) {
  case STATIC1CONFIG:
  case STATIC1CONFIG_TIED:
  case STATIC2CONFIG: {
    if (search_type==STATIC1CONFIG || search_type==STATIC1CONFIG_TIED) {
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE ){
        postream << "        STATIC1CONFIG search\n";
      }
#endif
      Process_Method = ContactNodeFaceInteraction::CLOSEST_POINT_PROJECTION_1;
    } else {
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE ){
        postream << "         STATIC2CONFIG search\n";
      }
#endif
      Process_Method = ContactNodeFaceInteraction::CLOSEST_POINT_PROJECTION_2;
    }
    face->FacetDecomposition(nfacets,ms_coordinates_c,ms_normals_c,POSITION);
    POSTCONDITION(nfacets<=4);
    int npairs = 1;
    Real* coords_node = node->Variable(POSITION);
    for (nn=0; nn<nfacets; ++nn) {
      pushback_dir_flag[nn] = 2; // normal pushback
    }
    for (nn=0; nn<nfacets*data_size; ++nn) {
      ctrcl_facets[nn] = 0.0;
    }
    for (nn=0; nn<nfacets; ++nn) {
      if (dimensionality==3) {
        FORTRAN(cnodetriangle_cpproj)(npairs, coords_node,
                                      &ms_coordinates_c[nn*9],
                                      &ms_normals_c[nn*3],
                                      &pushback_dir_flag[nn],
                                      ctrl,
                                      &ctrcl_facets[nn*data_size],
                                      REL_TANG_TOL_VAR,
                                      user_tang_tol,
                                      gap_tol);
      }
    }
    face->FacetStaticRestriction(nfacets, ms_coordinates_c, ms_normals_c,
                                 ctrcl_facets, ctrcl);
    }break;
  case DYNAMIC2CONFIG:{
#ifdef CONTACT_DEBUG_NODE
    if( PRINT_THIS_NODE ){
      postream << "        DYNAMIC2CONFIG search\n";
    }
#endif
    int npairs = 1;
    Process_Method = ContactNodeFaceInteraction::CLOSEST_POINT_PROJECTION_2;
#ifdef CONTACT_DEBUG_NODE
    if( PRINT_THIS_NODE ){
      switch (Process_Method) {
      case ContactNodeEntityInteraction::UNKNOWN_SOURCE:
        postream << "        Process_Method = UNKNOWN_SOURCE\n";
        break;
      case ContactNodeEntityInteraction::CLOSEST_POINT_PROJECTION_1:
        postream << "        Process_Method = CLOSEST_POINT_PROJECTION_1\n";
        break;
      case ContactNodeEntityInteraction::CLOSEST_POINT_PROJECTION_2:
        postream << "        Process_Method = CLOSEST_POINT_PROJECTION_2\n";
        break;
      case ContactNodeEntityInteraction::MOVING_INTERSECTION:
        postream << "        Process_Method = MOVING_INTERSECTION\n";
        break;
      case ContactNodeEntityInteraction::CLOSEST_POINT_PROJECTION_2_FROM_MOVING:
        postream << "        Process_Method = CLOSEST_POINT_PROJECTION_2_FROM_MOVING\n";
        break;
      case ContactNodeEntityInteraction::RETRIEVED_TIED:
        postream << "        Process_Method = RETRIEVED_TIED\n";
        break;
      case ContactNodeEntityInteraction::RETRIEVED_GLUED:
        postream << "        Process_Method = RETRIEVED_GLUED\n";
        break;
      }
    }
#endif
    if(Process_Method==ContactNodeFaceInteraction::MOVING_INTERSECTION) {
      face->FacetDecomposition(nfacets,
                         ms_coordinates_c, ms_normals_c, CURRENT_POSITION,
                         ms_coordinates_a, ms_normals_a, AUGMENTED_POSITION,
                         ms_coordinates_p, ms_normals_p, PREDICTED_POSITION );
    } else {
      face->FacetDecomposition(nfacets,
                         ms_coordinates_a, ms_normals_a, AUGMENTED_POSITION,
                         ms_coordinates_p, ms_normals_p, PREDICTED_POSITION);
    }
#ifdef CONTACT_DEBUG_NODE
    if( PRINT_THIS_NODE ){
      postream << "        processing "<<nfacets<<" facets\n";
    }
#endif
    Real* sn_coordinates_c = node->Variable(CURRENT_POSITION);
    Real* sn_coordinates_p = node->Variable(PREDICTED_POSITION);
    Real* sn_coordinates_a = node->Variable(AUGMENTED_POSITION);
    for (nn=0; nn<nfacets; ++nn) {
      pushback_dir_flag[nn] = 2; // normal pushback
    }
    for (nn=0; nn<nfacets*data_size; ++nn) {
      ctrcl_facets[nn] = 0.0;
    }
    if(Process_Method==ContactNodeFaceInteraction::MOVING_INTERSECTION){
      for (nn=0; nn<nfacets; ++nn) {
        if (dimensionality==3) {
          if (auto_tol) {
          FORTRAN(cnodetriangle_movsrch_aug)( npairs,
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
                                        REL_TANG_TOL_VAR,
                                        user_tang_tol );
          } else {
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
                                              user_search_tol,
                                              user_tang_tol );
          }
        }
#ifdef CONTACT_DEBUG_NODE
        if( PRINT_THIS_NODE ){
          if( ctrcl_facets[nn*data_size] == 0 )
            postream << "        Rejected with subtriangle " << nn << "\n";
          else if( ctrcl_facets[nn*data_size] > 0 )
            postream << "        Accepted with subtriangle " << nn << "\n";
          else
            postream << "        Deferred to CPP for subtriangle " << nn << "\n";
        }
#endif
      }
      face->FacetDynamicRestriction(nfacets, ctrcl_facets, ctrcl);
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE && ctrcl[0] > 0 ){
        Real local_coords[3];
        face->Compute_Local_Coordinates( POSITION, 
                                         &ctrcl[1], 
                                         local_coords );
        postream << "             Global Coords: (" 
                 << ctrcl[1] << "," << ctrcl[2] << "," 
                 << ctrcl[3] << ")\n";
        postream << "             Local Coords: (" 
                 << local_coords[0] << "," << local_coords[1] << "," 
                 << local_coords[2] << ")\n";
        postream << "             Gap:           " 
                 << ctrcl[4] << "\n";
        postream << "             Pushback Dir:  (" 
                 << ctrcl[5] << "," << ctrcl[6] << "," 
                 << ctrcl[7] << ")\n";
        postream << "             Normal Dir:    (" 
                 << ctrcl[8] << "," << ctrcl[9] << "," 
                 << ctrcl[10] << ")\n";
        postream << "             Location: "
                 << (int) ctrcl[11] << "\n";
        postream << "             Time - In/Out: "
                 << ctrcl[12] << "\n";
      }
#endif
    }
    // process interaction(s) flagged for CPP
    bool CPProj_Search = false;
    if( Process_Method ==
        ContactNodeFaceInteraction::CLOSEST_POINT_PROJECTION_2 ) 
      CPProj_Search = true;
    if( Process_Method == ContactNodeFaceInteraction::MOVING_INTERSECTION ){
      if( ctrcl[0] == -1 ){
        CPProj_Search = true;
        Process_Method = ContactNodeFaceInteraction::
          CLOSEST_POINT_PROJECTION_2_FROM_MOVING;
      }
    }
    if ( CPProj_Search ) {
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE ){
          postream << "        REL_TANG_TOL     = "<<REL_TANG_TOL_VAR<<"\n";
          postream << "        user_tang_tol    = "<<user_tang_tol<<"\n";
          postream << "        gap_tol          = "<<gap_tol<<"\n";
          postream << "        sn_coordinates_p = "<<sn_coordinates_p[0]<<", "<<sn_coordinates_p[1]<<", "<<sn_coordinates_p[2]<<"\n";
          postream << "        sn_coordinates_a = "<<sn_coordinates_a[0]<<", "<<sn_coordinates_a[1]<<", "<<sn_coordinates_a[2]<<"\n";
      }
#endif
      for (nn=0; nn<nfacets; ++nn) {
#ifdef CONTACT_DEBUG_NODE
        if( PRINT_THIS_NODE ){
            postream << "        Processing subtriangle " << nn << "\n";
            postream << "          ms_normals_p      = "<<ms_normals_p[0]<<", "<<ms_normals_p[1]<<", "<<ms_normals_p[2]<<"\n";
            postream << "          ms_normals_a      = "<<ms_normals_a[0]<<", "<<ms_normals_a[1]<<", "<<ms_normals_a[2]<<"\n";
            postream << "          ms_coordinates_p  = "<<ms_coordinates_p[0]<<", "<<ms_coordinates_p[1]<<", "<<ms_coordinates_p[2]<<"\n";
            postream << "                              "<<ms_coordinates_p[3]<<", "<<ms_coordinates_p[4]<<", "<<ms_coordinates_p[5]<<"\n";
            postream << "                              "<<ms_coordinates_p[6]<<", "<<ms_coordinates_p[7]<<", "<<ms_coordinates_p[8]<<"\n";
            postream << "          ms_coordinates_a  = "<<ms_coordinates_a[0]<<", "<<ms_coordinates_a[1]<<", "<<ms_coordinates_a[2]<<"\n";
            postream << "                              "<<ms_coordinates_a[3]<<", "<<ms_coordinates_a[4]<<", "<<ms_coordinates_a[5]<<"\n";
            postream << "                              "<<ms_coordinates_a[6]<<", "<<ms_coordinates_a[7]<<", "<<ms_coordinates_a[8]<<"\n";
            postream << "          pushback_dir_flag = "<<pushback_dir_flag[nn]<<"\n";
        }
#endif
        if (auto_tol) {
        FORTRAN(cnodetriangle_cpproj_aug)(npairs,
                                        sn_coordinates_p,
                                        &ms_coordinates_p[nn*9], 
                                        &ms_normals_p[nn*3],
                                        sn_coordinates_a,
                                        &ms_coordinates_a[nn*9], 
                                        &ms_normals_a[nn*3],
                                        &pushback_dir_flag[nn],
                                        ctrl,
                                        &ctrcl_facets[nn*data_size],
                                        REL_TANG_TOL_VAR,
                                        user_tang_tol,
                                        gap_tol);
        } else {
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
                                          REL_TANG_TOL_VAR,
                                          user_tang_tol,
                                          gap_tol);
        }
#ifdef CONTACT_DEBUG_NODE
        if( PRINT_THIS_NODE ){
          if( ctrcl_facets[nn*data_size] == 0 )
            postream << "          Rejected with subtriangle " << nn << "\n";
          else
            postream << "          Accepted with subtriangle " << nn << ", location = "<<(int)ctrcl_facets[11]<<"\n";
        }
#endif
      }
      face->FacetStaticRestriction(nfacets, ms_coordinates_p, 
                                   ms_normals_p, ctrcl_facets, ctrcl);
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE && ctrcl[0] > 0 ){
        Real local_coords[3];
        face->Compute_Local_Coordinates( POSITION, 
                                         &ctrcl[1], 
                                         local_coords );
        postream << "             Global Coords: (" 
                 << ctrcl[1] << "," << ctrcl[2] << "," 
                 << ctrcl[3] << ")\n";
        postream << "             Local Coords: (" 
                 << local_coords[0] << "," << local_coords[1] << "," 
                 << local_coords[2] << ")\n";
        postream << "             Gap:           " 
                 << ctrcl[4] << "\n";
        postream << "             Pushback Dir:  (" 
                 << ctrcl[5] << "," << ctrcl[6] << "," 
                 << ctrcl[7] << ")\n";
        postream << "             Normal Dir:    (" 
                 << ctrcl[8] << "," << ctrcl[9] << "," 
                 << ctrcl[10] << ")\n";
        postream << "             Location: "
                 << (int) ctrcl[11] << "\n";
        postream << "             Time - In/Out: "
                 << ctrcl[12] << "\n";
      }
#endif
    }

    }break;
  default:
    POSTCONDITION(0);
    break;
  }
  
  ContactNodeFaceInteraction* nfi_stored = static_cast<ContactNodeFaceInteraction*>
        (node->Get_NodeEntity_Interaction( new_interaction_num ));
        
  if( ctrcl[0] == 1 ){
#ifdef CONTACT_DEBUG_NODE
    if( PRINT_THIS_NODE ){
      postream << "        Valid Interaction for location "<<new_interaction_num<<"\n";
    }
#endif
    // The interaction is valid
    status = 1;
    
    Real* pf_norm = NULL;
#ifdef COMPUTE_FOR_SERIAL
    if( contact_number_of_processors(SearchComm)==1) {
      int physical_face = -1;
      int num_physical_faces = (int) *(NUMBER_PHYSICAL_FACES.Get_Scratch(node->ProcArrayIndex()));
      if( num_physical_faces>0 ){
        Real* pf_normals[3];
        Real most_opposed = 2.0;
        Real* face_normal = face->Variable(FACE_NORMAL);
        pf_normals[0] = PHYSICAL_FACE_NORMAL_1.Get_Scratch(node->ProcArrayIndex());
        pf_normals[1] = PHYSICAL_FACE_NORMAL_2.Get_Scratch(node->ProcArrayIndex());
        pf_normals[2] = PHYSICAL_FACE_NORMAL_3.Get_Scratch(node->ProcArrayIndex());
        Find_Physical_Face(num_physical_faces, pf_normals, face_normal, physical_face, most_opposed);
        pf_norm = pf_normals[physical_face];
      }
    } else {
      pf_norm = old_cnfi->Get_Physical_Face_Normal();
      //pf_norm = node->Variable(NODE_NORMAL);  // since no physical faces in the track search
    }
#else
    pf_norm = old_cnfi->Get_Physical_Face_Normal();
    //pf_norm = node->Variable(NODE_NORMAL);  // since no physical faces in the track search
#endif
    
    ContactNodeFaceInteraction* nfi =
      ContactNodeFaceInteraction::new_ContactNodeFaceInteraction(
          allocators[ALLOC_ContactNodeFaceInteraction],
          Process_Method, node, face, ctrcl, node_key,
          pf_norm, false, false, false, true, enable_off_face_tracking, POSITION );
    
    Real& slave_node_area = nfi->Scalar_Var(ContactNodeFaceInteraction::NODE_AREA);
    slave_node_area = node->Get_NodeEntity_Interaction(old_interaction_num,1)->Scalar_Var(ContactNodeFaceInteraction::NODE_AREA);
    POSTCONDITION(node==nfi->Node());
    if( !nfi_stored ){
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE ){
        postream << "        no interaction previously stored at this location, storing interaction\n";
      }
#endif
      node->Store_NodeEntity_Interaction( new_interaction_num, nfi);
    } else {
      if (Compete_Interactions_Connected(nfi_stored,nfi)==nfi_stored) {
#ifdef CONTACT_DEBUG_NODE
        if( PRINT_THIS_NODE ){
          postream << "        keeping the stored interaction\n";
        }
#endif
	nfi->Delete_Frag(allocators);
      } else {
#ifdef CONTACT_DEBUG_NODE
        if( PRINT_THIS_NODE ){
          postream << "        storing new interaction\n";
        }
#endif
	node->Store_NodeEntity_Interaction( new_interaction_num, nfi);
      }
    }
  } else {
    if( nfi_stored ) status = 1;
  }
    
#ifdef CONTACT_DEBUG_NODE
  if( PRINT_THIS_NODE ) {
    postream<<"        >>status = "<<status<<"\n";
    //postream.flush();
  }
#endif
  return status;
}

void ContactSearch::Compute_Remaining_Gap( Real* max_gap_mag, 
					   Real* max_gap,
                                           ContactSearch::Topology use_topology,
                                           ContactTopologyEntity<Real>::SearchContext status )
{
  ContactTopology *topology(0);
  switch(use_topology){
    case ContactSearch::PRIMARY:
      topology = primary_topology;
      break;
    case ContactSearch::SECONDARY :
      topology = secondary_topology;
      break;
    case ContactSearch::SEARCH :
      topology = search_topology;
      break;
    default:
      PRECONDITION(0);
      break;
  }

  max_gap[0] = 0.0;
  max_gap[1] = 0.0;
  max_gap[2] = 0.0;
  Real remaining_gap[3] = {0.0, 0.0, 0.0};

  Real shape_functions[MAX_NODES_PER_FACE];

  // Get the variable handles for everything I need
  VariableHandle CURRENT_POSITION = topology->Variable_Handle( ContactTopology::Current_Position );
  VariableHandle REMAINING_GAP    = topology->Variable_Handle( ContactTopology::Remaining_Gap );

  int number_of_nodes = topology->Number_of_Nodes();
  ContactNode<Real>** Nodes =  reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());

  for (int i=0; i<number_of_nodes; ++i) {
    ContactNode<Real>* node = Nodes[i];
    
    //if this node is not locally owned then skip it.
    if (node->Ownership() != ContactTopologyEntity<Real>::OWNED) continue;
    
    //if this node doesn't have the right status then skip it.
    //if(!node->CheckContext(status)) continue;
    
#ifdef CONTACT_DEBUG_NODE
    bool PRINT_THIS_NODE = primary_topology->Is_a_Debug_Node( node );
    if( node->Ownership() != ContactTopologyEntity<Real>::OWNED ) PRINT_THIS_NODE = false;
#endif

    Real* pgap = node->Variable( REMAINING_GAP );
    for( int j=0 ; j<dimensionality ; ++j ) pgap[j] = 0.0;
    
    // Loop over all the old interactions
    ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions(1);
//    postream << "Number of interactions = " << node->Number_NodeEntity_Interactions(1) << "\n";
    for (int j=0; j<node->Number_NodeEntity_Interactions(1); ++j) {
      ContactNodeEntityInteraction* cnei = interactions[j];
      // If the gap is negative, compute the gap left over
      if( cnei->Get_Gap_Cur() <= 0 ){
        // Get the global position of the contact point on the face
#ifdef CONTACT_DEBUG_NODE
        if( PRINT_THIS_NODE ){
          postream << "Computing Remaining Gap for Node ("
                   << cnei->Node()->Exodus_ID() << ") ";
          if(cnei->Get_Type() == ContactNodeEntityInteraction::NODE_FACE_INTERACTION) {
            ContactNodeFaceInteraction*  cnfi = 
              static_cast<ContactNodeFaceInteraction*>(cnei);
            ContactFace<Real>* face = cnfi->Face();
            POSTCONDITION(face);
            postream << "with Face " << face->Global_ID();
            break;
	  } else if (cnei->Get_Type() == ContactNodeEntityInteraction::NODE_SURFACE_INTERACTION) {
            ContactNodeSurfaceInteraction* cnsi = 
              static_cast<ContactNodeSurfaceInteraction*>(cnei);
            postream << "with Surface " << cnsi->SurfaceID();
            break;
          }
          postream << "\n";
        }
#endif
        Real* coordinates   = cnei->Get_Coordinates();
        Real* contact_point = cnei->Get_Contact_Point();
        Real  contact_point2[3];
        if(cnei->Get_Type()==ContactNodeEntityInteraction::NODE_FACE_INTERACTION){
          ContactNodeFaceInteraction* cnfi = static_cast<ContactNodeFaceInteraction*>(cnei);
          ContactFace<Real>* face = cnfi->Face();
          PRECONDITION( MAX_NODES_PER_FACE >= face->Nodes_Per_Face() );
          face->Evaluate_Shape_Functions( coordinates, shape_functions );
          contact_point2[0]=0.0; contact_point2[1]=0.0; contact_point2[2]=0.0;
          for( int k=0 ; k<face->Nodes_Per_Face() ; ++k ){
            ContactNode<Real>* Node = face->Node(k);
            Real* coord = Node->Variable(CURRENT_POSITION);
            contact_point2[0] += shape_functions[k]*coord[0];
            contact_point2[1] += shape_functions[k]*coord[1];
            contact_point2[2] += shape_functions[k]*coord[2];
          }
          //postream << "comp gap Coordinates   = " << coordinates[0] << ", "
          //                                        << coordinates[1] << ", "
          //                                        << coordinates[2] << "\n"
          //         << "Original contact_point = " << contact_point2[0]  << ", "
          //                                        << contact_point2[1]  << ", "
          //                                        << contact_point2[2]  << "  "
          //                                        << "For node " << cnei->Node()->Exodus_ID() << "\n"
          //         << "New contact point      = " << contact_point2[0] << ", "
          //                                        << contact_point2[1] << ", "
          //                                        << contact_point2[2] << "  "
          //                                        << "For node " << cnei->Node()->Exodus_ID() << "\n";
          contact_point[0] = contact_point2[0];
          contact_point[1] = contact_point2[1];
          contact_point[2] = contact_point2[2];
        }
#ifdef CONTACT_DEBUG_NODE
        if( PRINT_THIS_NODE )
          postream << "  Contact Point:  " << contact_point[0] << " "
                   << contact_point[1] << " " << contact_point[2] << "\n";
#endif
        // Compute a vector from the contact point to the node
        Real vec[3];
        Real* node_pos = node->Variable(CURRENT_POSITION);
        vec[0] = contact_point[0] - node_pos[0];
        vec[1] = contact_point[1] - node_pos[1];
        vec[2] = contact_point[2] - node_pos[2];
        // Retrieve the pushback direction
        Real* pbdir = cnei->Get_Pushback_Dir();
#ifdef CONTACT_DEBUG_NODE
        if( PRINT_THIS_NODE )
          postream << "  Pushback:  " << pbdir[0] << " "
                                      << pbdir[1] << " " 
                                      << pbdir[2] << "\n";
#endif
        // Compute the normal gap not enforced
        Real pmagl = pbdir[0]*vec[0] + pbdir[1]*vec[1] + pbdir[2]*vec[2]; 
#ifdef CONTACT_DEBUG_NODE
        if( PRINT_THIS_NODE )
          postream << "  Unenforced Normal Gap:  " << pmagl << "\n";
#endif  
        // Keep the maximum distance in each coordinate
        pgap[0] = std::max( pgap[0], std::fabs(pmagl*pbdir[0]) );
        pgap[1] = std::max( pgap[1], std::fabs(pmagl*pbdir[1]) );
        pgap[2] = std::max( pgap[2], std::fabs(pmagl*pbdir[2]) );
        // Keep the maximum over all nodes in each coordinate
        remaining_gap[0] = std::max( remaining_gap[0], pgap[0] );
        remaining_gap[1] = std::max( remaining_gap[1], pgap[1] );
        remaining_gap[2] = std::max( remaining_gap[2], pgap[2] );
#ifdef CONTACT_DEBUG_NODE
        if( primary_topology->Is_a_Debug_Node( node ) ){
          if( node->Owner() == contact_processor_number(SearchComm) ){
            postream << "Remaining Gap = " << pgap[0] << " " << pgap[1]
                     << " " << pgap[2] << "\n";
          }
        }
#endif
      }
    }
  }
  contact_global_maximum(remaining_gap, max_gap, 3, SearchComm);
  *max_gap_mag = Magnitude(max_gap);
}

void ContactSearch::Compute_Remaining_Gap( Real* max_gap_mag, 
					   Real* max_gap )
{
  ContactTopology *topology = primary_topology;

  max_gap[0] = 0.0;
  max_gap[1] = 0.0;
  max_gap[2] = 0.0;
  Real remaining_gap[3] = {0.0, 0.0, 0.0};

  VariableHandle REMAINING_GAP = topology->Variable_Handle( ContactTopology::Remaining_Gap );

  int number_of_nodes = topology->Number_of_Nodes();
  ContactNode<Real>** Nodes = reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());

  for (int i=0; i<number_of_nodes; ++i) {
    ContactNode<Real>* node = Nodes[i];
    
    //if this node is not locally owned then skip it.
    if (node->Ownership() != ContactTopologyEntity<Real>::OWNED) continue;
    
#ifdef CONTACT_DEBUG_NODE
    bool PRINT_THIS_NODE = primary_topology->Is_a_Debug_Node( node );
    if( node->Ownership() != ContactTopologyEntity<Real>::OWNED ) PRINT_THIS_NODE = false;
#endif

    Real* pgap = node->Variable( REMAINING_GAP );
    remaining_gap[0] = std::max( remaining_gap[0], pgap[0] );
    remaining_gap[1] = std::max( remaining_gap[1], pgap[1] );
    remaining_gap[2] = std::max( remaining_gap[2], pgap[2] );
  }
  contact_global_maximum(remaining_gap, max_gap, 3, SearchComm);
  *max_gap_mag = Magnitude(max_gap);
}

// void ContactSearch::Partition_Between_GapCur_and_GapOld( 
// 				           ContactNodeFaceInteraction* cnfi )
// {
//   // Upon entering this routine, all of the gap is in GAP_CUR and GAP_OLD
//   // has yet to be set.  The purpose of this routine is partition the gap
//   // into the part that occurred physically this time step (GAP_CUR) and
//   // the part of the gap that was pre-existing starting the step (GAP_OLD)
//   //
//   // There are six possible cases to handled for non-shell faces & nodes
//   // 1) Initially outside, moving towards the face but still outside
//   // 2) Initially outside, moving further outside
//   // 3) Initially outside, moving inside
//   // 4) Initially inside, moving further inside
//   // 5) Initially inside, moving outside
//   // 6) Initially inside, moving towards the face but still inside
//   //
//   
//   Real& gdotdt = cnfi->Scalar_Var(ContactNodeFaceInteraction::GAP_CUR);
//   Real& gapold = cnfi->Scalar_Var(ContactNodeFaceInteraction::GAP_OLD);
//   Real gtotal = gdotdt;
// 
//   // If the the total_gap is non-negative put everything in gap_cur
//   if( gtotal >= 0.0 ){
//     // Cases 1, 2 & 5
//     gapold = 0.0;
//     return;
//   }
// 
//   // Get the variable handles for everything I need
//   VariableHandle CURRENT_POSITION = 
//     search_topology->Variable_Handle( ContactTopology::Current_Position );
//   VariableHandle PREDICTED_POSITION =
//     search_topology->Variable_Handle( ContactTopology::Predicted_Position );
// 
//   // Get the pushback direction
//   Real* pb_dir = cnfi->Vector_Var(ContactNodeFaceInteraction::PUSHBACK_DIR);
// 
//   // Compute the amount of node motion along pushback_dir
//   Real* cur0 = cnfi->Node()->Variable(CURRENT_POSITION);
//   Real* curp = cnfi->Node()->Variable(PREDICTED_POSITION);
//   Real node_gdotdt = ((curp[0]-cur0[0])*pb_dir[0] +
// 		      (curp[1]-cur0[1])*pb_dir[1] +
// 		      (curp[2]-cur0[2])*pb_dir[2] );
//     
//   // Now compute the amount of face motion along pushback_dir
//   Real shape_functions[MAX_NODES_PER_FACE], contact_point_motion[3];
//   cnfi->Face()->Evaluate_Shape_Functions( 
// 	    cnfi->Vector_Var(ContactNodeFaceInteraction::COORDINATES),
// 	    shape_functions );
//   PRECONDITION( cnfi->Face()->Nodes_Per_Face() <= MAX_NODES_PER_FACE );
//   contact_point_motion[0]=contact_point_motion[1]=contact_point_motion[2] = 0.0;
//   for( int i=0 ; i<cnfi->Face()->Nodes_Per_Face() ; ++i ){
//     ContactNode<Real>* node = cnfi->Face()->Node(i);
//     Real* c0 = node->Variable(CURRENT_POSITION);
//     Real* cp = node->Variable(PREDICTED_POSITION);
//     contact_point_motion[0] += shape_functions[i]*(cp[0]-c0[0]);
//     contact_point_motion[1] += shape_functions[i]*(cp[1]-c0[1]);
//     contact_point_motion[2] += shape_functions[i]*(cp[2]-c0[2]);
//   }
//   Real face_gdotdt = (contact_point_motion[0]*pb_dir[0] +
// 		      contact_point_motion[1]*pb_dir[1] +
// 		      contact_point_motion[2]*pb_dir[2] );
//   
//   // Take the amount of the gap that was obtainable this step and store
//   // in gdotdt and leave the rest in gap_old;
//   Real physical_gap = node_gdotdt - face_gdotdt;
// 
//   if( physical_gap >= 0.0 ){
//     // case 6 (the other cases of moving "out" have been handled)
//     gapold = gtotal;
//     gdotdt = 0.0;
//   } else {
//     if( gtotal < physical_gap ){
//       // case 4
//       gdotdt = physical_gap;
//       gapold = gtotal - gdotdt;
//     } else {
//       // case 3
//       gapold = 0.0;
//     }
//   }
// }

void ContactSearch::Partition_Between_GapCur_and_GapOld( 
				           ContactNodeEntityInteraction* cnei )
{
  // Upon entering this routine, all of the gap is in GAP_CUR and GAP_OLD
  // has yet to be set.  The purpose of this routine is partition the gap
  // into the part that occurred physically this time step (GAP_CUR) and
  // the part of the gap that was pre-existing starting the step (GAP_OLD)
  //
  // There are six possible cases to handled for non-shell faces & nodes
  // 1) Initially outside, moving towards the face but still outside
  // 2) Initially outside, moving further outside
  // 3) Initially outside, moving inside
  // 4) Initially inside, moving further inside
  // 5) Initially inside, moving outside
  // 6) Initially inside, moving towards the face but still inside
  //
  
  Real& gdotdt = cnei->Get_Gap_Cur();
  Real& gapold = cnei->Get_Gap_Old();
  Real gtotal = gdotdt;

  // If the the total_gap is non-negative put everything in gap_cur
  if( gtotal >= 0.0 ){
    // Cases 1, 2 & 5
    gapold = 0.0;
    return;
  }

  // Get the variable handles for everything I need
  VariableHandle CURRENT_POSITION = 
    search_topology->Variable_Handle( ContactTopology::Current_Position );
  VariableHandle PREDICTED_POSITION =
    search_topology->Variable_Handle( ContactTopology::Predicted_Position );

  // Get the pushback direction
  Real* pb_dir = cnei->Get_Pushback_Dir();

  // Compute the amount of node motion along pushback_dir
  Real* cur0 = cnei->Node()->Variable(CURRENT_POSITION);
  Real* curp = cnei->Node()->Variable(PREDICTED_POSITION);
  Real node_gdotdt = ((curp[0]-cur0[0])*pb_dir[0] +
		      (curp[1]-cur0[1])*pb_dir[1] +
		      (curp[2]-cur0[2])*pb_dir[2] );

  Real physical_gap;
  if(cnei->Get_Type()==ContactNodeEntityInteraction::NODE_FACE_INTERACTION) {
    ContactNodeFaceInteraction *cnfi = static_cast<ContactNodeFaceInteraction*>(cnei);
    // Now compute the amount of face motion along pushback_dir
    Real shape_functions[MAX_NODES_PER_FACE], contact_point_motion[3];
    cnfi->Face()->Evaluate_Shape_Functions( 
	      cnfi->Vector_Var(ContactNodeFaceInteraction::COORDINATES),
	      shape_functions );
    PRECONDITION( cnfi->Face()->Nodes_Per_Face() <= MAX_NODES_PER_FACE );
    contact_point_motion[0]=contact_point_motion[1]=contact_point_motion[2] = 0.0;
    for( int i=0 ; i<cnfi->Face()->Nodes_Per_Face() ; ++i ){
      ContactNode<Real>* node = cnfi->Face()->Node(i);
      Real* c0 = node->Variable(CURRENT_POSITION);
      Real* cp = node->Variable(PREDICTED_POSITION);
      contact_point_motion[0] += shape_functions[i]*(cp[0]-c0[0]);
      contact_point_motion[1] += shape_functions[i]*(cp[1]-c0[1]);
      contact_point_motion[2] += shape_functions[i]*(cp[2]-c0[2]);
    }
    Real face_gdotdt = Dot(contact_point_motion,pb_dir);
  
    // Take the amount of the gap that was obtainable this step and store
    // in gdotdt and leave the rest in gap_old;
    physical_gap = node_gdotdt - face_gdotdt;
  } else {
    physical_gap = node_gdotdt;
  }

  if( physical_gap >= 0.0 ){
    // case 6 (the other cases of moving "out" have been handled)
    gapold = gtotal;
    gdotdt = 0.0;
  } else {
    if( gtotal < physical_gap ){
      // case 4
      gdotdt = physical_gap;
      gapold = gtotal - gdotdt;
    } else {
      // case 3
      gapold = 0.0;
    }
  }
}


// void ContactSearch::Partition_Between_GapCur_and_GapOld( 
// 				           ContactNodeSurfaceInteraction* cnsi )
// {
//   // There are six possible cases to handled
//   // 1) Initially outside, moving towards the face but still outside
//   // 2) Initially outside, moving further outside
//   // 3) Initially outside, moving inside
//   // 4) Initially inside, moving further inside
//   // 5) Initially inside, moving outside
//   // 6) Initially inside, moving towards the face but still inside
//   //
//   
//   // Get the variable handles for everything I need
//   VariableHandle CURRENT_POSITION = 
//     search_topology->Variable_Handle( ContactTopology::Current_Position );
//   VariableHandle PREDICTED_POSITION =
//     search_topology->Variable_Handle( ContactTopology::Predicted_Position );
// 
//   Real& gdotdt = cnsi->Scalar_Var(ContactNodeSurfaceInteraction::GAP_CUR);
//   Real& gapold = cnsi->Scalar_Var(ContactNodeSurfaceInteraction::GAP_OLD);
//   Real gtotal = gdotdt;
// 
//   // If the the total_gap is non-negative put everything in gap_cur
//   if( gtotal >= 0.0 ){
//     // Cases 1, 2 & 5
//     gapold = 0.0;
//     return;
//   }
// 
//   // Get the pushback direction
//   //bmp note: why is this the surface normal and not the pushback dir?
//   Real* n = cnsi->Get_Normal();
// 
//   // Compute the amount of node motion along pushback_dir
//   Real* cur0 = cnsi->Node()->Variable(CURRENT_POSITION);
//   Real* curp = cnsi->Node()->Variable(PREDICTED_POSITION);
//   Real node_gdotdt = ((curp[0]-cur0[0])*n[0] +
// 		      (curp[1]-cur0[1])*n[1] +
// 		      (curp[2]-cur0[2])*n[2] );
//   
//   // Take the amount of the gap that was obtainable this step and store
//   // in gdotdt and leave the rest in gap_old;
//   Real physical_gap = node_gdotdt;
// 
//   if( physical_gap >= 0.0 ){
//     // case 6 (the other cases of moving "out" have been handled)
//     gapold = gtotal;
//     gdotdt = 0.0;
//   } else {
//     if( gtotal < physical_gap ){
//       // case 4
//       gdotdt = physical_gap;
//       gapold = gtotal - gdotdt;
//     } else {
//       // case 3
//       gapold = 0.0;
//     }
//   }
// }

ContactNodeEntityInteraction*
ContactSearch::Compete_Interactions( VariableHandle NODE_COORDS,
				     ContactNodeEntityInteraction* cnei1,
				     ContactNodeEntityInteraction* cnei2 )
{
  //
  // Note: It is assumed that interaction 1 is the stored interaction and
  //       interaction 2 is the candidate that may replace interaction 1.
  //       This matters because in one branch we are averaging normal and
  //       we need to know who has the previously averaged values.
  //
  PRECONDITION( cnei1 );
  PRECONDITION( cnei2 );
  PRECONDITION( cnei1->Node() == cnei2->Node() );

  ContactNodeFaceInteraction *cnfi1 = dynamic_cast<ContactNodeFaceInteraction*>(cnei1);
  ContactNodeFaceInteraction *cnfi2 = dynamic_cast<ContactNodeFaceInteraction*>(cnei2);
  
#ifdef CONTACT_DEBUG_NODE
  ContactNode<Real>* dn = cnfi1->Node();
  bool PRINT_THIS_NODE = primary_topology->Is_a_Debug_Node( dn );
  if( PRINT_THIS_NODE ){
    postream << "    Compete_Interactions Competing faces for Node ("
             << dn->Exodus_ID() << ")\n";
  }
#endif

  bool connected_faces = false;
  if(cnfi1 != NULL && cnfi2 != NULL) {
    connected_faces = search_topology->Faces_Connected(cnfi1->Face(), cnfi2->Face());
  }

  ContactNodeEntityInteraction* cnei;

  if( connected_faces ) {
#ifdef CONTACT_DEBUG_NODE
    if( PRINT_THIS_NODE ) {
      postream << "      Stored and candidate face are connected \n";
    }
#endif
    cnei = Compete_Interactions_Connected( cnfi1, cnfi2 );
    return cnei;
  } else {
#ifdef CONTACT_DEBUG_NODE
    if( PRINT_THIS_NODE ) {
      postream << "      Stored and candidate face are unconnected \n";
    }
#endif
    cnei = Compete_Interactions_Unconnected( NODE_COORDS,cnei1, cnei2 );
    return cnei;
  }
}
//
//  Return the strength of contact for a single node-face interaction
//
Real Strength_Contact(ContactNodeFaceInteraction* cnfi,
                      VariableHandle &FACE_NORMAL) {
  Real *normal = cnfi->Face()->Variable(FACE_NORMAL);
  Real* sn_normal   = cnfi->Vector_Var(ContactNodeFaceInteraction::PHYSICAL_FACE_NORMAL);
  return -(normal[0]*sn_normal[0] +
	   normal[1]*sn_normal[1] +
	   normal[2]*sn_normal[2] ); // normal strength of contact
}
//
//  Return the concavity for a pair of node-face interactions
//
int Concavity(ContactNodeFaceInteraction* cnfi1,
	      ContactNodeFaceInteraction* cnfi2,
              VariableHandle &CENTROID,
              VariableHandle &FACE_NORMAL) {
  ContactFace<Real>* face_1 = cnfi1->Face(); // interaction face
  ContactFace<Real>* face_2 = cnfi2->Face(); // interaction face
  Real* centroid_1    = face_1->Variable(CENTROID);
  Real* centroid_2    = face_2->Variable(CENTROID);
  Real* normal_1      = face_1->Variable(FACE_NORMAL);
  Real* normal_2      = face_2->Variable(FACE_NORMAL);

  // vector from centroid 1 to centroid 2
  Real vec_1_2[3];
  vec_1_2[0] = centroid_2[0] - centroid_1[0];
  vec_1_2[1] = centroid_2[1] - centroid_1[1];
  vec_1_2[2] = centroid_2[2] - centroid_1[2];

  // projection of normal 1 onto the vector between centroids
  Real proj1 = Dot(normal_1, vec_1_2);

  // projection of normal 2 onto the vector between centroids
  Real proj2 = Dot(normal_2, vec_1_2);

  // determine concavity/convexity
  int iconcave;
  if( proj1 > proj2 )
    iconcave = 1;  // surfaces are concave
  else
    iconcave = -1; // surfaces are convex
  
  Real mag_v_1_2 = Dot(vec_1_2, vec_1_2);

  if( mag_v_1_2 > 0.0 ){
    if( std::fabs(proj2) < 1.0E-2*std::sqrt(mag_v_1_2) ) iconcave = 0;
  }
  return iconcave;
}

//
//  Determine if two faces are smoothly connected
//
bool Connected_Smooth(ContactNodeFaceInteraction* cnfi1,
	              ContactNodeFaceInteraction* cnfi2,
                      const Real sharp_smooth_curvature,
                      VariableHandle &FACE_NORMAL) {
  //
  //**********************************************************
  // Are contact faces connected smoothly?
  //**********************************************************
  //
  // if curvature_between_faces is less than sharp_smooth_curvature then
  // faces are connected_smooth
  //
  ContactFace<Real>* face_1          = cnfi1->Face(); // interaction face
  ContactFace<Real>* face_2          = cnfi2->Face(); // interaction face
  Real* normal_1               = face_1->Variable(FACE_NORMAL);
  Real* normal_2               = face_2->Variable(FACE_NORMAL);
  Real cos_angle_between_faces = Dot(normal_1, normal_2);
  Real face_curvature = 1.0-cos_angle_between_faces;
  if( face_curvature < sharp_smooth_curvature ) {
    return true;
  } else {
    return false;
  }
}


ContactNodeFaceInteraction*
ContactSearch::Compete_Interactions_Connected( ContactNodeFaceInteraction* cnfi1,
					       ContactNodeFaceInteraction* cnfi2 )
{
  // Note: It is assumed that interaction 1 is the stored interaction and
  //       interaction 2 is the candidate that may replace interaction 1.
  //       This matters because in one branch we are averaging normal and
  //       we need to know who has the previously averaged values.
  
  PRECONDITION( cnfi1 && cnfi2 );
  PRECONDITION( cnfi1->Node() == cnfi2->Node() );
  PRECONDITION( search_topology->Faces_Connected(cnfi1->Face(),cnfi2->Face()));
  
#ifdef CONTACT_DEBUG_NODE
  ContactNode<Real>* dn = cnfi1->Node();
  bool PRINT_THIS_NODE = primary_topology->Is_a_Debug_Node( dn );
  if( PRINT_THIS_NODE ){
    postream << "      Compete_Interactions_Connected\n";
  }
#endif

  VariableHandle CENTROID    = search_topology->Variable_Handle( ContactTopology::Centroid );
  VariableHandle FACE_NORMAL = search_topology->Variable_Handle( ContactTopology::Face_Normal );  

  const int INSIDE  = 1;
  const int OUTSIDE = 2;
  //
  //**********************************************************
  // Interaction location data
  //**********************************************************
  //
  int location_1      = cnfi1->Location();
  int location_2      = cnfi2->Location();

  // check that both locations are set and are in the proper range.
  PRECONDITION( (location_1 == INSIDE || location_1 == OUTSIDE ) &&
		(location_2 == INSIDE || location_2 == OUTSIDE ));
    
  //***********************************************************
  // Logic for interactions 1 and 2 INSIDE face
  //**********************************************************
  if( location_1 == INSIDE ) {
    if( location_2 == INSIDE ) {
      Real socn_1 = Strength_Contact(cnfi1, FACE_NORMAL);
      Real socn_2 = Strength_Contact(cnfi2, FACE_NORMAL);
      //
      // Check if the two are the same to within roundoff
      //
      Real max_mag = std::max( std::fabs(socn_1), std::fabs(socn_2) );
      if( max_mag ){
	Real relative_diff = std::fabs((socn_1 - socn_2)/max_mag);
	if( relative_diff < 1.E-14 ){
	  return cnfi1;
	}
      }
      if( socn_1 > socn_2 ) {
	return cnfi1;
      } else {
	return cnfi2;
      }
  //********************************************************************
  // Logic for interaction 1 INSIDE face and interaction 2 OUTSIDE face
  //********************************************************************
    } else {
      if( Connected_Smooth(cnfi1,cnfi2, sharp_smooth_curvature, FACE_NORMAL) ) {
 	return cnfi1;
      } else {
        Real socn_1 = Strength_Contact(cnfi1, FACE_NORMAL);
        Real socn_2 = Strength_Contact(cnfi2, FACE_NORMAL);
	if( socn_1 >= socn_2 ) {
	  return cnfi1;
	} else {
	  return cnfi2;
	}
      }
    }
  }
  if( location_1 == OUTSIDE ) {
    if( location_2 == INSIDE ) {
      //********************************************************************
      // Logic for interaction 2 INSIDE face and interaction 1 OUTSIDE face
      //********************************************************************
      if( Connected_Smooth(cnfi1, cnfi2, sharp_smooth_curvature, FACE_NORMAL) ) {
	return cnfi2;
      } else {
        Real socn_1 = Strength_Contact(cnfi1, FACE_NORMAL);
        Real socn_2 = Strength_Contact(cnfi2, FACE_NORMAL);
	if( socn_1 > socn_2 ) {
	  return cnfi1;
	} else {
	  return cnfi2;
	}
      }
    } else {
      //********************************************************************
      // Logic for interactions 1 and 2 OUTSIDE
      //********************************************************************
      Real* pbdir_1 = cnfi1->Get_Pushback_Dir();
      Real* pbdir_2 = cnfi2->Get_Pushback_Dir();
      //
      // indicator of sign on stored and candidate gaps
      //  + both penetrating or both non-penetrating
      //  - one penetrating one non-penetrating
      //
      Real gap_1          = cnfi1->Scalar_Var(ContactNodeFaceInteraction::GAP_CUR);
      Real gap_2          = cnfi2->Scalar_Var(ContactNodeFaceInteraction::GAP_CUR);
      Real pctps = gap_1*gap_2;
      //
      // indicators for check if contact points are same
      //
      Real discp = pbdir_1[0]*pbdir_2[0] +
  	           pbdir_1[1]*pbdir_2[1] +
	           pbdir_1[2]*pbdir_2[2];
      Real onemeps = 1.-1.e-6; // tolerance for same point check

      if( std::fabs(discp) > onemeps ) {
	//********************************************************************
	// Logic for interactions 1 and 2 going to same point
	//********************************************************************
	if( pctps == 0 && Concavity(cnfi1, cnfi2, CENTROID, FACE_NORMAL) <= 0 ) { // convex, one point on surface
	  // return interaction for smallest gap
	  if( gap_1 <= gap_2 ) {
	    return cnfi1; // gap_1 is smallest, return interaction 1
	  } else {
	    return cnfi2; // gap_2 is smallest, return interaction 2
	  }
	}
	if( pctps == 0 && Concavity(cnfi1, cnfi2, CENTROID, FACE_NORMAL) > 0 ) { // concave, one point on surface
	  // return interaction for smallest gap
	  if( gap_1 <= gap_2 ) {
	    return cnfi1; // gap_1 is smallest, return interaction 1
	  } else {
	    return cnfi2; // gap_2 is smallest, return interaction 2
	  }
	}
	if( pctps > 0 && Concavity(cnfi1, cnfi2, CENTROID, FACE_NORMAL) <= 0 ) { // convex, both pen. or non-pen.
	  if( Connected_Smooth(cnfi1, cnfi2, sharp_smooth_curvature, FACE_NORMAL) ) {
	    // if connected smooth, keep closest candidate
	    if( std::fabs(gap_1) <= std::fabs(gap_2) ) {
	      return cnfi1; // gap_1 is closest, return interaction 1
	    } else {
	      return cnfi2; // gap_2 is closest, return interaction 2
	    }
	  } else { 
	    // if not smooth connected chose "most-opposed"
            Real socn_1 = Strength_Contact(cnfi1, FACE_NORMAL);
            Real socn_2 = Strength_Contact(cnfi2, FACE_NORMAL);
	    if( socn_1 >= socn_2 ) {
	      return cnfi1; // interaction 1 has greatest normal soc
	    } else {
	      return cnfi2; // interaction 2 has greatest normal soc
	    }
	  }
	}
	if( pctps > 0 && Concavity(cnfi1, cnfi2, CENTROID, FACE_NORMAL) > 0 ) { // concave, both pen. or non-pen.
	  // return interaction for closest candidate
	  if( std::fabs(gap_1) <= std::fabs(gap_2) ) {
	    return cnfi1; // gap_1 is closest, return interaction 1
	  } else {
	    return cnfi2; // gap_2 is closest, return interaction 2
	  }
	}
	if( pctps < 0 && Concavity(cnfi1, cnfi2, CENTROID, FACE_NORMAL) <= 0 ) { // convex, one pen. one non-pen.
	  // chose penetrating constraint
	  if( gap_1 <= gap_2 ) {
	    return cnfi1; // gap_1 is penetrating, return interaction 1
	  } else {
	    return cnfi2; // gap_2 is penetrating, return interaction 2
	  }
	}
	if( pctps < 0 && Concavity(cnfi1, cnfi2, CENTROID, FACE_NORMAL) > 0 ) { // concave, one pen. one non-pen.
	  // chose penetrating constraint
	  if( gap_1 <= gap_2 ) {
	    return cnfi1; // gap_1 is penetrating, return interaction 1
	  } else {
	    return cnfi2; // gap_2 is penetrating, return interaction 2
	  }
	}   
      } else {
	//********************************************************************
	// Logic for interactions 1 and 2 going to different points
	//********************************************************************
	if( Connected_Smooth(cnfi1, cnfi2, sharp_smooth_curvature, FACE_NORMAL) ) {
	  // return interaction for closest candidate
	  if( std::fabs(gap_1) <= std::fabs(gap_2) ) {
	    return cnfi1; // gap_1 is closest, return interaction 1
	  } else {
	    return cnfi2; // gap_2 is closest, return interaction 2
	  }
	} else {
	  // if not smooth connected chose "most-opposed"
          Real socn_1 = Strength_Contact(cnfi1, FACE_NORMAL);
          Real socn_2 = Strength_Contact(cnfi2, FACE_NORMAL);
	  if( socn_1 >= socn_2 ) {
	    return cnfi1; // interaction 1 most-opposed, return 1
	  } else {
	    return cnfi2; // interaction 2 most-opposed, return 2
	  }
	}
      }
    }
  }
  return cnfi1;
}

ContactNodeEntityInteraction*
ContactSearch::Compete_Interactions_Unconnected( 
				VariableHandle NODE_COORDS,
				ContactNodeEntityInteraction* cefi1,
				ContactNodeEntityInteraction* cefi2 )
{
  // Note: It is assumed that interaction 1 is the stored interaction and
  //       interaction 2 is the candidate that may replace interaction 1.
  //       This matters because in one branch we are averaging normal and
  //       we need to know who has the previously averaged values.
  
  PRECONDITION( cefi1 && cefi2 );
  PRECONDITION( cefi1->Node() == cefi2->Node() );
  
#ifdef CONTACT_DEBUG_NODE
  ContactNode<Real>* dn = cefi1->Node();
  bool PRINT_THIS_NODE = primary_topology->Is_a_Debug_Node( dn );
  if( PRINT_THIS_NODE ){
    postream << "Compete_Interactions_Unconnected\n";
  }
#endif

  //
  //  Extract the current face normal for the interaction
  //
  Real normal_1[3];
  Real normal_2[3];
  cefi1->Get_Face_Normal(normal_1, search_topology);
  cefi2->Get_Face_Normal(normal_2, search_topology);

  VariableHandle NODE_NORMAL = search_topology->Variable_Handle( ContactTopology::Node_Normal );
  //
  //**********************************************************
  // Data for interaction 1
  //**********************************************************
  //
  Real* sn_normal_1 = cefi1->Vector_Var(ContactNodeFaceInteraction::PHYSICAL_FACE_NORMAL);
  Real* nd_normal_1 = cefi1->Node()->Variable(NODE_NORMAL);
  //
  //**********************************************************
  // Data for interaction 2
  //**********************************************************
  //
  Real* sn_normal_2 = cefi2->Vector_Var(ContactNodeFaceInteraction::PHYSICAL_FACE_NORMAL);
  Real* nd_normal_2 = cefi2->Node()->Variable(NODE_NORMAL);
  //
  //**********************************************************
  // calculate interpolated average face normals
  //**********************************************************
  //
  Real  int_face_norm_1[] = {0.0, 0.0, 0.0};
  Real  int_face_norm_2[] = {0.0, 0.0, 0.0};

  // interpolated (average) master face normals
  cefi1->Get_Avg_Face_Normal(int_face_norm_1, search_topology, NODE_COORDS);
  cefi2->Get_Avg_Face_Normal(int_face_norm_2, search_topology, NODE_COORDS);
  
  // Strength of Contact, normal   (socn) = <master_face_normal, slave_node_physical_face_normal>
  // Strength of Contact, tangent  (soct) = <master_face_average_normal, slave_node_normal>
  // Strength of Contact, combined (soc)  = socn + 0.5*soct

  // strength of contact for interaction 1
  Real socn_1         = -(normal_1[0]*sn_normal_1[0] +
			  normal_1[1]*sn_normal_1[1] +
			  normal_1[2]*sn_normal_1[2] );        // normal soc
  Real soct_1         =  (int_face_norm_1[0]*nd_normal_1[0] +
			  int_face_norm_1[1]*nd_normal_1[1] +
			  int_face_norm_1[2]*nd_normal_1[2] ); // tangential soc
  Real combined_soc_1 = socn_1 + 0.5*soct_1;                   // combined soc
  
  // strength of contact for interaction 2
  Real socn_2         = -(normal_2[0]*sn_normal_2[0] +
		  	  normal_2[1]*sn_normal_2[1] +
			  normal_2[2]*sn_normal_2[2] );        // normal soc
  Real soct_2         =  (int_face_norm_2[0]*nd_normal_2[0] +
			  int_face_norm_2[1]*nd_normal_2[1] +
			  int_face_norm_2[2]*nd_normal_2[2] ); // tangential soc
  Real combined_soc_2 = socn_2 + 0.5*soct_2;                   // combined soc

  if( combined_soc_2 > combined_soc_1 ) {
#ifdef CONTACT_DEBUG_NODE
    if( PRINT_THIS_NODE ){
      postream << "        Choosing candidate (largest combined soc)\n";
      postream << "           socn_1 =" << socn_1 << "\n";
      postream << "           soct_1 =" << soct_1 << "\n";
      postream << "             normal_1 =" << normal_1[0] 
               << " " << normal_1[1] << " " << normal_1[2] << "\n";
      postream << "             sn_normal_1 =" << sn_normal_1[0] 
               << " " << sn_normal_1[1] << " " << sn_normal_1[2] << "\n";
      postream << "             int_face_norm_1 =" << int_face_norm_1[0] 
             << " " << int_face_norm_1[1] << " " << int_face_norm_1[2] << "\n";
      postream << "             nd_normal_1 =" << nd_normal_1[0] 
               << " " << nd_normal_1[1] << " " << nd_normal_1[2] << "\n";
      postream << "           socn_2 =" << socn_2 << "\n";
      postream << "           soct_2 =" << soct_2 << "\n";
      postream << "             normal_2 =" << normal_2[0] 
               << " " << normal_2[1] << " " << normal_2[2] << "\n";
      postream << "             sn_normal_2 =" << sn_normal_2[0] 
               << " " << sn_normal_2[1] << " " << sn_normal_2[2] << "\n";
      postream << "             int_face_norm_2 =" << int_face_norm_2[0] 
             << " " << int_face_norm_2[1] << " " << int_face_norm_2[2] << "\n";
      postream << "             nd_normal_2 =" << nd_normal_2[0] 
               << " " << nd_normal_2[1] << " " << nd_normal_2[2] << "\n";
    }
#endif
    return cefi2;
  } else {
#ifdef CONTACT_DEBUG_NODE
    if( PRINT_THIS_NODE ){
      postream << "        Choosing stored (largest combined soc)\n";
      postream << "           socn_1 =" << socn_1 << "\n";
      postream << "           soct_1 =" << soct_1 << "\n";
      postream << "             normal_1 =" << normal_1[0] 
               << " " << normal_1[1] << " " << normal_1[2] << "\n";
      postream << "             sn_normal_1 =" << sn_normal_1[0] 
               << " " << sn_normal_1[1] << " " << sn_normal_1[2] << "\n";
      postream << "             int_face_norm_1 =" << int_face_norm_1[0] 
             << " " << int_face_norm_1[1] << " " << int_face_norm_1[2] << "\n";
      postream << "             nd_normal_1 =" << nd_normal_1[0] 
               << " " << nd_normal_1[1] << " " << nd_normal_1[2] << "\n";
      postream << "           socn_2 =" << socn_2 << "\n";
      postream << "           soct_2 =" << soct_2 << "\n";
      postream << "             normal_2 =" << normal_2[0] 
               << " " << normal_2[1] << " " << normal_2[2] << "\n";
      postream << "             sn_normal_2 =" << sn_normal_2[0] 
               << " " << sn_normal_2[1] << " " << sn_normal_2[2] << "\n";
      postream << "             int_face_norm_2 =" << int_face_norm_2[0] 
             << " " << int_face_norm_2[1] << " " << int_face_norm_2[2] << "\n";
      postream << "             nd_normal_2 =" << nd_normal_2[0] 
               << " " << nd_normal_2[1] << " " << nd_normal_2[2] << "\n";
    }
#endif
    return cefi1;
  }
}

void ContactSearch::Register_Timers()
{
#ifdef CONTACT_TIMINGS
  contact_time                         = timer.Register_Timer( CString("Total Contact") );
  search_time                          = timer.Register_Timer( CString("Total Search") );
  search_setup_time                    = timer.Register_Timer( CString("  Search Setup") );
  search_update_state_time             = timer.Register_Timer( CString("    Update State") );
  search_remaining_gap_time            = timer.Register_Timer( CString("    Remaining Gap") );
  augmented_config_time                = timer.Register_Timer( CString("  Augmented Configuration") );
  geom_time                            = timer.Register_Timer( CString("  Geometry Update") );
  context_time                         = timer.Register_Timer( CString("  Context Time") );
  
  face_charlen_geom_time               = timer.Register_Timer( CString("    Face Characteristic Len") );
  face_normal_geom_time                = timer.Register_Timer( CString("    Face Normals") );
  node_normal_geom_time                = timer.Register_Timer( CString("    Node Normals") );
  edge_curve_geom_time                 = timer.Register_Timer( CString("    Edge Curvature") );
  edge_curve_geom_time_serial          = timer.Register_Timer( CString("      Edge Curvature Serial") );
  edge_smooth_geom_time                = timer.Register_Timer( CString("    Edge Smoothing") );
  shell_loft_geom_time                 = timer.Register_Timer( CString("    Shell Lofting") );
  analytic_surf_geom_time              = timer.Register_Timer( CString("    Analytic Surfaces") );
  elem_geom_time                       = timer.Register_Timer( CString("    Element Volume/Centroids") );
  
  track_search_time                    = timer.Register_Timer( CString("  Tracked Search") );
  if( contact_number_of_processors( SearchComm ) > 1 ){
  track_ghost_time                     = timer.Register_Timer( CString("    TS Ghosting") );
  }
  track_scratch_time                   = timer.Register_Timer( CString("    TS Get Scratch") );
  track_gap_time                       = timer.Register_Timer( CString("    TS Compute Old Gaps") );
  track_retr_tracked_time              = timer.Register_Timer( CString("    TS Retrieve Tracked") );
  track_id_time                        = timer.Register_Timer( CString("    TS Interaction Definition") );
  track_release_scratch_time           = timer.Register_Timer( CString("    TS Release Scratch") );
  track_cleanup_time                   = timer.Register_Timer( CString("    TS Release Ghosts") );
  global_search_time                   = timer.Register_Timer( CString("  Global Search") );
  if( contact_number_of_processors( SearchComm ) > 1 ){
  create_secondary_time                = timer.Register_Timer( CString("    Create Secondary/Ghosting") );
#ifdef CONTACT_HSFC
  rcb_time                             = timer.Register_Timer( CString("      Zoltan HSFC") );
#else
  rcb_time                             = timer.Register_Timer( CString("      Zoltan RCB") );
#endif
  secondary_owner_migration_time       = timer.Register_Timer(CString("      Toplevel Copy") );
  owner_help_migrate_time              = timer.Register_Timer(CString("        TL Zoltan Migrate"));
  owner_help_migrate_size_time         = timer.Register_Timer(CString("          TL Migrate Size"));
  owner_help_migrate_pack_time         = timer.Register_Timer(CString("          TL Migrate Pack"));
  owner_help_migrate_unpack_time       = timer.Register_Timer(CString("          TL Migrate Unpack"));

  secondary_connect_mesh_time          = timer.Register_Timer( CString("      Finalize") );
  secondary_connect_mesh_phase2_time   = timer.Register_Timer( CString("        Build Lists") );
  secondary_connect_mesh_phase3_time   = timer.Register_Timer( CString("        Set Ownership") );
  secondary_connect_mesh_phase4_time   = timer.Register_Timer( CString("        Topology Connections") );
  secondary_connect_mesh_phase5_time   = timer.Register_Timer( CString("        Interaction Connections") );
  }
  on_processor_search_time             = timer.Register_Timer( CString("    On-Processor Search") );
  search_scratch_time                  = timer.Register_Timer( CString("      OPS Get Scratch") );
  search_gap_time                      = timer.Register_Timer( CString("      OPS Compute Old Gaps") );
  search_phys_face_time                = timer.Register_Timer( CString("      OPS Build Physical Faces") );
  search_retr_tied_time                = timer.Register_Timer( CString("      OPS Retrieve Tied") );
  search_serial_setup_time             = timer.Register_Timer( CString("      OPS Set Up") );
  search_main_time                     = timer.Register_Timer( CString("      OPS Detailed Searching") );
  search_id_time                       = timer.Register_Timer( CString("      OPS Interaction Definition") );
  search_release_scratch_time          = timer.Register_Timer( CString("      OPS Release Scratch") );
  if( contact_number_of_processors( SearchComm ) > 1 ){
  interaction_migration_time           = timer.Register_Timer( CString("    Interaction Copy") );
  interaction_onproc_copy_time         = timer.Register_Timer( CString("      On-processor Copy") );
  interaction_export_setup_time        = timer.Register_Timer( CString("      Export Setup") );
  interaction_help_migrate_time        = timer.Register_Timer( CString("      IC Zoltan Migrate") );
  interaction_help_migrate_size_time   = timer.Register_Timer( CString("        IC Migrate Size") );
  interaction_help_migrate_pack_time   = timer.Register_Timer( CString("        IC Migrate Pack") );
  interaction_help_migrate_unpack_time =  timer.Register_Timer( CString("        IC Migrate Unpack") );
  }
  restore_constraint_continuity_time   = timer.Register_Timer( CString("    Restore Constraint Continuity") );
  if( contact_number_of_processors( SearchComm ) > 1 ){
  cleanup_secondary_time               = timer.Register_Timer( CString("    Clean Up Secondary") );
  }
#ifdef CONTACT_TIMINGS1
  baseline_constructors_time = 
    timer.Register_Timer( CString("  Constructors") );
  baseline_node_time = 
    timer.Register_Timer( CString("    Nodes") );
  baseline_edge_time = 
    timer.Register_Timer( CString("    Edges") );
  baseline_face_time = 
    timer.Register_Timer( CString("    Faces") );
  baseline_size_time = 
    timer.Register_Timer( CString("  Size") );
  baseline_pack_time = 
    timer.Register_Timer( CString("  Pack") );
  baseline_unpack_time = 
    timer.Register_Timer( CString("  Unpack") );
  baseline_size_all_time = 
    timer.Register_Timer( CString("  Size_All") );
  baseline_pack_all_time = 
    timer.Register_Timer( CString("  Pack_All") );
  baseline_unpack_all_time = 
    timer.Register_Timer( CString("  Unpack_All") );
#endif
#endif
}

//
//  Loop through all analytic surfaces.  Detect and compete, and apply any interactions found.
//  Interactions will be stored as general node-entity interactions
//
void ContactSearch::Process_Analytic_Surfaces(
                                              ObjectBoundingBoxHierarchy *node_hierarchy,
                                              ContactNode<Real>** Nodes,
                                              int interaction_source,
                                              VariableHandle position_var,
                                              VariableHandle position_var_2
                                              ) {
  if(search_topology->Number_of_Analytic_Surfaces() == 0) return;
  VariableHandle NODE_NORMAL = search_topology->Variable_Handle( ContactTopology::Node_Normal );
  int *list = new int[search_topology->Number_of_Nodes()];

  Real xmin[MAX_DIMENSIONALITY];
  Real xmax[MAX_DIMENSIONALITY];

  // Process any analytic Surfaces
  for( int i=0 ; i<search_topology->Number_of_Analytic_Surfaces() ; ++i ){
    ContactAnalyticSurface* AnalyticSurface = 
      search_topology->Analytic_Surface( i );
      
#define NO_CONTACT_USE_ACTUAL_BB
#ifdef CONTACT_USE_ACTUAL_BB

    ContactBoundingBox* surface_box = AnalyticSurface->BoundingBox();
    
    //
    //  Compute the nodes that can interact with the surface
    //
    int list_size = 0;
    if(node_hierarchy != NULL) node_hierarchy->search_for_overlap_recurse(*surface_box, list, list_size);
#else

    AnalyticSurface->Bounding_Box( xmin,xmax );
    ContactBoundingBox surface_box;
    surface_box.add_point(xmin);
    surface_box.add_point(xmax);
    
    //
    //  Compute the nodes that can interact with the surface
    //
    int list_size = 0;
    if(node_hierarchy != NULL) ObjectBoundingBoxHierarchy::search_for_overlap_loop(node_hierarchy,surface_box, list, list_size);
#endif
    
    for( int j=0 ; j<list_size ; ++j ) {
      ContactNode<Real> *node = Nodes[list[j]];
      
      bool PRINT_THIS_NODE = false;
#ifdef CONTACT_DEBUG_NODE
      PRINT_THIS_NODE = primary_topology->Is_a_Debug_Node( node );
#endif

      Real* position  = node->Variable(position_var);
      Real* position2 = node->Variable(position_var_2);
      Real p_mag;
      Real cpoint[3],snorm[3];
      Real pushback_dir[3];
      Real time_to_contact;
      int  location;
      
      bool possible_interaction;
      
      if(position_var_2 == -1){
        //do the single position search
        possible_interaction = AnalyticSurface->Process( position, p_mag, cpoint, snorm, pushback_dir, time_to_contact, location, PRINT_THIS_NODE, &postream );
      }else{
        //do the two position search
        possible_interaction = AnalyticSurface->Process( position, position2, p_mag, cpoint, snorm, pushback_dir, time_to_contact, location, PRINT_THIS_NODE, &postream );
      }
      
      if( possible_interaction == true){
	int node_key = node->Entity_Key();
	int physical_face=-1;
        int num_physical_faces;
        Real* pf_normals[3];
	if( node_key == -1 ){
	  // Determine which node entity key and physical face to use
	  num_physical_faces = (int)*(NUMBER_PHYSICAL_FACES.Get_Scratch(node->ProcArrayIndex()));
	  POSTCONDITION( num_physical_faces > 0 );
	  Real most_opposed = 2.0;
	  pf_normals[0] = PHYSICAL_FACE_NORMAL_1.Get_Scratch(node->ProcArrayIndex());
	  pf_normals[1] = PHYSICAL_FACE_NORMAL_2.Get_Scratch(node->ProcArrayIndex());
	  pf_normals[2] = PHYSICAL_FACE_NORMAL_3.Get_Scratch(node->ProcArrayIndex());
          Find_Physical_Face(num_physical_faces, pf_normals, snorm, physical_face, most_opposed);
          node_key = node->GetFacePFEntityKey(physical_face);
	} else {
	  // Determine which physical face to use
	  num_physical_faces = (int)*(NUMBER_PHYSICAL_FACES.Get_Scratch(node->ProcArrayIndex()));
	  if( num_physical_faces > 0 ) {
	    Real most_opposed = 2.0;
	    pf_normals[0] = PHYSICAL_FACE_NORMAL_1.Get_Scratch(node->ProcArrayIndex());
	    pf_normals[1] = PHYSICAL_FACE_NORMAL_2.Get_Scratch(node->ProcArrayIndex());
	    pf_normals[2] = PHYSICAL_FACE_NORMAL_3.Get_Scratch(node->ProcArrayIndex());
            Find_Physical_Face(num_physical_faces, pf_normals, snorm, physical_face, most_opposed);
          } else {
	    physical_face = 0;
            pf_normals[0] = node->Variable(NODE_NORMAL);
          }
	}
        POSTCONDITION(physical_face>=0);
	int Interaction_Type = (int)
	  search_data->Get_Search_Data( ContactSearch::INTERACTION_TYPE,
					node_key, 
					AnalyticSurface->Entity_Key() );
	if( Interaction_Type != NO_INTERACTION ){
          //if this is a dynamic search with only one physical face
          if((position_var_2 != -1 && num_physical_faces == 1) && 
             (Interaction_Type == CPP_INTERACTION ||
              node->Get_NodeEntity_Interaction( 0, 1 )))
          {
           interaction_source = ContactNodeEntityInteraction::CLOSEST_POINT_PROJECTION_2;
          }
          int action = 0;
	  bool is_tied=(search_data->Get_Search_Data(INTERACTION_TYPE, 
					   node_key,
					   AnalyticSurface->Entity_Key())
		      == TIED_INTERACTION ) ? true : false;
	  bool is_infSlip=(search_data->Get_Search_Data(INTERACTION_TYPE, 
					   node_key,
					   AnalyticSurface->Entity_Key())
		      == INFINITESIMAL_SLIP_INTERACTION ) ? true : false;
	  bool is_glued=(search_data->Get_Search_Data(INTERACTION_TYPE, 
					   node_key,
					   AnalyticSurface->Entity_Key())
		      == GLUED_INTERACTION ) ? true : false;
          bool is_tracked=(tracking_type==LOCAL_TRACKING);
          ContactNodeSurfaceInteraction* cnsi = 
            ContactNodeSurfaceInteraction::new_ContactNodeSurfaceInteraction( 
                             allocators[ALLOC_ContactNodeSurfaceInteraction],
                             (ContactNodeEntityInteraction::InteractionSource)interaction_source,
                             node,AnalyticSurface,node_key,
                             pf_normals[physical_face],p_mag,cpoint,snorm,pushback_dir,time_to_contact,
                             location, is_tied, is_infSlip, is_glued, is_tracked );
          Process_Interaction(position_var, cnsi, action);
        }
      }
    }
  }

  delete [] list;
}



//
//  Loop through all analytic surfaces.  Detect and compete, and apply any interactions found.
//  Interactions will be stored as general node-entity interactions
//
void ContactSearch::Process_Analytic_Surfaces_New(
                                              ObjectBoundingBoxHierarchy **node_hierarchy_ptrs,
                                              ContactNode<Real> **Nodes,
                                              int interaction_source,
                                              VariableHandle position_var,
                                              VariableHandle position_var_2) {
  if(search_topology->Number_of_Analytic_Surfaces() == 0) return;
  VariableHandle NODE_NORMAL = search_topology->Variable_Handle( ContactTopology::Node_Normal );
  int *list = new int[search_topology->Number_of_Nodes()];

  Real xmin[MAX_DIMENSIONALITY];
  Real xmax[MAX_DIMENSIONALITY];

  int num_entities = search_data->Num_Search_Entities();
  ACME::ContactNode_Vector node_list;

  // Process any analytic Surfaces
  for( int i=0 ; i<search_topology->Number_of_Analytic_Surfaces() ; ++i ){
    ContactAnalyticSurface* AnalyticSurface = search_topology->Analytic_Surface( i );
      
#define NO_CONTACT_USE_ACTUAL_BB
#ifdef CONTACT_USE_ACTUAL_BB
    ContactBoundingBox* surface_box = AnalyticSurface->BoundingBox();
    
    //
    //  Compute the nodes that can interact with the surface
    //
    int list_size = 0;    
    for(int iblock = 0; iblock < num_entities; ++iblock) {
      int Interaction_Type = (int)search_data->Get_Search_Data( ContactSearch::INTERACTION_TYPE, iblock, 
					                        AnalyticSurface->Entity_Key() );
      if( Interaction_Type == NO_INTERACTION ) continue;

      ObjectBoundingBoxHierarchy *node_block_hierarchy = node_hierarchy_ptrs[iblock+1];
      if(node_block_hierarchy != NULL) ObjectBoundingBoxHierarchy::search_for_overlap_loop(node_block_hierarchy, surface_box, list, list_size);     
    }

    ObjectBoundingBoxHierarchy *node_block_hierarchy = node_hierarchy_ptrs[0];
    if(node_block_hierarchy != NULL) ObjectBoundingBoxHierarchy::search_for_overlap_loop(node_block_hierarchy, surface_box, list, list_size);     
#else
    AnalyticSurface->Bounding_Box( xmin,xmax );
    ContactBoundingBox surface_box;
    surface_box.add_point(xmin);
    surface_box.add_point(xmax);
    
    //
    //  Compute the nodes that can interact with the surface
    //
    //
    //  Compute the nodes that can interact with the surface
    //
    node_list.empty();
    for(int iblock = 0; iblock < num_entities; ++iblock) {
      int Interaction_Type = (int)search_data->Get_Search_Data( ContactSearch::INTERACTION_TYPE, iblock, 
					                        AnalyticSurface->Entity_Key() );
      if( Interaction_Type == NO_INTERACTION ) continue;

      ObjectBoundingBoxHierarchy *node_block_hierarchy = node_hierarchy_ptrs[iblock+1];
      if(node_block_hierarchy != NULL) {
        int list_size = 0;    
        ObjectBoundingBoxHierarchy::search_for_overlap_loop(node_block_hierarchy, surface_box, list, list_size);
        for(int j = 0; j < list_size; ++j) {
          node_list.push_back(Nodes[list[j]]);
	}
      }
    }

    ObjectBoundingBoxHierarchy *node_block_hierarchy = node_hierarchy_ptrs[0];
    if(node_block_hierarchy != NULL) {
      int list_size = 0;    
      ObjectBoundingBoxHierarchy::search_for_overlap_loop(node_block_hierarchy, surface_box, list, list_size);
      for(int j = 0; j < list_size; ++j) {
        node_list.push_back(Nodes[list[j]]);
      }
    }
#endif
     
    int num_nodes = node_list.size();
    for( int j=0 ; j<num_nodes ; ++j ) {
      ContactNode<Real> *node = node_list[j];

      bool PRINT_THIS_NODE = false;
#ifdef CONTACT_DEBUG_NODE
      PRINT_THIS_NODE = primary_topology->Is_a_Debug_Node( node );
#endif

      Real* position  = node->Variable(position_var);
      Real* position2 = node->Variable(position_var_2);
      Real p_mag;
      Real cpoint[3],snorm[3];
      Real pushback_dir[3];
      Real time_to_contact;
      int  location;
      bool possible_interaction;
      
      if(position_var_2 == -1){
        //do the single position search
        possible_interaction = AnalyticSurface->Process( position, p_mag, cpoint, snorm, pushback_dir, time_to_contact, location, PRINT_THIS_NODE, &postream );
      }else{
        //do the two position search
        possible_interaction = AnalyticSurface->Process( position, position2, p_mag, cpoint, snorm, pushback_dir, time_to_contact, location, PRINT_THIS_NODE, &postream );
      }

      if( possible_interaction == true){
	int node_key = node->Entity_Key();
	int physical_face=-1;
        int num_physical_faces;
        Real* pf_normals[3];
	if( node_key == -1 ){
	  // Determine which node entity key and physical face to use
	  num_physical_faces = (int)*(NUMBER_PHYSICAL_FACES.Get_Scratch(node->ProcArrayIndex()));
	  POSTCONDITION( num_physical_faces > 0 );
	  Real most_opposed = 2.0;
	  pf_normals[0] = PHYSICAL_FACE_NORMAL_1.Get_Scratch(node->ProcArrayIndex());
	  pf_normals[1] = PHYSICAL_FACE_NORMAL_2.Get_Scratch(node->ProcArrayIndex());
	  pf_normals[2] = PHYSICAL_FACE_NORMAL_3.Get_Scratch(node->ProcArrayIndex());
          Find_Physical_Face(num_physical_faces, pf_normals, snorm, physical_face, most_opposed);
          node_key = node->GetFacePFEntityKey(physical_face);
	} else {
	  // Determine which physical face to use
	  num_physical_faces = (int)*(NUMBER_PHYSICAL_FACES.Get_Scratch(node->ProcArrayIndex()));
	  if( num_physical_faces > 0 ) {
	    Real most_opposed = 2.0;
	    pf_normals[0] = PHYSICAL_FACE_NORMAL_1.Get_Scratch(node->ProcArrayIndex());
	    pf_normals[1] = PHYSICAL_FACE_NORMAL_2.Get_Scratch(node->ProcArrayIndex());
	    pf_normals[2] = PHYSICAL_FACE_NORMAL_3.Get_Scratch(node->ProcArrayIndex());
            Find_Physical_Face(num_physical_faces, pf_normals, snorm, physical_face, most_opposed);
          } else {
	    physical_face = 0;
            pf_normals[0] = node->Variable(NODE_NORMAL);
          }
	}
        POSTCONDITION(physical_face>=0);
	int Interaction_Type = (int)
	  search_data->Get_Search_Data( ContactSearch::INTERACTION_TYPE,
					node_key, 
					AnalyticSurface->Entity_Key() );
	if( Interaction_Type != NO_INTERACTION ){
          if((position_var_2 != -1 && num_physical_faces == 1) && 
             (Interaction_Type == CPP_INTERACTION ||
              node->Get_NodeEntity_Interaction( 0, 1 )))
          {
           interaction_source = ContactNodeEntityInteraction::CLOSEST_POINT_PROJECTION_2;
          }
          int action = 0;
	  bool is_tied=(search_data->Get_Search_Data(INTERACTION_TYPE, 
					   node_key,
					   AnalyticSurface->Entity_Key())
		      == TIED_INTERACTION ) ? true : false;
	  bool is_infSlip=(search_data->Get_Search_Data(INTERACTION_TYPE, 
					   node_key,
					   AnalyticSurface->Entity_Key())
		      == INFINITESIMAL_SLIP_INTERACTION ) ? true : false;
	  bool is_glued=(search_data->Get_Search_Data(INTERACTION_TYPE, 
					   node_key,
					   AnalyticSurface->Entity_Key())
		      == GLUED_INTERACTION ) ? true : false;
          bool is_tracked=(tracking_type==LOCAL_TRACKING);
          ContactNodeSurfaceInteraction* cnsi = 
            ContactNodeSurfaceInteraction::new_ContactNodeSurfaceInteraction( 
                             allocators[ALLOC_ContactNodeSurfaceInteraction],
                             (ContactNodeEntityInteraction::InteractionSource)interaction_source,
                             node,AnalyticSurface,node_key,
                             pf_normals[physical_face],p_mag,cpoint,snorm,pushback_dir,time_to_contact,
                             location,is_tied,is_infSlip,is_glued,is_tracked );
          Process_Interaction(position_var, cnsi, action);
        }
      }
    }
  }
  delete [] list;
}


//
//  Loop through all analytic surfaces.  Detect and compete, and apply any interactions found.
//  Interactions will be stored as general node-entity interactions
//
void ContactSearch::Process_Analytic_Surfaces(Real* position,
                                              int* index,
                                              int   number_of_nodes,
                                              Real* scratch,
                                              int* rank2, 
                                              int*  list,
                                              ContactNode<Real>** Nodes, int* node_map,
                                              int interaction_source,
                                              VariableHandle position_var_1,
                                              VariableHandle position_var_2) {
  Real xmin[MAX_DIMENSIONALITY];
  Real xmax[MAX_DIMENSIONALITY];
  VariableHandle NODE_NORMAL = search_topology->Variable_Handle( ContactTopology::Node_Normal );

  // Process any analytic Surfaces
  for( int i=0 ; i<search_topology->Number_of_Analytic_Surfaces() ; ++i ){
    ContactAnalyticSurface* AnalyticSurface = 
      search_topology->Analytic_Surface( i );
    AnalyticSurface->Bounding_Box( xmin,xmax );
    //
    //  Compute the nodes that can interact with the surface
    //
    int list_size;
    int ne = 1;
    Compute_Nodes_in_Box(ne, xmin, xmax, position, index, num_active_nodes, 
                         scratch, rank2, list_size, list);
    for( int j=0 ; j<list_size ; ++j ) {
      ContactNode<Real> *node = NULL;
      if (node_map) 
        node = Nodes[node_map[list[j]-1]];
      else
        node = Nodes[list[j]-1];
      if (node->Ownership() != ContactTopologyEntity<Real>::OWNED) continue;

      bool PRINT_THIS_NODE = false;
#ifdef CONTACT_DEBUG_NODE
      PRINT_THIS_NODE = primary_topology->Is_a_Debug_Node( node );
#endif

      Real* position1 = node->Variable(position_var_1);
      Real* position2 = node->Variable(position_var_2);
      Real p_mag;
      Real cpoint[3],snorm[3];
      Real pushback_dir[3];
      Real time_to_contact;
      int  location;
      bool possible_interaction;
      
      if(position_var_2 == -1){
        //do the single position search
        possible_interaction = AnalyticSurface->Process( position1, p_mag, cpoint, snorm, pushback_dir, time_to_contact, location, PRINT_THIS_NODE, &postream );
      }else{
        //do the two position search
        possible_interaction = AnalyticSurface->Process( position1, position2, p_mag, cpoint, snorm, pushback_dir, time_to_contact, location, PRINT_THIS_NODE, &postream );
      }
      
      if( possible_interaction == true ){
	int node_key = node->Entity_Key();
	int physical_face=-1;
        int num_physical_faces;
        Real* pf_normals[3];
	if( node_key == -1 ){
	  // Determine which node entity key and physical face to use
	  num_physical_faces = (int) *(NUMBER_PHYSICAL_FACES.Get_Scratch(node->ProcArrayIndex()));
	  POSTCONDITION( num_physical_faces > 0 );
	  POSTCONDITION( num_physical_faces < 4 );
	  Real most_opposed = 2.0;
	  pf_normals[0] = PHYSICAL_FACE_NORMAL_1.Get_Scratch(node->ProcArrayIndex());
	  pf_normals[1] = PHYSICAL_FACE_NORMAL_2.Get_Scratch(node->ProcArrayIndex());
	  pf_normals[2] = PHYSICAL_FACE_NORMAL_3.Get_Scratch(node->ProcArrayIndex());
          Find_Physical_Face(num_physical_faces, pf_normals, snorm, physical_face, most_opposed);
          PRECONDITION(physical_face>=0);
          PRECONDITION(physical_face<=3);
          node_key = node->GetFacePFEntityKey(physical_face);
	} else {
	  // Determine which physical face to use
	  num_physical_faces = (int)*(NUMBER_PHYSICAL_FACES.Get_Scratch(node->ProcArrayIndex()));
	  if( num_physical_faces > 0 ) {
	    Real most_opposed = 2.0;
	    pf_normals[0] = PHYSICAL_FACE_NORMAL_1.Get_Scratch(node->ProcArrayIndex());
	    pf_normals[1] = PHYSICAL_FACE_NORMAL_2.Get_Scratch(node->ProcArrayIndex());
	    pf_normals[2] = PHYSICAL_FACE_NORMAL_3.Get_Scratch(node->ProcArrayIndex());
            Find_Physical_Face(num_physical_faces, pf_normals, snorm, physical_face, most_opposed);
          } else {
	    physical_face = 0;
            pf_normals[0] = node->Variable(NODE_NORMAL);
          }
	}
        POSTCONDITION(physical_face>=0);
	int Interaction_Type = (int)
	  search_data->Get_Search_Data( ContactSearch::INTERACTION_TYPE,
					node_key, 
					AnalyticSurface->Entity_Key() );
	if( Interaction_Type != NO_INTERACTION ){
          if((position_var_2 != -1 && num_physical_faces == 1) && 
             (Interaction_Type == CPP_INTERACTION ||
              node->Get_NodeEntity_Interaction( 0, 1 )))
          {
           interaction_source = ContactNodeEntityInteraction::CLOSEST_POINT_PROJECTION_2;
          }
          int action = 0;
	  bool is_tied=(search_data->Get_Search_Data(INTERACTION_TYPE, 
					   node_key,
					   AnalyticSurface->Entity_Key())
		      == TIED_INTERACTION ) ? true : false;
	  bool is_infSlip=(search_data->Get_Search_Data(INTERACTION_TYPE, 
					   node_key,
					   AnalyticSurface->Entity_Key())
		      == INFINITESIMAL_SLIP_INTERACTION ) ? true : false;
	  bool is_glued=(search_data->Get_Search_Data(INTERACTION_TYPE, 
					   node_key,
					   AnalyticSurface->Entity_Key())
		      == GLUED_INTERACTION ) ? true : false;
          bool is_tracked=(tracking_type==LOCAL_TRACKING);
          ContactNodeSurfaceInteraction* cnsi = 
            ContactNodeSurfaceInteraction::new_ContactNodeSurfaceInteraction( 
                             allocators[ALLOC_ContactNodeSurfaceInteraction],
                             (ContactNodeEntityInteraction::InteractionSource)interaction_source,
                             node,AnalyticSurface,node_key,
                             pf_normals[physical_face],p_mag,cpoint,snorm,pushback_dir,time_to_contact,
                             location,is_tied,is_infSlip,is_glued,is_tracked );
          Process_Interaction(position_var_1, cnsi, action);
        }
      }
    }
  }
}

void ContactSearch::process_node_face_interactions(const int &list_size, 
                                                   ACME::ContactNode_Vector &list,
                                                   ContactFace<Real> *face, 
                                                   const VariableHandle &CALCULATION_POSITION,
                                                   ContactNode<Real> **Nodes,
                                                   const VariableHandle &NODE_NORMAL,
                                                   Real &REL_TANG_TOL_VAR,
                                                   Real &user_tang_tol,
                                                   Real &gap_tol,
                                                   int *physical_faces,
                                                   ACME::Int_Vector &node_keys,
                                                   int interaction_source){
  if( list_size ){
    int nfacets;
    face->FacetDecomposition(nfacets,ms_coordinates_c,
			     ms_normals_c,CALCULATION_POSITION);
    for( int k=0 ; k<list_size ; ++k ){
      ContactNode<Real> *node = list[k];
      int npairs = 1;
#ifdef CONTACT_DEBUG_NODE
      bool PRINT_THIS_NODE = primary_topology->Is_a_Debug_Node( node );
      if( PRINT_THIS_NODE ){
	postream << "  Processing Potential Interaction for Node "
		 << node->Exodus_ID() << " with Face "
		 << face->Global_ID() << "\n";
      }
#endif
      Real* coords_node = node->Variable(CALCULATION_POSITION);
      for (int nn=0; nn<nfacets; ++nn) {
	pushback_dir_flag[nn] = 2; // normal pushback
      }
      for (int nn=0; nn<nfacets*data_size; ++nn) {
	ctrcl_facets[nn] = 0.0;
      }
      if (dimensionality==2) {
        for (int nn=0; nn<nfacets; ++nn) {
	  FORTRAN(cnodeline_cpproj)(npairs, coords_node,
				    &ms_coordinates_c[nn*9],
				    &ms_normals_c[nn*3],
				    &pushback_dir_flag[nn],
				    ctrl,
				    &ctrcl_facets[nn*data_size],
                                    REL_TANG_TOL_VAR,user_tang_tol,gap_tol);
#ifdef CONTACT_DEBUG_NODE
	  if( PRINT_THIS_NODE ){
	    if( ctrcl_facets[nn*data_size] == 0 )
	      postream << "    Rejected with subtriangle " << nn << "\n";
	    else if( ctrcl_facets[nn*data_size] > 0 )
	      postream << "    Accepted with subtriangle " << nn << "\n";
	  }
#endif
        }
      } else if (dimensionality==3) {
        for (int nn=0; nn<nfacets; ++nn) {
	  FORTRAN(cnodetriangle_cpproj)(npairs, coords_node,
					&ms_coordinates_c[nn*9],
					&ms_normals_c[nn*3],
					&pushback_dir_flag[nn],
					ctrl,
					&ctrcl_facets[nn*data_size],
                                        REL_TANG_TOL_VAR,user_tang_tol,gap_tol);
#ifdef CONTACT_DEBUG_NODE
	  if( PRINT_THIS_NODE ){
	    if( ctrcl_facets[nn*data_size] == 0 )
	      postream << "    Rejected with subtriangle " << nn << "\n";
	    else if( ctrcl_facets[nn*data_size] > 0 )
	      postream << "    Accepted with subtriangle " << nn << "\n";
	  }
#endif
        }
      }
      face->FacetStaticRestriction(nfacets, ms_coordinates_c, ms_normals_c,
				   ctrcl_facets, ctrcl);
      if( ctrcl[0] == 1 ){
#ifdef CONTACT_DEBUG_NODE
	if( PRINT_THIS_NODE ){
	  postream << "	      Global Coords: (" 
		   << ctrcl[1] << "," << ctrcl[2] << "," 
		   << ctrcl[3] << ")\n";
	  postream << "	      Gap:	     " 
		   << ctrcl[4] << "\n";
	  postream << "	      Pushback Dir:  (" 
		   << ctrcl[5] << "," << ctrcl[6] << "," 
		   << ctrcl[7] << ")\n";
	  postream << "	      Normal Dir:    (" 
		   << ctrcl[8] << "," << ctrcl[9] << "," 
		   << ctrcl[10] << ")\n";
	  postream << "	      Location: "
		   << (int) ctrcl[11] << "\n";
	  postream << "	      Time - In/Out: "
		   << ctrcl[12] << "\n";
	  postream << "  Interaction Being Sent to Process Interaction\n";
	}
#endif
	// The interaction is valid
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
        POSTCONDITION(pf_norm!=0);
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
								     (ContactNodeFaceInteraction::InteractionSource)interaction_source,
								     node, face, ctrcl, node_keys[k],
								     pf_norm, is_tied, is_infSlip, is_glued, is_tracked,
                                                                     enable_off_face_tracking, CALCULATION_POSITION );
	POSTCONDITION(node==cnfi->Node());
        int action = 0;
	Process_Interaction( CALCULATION_POSITION, cnfi, action );
      }
    }
  }
}

void ContactSearch::Sort(int size, int* list)
{
  int i, j, k, n, entity;

  if (size>1) {
    k = (size>>1)+1;
    n = size;
    for (;;) {
      if (k>1) {
        entity = list[--k-1];
      } else {
        entity    = list[n-1];
        list[n-1] = list[0];
        if (--n == 1) {
          list[0] = entity;
          break;
        }
      }
      i = k;
      j = k<<1;
      while (j<=n) {
        if ((j<n) && (list[j-1]<list[j])) ++j;
        if (entity<list[j-1]) {
          list[i-1] = list[j-1];
          j += (i=j);
        } else {
          j = n+1;
        }
      }
      list[i-1] = entity;
    }
  }
}

void ContactSearch::Get_Scratch(const int number_of_nodes) {
  NUMBER_PHYSICAL_FACES.   Allocate_Scratch(number_of_nodes, 1);
  PHYSICAL_FACE_NORMAL_1.  Allocate_Scratch(number_of_nodes, 3);
  PHYSICAL_FACE_NORMAL_2.  Allocate_Scratch(number_of_nodes, 3);
  PHYSICAL_FACE_NORMAL_3.  Allocate_Scratch(number_of_nodes, 3);
}

void ContactSearch::Clear_Scratch() {
  NUMBER_PHYSICAL_FACES.Clear_Scratch();
  PHYSICAL_FACE_NORMAL_1.Clear_Scratch();
  PHYSICAL_FACE_NORMAL_2.Clear_Scratch();
  PHYSICAL_FACE_NORMAL_3.Clear_Scratch();
}

void Find_Physical_Face(int num_physical_faces, Real** PF_normals, Real * face_normal, int& physical_face, Real& most_opposed){
  most_opposed = 2.0;
  for( int i=0 ; i<num_physical_faces ; ++i ){
    Real dot = PF_normals[i][0]*face_normal[0] + 
               PF_normals[i][1]*face_normal[1] +
               PF_normals[i][2]*face_normal[2] ;
    if( dot < most_opposed ){
      physical_face = i;
      most_opposed = dot;
    }
  }
  
}

int ContactSearch::Number_Nodes_Per_Face( ContactFace_Type face_type) {
  return ContactFace<Real>::Nodes_Per_Face(face_type);
}

ContactSearch::ContactEdge_Type ContactSearch::Face_Edge_Type(ContactFace_Type face_type) {
  ContactSearch::ContactEdge_Type edge_type; 
  switch( face_type ){
  case ContactSearch::LINEFACEL2:
    edge_type = ContactSearch::NO_EDGES;
    break;
  case ContactSearch::LINEFACEQ3:
    edge_type = ContactSearch::NO_EDGES;
    break;
  case ContactSearch::QUADFACEL4:
    edge_type = ContactSearch::LINEEDGEL2;
    break;
  case ContactSearch::QUADFACEQ8:
    edge_type = ContactSearch::LINEEDGEQ3;
    break;
  case ContactSearch::QUADFACEQ9:
    edge_type = ContactSearch::LINEEDGEQ3;
    break;
  case ContactSearch::TRIFACEL3:
    edge_type = ContactSearch::LINEEDGEL2;
    break;
  case ContactSearch::TRIFACEQ6:
    edge_type = ContactSearch::LINEEDGEQ3;
    break;
  case ContactSearch::SHELLQUADFACEL4:
    edge_type = ContactSearch::LINEEDGEL2;
    break;
  case ContactSearch::SHELLTRIFACEL3:
    edge_type = ContactSearch::LINEEDGEL2;
    break;
  default:
    POSTCONDITION(0);
    edge_type = ContactSearch::NEDGE_TYPES;
    break;
  }
  return edge_type;
}

ContactFace<Real>* ContactSearch::New_ContactFace(ContactFace_Type face_type,   ContactFixedSizeAllocator* alloc) {
  ContactFace<Real>* face=NULL;
  switch( face_type ){
  case ContactSearch::QUADFACEL4 :
    face = ContactQuadFaceL4<Real>::new_ContactQuadFaceL4(alloc);
    break;
  case ContactSearch::QUADFACEQ8 :
    face = ContactQuadFaceQ8::new_ContactQuadFaceQ8(alloc);
    break;
  case ContactSearch::QUADFACEQ9 :
    face = ContactQuadFaceQ9::new_ContactQuadFaceQ9(alloc);
    break;
  case ContactSearch::TRIFACEL3 :
    face = ContactTriFaceL3<Real>::new_ContactTriFaceL3(alloc);
    break;
  case ContactSearch::TRIFACEQ6 :
    face = ContactTriFaceQ6::new_ContactTriFaceQ6(alloc);
    break;
  case ContactSearch::SHELLQUADFACEL4 :
    face = ContactShellQuadFaceL4<Real>::new_ContactShellQuadFaceL4(alloc);
    break;
  case ContactSearch::SHELLTRIFACEL3 :
    face = ContactShellTriFaceL3<Real>::new_ContactShellTriFaceL3(alloc);
    break;
  case ContactSearch::LINEFACEL2 :
    face = ContactLineFaceL2::new_ContactLineFaceL2(alloc);
    break;
  case ContactSearch::LINEFACEQ3 :
    face = ContactLineFaceQ3::new_ContactLineFaceQ3(alloc);
    break;
  default:
    POSTCONDITION(false);
    face = NULL;
    break;
  }
  return face;
}

ContactFixedSizeAllocator* ContactSearch::Get_ContactFaceAlloc(ContactFace_Type face_type, 
                                                               ContactFixedSizeAllocator *alloc) {
  ContactFixedSizeAllocator* a=NULL;
  switch (face_type){
  case ContactSearch::QUADFACEL4:
    a = &alloc[ContactSearch::ALLOC_ContactQuadFaceL4];
    break;
  case ContactSearch::QUADFACEQ8:
    a = &alloc[ContactSearch::ALLOC_ContactQuadFaceQ8];
    break;
  case ContactSearch::QUADFACEQ9:
    a = &alloc[ContactSearch::ALLOC_ContactQuadFaceQ9];
    break;
  case ContactSearch::TRIFACEL3:
    a = &alloc[ContactSearch::ALLOC_ContactTriFaceL3];
    break;
  case ContactSearch::TRIFACEQ6:
    a = &alloc[ContactSearch::ALLOC_ContactTriFaceQ6];
    break;
  case ContactSearch::SHELLQUADFACEL4:
    a = &alloc[ContactSearch::ALLOC_ContactShellQuadFaceL4];
    break;
  case ContactSearch::SHELLTRIFACEL3:
    a = &alloc[ContactSearch::ALLOC_ContactShellTriFaceL3];
    break;
  case ContactSearch::LINEFACEL2:
    a = &alloc[ContactSearch::ALLOC_ContactLineFaceL2];
    break;
  case ContactSearch::LINEFACEQ3:
    a = &alloc[ContactSearch::ALLOC_ContactLineFaceQ3];
    break;
  default:
    POSTCONDITION(false);
    a = NULL;
    break;
  }
  return a;
}

Real ContactSearch::ComputeCurvatureFromAngle(const Real angle)      {
  return (1.0-cos((PI*angle)/180.0));
}

void ContactSearch::Delete_ContactTopologyEntity(ContactTopologyEntity<Real>* entity,
                                                 ContactFixedSizeAllocator *alloc_array) {
  switch (entity->Base_Type()) {
  case (CT_NODE):
    {
      ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
      ContactFixedSizeAllocator* alloc = &alloc_array[ContactSearch::ALLOC_ContactNode];
      node->~ContactNode<Real>();
      alloc->Delete_Frag(node);
    }
    break;
  case (CT_SHELL_NODE):
    {
      ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
      ContactFixedSizeAllocator* alloc = &alloc_array[ContactSearch::ALLOC_ContactShellNode];
      (static_cast<ContactShellNode*> (node))->~ContactShellNode();
      alloc->Delete_Frag(node);
    }
    break;
  case (CT_EDGE):
    {
      ContactFixedSizeAllocator* alloc = NULL;
      ContactEdge<Real>* edge = static_cast<ContactEdge<Real>*>(entity);
      switch (edge->EdgeType()){
      case ContactSearch::LINEEDGEL2:
	alloc = &alloc_array[ContactSearch::ALLOC_ContactLineEdgeL2];
	break;
      case ContactSearch::LINEEDGEQ3:
	alloc = &alloc_array[ContactSearch::ALLOC_ContactLineEdgeQ3];
	break;
      default:
	POSTCONDITION(false);
	break;
      }
      edge->~ContactEdge<Real>();
      alloc->Delete_Frag(edge);
    }
    break;
  case (CT_FACE):
    {
      ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
      ContactFixedSizeAllocator* alloc = Get_ContactFaceAlloc(face->FaceType(), alloc_array);
      face->~ContactFace<Real>();
      alloc->Delete_Frag(face);
    }
    break;
  case (CT_ELEMENT):
    {
      ContactFixedSizeAllocator *alloc = NULL;
      ContactElement* element = static_cast<ContactElement*>(entity);
      switch (element->ElementType()){
      case ContactSearch::CARTESIANHEXELEMENTL8: 
	alloc =  &alloc_array[ContactSearch::ALLOC_ContactCartesianHexElementL8];
	break;
      case ContactSearch::HEXELEMENTL8: 
	alloc = &alloc_array[ContactSearch::ALLOC_ContactHexElementL8];
	break;
      default:
	POSTCONDITION(false);
	break;
      }
      element->~ContactElement();
      alloc->Delete_Frag(element);
    }
    break;
  default:
    POSTCONDITION(false);
    break;
  }
}
