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


#include "ContactSearch.h"
#include "Contact_Defines.h"

#ifndef CONTACT_NO_MPI
#include "mpi.h"
#endif

extern "C" {
  
  void FORTRAN(acme_version)( char* version );
  
  void FORTRAN(acme_versiondate)( char* versiondate );
  
  void FORTRAN(acme_mpi_compatibility)( int&, int& );

  void FORTRAN(build_search)( int& dimensionality,
			      int& num_states,
			      int& num_entity_keys,
			      int& num_node_blocks,
			      int* node_block_types,
			      int* num_nodes_per_block, 
			      int* node_exodus_ids,
			      int* node_global_ids,
			      double* coords,
			      int& number_face_blocks,
			      int* face_block_types,
			      int* number_faces_per_block,
			      int* face_global_ids, 
			      int* face_connectivity,
			      double* face_lofting_factors,
			      int& num_element_blocks,
			      int* element_block_types,
			      int* number_elements_per_block,
			      int* element_ids,
			      int* element_connectivity,
			      int& num_comm_partners,
			      int* comm_proc_id,
			      int* number_nodes_to_partner,
			      int* comm_node,
			      int& MPI_Comm,
			      int& error );
                     
  void FORTRAN(update_search)( int* num_node_deaths_per_block, 
			       int* node_deaths_global_ids,
                               int* num_face_deaths_per_block, 
			       int* face_deaths_global_ids,
                               int* num_element_deaths_per_block, 
			       int* element_deaths_global_ids,
                               
                               int* num_node_births_per_block, 
			       int* node_births_exodus_ids,
			       int* node_births_global_ids,
			       int* number_face_births_per_block, 
			       int* face_births_global_ids,
			       int* face_births_connectivity,
			       int* number_element_births_per_block, 
			       int* element_births_global_ids,
			       int* element_births_connectivity,
                               
			       int& num_node_exports,
                               int* node_export_id_list,
                               int* node_export_pid,
			       int& num_face_exports,
                               int* face_export_id_list,
                               int* face_export_pid,
			       int& num_element_exports,
                               int* element_export_id_list,
                               int* element_export_pid,
                               
                               int& num_nodes,
                               int* node_host_ids,
                               int& num_faces,
                               int* face_host_ids,
                               int& num_elements,
                               int* element_host_ids,
         
			       int& num_comm_partners,
			       int* comm_proc_id,
			       int* number_nodes_to_partner,
			       int* comm_node,
                               
			       int& error );

  void FORTRAN(add_analytic_surface)( int& surf_type, 
				      Real* data, 
				      int& error );
  
  void FORTRAN(set_analytic_surface_configuration)( int& id, 
				      Real* data, 
				      int& error );

  void FORTRAN(add_table)( int& ID, int& Num_Points, 
			   Real* abscissa, Real* ordinate, int& error );
  
  void FORTRAN(check_search_data_size)( int& size_data_per_pair,
					int& number_of_entity_keys,
					int& error );

  void FORTRAN(set_search_data)( Real* data );

  void FORTRAN(set_node_block_configuration)( int& config_type, 
					      int& node_blk_id,
					      Real* positions, 
					      int& error );
       
  void FORTRAN(set_node_block_remaining_gap)( int& node_blk_id, 
					      Real* gap, int& error );
       
  void FORTRAN(set_node_block_ghosting_gap)( int& node_blk_id, 
					     Real* gap, int& error );
       
  void FORTRAN(set_node_block_kin_cons)( int& node_blk_id, int* cons_per_node,
					 Real* cons_vec, int& error );

  void FORTRAN(set_node_block_attributes)( int& config_type, int& node_blk_id,
					   Real* attributes, int& error );

  void FORTRAN(set_face_block_attributes)( int& attr, int& face_blk_id,
					   Real* attributes, int& error );

  void FORTRAN(set_search_option)( int& option, int& status, Real* data,
				   int& ierror );
  
  void FORTRAN(static_search_1_configuration)( int& error );
  
  void FORTRAN(static_search_2_configuration)( int& error );
  
  void FORTRAN(dynamic_search_2_config)( Real& dt_old,
				         Real& dt,
					 int& error );
  
  void FORTRAN(dynamic_tracking_search_2_configuration)( int& error );

  void FORTRANNO_(display)();

  void FORTRAN(size_nodenode_interactions)(int& num_interactions,
					   int& data_size);

  void FORTRAN(get_nodenode_interactions)(int* slave_node_block_ids,
					  int* slave_node_indexes_in_block,
					  int* master_node_block_ids,
					  int* master_node_indexes_in_block,
					  int* master_node_proc,
					  Real* interaction_data);

  void FORTRAN(size_nodeface_interactions)(int& num_interactions,
					   int& data_size);

  void FORTRAN(get_nodeface_interactions)(int* node_block_ids,
					  int* node_indexes_in_block,
					  int* node_entity_keys,
					  int* face_block_ids,
					  int* face_indexes_in_block,
					  int* face_proc,
					  Real* interaction_data);

  void FORTRAN(size_nodesurface_interactions)(int& num_interactions,
					      int& data_size );

  void FORTRAN(get_nodesurface_interactions)( int* node_block_ids,
					      int* node_indexes_in_block,
					      int* surface_id,
					      Real* interaction_data );

  void FORTRAN(size_faceface_interactions)(int& num_interactions,
					   int& data_size);

  void FORTRAN(get_faceface_interactions)(int* slave_face_block_ids,
					  int* slave_face_indexes_in_block,
					  int* master_face_block_ids,
					  int* master_face_indexes_in_block,
					  int* master_face_proc,
					  int* interaction_index,
					  Real* interaction_data);

  void FORTRAN(size_facecoverage_interactions)(int& num_interactions,
					       int& data_size);

  void FORTRAN(get_facecoverage_interactions)(int* face_block_ids,
					      int* face_indexes_in_block,
					      int* interaction_index,
					      Real* interaction_data);

  void FORTRAN(size_elementelement_interactions)(int& num_interactions,
					         int& data_size);

  void FORTRAN(get_elementelement_interactions)(int* slave_element_block_ids,
					        int* slave_element_indexes_in_block,
					         int* master_element_block_ids,
					         int* master_element_indexes_in_block,
					         int* master_element_proc,
					         Real* interaction_data);

  void FORTRAN(number_of_search_errors)( int& num_errors );
  
  void FORTRAN(get_search_error_message)( int& i, char* message );

  void FORTRAN(cleanup_search)();

  void FORTRAN(contact_exodus_output)( int& exodus_id, Real& time, 
				       int& ierror );

  void FORTRAN(build_search_restart)( Real* restart_buffer, 
			              int* node_global_ids,
			              int* face_global_ids,
			              int* element_global_ids,
                                      int& MPI_Comm,
				      int& ierror );

  void FORTRAN(search_restart_size)( int& );

  void FORTRAN(search_extract_restart_data)( Real*, int& );

  void FORTRAN(search_num_node_rsvars)( int& num );

  void FORTRAN(search_num_edge_rsvars)( int& num );

  void FORTRAN(search_num_face_rsvars)( int& num );

  void FORTRAN(search_num_elem_rsvars)( int& num );

  void FORTRAN(search_num_general_rsvars)( int& num );

  void FORTRAN(search_extract_node_rsvars)( int& n, Real* data, int *node_ids, int& ierror );

  void FORTRAN(search_extract_edge_rsvars)( int& n, Real* data, int& ierror );

  void FORTRAN(search_extract_face_rsvars)( int& n, Real* data, int& ierror );

  void FORTRAN(search_extract_elem_rsvars)( int& n, Real* data, int& ierror );

  void FORTRAN(search_extract_general_rsvars)( Real* data, int& ierror );

  void FORTRAN(search_implant_node_rsvars)( int& n, Real* data, int& ierror );

  void FORTRAN(search_implant_edge_rsvars)( int& n, Real* data, int& ierror );

  void FORTRAN(search_implant_face_rsvars)( int& n, Real* data, int& ierror );

  void FORTRAN(search_implant_elem_rsvars)( int& n, Real* data, int& ierror );
 
  void FORTRAN(search_implant_general_rsvars)( Real* data, int& ierror );
 
  void FORTRAN(search_complete_restart)(int& ierror );

  void FORTRAN(add_debug_node)( int& node_id, int& ierror );

  void FORTRAN(overlapping_sphere_ghost)(int& num_sphere,
                                         Real *sphere_data,
                                         int& mpi_comm,
                                         int* ghost_indexes,
                                         int* ghost_procs);

  void FORTRAN(overlapping_sphere_search)(int& num_sphere1,
                                          Real *sphere_data1,
                                          int& num_sphere2,
                                          Real *sphere_data2,
                                          int *interaction_list,
                                          int *first_interaction,
                                          int *last_interaction);

  void FORTRAN(sphere_point_search)(int  &num_sphere,
                                    Real *sphere_data,
                                    int  &num_point,
                                    Real *point_data,
                                    int *interaction_list,
                                    int *first_interaction,
                                    int *last_interaction);
                                    
  void FORTRAN(print_search_summary)();

}

