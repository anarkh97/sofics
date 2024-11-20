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


// This file represents the FORTRAN/C interface for the ContactSearch
// object.  Only one search object is allowed for non-C++ codes.  These
// interfaces wrapper the ContactSearch object and return the error code
// in the arguement list.

#include "Search_Interface.h"
#include "ContactSearch.h"
#include "contact_assert.h"
#include "ContactVector.h"
#include "ContactRangeSearch.h"

#include <iostream>
#include <cstring>

ContactSearch* contact_search;

void FORTRAN(acme_version)( char* version )
{
  const char* vrsn = ACME_Version();
  int length = std::strlen(vrsn);
  std::memcpy( version, vrsn, length*sizeof(char) );
}


void FORTRAN(acme_versiondate)( char* versiondate )
{
  const char* vrsndate = ACME_VersionDate();
  int length = std::strlen(vrsndate);
  std::memcpy( versiondate, vrsndate, length*sizeof(char) );
}


void FORTRAN(acme_mpi_compatibility)( int& mpi_compile, int& error )
{
  error = ACME_MPI_Compatibility(mpi_compile);
}

void FORTRAN(print_search_summary)()
{
  contact_search->Print_Search_Summary();
}

void FORTRAN(build_search)( int& dimensionality, 
			    int& num_states,
			    int& num_analytic_surfaces,
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
			    int& number_element_blocks,
			    int* element_block_types,
			    int* number_elements_per_block,
			    int* element_global_ids,
			    int* element_connectivity,
			    int& num_comm_partners,
			    int* comm_proc_ids,
			    int* number_nodes_to_partner,
			    int* comm_node,
			    int& icomm,
			    int& ierror)
{
  int i;

  // I'm going to do an explicit cast of int to enum.  This works only if
  // the compiler doesn't use a smaller size for the enum (current none
  // do).  The following preconditions will fail if we ever get a 
  // compiler that does this.
  PRECONDITION( sizeof(int) == sizeof(ContactSearch::ContactNode_Type) );
  PRECONDITION( sizeof(int) == sizeof(ContactSearch::ContactFace_Type) );
  PRECONDITION( sizeof(int) == sizeof(ContactSearch::ContactEdge_Type) );
  PRECONDITION( sizeof(int) == sizeof(ContactSearch::AnalyticSurface_Type) );

#ifndef CONTACT_NO_MPI
  MPI_Comm comm = MPI_COMM_F2C((MPI_INT_TYPE)icomm);
#else
  int comm = icomm;
#endif
  int* face_block_ids = new int[number_face_blocks];
  for( i=0 ; i<number_face_blocks; ++i)
    face_block_ids[i] = i+1;
  ContactSearch::ContactErrorCode error;
  contact_search = new 
    ContactSearch( dimensionality, num_states, num_analytic_surfaces,
		   num_node_blocks,
		   (ContactSearch::ContactNode_Type*) node_block_types,
		   num_nodes_per_block, 
                   node_exodus_ids, 
                   node_global_ids, 
		   coords,
                   number_face_blocks,
		   (ContactSearch::ContactFace_Type*) face_block_types,
		   number_faces_per_block, 
                   face_global_ids, 
                   face_connectivity,
		   face_lofting_factors,
		   number_element_blocks, 
		   (ContactSearch::ContactElement_Type*) element_block_types,
		   number_elements_per_block, 
                   element_global_ids,
		   element_connectivity,
		   num_comm_partners, comm_proc_ids, 
		   number_nodes_to_partner, 
		   comm_node, comm, error);
  delete [] face_block_ids;

  ierror = ContactSearch::NO_ERROR;
}

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

			     int& ierror)
{
  ContactSearch::ContactErrorCode error;
  contact_search->UpdateSearch( num_node_deaths_per_block, 
			        node_deaths_global_ids,
                                num_face_deaths_per_block, 
			        face_deaths_global_ids,
                                num_element_deaths_per_block, 
			        element_deaths_global_ids,
                                
                                num_node_births_per_block, 
			        node_births_exodus_ids,
			        node_births_global_ids,
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
                                
			        error);
  ierror = ContactSearch::NO_ERROR;
}

void FORTRAN(add_table)( int& ID, int& Num_Points, 
			 Real* abscissa, Real* ordinate, int& error )
{
  error = contact_search->Add_Table( ID, Num_Points, abscissa, ordinate );
}

void FORTRAN(check_search_data_size)( int& size_data_per_pair,
				      int& number_of_entity_keys,
				      int& error )
{
  error = contact_search->Check_Search_Data_Size( size_data_per_pair,
						  number_of_entity_keys );
}

void FORTRAN(set_search_data)( Real* search_data )
{
  contact_search->Set_Search_Data( search_data );
}

void FORTRAN(add_analytic_surface)( int& surf_type, Real* data , int& error )
{
  error = contact_search->Add_Analytic_Surface( 
	  (ContactSearch::AnalyticSurface_Type) surf_type, data );
}

void FORTRAN(set_analytic_surface_configuration)( int& id, 
						  Real* data , 
						  int& error )
{
  error = contact_search->Set_Analytic_Surface_Configuration( id, data );
}

void FORTRAN(set_node_block_configuration)( int& config_type, int& node_blk_id,
					    Real* positions, int& error )
{
  ContactSearch::ContactNode_Configuration config = 
    (ContactSearch::ContactNode_Configuration) config_type;
  error = contact_search->Set_Node_Block_Configuration(config,node_blk_id,
						       positions );
}

void FORTRAN(set_node_block_remaining_gap)( int& node_block_id, 
				            Real* gap, int& error )
{
  error = contact_search->Set_Node_Block_Remaining_Gap( node_block_id,
							gap );
}

void FORTRAN(set_node_block_ghosting_gap)( int& node_block_id, 
				           Real* gap, int& error )
{
  error = contact_search->Set_Node_Block_Ghosting_Gap( node_block_id,
						       gap );
}

void FORTRAN(set_node_block_kin_cons)( int& node_block_id, int* cons_per_node, 
				       Real* cons_vec, int& error )
{
  error = contact_search->Set_Node_Block_Kinematic_Constraints( node_block_id,
								cons_per_node,
								cons_vec );
}

void FORTRAN(set_node_block_attributes)( int& attr_type, int& node_blk_id,
					 Real* attributes, int& error )
{
  ContactSearch::Node_Block_Attribute attr = 
    (ContactSearch::Node_Block_Attribute) attr_type;
  error = contact_search->Set_Node_Block_Attributes( attr, node_blk_id,
						     attributes );
}

void FORTRAN(set_face_block_attributes)( int& attr_type, int& face_blk_id,
					 Real* attributes, int& error )
{
  ContactSearch::Face_Block_Attribute attr = 
    (ContactSearch::Face_Block_Attribute) attr_type;
  error = contact_search->Set_Face_Block_Attributes( attr, face_blk_id,
						     attributes );
}

void FORTRAN(set_search_option)( int& option, int& status, Real* data,
				 int& error )
{
 error = contact_search->Set_Search_Option((ContactSearch::Search_Option) 
					   option,
					   (ContactSearch::Search_Option_Status)
					   status,
					   data );
}

void FORTRAN(static_search_1_configuration)(int& error)
{
  error = contact_search->Static_Search_1_Configuration();
}

void FORTRAN(static_search_2_configuration)(int& error)
{
  error = contact_search->Static_Search_2_Configuration();
}

void FORTRAN(dynamic_search_2_config)(Real& dt_old,
			              Real& dt,
				      int& error )
{
  error = contact_search->Dynamic_Search_2_Configuration( dt_old,
					                  dt );
}

void FORTRAN(dynamic_tracking_search_2_configuration)(int& error)
{
  error = contact_search->Dynamic_Tracking_Search_2_Configuration();
}

void FORTRAN(delete_all_interactions)(int& error)
{
  error = contact_search->Delete_All_Interactions();
}

void FORTRANNO_(display)()
{
  contact_search->Display();
}

void FORTRAN(size_nodenode_interactions)(int& num_interactions,
                                         int& data_size)
{
  contact_search->Size_NodeNode_Interactions(num_interactions,data_size);
}

void FORTRAN(get_nodenode_interactions)(int* slave_node_block_ids,
					int* slave_node_indexes_in_block,
					int* master_node_block_ids,
					int* master_node_indexes_in_block,
					int* master_node_proc,
					Real* interaction_data)
{
  contact_search->Get_NodeNode_Interactions(slave_node_block_ids,
					    slave_node_indexes_in_block,
					    master_node_block_ids,
					    master_node_indexes_in_block,
					    master_node_proc,interaction_data);
}

void FORTRAN(size_nodeface_interactions)(int& num_interactions,
                                         int& data_size)
{
  contact_search->Size_NodeFace_Interactions(num_interactions,data_size);
}

void FORTRAN(get_nodeface_interactions)(int* node_block_ids,
					int* node_indexes_in_block,
					int* node_entity_keys,
					int* face_block_ids,
					int* face_indexes_in_block,
					int* face_proc,
					Real* interaction_data)
{
  contact_search->Get_NodeFace_Interactions(node_block_ids,
					    node_indexes_in_block,
					    node_entity_keys,
					    face_block_ids,
					    face_indexes_in_block,
					    face_proc,interaction_data);
}

void FORTRAN(size_nodesurface_interactions)(int& num_interactions,
					    int& data_size)
{
  contact_search->Size_NodeSurface_Interactions(num_interactions,data_size);
}

void FORTRAN(get_nodesurface_interactions)(int* node_block_ids,
					   int* node_indexes_in_block,
					   int* surface_ids,
					   Real* interaction_data)
{
  contact_search->Get_NodeSurface_Interactions(node_block_ids,
					       node_indexes_in_block,
					       surface_ids, interaction_data);
}

void FORTRAN(size_faceface_interactions)(int& num_interactions,
                                         int& data_size)
{
  contact_search->Size_FaceFace_Interactions(num_interactions,data_size);
}

void FORTRAN(get_faceface_interactions)(int* slave_face_block_ids,
					int* slave_face_indexes_in_block,
					int* slave_face_proc,
					int* master_face_block_ids,
					int* master_face_indexes_in_block,
					int* master_face_proc,
                                        int* interaction_index,
					Real* interaction_data)
{
  contact_search->Get_FaceFace_Interactions(slave_face_block_ids,
					    slave_face_indexes_in_block,
						slave_face_proc,
					    master_face_block_ids,
					    master_face_indexes_in_block,
					    master_face_proc,
                                            interaction_index,
                                            interaction_data);
}

void FORTRAN(size_facecoverage_interactions)(int& num_interactions,
                                             int& data_size)
{
  contact_search->Size_FaceCoverage_Interactions(num_interactions,data_size);
}

void FORTRAN(get_facecoverage_interactions)(int* face_block_ids,
					    int* face_indexes_in_block,
                                        int* interaction_index,
					    Real* interaction_data)
{
  contact_search->Get_FaceCoverage_Interactions(face_block_ids,
					        face_indexes_in_block,
                                                interaction_index,
					        interaction_data);
}

void FORTRAN(size_elementelement_interactions)(int& num_interactions,
                                               int& data_size)
{
  contact_search->Size_ElementElement_Interactions(num_interactions,data_size);
}

void FORTRAN(get_elementelement_interactions)(int* slave_element_block_ids,
				   	      int* slave_element_indexes_in_block,
					      int* master_element_block_ids,
					      int* master_element_indexes_in_block,
					      int* master_element_proc,
					      Real* interaction_data)
{
  contact_search->Get_ElementElement_Interactions(slave_element_block_ids,
					          slave_element_indexes_in_block,
					          master_element_block_ids,
					          master_element_indexes_in_block,
					          master_element_proc,
                                                  interaction_data);
}

void FORTRAN(cleanup_search)()
{
  delete contact_search;
}

void FORTRAN(number_of_search_errors)( int& num_errors )
{
  num_errors = contact_search->Number_of_Errors();
}

void FORTRAN(get_search_error_message)( int& i, char* message )
{
  const char* msg = contact_search->Error_Message( i );
  std::memcpy( message, msg, 80*sizeof(char) );
}

void FORTRAN(contact_exodus_output)( int& exodus_id, Real& time, int& error )
{
  error = contact_search->Exodus_Output( exodus_id, time );
}

void FORTRAN(build_search_restart)( Real* restart_buffer, 
			            int* node_global_ids,
			            int* face_global_ids,
			            int* element_global_ids,
                                    int& icomm, 
				    int& error )
{
#ifndef CONTACT_NO_MPI
  MPI_Comm comm = MPI_COMM_F2C((MPI_INT_TYPE)icomm);
#else
  int comm = icomm;
#endif
  ContactSearch::ContactErrorCode ec;
  contact_search =  new ContactSearch( restart_buffer, 
                                       node_global_ids,
                                       face_global_ids, 
                                       element_global_ids, 
                                       comm, ec );
  error = ec;
}

void FORTRAN(search_restart_size)( int& size )
{
  size = contact_search->Restart_Size();
}

void FORTRAN(search_extract_restart_data)( Real* buffer, int& error )
{
  error = contact_search->Extract_Restart_Data( buffer );
}

void FORTRAN(search_num_general_rsvars)( int& num )
{
  num = contact_search->Number_General_Restart_Variables();
}

void FORTRAN(search_num_node_rsvars)( int& num )
{
  num = contact_search->Number_Nodal_Restart_Variables();
}

void FORTRAN(search_num_edge_rsvars)( int& num )
{
  num = contact_search->Number_Edge_Restart_Variables();
}

void FORTRAN(search_num_face_rsvars)( int& num )
{
  num = contact_search->Number_Face_Restart_Variables();
}

void FORTRAN(search_num_elem_rsvars)( int& num )
{
  num = contact_search->Number_Element_Restart_Variables();
}

void FORTRAN(search_extract_general_rsvars)( Real* data, int& ierror )
{
  ierror = contact_search->Extract_General_Restart_Variable( data );
}

void FORTRAN(search_extract_node_rsvars)( int& n, Real* data, int *node_ids, int& ierror )
{
  ierror = contact_search->Extract_Nodal_Restart_Variable( n, data, node_ids );
}

void FORTRAN(search_extract_edge_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_search->Extract_Edge_Restart_Variable( n, data );
}

void FORTRAN(search_extract_face_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_search->Extract_Face_Restart_Variable( n, data );
}

void FORTRAN(search_extract_elem_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_search->Extract_Element_Restart_Variable( n, data );
}

void FORTRAN(search_implant_general_rsvars)( Real* data, int& ierror )
{
  ierror = contact_search->Implant_General_Restart_Variable( data );
}

void FORTRAN(search_implant_node_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_search->Implant_Nodal_Restart_Variable( n, data );
}

void FORTRAN(search_implant_edge_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_search->Implant_Edge_Restart_Variable( n, data );
}

void FORTRAN(search_implant_face_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_search->Implant_Face_Restart_Variable( n, data );
}

void FORTRAN(search_implant_elem_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_search->Implant_Element_Restart_Variable( n, data );
}

void FORTRAN(search_complete_restart)(int& ierror)
{
  ierror = contact_search->Complete_Restart();
}

void FORTRAN(add_debug_node)(int& node_id, int& ierror )
{ 
  ierror = contact_search->Add_Debug_Node( node_id );
}
 

void FORTRAN(overlapping_sphere_ghost)(int& num_sphere,
                                       Real *sphere_data,
                                       int& icomm,
                                       int* ghost_indexes,
				       int* ghost_procs) {
  static ACME::Int_Vector ghost_index_vec;
  static ACME::Int_Vector ghost_procs_vec;

#ifndef CONTACT_NO_MPI
  MPI_Comm comm = MPI_COMM_F2C((MPI_INT_TYPE)icomm);
#else
  int comm = icomm;
#endif
  ACME::Overlapping_Sphere_Ghost(num_sphere, sphere_data, comm,
                                 ghost_index_vec, ghost_procs_vec);

  ghost_indexes = ghost_index_vec.get_buffer();
  ghost_procs = ghost_procs_vec.get_buffer();


}


void FORTRAN(overlapping_sphere_search)(int& num_sphere1,
                                        Real *sphere_data1,
                                        int& num_sphere2,
                                        Real *sphere_data2,
                                        int *interaction_list,
                                        int *first_interaction,
                                        int *last_interaction) {
  static ACME::Int_Vector interaction_list_vec;
  static ACME::Int_Vector first_interaction_vec;
  static ACME::Int_Vector last_interaction_vec;

  ACME::Overlapping_Sphere_Search(num_sphere1,
                                  sphere_data1,
                                  num_sphere2,
                                  sphere_data2,
                                  interaction_list_vec,
                                  first_interaction_vec,
                                  last_interaction_vec);

  interaction_list = interaction_list_vec.get_buffer();
  first_interaction = first_interaction_vec.get_buffer();
  last_interaction = last_interaction_vec.get_buffer();
}


void FORTRAN(sphere_point_search)(int  &num_sphere,
                                  Real *sphere_data,
                                  int  &num_point,
                                  Real *point_data,
                                  int *interaction_list,
                                  int *first_interaction,
                                  int *last_interaction) {
  static ACME::Int_Vector interaction_list_vec;
  static ACME::Int_Vector first_interaction_vec;
  static ACME::Int_Vector last_interaction_vec;
  ACME::Sphere_Point_Search(num_sphere,
                            sphere_data,
                            num_point,
                            point_data,
                            interaction_list_vec,
                            first_interaction_vec,
	                    last_interaction_vec);
  interaction_list = interaction_list_vec.get_buffer();
  first_interaction = first_interaction_vec.get_buffer();
  last_interaction = last_interaction_vec.get_buffer();
}

