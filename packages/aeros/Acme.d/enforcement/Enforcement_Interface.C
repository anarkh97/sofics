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


#include <cstring>
#include "Enforcement_Interface.h"
#include "contact_assert.h"
#include "ContactTDUserSubTypes.h"

//contact_td_enforcement  = NULL;
//contact_td_penalty      = NULL;
//contact_gap_removal     = NULL;
//contact_tied_kinematics = NULL;
//contact_volume_transfer = NULL;
//contact_mpc_eqns        = NULL;

extern ContactSearch* contact_search;

void FORTRAN(number_of_enforcement_errors)( int& num_errors )
{
  if (contact_td_enforcement!=NULL) 
    num_errors = contact_td_enforcement->Number_of_Errors();
  if (contact_td_penalty!=NULL) 
    num_errors = contact_td_penalty->Number_of_Errors();
  if (contact_gap_removal!=NULL) 
    num_errors = contact_gap_removal->Number_of_Errors();
  if (contact_tied_kinematics!=NULL) 
    num_errors = contact_tied_kinematics->Number_of_Errors();
  if (contact_volume_transfer!=NULL) 
    num_errors = contact_volume_transfer->Number_of_Errors();
  if (contact_mpc_eqns!=NULL) 
    num_errors = contact_mpc_eqns->Number_of_Errors();
}

void FORTRAN(get_enforcement_error_message)( int& i, char* message )
{
  if (contact_td_enforcement!=NULL) 
    std::memcpy( message, 
                 contact_td_enforcement->Error_Message( i ), 
                 80*sizeof(char) );
  if (contact_td_penalty!=NULL) 
    std::memcpy( message, 
                 contact_td_penalty->Error_Message( i ), 
                 80*sizeof(char) );
  if (contact_gap_removal!=NULL) 
    std::memcpy( message, 
                 contact_gap_removal->Error_Message( i ), 
                 80*sizeof(char) );
  if (contact_tied_kinematics!=NULL) 
    std::memcpy( message, 
                 contact_tied_kinematics->Error_Message( i ), 
                 80*sizeof(char) );
  if (contact_volume_transfer!=NULL) 
    std::memcpy( message, 
                 contact_volume_transfer->Error_Message( i ), 
                 80*sizeof(char) );
  if (contact_mpc_eqns!=NULL) 
    std::memcpy( message, 
                 contact_mpc_eqns->Error_Message( i ), 
                 80*sizeof(char) );
}

void FORTRAN(build_td_enforcement)( double* Enforcement_Data,
                                    int* get_cvars,
                                    int* plot_force,
				    int* error )
{
  ContactSearch::ContactErrorCode ec;
  contact_td_enforcement = new ContactTDEnforcement(Enforcement_Data,
						    contact_search,
						    ec,
                                                    *get_cvars,
                                                    *plot_force );
  *error = ec;
}


void FORTRAN(build_td_enf_restart)( double* restart_data,
                                    int* get_cvars,
                                    int* plot_force,
				    int* error )
{
  ContactSearch::ContactErrorCode ec;
  contact_td_enforcement = new ContactTDEnforcement(contact_search,
						    restart_data,
						    ec,
                                                    *get_cvars,
                                                    *plot_force );
  *error = ec;
}

void FORTRAN(td_enf_restart_size)( int& size )
{
  size = contact_td_enforcement->Restart_Size();
}

void FORTRAN( td_enf_extract_restart )( Real* buffer, int& error )
{
  error = contact_td_enforcement->Extract_Restart_Data( buffer );
}

void FORTRAN( set_td_iterations )( int* num_iter,
				   int* error )
{
  *error = contact_td_enforcement->Set_Number_of_Iterations( *num_iter );
}

void FORTRAN( set_td_convg_tol )( double* tol,
				                   int* error )
{
  *error = contact_td_enforcement->Set_Convergence_Tolerance( *tol );
}


void FORTRAN( set_td_enf_symm_nodes )( int& node_1, int& node_2, int& error )
{
  error = contact_td_enforcement->Enforce_Symmetry_on_Nodes( node_1, node_2 );
}

void FORTRAN( compute_td_contact_force)( double& dt_old,
					 double& dt,
					 double* mass,
					 double* density,
					 double* wavespeed,
					 double* forces,
					 int* error )
{
  *error = contact_td_enforcement->Compute_Contact_Force( dt_old,
							  dt,
							  mass,
							  density,
							  wavespeed,
							  forces );
}

void FORTRAN( get_td_plot_variable )( int& plot_var, double* data, int* error )
{
  ContactTDEnforcement::Contact_TDEnf_Plot_Vars Var =
    (ContactTDEnforcement::Contact_TDEnf_Plot_Vars) plot_var;
  *error = contact_td_enforcement->Get_Plot_Variable( Var, data );
}

void FORTRAN( get_td_global_variable )( int& plot_var, double data, int* error )
{
  ContactTDEnforcement::Contact_TDEnf_Global_Vars Var =
    (ContactTDEnforcement::Contact_TDEnf_Global_Vars) plot_var;
  *error = contact_td_enforcement->Get_Global_Variable( Var, data );
}

void FORTRAN(cleanup_td_enforcement)()
{
  delete contact_td_enforcement;
}

void FORTRAN(build_td_penalty)( double* Enforcement_Data,
                                   int* get_cvars,
                                   int* get_plot_force,
				   int* error )
{
  ContactSearch::ContactErrorCode ec;
  contact_td_penalty = new ContactTDEnfPenalty(Enforcement_Data,
						        contact_search,
						        ec,
                                                        *get_cvars,
                                                        *get_plot_force );
  *error = ec;
}

void FORTRAN(build_td_penalty_restart)( double* restart_data,
                                           int* get_cvars,
                                           int* get_plot_force,
				           int* error )
{
  ContactSearch::ContactErrorCode ec;
  contact_td_penalty = new ContactTDEnfPenalty(contact_search,
						    restart_data,
						    ec,
                                                    *get_cvars,
                                                    *get_plot_force );
  *error = ec;
}

void FORTRAN(td_penalty_restart_size)( int& size )
{
  size = contact_td_penalty->Restart_Size();
}

void FORTRAN( td_penalty_extract_restart )( Real* buffer, int& error )
{
  error = contact_td_penalty->Extract_Restart_Data( buffer );
}


void FORTRAN( compute_td_penalty_contact_force)( double& dt_old,
					 double& dt,
					 double* mass,
					 double* density,
					 double* wavespeed,
					 double* forces,
					 int* error )
{
  *error = contact_td_penalty->Compute_Contact_Force( dt_old,
							  dt,
							  mass,
							  density,
							  wavespeed,
							  forces );
}

void FORTRAN( set_td_penalty_iterations )( int* num_iter,
				   int* error )
{
  *error = contact_td_penalty->Set_Number_of_Iterations( *num_iter );
}

void FORTRAN( set_td_penalty_convg_tol )( double* tol,
				                   int* error )
{
  *error = contact_td_penalty->Set_Convergence_Tolerance( *tol );
}

void FORTRAN( set_td_penalty_pen_scale )( double* scale,
				                   int* error )
{
  *error = contact_td_penalty->Set_Penalty_Scale( *scale );
}

void FORTRAN(cleanup_td_penalty)()
{
  delete contact_td_penalty;
}

void FORTRAN( compute_vol_tran )( int& num_node_vars,
		                  int& num_elem_vars,
                                  double* donar_node_vars,
                                  double* donar_elem_vars,
			          double* receiver_node_vars,
				  double* receiver_elem_vars,
				  double* volume_fraction,
                                  int* error )
{
  ContactSearch::ContactErrorCode ec;
  ec = contact_volume_transfer->Compute_Volume_Transfer( num_node_vars,
							 num_elem_vars,
							 donar_node_vars,
							 donar_elem_vars,
							 receiver_node_vars,
							 receiver_elem_vars,
							 volume_fraction );
  *error = ec;
}

void FORTRAN( get_num_mpceqns )( int& nmpceq, int* error )
{
  ContactSearch::ContactErrorCode ec;
  ec = contact_mpc_eqns->Number_of_MPC_Equations(nmpceq);

  *error = ec;
}

void FORTRAN(cleanup_mpc_eqns)()
{
  delete contact_mpc_eqns;
}

void FORTRAN( compute_mpc_eqns )( int* error )
{
  ContactSearch::ContactErrorCode ec;
  ec = contact_mpc_eqns->Compute_MPCs( );
  *error = ec;
}

void FORTRAN( get_mpc_eqns )( int& num_mpcs,
			      int* snode_pid,
			      int* snode_lid,
			      int* mface_pid,
			      int* mface_lid,
			      int* nface_nodes,
			      int* fnode_pid,
			      int* fnode_lid,
			      Real* fnode_coefs,
			      int* error )
{
  ContactSearch::ContactErrorCode ec;

  int number_of_eqns = num_mpcs;

  ec = contact_mpc_eqns->Get_MPC_Equations( number_of_eqns,
					    snode_pid,
					    snode_lid,
					    mface_pid,
					    mface_lid,
					    nface_nodes,
					    fnode_pid,
					    fnode_lid,
					    fnode_coefs );

  *error = ec;
}


void FORTRAN(build_gap_removal)( double* Enforcement_Data,
				 int* error )
{
  ContactSearch::ContactErrorCode ec;
  contact_gap_removal = new ContactGapRemoval(Enforcement_Data,
					      contact_search,
					      ec );
  *error = ec;
}


void FORTRAN(build_gap_removal_restart)( double* restart_data,
				    int* error )
{
  ContactSearch::ContactErrorCode ec;
  contact_gap_removal = new ContactGapRemoval(contact_search,
					      restart_data,
					      ec );
  *error = ec;
}


void FORTRAN(gap_removal_restart_size)( int& size )
{
  size = contact_gap_removal->Restart_Size();
}

void FORTRAN( gap_removal_extract_restart )( Real* buffer, int& error )
{
  error = contact_gap_removal->Extract_Restart_Data( buffer );
}

void FORTRAN( compute_gap_removal)( int& max_iterations,
				    double& trivial_gap,
				    double* displc,
				    int* error )
{
  *error = contact_gap_removal->Compute_Gap_Removal( max_iterations,
						     trivial_gap, 
						     displc );
}

 

void FORTRAN(cleanup_gap_removal)()
{
  delete contact_gap_removal;
}


void FORTRAN(build_tied_kin)( double* Enforcement_Data,
				 int* error )
{
  ContactSearch::ContactErrorCode ec;
  contact_tied_kinematics = new ContactTiedKinematics(Enforcement_Data,
						      contact_search,
						      ec );
  *error = ec;
}


void FORTRAN(build_vol_tran)( double* Enforcement_Data,
			      int* error )
{
  ContactSearch::ContactErrorCode ec;
  contact_volume_transfer = new ContactVolumeTransfer(Enforcement_Data,
						      contact_search,
						      ec );
  *error = ec;
}

void FORTRAN(build_mpc_eqns)( double* Enforcement_Data,
			      int* error )
{
  ContactSearch::ContactErrorCode ec;
  contact_mpc_eqns = new ContactMPCs( Enforcement_Data,
				      contact_search,
				      ec );
  *error = ec;
}

void FORTRAN(build_tied_kin_restart)( double* restart_data,
				      int* error )
{
  ContactSearch::ContactErrorCode ec;
  contact_tied_kinematics = new ContactTiedKinematics(contact_search,
						      restart_data,
						      ec );
  *error = ec;
}

void FORTRAN(build_vol_tran_restart)( double* restart_data,
				      int* error )
{
  ContactSearch::ContactErrorCode ec;
  contact_volume_transfer = new ContactVolumeTransfer(contact_search,
						      restart_data,
						      ec );
  *error = ec;
}

void FORTRAN(build_mpc_eqns_restart)( double* restart_data,
				      int* error )
{
  ContactSearch::ContactErrorCode ec;
  contact_mpc_eqns = new ContactMPCs(contact_search,
				     restart_data,
				     ec );
  *error = ec;
}



void FORTRAN(tied_kin_restart_size)( int& size )
{
  size = contact_tied_kinematics->Restart_Size();
}

void FORTRAN(vol_tran_restart_size)( int& size )
{
  size = contact_volume_transfer->Restart_Size();
}

void FORTRAN(mpc_eqns_restart_size)( int& size )
{
  size = contact_mpc_eqns->Restart_Size();
}

void FORTRAN( tied_kin_extract_restart )( Real* buffer, int& error )
{
  error = contact_tied_kinematics->Extract_Restart_Data( buffer );
}

void FORTRAN( vol_tran_extract_restart )( Real* buffer, int& error )
{
  error = contact_volume_transfer->Extract_Restart_Data( buffer );
}

void FORTRAN( mpc_eqns_extract_restart )( Real* buffer, int& error )
{
  error = contact_mpc_eqns->Extract_Restart_Data( buffer );
}

void FORTRAN( compute_tied_kin_position)( double* pos,
					  int* error )
{
  *error = contact_tied_kinematics->Compute_Position( pos );
}

void FORTRAN(cleanup_tied_kinematics)()
{
  delete contact_tied_kinematics;
}

void FORTRAN(cleanup_vol_tran)()
{
  delete contact_volume_transfer;
}

void FORTRAN(number_of_td_errors)( int& num_errors )
{
  num_errors = contact_td_enforcement->Number_of_Errors();
}

void FORTRAN(get_td_error_message)( int& i, char* message )
{
  const char* msg = contact_td_enforcement->Error_Message( i );
  std::memcpy( message, msg, 80*sizeof(char) );
}
    
void FORTRAN(number_of_td_penalty_errors)( int& num_errors )
{
  num_errors = contact_td_penalty->Number_of_Errors();
}

void FORTRAN(get_td_penalty_error_message)( int& i, char* message )
{
  const char* msg = contact_td_penalty->Error_Message( i );
  std::memcpy( message, msg, 80*sizeof(char) );
}
    
void FORTRAN(number_of_gap_errors)( int& num_errors )
{
  num_errors = contact_gap_removal->Number_of_Errors();
}

void FORTRAN(get_gap_error_message)( int& i, char* message )
{
  const char* msg = contact_gap_removal->Error_Message( i );
  std::memcpy( message, msg, 80*sizeof(char) );
}

void FORTRAN(number_of_tied_errors)( int& num_errors )
{
  num_errors = contact_tied_kinematics->Number_of_Errors();
}

void FORTRAN(number_of_mpceqns_errors)( int& num_errors )
{
  num_errors = contact_mpc_eqns->Number_of_Errors();
}

void FORTRAN(number_of_voltrans_errors)( int& num_errors )
{
  num_errors = contact_volume_transfer->Number_of_Errors();
}

void FORTRAN(get_tied_kin_error_message)( int& i, char* message )
{
  const char* msg = contact_tied_kinematics->Error_Message( i );
  std::memcpy( message, msg, 80*sizeof(char) );
}

void FORTRAN(get_vol_tran_error_message)( int& i, char* message )
{
  const char* msg = contact_volume_transfer->Error_Message( i );
  std::memcpy( message, msg, 80*sizeof(char) );
}

void FORTRAN(get_mpc_eqns_error_message)( int& i, char* message )
{
  const char* msg = contact_mpc_eqns->Error_Message( i );
  std::memcpy( message, msg, 80*sizeof(char) );
}

void FORTRAN(gap_num_general_rsvars)( int& num )
{
  num = contact_gap_removal->Number_General_Restart_Variables();
}

void FORTRAN(td_num_general_rsvars)( int& num )
{
  num = contact_td_enforcement->Number_General_Restart_Variables();
}

void FORTRAN(td_penalty_num_general_rsvars)( int& num )
{
  num = contact_td_penalty->Number_General_Restart_Variables();
}

void FORTRAN(tied_kin_num_general_rsvars)( int& num )
{
  num = contact_tied_kinematics->Number_General_Restart_Variables();
}


void FORTRAN(vol_tran_num_general_rsvars)( int& num )
{
  num = contact_volume_transfer->Number_General_Restart_Variables();
}

void FORTRAN(mpc_eqns_num_general_rsvars)( int& num )
{
  num = contact_mpc_eqns->Number_General_Restart_Variables();
}

void FORTRAN(gap_num_node_rsvars)( int& num )
{
  num = contact_gap_removal->Number_Nodal_Restart_Variables();
}

void FORTRAN(td_num_node_rsvars)( int& num )
{
  num = contact_td_enforcement->Number_Nodal_Restart_Variables();
}

void FORTRAN(td_penalty_num_node_rsvars)( int& num )
{
  num = contact_td_penalty->Number_Nodal_Restart_Variables();
}


void FORTRAN(tied_kin_num_node_rsvars)( int& num )
{
  num = contact_tied_kinematics->Number_Nodal_Restart_Variables();
}

void FORTRAN(vol_tran_num_node_rsvars)( int& num )
{
  num = contact_volume_transfer->Number_Nodal_Restart_Variables();
}

void FORTRAN(mpc_eqns_num_node_rsvars)( int& num )
{
  num = contact_mpc_eqns->Number_Nodal_Restart_Variables();
}
void FORTRAN(gap_num_edge_rsvars)( int& num )
{
  num = contact_gap_removal->Number_Edge_Restart_Variables();
}

void FORTRAN(td_num_edge_rsvars)( int& num )
{
  num = contact_td_enforcement->Number_Edge_Restart_Variables();
}

void FORTRAN(td_penalty_num_edge_rsvars)( int& num )
{
  num = contact_td_penalty->Number_Edge_Restart_Variables();
}

void FORTRAN(tied_kin_num_edge_rsvars)( int& num )
{
  num = contact_tied_kinematics->Number_Edge_Restart_Variables();
}

void FORTRAN(vol_tran_num_edge_rsvars)( int& num )
{
  num = contact_volume_transfer->Number_Edge_Restart_Variables();
}

void FORTRAN(mpc_eqns_num_edge_rsvars)( int& num )
{
  num = contact_mpc_eqns->Number_Edge_Restart_Variables();
}

void FORTRAN(gap_num_face_rsvars)( int& num )
{
  num = contact_gap_removal->Number_Face_Restart_Variables();
}

void FORTRAN(td_num_face_rsvars)( int& num )
{
  num = contact_td_enforcement->Number_Face_Restart_Variables();
}

void FORTRAN(td_penalty_num_face_rsvars)( int& num )
{
  num = contact_td_penalty->Number_Face_Restart_Variables();
}

void FORTRAN(tied_kin_num_face_rsvars)( int& num )
{
  num = contact_tied_kinematics->Number_Face_Restart_Variables();
}

void FORTRAN(vol_tran_num_face_rsvars)( int& num )
{
  num = contact_volume_transfer->Number_Face_Restart_Variables();
}

void FORTRAN(mpc_eqns_num_face_rsvars)( int& num )
{
  num = contact_mpc_eqns->Number_Face_Restart_Variables();
}

void FORTRAN(gap_num_element_rsvars)( int& num )
{
  num = contact_gap_removal->Number_Element_Restart_Variables();
}

void FORTRAN(td_num_element_rsvars)( int& num )
{
  num = contact_td_enforcement->Number_Element_Restart_Variables();
}

void FORTRAN(td_penalty_num_element_rsvars)( int& num )
{
  num = contact_td_penalty->Number_Element_Restart_Variables();
}

void FORTRAN(tied_kin_num_element_rsvars)( int& num )
{
  num = contact_tied_kinematics->Number_Element_Restart_Variables();
}

void FORTRAN(vol_tran_num_element_rsvars)( int& num )
{
  num = contact_volume_transfer->Number_Element_Restart_Variables();
}

void FORTRAN(mpc_eqns_num_element_rsvars)( int& num )
{
  num = contact_mpc_eqns->Number_Element_Restart_Variables();
}

void FORTRAN(gap_extract_general_rsvars)( Real* data, int& ierror )
{
  ierror = contact_gap_removal->Extract_General_Restart_Variable( data );
}

void FORTRAN(td_extract_general_rsvars)( Real* data, int& ierror )
{
  ierror = contact_td_enforcement->Extract_General_Restart_Variable( data );
}

void FORTRAN(td_penalty_extract_general_rsvars)( Real* data, int& ierror )
{
  ierror = contact_td_penalty->Extract_General_Restart_Variable( data );
}

void FORTRAN(tied_kin_extract_general_rsvars)( Real* data, int& ierror )
{
  ierror = contact_tied_kinematics->Extract_General_Restart_Variable( data );
}


void FORTRAN(vol_tran_extract_general_rsvars)( Real* data, int& ierror )
{
  ierror = contact_volume_transfer->Extract_General_Restart_Variable( data );
}

void FORTRAN(mpc_eqns_extract_general_rsvars)( Real* data, int& ierror )
{
  ierror = contact_mpc_eqns->Extract_General_Restart_Variable( data );
}


void FORTRAN(gap_extract_node_rsvars)( int& n, Real* data, int* node_ids, int& ierror )
{
  ierror = contact_gap_removal->Extract_Nodal_Restart_Variable( n, data, node_ids );
}

void FORTRAN(td_extract_node_rsvars)( int& n, Real* data, int* node_ids, int& ierror )
{
  ierror = contact_td_enforcement->Extract_Nodal_Restart_Variable( n, data, node_ids );
}

void FORTRAN(td_penalty_extract_node_rsvars)( int& n, Real* data, int* node_ids, int& ierror )
{
  ierror = contact_td_penalty->Extract_Nodal_Restart_Variable( n, data, node_ids );
}

void FORTRAN(tied_kin_extract_node_rsvars)( int& n, Real* data, int* node_ids, int& ierror )
{
  ierror = contact_tied_kinematics->Extract_Nodal_Restart_Variable( n, data, node_ids );
}


void FORTRAN(vol_tran_extract_node_rsvars)( int& n, Real* data, int* node_ids, int& ierror )
{
  ierror = contact_volume_transfer->Extract_Nodal_Restart_Variable( n, data, node_ids );
}

void FORTRAN(mpc_eqns_extract_node_rsvars)( int& n, Real* data, int* node_ids, int& ierror )
{
  ierror = contact_mpc_eqns->Extract_Nodal_Restart_Variable( n, data, node_ids );
}

void FORTRAN(gap_extract_edge_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_gap_removal->Extract_Edge_Restart_Variable( n, data );
}

void FORTRAN(td_extract_edge_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_td_enforcement->Extract_Edge_Restart_Variable( n, data );
}

void FORTRAN(td_penalty_extract_edge_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_td_penalty->Extract_Edge_Restart_Variable( n, data );
}

void FORTRAN(tied_kin_extract_edge_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_tied_kinematics->Extract_Edge_Restart_Variable( n, data );
}

void FORTRAN(vol_tran_extract_edge_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_volume_transfer->Extract_Edge_Restart_Variable( n, data );
}

void FORTRAN(mpc_eqns_extract_edge_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_mpc_eqns->Extract_Edge_Restart_Variable( n, data );
}

void FORTRAN(gap_extract_face_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_gap_removal->Extract_Face_Restart_Variable( n, data );
}

void FORTRAN(td_extract_face_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_td_enforcement->Extract_Face_Restart_Variable( n, data );
}

void FORTRAN(td_penalty_extract_face_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_td_penalty->Extract_Face_Restart_Variable( n, data );
}

void FORTRAN(tied_kin_extract_face_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_tied_kinematics->Extract_Face_Restart_Variable( n, data );
}

void FORTRAN(vol_tran_extract_face_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_volume_transfer->Extract_Face_Restart_Variable( n, data );
}

void FORTRAN(mpc_eqns_extract_face_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_mpc_eqns->Extract_Face_Restart_Variable( n, data );
}

void FORTRAN(gap_extract_element_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_gap_removal->Extract_Element_Restart_Variable( n, data );
}

void FORTRAN(td_extract_element_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_td_enforcement->Extract_Element_Restart_Variable( n, data );
}

void FORTRAN(td_penalty_extract_element_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_td_penalty->Extract_Element_Restart_Variable( n, data );
}

void FORTRAN(tied_kin_extract_element_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_tied_kinematics->Extract_Element_Restart_Variable( n, data );
}

void FORTRAN(vol_tran_extract_element_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_volume_transfer->Extract_Element_Restart_Variable( n, data );
}

void FORTRAN(mpc_eqns_extract_element_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_mpc_eqns->Extract_Element_Restart_Variable( n, data );
}

void FORTRAN(gap_implant_general_rsvars)( Real* data, int& ierror )
{
  ierror = contact_gap_removal->Implant_General_Restart_Variable( data );
}

void FORTRAN(td_implant_general_rsvars)( Real* data, int& ierror )
{
  ierror = contact_td_enforcement->Implant_General_Restart_Variable( data );
}

void FORTRAN(td_penalty_implant_general_rsvars)( Real* data, int& ierror )
{
  ierror = contact_td_penalty->Implant_General_Restart_Variable( data );
}

void FORTRAN(tied_kin_implant_general_rsvars)( Real* data, int& ierror )
{
  ierror = contact_tied_kinematics->Implant_General_Restart_Variable( data );
}

void FORTRAN(vol_tran_implant_general_rsvars)( Real* data, int& ierror )
{
  ierror = contact_volume_transfer->Implant_General_Restart_Variable( data );
}

void FORTRAN(mpc_eqns_implant_general_rsvars)( Real* data, int& ierror )
{
  ierror = contact_mpc_eqns->Implant_General_Restart_Variable( data );
}

void FORTRAN(gap_implant_node_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_gap_removal->Implant_Nodal_Restart_Variable( n, data );
}

void FORTRAN(td_implant_node_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_td_enforcement->Implant_Nodal_Restart_Variable( n, data );
}

void FORTRAN(td_penalty_implant_node_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_td_penalty->Implant_Nodal_Restart_Variable( n, data );
}


void FORTRAN(tied_kin_implant_node_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_tied_kinematics->Implant_Nodal_Restart_Variable( n, data );
}

void FORTRAN(vol_tran_implant_node_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_volume_transfer->Implant_Nodal_Restart_Variable( n, data );
}

void FORTRAN(mpc_eqns_implant_node_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_mpc_eqns->Implant_Nodal_Restart_Variable( n, data );
}

void FORTRAN(gap_implant_edge_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_gap_removal->Implant_Edge_Restart_Variable( n, data );
}

void FORTRAN(td_implant_edge_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_td_enforcement->Implant_Edge_Restart_Variable( n, data );
}

void FORTRAN(td_penalty_implant_edge_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_td_penalty->Implant_Edge_Restart_Variable( n, data );
}


void FORTRAN(tied_kin_implant_edge_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_tied_kinematics->Implant_Edge_Restart_Variable( n, data );
}

void FORTRAN(vol_tran_implant_edge_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_volume_transfer->Implant_Edge_Restart_Variable( n, data );
}

void FORTRAN(mpc_eqns_implant_edge_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_mpc_eqns->Implant_Edge_Restart_Variable( n, data );
}

void FORTRAN(gap_implant_face_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_gap_removal->Implant_Face_Restart_Variable( n, data );
}

void FORTRAN(td_implant_face_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_td_enforcement->Implant_Face_Restart_Variable( n, data );
}

void FORTRAN(td_penalty_implant_face_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_td_penalty->Implant_Face_Restart_Variable( n, data );
}


void FORTRAN(tied_kin_implant_face_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_tied_kinematics->Implant_Face_Restart_Variable( n, data );
}

void FORTRAN(vol_tran_implant_face_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_volume_transfer->Implant_Face_Restart_Variable( n, data );
}

void FORTRAN(mpc_eqns_implant_face_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_mpc_eqns->Implant_Face_Restart_Variable( n, data );
}

void FORTRAN(gap_implant_element_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_gap_removal->Implant_Element_Restart_Variable( n, data );
}

void FORTRAN(td_implant_element_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_td_enforcement->Implant_Element_Restart_Variable( n, data );
}

void FORTRAN(td_penalty_implant_element_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_td_penalty->Implant_Element_Restart_Variable( n, data );
}

void FORTRAN(tied_kin_implant_element_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_tied_kinematics->Implant_Element_Restart_Variable( n, data );
}

void FORTRAN(vol_tran_implant_element_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_volume_transfer->Implant_Element_Restart_Variable( n, data );
}

void FORTRAN(mpc_eqns_implant_element_rsvars)( int& n, Real* data, int& ierror )
{
  ierror = contact_mpc_eqns->Implant_Element_Restart_Variable( n, data );
}

void FORTRAN(gap_complete_restart)(int& ierror)
{
  ierror = contact_gap_removal->Complete_Restart();
}

void FORTRAN(td_complete_restart)(int& ierror)
{
  ierror = contact_td_enforcement->Complete_Restart();
}

void FORTRAN(tied_kin_complete_restart)(int& ierror)
{
  ierror = contact_tied_kinematics->Complete_Restart();
}

void FORTRAN(vol_tran_complete_restart)(int& ierror)
{
  ierror = contact_volume_transfer->Complete_Restart();
}

void FORTRAN(mpc_eqns_complete_restart)(int& ierror)
{
  ierror = contact_mpc_eqns->Complete_Restart();
}

void FORTRAN(td_add_enf_model)(int& type, int& id, int* int_data,
			       Real* real_data, int& error )
{
  error = contact_td_enforcement->
    Add_Enforcement_Model( (ContactEnforcement::Enforcement_Model_Types) type,
			   id, int_data, real_data );
}

void FORTRAN(td_penalty_add_enf_model)(int& type, int& id, int* int_data,
			       Real* real_data, int& error )
{
  error = contact_td_penalty->
    Add_Enforcement_Model( (ContactEnforcement::Enforcement_Model_Types) type,
			   id, int_data, real_data );
}

void FORTRAN(user_fn_initialize_model)(int& id, CONTACT_INIT_MODEL_FN* fn)
{
  contact_td_enforcement->User_Initialize_Model_Fn(id, fn);
}
void FORTRAN(user_fn_initialize_time_step)(int& id, CONTACT_INIT_TIME_STEP_FN *fn)
{
  contact_td_enforcement->User_Initialize_Time_Step_Fn(id, fn);
}
void FORTRAN(user_fn_init_node_state_data)(int& id, CONTACT_INIT_NODE_STATE_DATA_FN *fn)
{
  contact_td_enforcement->User_Initialize_Node_State_Data_Fn(id, fn);
}
void FORTRAN(user_fn_limit_force)(int& id, CONTACT_LIMIT_FORCE_FN *fn)
{
  contact_td_enforcement->User_Limit_Force_Fn(id, fn);
}
void FORTRAN(user_fn_active)(int& id, CONTACT_INTERACTION_ACTIVE_FN *fn)
{
  contact_td_enforcement->User_Active_Fn(id, fn);
}
void FORTRAN(user_fn_interaction_type)(int& id, CONTACT_INTERACTION_TYPE_FN *fn)
{
  contact_td_enforcement->User_Interaction_Type_Fn(id, fn);
}
