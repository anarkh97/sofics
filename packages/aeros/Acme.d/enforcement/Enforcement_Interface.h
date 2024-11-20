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


#include "Search_Interface.h"
#include "ContactSearch.h"
#include "ContactTDEnforcement.h"
#include "ContactTDEnfPenalty.h"
#include "ContactGapRemoval.h"
#include "ContactTiedKinematics.h"
#include "ContactVolumeTransfer.h"
#include "ContactMPCs.h"


static ContactTDEnforcement* contact_td_enforcement = NULL;
static ContactTDEnfPenalty* contact_td_penalty = NULL;
static ContactGapRemoval* contact_gap_removal = NULL;
static ContactTiedKinematics* contact_tied_kinematics = NULL;
static ContactVolumeTransfer* contact_volume_transfer = NULL;
static ContactMPCs* contact_mpc_eqns = NULL;

extern "C" {

  void FORTRAN(number_of_enforcement_errors)( int& num_errors );
  
  void FORTRAN(get_enforcement_error_message)( int& i, char* message );
  
  void FORTRAN( build_td_enforcement )( double* Enforcement_Data,
                                        int* get_cvars,
                                        int* plot_force,
					int* error );
  void FORTRAN( build_td_enf_restart )( double* restart_data, 
                                        int* get_cvars,
                                        int* plot_force,
                                        int* error );
  void FORTRAN( td_enf_restart_size )( int& );
  void FORTRAN( td_enf_extract_restart )( double*, int& );
  void FORTRAN( set_td_iterations )( int*, int* );
  void FORTRAN( set_td_convg_tol )( double*, int* );
  void FORTRAN( set_td_enf_symm_nodes )( int&, int&, int& error );
  void FORTRAN( compute_td_contact_force )( double& dt_old,
					    double& dt,
					    double* mass,
					    double* density,
					    double* wavespeed,
					    double* forces,
					    int* error );
  void FORTRAN( get_td_plot_variable )( int& plot_var, double* data, 
					int* error );
  void FORTRAN( get_td_global_variable)(int& plot_var, double  data, 
					int* error );
  void FORTRAN( cleanup_td_enforcement )();
  
  void FORTRAN( build_td_penalty )( double* Enforcement_Data,
                                       int* get_cvars,
                                       int* get_plot_force,
				       int* error );
  void FORTRAN( build_td_penalty_restart )( double* restart_data,
                                               int* get_cvars,
                                               int* get_plot_force,
                                               int* error );
  void FORTRAN( td_penalty_restart_size )( int& );
  void FORTRAN( td_penalty_extract_restart )( double*, int& );
  void FORTRAN( set_td_penalty_iterations )( int*, int* );
  void FORTRAN( set_td_penalty_convg_tol )( double*, int* );
  void FORTRAN( compute_td_penalty_contact_force )( double& dt_old,
					    double& dt,
					    double* mass,
					    double* density,
					    double* wavespeed,
					    double* forces,
					    int* error );
  void FORTRAN( cleanup_td_penalty)();

  void FORTRAN( build_gap_removal )( double* Enforcement_Data,
				     int* error );
  void FORTRAN( build_gap_removal_restart )( double* restart_data, int* error );
  void FORTRAN( gap_removal_restart_size )( int& );
  void FORTRAN( gap_removal_extract_restart )( double*, int& );
  
  void FORTRAN( compute_gap_removal )( int& max_iterations,
				       double& trivial_gap,
				       double* displc,
				       int* error );
  
  void FORTRAN( cleanup_gap_removal )();
  
  void FORTRAN( build_tied_kin )( double* Enforcement_Data,
			          int* error );
  void FORTRAN( build_tied_kin_restart )( double* restart_data, int* error );
  void FORTRAN( tied_kin_restart_size )( int& );
  void FORTRAN( tied_kin_extract_restart )( double*, int& );
  
  void FORTRAN( compute_tied_kin_position )( double* position,
					     int* error );
  
  void FORTRAN( cleanup_tied_kinematics )();
  

  // build and cleanup of volume transfer object
  void FORTRAN( build_vol_tran )( double* Enforcement_Data,
				  int* error );
  void FORTRAN( compute_vol_tran )( int&, int&, double*, double*, double*, double*, double*, int* );
  void FORTRAN( build_vol_tran_restart )( double* restart_data, int* error );
  void FORTRAN( vol_tran_restart_size )( int& );
  void FORTRAN( vol_tran_extract_restart )( double*, int& );
  void FORTRAN( cleanup_vol_tran )();


  // mpc enforcement object utility functions
  void FORTRAN( build_mpc_eqns )( double* Enforcement_Data, int* error );
  void FORTRAN( compute_mpc_eqns )( int* error );
  void FORTRAN( build_mpc_eqns_restart )( double* restart_data, int* error );
  void FORTRAN( mpc_eqns_restart_size )( int& );
  void FORTRAN( mpc_eqns_extract_restart )( double*, int& );
  void FORTRAN( get_num_mpceqns )( int& nmpceqn, int* error );
  void FORTRAN( get_mpc_eqns )( int&, int*, int*, int*, int*, 
				int*, int*, int*, Real*, int* );
  
  void FORTRAN( cleanup_mpc_eqns )();

  void FORTRAN( number_of_td_errors )(int& num_errors );
  void FORTRAN( number_of_td_penalty_errors )(int& num_errors );
  void FORTRAN( number_of_gap_errors )( int& num_errors );
  void FORTRAN( number_of_tied_errors )( int& num_errors );
  void FORTRAN( number_of_voltrans_errors )( int& num_errors );
  void FORTRAN( number_of_mpc_eqns_errors )( int& num_errors );

  void FORTRAN( get_td_error_message )( int& num_errors, char* message );
  void FORTRAN( get_td_penalty_error_message )( int& num_errors, char* message );
  void FORTRAN( get_gap_error_message )( int& num_errors, char* message );
  void FORTRAN( get_tied_kin_error_message )( int& num_errors, char* message );
  void FORTRAN( get_vol_tran_error_message )( int& num_errors, char* message );
  void FORTRAN( get_mpc_eqns_error_message )( int& num_errors, char* message );

  void FORTRAN(gap_num_general_rsvars)( int& num );
  void FORTRAN(td_num_general_rsvars)( int& num );
  void FORTRAN(td_penalty_num_general_rsvars)( int& num );
  void FORTRAN(tied_kin_num_general_rsvars)( int& num );
  void FORTRAN(vol_tran_num_general_rsvars)( int& num );
  void FORTRAN(mpc_eqns_num_general_rsvars)( int& num );

  void FORTRAN(gap_num_node_rsvars)( int& num );
  void FORTRAN(td_num_node_rsvars)( int& num );
  void FORTRAN(td_penalty_num_node_rsvars)( int& num );
  void FORTRAN(tied_kin_num_node_rsvars)( int& num );
  void FORTRAN(vol_tran_num_node_rsvars)( int& num );
  void FORTRAN(mpc_eqns_num_node_rsvars)( int& num );

  void FORTRAN(gap_num_edge_rsvars)( int& num );
  void FORTRAN(td_num_edge_rsvars)( int& num );
  void FORTRAN(td_penalty_num_edge_rsvars)( int& num );
  void FORTRAN(tied_kin_num_edge_rsvars)( int& num );
  void FORTRAN(vol_tran_num_edge_rsvars)( int& num );
  void FORTRAN(mpc_eqns_num_edge_rsvars)( int& num );

  void FORTRAN(gap_num_face_rsvars)( int& num );
  void FORTRAN(td_num_face_rsvars)( int& num );
  void FORTRAN(td_penalty_num_face_rsvars)( int& num );
  void FORTRAN(tied_kin_num_face_rsvars)( int& num );
  void FORTRAN(vol_tran_num_face_rsvars)( int& num );
  void FORTRAN(mpc_eqns_num_face_rsvars)( int& num );

  void FORTRAN(gap_num_element_rsvars)( int& num );
  void FORTRAN(td_num_element_rsvars)( int& num );
  void FORTRAN(td_penalty_num_element_rsvars)( int& num );
  void FORTRAN(tied_kin_num_element_rsvars)( int& num );
  void FORTRAN(vol_tran_num_element_rsvars)( int& num );
  void FORTRAN(mpc_eqns_num_element_rsvars)( int& num );

  void FORTRAN(gap_extract_general_rsvars)( Real* data, int& ierror );
  void FORTRAN(td_extract_general_rsvars)( Real* data, int& ierror );
  void FORTRAN(td_penalty_extract_general_rsvars)( Real* data, int& ierror );
  void FORTRAN(tied_kin_extract_general_rsvars)( Real* data, int& ierror );
  void FORTRAN(vol_tran_extract_general_rsvars)( Real* data, int& ierror );
  void FORTRAN(mpc_eqns_extract_general_rsvars)( Real* data, int& ierror );

  void FORTRAN(gap_extract_node_rsvars)( int& n, Real* data, int* node_ids,
                                         int& ierror );
  void FORTRAN(td_extract_node_rsvars)( int& n, Real* data, int* node_ids,
                                        int& ierror );
  void FORTRAN(td_penalty_extract_node_rsvars)( int& n, Real* data, 
                                                   int* node_ids, int& ierror );
  void FORTRAN(tied_kin_extract_node_rsvars)( int& n, Real* data, 
					      int* node_ids, int& ierror );
  void FORTRAN(vol_tran_extract_node_rsvars)( int& n, Real* data,
					      int* node_ids, int& ierror );
  void FORTRAN(mpc_eqns_extract_node_rsvars)( int& n, Real* data,
					      int* node_ids, int& ierror );

  void FORTRAN(td_extract_edge_rsvars)( int& n, Real* data, int& ierror );
  void FORTRAN(td_penalty_extract_edge_rsvars)( int& n, Real* data, int& ierror );
  void FORTRAN(gap_extract_edge_rsvars)( int& n, Real* data, int& ierror );
  void FORTRAN(tied_kin_extract_edge_rsvars)( int& n, Real* data, 
					      int& ierror );
  void FORTRAN(vol_tran_extract_edge_rsvars)( int& n, Real* data, 
					       int& ierror );
  void FORTRAN(mpc_eqns_extract_edge_rsvars)( int& n, Real* data,
					       int& ierror );

  void FORTRAN(gap_extract_face_rsvars)( int& n, Real* data, int& ierror );
  void FORTRAN(td_extract_face_rsvars)( int& n, Real* data, int& ierror );
  void FORTRAN(td_penalty_extract_face_rsvars)( int& n, Real* data, int& ierror );
  void FORTRAN(tied_kin_extract_face_rsvars)( int& n, Real* data, 
					      int& ierror );
  void FORTRAN(vol_tran_extract_face_rsvars)( int& n, Real* data,
					      int& ierror );
  void FORTRAN(mpc_eqns_extract_face_rsvars)( int& n, Real* data,
					      int& ierror );

  void FORTRAN(gap_extract_element_rsvars)( int& n, Real* data, int& ierror );
  void FORTRAN(td_extract_element_rsvars)( int& n, Real* data, int& ierror );
  void FORTRAN(td_penalty_extract_element_rsvars)( int& n, Real* data, int& ierror );
  void FORTRAN(tied_kin_extract_element_rsvars)( int& n, Real* data, 
					      int& ierror );
  void FORTRAN(vol_tran_extract_element_rsvars)( int& n, Real* data,
						 int& ierror );
  void FORTRAN(mpc_eqns_extract_element_rsvars)( int& n, Real* data,
						 int& ierror );

  void FORTRAN(gap_implant_general_rsvars)( Real* data, int& ierror );
  void FORTRAN(td_implant_general_rsvars)( Real* data, int& ierror );
  void FORTRAN(td_penalty_implant_general_rsvars)( Real* data, int& ierror );
  void FORTRAN(tied_kin_implant_general_rsvars)( Real* data, int& ierror );
  void FORTRAN(vol_tran_implant_general_rsvars)( Real* data, int& ierror );
  void FORTRAN(mpc_eqns_implant_general_rsvars)( Real* data, int& ierror );

  void FORTRAN(gap_implant_node_rsvars)( int& n, Real* data, int& ierror );
  void FORTRAN(td_implant_node_rsvars)( int& n, Real* data, int& ierror );
  void FORTRAN(td_penalty_implant_node_rsvars)( int& n, Real* data, int& ierror );
  void FORTRAN(tied_kin_implant_node_rsvars)( int& n, Real* data, 
					      int& ierror );
  void FORTRAN(vol_tran_implant_node_rsvars)( int& n, Real* data,
					       int& ierror );
  void FORTRAN(mpc_eqns_implant_node_rsvars)( int& n, Real* data,
					       int& ierror );

  void FORTRAN(gap_implant_edge_rsvars)( int& n, Real* data, int& ierror );
  void FORTRAN(td_implant_edge_rsvars)( int& n, Real* data, int& ierror );
  void FORTRAN(td_penalty_implant_edge_rsvars)( int& n, Real* data, int& ierror );
  void FORTRAN(tied_kin_implant_edge_rsvars)( int& n, Real* data, 
					      int& ierror );
  void FORTRAN(vol_tran_implant_edge_rsvars)( int& n, Real* data, 
					       int& ierror );
  void FORTRAN(mpc_eqns_implant_edge_rsvars)( int& n, Real* data, 
					       int& ierror );

  void FORTRAN(gap_implant_face_rsvars)( int& n, Real* data, int& ierror );
  void FORTRAN(td_implant_face_rsvars)( int& n, Real* data, int& ierror );
  void FORTRAN(td_penalty_implant_face_rsvars)( int& n, Real* data, int& ierror );
  void FORTRAN(tied_kin_implant_face_rsvars)( int& n, Real* data, 
					      int& ierror );
  void FORTRAN(vol_tran_implant_face_rsvars)( int& n, Real* data, 
					       int& ierror );
  void FORTRAN(mpc_eqns_implant_face_rsvars)( int& n, Real* data, 
					       int& ierror );
  
  void FORTRAN(gap_implant_element_rsvars)( int& n, Real* data, int& ierror );
  void FORTRAN(td_implant_element_rsvars)( int& n, Real* data, int& ierror );
  void FORTRAN(td_penalty_implant_element_rsvars)( int& n, Real* data, int& ierror );
  void FORTRAN(tied_kin_implant_element_rsvars)( int& n, Real* data, 
					      int& ierror );
  void FORTRAN(vol_tran_implant_element_rsvars)( int& n, Real* data, 
					       int& ierror );
  void FORTRAN(mpc_eqns_implant_element_rsvars)( int& n, Real* data, 
					       int& ierror );
  
  void FORTRAN(gap_complete_restart)(int& ierror);
  void FORTRAN(td_complete_restart)(int& ierror);
  void FORTRAN(tied_kin_complete_restart)(int& ierror);
  void FORTRAN(vol_tran_complete_restart)(int& ierror);
  void FORTRAN(mpc_eqns_complete_restart)(int& ierror);

  void FORTRAN(td_add_enf_model)(int& type, int& id, int* int_data,
				 Real* real_data, int& error);
  void FORTRAN(td_penalty_add_enf_model)(int& type, int& id, int* int_data,
				 Real* real_data, int& error);

  void FORTRAN(user_fn_initialize_model)(int& id, CONTACT_INIT_MODEL_FN* fn);
  void FORTRAN(user_fn_initialize_time_step)(int& id, CONTACT_INIT_TIME_STEP_FN* fn);
  void FORTRAN(user_fn_init_node_state_data)(int& id, CONTACT_INIT_NODE_STATE_DATA_FN *fn);
  void FORTRAN(user_fn_limit_force)(int& id, CONTACT_LIMIT_FORCE_FN *fn);
  void FORTRAN(user_fn_active)(int& id, CONTACT_INTERACTION_ACTIVE_FN *fn);
  void FORTRAN(user_fn_interaction_type)(int& id, CONTACT_INTERACTION_TYPE_FN *fn);
}
