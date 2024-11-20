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


// NOTE to turn off debugging spew set LOCAL_PRINT_FLAG to 0 below

// To Do:
// * turn "ALL_ACTIVE" back on after time-to-contact in search is updated
// *     add zeroing of IDEAL_GAP and GAP_TO_ENFORCE
// * move 6 conditions in search into TDEnf (after head search is updated)
// * undo effects of ContactSearch::Partition_Between_GapCur_and_GapOld
// * try a time-to-contact based on a interaction surface rather than physical
// * make tolerances/parameters dimensionless
// * look at multiple with friction
// * why does using all interactions have problems(vs just currently contacting)

#include "Contact_Communication.h"
#include "ContactTDEnfModel.h"
#include "ContactErrors.h"
#include "ContactTDEnforcement.h"
#include "ContactSearch.h"
#include "ContactTopology.h"
#include "ContactShellHandler.h"
#include "ContactShellNode.h"
#include "ContactNodeFaceInteraction.h"
#include "ContactNodeSurfaceInteraction.h"
#include "ContactScratchManager.h"
#include "ContactTDTied.h"
#include "ContactTDFrictionless.h"
#include "ContactHexOverlap.h"
#include "cstring"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <algorithm>

using namespace std;

// Internal constants
#undef  NDIM
#define NDIM 3
#undef  ZERO_TOL
#define ZERO_TOL 1.0E-14
#undef  SMALL
#define SMALL 1.0e-6
#undef  PARTITION_ITERATIONS 

#define PARTITION_ITERATIONS 100
#undef  FRAC_CUR_GAP
#define FRAC_CUR_GAP 0.10
#undef  DET_TOL
#define DET_TOL 1.0e-30
#undef  COLLINEARITY_TOL
#define COLLINEARITY_TOL 0.99e0
#undef  CONVERGENCE_TOL
#define CONVERGENCE_TOL 1.0e-10
#undef  FAIL_TOL
#define FAIL_TOL 0.01

// Internal flags 
// LOCAL_PRINT_FLAG = 0 turns off all internal debugging prints
#undef  LOCAL_PRINT_FLAG 
#define LOCAL_PRINT_FLAG 0

#define PJSA_FIX
#ifdef USE_EIGEN3
#include <Eigen/Dense>
#endif
extern int contactPrintFlag;

#ifdef CONTACT_DEBUG_NODE
#include "ContactParOStream.h"
#elif LOCAL_PRINT_FLAG > 0
#include "ContactParOStream.h"
#endif

/* SYMBOLS ______________________________________________________
 * x = position
 * u = displacement
 * v = velocity
 * a = acceleration
 * f = force
 * m = mass
 * nm = outward unit normal vector ( f. nm < 0 for compression )
 *      (nm is typically outward for the master face and the
 *       force is typically stored with the slave node)
 * uj = jump in displacement i.e. relative displacement
 * vj = jump in velocity i.e. relative velocity
 * gap = uj . nm  , where  gap > 0 for separation
 * slip = vj * dt
 * kp ? = kinematic partition
 *     ... subscripts ....
 * f_n = the normal component of force
 * f_t = the tangential force
 * x_sn = slave node position
 * x_mn = master node position
 * x_cp = contact point position
 * x_cur = "CURRENT" position ?? from last time step after contact correction
 * x_new = "NEW" position ?? from current iteration
 * x_pre = "PREDICTED" position ?? from central difference predictor w/o contact *         force.
 * _____________________________________________________________*/


/* ASSUMPTIONS__________________________________________________
 *  for node-surface interactions the kinematic partition is 
 *  assumed to be 1.0
 * _____________________________________________________________*/



ContactTDEnforcement::ContactTDEnforcement( const Real* enf_data, 
					    ContactSearch* Search,
					    ContactSearch::ContactErrorCode& error,
                                            bool get_cvars,
					    bool calc_plot_force)
  : ContactEnforcement( error, Search, ContactEnforcement::TDEnforcement,
			NSIZED, enf_data, true, FRICTION_MODEL_ID ),
    SAVE_CVARS(get_cvars),
    SAVE_CGVARS(get_cvars),
    CALC_PLOT_FORCE(calc_plot_force)
{  
  if (error != ContactSearch::NO_ERROR) return;
  only_frictionless_or_tied = false;

  search->Get_Search_Option(ContactSearch::NEW_TIED_ENFORCEMENT, new_tied_enforcement, NULL);

  // plot vars for regression
  number_global_plot_vars  = 2;
  number_nodal_plot_vars   = 6;
  number_element_plot_vars = 0;

  max_interactions = search->Max_Interactions_Per_Node();
  
  plot_force = NULL;
  if(CALC_PLOT_FORCE) {
    plot_force = new Real[number_of_nodes*number_nodal_plot_vars];
    std::memset( plot_force, 0, number_of_nodes*number_nodal_plot_vars*sizeof(Real) );
  }

  disp_old_aug = new Real[number_of_nodes*3];
  std::memset( disp_old_aug, 0, number_of_nodes*3*sizeof(Real) );
  num_iterations = 1;
  time_mult_old_aug = 0.0;
  if( SAVE_CVARS ){
    conface                   = new Real[number_of_nodes];

    normal_force_mag          = new Real[number_of_nodes];
    normal_traction_mag       = new Real[number_of_nodes];
    tangential_force_mag      = new Real[number_of_nodes];
    tangential_traction_mag   = new Real[number_of_nodes];

    normal_dir_x              = new Real[number_of_nodes];
    normal_dir_y              = new Real[number_of_nodes];
    normal_dir_z              = new Real[number_of_nodes];
  
    tangential_dir_x          = new Real[number_of_nodes];
    tangential_dir_y          = new Real[number_of_nodes];
    tangential_dir_z          = new Real[number_of_nodes];

    slipmag                   = new Real[number_of_nodes];
    nodal_dissipation         = new Real[number_of_nodes];
    nodal_dissipation_density = new Real[number_of_nodes];
    contact_area              = new Real[number_of_nodes];
    plot_gap_cur              = new Real[number_of_nodes];
    plot_gap_old              = new Real[number_of_nodes];
    kinematic_partition       = new Real[number_of_nodes];
  } else {
    conface                   = NULL;
    normal_force_mag          = NULL;
    normal_traction_mag       = NULL;
    tangential_force_mag      = NULL;
    tangential_traction_mag   = NULL;
    
    normal_dir_x              = NULL;
    normal_dir_y              = NULL;
    normal_dir_z              = NULL;

    tangential_dir_x          = NULL;
    tangential_dir_y          = NULL;
    tangential_dir_z          = NULL;

    slipmag                   = NULL;
    nodal_dissipation         = NULL;
    nodal_dissipation_density = NULL;
    contact_area              = NULL;
    plot_gap_cur              = NULL;
    plot_gap_old              = NULL;
    kinematic_partition       = NULL;
  }

  have_symmetry_nodes = false;
  symmetry_node_1_exodus_id = -1;
  symmetry_node_2_exodus_id = -1;

  convergence_tolerance = CONVERGENCE_TOL;

#ifdef CONTACT_TIMINGS
  enforcement_time               = Timer().Register_Timer( CString("TD Enforcement Time") );
  initialization_time            = Timer().Register_Timer( CString("  TD Initialization") );
  Register_Enforcement_Timers();
  scratch_allocation_time        = Timer().Register_Timer( CString("    Scratch Allocation") );
  copy_nodal_to_scratch          = Timer().Register_Timer( CString("    Copy Nodal Var to Scratch Var") );
  init_host_code_vars            = Timer().Register_Timer( CString("    Init Host Code Vars") );
  setup_interaction_list         = Timer().Register_Timer( CString("    Setup Entity Interaction List") );
  td_enforcement_setup           = Timer().Register_Timer( CString("    TD Enforcement Model Setup") );


  iteration_time                 = Timer().Register_Timer( CString("  TD Iteration Loop") );
  compute_gap_time               = Timer().Register_Timer( CString("    Compute Gap") );

  accel_cor_time                 = Timer().Register_Timer( CString("    Accleration Correction") );
  response_prediction_time       = Timer().Register_Timer( CString("      Response Prediction") );
  response_correction_time       = Timer().Register_Timer( CString("      Response Correction") );
  swap_add_time                  = Timer().Register_Timer( CString("      Swap Add") );
  contact_force_time             = Timer().Register_Timer( CString("      Contact Force") );
  response_constitutive_correction_time = 
                                   Timer().Register_Timer( CString("      Response Constitutive Correction") );
  force_adjust_time              = Timer().Register_Timer( CString("      Force Adjust") );

  penalty_time                   = Timer().Register_Timer( CString("    Penalty Forces") );
  assemble_forces_time           = Timer().Register_Timer( CString("    Assemble Forces") );
  kinematic_consistent_time      = Timer().Register_Timer( CString("    Kinematic Consistent") );
  update_config_time             = Timer().Register_Timer( CString("    Update Configuration") );
  extra_ghost_len_time           = Timer().Register_Timer( CString("  TD Extra Ghosting Length") );
  cleanup_time                   = Timer().Register_Timer( CString("  TD Clean Up") );
#endif

  Verify_Enforcement_Data();
  
  error = ContactSearch::NO_ERROR;
}

ContactTDEnforcement::ContactTDEnforcement( ContactSearch* Search,
					    const Real* restart_data,
			     ContactSearch::ContactErrorCode& error,
                             bool get_cvars,
                             bool calc_plot_force)
  : ContactEnforcement( error, Search, ContactEnforcement::TDEnforcement,
			restart_data ),
    SAVE_CVARS(get_cvars),
    SAVE_CGVARS(get_cvars),
    CALC_PLOT_FORCE(calc_plot_force)
{  
  if (error != ContactSearch::NO_ERROR) return;
  only_frictionless_or_tied = false;

  search->Get_Search_Option(ContactSearch::NEW_TIED_ENFORCEMENT, new_tied_enforcement, NULL);

  // plot vars for regression
  number_global_plot_vars  = 2;
  number_nodal_plot_vars   = 6;
  number_element_plot_vars = 0;

  if( SAVE_CVARS ){
    conface                   = new Real[number_of_nodes];

    normal_force_mag          = new Real[number_of_nodes];
    normal_traction_mag       = new Real[number_of_nodes];
    tangential_force_mag      = new Real[number_of_nodes];
    tangential_traction_mag   = new Real[number_of_nodes];

    normal_dir_x              = new Real[number_of_nodes];
    normal_dir_y              = new Real[number_of_nodes];
    normal_dir_z              = new Real[number_of_nodes];

    tangential_dir_x          = new Real[number_of_nodes];
    tangential_dir_y          = new Real[number_of_nodes];
    tangential_dir_z          = new Real[number_of_nodes];


    slipmag                   = new Real[number_of_nodes];
    nodal_dissipation         = new Real[number_of_nodes];
    nodal_dissipation_density = new Real[number_of_nodes];
    contact_area              = new Real[number_of_nodes];
    plot_gap_cur              = new Real[number_of_nodes];
    plot_gap_old              = new Real[number_of_nodes];
    kinematic_partition       = new Real[number_of_nodes];
  } else {
    conface                   = NULL;

    normal_force_mag          = NULL;
    normal_traction_mag       = NULL;
    tangential_force_mag      = NULL;
    tangential_traction_mag   = NULL;

    normal_dir_x              = NULL;
    normal_dir_y              = NULL;
    normal_dir_z              = NULL;

    tangential_dir_x          = NULL;
    tangential_dir_y          = NULL;
    tangential_dir_z          = NULL;

    slipmag                   = NULL;
    nodal_dissipation         = NULL;
    nodal_dissipation_density = NULL;
    contact_area              = NULL;
    plot_gap_cur              = NULL;
    plot_gap_old              = NULL;
    kinematic_partition       = NULL;
  }
  num_iterations = 1;

  max_interactions = search->Max_Interactions_Per_Node(); 

  plot_force = NULL;
  if(CALC_PLOT_FORCE) {
    plot_force = new Real[number_of_nodes*number_nodal_plot_vars];
  }
  // Get restart data. First need to read base class data since
  // it was written to the binary file first.
  const Real* td_restart_data = restart_data + 
    ContactEnforcement::Restart_Size();

  // Get time multiplier for last step. Needed for augmented search.
  time_mult_old_aug = *td_restart_data++;

  // Get enforcement interations
  num_iterations = static_cast<int>(*td_restart_data);
  ++td_restart_data;

  // Get old displacements. Needed for augmented search.
  int number_of_data = static_cast<int>(*td_restart_data);
  ++td_restart_data;
  disp_old_aug = new Real[number_of_data];
  std::memcpy( disp_old_aug, td_restart_data, number_of_data*sizeof(Real) );
  td_restart_data += number_of_data;

  // Get the symmetry node information
  if( static_cast<int>(*td_restart_data) == 0 )
    have_symmetry_nodes = false;
  else
    have_symmetry_nodes = true;
  ++td_restart_data;
  symmetry_node_1_exodus_id = static_cast<int>(*td_restart_data);
  ++td_restart_data;
  symmetry_node_2_exodus_id = static_cast<int>(*td_restart_data);
  ++td_restart_data;

  convergence_tolerance = CONVERGENCE_TOL;

#ifdef CONTACT_TIMINGS
  enforcement_time               = Timer().Register_Timer( CString("TD Enforcement Time") );
  initialization_time            = Timer().Register_Timer( CString("  TD Initialization") );
  Register_Enforcement_Timers();
  scratch_allocation_time        = Timer().Register_Timer( CString("    Scratch Allocation") );
  copy_nodal_to_scratch          = Timer().Register_Timer( CString("    Copy Nodal Var to Scratch Var") );
  init_host_code_vars            = Timer().Register_Timer( CString("    Init Host Code Vars") );
  setup_interaction_list         = Timer().Register_Timer( CString("    Setup Entity Interaction List") );
  td_enforcement_setup           = Timer().Register_Timer( CString("    TD Enforcement Model Setup") );
  iteration_time                 = Timer().Register_Timer( CString("  TD Iteration Loop") );
  compute_gap_time               = Timer().Register_Timer( CString("    Compute Gap") );
  accel_cor_time                 = Timer().Register_Timer( CString("    Accleration Correction") );
  response_prediction_time       = Timer().Register_Timer( CString("      Response Prediction") );
  response_correction_time       = Timer().Register_Timer( CString("      Response Correction") );
  swap_add_time                  = Timer().Register_Timer( CString("      Swap Add") );
  contact_force_time             = Timer().Register_Timer( CString("      Contact Force") );
  response_constitutive_correction_time = 
                                   Timer().Register_Timer( CString("      Response Constitutive Correction") );
  force_adjust_time              = Timer().Register_Timer( CString("      Force Adjust") );
  penalty_time                   = Timer().Register_Timer( CString("    Penalty Forces") );
  assemble_forces_time           = Timer().Register_Timer( CString("    Assemble Forces") );
  kinematic_consistent_time      = Timer().Register_Timer( CString("    Kinematic Consistent") );
  update_config_time             = Timer().Register_Timer( CString("    Update Configuration") );
  extra_ghost_len_time           = Timer().Register_Timer( CString("  TD Extra Ghosting Length") );
  cleanup_time                   = Timer().Register_Timer( CString("  TD Clean Up") );
#endif

  error = ContactSearch::NO_ERROR;
}

ContactTDEnforcement::~ContactTDEnforcement()
{
  if (CALC_PLOT_FORCE) delete [] plot_force;
  if (disp_old_aug)    delete [] disp_old_aug;
  if (SAVE_CVARS){
    delete [] conface;

    delete [] normal_force_mag;
    delete [] normal_traction_mag;
    delete [] tangential_force_mag;
    delete [] tangential_traction_mag;

    delete [] normal_dir_x;
    delete [] normal_dir_y;
    delete [] normal_dir_z;

    delete [] tangential_dir_x;
    delete [] tangential_dir_y;
    delete [] tangential_dir_z;

    delete [] slipmag;
    delete [] nodal_dissipation;
    delete [] nodal_dissipation_density;
    delete [] contact_area;
    delete [] plot_gap_cur;
    delete [] plot_gap_old;
    delete [] kinematic_partition;
  }
}


ContactSearch::ContactErrorCode 
ContactTDEnforcement::Get_Plot_Variable( Contact_TDEnf_Plot_Vars var,
					       Real * buffer )
{
  int number_of_nodes = number_host_code_nodes; // PJSA size of buffer is number_host_code_nodes (host doesn't know about extra nodes added for shell lofting)
  if( !SAVE_CVARS ){
    errors->Add_Error_Message( "CVARS weren't requested" );
    return ContactSearch::INVALID_DATA;
  }
  switch( var ){
  case( CONFACE ):
    std::memcpy( buffer, conface, number_of_nodes*sizeof(Real) );
    break;
  case( NORMAL_FORCE_MAG ):
    std::memcpy( buffer, normal_force_mag, number_of_nodes*sizeof(Real) );
    break;
  case( NORMAL_TRACTION_MAG ):
    std::memcpy( buffer, normal_traction_mag, number_of_nodes*sizeof(Real) );
    break;
  case( TANGENTIAL_FORCE_MAG ):
    std::memcpy( buffer, tangential_force_mag, number_of_nodes*sizeof(Real) );
    break;
  case( TANGENTIAL_TRACTION_MAG ):
    std::memcpy( buffer, tangential_traction_mag, number_of_nodes*sizeof(Real) );
    break;
  case( CDIRNORX ):
    std::memcpy( buffer, normal_dir_x, number_of_nodes*sizeof(Real) );
    break;
  case( CDIRNORY ):
    std::memcpy( buffer, normal_dir_y, number_of_nodes*sizeof(Real) );
    break;
  case( CDIRNORZ ):
    std::memcpy( buffer, normal_dir_z, number_of_nodes*sizeof(Real) );
    break;
  case( CDIRTANX ):
    std::memcpy( buffer, tangential_dir_x, number_of_nodes*sizeof(Real) );
    break;
  case( CDIRTANY ):
    std::memcpy( buffer, tangential_dir_y, number_of_nodes*sizeof(Real) );
    break;
  case( CDIRTANZ ):
    std::memcpy( buffer, tangential_dir_z, number_of_nodes*sizeof(Real) );
    break;
  case( SLIP_MAG ):
    std::memcpy( buffer, slipmag, number_of_nodes*sizeof(Real) );
    break;
  case( NODAL_DISSIPATION ):
    std::memcpy( buffer, nodal_dissipation, number_of_nodes*sizeof(Real) );
    break;
  case( NODAL_DISSIPATION_DENSITY ):
    std::memcpy( buffer, nodal_dissipation_density, number_of_nodes*sizeof(Real) );
    break;
  case (CONTACT_AREA):
    std::memcpy(buffer, contact_area, number_of_nodes*sizeof(Real));
    break;
  case (GAP_CUR):
    std::memcpy(buffer, plot_gap_cur, number_of_nodes*sizeof(Real));
    break;
  case (GAP_OLD):
    std::memcpy(buffer, plot_gap_old, number_of_nodes*sizeof(Real));
    break;
  case (KINEMATIC_PARTITION_VALUE):
    std::memcpy(buffer, kinematic_partition, number_of_nodes*sizeof(Real));
    break;
  };
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode 
ContactTDEnforcement::Get_Global_Variable( Contact_TDEnf_Global_Vars var,
					       Real& value )
{
  if( !SAVE_CGVARS ){
    errors->Add_Error_Message( "CGVARS weren't requested" );
    return ContactSearch::INVALID_DATA;
  }
  switch( var ){
  case( FORCE_X ):
    value = global_force_x;
    break;
  case( FORCE_Y ):
    value = global_force_y;
    break;
  case( FORCE_Z ):
    value = global_force_z;
    break;
  case( FORCE_NORM ):
    value = global_force_norm;
    break;
  case( DISSIPATION ):
    value = global_dissipation;
    break;
  case( CONSTRAINT_NORM ):
    value = global_constraint_norm;
    break;
  case( INC_FORCE_NORM ):
    value = global_inc_force_norm;
    break;
  };
  return ContactSearch::NO_ERROR;
}

int ContactTDEnforcement::Restart_Size()
{
  int size = 0;

  // Get base class restart size first
  size = ContactEnforcement::Restart_Size();

  // Now add size for additional td_enforcement info.
  size += 1; //time multiplier

  size += 1; // enforcement iterations

  size += 1; // number of old augmented displacements

  size += 3*number_of_nodes; // displacements

  size += 3; // symmetry nodes information

  return size;
}

ContactSearch::ContactErrorCode
ContactTDEnforcement::Extract_Restart_Data( Real* restart_data )
{
  ContactSearch::ContactErrorCode ec;
  // add base class restart data first
  ec = ContactEnforcement::Extract_Restart_Data( restart_data );
  int words_added = ContactEnforcement::Restart_Size(); // base class size
  // now add td enforcement specific data to restart
  restart_data[words_added++] = time_mult_old_aug; // time step multiplier
  restart_data[words_added++] = num_iterations;    // num of enforcement iters.
  restart_data[words_added++] = 3*number_of_nodes; // number of displacements
  for(int i=0; i<3*number_of_nodes; ++i) {
    restart_data[words_added++] = disp_old_aug[i];
  }
  restart_data[words_added++] = have_symmetry_nodes;
  restart_data[words_added++] = symmetry_node_1_exodus_id;
  restart_data[words_added++] = symmetry_node_2_exodus_id;

  POSTCONDITION( words_added == Restart_Size() );
  return ec;
}
  
ContactSearch::ContactErrorCode 
ContactTDEnforcement::Set_Number_of_Iterations( int iter )
{
  if( iter > 0 ){
    num_iterations = iter;
    return ContactSearch::NO_ERROR;
  }
  errors->Add_Error_Message( "Number of iterations must be positive" );
  return ContactSearch::INVALID_DATA;
}

ContactSearch::ContactErrorCode 
ContactTDEnforcement::Set_Convergence_Tolerance( Real tol )
{
  if( tol > 0.0 ){
    convergence_tolerance = tol;
    return ContactSearch::NO_ERROR;
  }
  errors->Add_Error_Message( "Convergence tolerance must be positive" );
  return ContactSearch::INVALID_DATA;
}

ContactSearch::ContactErrorCode
ContactTDEnforcement::Enforce_Symmetry_on_Nodes( int n1, int n2 )
{
  if( have_symmetry_nodes ){
    errors->Add_Error_Message( "Can only specify one pair of symmetry nodes" );
    return ContactSearch::INVALID_DATA;
  }
  if( n1 <= 0  || n2 <= 0 ){
    errors->Add_Error_Message( "Enforce_Symmetry_on_Nodes: ids must be positive" );
    return ContactSearch::INVALID_DATA;
  }
  // Verify someone has this node even if I don't
  int have_node_n1 = 0;
  int num_nodes = topology->Number_of_Nodes();
  ContactNode<Real>** Nodes = reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());
  for (int i=0; i<num_nodes; ++i) {
    if( Nodes[i]->Exodus_ID() == n1 ){
      have_node_n1 = true;
      break;
    }
  }

  int have_node_n2 = 0;
  for (int i=0; i<num_nodes; ++i) {
    if( Nodes[i]->Exodus_ID() == n2 ){
      have_node_n2 = true;
      break;
    }
  }

  int sum_input[2], sum_output[2];
  sum_input[0] = have_node_n1;
  sum_input[1] = have_node_n2;

  contact_global_sum( sum_input, sum_output, 2, communicator );

  have_node_n1 = sum_output[0];
  have_node_n2 = sum_output[1];

  if( have_node_n1 == 0 ){
    std::sprintf( message, "Node %d is not part of the contact surface", n1 );
    return ContactSearch::INVALID_DATA;
  }
  symmetry_node_1_exodus_id = n1;


  if( have_node_n2 == 0 ){
    symmetry_node_1_exodus_id = -1;
    std::sprintf( message, "Node %d is not part of the contact surface", n2 );
    return ContactSearch::INVALID_DATA;
  }
  symmetry_node_2_exodus_id = n2;
  have_symmetry_nodes = true;

  return ContactSearch::NO_ERROR;
}

Real ContactTDEnforcement::Get_Global_Plot_Variable( int var_num )
{
  PRECONDITION( var_num < number_global_plot_vars );
  if( var_num == 0 )
    return dt_old;
  else
    return dt;
}

void ContactTDEnforcement::Get_Nodal_Plot_Variable( int var_num, Real* data )
{
  PRECONDITION( var_num < number_nodal_plot_vars );
  PRECONDITION(CALC_PLOT_FORCE);
  std::memcpy( data, plot_force+var_num*number_of_nodes, 
               number_of_nodes*sizeof(Real) );
}

void ContactTDEnforcement::Get_Old_Displacement_Ptrs( int var_num, Real** data )
{
  PRECONDITION( var_num < 3 );
  *data = disp_old_aug+var_num*number_of_nodes;
}

void ContactTDEnforcement::Get_Old_Displacement( int var_num, Real* data )
{
  PRECONDITION( var_num < 3 );
  std::memcpy( data, disp_old_aug+var_num*number_of_nodes, 
               number_of_nodes*sizeof(Real) );
}

void ContactTDEnforcement::Get_Old_Disp_Restart( int var_num, Real* data )
{
  if (search->Get_Primary_Topology()->Have_Shells()) {
    ContactShellHandler * shell_handler = 
      search->Get_Primary_Topology()->Shell_Handler();
    int host_to_acme_node = 0;
    int var_idx = var_num;
    while (var_idx >= 3){
      var_idx -= 3;
      ++host_to_acme_node;
    }
 
    for ( int i = 0; i < shell_handler->Number_Host_Code_Nodes(); ++i) {
      if (host_to_acme_node < shell_handler->Num_Acme_Nodes_for_Host_Node(i)){
	int node_idx = 
	  shell_handler->Acme_Node_for_Host_Node(i,host_to_acme_node);
	data[i] = disp_old_aug[var_idx*number_of_nodes + node_idx];
      }
      else data[i] = 0.0;
    }
  }
  else Get_Old_Displacement(var_num, data);
}

void ContactTDEnforcement::Set_Old_Disp_Restart( int var_num, Real* data )
{
  if (search->Get_Primary_Topology()->Have_Shells()) {
    ContactShellHandler * shell_handler = 
      search->Get_Primary_Topology()->Shell_Handler();
    int host_to_acme_node = 0;
    int var_idx = var_num;
    while (var_idx >= 3){
      var_idx -= 3;
      ++host_to_acme_node;
    }
    
    for ( int i = 0; i < shell_handler->Number_Host_Code_Nodes(); ++i) {
      if (host_to_acme_node < shell_handler->Num_Acme_Nodes_for_Host_Node(i)){
	int node_idx = 
	  shell_handler->Acme_Node_for_Host_Node(i,host_to_acme_node);
	disp_old_aug[var_idx*number_of_nodes + node_idx] = data[i];
      }
      else continue;
    }
    
  }
  else {
    PRECONDITION( var_num < 3 );
    std::memcpy( disp_old_aug+var_num*number_of_nodes, data,
	    number_of_nodes*sizeof(Real) );
  }
}

ContactSearch::ContactErrorCode
ContactTDEnforcement::Extract_General_Restart_Variable( Real* data )
{
  Real* buffer = data;
  Extract_Base_General_Restart_Variable( buffer );
  buffer += Number_Base_General_Restart_Variables();
  buffer[0] = time_mult_old_aug;
  
  return ContactSearch::NO_ERROR;           
}

int
ContactTDEnforcement::Number_Nodal_Restart_Variables()
{
  int num_vars = 3; // three terms in old displacement

  // now loop over the enforcement models and count up their state data
  for (int i = 0; i < num_enforcement_models; ++i)
    num_vars += enforcement_models[i]->Num_Node_State_Variables();

  if (search->Get_Primary_Topology()->Have_Shells())
    num_vars *= 
      search->Get_Primary_Topology()->Shell_Handler()->
      Max_Num_Acme_Nodes_for_Host_Node();

  return num_vars;
}

ContactSearch::ContactErrorCode
ContactTDEnforcement::Extract_Nodal_Restart_Variable( int n, Real* data, int* node_id )
{
  int num_restart_vars = Number_Nodal_Restart_Variables();
  int max_host_to_acme_size;
  bool have_shells =  search->Get_Primary_Topology()->Have_Shells();
  if (have_shells) 
    max_host_to_acme_size = 
      search->Get_Primary_Topology()->Shell_Handler()->
      Max_Num_Acme_Nodes_for_Host_Node();
  else max_host_to_acme_size = 1;
  int num_restart_vars_per_node = 
    Number_Nodal_Restart_Variables()/max_host_to_acme_size;
  POSTCONDITION ( max_host_to_acme_size * num_restart_vars_per_node 
		  == Number_Nodal_Restart_Variables());
  if( n>num_restart_vars || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  
  // get node ids to return to host code
  int nnod;
  if ( have_shells )
    nnod = topology->Shell_Handler()->Number_Host_Code_Nodes();
  else
    nnod = topology->Number_of_Nodes();
  for( int i=0 ; i<nnod ; ++i){
    if (have_shells) {
      int hi, lo;
      topology->Shell_Handler()->Acme_NodeGID_for_Host_Node(i,0,hi,lo);
      ContactHostGlobalID first_gid(hi,lo);
      ContactNode<Real> * node = static_cast<ContactNode<Real>*>
        (topology->NodeList()->Find(first_gid));
      ContactShellNode * shell_node = dynamic_cast<ContactShellNode*>(node);
      if (NULL != shell_node) {
        node_id[i]=shell_node->Shell_Node_Base_ID();
      } else {
        node_id[i]=node->Exodus_ID();
      }
    } else {
      ContactNode<Real> * node =
        static_cast<ContactNode<Real>*>(topology->NodeList()->Find(i));
      node_id[i]=node->Exodus_ID();
    }
  }

  // figure which shell node on the host code node we are processing
  // (this is always zero if there are no shells in the analysis)
  int host_to_acme_node = 0;
  int var_num = n;
  while (var_num > num_restart_vars_per_node){
    var_num -= num_restart_vars_per_node;
    ++host_to_acme_node;
  }
  
  PRECONDITION ( have_shells || host_to_acme_node == 0 );

  // copy displacements into data
  if ( var_num <= 3 ) {
    int var_disp = var_num - 1 + host_to_acme_node*3;
    Get_Old_Disp_Restart(var_disp,data);
  }
  else {
    int this_enf_var_num = var_num - 3;
    int enf_var_count(0);
    for (int i = 0; i < num_enforcement_models; ++i) {
      int num_enf_model_vars = 
	enforcement_models[i]->Num_Node_State_Variables();
      if (num_enf_model_vars == 0 ) continue;
      if (enf_var_count+num_enf_model_vars >= this_enf_var_num ) {
	int local_enf_var = this_enf_var_num - enf_var_count;
	PRECONDITION( local_enf_var > 0 && local_enf_var < num_enf_model_vars);
	local_enf_var += host_to_acme_node * num_enf_model_vars;
	ContactSearch::ContactErrorCode err = 
	  enforcement_models[i]->Extract_Nodal_Restart_Variable(local_enf_var,
								data);
	if ( err != ContactSearch::NO_ERROR ) return err;
	break;
      }
      enf_var_count += num_enf_model_vars;
    }
  }
  
  return ContactSearch::NO_ERROR;           
}

ContactSearch::ContactErrorCode
ContactTDEnforcement::Implant_Nodal_Restart_Variable( int n, Real* data )
{
  int num_restart_vars = Number_Nodal_Restart_Variables();
  int max_host_to_acme_size;
  if (search->Get_Primary_Topology()->Have_Shells()) 
    max_host_to_acme_size = 
      search->Get_Primary_Topology()->Shell_Handler()->
      Max_Num_Acme_Nodes_for_Host_Node();
  else max_host_to_acme_size = 1;
  int num_restart_vars_per_node = 
    Number_Nodal_Restart_Variables()/max_host_to_acme_size;
  POSTCONDITION ( max_host_to_acme_size * num_restart_vars_per_node 
		  == Number_Nodal_Restart_Variables());
  if( n>num_restart_vars  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  
  int host_to_acme_node = 0;
  int var_num = n;
  while (var_num > num_restart_vars_per_node){
    var_num -= num_restart_vars_per_node;
    ++host_to_acme_node;
  }
  
  PRECONDITION ( search->Get_Primary_Topology()->Have_Shells() || 
		 host_to_acme_node == 0 );
  
  // copy displacements into data
  if ( var_num <= 3 ) {
    int var_disp = var_num - 1 + host_to_acme_node*3;
    Set_Old_Disp_Restart(var_disp,data);
  }
  else {
    int this_enf_var_num = n - 3;
    int enf_var_count(0);
    for (int i = 0; i < num_enforcement_models; ++i) {
      int num_enf_model_vars = 
	enforcement_models[i]->Num_Node_State_Variables();
      if (num_enf_model_vars == 0 ) continue;
      if (enf_var_count+num_enf_model_vars >= this_enf_var_num ) {
	int local_enf_var = this_enf_var_num - enf_var_count;
	PRECONDITION( local_enf_var > 0 && local_enf_var < num_enf_model_vars);
	local_enf_var += host_to_acme_node * num_enf_model_vars;
	ContactSearch::ContactErrorCode err = 
	  enforcement_models[i]->Implant_Nodal_Restart_Variable(local_enf_var,
								data);
	if ( err != ContactSearch::NO_ERROR ) return err;
	break;
      }
      enf_var_count += num_enf_model_vars;
    }
  }

  // if all the enforcement model data has been set up, then initialize
  // the enforcement models
  if ( n == Number_Nodal_Restart_Variables() ) {
    for( int i=0 ; i<num_enforcement_models ; ++i ){
      enforcement_models[i]->Initialize_Model( num_enforcement_models,
					       enforcement_models,
					       search->Num_Tables(),
					       search->Tables() );
    }
  }

  return ContactSearch::NO_ERROR;       
}

ContactSearch::ContactErrorCode
ContactTDEnforcement::Extract_Edge_Restart_Variable( int n, Real* data )
{
  if( n>Number_Edge_Restart_Variables() || n<=0  ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR; 
}

ContactSearch::ContactErrorCode
ContactTDEnforcement::Implant_General_Restart_Variable( Real* data )
{
  Real* buffer = data;
  Implant_Base_General_Restart_Variable( buffer );
  buffer += Number_Base_General_Restart_Variables();
  time_mult_old_aug = buffer[0];
  
  return ContactSearch::NO_ERROR;           
}

ContactSearch::ContactErrorCode
ContactTDEnforcement::Implant_Edge_Restart_Variable( int n, Real* data )
{
  if( n>Number_Edge_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;   
}

ContactSearch::ContactErrorCode
ContactTDEnforcement::Extract_Face_Restart_Variable( int n, Real* data )
{
  if( n>Number_Face_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  PRECONDITION(CALC_PLOT_FORCE);
  std::memcpy( data, plot_force+(n-1), number_of_nodes*sizeof(Real) );
  return ContactSearch::NO_ERROR;         
}

ContactSearch::ContactErrorCode
ContactTDEnforcement::Implant_Face_Restart_Variable( int n, Real* data )
{
  if( n>Number_Face_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  PRECONDITION(CALC_PLOT_FORCE);
  std::memcpy( plot_force+(n-1), data, number_of_nodes*sizeof(Real) );
  return ContactSearch::NO_ERROR;      
}

ContactSearch::ContactErrorCode
ContactTDEnforcement::Extract_Element_Restart_Variable( int n, Real* data )
{
  if( n>Number_Element_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;         
}

ContactSearch::ContactErrorCode
ContactTDEnforcement::Implant_Element_Restart_Variable( int n, Real* data )
{
  if( n>Number_Element_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;      
}


ContactSearch::ContactErrorCode 
ContactTDEnforcement::Compute_Contact_Force( Real DT_old, Real DT,
					     Real* Mass,
					     Real* Density,
					     Real* Wavespeed,
					     Real* Force  )
{
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
  Timer().Start_Timer( search->Contact_time() );
  Timer().Start_Timer( enforcement_time );
  Timer().Start_Timer( initialization_time );
#endif

#if defined (CONTACT_DEBUG_NODE) || (LOCAL_PRINT_FLAG>0) || (CONTACT_DEBUG_PRINT_LEVEL>=2)
  ContactParOStream& postream = ParOStream();
#endif

#if CONTACT_DEBUG_PRINT_LEVEL>=2
  postream<<"TD Enforcement\n";
  postream<<"  max iterations    = "<<num_iterations<<"\n";
  if ( convergence_tolerance > 0.0) {
    postream<<"  tolerence         = "<<convergence_tolerance<<"\n";
  }
#endif
#ifdef CONTACT_HEARTBEAT 
  if (contact_processor_number(communicator)==0) {
    std::cout<<"TD Enforcement\n";
    std::cout<<"  max iterations    = "<<num_iterations<<"\n";
    if ( convergence_tolerance > 0.0) {
      std::cout<<"  tolerence         = "<<convergence_tolerance<<"\n";
    }
  }
#endif

  ContactSearch::ContactErrorCode error_code = ContactSearch::NO_ERROR;

  dt = DT;
  dt_old = DT_old;
  dt2 = 1.0/(0.5*(dt+dt_old)*dt);

#ifdef CONTACT_DEBUG_NODE
  if( Number_Debug_Nodes() ){
    if( contact_processor_number(communicator) == 0 ){
      postream << "dt = " << dt << "  dt_old = " << dt_old << "\n";
    }
  }  
  for(int i=0 ; i<Number_Debug_Nodes() ; ++i ){
    ContactNode<Real>* dn = Debug_Node( i );
    if( dn ){
      int index = Get_Node_Host_Array_Index(dn);
      postream << "TD Enforcement (m, rho, c) for Node "
               << dn->Exodus_ID() << " ("
               << Mass[index] << " , "
               << Density[index] << " , "
               << Wavespeed[index] << ")\n";
    }
  }
  postream.flush();
#endif

  int have_interactions = contact_global_sum(topology->Number_NodeEntity_Interactions(),communicator);
  if (have_interactions==0) {
    std::memset( Force,0, number_host_code_nodes*sizeof(Real)*NDIM );
    if(SAVE_CVARS) Initialize_CVARS();
#ifdef CONTACT_HEARTBEAT 
    if (contact_processor_number(communicator)==0) {
      std::cout << "End of enforcement\n"<<flush;
    }
#endif
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Stop_Timer( initialization_time );
    Timer().Stop_Timer( enforcement_time );
    Timer().Stop_Timer( search->Contact_time() );
#endif
    return error_code;
  }


  // Call base class Set_Up function to prepare for enforcement
  error_code = Set_Up();

  error_code = (ContactSearch::ContactErrorCode) contact_global_error_check( error_code, communicator );

  if( error_code ){
    Clean_Up();
    return error_code;
  }

  //
  //  Determine if automatic kinematic partionals are every used
  //
  have_auto_kin_part = 0;
  for( int iface_entity=0 ; iface_entity<number_entity_keys ; ++iface_entity ){
    for( int inode_entity=0 ; inode_entity<number_entity_keys ; ++inode_entity ){
      Real kin_part = Enforcement_Data( KINEMATIC_PARTITION, iface_entity, inode_entity );
      if( kin_part+0.25 >= 2.0 && kin_part+0.25 < 3.0) {
        have_auto_kin_part = 1;
      }
    }
  }
  
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
  Timer().Start_Timer( scratch_allocation_time );
#endif
  
  // Get the scratch memory 
  NODAL_MASS.       Allocate_Scratch(number_of_total_nodes, 1, ScratchVariable::ZERO_SCRATCH);
  SN_ACCELERATION.  Allocate_Scratch(number_of_total_nodes, 3, ScratchVariable::ZERO_SCRATCH);
  NEW_POSITION.     Allocate_Scratch(number_of_total_nodes, 3, ScratchVariable::ZERO_SCRATCH);
  TOTAL_FORCE.      Allocate_Scratch(number_of_total_nodes, 3, ScratchVariable::ZERO_SCRATCH);

  NEI_TOTAL_SLIP.         Allocate_Scratch(number_node_entity_constraints, 3, ScratchVariable::ZERO_SCRATCH);
  NEI_REL_DISP.           Allocate_Scratch(number_node_entity_constraints, 3, ScratchVariable::ZERO_SCRATCH);
  NEI_FORCE.              Allocate_Scratch(number_node_entity_constraints, 3, ScratchVariable::ZERO_SCRATCH);
  NEI_TAN_SLIP.           Allocate_Scratch(number_node_entity_constraints, 3, ScratchVariable::ZERO_SCRATCH);

  NEI_KINEMATIC_PARTITION.Allocate_Scratch(number_node_entity_constraints, 1);
  NEI_GAP_TO_ENFORCE.     Allocate_Scratch(number_node_entity_constraints, 1);
  NEI_TARGET_GAP.         Allocate_Scratch(number_node_entity_constraints, 1);
  //
  //  Various optional scratch memory
  //
  if(have_auto_kin_part || plot_force) {
    NODAL_DENSITY.    Allocate_Scratch(number_of_total_nodes, 1, ScratchVariable::ZERO_SCRATCH);
    NODAL_WAVESPEED.  Allocate_Scratch(number_of_total_nodes, 1, ScratchVariable::ZERO_SCRATCH);
  }


  PREDICTED_POS.    Allocate_Scratch(number_of_total_nodes, 3);
  CURRENT_POS.      Allocate_Scratch(number_of_total_nodes, 3);
  KIN_CONSTR_VECTOR.Allocate_Scratch(number_of_total_nodes, 3);
  NUM_KIN_CONSTR.   Allocate_Scratch(number_of_total_nodes, 1);

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
  Timer().Stop_Timer( scratch_allocation_time );
  Timer().Start_Timer( copy_nodal_to_scratch );
#endif
  //
  //  Copy permanante nodal variables to scratch variables
  //
  {
    VariableHandle var_handles[4];
    var_handles[0] = topology->Variable_Handle( ContactTopology::Current_Position );
    var_handles[1] = topology->Variable_Handle( ContactTopology::Predicted_Position );
    var_handles[2] = topology->Variable_Handle( ContactTopology::Kin_Constr_Vector );
    var_handles[3] = topology->Variable_Handle( ContactTopology::Num_Kin_Constr );
    ScratchVariable* scratch_arrays[4];  
    scratch_arrays[0] = &CURRENT_POS;
    scratch_arrays[1] = &PREDICTED_POS;
    scratch_arrays[2] = &KIN_CONSTR_VECTOR;
    scratch_arrays[3] = &NUM_KIN_CONSTR;
    Copy_Variable_to_Scratch(4, var_handles, scratch_arrays);
  }
  //
  //  Set initial updated position to the current predicted positiion
  //
  NEW_POSITION.Duplicate_Scratch(PREDICTED_POS);

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
  Timer().Stop_Timer( copy_nodal_to_scratch );
  Timer().Start_Timer( init_host_code_vars );
#endif
  //
  // Shared node-face, node-entity scratch data
  //
  // Initialize the scratch memory 
  //
  std::memset( Force,0, number_host_code_nodes*sizeof(Real)*NDIM );
    
  if(have_auto_kin_part || plot_force) {
    ScratchVariable* scratch_arrays[3];
    scratch_arrays[0] = &NODAL_MASS;
    scratch_arrays[1] = &NODAL_DENSITY;
    scratch_arrays[2] = &NODAL_WAVESPEED;
    Real* host_arrays[3];
    host_arrays[0] = Mass;
    host_arrays[1] = Density;
    host_arrays[2] = Wavespeed;
    Copy_Host_Scalar_Arrays_to_Node_Scratch( 3, host_arrays, scratch_arrays );
#ifndef CONTACT_NO_MPI
    if(contact_number_of_processors( communicator ) > 1 ){
      contact_import_scratch_vars(communicator, 
                                  *Node_AsymComm, 
                                  *(search->Get_Comm_Buffer()),
                                  3,
                                  scratch_arrays);
    }
#endif
  } else {
    ScratchVariable* scratch_arrays[1];
    scratch_arrays[0] = &NODAL_MASS;
    Real* host_arrays[1];
    host_arrays[0] = Mass;
    Copy_Host_Scalar_Arrays_to_Node_Scratch( 1, host_arrays, scratch_arrays );
#ifndef CONTACT_NO_MPI
    if(contact_number_of_processors( communicator ) > 1 ){
      contact_import_scratch_var(communicator, 
                                 *Node_AsymComm, 
                                  *(search->Get_Comm_Buffer()),
                                  NODAL_MASS);
    }
#endif
  }

#ifdef CONTACT_DEBUG_NODE
    for(int cdn_i=0 ; cdn_i<Number_Debug_Nodes() ; cdn_i++ ){
      ContactNode<Real>* dn = Debug_Node( cdn_i );
      if( dn ){
        Real* pred_coords = NEW_POSITION.Get_Scratch(dn);
        postream << "New Predicted Configuration for Node "
	       << dn->Exodus_ID() << " = " 
	       << pred_coords[0] << " " << pred_coords[1] << " " 
	       << pred_coords[2] << "\n";
      }
    }
    postream.flush();
#endif

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
  Timer().Stop_Timer( init_host_code_vars );
  Timer().Start_Timer( td_enforcement_setup );
#endif


  error_code = TD_Enforcement_Set_Up();


#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
  Timer().Stop_Timer( td_enforcement_setup );
  Timer().Start_Timer( setup_interaction_list );
#endif
  //
  //  Compact out duplicate or invalid interactions
  //
  Compact_Node_Entity_List();
  //
  //  Set the fast lookup methods for the entity list
  //    
  node_entity_list.Finalize();
  
#ifdef CONTACT_HEARTBEAT 
  int cnt_global = contact_global_sum(node_entity_list.Size(),communicator);
  if (contact_processor_number(communicator)==0) {
    std::cout << "  # of interactions = " << cnt_global << "\n";
    std::cout << "  have_multiple     = " << have_multiple <<"\n";
  }
#endif

  error_code = (ContactSearch::ContactErrorCode) contact_global_error_check( error_code, communicator );

  if(  error_code ){
    TD_Enforcement_Clean_Up();
    Clean_Up();
    return error_code;
  }

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
  Timer().Stop_Timer( setup_interaction_list );
  Timer().Stop_Timer( initialization_time );
  Timer().Start_Timer( iteration_time );
#endif

  //
  //  Create iteration scratch variables
  //
  ASSEMBLED_MASS.   Allocate_Scratch(number_of_total_nodes, 1);
  ASSEMBLED_FORCE.  Allocate_Scratch(number_of_total_nodes, 3);
  INC_FORCE.        Allocate_Scratch(number_of_total_nodes, 3);
  A_CORRECTION.     Allocate_Scratch(number_of_total_nodes, 3);
  A_CONSTITUTIVE.   Allocate_Scratch(number_of_total_nodes, 3);
  A_PREDICTION.     Allocate_Scratch(number_of_total_nodes, 3);
  //
  //  Calculate the initial assembled masses (face mass)
  //
  if(!have_multiple) {
    //cout<<"assembling masses"<<endl;
    Assemble_Masses();
  }

  Real mult = 0.0;

  int iteration;
  for(iteration=0 ; iteration<num_iterations ; ++iteration ){
    mult = mult + 1.0;

    if (contactPrintFlag > 1 && contact_processor_number(communicator)==0) {
      std::cout << "  Iteration "<<iteration<<"\n";
    }

    //
    // Zero the scratch memory
    //
    if(have_multiple) {
      //cout<<"zero mass"<<endl;

      ASSEMBLED_MASS. Zero_Scratch();
    }
    ASSEMBLED_FORCE.Zero_Scratch();
    INC_FORCE      .Zero_Scratch();
    A_PREDICTION   .Zero_Scratch();  // MWG: necessary???
    A_CORRECTION   .Zero_Scratch();
    A_CONSTITUTIVE .Zero_Scratch();

    // Compute the gap to enforce and reset the pushback directoin
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Start_Timer( compute_gap_time );
#endif

    Compute_Kinematic_Quantities( iteration );
    
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Stop_Timer( compute_gap_time );
#endif

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Start_Timer( penalty_time );
#endif

    // Compute penalty forces
    Compute_Penalty_Forces();
    
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Stop_Timer( penalty_time );
#endif

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Start_Timer( kinematic_consistent_time );
#endif

    Make_Scratch_Vector_Consistent_with_BCs( SN_ACCELERATION );

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Stop_Timer( kinematic_consistent_time );
#endif

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Start_Timer( assemble_forces_time );
#endif
    //
    // Assemble nodal forces
    // SWAPADD the assembled_force and if necessary the assembled_mass
    //
    if(have_multiple){
      //cout<<"assemble mult"<<endl;
      Assemble_Nodal_Forces_and_Masses();
      ScratchVariable *vars[2];
      vars[0] = &ASSEMBLED_FORCE;
      vars[1] = &ASSEMBLED_MASS;
      swapadd_node_scratch_vars( 2, vars );
    } else {
      //cout<<"assemble no mult"<<endl;
      Assemble_Nodal_Forces_No_Multiple();
      swapadd_node_scratch(ASSEMBLED_FORCE);
    }

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Stop_Timer( assemble_forces_time );
#endif

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Start_Timer( kinematic_consistent_time );
#endif

    Make_Scratch_Vector_Consistent_with_BCs( ASSEMBLED_FORCE );

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Stop_Timer( kinematic_consistent_time );
#endif

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Start_Timer( accel_cor_time );
    Timer().Start_Timer( response_prediction_time );
#endif

    Response_Prediction(mult);

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Stop_Timer( response_prediction_time );
    Timer().Start_Timer( response_correction_time );
#endif

    Response_Correction();

#ifndef CONTACT_NO_MPI

#if defined (CONTACT_TIMINGS)
    Timer().Stop_Timer( response_correction_time );
    Timer().Start_Timer( swap_add_time );
#endif
    //
    // Import the node acceleration correction
    //  Note, as constraints are owned by the processor that owns the constraint node,
    //  only an import is required here, not a swap add. 
    //
    contact_import_scratch_var(communicator, 
                               *Node_AsymComm_Full, 
                               *(search->Get_Comm_Buffer()),
                               A_CORRECTION);
    
#if defined (CONTACT_TIMINGS)
    Timer().Stop_Timer( swap_add_time );
    Timer().Start_Timer( contact_force_time );
#endif

#endif

    //
    // Compute the contact force (by nodes)
    //
    Real *a_pre = A_PREDICTION.Get_Data();
    Real *a_cor = A_CORRECTION.Get_Data();
    Real *m_sn  = NODAL_MASS.Get_Data();
    Real *f_inc = INC_FORCE.Get_Data();    

    for(int inode=0 ; inode<number_of_total_nodes ; ++inode ){
      f_inc[0] = m_sn[0]*(a_pre[0]+a_cor[0]);
      f_inc[1] = m_sn[0]*(a_pre[1]+a_cor[1]);
      f_inc[2] = m_sn[0]*(a_pre[2]+a_cor[2]);
      f_inc += 3;
      a_pre += 3;
      a_cor += 3;
      m_sn += 1;
    }    

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Stop_Timer( contact_force_time );
#endif


    if (! only_frictionless_or_tied) {
      //
      // correct for constitutive behavior
      //
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
      Timer().Start_Timer( response_constitutive_correction_time );
#endif

      Response_Constitutive_Correction();

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
      Timer().Stop_Timer( response_constitutive_correction_time );
      Timer().Start_Timer( swap_add_time );
#endif

      // SWAPADD the node acceleration correction
      swapadd_node_scratch( A_CONSTITUTIVE );

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
      Timer().Stop_Timer( swap_add_time );
      Timer().Start_Timer( force_adjust_time );
#endif

      f_inc = INC_FORCE.Get_Data();    
      m_sn  = NODAL_MASS.Get_Data();
      Real *a_con = A_CONSTITUTIVE.Get_Data();
      for(int i=0 ; i<number_of_total_nodes ; ++i ){
        f_inc[0] += m_sn[0]*a_con[0];
        f_inc[1] += m_sn[0]*a_con[1];
        f_inc[2] += m_sn[0]*a_con[2];
        f_inc += 3;
        a_con += 3;
        m_sn  += 1;
      }

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
      Timer().Stop_Timer( force_adjust_time );
#endif
    }

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Stop_Timer( accel_cor_time );
#endif

    if( topology->Have_Shells() ) Assemble_Shell_Forces();

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Start_Timer( kinematic_consistent_time );
#endif
    // make force compatible with prescibed boundary conditions

    Make_Scratch_Vector_Consistent_with_BCs( INC_FORCE );

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Stop_Timer( kinematic_consistent_time );
#endif

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Start_Timer( update_config_time );
#endif
    //
    // Compute the new coordinates this implies
    //
    Real *x_new = NEW_POSITION.Get_Data();
    f_inc = INC_FORCE.Get_Data();
    Real *f_tot = TOTAL_FORCE.Get_Data();
    m_sn  = NODAL_MASS.Get_Data();

    for(int i=0 ; i<number_of_total_nodes ; ++i ){
      Real  scale = 1.0/(m_sn[0]*dt2);
      x_new[0] += f_inc[0]*scale;
      x_new[1] += f_inc[1]*scale;
      x_new[2] += f_inc[2]*scale;
      f_tot[0] += f_inc[0];
      f_tot[1] += f_inc[1];
      f_tot[2] += f_inc[2];
      x_new += 3;
      f_inc += 3;
      f_tot += 3;
      m_sn  += 1;
    }

    // We need to get the predicted positions of the phantom nodes

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Stop_Timer( update_config_time );
#endif


#ifdef CONTACT_DEBUG_NODE
    for(int i=0 ; i<Number_Debug_Nodes() ; ++i ){
      ContactNode<Real>* dn = Debug_Node( i );
      if( dn ){
        const int node_index = dn->EnfArrayIndex();
	Real m_sn_dbg   = *(NODAL_MASS.  Get_Scratch(node_index));
	Real* inc_force =  INC_FORCE.    Get_Scratch(node_index);
	Real* pcf       =  NEW_POSITION. Get_Scratch(node_index);
	Real* pci       =  PREDICTED_POS.Get_Scratch(node_index);
	Real Mult = m_sn_dbg*dt2;
	postream << "Iterative Force  for Node " << dn->Exodus_ID() << " = "
		 << inc_force[0] << " " << inc_force[1] << " " << inc_force[2] 
		 << ", for iteration : " << iteration << "\n";
	postream << "Total Force     for Node " << dn->Exodus_ID() << " = "
		 << Mult*(pcf[0]-pci[0]) << " "
		 << Mult*(pcf[1]-pci[1]) << " "
		 << Mult*(pcf[2]-pci[2]) 
		 << ", for iteration : " << iteration << "\n";
      }
    }
    postream.flush();
#endif

    // calculate norm of the incremental force (multiplier)
    if ( convergence_tolerance > 0.0) {
      Real inc_force_norm = 0.0;
      
      f_inc = INC_FORCE.Get_Data();
      ContactTopologyEntity<Real> **e = topology->NodeList()->EntityList();
      for(int i=0 ; i<number_of_nodes ; ++i ){
	if ( e[i]->Ownership() == ContactTopologyEntity<Real>::OWNED ) {
	  inc_force_norm += f_inc[0]*f_inc[0];
	  inc_force_norm += f_inc[1]*f_inc[1];
	  inc_force_norm += f_inc[2]*f_inc[2];
	}
        f_inc += 3;
      }
      global_inc_force_norm = std::sqrt(contact_global_sum(inc_force_norm, communicator));
      if (iteration == 0) {
        initial_global_inc_force_norm = global_inc_force_norm;

        if (contactPrintFlag > 1 && contact_processor_number(communicator)==0) {
          std::cout << "    alt convergence tol = "<<convergence_tolerance*initial_global_inc_force_norm<<"\n";
        }

      }

      if (contactPrintFlag > 1 && contact_processor_number(communicator)==0) {
        std::cout << "    convergence norm    = "<<global_inc_force_norm<<"\n";
      }

      if( global_inc_force_norm
         < convergence_tolerance*initial_global_inc_force_norm
       || global_inc_force_norm < convergence_tolerance ) {
#if LOCAL_PRINT_FLAG > 0 ||  CONTACT_DEBUG_PRINT_LEVEL>9
        Debug_Dump(-iteration);
#endif
        break;
      } else {
#if LOCAL_PRINT_FLAG > 0 ||  CONTACT_DEBUG_PRINT_LEVEL>9
        Debug_Dump(iteration);
#endif
      }
    }

  } // end of iteration loop

#if CONTACT_DEBUG_PRINT_LEVEL>=2
  if(iteration == num_iterations) std::cerr<<"iteration loop did not converge: rel. residual = " << global_inc_force_norm/initial_global_inc_force_norm << "\n";
#endif
  if(iteration == num_iterations && num_iterations > 1 && global_inc_force_norm > std::max(initial_global_inc_force_norm,FAIL_TOL)) {
    if(contact_processor_number(communicator)==0)
      std::cerr << "\rerror: contact enforcement failed (inc. force = " << global_inc_force_norm 
                << ", target = " << std::max(convergence_tolerance*initial_global_inc_force_norm, convergence_tolerance) << ")\n";
    if(global_inc_force_norm > initial_global_inc_force_norm) {
      // Set final updated position to the current predicted positiion
      NEW_POSITION.Duplicate_Scratch(PREDICTED_POS);
    }
  }

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Stop_Timer( iteration_time );
#endif

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
  Timer().Start_Timer( extra_ghost_len_time );
#endif

  Compute_Extra_Ghosting_Length();

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
  Timer().Stop_Timer( extra_ghost_len_time );
  Timer().Start_Timer( cleanup_time );
#endif
  //
  // KHB: 
  // I'm going to use the A_PREDICTION scratch to hold the total force on the
  // node which will then be passed through the base class functions to deal
  // with shells
  //
  Real* x_i  = PREDICTED_POS.Get_Data();
  Real* x_f  = NEW_POSITION. Get_Data();
  Real* f_t  = TOTAL_FORCE.  Get_Data();
  Real* m_sn = NODAL_MASS.   Get_Data();

  Real* disp_old_aug_0 = disp_old_aug;
  Real* disp_old_aug_1 = disp_old_aug + number_of_nodes;
  Real* disp_old_aug_2 = disp_old_aug + 2*number_of_nodes;

  for(int i=0 ; i<number_of_nodes ; ++i ){
    Real pen = m_sn[0]*dt2;
    
    Real x_diff       = x_f[0]-x_i[0];
    disp_old_aug_0[i] = x_diff;
    f_t[0]            = pen*x_diff;
    
    x_diff            = x_f[1]-x_i[1];
    disp_old_aug_1[i] = x_diff;
    f_t[1]            = pen*x_diff;
    
    x_diff            = x_f[2]-x_i[2];
    disp_old_aug_2[i] = x_diff;
    f_t[2]            = pen*x_diff;
    
    m_sn += 1;
    x_f  += 3;
    f_t  += 3;
    x_i  += 3;
  }
  
  {
    ScratchVariable* scratch_arrays[1];
    scratch_arrays[0] = &TOTAL_FORCE;
    Copy_Node_Vector_Scratch_to_Host_Arrays( 1, &Force, scratch_arrays );
  }
  if( topology->Have_Shells() ) Store_Shell_Final_Lofted_Positions();
  // Save the force in plot_force, only do this if saving CVARS
  // KHB:
  // I'm using the scratch mass, density & wavespeed to handle shells

  if(CALC_PLOT_FORCE) {
          m_sn = NODAL_MASS.     Get_Data();
    Real *d_sn = NODAL_DENSITY.  Get_Data();
    Real *w_sn = NODAL_WAVESPEED.Get_Data();
          f_t  = TOTAL_FORCE    .Get_Data();

    Real *plot_force_0 = plot_force;
    Real *plot_force_1 = plot_force + number_of_nodes;
    Real *plot_force_2 = plot_force + number_of_nodes*2;
    Real *plot_force_3 = plot_force + number_of_nodes*3;
    Real *plot_force_4 = plot_force + number_of_nodes*4;
    Real *plot_force_5 = plot_force + number_of_nodes*5;

    for(int i=0 ; i<number_of_nodes; ++i ){
      plot_force_0[i] = f_t[0];
      plot_force_1[i] = f_t[1];
      plot_force_2[i] = f_t[2];
      plot_force_3[i] = m_sn[0];
      plot_force_4[i] = d_sn[0];
      plot_force_5[i] = w_sn[0];
      f_t  += 3;
      m_sn += 1;
      d_sn += 1;
      w_sn += 1;
    }
  }
  // Compute Old timestep multiplier for augmented search
  time_mult_old_aug = dt*(0.5*(dt+dt_old));

  // Compute & set the CVARS variables if requested
  if( SAVE_CVARS ) Set_CVARS();

  // Compute & set the CVARS variables if requested
  if( SAVE_CGVARS ) Set_CGVARS();
  // Clean up this class
  TD_Enforcement_Clean_Up();
  // Call the base class clean up function
  Clean_Up();

#ifdef CONTACT_DEBUG_NODE
  for(int i=0 ; i<Number_Debug_Nodes() ; ++i ){
    ContactNode<Real>* dn = Debug_Node( i );
    if( dn ){
      int index = Get_Node_Host_Array_Index(dn);
      postream << "TD Enforcement Contact Force for Node = "
               << dn->Exodus_ID() << " = "
               << Force[index*NDIM+0]
               << " "
               << Force[index*NDIM+1]
               << " "
               << Force[index*NDIM+2]
               << "\n";
    }
  }
  postream.flush();
#endif
 
#ifdef CONTACT_HEARTBEAT 
    if (contact_processor_number(communicator)==0) {
      std::cout << "End of enforcement\n"<<flush;
    }
#endif

#if !defined (CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  Timer().Stop_Timer( cleanup_time );
  Timer().Stop_Timer( enforcement_time );
  Timer().Stop_Timer( search->Contact_time() );
#endif

  return (ContactSearch::ContactErrorCode) 
    contact_global_error_check( error_code, communicator );
}

void ContactTDEnforcement::Debug_Dump(int iteration) 
{
#if LOCAL_PRINT_FLAG > 0 || CONTACT_DEBUG_PRINT_LEVEL>9
  ContactParOStream& postream = ParOStream();

  Real *m     = NODAL_MASS.     Get_Data();
  Real *a_sn  = SN_ACCELERATION.Get_Data();
  Real *a_pre = A_PREDICTION.   Get_Data();
  Real *a_con = A_CONSTITUTIVE. Get_Data();
  Real *a_cor = A_CORRECTION.   Get_Data();
  Real *f_inc = INC_FORCE.      Get_Data();
  Real *f_tot = TOTAL_FORCE.    Get_Data();

  // Print all forces
  Real total_pre[3] = {0.0,0.0,0.0};
  Real total_cor[3] = {0.0,0.0,0.0};
  Real total_con[3] = {0.0,0.0,0.0};
  Real total_inc[3] = {0.0,0.0,0.0};
  Real total_tot[3] = {0.0,0.0,0.0};
  postream << "\n FORCES____iteration = " << iteration << "\n";
  for(int i=0 ; i<number_of_nodes ; ++i ){
	int ExoID = enforcement_node_list[i]->Exodus_ID();
        postream << ExoID << " mass  " <<m[0]<<"\n";
        postream << ExoID << " f_sn  " <<m[0]*a_sn[0] <<","
                                       <<m[0]*a_sn[1] <<","
                                       <<m[0]*a_sn[2] <<"\n";
        postream << ExoID << " f_pre " <<m[0]*a_pre[0]<<","
                                       <<m[0]*a_pre[1]<<","
                                       <<m[0]*a_pre[2]<<"\n";
        postream << ExoID << " f_cor " <<m[0]*a_cor[0]<<","
                                       <<m[0]*a_cor[1]<<","
                                       <<m[0]*a_cor[2]<<"\n";
        postream << ExoID << " f_con " <<m[0]*a_con[0]<<","
                                       <<m[0]*a_con[1]<<","
                                       <<m[0]*a_con[2]<<"\n";
        postream << ExoID << " f_inc " <<f_inc[0]<<","
                                       <<f_inc[1]<<","
                                       <<f_inc[2]<< "\n";
        postream << ExoID << " f_tot " <<f_tot[0]<<","
                                       <<f_tot[1]<<","
                                       <<f_tot[2]<< "\n";
        postream << "\n";
        for(int j=0 ; j<NDIM ; ++j ){ 
          total_pre[j] += m[0]*a_pre[j]; 
          total_cor[j] += m[0]*a_cor[j];
          total_con[j] += m[0]*a_con[j]; 
          total_inc[j] += f_inc[j]; 
          total_tot[j] += f_tot[j]; 
        }
        m += 1;
        a_sn += 3;
        a_pre += 3;
        a_con += 3;
        a_cor += 3;
        f_inc += 3;
        f_tot += 3;
  }
  postream << "_______________________________________\n";
  postream.flush();

  // Print momentum balance
  postream << "\n Total forces----iteration = " << iteration << "\n";
  postream << "-> Total predictor force " << total_pre[0] << ", " 
                                          << total_pre[1] << ", " 
                                          << total_pre[2] << "\n" ;
  postream << "-> Total accl. corrector " << total_cor[0] << ", " 
                                          << total_cor[1] << ", " 
                                          << total_cor[2] << "\n" ;
  postream << "-> Total cons. corrector " << total_con[0] << ", " 
                                          << total_con[1] << ", " 
                                          << total_con[2] << "\n" ;
  postream << "-> Total increm. force   " << total_inc[0] << ", " 
                                          << total_inc[1] << ", " 
                                          << total_inc[2] << "\n" ;
  postream << "-> Total contact force   " << total_tot[0] << ", " 
                                          << total_tot[1] << ", " 
                                          << total_tot[2] << "\n" ;
  postream << "_______________________________________\n";
  postream.flush();

  // Calculate : energy, linear and angular momentum conservation
  if (SAVE_CGVARS) {
  postream << "\nQUALITY MEASURES________________________\n";
  postream << "TOTAL Force_X   = "     << global_force_x << "\n";
  postream << "      Force_Y   = "     << global_force_y << "\n";
  postream << "      Force_Z   = "     << global_force_z << "\n";
  postream << "    ||Force||_2 = "     << global_force_norm << "\n";
  postream << "    Dissipation = "     << global_dissipation << "\n";
  postream << "||constraints||_max : " << global_constraint_norm << "\n";
  postream << "||Force_i||_2       : " << global_inc_force_norm << "\n";
  postream << "________________________________________\n";
  postream.flush();
  }

  // check dissipation at nodes
  // REJ : use CVARS instead?
  f_tot = TOTAL_FORCE.    Get_Data();
  for(Node_Constraint_Group cnei_group = node_entity_list.Node_Group_Start(); 
      cnei_group.Valid_Group();
      cnei_group = node_entity_list.Node_Group_Next()) {
    // REJ : this isn't quite right for multiple interactions
    for(int i = 0; i<cnei_group.Num_Interactions(); ++i) {
      ContactNodeEntityInteraction* cnei  = cnei_group.Get_Interaction(i);
      int interaction_index = cnei_group.Get_Index(i);
      Real* slip     = NEI_TOTAL_SLIP.Get_Scratch(interaction_index);
      ContactNode<Real>* node = cnei->Node();
      Real dissipation = 0.0;
      for (int k=0; k< NDIM; ++k) { dissipation += f_tot[k]*slip[k];}
      if (dissipation > ZERO_TOL) {
        postream << "WARNING: energy creation @ " << node->Exodus_ID()
      	         <<  "("<<i<<")"<<",  " << dissipation << "\n";
      }
    }
  }
  postream.flush();

  // Print convergence 
  f_tot = TOTAL_FORCE.  Get_Data();
  f_inc = INC_FORCE.    Get_Data();
  Real inc_norm = 0.0, tot_norm = 0.0;
  for(int i=0 ; i<number_of_nodes ; ++i ){
	Real inc_mag = 0.0, tot_mag = 0.0;
        for(int j=0 ; j<NDIM ; ++j ){ inc_mag += f_inc[j]*f_inc[j]; }
        for(int j=0 ; j<NDIM ; ++j ){ tot_mag += f_tot[j]*f_tot[j]; }
	inc_norm = std::max(inc_mag,inc_norm);
	tot_norm = std::max(tot_mag,tot_norm);
        f_inc += 3;
        f_tot += 3;
  }
  inc_norm = contact_global_maximum(inc_norm, communicator);
  tot_norm = contact_global_maximum(tot_norm, communicator);
  postream           << " iteration = " << iteration
                     << ", || f_inc ||_max = " << inc_norm
                     << ", || f_tot ||_max = " << tot_norm << "\n";
  postream.flush();
#endif
}

void ContactTDEnforcement::Verify_Enforcement_Data()
{
  //Making certain that the kinmatic partition for surfaces is entirely on the other entity.
  int surface_count = topology->Number_of_Analytic_Surfaces();
  for(int i = 0; i < surface_count; ++i)
  {
    int surface_key = topology->Analytic_Surface(i)->Entity_Key();
    for(int j = 0; j < enforcement_data.Number_Entity_Keys(); ++j)
    {
      //if this entity is the surface then we should do nothing.
      if(surface_key != j)
      {
        enforcement_data.Set_Data(KINEMATIC_PARTITION, surface_key, j, 1.0);
      }
    }
  }
}

ContactSearch::ContactErrorCode ContactTDEnforcement::TD_Enforcement_Set_Up()
{
  only_frictionless_or_tied = Only_Frictionless_or_Tied();
  //
  // Initialize the enforcement models for a new time step
  //
  for(int i=0 ; i<num_enforcement_models ; ++i ) {
    ((ContactTDEnfModel*)enforcement_models[i])->Initialize_for_Time_Step();
  }
  //
  // Modify interactions to enforce symmetry if requested
  //
  ContactSearch::ContactErrorCode error_code = ContactSearch::NO_ERROR;
  if( have_symmetry_nodes ) error_code = Symmetry_Node_Correction();
  //
  // Globally syncronize the error code
  //
  error_code = (ContactSearch::ContactErrorCode) contact_global_error_check(  error_code, communicator );
  return error_code;
}

void ContactTDEnforcement::Assemble_Masses() {
  ASSEMBLED_MASS.Zero_Scratch();
  for(Node_Constraint_Group cnei_group = node_entity_list.Node_Group_Start(); 
      cnei_group.Valid_Group();
      cnei_group = node_entity_list.Node_Group_Next()) {
    PRECONDITION(cnei_group.Num_Interactions() == 1);
    int node_index = cnei_group.Get_Constrained_Node_Index();
    Real m         = *(NODAL_MASS.     Get_Scratch(node_index));
    //	
    // Calculate the mass and acceleration contributions
    //
    const int num_face_nodes = cnei_group.Get_Num_Face_Nodes(0);
    for( int k=0 ; k < num_face_nodes ; ++k ){
      int face_node_index = cnei_group.Get_Face_Node_Index(0, k);
      Real face_shape_function = cnei_group.Get_Face_Shape_Function(0,k);
      *(ASSEMBLED_MASS.Get_Scratch(face_node_index)) += face_shape_function*m;
    }
  }
  swapadd_node_scratch(ASSEMBLED_MASS);
}

ContactSearch::ContactErrorCode ContactTDEnforcement::Symmetry_Node_Correction()
{
#if defined (CONTACT_DEBUG_NODE) || (LOCAL_PRINT_FLAG>0)
  ContactParOStream& postream = ParOStream();
#endif
    bool have_an_error = false;
    for(ContactNodeEntityInteraction *cnei = node_entity_list.Iterator_Start();
        cnei != NULL;
        cnei = node_entity_list.Iterator_Next()) {
      ContactNodeFaceInteraction *cnfi = dynamic_cast<ContactNodeFaceInteraction*>(cnei);
      if(cnfi){
#ifdef CONTACT_DEBUG_NODE
        bool PRINT_THIS_NODE = topology->Is_a_Debug_Node( cnfi->Node() );
        int node_id = -1;
#endif
        ContactFace<Real>* face = NULL;
        ContactNode<Real>* paired_node = NULL;

        int node_face_index = -1;
        bool is_a_symmetry_node = false;
        if( cnfi->Node()->Exodus_ID() == symmetry_node_1_exodus_id ){
	  is_a_symmetry_node = true;
	  // See if node_2 is part of the face
	  face = cnfi->Face();
	  for(int j=0 ; j<cnfi->Face()->Nodes_Per_Face() ; ++j ){
	    if( face->Node(j)->Exodus_ID() == symmetry_node_2_exodus_id ){
	      paired_node = face->Node(j);
	      node_face_index = j;
#ifdef CONTACT_DEBUG_NODE
	      node_id = symmetry_node_1_exodus_id;
#endif
	      break;
	    }
	  }
        }
        if( cnfi->Node()->Exodus_ID() == symmetry_node_2_exodus_id ){
	  is_a_symmetry_node = true;
	  // See if node_1 is part of the face
	  face = cnfi->Face();
	  for(int j=0 ; j<cnfi->Face()->Nodes_Per_Face() ; ++j ){
	    if( face->Node(j)->Exodus_ID() == symmetry_node_1_exodus_id ){
	      paired_node = face->Node(j);
	      node_face_index = j;
#ifdef CONTACT_DEBUG_NODE
	      node_id = symmetry_node_2_exodus_id;
#endif
	      break;
	    }
	  }
        }
        if( paired_node ){
          POSTCONDITION(face!=NULL);
          POSTCONDITION(node_face_index>=0);
	  ContactNode<Real>* slave_node = cnfi->Node();

          VariableHandle PREDICTED_POSITION = topology->Variable_Handle( ContactTopology::Predicted_Position );
	  Real* sn_pos = slave_node->Variable(PREDICTED_POSITION);
	  Real* ms_pos = paired_node->Variable(PREDICTED_POSITION);

	  Real vec[3];
	  vec[0] = ms_pos[0] - sn_pos[0];
	  vec[1] = ms_pos[1] - sn_pos[1];
	  vec[2] = ms_pos[2] - sn_pos[2];
	  Real mag = std::sqrt( vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2] );
	  Real* n_dir = cnfi->Vector_Var(ContactNodeFaceInteraction::NORMAL_DIR);
	  Real* pb_dir = cnfi->Vector_Var(ContactNodeFaceInteraction::PUSHBACK_DIR);
	  Real* LC = cnfi->Vector_Var(ContactNodeFaceInteraction::COORDINATES);
	  Real vec_dot_n = vec[0]*n_dir[0] + vec[1]*n_dir[1] + vec[2]*n_dir[2];
	  if( mag ){
	    vec[0] /= mag;
	    vec[1] /= mag;
	    vec[2] /= mag;
	    if( vec_dot_n > 0.0 ){
	      mag *= -1.0;
	    } else {
	      vec[0] *= -1.0;
	      vec[1] *= -1.0;
	      vec[2] *= -1.0;
	    }
	  } else {
	    vec[0] = pb_dir[0];
	    vec[1] = pb_dir[1];
	    vec[2] = pb_dir[2];
	  }
	  // following only works for 4 node quads. Error if not.
          if( face->FaceType() != ContactSearch::QUADFACEL4 ) {
	    std::sprintf( message, "Face paired with symmetry node %d is not QUADFACE4L. Results are suspect.", 
		  slave_node->Exodus_ID() );
	    errors->Add_Error_Message(message);
	    have_an_error = true;
	  }
	  Real lc[2] = {0.0, 0.0};
	  switch( node_face_index ){
	  case 0:
	    lc[0] = -1.0;
	    lc[1] = -1.0;
	    break;
	  case 1:
	    lc[0] =  1.0;
	    lc[1] = -1.0;
	    break;
	  case 2:
	    lc[0] =  1.0;
	    lc[1] =  1.0;
	    break;
	  case 3:
	    lc[0] = -1.0;
	    lc[1] =  1.0;
	    break;
	  }
#ifdef CONTACT_DEBUG_NODE
	  if( PRINT_THIS_NODE ){
	    postream << "Resetting the Interaction for Node " 
		     << node_id << " due to symmtery" << "\n";
	    postream << "   Modifying local coords from (" << LC[0] << "," 
		     << LC[1] << ") to (" << lc[0] << "," << lc[1] << ")\n";
	    postream << "   Modifying gap_cur from " 
		     << cnfi->Scalar_Var(ContactNodeFaceInteraction::GAP_CUR) 
		     << " to " 
		     << mag << "\n";
	    postream << "   Modifying gap_old from "
		     << cnfi->Scalar_Var(ContactNodeFaceInteraction::GAP_OLD)
		     << " to 0.0\n";
	    postream << "   Modifying pb_dir from (" << pb_dir[0] << ","
		     << pb_dir[1] << "," << pb_dir[2] << ") to (" 
		     << vec[0] << "," << vec[1] << "," << vec[2] << ")\n";
	    postream << "   Modifying n_dir from  (" << n_dir[0] << ","
		     << n_dir[1] << "," << n_dir[2] << ") to ("
		     << vec[0] << "," << vec[1] << "," << vec[2] << ")\n";
          }
#endif
	  pb_dir[0] = vec[0];
	  pb_dir[1] = vec[1];
	  pb_dir[2] = vec[2];
	  n_dir[0] = vec[0];
	  n_dir[1] = vec[1];
	  n_dir[2] = vec[2];
	  cnfi->Scalar_Var(ContactNodeFaceInteraction::GAP_CUR) = mag;
	  cnfi->Scalar_Var(ContactNodeFaceInteraction::GAP_OLD) = 0.0;
	  LC[0] = lc[0];
	  LC[1] = lc[1];
        } else if( is_a_symmetry_node ){
	  std::sprintf( message, "Symmetry Node %d is paired with a face that does not contain its paired node.", cnfi->Node()->Exodus_ID() );
	  errors->Add_Error_Message(message);
	  have_an_error = true;
        }
      }
    }
    if( contact_global_error_check( have_an_error, communicator ) )  {
      return ContactSearch::INTERNAL_ERROR;
    } else {
      return ContactSearch::NO_ERROR;
    }
}

void ContactTDEnforcement::Compute_Penalty_Forces()
{
#if defined (CONTACT_DEBUG_NODE) || (LOCAL_PRINT_FLAG)
  ContactParOStream& postream = ParOStream();
#endif
  Real f_sn[3]; 
  // Now pushback directions are fixed . . . . . . . .
  // Form unique gap vector and force prediction per node
  // f_predictor @ sn = mass@node * gap@node / (1/2(dt+dt_old)dt) 
  // 
  Real* normals[3];
  Real gaps[3], g[3];

  for(Node_Constraint_Group cnei_group = node_entity_list.Node_Group_Start();
      cnei_group.Valid_Group();
      cnei_group = node_entity_list.Node_Group_Next()) {
    int ncc = cnei_group.Num_Interactions();
    int node_index = cnei_group.Get_Constrained_Node_Index();
    if(ncc == 1) {
      ContactNodeEntityInteraction  *cnei = cnei_group.Get_Interaction(0);
      int interaction_index = cnei_group.Get_Index(0);
      Real kin_partition = *(NEI_KINEMATIC_PARTITION.Get_Scratch(interaction_index));  
      if (kin_partition > 0.0) {
        Real* a_sn =   SN_ACCELERATION.Get_Scratch(node_index);
        Real  m_sn = *(NODAL_MASS.     Get_Scratch(node_index));
        Real* n_dir = cnei->Get_Normal();
        Real  gap   = *(NEI_GAP_TO_ENFORCE.Get_Scratch(interaction_index));
	Real pen = kin_partition*m_sn*dt2;
        Real f_n = (gap < 0.0) ? pen*gap : 0.0;
        for(int j=0 ; j<NDIM ; ++j ){ 
          f_sn[j] = -f_n*n_dir[j]; 
        } 
	const int model_type = TDEnfModel(cnei)->Interaction_Type(cnei);  
        if (model_type ==ContactTDEnfModel::TDEM_SPRINGY){
          Real* slip     = NEI_TOTAL_SLIP.Get_Scratch(interaction_index);
          Real* rel_disp = NEI_REL_DISP  .Get_Scratch(interaction_index);
          Real s_n = 0.0;
          for(int j=0 ; j<NDIM ; ++j ){ s_n += slip[j]*n_dir[j];}          
          Real d_n = 0.0;
          for(int j=0 ; j<NDIM ; ++j ){ d_n += rel_disp[j]*n_dir[j];}          
          f_n = (s_n < 0.0 && d_n < 0.0) ? pen*s_n : 0.0;
          for(int j=0 ; j<NDIM ; ++j ){ f_sn[j] = -f_n*n_dir[j];}   
        }
	else if (model_type == ContactTDEnfModel::TDEM_TIED ) {
// REJ this "if" will probably be replaced by the TIED/GLUED distinction
          ContactEnfModel *model      = TDEnfModel(cnei);
          ContactTDTied   *tied_model = dynamic_cast<ContactTDTied*>(model);
          if (tied_model != NULL && new_tied_enforcement == ContactSearch::ACTIVE) {
            Real* rel_disp = NEI_REL_DISP  .Get_Scratch(interaction_index);
	    for(int j=0 ; j<NDIM ; ++j ){ f_sn[j] = -pen*rel_disp[j];}
          } else {
            Real* slip = NEI_TOTAL_SLIP.Get_Scratch(interaction_index);
            for(int j=0 ; j<NDIM ; ++j ){ f_sn[j] = -pen*slip[j];}
          }
	}
        Real m_sn_inv = 1.0/m_sn;
	for(int j=0 ; j<NDIM ; ++j ){ a_sn[j] = f_sn[j] * m_sn_inv; }
      }
    } else if ( ncc > 1) {
      Get_Normals_and_Gaps_at_Node(cnei_group, normals, gaps);
      Unify( normals, gaps, g );
      Real* a_sn = SN_ACCELERATION.Get_Scratch(node_index);
      for (int j=0 ; j<NDIM ; ++j ){
	//
  	// The minus sign here is because gap is negative for penetration
	//
        a_sn[j] = -g[j]*dt2;
      }
    }
  }
#ifdef CONTACT_DEBUG_NODE
  postream.flush();
  for( int i=0 ; i<Number_Debug_Nodes() ; ++i ){
    ContactNode<Real>* dn = Debug_Node( i );
    if( dn ){
      int node_index = dn->EnfArrayIndex();
      Real m_sn  = *(NODAL_MASS.     Get_Scratch(node_index));
      Real* a_sn =   SN_ACCELERATION.Get_Scratch(node_index);
      postream << "Penalty Force for node " << dn->Exodus_ID() << " = " 
               << m_sn*a_sn[0] << " " << m_sn*a_sn[1] << " " << m_sn*a_sn[2] 
               << "\n";
    }
  }
  postream.flush();
#endif
}

void ContactTDEnforcement::Apply_Face_Forces(Real face_force[3], 
                                             Node_Constraint_Group &cnei_group,
                                             int constraint_num,
                                             Real m) {


  ContactNodeEntityInteraction * cnei = cnei_group.Get_Interaction(constraint_num);
  const int num_face_nodes = cnei_group.Get_Num_Face_Nodes(constraint_num);
  ContactNodeFaceInteraction * cnfi = dynamic_cast<ContactNodeFaceInteraction*>(cnei);
  ContactEnfModel *model      = TDEnfModel(cnei);
  ContactTDTied   *tied_model = dynamic_cast<ContactTDTied*>(model);
  //
  //  Use special processing for tied node face contact when lofting of the node is present.  For all other cases, simply weight the node
  //  forces onto the face nodes by the face shape functions.
  //
  if(cnfi != NULL && tied_model != NULL && new_tied_enforcement == ContactSearch::ACTIVE) {
    //
    //  Extract face vital information.  Need to get enough data to predict how the lofted tied constraint will move as the face
    //  translates and rotates.
    //
    Real* x_face[MAX_NODES_PER_FACE];
    ContactFace<Real>* face = cnfi->Face();
    Real node_positions[MAX_NODES_PER_FACE][3];
    Real n_dir[3] = {0,0,0}; 
    Real gap0 = cnei->Get_Gap_Initial();
    Real xp_face[3] = {0.0, 0.0, 0.0};
    for(int ii=0 ; ii<num_face_nodes ; ++ii ){
      const int node_face_index = cnei_group.Get_Face_Node_Index(constraint_num, ii);
      const Real face_shape_function = cnei_group.Get_Face_Shape_Function(constraint_num ,ii);
      Real* x_mn_new = NEW_POSITION.Get_Scratch(node_face_index);
      for(int j=0 ; j<NDIM ; ++j ){
	xp_face[j] +=  face_shape_function*x_mn_new[j] ;
      }
      x_face[ii] = NEW_POSITION.Get_Scratch(node_face_index);
      node_positions[ii][0] = x_face[ii][0];
      node_positions[ii][1] = x_face[ii][1];
      node_positions[ii][2] = x_face[ii][2];
    }
    face->Compute_Normal(x_face,cnfi->Get_Coordinates(),n_dir);

    //
    //  Compute the moment applied to the face by the force 'face_force' acting at the distance of gap0.  
    //
    Real in_plane_force[3];  
    Real f_dot_n = dot_vectors(face_force, n_dir);
    in_plane_force[0] = face_force[0] - f_dot_n * n_dir[0];
    in_plane_force[1] = face_force[1] - f_dot_n * n_dir[1];
    in_plane_force[2] = face_force[2] - f_dot_n * n_dir[2];

    Real force_mag = norm_vector(face_force);
    Real in_plane_f_mag = norm_vector(in_plane_force);
    Real moment_mag = in_plane_f_mag * gap0;
    //
    //  Pick two points on the face that lie in the same line as the force.  Equal and opposite
    //  force will be applied to these two points to add in a resultant face moment.
    //
    Real p1[3], p2[3];
    Real in_plane_dir[3] = {node_positions[0][0] - node_positions[2][0], 
                            node_positions[0][1] - node_positions[2][1], 
			    node_positions[0][2] - node_positions[2][2]};
    if(in_plane_f_mag > 0.0) {
      in_plane_dir[0] = in_plane_force[0] / in_plane_f_mag;
      in_plane_dir[1] = in_plane_force[1] / in_plane_f_mag;
      in_plane_dir[2] = in_plane_force[2] / in_plane_f_mag;
    }
    Real max_coord_in_plane = node_positions[0][0] * in_plane_dir[0] + 
      node_positions[0][1] * in_plane_dir[1] +
      node_positions[0][2] * in_plane_dir[2];
    Real min_coord_in_plane = max_coord_in_plane;
    for(int i = 0; i < num_face_nodes; ++i) {
      Real cur_coord = node_positions[i][0] * in_plane_dir[0] + 
	node_positions[i][1] * in_plane_dir[1] +
	node_positions[i][2] * in_plane_dir[2];

      if(cur_coord > max_coord_in_plane) max_coord_in_plane = cur_coord;
      if(cur_coord < min_coord_in_plane) min_coord_in_plane = cur_coord;
    }
    Real face_size = max_coord_in_plane - min_coord_in_plane;

    Real torque_arm = face_size * 0.1;
    p1[0] = xp_face[0] + torque_arm*in_plane_dir[0];
    p1[1] = xp_face[1] + torque_arm*in_plane_dir[1];
    p1[2] = xp_face[2] + torque_arm*in_plane_dir[2];
    p2[0] = xp_face[0] - torque_arm*in_plane_dir[0];
    p2[1] = xp_face[1] - torque_arm*in_plane_dir[1];
    p2[2] = xp_face[2] - torque_arm*in_plane_dir[2];

    Real M_force = 0.0;
    if(torque_arm > 0.0) {
      M_force = moment_mag/(2.0*torque_arm);
    }
    Real local_coords1[3], local_coords2[3];
    face->Compute_Local_Coords(node_positions, p1, local_coords1);
    face->Compute_Local_Coords(node_positions, p2, local_coords2);
    Real shape1[MAX_NODES_PER_FACE], shape2[MAX_NODES_PER_FACE];
    face->Evaluate_Shape_Functions(local_coords1, shape1);
    face->Evaluate_Shape_Functions(local_coords2, shape2);
    //
    //  Weight the moment contribution vs, standard force.  Both forces will move nodes, need to avoid moving the nodes
    //  too far as this can yield an instability.  Some portion of the force will be used as moment, and some as direct
    //  force.  The actuall weighting factor is a huerstic guess.
    //
    Real moment_weight = 0.0;
    if(gap0 * in_plane_f_mag + face_size*force_mag > 0.0) {
      moment_weight = gap0*in_plane_f_mag / ( gap0 * in_plane_f_mag + face_size * force_mag);
    }
    //
    //  saftey factor, when using lofted constraints a displacement of node may be magnified.  Apply
    //  the saftey factor to the final forces to avoid large overshooting of the solution.  A force to the underlying
    //  nodes to statisfy that constraint may overshoot the solution by an excessive ammount.  The saftey factor
    //  scales down resultant stresses to avoid the overshoot.  Note, iterative enforcement of tied contstraints with
    //  gaps is a nessecity even for simple pure master/slave cases.  As the face is 2D, a motion in one node can cause
    //  the lofted node to move in two directions, thus need to do a vector sum of the moment multiplication factors.
    //
    static Real sqrt_2 = sqrt(2.0);
    Real saftey_factor = (face_size / (gap0 * sqrt_2 + face_size));
    for( int k=0 ; k < num_face_nodes ; ++k ){
      Real contact_point_shape_function = cnei_group.Get_Face_Shape_Function(constraint_num,k);
      Real moment_point_1_shape_function = shape1[k]; 
      Real moment_point_2_shape_function = shape2[k];
      int face_node_index = cnei_group.Get_Face_Node_Index(constraint_num, k);
      Real* f_K = ASSEMBLED_FORCE.Get_Scratch(face_node_index);
      *(ASSEMBLED_MASS.Get_Scratch(face_node_index)) += contact_point_shape_function*m;
      for(int j=0 ; j<NDIM ; ++j ){ 
        Real direct_force = contact_point_shape_function * face_force[j]; 
        Real moment_force = (moment_point_2_shape_function -  moment_point_1_shape_function) * M_force * n_dir[j];
        f_K[j] += (direct_force * (1.0-moment_weight) + moment_force * moment_weight)*saftey_factor;
      }
    }
  } else {
    for( int k=0 ; k < num_face_nodes ; ++k ){
      const int face_node_index = cnei_group.Get_Face_Node_Index(constraint_num, k);
      const Real face_shape_function = cnei_group.Get_Face_Shape_Function(constraint_num, k);
      Real* f_K = ASSEMBLED_FORCE.Get_Scratch(face_node_index);
      *(ASSEMBLED_MASS.Get_Scratch(face_node_index)) += face_shape_function*m;
      //
      // force on the face is equal and opposite to the slave node force
      //
      for(int j=0 ; j<NDIM ; ++j ) {
        f_K[j] += face_shape_function * face_force[j];
      }
    }
  }
}



void ContactTDEnforcement::Assemble_Nodal_Forces_and_Masses()
{ 
#ifdef CONTACT_DEBUG_NODE
  ContactParOStream& postream = ParOStream();
#endif

  for(Node_Constraint_Group cnei_group = node_entity_list.Node_Group_Start(); 
      cnei_group.Valid_Group();
      cnei_group = node_entity_list.Node_Group_Next()) {
    int ncc = cnei_group.Num_Interactions();
    if( ncc == 1) {
      int node_index = cnei_group.Get_Constrained_Node_Index();
      Real m       = *(NODAL_MASS.     Get_Scratch(node_index));
      Real* a_sn   =   SN_ACCELERATION.Get_Scratch(node_index);
      Real face_force[3];
      face_force[0] = -m * a_sn[0];      
      face_force[1] = -m * a_sn[1];      
      face_force[2] = -m * a_sn[2];      
      Apply_Face_Forces(face_force, cnei_group, 0, m);
    } else if ( ncc > 1) {

      Real f_sn[3];
      Real* normals[3];
      Real gaps[3], factors[3], f_constr[3][3];
      //
      // (1) Partition 
      //
      Get_Normals_and_Gaps_at_Node(cnei_group, normals, gaps);
      int node_index = cnei_group.Get_Constrained_Node_Index();
      Real m     = *(NODAL_MASS     .Get_Scratch(node_index));
      Real* a_sn =   SN_ACCELERATION.Get_Scratch(node_index);
#ifdef CONTACT_DEBUG_NODE
      for( int i=0 ; i<Number_Debug_Nodes() ; i++ ){
	ContactNode<Real>* dn = Debug_Node( i );
	if (!dn) continue;
	ContactNodeEntityInteraction * cnei = 
	  cnei_group.Get_Interaction(ncc-1);
	if (!cnei) continue;
	ContactNode<Real>* node = cnei->Node();
	if( dn->Exodus_ID() == node->Exodus_ID() ){
	  postream << "Node " << dn->Exodus_ID() << " in assemble loop, mult constr.\n";
	}
      }
#endif
      for(int j=0 ; j<NDIM ; ++j ){ f_sn[j] = m * a_sn[j]; }
      Real surface_stiffness_inv[3] = {0.0,0.0,0.0};
      //
      // get surface stiffness of NF constraints
      //
      for(int j=0 ; j< ncc ; ++j ){
	//
	// Strictly speaking, we need the delta_t's.  However, since all we
        // are going to use this for is ratio's we will simply neglect it.
	//
        const int num_face_nodes = cnei_group.Get_Num_Face_Nodes(j);
        for(int k=0 ; k < num_face_nodes ; ++k ){
          int face_node_index = cnei_group.Get_Face_Node_Index(j, k);
          Real face_shape_function = cnei_group.Get_Face_Shape_Function(j, k);
          surface_stiffness_inv[j] += face_shape_function*(*(NODAL_MASS.Get_Scratch(face_node_index)));
	}
        //
	//  Analytic surfaces are assumed to be inifinitly stiff.  However, they have no nodes so the above calculation
	//  will return zero stiffness.  Correct that here.
	//
        if(num_face_nodes == 0) {
  	  surface_stiffness_inv[j] = 0.0;
	} else {
  	  surface_stiffness_inv[j] = 1.0/surface_stiffness_inv[j];
	}
      }
#ifdef PJSA_FIX
      if(ncc == 2) { // project f_sn onto plane
        double toto = f_sn[0]*normals[2][0]+f_sn[1]*normals[2][1]+f_sn[2]*normals[2][2];
        for(int j=0; j<3; ++j) f_sn[j] -= normals[2][j]*toto;
      }
#endif
      Partition_Force(ncc, normals, f_sn, surface_stiffness_inv, f_constr );
      Real m_sn  = *(NODAL_MASS.Get_Scratch(node_index));
      Partition_Mass(ncc, normals, f_sn, factors);
      //
      // (2) Assemble to faces
      // NOTE : this step will lose all force and mass to be assembled 
      // to the analytic surfaces.
      //
      for( int j=0 ; j<ncc ; ++j ) {
	m = m_sn*factors[j];
        Real face_forces[3];
        face_forces[0] = -f_constr[j][0];
        face_forces[1] = -f_constr[j][1];
        face_forces[2] = -f_constr[j][2];
	Apply_Face_Forces(face_forces, cnei_group, j, m);
      }
    }
  }
#ifdef CONTACT_DEBUG_NODE
  postream.flush();
  for( int i=0 ; i<Number_Debug_Nodes() ; i++ ){
    ContactNode<Real>* dn = Debug_Node( i );
    if( dn ){
      Real* af = ASSEMBLED_FORCE.Get_Scratch(dn);
      postream << "Interface at Node " << dn->Exodus_ID() << "\n";
      postream << "  Assembled Force = " << af[0] << " " << af[1] << " " 
	       << af[2] << "\n";
    }
  }
  postream.flush();

#endif
}


void ContactTDEnforcement::Assemble_Nodal_Forces_No_Multiple()
{ 
  for(Node_Constraint_Group cnei_group = node_entity_list.Node_Group_Start(); 
      cnei_group.Valid_Group();
      cnei_group = node_entity_list.Node_Group_Next()) {
    PRECONDITION(cnei_group.Num_Interactions() == 1);
    int node_index = cnei_group.Get_Constrained_Node_Index();
    Real m       = *(NODAL_MASS.     Get_Scratch(node_index));
    Real* a_sn   =   SN_ACCELERATION.Get_Scratch(node_index);
    Real face_force[3];
    face_force[0] = -m * a_sn[0];      
    face_force[1] = -m * a_sn[1];      
    face_force[2] = -m * a_sn[2];      
    Apply_Face_Forces(face_force, cnei_group, 0, 0.0);
  }
}

void ContactTDEnforcement::Response_Prediction(Real mult)
{
  mult = 1.0;
  // Calculate the acceleration response of the surface (by nodes)
  // This is the PREDICTOR
#ifdef CONTACT_DEBUG_NODE
  ContactParOStream& postream = ParOStream();
#endif

  Real *a_pre   = A_PREDICTION.   Get_Data();
  Real *f       = ASSEMBLED_FORCE.Get_Data();
  Real *m_sn    = NODAL_MASS.     Get_Data();
  Real *m_assem = ASSEMBLED_MASS. Get_Data();

  for(int i=0 ; i<number_of_total_nodes ; ++i ){
    Real mass_sum_inv = 1.0 / (m_sn[0] + m_assem[0] * mult);
    a_pre[0] = f[0] * mass_sum_inv;
    a_pre[1] = f[1] * mass_sum_inv;
    a_pre[2] = f[2] * mass_sum_inv;
    a_pre   += 3;
    f       += 3;
    m_sn    += 1;
    m_assem += 1;
  }
#ifdef CONTACT_DEBUG_NODE
  for(int i=0 ; i<Number_Debug_Nodes() ; ++i ){
    ContactNode<Real>* dn = Debug_Node( i );
    if( dn ){
      const int node_index = dn->EnfArrayIndex();
      Real m = *(NODAL_MASS.  Get_Scratch(node_index));
      Real* ap = A_PREDICTION.Get_Scratch(node_index);
      postream << "Predictor Force for Node "
	       << dn->Exodus_ID() << " = " 
	       << m*ap[0] << " " << m*ap[1] << " " << m*ap[2] << "\n";
    }
  }
  postream.flush();
#endif
}

void ContactTDEnforcement::Response_Correction()
{
  // Compute a normal correction for nodes to account for how the master surface
  // responded as a result of "interface motion".  
  //
  // Note: There is no correction for analytic surfaces since they aren't 
  // currently allowed to move.  All we do for them is assemble the force to 
  // the constraint for friction.

#ifdef CONTACT_DEBUG_NODE
  ContactParOStream& postream = ParOStream();
#endif
  // used in force estimate for mult. interactions 
  Real a_cor_mag;
  // Kinematically admissible correction for surface motion (by interaction)
  // This is the CORRECTOR, based on acceleration mismatch
  Real* normals[3];
  Real gaps[3];
  Real a_ms[3], mag_a_sn[3];
  Real a_pre_mag[3];

  for(Node_Constraint_Group cnei_group = node_entity_list.Node_Group_Start();
      cnei_group.Valid_Group();
      cnei_group = node_entity_list.Node_Group_Next()) {
    int ncc = cnei_group.Num_Interactions();
    if(ncc == 1) {
      ContactNodeEntityInteraction  *cnei = cnei_group.Get_Interaction(0);
      const int model_type = TDEnfModel(cnei)->Interaction_Type(cnei);
      Real a_mf[3] = {0.0, 0.0, 0.0};
      //
      // Compute the acceleration of the master surface at the contact pt.
      //
      const int num_face_nodes = cnei_group.Get_Num_Face_Nodes(0);
      for(int k=0 ; k < num_face_nodes ; ++k ){
        const int face_node_index = cnei_group.Get_Face_Node_Index(0, k);
        const Real face_shape_function = cnei_group.Get_Face_Shape_Function(0, k);
        Real* a_mf_K = A_PREDICTION.Get_Scratch(face_node_index);
        for(int j=0 ; j<NDIM ; ++j ){ a_mf[j] += a_mf_K[j]*face_shape_function; }
      }
      //
      // acceleration due to the penalty force on the slave node
      //
      const int node_index = cnei_group.Get_Constrained_Node_Index();
      Real* a_sn  =   SN_ACCELERATION.Get_Scratch(node_index);
      Real* a_cor =   A_CORRECTION   .Get_Scratch(node_index);
      // corrector relase 
      // alt: use  a_cor[j] = a_mf[j] + a_sn[j]; for frictionless && a_n >0
      // alt: no release
      if( model_type == ContactTDEnfModel::TDEM_FRICTIONLESS || model_type == ContactTDEnfModel::TDEM_FRICTIONAL ){
	Real a_n = 0.0;
	Real* nm = cnei->Get_Normal();
	for(int j=0 ; j<NDIM ; ++j){ a_n += (a_mf[j] + a_sn[j])*nm[j]; }
	if (a_n > 0.0) { 
	  for(int j=0 ; j<NDIM ; ++j){ a_cor[j] = a_n*nm[j]; }
	}
      } else {
	for(int j=0 ; j<NDIM ; ++j){ a_cor[j] = a_mf[j] + a_sn[j]; }
      }
      
      int interaction_index = cnei_group.Get_Index(0);
      Real* f_tot = TOTAL_FORCE.Get_Scratch(node_index);
      Real* f_est = NEI_FORCE.Get_Scratch(interaction_index);
      f_est[0] = f_tot[0];
      f_est[1] = f_tot[1];
      f_est[2] = f_tot[2];
    } else if( ncc>1 ) {
      //
      // (1) Partition 
      //
      Get_Normals_and_Gaps_at_Node(cnei_group, normals, gaps);
      int node_index = cnei_group.Get_Constrained_Node_Index();
      //compute master surface accelerations
      Real* a_cor  = A_CORRECTION.Get_Scratch(node_index);
      Real mag_a_ms[3] = {0.0, 0.0, 0.0};
      for(int j=0 ; j<ncc ; ++j ){
	Real* nm = cnei_group.Get_Interaction(j)->Get_Normal();
	//
        // acceleration at contact point
        //
        // the normal acceleration of the master surface at the contact pt.
        Real nm_a_mf = 0.0 ;
        const int num_face_nodes = cnei_group.Get_Num_Face_Nodes(j);
        for(int k=0 ; k< num_face_nodes ; ++k ){
          const int face_node_index = cnei_group.Get_Face_Node_Index(j, k);
          const Real face_shape_function = cnei_group.Get_Face_Shape_Function(j,k);
          Real* a_mf = A_PREDICTION.Get_Scratch(face_node_index);
	  for(int l=0 ; l<NDIM ; ++l ) nm_a_mf += nm[l]*a_mf[l] * face_shape_function;
	}
	mag_a_ms[j] = nm_a_mf;
      }
#ifdef CONTACT_DEBUG_NODE
      for( int i=0 ; i<Number_Debug_Nodes() ; i++ ){
	ContactNode<Real>* dn = Debug_Node( i );
	if (!dn) continue;
	ContactNodeEntityInteraction  *cnei = cnei_group.Get_Interaction(0);
	if (!cnei) continue;
	ContactNode<Real>* node = cnei->Node();
	if( dn->Exodus_ID() == node->Exodus_ID() ){
	  postream << "Master surface acceleration for node : "
		   << node->Exodus_ID()
		   << " " << mag_a_ms[0] << ","
		   << mag_a_ms[1] << ","
		   << mag_a_ms[2] << "\n";
	}
      }
#endif
      // (1) Unify then Partition accelerations for master surfaces
      Unify (normals, mag_a_ms, a_ms);
      Partition_Gap( normals, mag_a_ms, a_ms );
      // (2) Parition acceleration for the slave node
      Real* a_sn =   SN_ACCELERATION.Get_Scratch(node_index);
      Real m_sn  = *(NODAL_MASS.     Get_Scratch(node_index));
      Partition_Gap(normals, mag_a_sn, a_sn );
      Real* a_pre =  A_PREDICTION.   Get_Scratch(node_index);
      Partition_Gap(normals,a_pre_mag, a_pre);
      int k = 0;
      for (int ii=0 ; ii<ncc ; ++ii ){
        int nei_index = cnei_group.Get_Index(ii);
	Real* f_est   = NEI_FORCE.Get_Scratch(nei_index);
        Real* nm      = normals[k];
        a_cor_mag     = mag_a_ms[k] + mag_a_sn[k];
        for(int j=0 ; j<NDIM ; ++j ){
	  f_est[j] += m_sn*(a_pre_mag[k] + a_cor_mag)*nm[j]; 
	  a_cor[j] +=       a_cor_mag*nm[j];
	}
        ++k; 
      }
    }
  }

    
#ifdef CONTACT_DEBUG_NODE
  for(int i=0 ; i<Number_Debug_Nodes() ; ++i ){
    ContactNode<Real>* dn = Debug_Node( i );
    if( dn ){
      const int node_index = dn->EnfArrayIndex();
      Real m = *(NODAL_MASS.  Get_Scratch(node_index));
      Real* ac = A_CORRECTION.Get_Scratch(node_index);
      postream << "Corrector Force for Node "
	       << dn->Exodus_ID() << " = " 
	       << m*ac[0] << " " << m*ac[1] << " " << m*ac[2] << "\n";
    }
  }
  postream.flush();
#endif
  
}

void ContactTDEnforcement::Response_Constitutive_Correction()
{
  // This is a 2nd CORRECTOR, based on constitutive behavior
#if defined (CONTACT_DEBUG_NODE) || LOCAL_PRINT_FLAG>0
  ContactParOStream& postream = ParOStream();
#endif

  Real f_trial[3]; 
  Real* normals[3];
  Real gaps[3];

  for(Node_Constraint_Group cnei_group = node_entity_list.Node_Group_Start();
      cnei_group.Valid_Group();
      cnei_group = node_entity_list.Node_Group_Next()) {
    int ncc = cnei_group.Num_Interactions();
    if(ncc == 1) {
      ContactNodeEntityInteraction  *cnei = cnei_group.Get_Interaction(0);
      int interaction_index = cnei_group.Get_Index(0);
      int node_index = cnei_group.Get_Constrained_Node_Index();
      //
      //  Update the current force by the force increment
      //
      Real* f_inc = INC_FORCE.  Get_Scratch(node_index);
      Real* f_est = NEI_FORCE.Get_Scratch(interaction_index);
      f_est[0] += f_inc[0];
      f_est[1] += f_inc[1];
      f_est[2] += f_inc[2];
      
      Real kin_partition = *(NEI_KINEMATIC_PARTITION.Get_Scratch(interaction_index));

      ContactEnfModel       *model              = TDEnfModel(cnei);
      ContactTDTied         *tied_model         = dynamic_cast<ContactTDTied*>(model);
      ContactTDFrictionless *frictionless_model = dynamic_cast<ContactTDFrictionless*>(model);
      // fine-grained test consistent with the global Is_Frictionless_or_Tied
      if ((kin_partition>0.0) && (frictionless_model==NULL) && (tied_model==NULL)){
	Real  m_sn     = *(NODAL_MASS.Get_Scratch(node_index));
	Real* slip     = NEI_TOTAL_SLIP.Get_Scratch(interaction_index);
	Real* rel_disp = NEI_REL_DISP  .Get_Scratch(interaction_index);
	Real* nm       = cnei->Get_Normal();
	Real  gap      = *(NEI_GAP_TO_ENFORCE.Get_Scratch(interaction_index));
	Real area_sn   = cnei->Get_Node_Area();
	if (area_sn <= 0.0) area_sn = 1.0; // protect against zero/neg. areas
	const int model_type = TDEnfModel(cnei)->Interaction_Type(cnei);
	for(int k=0 ; k<NDIM ; ++k) { f_trial[k] = f_est[k] ;}
	Real tot_eff_mass = 0.0;
	const int num_face_nodes = cnei_group.Get_Num_Face_Nodes(0);
	Real face_nodal_partition[MAX_NODES_PER_FACE];
	Real *a_increm_master[MAX_NODES_PER_FACE];
	//
	// master surface
	//
	for(int inode = 0; inode < num_face_nodes; ++inode) {
          const Real face_shape_function = cnei_group.Get_Face_Shape_Function(0, inode);
          const int face_node_index = cnei_group.Get_Face_Node_Index(0, inode);
	  face_nodal_partition[inode] = face_shape_function*m_sn / *(NODAL_MASS.Get_Scratch(face_node_index));
	  a_increm_master     [inode] = A_CONSTITUTIVE.Get_Scratch(face_node_index);
	  tot_eff_mass += face_nodal_partition[inode];
	}
	if (model_type == ContactTDEnfModel::TDEM_FRICTIONAL ||  
	    model_type == ContactTDEnfModel::TDEM_ADHESIVE ||  
	    model_type == ContactTDEnfModel::TDEM_ADHESIVE_FRICTION ||  
	    model_type == ContactTDEnfModel::TDEM_SPRINGY ) {
	  //
	  // to emulate the predictor phase
	  //
	  Real pen =  m_sn*dt2/(1.0+tot_eff_mass);
	  if (model_type == ContactTDEnfModel::TDEM_ADHESIVE_FRICTION || model_type == ContactTDEnfModel::TDEM_FRICTIONAL ) {
	    Real* slip_t = NEI_TAN_SLIP.Get_Scratch(interaction_index);
	    for(int j=0 ; j<NDIM ; ++j ){ f_trial[j] += -pen*slip_t[j]; }
	  }
	  if ((model_type == ContactTDEnfModel::TDEM_ADHESIVE || 
               model_type == ContactTDEnfModel::TDEM_ADHESIVE_FRICTION) && 
	      TDEnfModel(cnei)->Active_Interaction(cnei,gap) 
	      && gap > 0.0) { // for adhesion
	    for(int j=0 ; j<NDIM ; ++j ){ f_trial[j] += -pen*gap*nm[j]; }
	  }
	  if (model_type == ContactTDEnfModel::TDEM_SPRINGY ) {
	    Real* slip_t = NEI_TAN_SLIP.Get_Scratch(interaction_index);
	    for(int j=0 ; j<NDIM ; ++j ){ f_trial[j] += -pen*slip_t[j]; }
	    if (gap > 0.0) 
	      for(int j=0 ; j<NDIM ; ++j ){ f_trial[j] += -pen*gap*nm[j]; }
	  }
	}

	TDEnfModel(cnei)->Limit_Force(cnei, gap, rel_disp, slip, nm, dt, area_sn, f_trial);
	Real* a_con = A_CONSTITUTIVE.Get_Scratch(node_index);
        
	//
	// add (increm. force) correction to both surfaces maintain mom. 
        // balance slave node
	//
        Real m_sn_inv = 1.0/m_sn;
	for(int k=0 ; k<NDIM ; ++k ){
	  //
	  // difference with actual forces to form corrector
	  //
	  const Real a_inc = kin_partition*(f_trial[k] - f_est[k]) * m_sn_inv;
	  a_con[k] += a_inc;
	  f_est[k] = f_trial[k]; // reset NEI_FORCE, ~ adding m*a_inc

	  for(int j=0 ;  j < num_face_nodes ; ++j ){
	    // the portion of the total force for the master node j
	    // i.e. f_j = N_j f_sn
	    a_increm_master[j][k] -= face_nodal_partition[j]*a_inc;
	  }
	}
	
      }
    } else if ( ncc > 1) {
      Get_Normals_and_Gaps_at_Node(cnei_group, normals, gaps);
      const int node_index = cnei_group.Get_Constrained_Node_Index();
      Real  m_sn   = *(NODAL_MASS.Get_Scratch(node_index));
      if (ncc == 2) { 
	Real* free_dir = normals[2];
	ContactNodeEntityInteraction* cneis[2] ={NULL,NULL};
	int ncnei= 0;
        ContactNodeEntityInteraction *cnei0 = cnei_group.Get_Interaction(0);
        ContactNodeEntityInteraction *cnei1 = cnei_group.Get_Interaction(1);
	//
	// the tied interaction trumps all others
	//
	const int model_type1 = TDEnfModel(cnei0)->Interaction_Type(cnei0);
	const int model_type2 = TDEnfModel(cnei1)->Interaction_Type(cnei1);
        int indexes[2];
	if   ( (model_type1 == ContactTDEnfModel::TDEM_TIED ) && !(model_type2 == ContactTDEnfModel::TDEM_TIED ) ) {
	  cneis[0] = cnei0;
          indexes[0] = cnei_group.Get_Index(0);
	  ncnei = 1;
	} else if ( (model_type2 == ContactTDEnfModel::TDEM_TIED ) && !(model_type1 == ContactTDEnfModel::TDEM_TIED ) ) {
	  cneis[0] = cnei1;
          indexes[0] = cnei_group.Get_Index(1);
          ncnei =1;
	} else {
	  cneis[0] = cnei0; 
          cneis[1] = cnei1; 
          indexes[0] = cnei_group.Get_Index(0);
          indexes[1] = cnei_group.Get_Index(1);
          ncnei =2;
	}
	for(int n = 0; n < ncnei ; ++n) {
	  ContactNodeEntityInteraction *cnei = cneis[n];
          int interaction_index = indexes[n];
	  Real* tang_slip = NEI_TAN_SLIP.Get_Scratch(interaction_index);
	  Real free_slip[3];
	  for (int j=0;j<NDIM;  ++j ){free_slip[j] = tang_slip[j];}
	  //
	  // (1) Project on to unconstrained direction
	  //
	  Real dot = 0.0;
	  for (int j=0;j<NDIM;  ++j ){dot += free_slip[j]*free_dir[j];}
	  for (int j=0;j<NDIM;  ++j ){free_slip[j] = dot*free_dir[j];}
	  //
	  // (2) Constitutive Limit & Partition to the Face
	  //
          Real tot_eff_mass = 0.0;
	  //
	  // to emulate the predictor phase
	  //
          const int num_nodes_per_face = cnei_group.Get_Num_Face_Nodes(n);
          for(int j=0 ; j<num_nodes_per_face ; ++j ){
            const int node_face_index = cnei_group.Get_Face_Node_Index(n, j);
            const Real face_shape_function = cnei_group.Get_Face_Shape_Function(n, j);
            tot_eff_mass += face_shape_function*m_sn / *(NODAL_MASS.Get_Scratch(node_face_index));
          }
          Real* slip      = NEI_TOTAL_SLIP.Get_Scratch(interaction_index);
	  Real* rel_disp  = NEI_REL_DISP.  Get_Scratch(interaction_index);
	  Real* f_est     = NEI_FORCE.     Get_Scratch(interaction_index);
	  Real kin_partition = *(NEI_KINEMATIC_PARTITION.Get_Scratch(interaction_index));
	  for(int k=0 ; k<NDIM ; ++k) { f_trial[k]  = f_est[k] ;}
	  Real pen =  m_sn*dt2/(1.0+tot_eff_mass);
	  for(int k=0 ; k<NDIM ; ++k) { f_trial[k] += -pen*free_slip[k];}
	  Real* nm      = cnei->Get_Normal();
	  Real  gap     = *(NEI_GAP_TO_ENFORCE.Get_Scratch(interaction_index));
	  Real  area_sn = cnei->Get_Node_Area();
	  if (area_sn <= 0.0) area_sn = 1.0; //protect against zero/neg.areas
	  TDEnfModel(cnei)->Limit_Force(cnei,gap,rel_disp,slip,nm,dt,area_sn,f_trial);
	  // take only the correction in the free direction
	  Real f_diff = 0.0;
	  for(int k=0 ; k<NDIM ; ++k) {
            f_diff += kin_partition*(f_trial[k]-f_est[k])*free_dir[k];
          }
	  // add  correction to both surfaces maintain mom. balance
	  // slave node
	  Real* a_increm = A_CONSTITUTIVE.Get_Scratch(node_index);
	  
	  for(int k=0 ; k<NDIM ; ++k ){
	    // difference with actual forces to form corrector
	    const Real a_inc = f_diff*free_dir[k]/m_sn;
	    a_increm[k] += a_inc;
	    f_est[k] += m_sn*a_inc;
	    //
	    // master surface
	    //
            const int num_face_nodes = cnei_group.Get_Num_Face_Nodes(n);
  	    for(int j=0 ;  j< num_face_nodes ; ++j ){
              const int face_node_index = cnei_group.Get_Face_Node_Index(n, j);
              const Real face_shape_function = cnei_group.Get_Face_Shape_Function(n, j);
	      Real* a_con = A_CONSTITUTIVE.Get_Scratch(face_node_index);
              //
	      // the portion of the total force for the master node j
              //
	      Real partition = face_shape_function*m_sn / *(NODAL_MASS.Get_Scratch(face_node_index));
	      a_con[k] -= partition*a_inc;
      	    }
	  }
	}
      }
    }
  }  
#ifdef CONTACT_DEBUG_NODE
  for(int i=0 ; i<Number_Debug_Nodes() ; ++i ){
    ContactNode<Real>* dn = Debug_Node( i );
    if( dn ){
      const int node_index = dn->EnfArrayIndex();
      Real  m   = *(NODAL_MASS.    Get_Scratch(node_index));
      Real* ai =    A_CONSTITUTIVE.Get_Scratch(node_index);
      postream << "Constitut.Force for Node "
	       << dn->Exodus_ID() << " = " 
	       << m*ai[0] << " " << m*ai[1] << " " << m*ai[2] << "\n";
    }
  }
  postream.flush();
#endif
}



//
//  NKC, perf, switch the variable handle to scratch vars to avoid the nodal look ups, this routine called many a time...
//
void ContactTDEnforcement::Make_Scratch_Vector_Consistent_with_BCs(ScratchVariable &var )
{
  for( int i=0 ; i<number_of_total_nodes ; ++i ){
    int   nkc    = static_cast<int>(*(NUM_KIN_CONSTR.Get_Scratch(i)));
    Real* c_vec  = KIN_CONSTR_VECTOR.Get_Scratch(i);
    Real* vec    = var.Get_Scratch(i);
    Real dot;
    switch (nkc) {
      //case 0 : // no modification
    case 1 : // here 'c_vec' is the constrained direction
      dot = c_vec[0]*vec[0] + c_vec[1]*vec[1] + c_vec[2]*vec[2];
      vec[0] -= dot * c_vec[0] ;
      vec[1] -= dot * c_vec[1] ;
      vec[2] -= dot * c_vec[2] ;
      break;
    case 2 : // here 'c_vec' is the unconstrained direction
      dot = c_vec[0]*vec[0] + c_vec[1]*vec[1] + c_vec[2]*vec[2] ;
      vec[0] = dot * c_vec[0] ;
      vec[1] = dot * c_vec[1] ;
      vec[2] = dot * c_vec[2] ;
      break;
    case 3 : // no contact enforcement possible
      vec[0] = 0.0 ;
      vec[1] = 0.0 ;
      vec[2] = 0.0 ;
      break;
    }
  }
}

void ContactTDEnforcement::Get_Normals_and_Gaps_at_Node(const Node_Constraint_Group cnei_group, Real** normals, Real* gaps)
{
  int ncc = cnei_group.Num_Interactions();
  for(int j=0; j < ncc ; ++j) {
    ContactNodeEntityInteraction *cnei = cnei_group.Get_Interaction(j);
    int index = cnei_group.Get_Index(j);

    normals[j]  = cnei->Get_Normal();
    gaps[j] = *(NEI_GAP_TO_ENFORCE.Get_Scratch(index));

    if (gaps[j] > 0.0) gaps[j] = 0.0 ;
    Real* n_dir = normals[j];

    const int model_type = TDEnfModel(cnei)->Interaction_Type(cnei);
    if (model_type == ContactTDEnfModel::TDEM_SPRINGY){
      Real* slip     = NEI_TOTAL_SLIP.Get_Scratch(index);
      Real* rel_disp = NEI_REL_DISP  .Get_Scratch(index);
      Real s_n = 0.0;
      for( int k=0 ; k<NDIM ; ++k ){ s_n += slip[k]*n_dir[k];}
      Real d_n = 0.0;
      for( int k=0 ; k<NDIM ; ++k ){ d_n += rel_disp[k]*n_dir[k];}
      gaps[j]  = (s_n < 0.0 && d_n < 0.0) ? s_n : 0.0;
    }
    else if (model_type == ContactTDEnfModel::TDEM_TIED ) {
      Real* slip = NEI_TOTAL_SLIP.Get_Scratch(index);
      Real s_n = 0.0;
      for( int k=0 ; k<NDIM ; ++k ){ s_n += slip[k]*n_dir[k];}
      gaps[j] = s_n;
    }

    // multiply gaps by kp to get amount to penalize
    gaps[j] *= *(NEI_KINEMATIC_PARTITION.Get_Scratch(index));
  }
  if ( ncc == 2) { // number of contact constraints
     Real* n_1 = normals[0];
     Real* n_2 = normals[1];
     normals[2] = third_vector; // this is scratch local to CTDE
     Real* n_3 =  normals[2];
     n_3[0] = n_1[1]*n_2[2] - n_1[2]*n_2[1] ;
     n_3[1] = n_1[2]*n_2[0] - n_1[0]*n_2[2] ;
     n_3[2] = n_1[0]*n_2[1] - n_1[1]*n_2[0] ;
     Real scale = 1.0/std::sqrt(n_3[0]*n_3[0]+n_3[1]*n_3[1]+n_3[2]*n_3[2]);
     n_3[0] = n_3[0]*scale ;
     n_3[1] = n_3[1]*scale ;
     n_3[2] = n_3[2]*scale ;
     gaps[2] = 0.0 ;
  }
}




   
//
//  NKC, tag for performace!  
//  There are better, much faster, ways to do this solve than using matrix inversion inversion!
//
//  Should make a generic 3x3 Ax = b matrix solver
//
void ContactTDEnforcement::Unify(Real** normals, Real* magnitudes, Real* vector){ // Solve  g . n_i = g_i for g
#if defined(USE_EIGEN3) && defined(PJSA_FIX)
  Real* n_1 = normals[0];
  Real* n_2 = normals[1];
  Real* n_3 = normals[2];
  Eigen::Matrix<Real,3,3> A;
  A << n_1[0], n_1[1], n_1[2],
       n_2[0], n_2[1], n_2[2],
       n_3[0], n_3[1], n_3[2];
  Eigen::Map<Eigen::Matrix<Real,3,1> > b(magnitudes);
  Eigen::Map<Eigen::Matrix<Real,3,1> > x(vector);
  x = A.colPivHouseholderQr().solve(b);
#else
  Real* n_1 = normals[0];
  Real* n_2 = normals[1];
  Real* n_3 = normals[2]; 
  Real unify_matrix [3] [3] ;
  Real partition_matrix  [3] [3] ;
  partition_matrix[0] [0] = n_1[0] ;
  partition_matrix[0] [1] = n_1[1] ;
  partition_matrix[0] [2] = n_1[2] ;
  partition_matrix[1] [0] = n_2[0] ;
  partition_matrix[1] [1] = n_2[1] ;
  partition_matrix[1] [2] = n_2[2] ;
  partition_matrix[2] [0] = n_3[0] ; 
  partition_matrix[2] [1] = n_3[1] ;
  partition_matrix[2] [2] = n_3[2] ; 

  Invert_3x3Matrix(partition_matrix,unify_matrix);

  vector[0] = unify_matrix[0] [0]*magnitudes [0] 
            + unify_matrix[0] [1]*magnitudes [1]
            + unify_matrix[0] [2]*magnitudes [2] ;
  vector[1] = unify_matrix[1] [0]*magnitudes [0] 
            + unify_matrix[1] [1]*magnitudes [1]
            + unify_matrix[1] [2]*magnitudes [2] ;
  vector[2] = unify_matrix[2] [0]*magnitudes [0] 
            + unify_matrix[2] [1]*magnitudes [1]
            + unify_matrix[2] [2]*magnitudes [2] ;
#endif
}

//
//  What is this routine doing, can it be combined with partition Force?
//
void ContactTDEnforcement::Partition_Mass(int ncc, Real** normals, Real* f_sn, Real* partition_factors) {
#ifdef CONTACT_DEBUG_NODE
  ContactParOStream& postream = ParOStream();
  bool PRINT_THIS_NODE = false;
#endif
  PRECONDITION( ncc == 2 || ncc == 3 );
  switch( ncc ){
  case 2: {
    Real* n0 = normals[0];
    Real* n1 = normals[1];
    Real dot = n0[0]*n1[0] + n0[1]*n1[1] + n0[2]*n1[2];
    Real denom = 1.0 - dot*dot;
    // Recover the force at the slave node from the penalty phase
    //   * This give the direction for the momentum balance
    Real  f_sn_mag = std::sqrt( f_sn[0]*f_sn[0] + f_sn[1]*f_sn[1] +
				f_sn[2]*f_sn[2] );
    //
    // Get the amount of f_sn along each of the pushback directions
    //
    Real f_n0 = f_sn[0]*n0[0] + f_sn[1]*n0[1] + f_sn[2]*n0[2] ; 
    Real f_n1 = f_sn[0]*n1[0] + f_sn[1]*n1[1] + f_sn[2]*n1[2] ; 
    Real a0,a1; // The distribution "factors" for the mass
#ifdef PJSA_FIX
    // Shouldn't the factors a0 and a1 always add up to 1?
    a0 = (     f_n0 - dot*f_n1 );
    a1 = (-dot*f_n0 +     f_n1 );
    if(a0 != 0 && a1 != 0) {
      Real inv_factors = 1/sqrt(a0*a0+a1*a1);
      a0 *= inv_factors;
      a1 *= inv_factors;
#else
    if( denom > 0.0  &&  f_sn_mag > 0.0 ){
      Real inv_factors = (1.0/denom) * (1.0/f_sn_mag);
      a0 = (     f_n0 - dot*f_n1 ) * inv_factors;
      a1 = (-dot*f_n0 +     f_n1 ) * inv_factors;
#endif
      a0 *= a0;
      a1 *= a1;
    } else {
      a0 = 0.5;
      a1 = 0.5;
    }
#ifdef CONTACT_DEBUG_NODE
    if( PRINT_THIS_NODE )
      postream << "Mass partitioning for node " 
	//		 << node->Exodus_ID() 
	       << " a0=" << a0 << " a1=" << a1 << " sum= " << a0+a1 << "\n";
#endif
    partition_factors[0] = a0;
    partition_factors[1] = a1;
    partition_factors[2] = 0.0;
    break;
  }
  case 3: {
    Real* n0 = normals[0];
    Real* n1 = normals[1];
    Real* n2 = normals[2];
    // Recover the force at the slave node from the penalty phase
    //   * This give the direction for the momentum balance
    Real  f_sn_mag = std::sqrt( f_sn[0]*f_sn[0] + f_sn[1]*f_sn[1] +
				f_sn[2]*f_sn[2] );
    // Get the amount of f_sn along each of the pushback directions
    Real f_n0 = f_sn[0]*n0[0] + f_sn[1]*n0[1] + f_sn[2]*n0[2] ; 
    Real f_n1 = f_sn[0]*n1[0] + f_sn[1]*n1[1] + f_sn[2]*n1[2];
    Real f_n2 = f_sn[0]*n2[0] + f_sn[1]*n2[1] + f_sn[2]*n2[2];
    // Compute the angle between the interaction directions
    Real n0dn1 = n0[0]*n1[0] + n0[1]*n1[1] + n0[2]*n1[2];
    Real n0dn2 = n0[0]*n2[0] + n0[1]*n2[1] + n0[2]*n2[2];
    Real n1dn2 = n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2];
    Real A[3][3],A_inv[3][3];
    A[0][0] = 1.0;
    A[0][1] = n0dn1;
    A[0][2] = n0dn2;
    A[1][0] = n0dn1;
    A[1][1] = 1.0;
    A[1][2] = n1dn2;
    A[2][0] = n0dn2;
    A[2][1] = n1dn2;
    A[2][2] = 1.0;
    Real determinant = Invert_3x3Matrix( A, A_inv );
    Real a0,a1,a2;
    if( determinant <= ZERO_TOL ){
      // The constraints are co-aligned
      a0 = 1.0/3.0;
      a1 = 1.0/3.0;
      a2 = 1.0/3.0;
    } else {
      a0 = A_inv[0][0]*f_n0 + A_inv[0][1]*f_n1 + A_inv[0][2]*f_n2; 
      a1 = A_inv[1][0]*f_n0 + A_inv[1][1]*f_n1 + A_inv[1][2]*f_n2; 
      a2 = A_inv[2][0]*f_n0 + A_inv[2][1]*f_n1 + A_inv[2][2]*f_n2; 
#ifdef PJSA_FIX
      if( a0 != 0 && a1 !=0 && a2 != 0 ){
        // Shouldn't the factors a0 and a1 always add up to 1?
        Real f_sn_mag_inv = 1.0/sqrt(a0*a0+a1*a1+a2*a2);
#else
      if( f_sn_mag > 0.0 ){
        Real f_sn_mag_inv = 1.0/f_sn_mag;
#endif
	a0 *= f_sn_mag_inv;
	a1 *= f_sn_mag_inv;
	a2 *= f_sn_mag_inv;
	a0 = a0*a0;
	a1 = a1*a1;
	a2 = a2*a2;
      } else {
	a0 = 1.0/3.0;
	a1 = 1.0/3.0;
	a2 = 1.0/3.0;
      }
    }
    partition_factors[0] = a0;
    partition_factors[1] = a1;
    partition_factors[2] = a2;
    break;
  }
  }
}

//
// Solve (p n)_i = F
//
void ContactTDEnforcement::Partition_Force (int ncc, 
                                            Real** normals,  
                                            Real* total_force, 
                                            Real* surface_stiffness_inv,
                                            Real f_constr[3][3] )
{ 

#if defined(USE_EIGEN3) && defined(PJSA_FIX)
  Real* n_1 = normals[0];
  Real* n_2 = normals[1];
  Real* n_3 = normals[2];
  Eigen::Matrix<Real,3,3> N;
  N << n_1[0], n_2[0], n_3[0],
       n_1[1], n_2[1], n_3[1],
       n_1[2], n_2[2], n_3[2];
  Eigen::Map<Eigen::Matrix<Real,3,1> > F(total_force);
  Eigen::Map<Eigen::Matrix<Real,3,1> > f_constr0(f_constr[0]), f_constr1(f_constr[1]), f_constr2(f_constr[2]);

  if(ncc == 2) {
    Eigen::Matrix<Real,3,3> J = N.transpose()*(N.col(0)*N.col(0).transpose()-N.col(1)*N.col(1).transpose());
    Eigen::Matrix<Real,3,1> f_dot_mags = N.transpose()*(F - F.dot(N.col(0))*N.col(0));

    f_constr2.setZero();
    f_constr1 = (-J.colPivHouseholderQr().solve(f_dot_mags)).dot(N.col(1))*N.col(1);
    f_constr0 = (F-f_constr1).dot(N.col(0))*N.col(0);
  }
  else if(ncc == 3) {

    // first, partition total_force into f_constr2 and f_constr01
    Eigen::Matrix<Real,3,3> P01 = N.leftCols(2)*(N.leftCols(2).transpose()*N.leftCols(2)).colPivHouseholderQr().solve(N.leftCols(2).transpose());
    Eigen::Matrix<Real,3,3> J = N.transpose()*(P01-N.col(2)*N.col(2).transpose());
    Eigen::Matrix<Real,3,1> f_dot_mags = N.transpose()*(Eigen::Matrix3d::Identity()-P01)*F;

    f_constr2 = (-J.colPivHouseholderQr().solve(f_dot_mags)).dot(N.col(2))*N.col(2);
    Eigen::Matrix<Real,3,1> f_constr01 = P01*(F-f_constr2);

    {
      // second, partition f_constr01 into f_constr0 and f_constr1
      Eigen::Matrix<Real,3,3> J = N.transpose()*(N.col(0)*N.col(0).transpose()-N.col(1)*N.col(1).transpose());
      Eigen::Matrix<Real,3,1> f_dot_mags = N.transpose()*(f_constr01 - f_constr01.dot(N.col(0))*N.col(0));

      f_constr1 = (-J.colPivHouseholderQr().solve(f_dot_mags)).dot(N.col(1))*N.col(1);
      f_constr0 = (f_constr01-f_constr1).dot(N.col(0))*N.col(0);
    }
  }
#else
#ifdef CONTACT_DEBUG_NODE
  ContactParOStream& postream = ParOStream();
  bool PRINT_THIS_NODE = 0;
  //  if( PRINT_THIS_NODE )
  //    postream << "Partitioning Force for node " << node->Exodus_ID() << "\n";
#endif

  int num_constraints = ncc;
  PRECONDITION( num_constraints == 2 || num_constraints == 3 );
  
#ifdef CONTACT_DEBUG_NODE
  if( PRINT_THIS_NODE )
    postream << "    Inverse Surface Stiffnesses = " << surface_stiffness_inv[0]
	     << " " << surface_stiffness_inv[1] << " " 
	     << surface_stiffness_inv[2] << "\n";
#endif

  Real k_bar  =1.0;
  Real k_tilda=1.0;
  // the constraints are ordered such that NF come before NS
  // so if the first constraint has a zero inverse stiffness they all do
  // this allows us to avoid pivoting the constraints
  if(surface_stiffness_inv[0] != 0.0 ){
    k_bar = 1.0 + surface_stiffness_inv[1]/surface_stiffness_inv[0];
    if( num_constraints == 3 && surface_stiffness_inv[1] != 0.0 ){
      k_tilda = 1.0 + (surface_stiffness_inv[1]*surface_stiffness_inv[2] + 
		       surface_stiffness_inv[0]*surface_stiffness_inv[2])/
	(surface_stiffness_inv[0]*surface_stiffness_inv[1]);
    }
#ifdef CONTACT_DEBUG_NODE
    if( PRINT_THIS_NODE )
      postream << "    K_Bar = " << k_bar << "   K_Tilda = " << k_tilda << "\n";
#endif
  }
  Real k_bar_inv = 1.0/k_bar;
  Real k_tilda_inv = 1.0/k_tilda;

  Real* n_0 = normals[0];
  Real* n_1 = normals[1];
  Real* n_2 = normals[2];

  Real partition_matrix  [3] [3] ;
  partition_matrix[0] [0] = n_0[0] ;
  partition_matrix[0] [1] = n_0[1] ;
  partition_matrix[0] [2] = n_0[2] ;
  partition_matrix[1] [0] = n_1[0] ;
  partition_matrix[1] [1] = n_1[1] ;
  partition_matrix[1] [2] = n_1[2] ;
  partition_matrix[2] [0] = n_2[0] ;      
  partition_matrix[2] [1] = n_2[1] ;
  partition_matrix[2] [2] = n_2[2] ; 

  Real partition_inverse [3] [3];
  Invert_3x3Matrix(partition_matrix,partition_inverse);

  Real f_dot_mags[3];
  f_dot_mags[0] = n_0[0]*total_force[0] + n_0[1]*total_force[1] + 
                  n_0[2]*total_force[2];
  f_dot_mags[1] = n_1[0]*total_force[0] + n_1[1]*total_force[1] + 
                  n_1[2]*total_force[2];
  f_dot_mags[2] = n_2[0]*total_force[0] + n_2[1]*total_force[1] + 
                  n_2[2]*total_force[2];

  if( num_constraints == 2 ){
    f_constr[1][0] = partition_inverse[0] [0]* f_dot_mags[0] 
                   + partition_inverse[0] [1]* f_dot_mags[1] 
                   + partition_inverse[0] [2]* f_dot_mags[2] ;
    f_constr[1][1] = partition_inverse[1] [0]* f_dot_mags[0] 
                   + partition_inverse[1] [1]* f_dot_mags[1] 
                   + partition_inverse[1] [2]* f_dot_mags[2] ;
    f_constr[1][2] = partition_inverse[2] [0]* f_dot_mags[0] 
                   + partition_inverse[2] [1]* f_dot_mags[1] 
                   + partition_inverse[2] [2]* f_dot_mags[2] ;
    f_constr[1][0] *= k_bar_inv;
    f_constr[1][1] *= k_bar_inv;
    f_constr[1][2] *= k_bar_inv;
    f_constr[0][0] = total_force[0] - f_constr[1][0];
    f_constr[0][1] = total_force[1] - f_constr[1][1];
    f_constr[0][2] = total_force[2] - f_constr[1][2];
    f_constr[2][0] = 0.0;
    f_constr[2][1] = 0.0;
    f_constr[2][2] = 0.0;
    Real total_force_magsqr = total_force[0]*total_force[0] +
			      total_force[1]*total_force[1] +
			      total_force[2]*total_force[2] ;
    Real reference_force = std::fabs(f_dot_mags[0]) + std::fabs(f_dot_mags[1]);
#ifdef CONTACT_DEBUG_NODE
    if( PRINT_THIS_NODE ){
      postream << "    Beginning iterative Solve with 2 constraints" << "\n";
      postream << "      Total Force     = " << total_force[0] << " "
	       << total_force[1] << " " << total_force[2] << "\n";
      postream << "      Total Force Mag**2 = " << total_force_magsqr << "\n";
      postream << "      F_Constr_0      = " << f_constr[0][0] << " "
	       << f_constr[0][1] << " " << f_constr[0][2] << "\n";
      postream << "      F_Constr_1      = " << f_constr[1][0] << " "
	       << f_constr[1][1] << " " << f_constr[1][2] << "\n";
    }
#endif

    int iter = 1;
    bool converged = false;
    while( !converged ) {
      Real f_constr1_dot_n_1 = f_constr[1][0]*n_1[0] + f_constr[1][1]*n_1[1]
	                     + f_constr[1][2]*n_1[2];
      f_constr[1][0] = f_constr1_dot_n_1*n_1[0];
      f_constr[1][1] = f_constr1_dot_n_1*n_1[1];
      f_constr[1][2] = f_constr1_dot_n_1*n_1[2];
      Real f_constr_temp[3];
      f_constr_temp[0] = total_force[0] - f_constr[1][0];
      f_constr_temp[1] = total_force[1] - f_constr[1][1];
      f_constr_temp[2] = total_force[2] - f_constr[1][2];
      Real f_constr0_dot_n_0 = f_constr_temp[0]*n_0[0] + f_constr_temp[1]*n_0[1]
	                     + f_constr_temp[2]*n_0[2];
      f_constr[0][0] = f_constr0_dot_n_0*n_0[0];
      f_constr[0][1] = f_constr0_dot_n_0*n_0[1];
      f_constr[0][2] = f_constr0_dot_n_0*n_0[2];
      // MWH : Kevin lets re-test this (I corrected a problem with it)
      Real denominator = total_force[0]*(f_constr[0][0] + f_constr[1][0]) + 
	                 total_force[1]*(f_constr[0][1] + f_constr[1][1]) +
                         total_force[2]*(f_constr[0][2] + f_constr[1][2]) ;
      if( denominator > ZERO_TOL ){
	Real alpha = total_force_magsqr/denominator;
	f_constr[0][0] *= alpha;
	f_constr[0][1] *= alpha;
	f_constr[0][2] *= alpha;
	f_constr[1][0] *= alpha;
	f_constr[1][1] *= alpha;
	f_constr[1][2] *= alpha;
      }
      Real unpartitioned_force[3];
      unpartitioned_force[0] = total_force[0] - f_constr[0][0] - f_constr[1][0];
      unpartitioned_force[1] = total_force[1] - f_constr[0][1] - f_constr[1][1];
      unpartitioned_force[2] = total_force[2] - f_constr[0][2] - f_constr[1][2];
      f_dot_mags[0] = n_0[0]*unpartitioned_force[0] + 
	              n_0[1]*unpartitioned_force[1] + 
	              n_0[2]*unpartitioned_force[2];
      f_dot_mags[1] = n_1[0]*unpartitioned_force[0] + 
                      n_1[1]*unpartitioned_force[1] +  
	              n_1[2]*unpartitioned_force[2];
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE ){
	postream << "      Data on Iteration " << iter << "\n";
	postream << "         F_Constr_0      = " << f_constr[0][0] << " "
		 << f_constr[0][1] << " " << f_constr[0][2] << "\n";
	postream << "         F_Constr_1      = " << f_constr[1][0] << " "
		 << f_constr[1][1] << " " << f_constr[1][2] << "\n";
	postream << "         Unpartitioned F = " << unpartitioned_force[0] 
		 << " " << unpartitioned_force[1] << " " 
		 << unpartitioned_force[2] << "\n";
	postream << "         F_dot_Mags      = " << f_dot_mags[0] << " "
		 << f_dot_mags[1] << "\n";
      }
#endif
      if( iter != PARTITION_ITERATIONS && 
	  (std::fabs(f_dot_mags[0]) + std::fabs(f_dot_mags[1]) > 
           std::max(ZERO_TOL,ZERO_TOL*reference_force) ) ){
	++iter;
	Real f_constr1_iter[3];
	f_constr1_iter[0] = partition_inverse[0] [0]* f_dot_mags[0] 
	  + partition_inverse[0] [1]* f_dot_mags[1] 
	  + partition_inverse[0] [2]* f_dot_mags[2] ;
	f_constr1_iter[1] = partition_inverse[1] [0]* f_dot_mags[0] 
	  + partition_inverse[1] [1]* f_dot_mags[1] 
	  + partition_inverse[1] [2]* f_dot_mags[2] ;
	f_constr1_iter[2] = partition_inverse[2] [0]* f_dot_mags[0] 
	  + partition_inverse[2] [1]* f_dot_mags[1] 
	  + partition_inverse[2] [2]* f_dot_mags[2] ;
	f_constr[1][0] += f_constr1_iter[0];
	f_constr[1][1] += f_constr1_iter[1];
	f_constr[1][2] += f_constr1_iter[2];
	f_constr[0][0] += (unpartitioned_force[0] - f_constr1_iter[0]);
	f_constr[0][1] += (unpartitioned_force[1] - f_constr1_iter[1]);
	f_constr[0][2] += (unpartitioned_force[2] - f_constr1_iter[2]);
      } else {
	converged = true;
      }
    }
#if CONTACT_DEBUG_PRINT_LEVEL>=2
    if (iter == PARTITION_ITERATIONS) {
      std::cerr << "ACME did not converge in Partition_Force" << std::endl;
      std::cerr << "  f_dot_mags[0]     = " << f_dot_mags[0] << std::endl;
      std::cerr << "  f_dot_mags[1]     = " << f_dot_mags[1] << std::endl;
      std::cerr << "  total_force_mag^2 = " << total_force_magsqr << std::endl;
    }
#endif

  } else if( num_constraints == 3 ){
    f_constr[2][0] = partition_inverse[0][0] * f_dot_mags[0] 
                   + partition_inverse[0][1] * f_dot_mags[1] 
                   + partition_inverse[0][2] * f_dot_mags[2] ;
    f_constr[2][1] = partition_inverse[1][0] * f_dot_mags[0] 
                   + partition_inverse[1][1] * f_dot_mags[1] 
                   + partition_inverse[1][2] * f_dot_mags[2] ;
    f_constr[2][2] = partition_inverse[2][0] * f_dot_mags[0] 
                   + partition_inverse[2][1] * f_dot_mags[1] 
                   + partition_inverse[2][2] * f_dot_mags[2] ;
    f_constr[2][0] *= k_tilda_inv;
    f_constr[2][1] *= k_tilda_inv;
    f_constr[2][2] *= k_tilda_inv;

    Real f_constr_remainder[3];
    f_constr_remainder[0] = f_dot_mags[0] - (f_constr[2][0]*n_0[0] + 
                                             f_constr[2][1]*n_0[1] + 
                                             f_constr[2][2]*n_0[2] ); 
    f_constr_remainder[1] = f_dot_mags[1] - (f_constr[2][0]*n_1[0] + 
                                             f_constr[2][1]*n_1[1] + 
                                             f_constr[2][2]*n_1[2] ); 
    f_constr_remainder[2] = f_dot_mags[2] - (f_constr[2][0]*n_2[0] + 
                                             f_constr[2][1]*n_2[1] + 
                                             f_constr[2][2]*n_2[2] ); 

    f_constr[1][0] = partition_inverse[0][0] * f_constr_remainder[0] 
                   + partition_inverse[0][1] * f_constr_remainder[1] 
                   + partition_inverse[0][2] * f_constr_remainder[2];
    f_constr[1][1] = partition_inverse[1][0] * f_constr_remainder[0] 
                   + partition_inverse[1][1] * f_constr_remainder[1] 
                   + partition_inverse[1][2] * f_constr_remainder[2];
    f_constr[1][2] = partition_inverse[2][0] * f_constr_remainder[0] 
                   + partition_inverse[2][1] * f_constr_remainder[1] 
                   + partition_inverse[2][2] * f_constr_remainder[2];
    f_constr[1][0] *= k_bar_inv;
    f_constr[1][1] *= k_bar_inv;
    f_constr[1][2] *= k_bar_inv;

    f_constr[0][0] = total_force[0] - f_constr[2][0] - f_constr[1][0];
    f_constr[0][1] = total_force[1] - f_constr[2][1] - f_constr[1][1];
    f_constr[0][2] = total_force[2] - f_constr[2][2] - f_constr[1][2];

    Real reference_force = std::fabs(f_dot_mags[0]) + 
                           std::fabs(f_dot_mags[1]) + 
                           std::fabs(f_dot_mags[2]);
    Real total_force_magsqr = total_force[0]*total_force[0] +
	                      total_force[1]*total_force[1] + 
	                      total_force[2]*total_force[2] ;
#ifdef CONTACT_DEBUG_NODE
    if( PRINT_THIS_NODE ){
      postream << "    Beginning iterative Solve with 3 constraints" << "\n";
      postream << "      Total Force     = " << total_force[0] << " "
	            << total_force[1] << " " << total_force[2] << "\n";
      postream << "      Total Force Mag**2 = " << total_force_magsqr << "\n";
      postream << "      F_Constr_0      = " << f_constr[0][0] << " "
	       << f_constr[0][1] << " " << f_constr[0][2] << "\n";
      postream << "      F_Constr_1      = " << f_constr[1][0] << " "
	       << f_constr[1][1] << " " << f_constr[1][2] << "\n";
      postream << "      F_Constr_2      = " << f_constr[2][0] << " "
               << f_constr[2][1] << " " << f_constr[2][2] << "\n";
    }
#endif

    int iter = 1;
    bool converged = false;
    while( !converged ) {
      Real f_constr2_dot_n_2 = f_constr[2][0]*n_2[0]
                             + f_constr[2][1]*n_2[1]
	                     + f_constr[2][2]*n_2[2];      
      Real f_constr_temp[3];
      //amount of f_constr[2] that is not supported by normal_2 force
      //later applied to normal_1
      f_constr_temp[0] = f_constr[2][0] - f_constr2_dot_n_2*n_2[0];
      f_constr_temp[1] = f_constr[2][1] - f_constr2_dot_n_2*n_2[1];
      f_constr_temp[2] = f_constr[2][2] - f_constr2_dot_n_2*n_2[2];
      //normal force on surface n_2
      f_constr[2][0] = f_constr2_dot_n_2*n_2[0];
      f_constr[2][1] = f_constr2_dot_n_2*n_2[1];
      f_constr[2][2] = f_constr2_dot_n_2*n_2[2];

      Real f_constr_temp2[3];
      f_constr_temp2[0] = f_constr_temp[0]*n_0[0] + 
                          f_constr_temp[1]*n_0[1] +
      	                  f_constr_temp[2]*n_0[2] ;
      f_constr_temp2[1] = f_constr_temp[0]*n_1[0] + 
                          f_constr_temp[1]*n_1[1] +
	                  f_constr_temp[2]*n_1[2] ;
      f_constr_temp2[2] = f_constr_temp[0]*n_2[0] + 
                          f_constr_temp[1]*n_2[1] +
	                  f_constr_temp[2]*n_2[2] ;
      //re-calculation of f_constr[1] (adding in part of f_constr[2] that
      //could not be supported by nomal force on n_2)
      f_constr[1][0] += (partition_inverse[0][0] * f_constr_temp2[0] 
                       + partition_inverse[0][1] * f_constr_temp2[1] 
                       + partition_inverse[0][2] * f_constr_temp2[2])
	               * k_bar_inv ;
      f_constr[1][1] += (partition_inverse[1][0] * f_constr_temp2[0] 
                       + partition_inverse[1][1] * f_constr_temp2[1] 
                       + partition_inverse[1][2] * f_constr_temp2[2])
                       * k_bar_inv ;
      f_constr[1][2] += (partition_inverse[2][0] * f_constr_temp2[0] 
                       + partition_inverse[2][1] * f_constr_temp2[1] 
                       + partition_inverse[2][2] * f_constr_temp2[2])
                       * k_bar_inv ;

      Real f_constr1_dot_n_1 = f_constr[1][0]*n_1[0]
                             + f_constr[1][1]*n_1[1]
	                     + f_constr[1][2]*n_1[2];
      f_constr[1][0] = f_constr1_dot_n_1*n_1[0];
      f_constr[1][1] = f_constr1_dot_n_1*n_1[1];
      f_constr[1][2] = f_constr1_dot_n_1*n_1[2];

      //re-calculation of f_constr[0]
      f_constr[0][0] = total_force[0] - f_constr[2][0] - f_constr[1][0];
      f_constr[0][1] = total_force[1] - f_constr[2][1] - f_constr[1][1];
      f_constr[0][2] = total_force[2] - f_constr[2][2] - f_constr[1][2];

      Real f_constr0_dot_n_0 = f_constr[0][0]*n_0[0]
                             + f_constr[0][1]*n_0[1]
	                     + f_constr[0][2]*n_0[2];
      f_constr[0][0] = f_constr0_dot_n_0*n_0[0];
      f_constr[0][1] = f_constr0_dot_n_0*n_0[1];
      f_constr[0][2] = f_constr0_dot_n_0*n_0[2];

      Real denominator = total_force[0]*(f_constr[0][0]+
                                         f_constr[1][0]+
                                         f_constr[2][0]) + 
	                 total_force[1]*(f_constr[0][1]+
                                         f_constr[1][1]+
                                         f_constr[2][1]) +
                         total_force[2]*(f_constr[0][2]+
                                         f_constr[1][2]+
                                         f_constr[2][2]) ;

      if( denominator > ZERO_TOL ){
        Real alpha = total_force_magsqr/denominator;
        f_constr[0][0] *= alpha;
        f_constr[0][1] *= alpha;
        f_constr[0][2] *= alpha;
        f_constr[1][0] *= alpha;
        f_constr[1][1] *= alpha;
        f_constr[1][2] *= alpha;
        f_constr[2][0] *= alpha;
        f_constr[2][1] *= alpha;
        f_constr[2][2] *= alpha;
      }

      Real unpartitioned_force[3];
      unpartitioned_force[0] = total_force[0] - f_constr[0][0]
                                              - f_constr[1][0]
                                              - f_constr[2][0];
      unpartitioned_force[1] = total_force[1] - f_constr[0][1]
                                              - f_constr[1][1]
                                              - f_constr[2][1];
      unpartitioned_force[2] = total_force[2] - f_constr[0][2]
                                              - f_constr[1][2]
                                              - f_constr[2][2];
      f_dot_mags[0] = n_0[0]*unpartitioned_force[0] + 
	              n_0[1]*unpartitioned_force[1] + 
	              n_0[2]*unpartitioned_force[2];
      f_dot_mags[1] = n_1[0]*unpartitioned_force[0] + 
                      n_1[1]*unpartitioned_force[1] +  
	              n_1[2]*unpartitioned_force[2];
      f_dot_mags[2] = n_2[0]*unpartitioned_force[0] + 
                      n_2[1]*unpartitioned_force[1] +  
	              n_2[2]*unpartitioned_force[2];
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE ){
	postream << "      Data on Iteration " << iter << "\n";
	postream << "         F_Constr_0      = " << f_constr[0][0] << " "
		 << f_constr[0][1] << " " << f_constr[0][2] << "\n";
	postream << "         F_Constr_1      = " << f_constr[1][0] << " "
		 << f_constr[1][1] << " " << f_constr[1][2] << "\n";
	postream << "         F_Constr_2      = " << f_constr[2][0] << " "
		 << f_constr[2][1] << " " << f_constr[2][2] << "\n";
	postream << "         Unpartitioned F = " << unpartitioned_force[0] 
		 << " " << unpartitioned_force[1] << " " 
		 << unpartitioned_force[2] << "\n";
	postream << "         F_dot_Mags      = " << f_dot_mags[0] << " "
		 << f_dot_mags[1] << " " << f_dot_mags[2] << "\n";
      }
#endif
      if( iter != PARTITION_ITERATIONS && 
	  (std::fabs(f_dot_mags[0])+std::fabs(f_dot_mags[1])+std::fabs(f_dot_mags[2]) > 
           std::max(ZERO_TOL,ZERO_TOL*reference_force) ) ){
	++iter;
	Real f_constr2_iter[3];
	f_constr2_iter[0] = (partition_inverse[0][0] * f_dot_mags[0] 
	                   + partition_inverse[0][1] * f_dot_mags[1] 
	                   + partition_inverse[0][2] * f_dot_mags[2])*k_tilda_inv ;
	f_constr2_iter[1] = (partition_inverse[1][0] * f_dot_mags[0] 
                           + partition_inverse[1][1] * f_dot_mags[1] 
	                   + partition_inverse[1][2] * f_dot_mags[2])*k_tilda_inv ;
	f_constr2_iter[2] = (partition_inverse[2][0] * f_dot_mags[0] 
                           + partition_inverse[2][1] * f_dot_mags[1] 
                           + partition_inverse[2][2] * f_dot_mags[2])*k_tilda_inv ;
	f_constr[2][0] += f_constr2_iter[0];
	f_constr[2][1] += f_constr2_iter[1];
	f_constr[2][2] += f_constr2_iter[2];


        f_constr_temp[0] = f_dot_mags[0] - (f_constr2_iter[0]*n_0[0] + 
                                            f_constr2_iter[1]*n_0[1] + 
                                            f_constr2_iter[2]*n_0[2] ); 
        f_constr_temp[1] = f_dot_mags[1] - (f_constr2_iter[0]*n_1[0] + 
                                            f_constr2_iter[1]*n_1[1] + 
                                            f_constr2_iter[2]*n_1[2] ); 
        f_constr_temp[2] = f_dot_mags[2] - (f_constr2_iter[0]*n_2[0] + 
                                            f_constr2_iter[1]*n_2[1] + 
                                            f_constr2_iter[2]*n_2[2] ); 

	Real f_constr1_iter[3];
        f_constr1_iter[0] = (partition_inverse[0][0]* f_constr_temp[0] 
                           + partition_inverse[0][1]* f_constr_temp[1] 
                           + partition_inverse[0][2]* f_constr_temp[2])*k_bar_inv ;
        f_constr1_iter[1] = (partition_inverse[1][0]* f_constr_temp[0] 
                           + partition_inverse[1][1]* f_constr_temp[1] 
                           + partition_inverse[1][2]* f_constr_temp[2])*k_bar_inv ;
        f_constr1_iter[2] = (partition_inverse[2][0]* f_constr_temp[0] 
                           + partition_inverse[2][1]* f_constr_temp[1] 
                           + partition_inverse[2][2]* f_constr_temp[2])*k_bar_inv ;
        f_constr[1][0] += f_constr1_iter[0];
        f_constr[1][1] += f_constr1_iter[1];
	f_constr[1][2] += f_constr1_iter[2];

	Real f_constr0_iter[3];
        f_constr0_iter[0] = unpartitioned_force[0] - f_constr2_iter[0]
                                                   - f_constr1_iter[0];
        f_constr0_iter[1] = unpartitioned_force[1] - f_constr2_iter[1]
                                                   - f_constr1_iter[1];
        f_constr0_iter[2] = unpartitioned_force[2] - f_constr2_iter[2]
                                                   - f_constr1_iter[2];
        f_constr[0][0] += f_constr0_iter[0];
        f_constr[0][1] += f_constr0_iter[1];
	f_constr[0][2] += f_constr0_iter[2];


      } else {
	converged = true;
      }
    }
#if CONTACT_DEBUG_PRINT_LEVEL>=2
    if (iter == PARTITION_ITERATIONS) {
      std::cerr << "ACME did not converge in Partition_Force" << std::endl;
      std::cerr << "  f_dot_mags[0]     = " << f_dot_mags[0] << std::endl;
      std::cerr << "  f_dot_mags[1]     = " << f_dot_mags[1] << std::endl;
      std::cerr << "  f_dot_mags[2]     = " << f_dot_mags[2] << std::endl;
      std::cerr << "  total_force_mag^2 = " << total_force_magsqr << std::endl;
    }
#endif
  }
#endif
}


void 
ContactTDEnforcement::Partition_Gap
( Real** normals, Real* magnitudes, Real* vector)
{ // Solve  [n_i]T {g_i} = {g}  for g_i by least squares
  // leastsquares_matrix = [n_i][n_i]T (formed directly)
#if defined(USE_EIGEN3) && defined(PJSA_FIX)
  Real* n_1 = normals[0];
  Real* n_2 = normals[1];
  Real* n_3 = normals[2];
  Eigen::Matrix<Real,3,3> A;
  A << n_1[0], n_1[1], n_1[2],
       n_2[0], n_2[1], n_2[2],
       n_3[0], n_3[1], n_3[2];
  Eigen::Map<Eigen::Matrix<Real,3,1> > x(magnitudes);
  Eigen::Map<Eigen::Matrix<Real,3,1> > b(vector);
  x = A.transpose().colPivHouseholderQr().solve(b);
#else
  Real leastsquares_matrix [3] [3] ;
  Real partition_matrix  [3] [3] ;
  Real* n_1 = normals[0];
  Real* n_2 = normals[1];
  Real* n_3 = normals[2];
  Real n1dn2 = 0.0;
  Real n1dn3 = 0.0;
  Real n2dn3 = 0.0;
  for(int j=0 ; j<NDIM ; ++j ){
    n1dn2 += n_1[j]*n_2[j];
    n1dn3 += n_1[j]*n_3[j];
    n2dn3 += n_2[j]*n_3[j];
  }

  leastsquares_matrix[0] [0] = 1.0 ;
  leastsquares_matrix[0] [1] = n1dn2 ;
  leastsquares_matrix[0] [2] = n1dn3 ;
  leastsquares_matrix[1] [0] = n1dn2 ;
  leastsquares_matrix[1] [1] = 1.0 ;
  leastsquares_matrix[1] [2] = n2dn3 ;
  leastsquares_matrix[2] [0] = n1dn3 ; 
  leastsquares_matrix[2] [1] = n2dn3 ;
  leastsquares_matrix[2] [2] = 1.0 ; 

  Invert_3x3Matrix(leastsquares_matrix,partition_matrix);

  Real vdn1 = 0.0;
  Real vdn2 = 0.0;
  Real vdn3 = 0.0;
  for(int j=0 ; j<NDIM ; ++j ){
     vdn1 += n_1[j]*vector[j];
     vdn2 += n_2[j]*vector[j];
     vdn3 += n_3[j]*vector[j];
  }

  magnitudes[0] = partition_matrix[0] [0]*vdn1
                + partition_matrix[0] [1]*vdn2
                + partition_matrix[0] [2]*vdn3;
  magnitudes[1] = partition_matrix[1] [0]*vdn1
                + partition_matrix[1] [1]*vdn2
                + partition_matrix[1] [2]*vdn3;
  magnitudes[2] = partition_matrix[2] [0]*vdn1
                + partition_matrix[2] [1]*vdn2
                + partition_matrix[2] [2]*vdn3;
#endif
}



Real ContactTDEnforcement::Invert_3x3Matrix(Real M [3][3], Real M_inv [3][3])
{ 
  //
  // Cofactors
  //
  Real cofM [3] [3];
  cofM[0][0] =   M[1][1]*M[2][2] - M[1][2]*M[2][1] ;
  cofM[1][0] = -(M[0][1]*M[2][2] - M[0][2]*M[2][1]) ;
  cofM[2][0] =   M[0][1]*M[1][2] - M[0][2]*M[1][1] ;
  cofM[0][1] = -(M[1][0]*M[2][2] - M[1][2]*M[2][0]) ;
  cofM[1][1] =   M[0][0]*M[2][2] - M[0][2]*M[2][0] ;
  cofM[2][1] = -(M[0][0]*M[1][2] - M[0][2]*M[1][0]) ;
  cofM[0][2] =   M[1][0]*M[2][1] - M[1][1]*M[2][0] ;
  cofM[1][2] = -(M[0][0]*M[2][1] - M[0][1]*M[2][0]) ;
  cofM[2][2] =   M[0][0]*M[1][1] - M[0][1]*M[1][0] ;
  //
  // Determinant
  //
  Real det = M[0][0]*cofM[0][0] + M[0][1]*cofM[0][1] + M[0][2]*cofM[0][2] ;
  if( std::fabs(det) < DET_TOL ) {
    return 0;
  }
  Real det_inv = 1.0/det  ;
  M_inv[0][0] = cofM[0][0]*det_inv ;
  M_inv[1][0] = cofM[0][1]*det_inv ;
  M_inv[2][0] = cofM[0][2]*det_inv ;
  M_inv[0][1] = cofM[1][0]*det_inv ;
  M_inv[1][1] = cofM[1][1]*det_inv ;
  M_inv[2][1] = cofM[1][2]*det_inv ;
  M_inv[0][2] = cofM[2][0]*det_inv ;
  M_inv[1][2] = cofM[2][1]*det_inv ;
  M_inv[2][2] = cofM[2][2]*det_inv ;
  return det;
}








void ContactTDEnforcement::Display_Enforcement_Data()
{
#if CONTACT_DEBUG_PRINT_LEVEL>=2
  ContactParOStream& postream = ParOStream();

  const char* model[15] = { "TD_FRICTIONLESS",    
                      "TD_CONSTANT_FRICTION", 
                      "TD_TIED",              
                      "TD_SPOT_WELD",         
                      "TD_PRESSURE_DEPENDENT",
                      "TD_VELOCITY_DEPENDENT",
                      "TD_SPRING_WELD",       
                      "TD_ADHESION",          
                      "TD_COHESIVE_ZONE",     
                      "TD_JUNCTION",          
                      "TD_THREADED",          
                      "TD_PV_DEPENDENT",      
                      "TD_AREA_WELD",         
                      "TD_USER",              
                      "TD_SHARED" };
  postream << "TD Enforcement Models:\n";
  for( int i=0 ; i<num_enforcement_models ; ++i ){
    int id    = enforcement_models[i]->ID();
    int index = enforcement_models[i]->Type();
    postream << "  "<<i<<": Friction_Model[ID="<<id<<"] = "<<model[index-1]<<"\n";
  }
  postream << "TD Enforcement Data:\n";
  postream << "  Number of Entity Keys = " << number_entity_keys << "\n";
  for( int ifaces=0 ; ifaces<number_entity_keys ; ++ifaces ){
    for( int inodes=0 ; inodes<number_entity_keys ; ++inodes ){
      int ID = (int) Enforcement_Data(FRICTION_MODEL_ID,ifaces,inodes);
      ContactEnfModel* enf_model = NULL;
      if (Initialized()) {
        enf_model = enforcement_models[ID];
      } else {
        for( int i=0 ; i<num_enforcement_models ; ++i ){
          if( ID == enforcement_models[i]->ID() ){
            enf_model = enforcement_models[i];
            break;
          }
        }
      }
      postream << "    Data for Nodes of Entity " << inodes+1
               << " against Faces of Entity " << ifaces+1 << "\n";
      postream << "        Kinematic Partition  = "
               << Enforcement_Data( KINEMATIC_PARTITION,
                               ifaces, inodes ) << "\n";
      if (enf_model!=NULL) {
      postream << "        Friction Model ID    = "
               << enf_model->ID() << ", "<<model[(int)(enf_model->Type())-1]<<"\n";
      } else {
      postream << "        Friction Model ID    = "
               << "??, Unknown\n";
      }
    }
  }
#endif
}

void ContactTDEnforcement::Display0_Enforcement_Data()
{
  const char* model[15] = { "TD_FRICTIONLESS",    
                      "TD_CONSTANT_FRICTION", 
                      "TD_TIED",              
                      "TD_SPOT_WELD",         
                      "TD_PRESSURE_DEPENDENT",
                      "TD_VELOCITY_DEPENDENT",
                      "TD_SPRING_WELD",       
                      "TD_ADHESION",          
                      "TD_COHESIVE_ZONE",     
                      "TD_JUNCTION",          
                      "TD_THREADED",          
                      "TD_PV_DEPENDENT",      
                      "TD_AREA_WELD",         
                      "TD_USER",              
                      "TD_SHARED" };
  std::cout << "TD Enforcement Models:\n";
  for( int i=0 ; i<num_enforcement_models ; ++i ){
    int id    = enforcement_models[i]->ID();
    int index = enforcement_models[i]->Type();
    std::cout << "  "<<i<<": Friction_Model[ID="<<id<<"] = "<<model[index-1]<<"\n";
  }
  std::cout << "TD Enforcement Data:\n";
  std::cout << "  Number of Entity Keys = " << number_entity_keys << "\n";
  for( int ifaces=0 ; ifaces<number_entity_keys ; ++ifaces ){
    for( int inodes=0 ; inodes<number_entity_keys ; ++inodes ){
      int ID = (int) Enforcement_Data(FRICTION_MODEL_ID,ifaces,inodes);
      ContactEnfModel* enf_model = NULL;
      if (Initialized()) {
        enf_model = enforcement_models[ID];
      } else {
        for( int i=0 ; i<num_enforcement_models ; ++i ){
          if( ID == enforcement_models[i]->ID() ){
            enf_model = enforcement_models[i];
            break;
          }
        }
      }
      std::cout << "    Data for Nodes of Entity " << inodes+1
                << " against Faces of Entity " << ifaces+1 << "\n";
      std::cout << "        Kinematic Partition  = "
                << Enforcement_Data(KINEMATIC_PARTITION,ifaces,inodes) << "\n";
      if (enf_model!=NULL) {
      std::cout << "        Friction Model ID    = "
                << enf_model->ID() << ", "<<model[(int)(enf_model->Type())-1]<<"\n";
      } else {
      std::cout << "        Friction Model ID    = "
                << "??, Unknown\n";
      }
    }
  }
}

void ContactTDEnforcement::TD_Enforcement_Clean_Up()
{
  ASSEMBLED_MASS.         Clear_Scratch();
  NODAL_MASS.             Clear_Scratch();
  NODAL_DENSITY.          Clear_Scratch();
  NODAL_WAVESPEED.        Clear_Scratch();
  ASSEMBLED_FORCE.        Clear_Scratch();
  SN_ACCELERATION.        Clear_Scratch();
  INC_FORCE.              Clear_Scratch();
  A_PREDICTION.           Clear_Scratch();
  A_CORRECTION.           Clear_Scratch();
  A_CONSTITUTIVE.         Clear_Scratch();  
  TOTAL_FORCE.            Clear_Scratch();  
  NEW_POSITION.           Clear_Scratch();  

  PREDICTED_POS.          Clear_Scratch();  
  CURRENT_POS.            Clear_Scratch();  
  NUM_KIN_CONSTR.         Clear_Scratch();
  KIN_CONSTR_VECTOR.      Clear_Scratch();

  NEI_TOTAL_SLIP.         Clear_Scratch();
  NEI_GAP_TO_ENFORCE.     Clear_Scratch();
  NEI_REL_DISP.           Clear_Scratch();
  NEI_FORCE.              Clear_Scratch();
  NEI_TAN_SLIP.           Clear_Scratch();
  NEI_TARGET_GAP.         Clear_Scratch();
  NEI_KINEMATIC_PARTITION.Clear_Scratch();

  node_entity_list.Clear();
}


Real 
ContactTDEnforcement::Auto_Kinematic_Partition( ContactNodeFaceInteraction *cnfi) {
  Real kin_partition = 0.5;

  // Get the Node and Face for this interaction
  ContactNode<Real>* node = cnfi->Node();
  ContactFace<Real>* face = cnfi->Face();

  int node_index = node->EnfArrayIndex();
  Real density_sn   = *(NODAL_DENSITY.Get_Scratch(node_index));
  Real wavespeed_sn = *(NODAL_WAVESPEED.Get_Scratch(node_index));
  
  // Calculate the face distribution factors
  Real shape_functions[MAX_NODES_PER_FACE];
  PRECONDITION( MAX_NODES_PER_FACE>=face->Nodes_Per_Face() );
  Real* coordinates = cnfi->Vector_Var(ContactNodeFaceInteraction::COORDINATES);
  face->Evaluate_Shape_Functions( coordinates, shape_functions );      
  
  Real face_density   = 0.;
  Real face_wavespeed = 0.;
  for(int j=0; j<face->Nodes_Per_Face() ; ++j ){
    face_density   += shape_functions[j] * (*(NODAL_DENSITY.  Get_Scratch(face->Node(j))));
    face_wavespeed += shape_functions[j] * (*(NODAL_WAVESPEED.Get_Scratch(face->Node(j))));
  }
  kin_partition = (face_density*face_wavespeed)/
    (density_sn*wavespeed_sn + face_density*face_wavespeed);
#ifdef CONTACT_DEBUG_NODE
      ContactParOStream& postream = ParOStream();
      if( topology->Is_a_Debug_Node(cnfi->Node()) ){
	postream << "Information for Debug Node: " 
		 << cnfi->Node()->Exodus_ID() << "\n"  
		 << "  face density   = " << face_density << "\n"
		 << "  face wavespeed = " << face_wavespeed << "\n"
		 << "  node density   = " << density_sn << "\n"
		 << "  node wavespeed = " << wavespeed_sn << "\n"
		 << "  Automatic Kinematic Partition =  " 
		 << kin_partition << "\n";
      }
#endif
  return ( kin_partition );
}


void ContactTDEnforcement::Compute_Kinematic_Quantities( int iteration )
{

#ifdef CONTACT_DEBUG_NODE
  ContactParOStream& postream = ParOStream();
#endif

  // (1) Compute the target gap. 
  // This is based on assuming we want all of g_dot*dt and a fraction
  // of gap_old based on an energy consideration relative to g_dot*dt.
  // We are going to only consider the "normal" direction.  Therefore
  // we will project gap_cur & old onto the normal direction to compute
  // the target normal_gap
  //
  // NOTE at the first iteration NEW = PREDICTED != CURRENT
  // so relative displacment != current gap
  if( iteration == 0 ){
    //
    // Get the variable handles for everything I might need
    //
    for(Node_Constraint_Group cnei_group = node_entity_list.Node_Group_Start(); 
        cnei_group.Valid_Group();
        cnei_group = node_entity_list.Node_Group_Next()) {

      int node_index = cnei_group.Get_Constrained_Node_Index();

      for(int i = 0; i<cnei_group.Num_Interactions(); ++i) {

        int index = cnei_group.Get_Index(i);
        ContactNodeEntityInteraction *cnei = cnei_group.Get_Interaction(i);

        Real& gap_cur=cnei->Get_Gap_Cur();
        Real& gap_old=cnei->Get_Gap_Old();
        Real gap_old_to_enforce = 0.0;
        
        const int model_type = TDEnfModel(cnei)->Interaction_Type(cnei);

        if ( model_type == ContactTDEnfModel::TDEM_ADHESIVE
             ||  model_type == ContactTDEnfModel::TDEM_ADHESIVE_FRICTION
             ||  model_type == ContactTDEnfModel::TDEM_SPRINGY ) {
          if (TDEnfModel(cnei)->Active_Interaction(cnei,gap_cur)) {
	    gap_old_to_enforce = gap_old;
          }
        } else {

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
  
	  Real gtotal = gap_cur;

	  // If the the total_gap is non-negative put everything in gap_cur
	  if( gtotal >= 0.0 ){
	    // Cases 1, 2 & 5
	    gap_old = 0.0;
	  } else {
	    Real* normal = cnei->Get_Pushback_Dir();
	    //
	    // Compute the amount of node motion along the face normal
	    //
	    Real* cur0 = CURRENT_POS.Get_Scratch(node_index);
	    Real* curp = PREDICTED_POS.Get_Scratch(node_index);
	    Real node_gap_cur = ((curp[0]-cur0[0])*normal[0] +
				 (curp[1]-cur0[1])*normal[1] +
				 (curp[2]-cur0[2])*normal[2] );

	    Real physical_gap;
	    //
	    // Now compute the amount of face motion along the face normal
	    //
            Real contact_point_motion[3] = {0.0, 0.0, 0.0};
            const int num_face_nodes = cnei_group.Get_Num_Face_Nodes(i);
	    for( int j=0 ; j < num_face_nodes ; ++j ){
              const int face_node_index = cnei_group.Get_Face_Node_Index(i,j);
              Real face_shape_function = cnei_group.Get_Face_Shape_Function(i, j);
	      Real* c0 = CURRENT_POS.Get_Scratch(face_node_index);
	      Real* cp = PREDICTED_POS.Get_Scratch(face_node_index);
	      contact_point_motion[0] += face_shape_function*(cp[0]-c0[0]);
	      contact_point_motion[1] += face_shape_function*(cp[1]-c0[1]);
	      contact_point_motion[2] += face_shape_function*(cp[2]-c0[2]);
	    }
	    //
	    //  NKC, once again spurious skip of this calculation for analytic surfaces
	    //
            if(num_face_nodes > 0) {
	      Real face_gap_cur = (contact_point_motion[0]*normal[0] +
				   contact_point_motion[1]*normal[1] +
				   contact_point_motion[2]*normal[2] );
  
	      // Take the amount of the gap that was obtainable this step and store
	      // in gap_cur and leave the rest in gap_old;
	      physical_gap = node_gap_cur - face_gap_cur;
	    } else {
	      physical_gap = node_gap_cur;
	    }

	    if( physical_gap >= 0.0 ){
	      // case 6 (the other cases of moving "out" have been handled)
	      gap_old = gtotal;
	      gap_cur = 0.0;
	    } else {
	      if( gtotal < physical_gap ){
		// case 4
		gap_cur = physical_gap;
		gap_old = gtotal - gap_cur;
	      } else {
		// case 3
		gap_old = 0.0;
	      }
	    }
	  }
	  if( gap_cur < 0.0  && gap_old < 0.0){
	    Real limit_old_gap_to_enforce = std::max( FRAC_CUR_GAP*gap_cur, gap_old );
	    if( gap_old < limit_old_gap_to_enforce ) {
	      gap_old_to_enforce = limit_old_gap_to_enforce;
	    } else {
	      gap_old_to_enforce = gap_old;
 	    }
	  }
	}

        Real* nor_dir = cnei->Get_Normal();
        Real* pb_dir  = cnei->Get_Pushback_Dir();
        Real  dot = 0.0;
        for( int j=0 ; j<NDIM ; j++ ){ dot += pb_dir[j]*nor_dir[j];}
        *(NEI_TARGET_GAP.Get_Scratch(index))     = dot*(gap_old-gap_old_to_enforce);
        *(NEI_GAP_TO_ENFORCE.Get_Scratch(index)) = dot*(gap_cur+gap_old_to_enforce);

#ifdef CONTACT_DEBUG_NODE
        if( topology->Is_a_Debug_Node(cnei->Node()) ){
          postream << "Information for Debug Node: " 
      	           << cnei->Node()->Exodus_ID()  
                   << " for iteration " << iteration << "\n";
          postream << "   Gap to enforce     = " 
		   <<  *(NEI_GAP_TO_ENFORCE.Get_Scratch(index)) << "\n";
          postream << "   Target gap         = " 
		   <<  *(NEI_TARGET_GAP.Get_Scratch(index)) << "\n";
          postream << "   Cur gap         = " 
		   <<  cnei->Get_Gap_Cur() << "\n";
          postream << "   Old gap         = " 
		   <<  cnei->Get_Gap_Old() << "\n";
        }
#endif
      }
    }
  } else {
    for(Node_Constraint_Group cnei_group = node_entity_list.Node_Group_Start(); 
        cnei_group.Valid_Group();
        cnei_group = node_entity_list.Node_Group_Next()) {

      const int node_index = cnei_group.Get_Constrained_Node_Index();
      Real *x_sn = NEW_POSITION.Get_Scratch(node_index);

      for(int i = 0; i<cnei_group.Num_Interactions(); ++i) {

        int interaction_index = cnei_group.Get_Index(i);
        ContactNodeEntityInteraction *cnei = cnei_group.Get_Interaction(i);

        Real target_gap = *(NEI_TARGET_GAP.Get_Scratch(interaction_index));
        
        const int model_type = TDEnfModel(cnei)->Interaction_Type(cnei);

        Real  remaining_gap = 0.0;  // NOTE this is different than the NFI case
        Real  x_cp[3] = {0.,0.,0.}; // Here for debug node
	Real* normal = cnei->Get_Normal();

	//
	// Compute the total gap in the normal direction
	//
        const int num_face_nodes = cnei_group.Get_Num_Face_Nodes(i);
  	for(int k=0 ; k< num_face_nodes ; ++k ){
          const int face_node_index = cnei_group.Get_Face_Node_Index(i, k);
          const Real face_shape_function = cnei_group.Get_Face_Shape_Function(i, k);
	  Real* x_mn_K = NEW_POSITION.Get_Scratch(face_node_index);
	  for(int j=0 ; j<NDIM ; ++j ) x_cp[j] += face_shape_function*x_mn_K[j];
	}

	//
	//  NKC, This will only compute remaining gap for node face interactions.  For node surface leave at 0.0
	//  NKC, Note:  This is almost certainly wrong for node-surface interactions!!!
	//  Need to check and debug to get Node-surface interactions working!!!!
	//
        if(num_face_nodes > 0) {
	  for (int j=0 ; j<NDIM ; ++j ) remaining_gap += (x_sn[j]-x_cp[j])*normal[j];
        }
        if ( model_type == ContactTDEnfModel::TDEM_ADHESIVE
	     || model_type == ContactTDEnfModel::TDEM_ADHESIVE_FRICTION
             || model_type == ContactTDEnfModel::TDEM_SPRINGY ) {
          if (TDEnfModel(cnei)->Active_Interaction(cnei,remaining_gap)) {
	    *(NEI_GAP_TO_ENFORCE.Get_Scratch(interaction_index)) = remaining_gap;
	  } else {
            *(NEI_GAP_TO_ENFORCE.Get_Scratch(interaction_index)) = 0.0;
          }
        } else {
          if( remaining_gap < target_gap ){
            *(NEI_GAP_TO_ENFORCE.Get_Scratch(interaction_index)) = remaining_gap - target_gap;
          } else {
            *(NEI_GAP_TO_ENFORCE.Get_Scratch(interaction_index)) = 0.0;
          }
        }

#ifdef CONTACT_DEBUG_NODE
        if( topology->Is_a_Debug_Node(cnei->Node()) ){
          postream << "Information for Debug Node: " 
      	           << cnei->Node()->Exodus_ID()  
                   << " for iteration " << iteration << "\n";
          postream << "   Node Predicted Coords = " << x_sn[0]
          	   << " " << x_sn[1] << " " << x_sn[2] << "\n";
          postream << "   Contact Point         = " << x_cp[0]
                   << " " << x_cp[1] << " " << x_cp[2]
                   << "\n";
          postream << "   Remaining Nor Gap     = " << remaining_gap << "\n";
          postream << "   Gap to enforce     = " 
		   <<  *(NEI_GAP_TO_ENFORCE.Get_Scratch(interaction_index)) << "\n";
          postream << "   Target gap         = " 
		   <<  *(NEI_TARGET_GAP.Get_Scratch(interaction_index)) << "\n";
        }
#endif
      }
    }
  }
  //
  // (2) Compute slip, relative position, tangential slip
  //
  for(Node_Constraint_Group cnei_group = node_entity_list.Node_Group_Start(); 
      cnei_group.Valid_Group();
      cnei_group = node_entity_list.Node_Group_Next()) {

    const int node_index = cnei_group.Get_Constrained_Node_Index();
    Real *x_sn_new = NEW_POSITION.Get_Scratch(node_index);

    for(int i = 0; i<cnei_group.Num_Interactions(); ++i) {

      int interaction_index = cnei_group.Get_Index(i);
      ContactNodeEntityInteraction *cnei = cnei_group.Get_Interaction(i);

      Real* x_sn_cur = CURRENT_POS.Get_Scratch(node_index);
      Real* rel_disp = NEI_REL_DISP.Get_Scratch(interaction_index);
      Real* slip     = NEI_TOTAL_SLIP.Get_Scratch(interaction_index);
      Real* slip_t   = NEI_TAN_SLIP.Get_Scratch(interaction_index);
      Real* nm       = cnei->Get_Normal();


      Real xp_new[3] = {0.0,0.0,0.0};
      const int num_face_nodes = cnei_group.Get_Num_Face_Nodes(i);
      if(num_face_nodes > 0) {
        Real xp_cur[3] = {0.0,0.0,0.0};
        for(int ii=0 ; ii<num_face_nodes ; ++ii ){
          const int node_face_index = cnei_group.Get_Face_Node_Index(i, ii);
          const Real face_shape_function = cnei_group.Get_Face_Shape_Function(i ,ii);
	  Real* x_mn_cur = CURRENT_POS.Get_Scratch(node_face_index);
	  Real* x_mn_new = NEW_POSITION.Get_Scratch(node_face_index);
	  for(int j=0 ; j<NDIM ; ++j ){
	    xp_cur[j] +=  face_shape_function*x_mn_cur[j] ;
  	    xp_new[j] +=  face_shape_function*x_mn_new[j] ;
	  }
        }
        // subtract initial gap in the direction of the current normal
        ContactNodeFaceInteraction * cnfi = dynamic_cast<ContactNodeFaceInteraction*>(cnei);
        ContactEnfModel *model      = TDEnfModel(cnei);
        ContactTDTied   *tied_model = dynamic_cast<ContactTDTied*>(model);
        if (cnfi != NULL && tied_model != NULL) {
          Real gap0 = 0.0;
          if(new_tied_enforcement == ContactSearch::ACTIVE) gap0 = cnei->Get_Gap_Initial();
          Real n_dir[3] = {0,0,0}; 
          Real* x_face[MAX_NODES_PER_FACE];
          PRECONDITION( MAX_NODES_PER_FACE>= num_face_nodes);
          // current (read "old") gap vector
          for(int ii=0 ; ii<num_face_nodes ; ++ii ){
            const int node_face_index = cnei_group.Get_Face_Node_Index(i, ii);
	    x_face[ii] = CURRENT_POS.Get_Scratch(node_face_index);
          }
          cnfi->Face()->Compute_Normal(x_face,cnfi->Get_Coordinates(),n_dir);
          for(int j=0 ; j<NDIM ; ++j ){ xp_cur[j] +=  gap0*n_dir[j]; }
          // new  gap vector
          for(int ii=0 ; ii<num_face_nodes ; ++ii ){
            const int node_face_index = cnei_group.Get_Face_Node_Index(i, ii);
	    x_face[ii] = NEW_POSITION.Get_Scratch(node_face_index);
          }
          cnfi->Face()->Compute_Normal(x_face,cnfi->Get_Coordinates(),n_dir);
          for(int j=0 ; j<NDIM ; ++j ){ xp_new[j] +=  gap0*n_dir[j]; }
        }

        for(int j = 0; j < NDIM; ++j) {
	  slip[j] = (x_sn_new[j] - xp_new[j]) - (x_sn_cur[j] - xp_cur[j]);
        }
      } else {
        ContactNodeSurfaceInteraction* cnsi = dynamic_cast<ContactNodeSurfaceInteraction*>(cnei);
        if(cnsi != NULL){
          ContactAnalyticSurface* surface = cnsi->Surface();
          Real gap;
          Real time_to_contact;
          Real pushback_dir[3];
          int  location;
          if ( ! surface->Process( x_sn_new, gap, xp_new, nm, pushback_dir, time_to_contact, location ) ) continue;

          for(int j = 0; j < NDIM; ++j) {
            slip[j] = x_sn_new[j] - x_sn_cur[j];
	  }
        }
      }

      Real dot = slip[0]*nm[0] + slip[1]*nm[1] + slip[2]*nm[2];
      for(int j=0 ; j<NDIM ; ++j ) { 
        rel_disp[j] = x_sn_new[j] - xp_new[j];
        slip_t[j] = slip[j] - dot*nm[j]; 
      }
#ifdef CONTACT_DEBUG_NODE
      if( topology->Is_a_Debug_Node(cnei->Node()) ){ 
        postream << "Information for Debug Node: " 
                 << cnei->Node()->Exodus_ID() << "\n"; 
        postream << "   Relative Displacement = " << rel_disp[0]
         	 << " " << rel_disp[1] << " " << rel_disp[2] << "\n";
        postream << "   Total Slip           = " << slip[0]
                 << " " << slip[1] << " " << slip[2] << "\n";
        postream << "   Tangential Slip      = " << slip_t[0]
                 << " " << slip_t[1] << " " << slip_t[2] << "\n";
        dot = rel_disp[0]*nm[0] + rel_disp[1]*nm[1] + rel_disp[2]*nm[2];//HACK
        postream <<   "   Initial gap = " << cnei->Get_Gap_Initial()
		 << "\n   Old gap     = " << cnei->Get_Gap_Old()
		 << "\n   Current gap = " << cnei->Get_Gap_Cur()
		 << "\n   New gap     = " << dot << "\n";
      }
#endif
    }
  }
#ifdef CONTACT_DEBUG_NODE
  postream.flush();
#endif
}


ContactSearch::ContactErrorCode
ContactTDEnforcement::Update_For_Topology_Change( int new_number_of_nodes,
						  int* old_to_new_map )
{
  if(CALC_PLOT_FORCE) {
    Real* old_plot_force = plot_force;
    if( new_number_of_nodes ) {
      plot_force = new Real[3*new_number_of_nodes];
    } else { 
      plot_force = NULL;
    }
    Map_Array_to_New_Topology( new_number_of_nodes,
			       number_of_nodes,
			       old_to_new_map,
			       old_plot_force,
			       plot_force,
			       3,
			       BY_VARIABLE );
  }
  Real* old_disp_old_aug = disp_old_aug;  
  if( new_number_of_nodes ) {
    disp_old_aug = new Real[3*new_number_of_nodes];
  } else { 
    disp_old_aug = NULL;
  }
  Map_Array_to_New_Topology( new_number_of_nodes,
			     number_of_nodes,
			     old_to_new_map,
			     old_disp_old_aug,
			     disp_old_aug,
			     3,
			     BY_VARIABLE );

  // Update the friction models
  for(int i=0 ; i<num_enforcement_models ; ++i )
    enforcement_models[i]->Update_For_Topology_Change( new_number_of_nodes,
						       old_to_new_map );

  number_of_nodes = new_number_of_nodes;
  return ContactSearch::NO_ERROR;
}

ContactTDEnfModel*
ContactTDEnforcement::TDEnfModel( ContactNodeEntityInteraction* cnei)
{
  int e_key = cnei->Get_Entity_Key();
  int n_key = cnei->Get_Node_Key();
  int index = (int) Enforcement_Data(FRICTION_MODEL_ID, e_key, n_key );
  return (ContactTDEnfModel*) enforcement_models[index];
}

Real ContactTDEnforcement::Enforcement_Data( Enforcement_Data_Index EDI, ContactNodeEntityInteraction* cnei)
{
  int e_key = cnei->Get_Entity_Key();
  int n_key = cnei->Get_Node_Key();
  return enforcement_data.Get_Data( EDI, e_key, n_key );
}

Real ContactTDEnforcement::Enforcement_Data( Enforcement_Data_Index EDI,
					     int key_1, int key_2 )
{
  return enforcement_data.Get_Data( EDI, key_1, key_2 );
}

void ContactTDEnforcement::Set_Save_Nodal_Vars(bool flag)
{
  if(SAVE_CVARS == false && flag == true)
  {
    SAVE_CVARS = flag;

    conface                   = new Real[number_of_nodes];

    normal_force_mag          = new Real[number_of_nodes];
    normal_traction_mag       = new Real[number_of_nodes];
    tangential_force_mag      = new Real[number_of_nodes];
    tangential_traction_mag   = new Real[number_of_nodes];

    normal_dir_x              = new Real[number_of_nodes];
    normal_dir_y              = new Real[number_of_nodes];
    normal_dir_z              = new Real[number_of_nodes];

    tangential_dir_x          = new Real[number_of_nodes];
    tangential_dir_y          = new Real[number_of_nodes];
    tangential_dir_z          = new Real[number_of_nodes];

    slipmag                   = new Real[number_of_nodes];
    nodal_dissipation         = new Real[number_of_nodes];
    nodal_dissipation_density = new Real[number_of_nodes];
    contact_area              = new Real[number_of_nodes];
    plot_gap_cur              = new Real[number_of_nodes];
    plot_gap_old              = new Real[number_of_nodes];
    kinematic_partition       = new Real[number_of_nodes];
    Set_CVARS();
  }
  else if(SAVE_CVARS == true && flag == false)
  {
    SAVE_CVARS = flag;

    delete conface;

    delete normal_force_mag;
    delete normal_traction_mag;
    delete tangential_force_mag;
    delete tangential_traction_mag;

    delete normal_dir_x;
    delete normal_dir_y;
    delete normal_dir_z;

    delete tangential_dir_x;
    delete tangential_dir_y;
    delete tangential_dir_z;

    delete slipmag;
    delete nodal_dissipation;
    delete nodal_dissipation_density;
    delete contact_area;
    delete plot_gap_cur;
    delete plot_gap_old;
    delete kinematic_partition;

    conface                   = NULL;

    normal_force_mag          = NULL;
    normal_traction_mag       = NULL;
    tangential_force_mag      = NULL;
    tangential_traction_mag   = NULL;

    normal_dir_x              = NULL;
    normal_dir_y              = NULL;
    normal_dir_z              = NULL;

    tangential_dir_x          = NULL;
    tangential_dir_y          = NULL;
    tangential_dir_z          = NULL;

    slipmag                   = NULL;
    nodal_dissipation         = NULL;
    nodal_dissipation_density = NULL;
    contact_area              = NULL;
    plot_gap_cur              = NULL;
    plot_gap_old              = NULL;
    kinematic_partition       = NULL;
  }
  
}

void ContactTDEnforcement::Initialize_CVARS()
{
  PRECONDITION( SAVE_CVARS );

  // Initialize CVARS
  std::memset( normal_force_mag,          0, number_of_nodes*sizeof(Real) );
  std::memset( normal_traction_mag,       0, number_of_nodes*sizeof(Real) );
  std::memset( tangential_force_mag,      0, number_of_nodes*sizeof(Real) );
  std::memset( tangential_traction_mag,   0, number_of_nodes*sizeof(Real) );

  std::memset( normal_dir_x,              0, number_of_nodes*sizeof(Real) );
  std::memset( normal_dir_y,              0, number_of_nodes*sizeof(Real) );
  std::memset( normal_dir_z,              0, number_of_nodes*sizeof(Real) );

  std::memset( tangential_dir_x,          0, number_of_nodes*sizeof(Real) );
  std::memset( tangential_dir_y,          0, number_of_nodes*sizeof(Real) );
  std::memset( tangential_dir_z,          0, number_of_nodes*sizeof(Real) );

  std::memset( slipmag,                   0, number_of_nodes*sizeof(Real) );
  std::memset( nodal_dissipation,         0, number_of_nodes*sizeof(Real) );
  std::memset( nodal_dissipation_density, 0, number_of_nodes*sizeof(Real) );
  std::memset( contact_area,              0, number_of_nodes*sizeof(Real) );
  std::memset( plot_gap_cur,              0, number_of_nodes*sizeof(Real) );
  std::memset( plot_gap_old,              0, number_of_nodes*sizeof(Real) );
  std::memset( kinematic_partition,       0, number_of_nodes*sizeof(Real) );

  for(int i=0 ; i<number_of_nodes ; ++i ) conface[i] = 0.5;
}

void ContactTDEnforcement::Set_CVARS()
{
  PRECONDITION( SAVE_CVARS );

  // Initialize CVARS
  std::memset( normal_force_mag,          0, number_of_nodes*sizeof(Real) );
  std::memset( normal_traction_mag,       0, number_of_nodes*sizeof(Real) );
  std::memset( tangential_force_mag,      0, number_of_nodes*sizeof(Real) );
  std::memset( tangential_traction_mag,   0, number_of_nodes*sizeof(Real) );

  std::memset( normal_dir_x,              0, number_of_nodes*sizeof(Real) );
  std::memset( normal_dir_y,              0, number_of_nodes*sizeof(Real) );
  std::memset( normal_dir_z,              0, number_of_nodes*sizeof(Real) );

  std::memset( tangential_dir_x,          0, number_of_nodes*sizeof(Real) );
  std::memset( tangential_dir_y,          0, number_of_nodes*sizeof(Real) );
  std::memset( tangential_dir_z,          0, number_of_nodes*sizeof(Real) );

  std::memset( slipmag,                   0, number_of_nodes*sizeof(Real) );
  std::memset( nodal_dissipation,         0, number_of_nodes*sizeof(Real) );
  std::memset( nodal_dissipation_density, 0, number_of_nodes*sizeof(Real) );
  std::memset( contact_area,              0, number_of_nodes*sizeof(Real) );
  std::memset( plot_gap_cur,              0, number_of_nodes*sizeof(Real) );
  std::memset( plot_gap_old,              0, number_of_nodes*sizeof(Real) );
  std::memset( kinematic_partition,       0, number_of_nodes*sizeof(Real) );

  for(int i=0 ; i<number_of_nodes ; ++i ) conface[i] = 0.5;

  Compute_Kinematic_Quantities( 1 );

  for(Node_Constraint_Group cnei_group = node_entity_list.Node_Group_Start();
      cnei_group.Valid_Group();
      cnei_group = node_entity_list.Node_Group_Next()) {
    int ncc = cnei_group.Num_Interactions();
    ContactNodeEntityInteraction* cnei = cnei_group.Get_Interaction(0);
    int cnei_index = cnei_group.Get_Index(0);
    const int node_enf_index = cnei_group.Get_Constrained_Node_Index();
    int node_host_index = Get_Node_Host_Array_Index(cnei->Node());

    Real* slip = NEI_TOTAL_SLIP.Get_Scratch(cnei_index);
    Real* f_t  = TOTAL_FORCE.Get_Scratch(node_enf_index);
    for (int k=0; k< NDIM; ++k) { nodal_dissipation[node_host_index] += f_t[k]*slip[k];}
    if( conface[node_host_index] == 0.5 )
      conface[node_host_index] = 1.0;
    else
      conface[node_host_index] += 1.0;

    Real* f     = NEI_FORCE.   Get_Scratch(cnei_index);
    Real* t_dir = NEI_TAN_SLIP.Get_Scratch(cnei_index);

    contact_area[node_host_index] = cnei->Get_Node_Area();
    plot_gap_cur[node_host_index] = cnei->Get_Gap_Cur();
    plot_gap_old[node_host_index] = cnei->Get_Gap_Old();
    kinematic_partition[node_host_index] = 
      *(NEI_KINEMATIC_PARTITION.Get_Scratch(cnei_index));

    Real* n_dir = cnei->Get_Normal();
    normal_dir_x[node_host_index] = n_dir[0];
    normal_dir_y[node_host_index] = n_dir[1];
    normal_dir_z[node_host_index] = n_dir[2];

    tangential_dir_x[node_host_index] = t_dir[0];
    tangential_dir_y[node_host_index] = t_dir[1];
    tangential_dir_z[node_host_index] = t_dir[2];

    Real f_tan[3];
    Real f_n = 0.0;
    Real mag_f_tan = 0.0;
    for (int k=0;k<dimensionality;++k){f_n += n_dir[k]*f[k] ;}
    for (int k=0;k<dimensionality;++k){
      f_tan[k] = f[k] - f_n * n_dir[k];
      mag_f_tan += f_tan[k]*f_tan[k] ;
    }
    mag_f_tan = std::sqrt(mag_f_tan);
    normal_force_mag[node_host_index] = f_n;
    tangential_force_mag[node_host_index] = mag_f_tan;

    //calculate traction
    VariableHandle NODE_NORMAL = topology->Variable_Handle(ContactTopology::Node_Normal);
    Real *node_normal = cnei->Node()->Variable(NODE_NORMAL);
    Real *entity_normal = cnei->Get_Normal();
    Real node_entity_dot_product = 0.0;
    for(int k = 0; k < dimensionality; ++k){
      node_entity_dot_product += node_normal[k] * entity_normal[k];
    }

//     normal_traction_mag[node_host_index]     = contact_area[node_host_index] * node_entity_dot_product * f_n;
//     tangential_traction_mag[node_host_index] = contact_area[node_host_index] * node_entity_dot_product * mag_f_tan;
    normal_traction_mag[node_host_index]     =  (node_entity_dot_product * f_n)       / contact_area[node_host_index];
    tangential_traction_mag[node_host_index] =  (node_entity_dot_product * mag_f_tan) / contact_area[node_host_index];

    Real tmp = std::sqrt( t_dir[0]*t_dir[0]+t_dir[1]*t_dir[1]+ t_dir[2]*t_dir[2]);
    slipmag[node_host_index] = tmp;
    if(ncc > 1) {
      for(int i=0; i < ncc ; ++i) {
        cnei_index = cnei_group.Get_Index(i);
        f_t  = TOTAL_FORCE.Get_Scratch(node_enf_index);
        slip = NEI_TOTAL_SLIP.Get_Scratch(cnei_index);
        for (int k=0; k< NDIM; ++k) { 
          nodal_dissipation[node_host_index] += f_t[k]*slip[k];
        }
      }
    }
    nodal_dissipation_density[node_host_index] =
      nodal_dissipation[node_host_index] / contact_area[node_host_index];
  }
}

void ContactTDEnforcement::Set_CGVARS()
{
  PRECONDITION( SAVE_CGVARS );
  Real force_x     = 0.0;
  Real force_y     = 0.0;
  Real force_z     = 0.0;
  Real dissipation = 0.0;
  Real max_pen     = 0.0;
  
  Real* x_f  = NEW_POSITION.Get_Data();
  Real* f_t  = TOTAL_FORCE.Get_Data();
  Real* a_sn = SN_ACCELERATION.Get_Data();
  Real* m    = NODAL_MASS.Get_Data();
  Real* x_i  = CURRENT_POS.Get_Data();
  for(int i=0 ; i<number_of_nodes ; ++i ){
    force_x += f_t[0];
    force_y += f_t[1];
    force_z += f_t[2];
    dissipation += ( f_t[0]*( x_f[0] - x_i[0] )
                   + f_t[1]*( x_f[1] - x_i[1] )
                   + f_t[2]*( x_f[2] - x_i[2] )  )/dt ;
    Real mag = 0.0;
    for(int j=0 ; j<NDIM ; ++j ){ mag += m[0]*a_sn[j]*m[0]*a_sn[j]; }
    max_pen = std::max(max_pen,std::sqrt(mag));
    x_f += 3;
    f_t += 3;
    a_sn += 3;
    m += 1;
    x_i += 3;
  }
  Real sum_input[4], sum_output[4];
  sum_input[0] = force_x;
  sum_input[1] = force_y;
  sum_input[2] = force_z;
  sum_input[3] = dissipation;
  contact_global_sum(sum_input, sum_output, 4, communicator);
  global_force_x = sum_output[0];
  global_force_y = sum_output[1];
  global_force_z = sum_output[2];
  global_dissipation = sum_output[3];


  // l2 norm of sum of the contact force
  global_force_norm =  std::sqrt(global_force_x*global_force_x
                           +global_force_y*global_force_y
                           +global_force_z*global_force_z);
  // max norm of the constraint/penalty force  
  Real global_max_pen = contact_global_maximum(max_pen, communicator);
  global_constraint_norm = global_max_pen;
}


void ContactTDEnforcement::Store_Shell_Final_Lofted_Positions()
{
#ifdef CONTACT_DEBUG_NODE
  ContactParOStream& postream = ParOStream();
#endif
  // This routines computes and stores the final "compressed" configuration
  // for the shells.  The previous lofting is set for the shell nodes so the
  // CURRENT configuration can be set to the final compressed location.  This
  // will allow the partitioning between g_dot*dt and g_old.

  Real dt_sq = (0.5*(dt_old+dt))*dt;

  // Now compute the final lofting.
  Real* xp = PREDICTED_POS.Get_Data();
  Real* xf = NEW_POSITION. Get_Data();
  Real*  m = NODAL_MASS.   Get_Data();
  Real*  f = TOTAL_FORCE.  Get_Data();


  for(int i=0 ; i<number_of_nodes ; ++i, xp += 3, xf += 3, m += 1, f+= 3 ){
    ContactNode<Real>* node = enforcement_node_list[i];
    if( !node->Is_a_Shell_Node() ) continue;
    if( node->Physical_Type() == ContactNode<Real>::SHELL_NODE ){
      ContactShellNode* snode = static_cast<ContactShellNode*> (node);
      Real* loft = snode->Previous_Lofting();
#ifdef CONTACT_DEBUG_NODE
      bool PRINT_THIS_NODE = topology->Is_a_Debug_Node( node );
      if( PRINT_THIS_NODE ){
	postream << "Storing Final Lofting for Node " << node->Exodus_ID()
		 << "\n    m = " << m[0]
		 << "\n   xf = " << xf[0] << " " << xf[1] << " " << xf[2]
		 << "\n   xp = " << xp[0] << " " << xp[1] << " " << xp[2]
		 << "\n    f = " << f[0] << " " << f[1] << " " << f[2] 
		 << "\n   lp = " << loft[0] << " " << loft[1] << " " 
		 << loft[2] << "\n";
      }
#endif
      Real m_inv_dt = 1.0/(m[0]) * dt_sq;

      loft[0] = xf[0] - xp[0] + loft[0] - f[0]*m_inv_dt;
      loft[1] = xf[1] - xp[1] + loft[1] - f[1]*m_inv_dt;
      loft[2] = xf[2] - xp[2] + loft[2] - f[2]*m_inv_dt;
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE ){
	postream<< "   lf = " << loft[0] << " " << loft[1] << " " 
		 << loft[2] << "\n";
      }
#endif
    }
  }
#ifdef CONTACT_DEBUG_NODE
    postream.flush();
#endif
}

void ContactTDEnforcement::Assemble_Shell_Forces()
// KHB: 3/17/03
// We've decided that we can not treat the two sides of the shells
// independently (tied cases in tension were unconditionally 
// unstable).  Our solution is to let the two sides of the shells
// know about one another during the iteration process.  To do that,
// I'm going to "assemble" the force from all the nodes that were
// lofted from the same shell node to each of the lofted nodes.
{

#ifdef CONTACT_DEBUG_NODE
  ContactParOStream& postream = ParOStream();
#endif

  ContactShellHandler* sh = topology->Shell_Handler();
  int  num_host_code_nodes = sh->Number_Host_Code_Nodes();
  
  for(int i=0 ; i<num_host_code_nodes ; ++i ){
    if( sh->Num_Acme_Nodes_for_Host_Node(i) > 1 ){
      Real assembled_force[3];
      assembled_force[0] = 0.0;
      assembled_force[1] = 0.0;
      assembled_force[2] = 0.0;
#ifdef CONTACT_DEBUG_NODE
      int hi_cdn, lo_cdn;
      sh->Acme_NodeGID_for_Host_Node(i,0,hi_cdn,lo_cdn);
      ContactHostGlobalID gid_cdn(hi_cdn,lo_cdn);
      ContactNode<Real>* node_cdn = static_cast<ContactNode<Real>*>
         (topology->NodeList()->Find(gid_cdn));
      POSTCONDITION(node_cdn);
      bool PRINT_THIS_NODE = topology->Is_a_Debug_Node(node_cdn);
      if( PRINT_THIS_NODE ){
	postream << "Assembling Shell Forces for Node "
		 << node_cdn->Exodus_ID() << " (" 
		 << sh->Num_Acme_Nodes_for_Host_Node(i) 
		 << " lofted nodes)\n";
      }
#endif
      for(int j=0 ; j<sh->Num_Acme_Nodes_for_Host_Node(i) ; ++j ){
        int hi, lo;
        sh->Acme_NodeGID_for_Host_Node(i,j,hi,lo);
        ContactHostGlobalID gid(hi,lo);
        ContactNode<Real>* node = static_cast<ContactNode<Real>*>
          (topology->NodeList()->Find(gid));
        POSTCONDITION(node);
	Real* f = INC_FORCE.Get_Scratch(node);
#ifdef CONTACT_DEBUG_NODE
	if( PRINT_THIS_NODE ){
	  postream << "  Force for node " << node->Global_ID() << " = "
		   << f[0] << " " << f[1] << " " << f[2] << "\n";
	}
#endif
	assembled_force[0] += f[0];
	assembled_force[1] += f[1];
	assembled_force[2] += f[2];
      }
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE ){
	postream << "  Total Force = " << assembled_force[0] << " "
		 << assembled_force[1] << " " << assembled_force[2] << "\n";
      }
#endif
      for(int j=0 ; j<sh->Num_Acme_Nodes_for_Host_Node(i) ; ++j ){
        int hi, lo;
        sh->Acme_NodeGID_for_Host_Node(i,j,hi,lo);
        ContactHostGlobalID gid(hi,lo);
        ContactNode<Real>* node = static_cast<ContactNode<Real>*>
          (topology->NodeList()->Find(gid));
        POSTCONDITION(node);
	Real* f = INC_FORCE.Get_Scratch(node);
	f[0] = assembled_force[0];
	f[1] = assembled_force[1];
	f[2] = assembled_force[2];
      }
    }
  }
#ifdef CONTACT_DEBUG_NODE
  postream.flush();
#endif
}

void ContactTDEnforcement::Compute_Extra_Ghosting_Length()
{
#ifdef CONTACT_DEBUG_NODE
  ContactParOStream& postream = ParOStream();
#endif

#ifndef CONTACT_CALC_REMAINING_GAP_IN_SEARCH
  VariableHandle REMAINING_GAP  = topology->Variable_Handle( ContactTopology::Remaining_Gap );
  VariableHandle NODE_GHOST_GAP = topology->Variable_Handle( ContactTopology::Node_Ghost_Gap );
  for(int i=0 ; i<number_of_nodes ; ++i ){
    ContactNode<Real>* node = enforcement_node_list[i];
    Real* remaining_gap = node->Variable( REMAINING_GAP  );
    Real* ghost_gap     = node->Variable( NODE_GHOST_GAP );
    for( int j=0 ; j<dimensionality ; ++j ) {
      remaining_gap[j] = 0.0;
      ghost_gap[j]     = 0.0;
    }
  }
#endif
  for(Node_Constraint_Group cnei_group = node_entity_list.Node_Group_Start(); 
      cnei_group.Valid_Group();
      cnei_group = node_entity_list.Node_Group_Next()) {
    const int node_index = cnei_group.Get_Constrained_Node_Index();
    Real* node_pos = NEW_POSITION.Get_Scratch(node_index);
    for(int i = 0; i<cnei_group.Num_Interactions(); ++i) {
      Real contact_point[3];
      contact_point[0] = 0.0;
      contact_point[1] = 0.0;
      contact_point[2] = 0.0;
      const int num_face_nodes = cnei_group.Get_Num_Face_Nodes(i);
      for(int j=0 ; j<num_face_nodes ; ++j ){
        const int face_node_index = cnei_group.Get_Face_Node_Index(i, j);
        const Real face_shape_function = cnei_group.Get_Face_Shape_Function(i, j);
        Real* coord = NEW_POSITION.Get_Scratch(face_node_index);
        contact_point[0] += face_shape_function*coord[0];
        contact_point[1] += face_shape_function*coord[1];
        contact_point[2] += face_shape_function*coord[2];
      }
      //
      //  NKC, again this calculation only performed for node-face, will be skipped by node surface.
      //  again probally wrong
      //
      if(num_face_nodes > 0) {
        Real vec[3];
        vec[0] = contact_point[0] - node_pos[0];
        vec[1] = contact_point[1] - node_pos[1];
        vec[2] = contact_point[2] - node_pos[2];
#ifdef CONTACT_DEBUG_NODE
	{ 
          ContactNodeFaceInteraction * cnfi = dynamic_cast<ContactNodeFaceInteraction*>(cnei_group.Get_Interaction(i));
	  if (!cnfi) continue;
	  ContactNode<Real>* node = cnfi->Node();
	  ContactFace<Real>* face = cnfi->Face();
	  if( topology->Is_a_Debug_Node( node ) ){
	    postream << "Compute Extra Ghost Len for Node " 
	  	   << node->Exodus_ID() << " and face " 
		   << face->Global_ID() << "\n";
	    postream << "   Node Position = " << node_pos[0] << " " 
		   << node_pos[1] << " " << node_pos[2] << "\n";
	    postream << "   Contact Point = " << contact_point[0] << " "
		   << contact_point[1] << " " << contact_point[2] << "\n";
    	  }
        }
#endif
#ifndef CONTACT_CALC_REMAINING_GAP_IN_SEARCH
        ContactNodeFaceInteraction *cnfi = dynamic_cast<ContactNodeFaceInteraction*>(cnei_group.Get_Interaction(i));
        if( cnfi->Scalar_Var(ContactNodeFaceInteraction::GAP_CUR) <= 0 ){
          Real* remaining_gap = cnfi->Node()->Variable( REMAINING_GAP );
          Real* n_dir         = cnfi->Vector_Var(ContactNodeFaceInteraction::PUSHBACK_DIR);
          // Compute the normal gap not enforced
          Real pmagl = n_dir[0]*vec[0] + n_dir[1]*vec[1] + n_dir[2]*vec[2]; 
          // Keep the maximum distance in each coordinate
          remaining_gap[0] = std::max( remaining_gap[0], std::fabs(pmagl*n_dir[0]) );
          remaining_gap[1] = std::max( remaining_gap[1], std::fabs(pmagl*n_dir[1]) );
          remaining_gap[2] = std::max( remaining_gap[2], std::fabs(pmagl*n_dir[2]) );
        }
        Real* ghost_gap = cnfi->Node()->Variable( NODE_GHOST_GAP );
        ghost_gap[0]    = std::max(ghost_gap[0], std::fabs(vec[0]));
        ghost_gap[1]    = std::max(ghost_gap[1], std::fabs(vec[1]));
        ghost_gap[2]    = std::max(ghost_gap[2], std::fabs(vec[2]));
#endif
      }
    }
  }
  //
  //  Compute remaining gap on all extra nodes
  //
  for(int iconst = 0; iconst < extra_constraints.size(); ++iconst) { 
    ContactNodeEntityInteraction *cnei = extra_constraints[iconst];
    ContactNodeFaceInteraction *cnfi = dynamic_cast<ContactNodeFaceInteraction*>(cnei);
    if(cnfi == NULL) continue;
    ContactNode<Real> *node = cnei->Node();
    ContactFace<Real> *face = cnfi->Face();
    const int node_index = node->EnfArrayIndex();
    Real* node_pos = NEW_POSITION.Get_Scratch(node_index);
    Real contact_point[3];
    contact_point[0] = 0.0;
    contact_point[1] = 0.0;
    contact_point[2] = 0.0;
    const int num_face_nodes = face->Nodes_Per_Face();
    Real face_shape_function[MAX_NODES_PER_FACE];
    Real *coordinates = cnfi->Vector_Var(ContactNodeFaceInteraction::COORDINATES);
    face->Evaluate_Shape_Functions(coordinates, face_shape_function);    
    for(int j=0 ; j<num_face_nodes ; ++j ){
      ContactNode<Real> *face_node = face->Node(j);
      Real *coord = NEW_POSITION.Get_Scratch(face_node->EnfArrayIndex());
      contact_point[0] += face_shape_function[j]*coord[0];
      contact_point[1] += face_shape_function[j]*coord[1];
      contact_point[2] += face_shape_function[j]*coord[2];
    }
    //
    //  NKC, again this calculation only performed for node-face,
    //  will be skipped by node surface. again probally wrong
    //
    if(num_face_nodes > 0) {
      Real vec[3];
      vec[0] = contact_point[0] - node_pos[0];
      vec[1] = contact_point[1] - node_pos[1];
      vec[2] = contact_point[2] - node_pos[2];
#ifdef CONTACT_DEBUG_NODE
      if( topology->Is_a_Debug_Node( node ) ){
	cout << "Compute Extra Ghost Len for Node " 
	     << node->Exodus_ID() << " and face " 
	     << face->Global_ID() << "\n";
        cout << "   Node Position = " << node_pos[0] << " " 
	     << node_pos[1] << " " << node_pos[2] << "\n";
        cout << "   Contact Point = " << contact_point[0] << " "
	     << contact_point[1] << " " << contact_point[2] << "\n";
      }
#endif
#ifndef CONTACT_CALC_REMAINING_GAP_IN_SEARCH
      Real* ghost_gap = cnfi->Node()->Variable( NODE_GHOST_GAP );
      ghost_gap[0]    = std::max(ghost_gap[0], std::fabs(vec[0]));
      ghost_gap[1]    = std::max(ghost_gap[1], std::fabs(vec[1]));
      ghost_gap[2]    = std::max(ghost_gap[2], std::fabs(vec[2]));
#endif
    }
  }
#ifdef CONTACT_DEBUG_NODE
  postream.flush();
#endif
}

void ContactTDEnforcement::Compact_Node_Entity_List() {
  //
  //  Compact the global constraint lists.  Combine any colinear constraints, eliminate multiple tied interactions,
  //  Remove interactions with a zero kinematic partition.
  //
  have_multiple = 0;
  //
  //  Extra constraints will hold any defined interactions that are thrown out by this routine.  Remainig gaps must
  //  still be computed on this interactions to ensure correct ghosting.
  //
  extra_constraints.clear();

  for(Node_Constraint_Group cnei_list = node_entity_list.Node_Group_Start();
      cnei_list.Valid_Group();
      cnei_list = node_entity_list.Node_Group_Next()) {
    int ncc = cnei_list.Num_Interactions();
    Real *constraint_normals[3];
    Real constraint_gaps[3];
    Real kin_part[3];
    ContactNodeEntityInteraction *constraint_list[3];
    int num_constraints = 0;
    //
    //  Extract information about each constraint
    //
    for(int j=0 ; j<ncc ; ++j ){
      ContactNodeEntityInteraction  *cnei = cnei_list.Get_Interaction(j);
      //
      // set kinematic partition based on value
      // if two calculate using wavespeed/density method
      // else use fixed value
      //  
      Real kp = Enforcement_Data(KINEMATIC_PARTITION, cnei);
      if( kp+0.25 >= 2.0 && kp+0.25 < 3.0) {
        //we should never get here if we have a nsi, so only doing a static
        //cast. However, Just in case something is broken make a check during debug.
        PRECONDITION(dynamic_cast<ContactNodeFaceInteraction*>(cnei) != NULL);
        kp = Auto_Kinematic_Partition(static_cast<ContactNodeFaceInteraction*>(cnei));
      }

      if(kp > 0.0 && TDEnfModel(cnei)->Active_Interaction(cnei, cnei->Get_Gap_Cur())) {
        constraint_normals[num_constraints] = cnei->Get_Normal();
        constraint_gaps[num_constraints] = cnei->Get_Gap_Cur();
        kin_part[num_constraints] = kp;
        constraint_list[num_constraints] = cnei;
        ++num_constraints;
      } else {
        extra_constraints.push_back(cnei);
      }
    }
    //
    //  Determine which constraints to keep
    //
    bool keep[3] = {true, true, true};
    if(num_constraints > 1) {
      Real dot = 0.0;
      for(int j=0; j<NDIM; ++j) {dot += constraint_normals[0][j] * constraint_normals[1][j]; }
      if (std::fabs(dot) > COLLINEARITY_TOL) {
        if (constraint_gaps[0] > constraint_gaps[1]) {
          keep[1] = false;
        } else {
          keep[0] = false;
        }
      }
      if(num_constraints == 3) {
        dot = 0.0;
        for(int j=0 ; j<NDIM ; ++j ){dot += constraint_normals[2][j]*constraint_normals[1][j]; }
        if (std::fabs(dot) > COLLINEARITY_TOL) {
          if (constraint_gaps[2] > constraint_gaps[1]) {
            keep[1] = false;
          } else {
            keep[2] = false;
          }
        }
        dot = 0.0;
        for(int j=0 ; j<NDIM ; ++j ){dot += constraint_normals[0][j]*constraint_normals[2][j]; }
        if (std::fabs(dot) > COLLINEARITY_TOL) {
          if (constraint_gaps[0] > constraint_gaps[2]) {
            keep[2] = false;
          } else {
            keep[0] = false;
          }
        }
      }
    }
    bool have_tied = false;
    for(int j = 0; j < num_constraints; ++j) {
      if(keep[j]) {
        //
        //  Eliminate multiple tied constraints if they occur
        //
        ContactNodeEntityInteraction  *cnei = constraint_list[j];
        if(TDEnfModel(cnei)->Needs_Glued_Search(TDEnfModel(cnei)->Interaction_Type(cnei)) )  {
          if (!(cnei->Is_Tied()||cnei->Is_InfSlip())) cnei->Is_Glued( true );     
          if(have_tied) {
            keep[j] = false;
	  } else {
            have_tied = true;
	  }
	}
      }
    }
    //
    //  Compact the lists (set any invalid interaction to null, and then compact), 
    //  set the values for GAP_TO_ENFORCE and the KINEMATIC_PARTITION
    //  variables
    //
    int num_kept = 0;
    for(int j = 0; j < num_constraints; ++j) {
      if(keep[j]) {
        cnei_list.Set_Interaction(num_kept,constraint_list[j]);
        int index = cnei_list.Get_Index(num_kept);
        *(NEI_GAP_TO_ENFORCE.     Get_Scratch(index)) = constraint_gaps[j] ;
        *(NEI_KINEMATIC_PARTITION.Get_Scratch(index)) = kin_part[j];
        ++num_kept;
      } else {
        extra_constraints.push_back(constraint_list[j]);
      }
    }

    for(int j = num_kept; j < ncc; ++j) {
      cnei_list.Set_Interaction(j, NULL);
    }
    if(num_kept > 1) have_multiple = 1; 
  }
  have_multiple = contact_global_sum(have_multiple, communicator);
}

bool
ContactTDEnforcement::Only_Frictionless_or_Tied(void) {
  int n = enforcement_data.Number_Entity_Keys();
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      int fric_id = (int) enforcement_data.Get_Data( FRICTION_MODEL_ID, i, j );
      ContactEnfModel *model = enforcement_models[fric_id];
      ContactTDTied         *tied_model         = dynamic_cast<ContactTDTied*>(model);      
      ContactTDFrictionless *frictionless_model = dynamic_cast<ContactTDFrictionless*>(model);      
      if(tied_model != NULL || frictionless_model != NULL) {
        continue;
      } else {
        return false;
      }
    }
  }
  return true;
}
