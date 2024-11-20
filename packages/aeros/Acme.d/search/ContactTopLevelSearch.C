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


#include <algorithm>

#include "allocators.h"
#include "contact_assert.h"
#include "ContactAnalyticCylinderInside.h"
#include "ContactAnalyticCylinderOutside.h"
#include "ContactAnalyticPlane.h"
#include "ContactAnalyticSphere.h"
#include "ContactAnalyticSurface.h"
#include "ContactCommBuffer.h"
#include "Contact_Communication.h"
#include "ContactEdgeBlock.h"
#include "ContactEnforcement.h"
#include "ContactErrors.h"
#include "ContactFaceBlock.h"
#include "ContactNodeNodeInteraction.h"
#include "ContactNodeFaceInteraction.h"
#include "ContactNodeSurfaceInteraction.h"
#include "ContactFaceFaceInteraction.h"
#include "ContactElementElementInteraction.h"
#include "ContactElement.h"
#include "ContactLineEdgeL2.h"
#include "ContactTriFaceL3.h"
#include "ContactQuadFaceL4.h"
#include "ContactHexElementL8.h"
#include "ContactWedgeElementL6.h"
#include "ContactShellNode.h"
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
#include "ContactRangeSearch.h"
#include "contact_tolerances.h"
#include "ContactEnforcement.h"
#include "ContactTDEnforcement.h"

#ifndef CONTACT_NO_MPI
#include "mpi.h"
#include "ContactZoltan.h"
#include "ContactZoltanCommUtils.h"
#endif

#include <iostream>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>

using namespace std;

ContactSearch::ContactErrorCode ContactSearch::Static_Search_1_Configuration()
{
  Real dt           = 0.0;
  Real dt_old       = 0.0;
  if ( initialized_tied ) {
    initializing_tied = false;
  } else {
    initializing_tied = true;
  }
  Initialize_Context();
  TopLevel_Search(STATIC1CONFIG, dt_old, dt);
  return (ContactErrorCode) contact_global_error_check( error_code, SearchComm );
}

ContactSearch::ContactErrorCode ContactSearch::Static_Search_2_Configuration()
{
  Real dt     = 0.0;
  Real dt_old = 0.0;
  if (search_data->Have_Tied_Interactions() && !initialized_tied) {
    // initialize all tied interactions
    ContactErrorCode err = Static_Search_1_Configuration_Tied();
    if (err != NO_ERROR) return err;
  } else {
    Initialize_Context();
  }
  // go about my normal business...
  TopLevel_Search(STATIC2CONFIG, dt_old, dt);
  return (ContactErrorCode) contact_global_error_check( error_code, SearchComm );
}

ContactSearch::ContactErrorCode ContactSearch::Dynamic_Search_2_Configuration( Real& dt_old, Real& dt )
{
  if (search_data->Have_Tied_Interactions() && !initialized_tied) {
    // initialize all tied interactions
    ContactErrorCode err = Static_Search_1_Configuration_Tied();
    if (err != NO_ERROR) return err;
  } else {
    Initialize_Context();
  }
  // go about my normal business...
  TopLevel_Search(DYNAMIC2CONFIG, dt_old, dt);
  return (ContactErrorCode) contact_global_error_check( error_code, SearchComm );
}

ContactSearch::ContactErrorCode ContactSearch::Static_Search_1_Configuration_Tied()
{
  error_code = NO_ERROR;
  Initialize_Context();
  // all this is for initializing the gap based on a static configuration
  if (search_data->Have_Tied_Interactions() && !initialized_tied) {
    initializing_tied = true;
    search_data->SetOnlyTied();
    // lofted shells will need to have their positions stored and then restored
    Real* cur_pos = NULL;
    Real* pre_pos = NULL;
    if (primary_topology->Have_Shells()) {
      int num_primary_nodes             = primary_topology->Number_of_Primary_Nodes();
      ContactNode<Real>** Nodes               = reinterpret_cast<ContactNode<Real>**>(primary_topology->PrimaryNodeList()->EntityList());
      VariableHandle CURRENT_POSITION   = primary_topology->Variable_Handle( ContactTopology::Current_Position );
      VariableHandle PREDICTED_POSITION = primary_topology->Variable_Handle( ContactTopology::Predicted_Position );
      cur_pos = new Real[3*num_primary_nodes];
      pre_pos = new Real[3*num_primary_nodes];
      for (int i=0; i<num_primary_nodes; ++i) {
        Real* pos_c    = Nodes[i]->Variable( CURRENT_POSITION );
        Real* pos_p    = Nodes[i]->Variable( PREDICTED_POSITION );
        cur_pos[3*i+0] = pos_c[0];
        cur_pos[3*i+1] = pos_c[1];
        cur_pos[3*i+2] = pos_c[2];
        pre_pos[3*i+0] = pos_p[0];
        pre_pos[3*i+1] = pos_p[1];
        pre_pos[3*i+2] = pos_p[2];
      }
    }
    // saving options & settings
    Search_Option_Status save_multiple_interaction_status = multiple_interaction_status;
    Search_Option_Status save_enable_tracking             = enable_tracking;
    Track_Type           save_tracking_type               = tracking_type;
    multiple_interaction_status = INACTIVE;
    enable_tracking             = INACTIVE;
    tracking_type               = NO_TRACKING;
    Real tmp_dt                 = 0.0;
    Real tmp_dt_old             = 0.0;
    // a static 1 configuration search to set initial gaps
    TopLevel_Search(STATIC1CONFIG_TIED, tmp_dt_old, tmp_dt);
    ContactErrorCode err = 
      (ContactErrorCode) contact_global_error_check( error_code, SearchComm );
    if (err != NO_ERROR) return err;
    Set_Initial_Gap();
    // restore options & settings
    multiple_interaction_status = save_multiple_interaction_status;
    enable_tracking             = save_enable_tracking;
    tracking_type               = save_tracking_type;
    step_number                 = 0;
    tracking_step               = 0;
    // restore lofted shell positions
    if (primary_topology->Have_Shells()) {
      int num_primary_nodes = primary_topology->Number_of_Primary_Nodes();
      ContactNode<Real>** Nodes   = reinterpret_cast<ContactNode<Real>**>(primary_topology->PrimaryNodeList()->EntityList());
      VariableHandle CURRENT_POSITION   = primary_topology->Variable_Handle( ContactTopology::Current_Position );
      VariableHandle PREDICTED_POSITION = primary_topology->Variable_Handle( ContactTopology::Predicted_Position );
      for (int i=0; i<num_primary_nodes; ++i) {
        Real* pos_c = Nodes[i]->Variable( CURRENT_POSITION );
        Real* pos_p = Nodes[i]->Variable( PREDICTED_POSITION );
        pos_c[0]    = cur_pos[3*i+0];
        pos_c[1]    = cur_pos[3*i+1];
        pos_c[2]    = cur_pos[3*i+2];
        pos_p[0]    = pre_pos[3*i+0];
        pos_p[1]    = pre_pos[3*i+1];
        pos_p[2]    = pre_pos[3*i+2];
      }
      delete [] cur_pos;
      delete [] pre_pos;
      cur_pos = NULL;
      pre_pos = NULL;
    }
    search_data->SetAllButTied();
    #ifndef CONTACT_NO_MPI
    primary_topology->GhostTiedFaces();
    #endif
    initialized       = false;
    initializing_tied = false;
    initialized_tied  = true;
  }
  return (ContactErrorCode) contact_global_error_check( error_code, SearchComm );
}
 
void
ContactSearch::TopLevel_Search( SearchType search_type, Real& dt_old, Real& dt )
{
  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Start_Timer( contact_time );
  timer.Start_Timer( search_time );
  timer.Start_Timer( search_setup_time );
  #endif
  
  if (!initialized || restart) {
    Initialize_Search();
  }
  
  error_code = NO_ERROR;
  
  //------------------------------------------------------------------------------
  // Determine the number of configurations
  //------------------------------------------------------------------------------
  int num_configs = 0;
  ContactTDEnforcement* td_enf = NULL;
  switch (search_type) {
  case STATIC1CONFIG:
  case STATIC1CONFIG_TIED:
    num_configs = 1;
    break;
  case STATIC2CONFIG:
    num_configs = 2;
    break;
  case DYNAMIC2CONFIG:
    if( old_dynamic_search==INACTIVE ) {
      num_configs = 3;
      if( !enforcement ){
        errors->Add_Error_Message("Dynamic_Search_2_Configuration requires that ");
        errors->Add_Error_Message("an ACME (Dynamic) Enforcement be instantiated when the");
        errors->Add_Error_Message("AUGMENTED option is used.");
        error_code = UNKNOWN_TYPE;
        return;
      }
      
      // Look for a Transient Dynamics Enforcement
      
      bool found_enforcement = false;
      for( int i=0 ; i<number_registered_enforcements ; ++i ){
        if( enforcement[i]->Type() == ContactEnforcement::TDEnforcement 
         || enforcement[i]->Type() == ContactEnforcement::TDEnfPenalty ){
          td_enf = ( ContactTDEnforcement* ) enforcement[i];
          found_enforcement = true;
          break;
        }
      }
      
      if( !found_enforcement ){
        errors->Add_Error_Message("Dynamic_Search_2_Configuration reguires that ");
        errors->Add_Error_Message("an ACME Dynamic Enforcement be instantiated when the");
        errors->Add_Error_Message("AUGMENTED option is used.");
        error_code = UNKNOWN_TYPE;
        return;
      }
    } else {
      num_configs = 2;
    }
  }
  nconfigs = num_configs;
  
  //------------------------------------------------------------------------------
  // Print some info at the start of the search (if requested)
  //------------------------------------------------------------------------------
  Print_Search_Header(search_type, dt, dt_old, td_enf);
  
  //------------------------------------------------------------------------------
  // Advance the state of the interactions
  //------------------------------------------------------------------------------
  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Start_Timer( search_update_state_time );
  #endif
  primary_topology->Update_State();
  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Stop_Timer( search_update_state_time );
  #endif

  //------------------------------------------------------------------------------
  // Set the configuration variable handles
  //------------------------------------------------------------------------------
  // The variable handles are the same in either decomposition so use primary
  VariableHandle CURRENT_POSITION   = -1;
  VariableHandle PREDICTED_POSITION = -1;
  VariableHandle AUGMENTED_POSITION = -1;
  switch (search_type) {
  case STATIC1CONFIG:
  case STATIC1CONFIG_TIED:
    CURRENT_POSITION   = primary_topology->Variable_Handle( ContactTopology::Current_Position );
    PREDICTED_POSITION = primary_topology->Variable_Handle( ContactTopology::Current_Position );
    AUGMENTED_POSITION = primary_topology->Variable_Handle( ContactTopology::Current_Position );
    break;
  case STATIC2CONFIG:
    CURRENT_POSITION   = primary_topology->Variable_Handle( ContactTopology::Current_Position   );
    PREDICTED_POSITION = primary_topology->Variable_Handle( ContactTopology::Predicted_Position );
    AUGMENTED_POSITION = primary_topology->Variable_Handle( ContactTopology::Predicted_Position );
    break;
  case DYNAMIC2CONFIG:
    CURRENT_POSITION   = primary_topology->Variable_Handle( ContactTopology::Current_Position   );
    PREDICTED_POSITION = primary_topology->Variable_Handle( ContactTopology::Predicted_Position );
    AUGMENTED_POSITION = primary_topology->Variable_Handle( ContactTopology::Augmented_Position );
    break;
  default:
    POSTCONDITION(0);
    break;  
  }
  POSTCONDITION(CURRENT_POSITION>=0);
  POSTCONDITION(PREDICTED_POSITION>=0);
  POSTCONDITION(AUGMENTED_POSITION>=0);
  //------------------------------------------------------------------------------
  // Set some misc data
  //------------------------------------------------------------------------------
  Real min_characteristic_length = primary_topology->Min_Characteristic_Length();
  
  //------------------------------------------------------------------------------
  // Set the search options in the topology
  //------------------------------------------------------------------------------
  primary_topology->Set_Search_Options( multiple_interaction_status,
					normal_smoothing_status,
					sharp_smooth_curvature,
                                        no_parallel_consistency,
                                        auto_tol );
  
  //------------------------------------------------------------------------------
  // Set the global tolerances
  //------------------------------------------------------------------------------
  Real user_search_tol = search_data->Max_Search_Tolerance();
  Real user_tang_tol   = search_data->Max_Search_Tangential_Tolerance();
  Real gap_tol         = 0.0;
  Real gap_tol_auto    = 0.0;
  Real gap_tol_user    = 0.0;
  Real capture_fudge   = 1.1;
  
  switch (search_type) {
  case STATIC1CONFIG:
  case STATIC1CONFIG_TIED:
    max_node_motion[0]    = 0.0;
    max_node_motion[1]    = 0.0;
    max_node_motion[2]    = 0.0;
    max_node_displacement = 0.0;
    max_remaining_gaps[0] = 0.0;
    max_remaining_gaps[1] = 0.0;
    max_remaining_gaps[2] = 0.0;
    max_remaining_gap_mag = 0.0;
    reasonable_gap        = 0.0;
    capture_motion        = 0.0;
    gap_tol               = user_search_tol;
    break;
  case STATIC2CONFIG:
    primary_topology->Compute_Max_Relative_Node_Motion( max_node_motion );
    max_node_displacement = std::sqrt( max_node_motion[0]*max_node_motion[0]+ 
  				       max_node_motion[1]*max_node_motion[1]+ 
  				       max_node_motion[2]*max_node_motion[2]);
    //max_node_displacement *= (Real)(std::max(1,(global_tracking_interval-1))); 
    capture_motion = capture_fudge*max_node_displacement*global_tracking_interval;
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    timer.Start_Timer( search_remaining_gap_time );
    #endif
    Compute_Remaining_Gap( &max_remaining_gap_mag, &(max_remaining_gaps[0]) );
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    timer.Stop_Timer( search_remaining_gap_time );
    #endif
    if (auto_tol) {
      reasonable_gap  = max_remaining_gap_mag + max_node_displacement;
      reasonable_gap *= gap_inflation;
      gap_tol         = std::max(reasonable_gap,user_search_tol);
    } else {
      reasonable_gap = 0.0;
      gap_tol        = max_node_displacement + max_remaining_gap_mag + user_search_tol;
    }
    break;
  case DYNAMIC2CONFIG:
    Real auto_norm_tol = 0.0;
    Real auto_tang_tol = 0.0;
    primary_topology->Compute_Max_Relative_Node_Motion( max_node_motion );
    max_node_displacement = std::sqrt( max_node_motion[0]*max_node_motion[0]+ 
  				       max_node_motion[1]*max_node_motion[1]+ 
  				       max_node_motion[2]*max_node_motion[2]);
    capture_motion = capture_fudge*max_node_displacement*global_tracking_interval;
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    timer.Start_Timer( search_remaining_gap_time );
    #endif
    Compute_Remaining_Gap( &max_remaining_gap_mag, &(max_remaining_gaps[0]) );
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    timer.Stop_Timer( search_remaining_gap_time );
    #endif
    if (auto_tol) {

      // restrict the value of max_node_displacement (in case of errant deformation)
      max_node_displacement = std::min( max_node_displacement,min_characteristic_length );

      // automatically computed 'normal search tolerance'
      reasonable_gap  = max_remaining_gap_mag + max_node_displacement;
            //    reasonable_gap  += max_node_displacement;
      auto_norm_tol = gap_inflation * reasonable_gap;
      // automatically computed 'tangential search tolerance'
      auto_tang_tol = 0.01 * min_characteristic_length;
      if( user_tang_tol>1.0e-20 ){
        // the user has specified a tangential tolerance, so use it
        auto_tang_tol = user_tang_tol;
      }

      gap_tol_auto    = std::max(auto_norm_tol,auto_tang_tol);
      // finally consider consider the larger of the user vs. auto search tolerances
      gap_tol_auto = std::max(gap_tol_auto,user_search_tol);
      
      gap_tol_user = max_node_displacement + max_remaining_gap_mag + user_search_tol;
      gap_tol = gap_tol_auto;
    } else {
      // automatically computed 'normal search tolerance' with restricted value of
      // max_node_displacement (in case of errant deformation)
      reasonable_gap  = max_remaining_gap_mag + std::min( max_node_displacement,min_characteristic_length );
          //    reasonable_gap  += max_node_displacement;
      auto_norm_tol = gap_inflation * reasonable_gap;
      // automatically computed 'tangential search tolerance'
      auto_tang_tol = 0.01 * min_characteristic_length;
      if( user_tang_tol>1.0e-20 ){
        // the user has specified a tangential tolerance, so use it
        auto_tang_tol = user_tang_tol;
      }

      gap_tol_auto = std::max(auto_norm_tol,auto_tang_tol);
      // finally consider consider the larger of the user vs. auto search tolerances
      gap_tol_auto = std::max(gap_tol_auto,user_search_tol);

      gap_tol_user = max_node_displacement + max_remaining_gap_mag + user_search_tol;
      gap_tol      = gap_tol_user;
    }
    break;
  }
  #if CONTACT_DEBUG_PRINT_LEVEL>=2
  postream<<"    max_remaining_gap     = "<<max_remaining_gaps[0]<<", "<<max_remaining_gaps[1]<<", "<<max_remaining_gaps[2]<<"\n";
  postream<<"    max_remaining_gap_mag = "<<max_remaining_gap_mag<<"\n";
  postream<<"    reasonable_gap        = "<<reasonable_gap<<"\n";
  postream<<"    user tolerance gap    = "<<user_search_tol<<"\n";
  postream<<"    gap tolerance         = "<<gap_tol<<"\n";
  postream<<"    gap tolerance user    = "<<gap_tol_user<<"\n";
  postream<<"    gap tolerance auto    = "<<gap_tol_auto<<"\n";
  postream<<"    min char length       = "<<min_characteristic_length<<"\n";
  postream<<"    max char length       = "<<primary_topology->Max_Characteristic_Length()<<"\n";
  postream<<"    max node displacement = "<<max_node_displacement<<"\n";
  postream<<"    capture motion        = "<<capture_motion<<"\n";
  #else
  #ifdef CONTACT_HEARTBEAT 
  if (contact_processor_number(SearchComm)==0) {
    std::cout<<"    max_remaining_gap	  = "<<max_remaining_gaps[0]<<", "<<max_remaining_gaps[1]<<", "<<max_remaining_gaps[2]<<"\n";
    std::cout<<"    max_remaining_gap_mag = "<<max_remaining_gap_mag<<"\n";
    std::cout<<"    reasonable_gap	  = "<<reasonable_gap<<"\n";
    std::cout<<"    user tolerance gap    = "<<user_search_tol<<"\n";
    std::cout<<"    gap tolerance	  = "<<gap_tol<<"\n";
    std::cout<<"    gap tolerance user    = "<<gap_tol_user<<"\n";
    std::cout<<"    gap tolerance auto    = "<<gap_tol_auto<<"\n";
    std::cout<<"    min char length       = "<<min_characteristic_length<<"\n";
    std::cout<<"    max char length       = "<<primary_topology->Max_Characteristic_Length()<<"\n";
    std::cout<<"    max node displacement = "<<max_node_displacement<<"\n";
    std::cout<<"    capture motion        = "<<capture_motion<<"\n";
  }
  #endif
  #endif
    
  #if CONTACT_DEBUG_PRINT_LEVEL>=6
  primary_topology->Display( postream );
  #endif

  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Stop_Timer( search_setup_time );
  #endif
  
  //------------------------------------------------------------------------------
  // Set the augmented position for dynamic search
  //------------------------------------------------------------------------------
  if (search_type==DYNAMIC2CONFIG) {
  
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    timer.Start_Timer( augmented_config_time );
    #endif
    
    // Construct the augmented configuration by using the old displacements
    int num_primary_nodes = primary_topology->Number_of_Primary_Nodes();
    ContactNode<Real>** Nodes = 
      reinterpret_cast<ContactNode<Real>**>(primary_topology->PrimaryNodeList()->EntityList());

    //
    // Asked to use augmented search.
    if( old_dynamic_search==INACTIVE ) {
      Real old_time_multiplier = td_enf->Get_Old_Time_Multiplier();
      if( old_time_multiplier == 0 ) {
        for (int i=0; i<num_primary_nodes; ++i) {
          ContactNode<Real>* node = Nodes[i];
          Real* pos_p = node->Variable( PREDICTED_POSITION );
          Real* pos_a = node->Variable( AUGMENTED_POSITION );
          pos_a[0] = pos_p[0]; 
          pos_a[1] = pos_p[1]; 
          pos_a[2] = pos_p[2]; 
          #ifdef CONTACT_DEBUG_NODE
	  if( primary_topology->Is_a_Debug_Node( node ) ){
            Real* pos_c = node->Variable( CURRENT_POSITION );
	    postream << "1> Positions for Node (" << node->Exodus_ID() << ")\n" 
		     << "   Current:    " << pos_c[0] << " " << pos_c[1]
		     << " " << pos_c[2] << "\n"
		     << "   Predicted:  " << pos_p[0] << " " << pos_p[1]
		     << " " << pos_p[2] << "\n"
		     << "   Augmented:  " << pos_a[0] << " " << pos_a[1]
		     << " " << pos_a[2] << "\n"
		     << "   Aug Displ:  " << 0. 
		     << " " << 0. << " " 
		     << 0. << "\n";
	  }
          #endif
        }
      } else {
        if( num_primary_nodes ) {
          Real *disp_x_old = NULL, *disp_y_old = NULL, *disp_z_old = NULL;
          td_enf->Get_Old_Displacement_Ptrs( 0, &disp_x_old);
          td_enf->Get_Old_Displacement_Ptrs( 1, &disp_y_old);
          td_enf->Get_Old_Displacement_Ptrs( 2, &disp_z_old);
          Real ratio_of_multipliers = (dt*(0.5*(dt+dt_old)))/old_time_multiplier;
          // Assume acceleration is the same as last time-step
          // to predict augmented configuration
          //postream<<"  updating augmented position for "<<num_primary_nodes<<" nodes\n";
          //postream.flush();
          for (int i=0; i<num_primary_nodes; ++i) {
            ContactNode<Real>* node = Nodes[i];
            Real* pos_p = node->Variable( PREDICTED_POSITION );
            Real* pos_a = node->Variable( AUGMENTED_POSITION );
            Real old_disp_x = disp_x_old[i];
            Real old_disp_y = disp_y_old[i];
            Real old_disp_z = disp_z_old[i];
            pos_a[0]    = pos_p[0] + ratio_of_multipliers*old_disp_x; 
            pos_a[1]    = pos_p[1] + ratio_of_multipliers*old_disp_y; 
            pos_a[2]    = pos_p[2] + ratio_of_multipliers*old_disp_z; 
            #ifdef CONTACT_DEBUG_NODE
            if( primary_topology->Is_a_Debug_Node( node ) ){
              Real* pos_c = node->Variable( CURRENT_POSITION );
              postream << "2> Positions for Node (" << node->Exodus_ID() << ")\n" 
                       << "   Current:    " << pos_c[0] << " " << pos_c[1]
                       << " " << pos_c[2] << "\n"
                       << "   Predicted:  " << pos_p[0] << " " << pos_p[1]
                       << " " << pos_p[2] << "\n"
                       << "   Augmented:  " << pos_a[0] << " " << pos_a[1]
                       << " " << pos_a[2] << "\n"
                       << "   Aug Displ:  " << ratio_of_multipliers*disp_x_old[i] 
                       << " " << ratio_of_multipliers*disp_y_old[i] << " " 
                       << ratio_of_multipliers*disp_z_old[i] << "\n";
            }
            #endif
          }
        }
      }
    } else { 
      // use old dynamic search method 
      for (int i=0; i<num_primary_nodes; ++i) {
        ContactNode<Real>* node = Nodes[i];
        Real* pos_p = node->Variable( PREDICTED_POSITION );
        Real* pos_a = node->Variable( AUGMENTED_POSITION );
        pos_a[0] = pos_p[0]; 
        pos_a[1] = pos_p[1]; 
        pos_a[2] = pos_p[2];
        #ifdef CONTACT_DEBUG_NODE
        if( primary_topology->Is_a_Debug_Node( node ) ){
          Real* pos_c = node->Variable( CURRENT_POSITION );
          postream << "3> Positions for Node (" << node->Exodus_ID() << ")\n" 
                   << "   Current:    " << pos_c[0] << " " << pos_c[1]
                   << " " << pos_c[2] << "\n"
                   << "   Predicted:  " << pos_p[0] << " " << pos_p[1]
                   << " " << pos_p[2] << "\n"
                   << "   Augmented:  " << pos_a[0] << " " << pos_a[1]
                   << " " << pos_a[2] << "\n"
                   << "   Aug Displ:  " << 0.0 << " " << 0.0 << " " 
                   << 0.0 << "\n";
        }
        #endif
      }
    }

    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    timer.Stop_Timer( augmented_config_time );
    #endif
  }
  
  #if defined(CONTACT_TIMINGS)
  timer.Start_Timer( context_time );
  #endif
  
  SetContextForSearchSlaves();
    
  SetContextForGeometryUpdate();
  
  SetContextForGhostingRCB();
  
  #if defined(CONTACT_TIMINGS)
  timer.Stop_Timer( context_time );
  #endif

  //------------------------------------------------------------------------------
  // Set some geometry values
  //------------------------------------------------------------------------------
  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Start_Timer( geom_time );
  #endif
  if (do_elem_elem_search) {
    primary_topology->Compute_Element_Geometry( AUGMENTED_POSITION );
  }
  if (do_face_face_search && !do_node_face_search) {
    primary_topology->Compute_Face_Geometry( AUGMENTED_POSITION, false, true );
  }
  if (do_node_face_search) {
    if (tracking_type==GLOBAL_TRACKING && tracking_step>0) {
      primary_topology->Compute_Surface_Geometry( AUGMENTED_POSITION, 3, true );
    } else {
      primary_topology->Compute_Surface_Geometry( AUGMENTED_POSITION, 3, false );
    }
  }
  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Stop_Timer( geom_time );
  #endif
  
  if (tracking_type==GLOBAL_TRACKING) {
    //------------------------------------------------------------------------------
    // Global tracked search step
    //------------------------------------------------------------------------------
    NewGlobalSearch(search_type, 
                    gap_tol, num_configs,
                    CURRENT_POSITION, 
                    PREDICTED_POSITION,
                    AUGMENTED_POSITION);
      if (error_code != NO_ERROR) return;
  } else {
    if(tracking_step > 0){
      //------------------------------------------------------------------------------
      // Local tracked search step
      //------------------------------------------------------------------------------
      TrackedSearch(search_type, 
                    gap_tol, num_configs,
                    CURRENT_POSITION, 
                    PREDICTED_POSITION,
                    AUGMENTED_POSITION);
      if (error_code != NO_ERROR) return;

    }
    if(tracking_step==0){
      //------------------------------------------------------------------------------
      // global search step
      //------------------------------------------------------------------------------
      GlobalSearch(search_type, 
                   gap_tol, num_configs,
                   CURRENT_POSITION, 
                   PREDICTED_POSITION,
                   AUGMENTED_POSITION);
      if (error_code != NO_ERROR) return;
    }
  }
  if (do_coverage_search) {
    Process_Face_Coverage();
  }
  
  #ifdef CONTACT_DEBUG
  {
    // If in debug mode, do a sanity check to make sure there
    // are not too many interactions defined at each node.
    int number_of_nodes = primary_topology->Number_of_Nodes();
    ContactNode<Real>** Nodes = 
      reinterpret_cast<ContactNode<Real>**>(primary_topology->NodeList()->EntityList());
    for (int i=0; i<number_of_nodes; ++i) {
      ContactNode<Real>* node = Nodes[i];
      int ncc = node->Number_NodeEntity_Interactions();
      PRECONDITION(ncc<=MAX_NODE_ENTITY_INTERACTIONS_PER_NODE);
    }
  }
  #endif
  
  //------------------------------------------------------------------------------
  // Print some info at the start of the search (if requested)
  //------------------------------------------------------------------------------
  Print_Search_Footer();

  //------------------------------------------------------------------------------
  // All done!
  //------------------------------------------------------------------------------
  ++step_number;
  if (tracking_type!=NO_TRACKING) {
    tracking_step++;
    if (global_tracking_interval>1 && tracking_step>=global_tracking_interval) tracking_step=0;
  }
  restart = false;
  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  timer.Stop_Timer( search_time );
  timer.Stop_Timer( contact_time );
  #endif
}

int
ContactSearch::Dynamic_Process_Method( ContactNode<Real>* node, ContactFace<Real>* face )
{     
#ifdef CONTACT_DEBUG_NODE
  bool PRINT_THIS_NODE = primary_topology->Is_a_Debug_Node( node );
  if( PRINT_THIS_NODE ){
    postream << "In Dynamic_Process_Method for Node ("
	     << node->Exodus_ID() << ") with Face "
	     << face->Global_ID() << "\n";
  }
#endif
  // Decide which detailed check to use.  The logic is as follows
  //  1) If no previous interaction, use the moving search
  //  2) If previous interaction(s), check each old interaction and
  //     if any are connected to the face we're processing, use the
  //     static search.
  //  3) If previous interaction(s), and this face isn't connected
  //     to any of them, use the dynamic search. 
  int Process_Method = ContactNodeFaceInteraction::MOVING_INTERSECTION;
  int num_physical_faces = (int) *(NUMBER_PHYSICAL_FACES.Get_Scratch(node->ProcArrayIndex()));
  VariableHandle FACE_NORMAL = search_topology->Variable_Handle( ContactTopology::Face_Normal );
  // MWG: what about when num_physical_faces=0 (e.g. sph node)
  if( num_physical_faces == 1 ){
    // Only one physical face the processing is easy
    //   a) no interaction so use moving
    //   b) an interaction so use cpproj
    int node_key = node->Entity_Key();
    if( node_key == -1) { node_key = node->GetFace(0)->Entity_Key(); }
    int  type = (int) search_data->Get_Search_Data(INTERACTION_TYPE,
                                                   node_key,
                                                   face->Entity_Key());
    if( (type == CPP_INTERACTION) ||
     ((step_number != 0 || restart) && node->Get_NodeEntity_Interaction(0,STATE1))){
      Process_Method = ContactNodeFaceInteraction::CLOSEST_POINT_PROJECTION_2;
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE )
	postream << 
	  "  Choosing CPP2 with 1 physical face \n";
#endif
    }
  } else {
   if (step_number != 0 || restart) {
#ifdef CONTACT_DEBUG_NODE
     if( PRINT_THIS_NODE )
	postream << 
	  "  Processing with " << num_physical_faces << " physical faces\n";
#endif
     Real *normals[3];
     normals[0] = PHYSICAL_FACE_NORMAL_1.Get_Scratch(node->ProcArrayIndex());
     normals[1] = PHYSICAL_FACE_NORMAL_2.Get_Scratch(node->ProcArrayIndex());
     normals[2] = PHYSICAL_FACE_NORMAL_3.Get_Scratch(node->ProcArrayIndex());
     Real* face_normal = face->Variable(FACE_NORMAL);

     // First determine which physical face we will be contacting this step
     Real most_opposed = 2.0;
     int physical_face = 0;
     Find_Physical_Face(num_physical_faces, normals, face_normal, physical_face, most_opposed);
#ifdef CONTACT_DEBUG_NODE
     if( PRINT_THIS_NODE )
      postream << "  Most Opposed physical face: " << physical_face << "\n";
#endif
    
    // Second determine if one of the old interactions is most opposed to
    // the same physical face
    ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions(STATE1);
    for( int j=0 ; j<node->Number_NodeEntity_Interactions(STATE1) ; ++j ){
      ContactNodeEntityInteraction* cnei_old = interactions[j];
      if( cnei_old ){
#ifdef CONTACT_DEBUG_NODE
         if( PRINT_THIS_NODE )
	  postream << "  Processing old interaction: " << j << "\n";
#endif
         // If the old interaction was against the same face as this
	 // new pairing, use closest point projection
	 if ( cnei_old->Get_Type() == ContactNodeEntityInteraction::NODE_FACE_INTERACTION ) {
	   ContactNodeFaceInteraction* cnfi_old = dynamic_cast<ContactNodeFaceInteraction*>(cnei_old);
	   POSTCONDITION(cnfi_old);
           ContactHostGlobalID old_face_gid( cnfi_old->Get_Entity_Data()->host_gid[0], 
                                            cnfi_old->Get_Entity_Data()->host_gid[1] );
	   if (old_face_gid == face->Global_ID()) {
	   //if (cnfi_old->Face()->Global_ID() == face->Global_ID()) {
	     if( cnei_old ){
#ifdef CONTACT_DEBUG_NODE
	       if( PRINT_THIS_NODE )
	 	postream << "  Choosing CPP2 since old interaction is against same face\n";
#endif
	       Process_Method = 
	 	ContactNodeFaceInteraction::CLOSEST_POINT_PROJECTION_2;
	       break;
	     }
	   }
         }
         // find out which physical face this interactions is most opposed to
         // NOTE: we are assuming rotations are small each step
         Real* pf_normal = cnei_old->Get_Physical_Face_Normal();
#ifdef CONTACT_DEBUG_NODE
         if( PRINT_THIS_NODE ){
           postream << "    pf_normal = (" << pf_normal[0] << ","
                   << pf_normal[1] << "," << pf_normal[2] << ")\n";
         }
#endif
         Real most_aligned_old = normals[0][0]*pf_normal[0] +
                                 normals[0][1]*pf_normal[1] +
                                 normals[0][2]*pf_normal[2] ;
         int physical_face_old = 0;
         for( int i=1 ; i<num_physical_faces ; ++i ){
           Real dot = normals[i][0]*pf_normal[0] +
                      normals[i][1]*pf_normal[1] +
                      normals[i][2]*pf_normal[2] ;
           if( dot > most_aligned_old ){
             most_aligned_old = dot;
             physical_face_old = i;
           }
         }
#ifdef CONTACT_DEBUG_NODE
         if( PRINT_THIS_NODE )
           postream << "    Most Opposed: " << physical_face_old << "\n";
#endif
          // If the current physical face matches the old physical face,
          // we are done and can use CPPROJ
         if( physical_face_old == physical_face ){
#ifdef CONTACT_DEBUG_NODE
           if( PRINT_THIS_NODE )
              postream << "  Choosing CPP2\n";
#endif
           Process_Method = 
              ContactNodeFaceInteraction::CLOSEST_POINT_PROJECTION_2;
           break;
         }
       }
     }
   }
 }
 return Process_Method;
}

void
ContactSearch::Print_Search_Summary()
{
  if (contact_processor_number(SearchComm)==0) {
    do_node_node_search = 0;
    do_node_face_search = 0;
    do_face_face_search = 0;
    do_coverage_search  = 0;
    do_elem_elem_search = 0;
    if (search_data->Have_Node_Node_Interactions()      ) do_node_node_search = 1;
    if (search_data->Have_Node_Face_Interactions()      ) do_node_face_search = 1;
    if (search_data->Have_Face_Face_Interactions()      ) do_face_face_search = 1;
    if (search_data->Have_Face_Coverage_Interactions()  ) do_coverage_search  = 1;
    if (search_data->Have_Element_Element_Interactions()) do_elem_elem_search = 1;
    const char* offon[2]  = { "off", "on" };
    const char* sm_res[2] = { "NODE_BASED", "EDGE_BASED" };
    const char* pf_alg[3] = { "NONE", "FACE_BASED", "EDGE_BASED" };
    const char* cull[2]   = { "NONE", "SLAVE" };
    const char* track[3]  = { "NONE", "LOCAL", "GLOBAL" };
    std::cout << "ACME Options ----------------------------------------\n";
    
    std::cout << "  Use Secondary             = " << offon[!no_secondary] << "\n";
    std::cout << "  Search Cull               = " << cull[global_search_cull] << "\n";
    std::cout << "  Parallel Consistency      = " << offon[!no_parallel_consistency] << "\n";
    
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
    std::cout << "    Sharp/Smooth Curvature  = " << sharp_smooth_curvature << "\n";
    }
    std::cout << "  Normal smoothing          = " << offon[normal_smoothing_status] << "\n";
    if (normal_smoothing_status){
    std::cout << "  Sharp/Smooth Curvature    = " << sharp_smooth_curvature << "\n";
    std::cout << "  Normal Smoothing Distance = " << normal_smoothing_distance << "\n";
    std::cout << "  Smoothing Resolution      = " << sm_res[smoothing_resolution] << "\n";
    }
    std::cout << "  Partition Gap Status      = " << offon[partition_gap_status] << "\n";
    std::cout << "  Compute Node Areas        = " << offon[compute_node_areas] << "\n";
    std::cout << "  Auto Tolerance            = " << offon[auto_tol] << "\n";
    if (auto_tol) {
    std::cout << "    aggressive_tolerances   = " << offon[aggressive_tolerances] << "\n";
    }
    std::cout << "    box_inflation           = " << box_inflation << "\n";
    std::cout << "    gap_inflation           = " << gap_inflation << "\n";
    if (primary_topology->Have_Shells()) {
    std::cout << "  Simple shell lofting      = " << offon[compute_node_areas] << "\n";
    }
    }
    if (do_elem_elem_search) {
    std::cout << "  No warped volume          = " << offon[no_warped_volume] << "\n";
    }
    std::cout << flush;
  }
}

void
ContactSearch::Print_Search_Header(SearchType search_type, Real dt, Real dt_old, ContactTDEnforcement* td_enf)
{
  #if CONTACT_DEBUG_PRINT_LEVEL>=2
  if (step_number==0) {
    const char* offon[2]  = { "off", "on" };
    const char* sm_res[2] = { "NODE_BASED", "EDGE_BASED" };
    const char* pf_alg[3] = { "NONE", "FACE_BASED", "EDGE_BASED" };
    const char* cull[2]   = { "NONE", "SLAVE" };
    const char* track[3]  = { "NONE", "LOCAL", "GLOBAL" };
    postream << "ACME Options/Data ----------------------------------------\n";
    
    postream << "  No Secondary              = " << offon[no_secondary] << "\n";
    postream << "  Search Cull               = " << cull[global_search_cull] << "\n";
    postream << "  Keep Ghosting             = " << offon[keep_ghosting] << "\n";
    postream << "  No Parallel Consistency   = " << offon[no_parallel_consistency] << "\n";
    
    if (do_node_face_search) {
    postream << "  Tracking Type             = " << track[tracking_type] << "\n";
    if (tracking_type==GLOBAL_TRACKING) {
    postream << "    Tracking Interval       = " << global_tracking_interval << "\n";
    }
    if (tracking_type==LOCAL_TRACKING) {
    postream << "    Tracking Interval       = " << global_tracking_interval << "\n";
    postream << "    Off Face Tracking       = " << offon[enable_off_face_tracking] << "\n";
    }
    postream << "  Old Dynamic Search        = " << offon[old_dynamic_search] << "\n";
    postream << "  Physical Face Algorithm   = " << pf_alg[physical_face_algorithm]     << "\n";
    postream << "  Multiple interactions     = " << offon[multiple_interaction_status] << "\n";
    if (multiple_interaction_status){
    postream << "    Sharp/Smooth Curvature  = " << sharp_smooth_curvature << "\n";
    }
    postream << "  Normal smoothing          = " << offon[normal_smoothing_status] << "\n";
    if (normal_smoothing_status){
    postream << "  Sharp/Smooth Curvature    = " << sharp_smooth_curvature << "\n";
    postream << "  Normal Smoothing Distance = " << normal_smoothing_distance << "\n";
    postream << "  Smoothing Resolution      = " << sm_res[smoothing_resolution] << "\n";
    }
    postream << "  Partition Gap Status      = " << offon[partition_gap_status] << "\n";
    postream << "  Compute Node Areas        = " << offon[compute_node_areas] << "\n";
    postream << "  Auto Tolerance            = " << offon[auto_tol] << "\n";
    if (auto_tol) {
    postream << "    aggressive_tolerances   = " << offon[aggressive_tolerances] << "\n";
    }
    postream << "    box_inflation           = " << box_inflation << "\n";
    postream << "    gap_inflation           = " << gap_inflation << "\n";
    if (primary_topology->Have_Shells())
       postream << "  Simple shell lofting      = " << offon[compute_node_areas] << "\n";
    }
    if (do_elem_elem_search) {
      postream << "  No warped volume          = " << offon[no_warped_volume] << "\n";
    }
    search_data->Display(postream);
    if( td_enf ) td_enf->Display_Enforcement_Data();
  }
  postream<<"\nContactSearch::TopLevelSearch(), step = "<<step_number<<"\n";
  switch (search_type) {
  case STATIC1CONFIG:
    postream<<"  Static_Search_1_Configuration()\n";
    break;
  case STATIC1CONFIG_TIED:
    postream<<"  Static_Search_1_Configuration_Tied()\n";
    break;
  case STATIC2CONFIG:
    postream<<"  Static_Search_2_Configuration()\n";
    break;
  case DYNAMIC2CONFIG:
    postream<<"  Dynamic_Search_2_Configuration()\n";
    if( old_dynamic_search==ContactSearch::INACTIVE ) 
      postream << "    augmented search\n";
    else
      postream << "    non-augmented search\n";
    postream << "    dt_old = " << dt_old << "\n";
    postream << "    dt     = " << dt << "\n";
    break;
  }
  #else
  #ifdef CONTACT_HEARTBEAT 
  if (contact_processor_number(SearchComm)==0) {
    if (step_number==0) {
      const char* offon[2]  = { "off", "on" };
      const char* sm_res[2] = { "NODE_BASED", "EDGE_BASED" };
      const char* pf_alg[3] = { "NONE", "FACE_BASED", "EDGE_BASED" };
      const char* cull[2]   = { "NONE", "SLAVE" };
      const char* track[3]  = { "NONE", "LOCAL", "GLOBAL" };
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
      std::cout << "    Sharp/Smooth Curvature  = " << sharp_smooth_curvture << "\n";
      }
      std::cout << "  Normal smoothing          = " << offon[normal_smoothing_status] << "\n";
      if (normal_smoothing_status){
      std::cout << "  Sharp/Smooth Curvature    = " << sharp_smooth_curvature << "\n";
      std::cout << "  Normal Smoothing Distance = " << normal_smoothing_distance << "\n";
      std::cout << "  Smoothing Resolution      = " << sm_res[smoothing_resolution] << "\n";
      }
      std::cout << "  Partition Gap Status      = " << offon[partition_gap_status] << "\n";
      std::cout << "  Compute Node Areas        = " << offon[compute_node_areas] << "\n";
      std::cout << "  Auto Tolerance            = " << offon[auto_tol] << "\n";
      if (auto_tol) {
      std::cout << "    aggressive_tolerances   = " << offon[aggressive_tolerances] << "\n";
      }
      std::cout << "    box_inflation           = " << box_inflation << "\n";
      std::cout << "    gap_inflation           = " << gap_inflation << "\n";
      if (primary_topology->Have_Shells())
         std::cout << "  Simple shell lofting      = " << offon[compute_node_areas] << "\n";
      }
      if (do_elem_elem_search) {
        std::cout << "  No warped volume          = " << offon[no_warped_volume] << "\n";
      }
      search_data->Display0();
      if( td_enf ) td_enf->Display0_Enforcement_Data();
    }
    std::cout<<"\nContactSearch::TopLevelSearch(), step = "<<step_number<<"\n";
    switch (search_type) {
    case STATIC1CONFIG:
      std::cout<<"  Static_Search_1_Configuration()\n";
      break;
    case STATIC1CONFIG_TIED:
      std::cout<<"  Static_Search_1_Configuration_Tied()\n";
      break;
    case STATIC2CONFIG:
      std::cout<<"  Static_Search_2_Configuration()\n";
      break;
    case DYNAMIC2CONFIG:
      std::cout<<"  Dynamic_Search_2_Configuration()\n";
      if( old_dynamic_search==ContactSearch::INACTIVE ) 
        std::cout << "    augmented search\n";
      else
        std::cout << "    non-augmented search\n";
      std::cout << "    dt_old = " << dt_old << "\n";
      std::cout << "    dt     = " << dt << "\n";
      break;
    }
    std::cout << std::flush;
  }
  contact_global_sync(SearchComm);
  #else
  #ifdef CONTACT_ANALYZE_DATA_XFER
  postream<<"\nContactSearch::TopLevelSearch(), step = "<<step_number<<"\n";
  #endif
  #endif
  #endif

  #ifdef CONTACT_DEBUG_NODE
  primary_topology->Display_Debug_Node_IDs( postream );
  #endif
}

void
ContactSearch::Print_Search_Footer()
{
  #if CONTACT_DEBUG_PRINT_LEVEL>=2
  if (do_node_node_search) {
    primary_topology->Display_NodeNode_Interactions_Summary(postream, (char*)"  ");
  }
  if (do_node_face_search) {
    unsigned int all_search = ContactTopologyEntity<Real>::TRACK_SEARCH_SLAVE | ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE;
    primary_topology->Display_NodeEntity_Interactions_Summary(postream, all_search, (char*)"  ");
  }
  if (do_face_face_search) {
    primary_topology->Display_FaceFace_Interactions_Summary(postream, (char*)"  ");
  }
  if (do_coverage_search) {
    primary_topology->Display_FaceCoverage_Interactions_Summary(postream, (char*)"  ");
  }
  if (do_elem_elem_search) {
    primary_topology->Display_ElementElement_Interactions_Summary(postream, (char*)"  ");
  }
  #if CONTACT_DEBUG_PRINT_LEVEL>=4
  if (do_node_node_search) {
    postream << "\n\n\n***** Node-Node Interactions *****\n";
    primary_topology->Display_NodeNode_Interactions(postream);
  }
  if (do_node_face_search) {
    if( contact_processor_number(SearchComm) == 0 )
      std::cout << "\n\n\n***** Node-Entity Interactions *****" << std::endl;
    primary_topology->Display_NodeEntity_Interactions(postream);
  }
  if (do_face_face_search) {
    postream << "\n\n\n***** Face-Face Interactions *****\n";
    primary_topology->Display_FaceFace_Interactions(postream);
  }
  if (do_coverage_search) {
    postream << "\n\n\n***** Face-Coverage Interactions *****\n";
    primary_topology->Display_FaceCoverage_Interactions(postream);
  }
  if (do_elem_elem_search) {
    postream << "\n\n\n***** ElementElement Interactions *****\n";
    primary_topology->Display_ElementElement_Interactions(postream);
  }
  #endif
  #endif

  #ifdef CONTACT_DEBUG_NODE
  primary_topology->Display_Debug_Nodes( postream );
  #endif
  
  #if CONTACT_DEBUG_PRINT_LEVEL>=2 || defined(CONTACT_DEBUG_NODE) || defined(CONTACT_ANALYZE_DATA_XFER)
  #if CONTACT_DEBUG_PRINT_LEVEL>=2
  postream << "End of search\n";
  #endif
  postream.flush();
  #else
  #ifdef CONTACT_HEARTBEAT 
  if (contact_processor_number(SearchComm)==0) {
    std::cout<<"End of search\n";
    std::cout<<flush;
  }
  #endif
  #endif
}

void
ContactSearch::Initialize_Search()
{
  initialized         = true;
  do_node_node_search = 0;
  do_node_face_search = 0;
  do_face_face_search = 0;
  do_coverage_search  = 0;
  do_elem_elem_search = 0;
  if (search_data->Have_Node_Node_Interactions()      ) do_node_node_search = 1;
  if (search_data->Have_Node_Face_Interactions()      ) do_node_face_search = 1;
  if (search_data->Have_Face_Face_Interactions()      ) do_face_face_search = 1;
  if (search_data->Have_Face_Coverage_Interactions()  ) do_coverage_search  = 1;
  if (search_data->Have_Element_Element_Interactions()) do_elem_elem_search = 1;
  SetInitialContext();
  if (global_tracking_interval==0) {
    keep_ghosting   = INACTIVE;
    enable_tracking = INACTIVE;
    tracking_type   = NO_TRACKING;
  }
  if (tracking_type==GLOBAL_TRACKING) {
    if (global_tracking_interval>1) {
      keep_ghosting   = ACTIVE;
      no_secondary    = ACTIVE;
    } else {
      keep_ghosting   = INACTIVE;
      enable_tracking = INACTIVE;
      tracking_type   = NO_TRACKING;
    }
  } else {
    keep_ghosting = INACTIVE;
  }
  if (search_data->Have_Only_Tied_Interactions()) {
    keep_ghosting   = INACTIVE;
    enable_tracking = INACTIVE;
    tracking_type   = NO_TRACKING;
  }
  primary_topology->Compute_Characteristic_Length();
}

void
ContactSearch::Initialize_Context()
{ 
  if (!initialized_context) {
    initialized_context   = true;
    int num_primary_nodes = primary_topology->Number_of_Nodes();
    ContactNode<Real>** primary_nodes = 
      reinterpret_cast<ContactNode<Real>**>(primary_topology->NodeList()->EntityList());
    for (int i=0; i<num_primary_nodes; ++i) {
      primary_nodes[i]->ClearContext();
    }
        
    int num_primary_faces = primary_topology->Number_of_Faces();
    ContactFace<Real>** primary_faces = 
      reinterpret_cast<ContactFace<Real>**>(primary_topology->FaceList()->EntityList());
    for (int i=0; i<num_primary_faces; ++i) {
      primary_faces[i]->ClearContext();
    }
    
    int num_primary_elems = primary_topology->Number_of_Elements();
    ContactElement** primary_elems = 
      reinterpret_cast<ContactElement**>(primary_topology->ElemList()->EntityList());
    for (int i=0; i<num_primary_elems; ++i) {
      primary_elems[i]->ClearContext();
    }
  }
}

//
// Tag all nodes, faces, and elements as being potential 
// masters or slaves based on the search data matrix.
// It is meant to only be called once as an initialization.
//
void
ContactSearch::SetInitialContext()
{ 
  int num_primary_nodes = primary_topology->Number_of_Nodes();
  ContactNode<Real>** primary_nodes = 
    reinterpret_cast<ContactNode<Real>**>(primary_topology->NodeList()->EntityList());
  for (int i=0; i<num_primary_nodes; ++i) {
    primary_nodes[i]->ClearNonTiedContext();
  }
      
  int num_primary_faces = primary_topology->Number_of_Faces();
  ContactFace<Real>** primary_faces = 
    reinterpret_cast<ContactFace<Real>**>(primary_topology->FaceList()->EntityList());
  for (int i=0; i<num_primary_faces; ++i) {
    primary_faces[i]->ClearNonTiedContext();
  }
  
  int num_primary_elems = primary_topology->Number_of_Elements();
  ContactElement** primary_elems = 
    reinterpret_cast<ContactElement**>(primary_topology->ElemList()->EntityList());
  for (int i=0; i<num_primary_elems; ++i) {
    primary_elems[i]->ClearNonTiedContext();
  }
  
  if (do_node_node_search) {
    int first_block = 0;
    if (do_node_face_search) first_block = 1;
    ContactTopologyEntityList* node_list = primary_topology->NodeList();
    for (int i=first_block; i<primary_topology->Number_of_Node_Blocks(); ++i) {
      if (search_data->Is_NodeBlock_NodeNodeSlave(i)) {
        int nnodes = node_list->BlockNumEntities(i);
        ContactNode<Real>** Nodes = reinterpret_cast<ContactNode<Real>**>(node_list->BlockEntityList(i));
        for (int j=0; j<nnodes; ++j){
          Nodes[j]->SetSlave();
        } 
        primary_topology->Node_Block(i)->SetSlave();
      }
      if (search_data->Is_NodeBlock_NodeNodeMaster(i)) {
        int nnodes = node_list->BlockNumEntities(i);
        ContactNode<Real>** Nodes = reinterpret_cast<ContactNode<Real>**>(node_list->BlockEntityList(i));
        for (int j=0; j<nnodes; ++j){
          Nodes[j]->SetMaster();
        } 
        primary_topology->Node_Block(i)->SetMaster();
      }
    }
  }
    
  if (do_node_face_search) {
    ContactTopologyEntityList* face_list = primary_topology->FaceList();
    for (int i=0; i<primary_topology->Number_of_Face_Blocks(); ++i) {
      if (search_data->Is_FaceBlock_NodeFaceSlave(i)) {
        int nfaces = face_list->BlockNumEntities(i);
        ContactFace<Real>** Faces = 
          reinterpret_cast<ContactFace<Real>**>(face_list->BlockEntityList(i));
        for (int j=0; j<nfaces; ++j){
          ContactFace<Real>* face = Faces[j];
          int num_nodes = face->Nodes_Per_Face();
          for (int k=0; k<num_nodes; ++k) {
            ContactNode<Real>* node = face->Node(k);
            //if (node->Ownership()     != ContactTopologyEntity<Real>::OWNED) continue;
            if (node->Physical_Type() == ContactNode<Real>::SHELL_TAB_NODE) continue;
            if (node->Physical_Type() == ContactNode<Real>::MIXED_TAB_NODE) continue;
            node->SetSlave();
          }
        }
      }
      if (search_data->Is_FaceBlock_NodeFaceMaster(i)) {
        int nfaces = face_list->BlockNumEntities(i);
        ContactFace<Real>** Faces = 
          reinterpret_cast<ContactFace<Real>**>(face_list->BlockEntityList(i));
        for (int j=0; j<nfaces; ++j){
          Faces[j]->SetMaster();
        }
        primary_topology->Face_Block(i)->SetMaster();
      }
    }
    ContactTopologyEntityList* node_list = primary_topology->NodeList();
    for (int i=1; i<primary_topology->Number_of_Node_Blocks(); ++i) {
      if (search_data->Is_NodeBlock_NodeFaceSlave(i)) {
        int nnodes = node_list->BlockNumEntities(i);
        ContactNode<Real>** Nodes = reinterpret_cast<ContactNode<Real>**>(node_list->BlockEntityList(i));
        for (int j=0; j<nnodes; ++j){
          Nodes[j]->SetSlave();
        }
        primary_topology->Node_Block(i)->SetSlave();
      }
    }
  }
  
  if (do_face_face_search || do_coverage_search) {
    ContactTopologyEntityList* face_list = primary_topology->FaceList();
    for (int i=0; i<primary_topology->Number_of_Face_Blocks(); ++i) {
      if (search_data->Is_FaceBlock_FaceFaceSlave(i)) {
        int nfaces = face_list->BlockNumEntities(i);
        ContactFace<Real>** Faces = 
          reinterpret_cast<ContactFace<Real>**>(face_list->BlockEntityList(i));
        for (int j=0; j<nfaces; ++j){
          Faces[j]->SetSlave();
        }
        primary_topology->Face_Block(i)->SetSlave();
      }
      if (search_data->Is_FaceBlock_FaceFaceMaster(i)) {
        int nfaces = face_list->BlockNumEntities(i);
        ContactFace<Real>** Faces = 
          reinterpret_cast<ContactFace<Real>**>(face_list->BlockEntityList(i));
        for (int j=0; j<nfaces; ++j){
          Faces[j]->SetMaster();
        }
        primary_topology->Face_Block(i)->SetMaster();
      }
    }
  }
  
  if (do_elem_elem_search) {
    ContactTopologyEntityList* elem_list = primary_topology->ElemList();
    for (int i=0; i<primary_topology->Number_of_Element_Blocks(); ++i) {
      if (search_data->Is_ElemBlock_ElemElemSlave(i)) {
        int nelems = elem_list->BlockNumEntities(i);
        ContactElement** Elems = 
          reinterpret_cast<ContactElement**>(elem_list->BlockEntityList(i));
        for (int j=0; j<nelems; ++j){
          Elems[j]->SetSlave();
        }
        primary_topology->Element_Block(i)->SetSlave();
      }
      if (search_data->Is_ElemBlock_ElemElemMaster(i)) {
        int nelems = elem_list->BlockNumEntities(i);
        ContactElement** Elems = 
          reinterpret_cast<ContactElement**>(elem_list->BlockEntityList(i));
        for (int j=0; j<nelems; ++j){
          Elems[j]->SetMaster();
        }
        primary_topology->Element_Block(i)->SetMaster();
      }
    }
  }
}

void ContactSearch::ResetContextForSearchSlaves() {
  int num_primary_nodes = primary_topology->Number_of_Nodes();
  ContactNode<Real>** primary_nodes = 
    reinterpret_cast<ContactNode<Real>**>(primary_topology->NodeList()->EntityList());
  for (int i=0; i<num_primary_nodes; ++i) {
    primary_nodes[i]->ClearContextBit(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE);
  }
      
  int num_primary_faces = primary_topology->Number_of_Faces();
  ContactFace<Real>** primary_faces = 
    reinterpret_cast<ContactFace<Real>**>(primary_topology->FaceList()->EntityList());
  for (int i=0; i<num_primary_faces; ++i) {
    primary_faces[i]->ClearContextBit(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE);
  }
  
  int num_primary_elems = primary_topology->Number_of_Elements();
  ContactElement** primary_elems = 
    reinterpret_cast<ContactElement**>(primary_topology->ElemList()->EntityList());
  for (int i=0; i<num_primary_elems; ++i) {
    primary_elems[i]->ClearContextBit(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE);
  }
}

//
// Tag all nodes, faces, and elements as being active in the  
// slaves in global search based either: (1) just the search 
// data matrix or; or (2) the search data matrix and entity
// block overlap.
//
void
ContactSearch::SetContextForSearchSlaves()
{
  if (step_number==0 || (tracking_type==LOCAL_TRACKING  && tracking_step <= 1) || restart) {
    bool context_set = false;

    int num_primary_nodes = primary_topology->Number_of_Nodes();
    int num_primary_faces = primary_topology->Number_of_Faces();  
    int num_primary_elems = primary_topology->Number_of_Elements();

    ContactNode<Real>** primary_nodes = 
      reinterpret_cast<ContactNode<Real>**>(primary_topology->NodeList()->EntityList());
    ContactFace<Real>** primary_faces = 
      reinterpret_cast<ContactFace<Real>**>(primary_topology->FaceList()->EntityList());
    ContactElement** primary_elems = 
      reinterpret_cast<ContactElement**>(primary_topology->ElemList()->EntityList());

    if (global_search_cull==SLAVE_CULL) {
      context_set = true;
      // if culling by slave specification in the search data, then the status
      // won't change from step-to-step so just set the contact once at step zero
      ResetContextForSearchSlaves();
      //
      // Master/Slave Culling
      //
      if (do_face_face_search || do_coverage_search) {
        ContactTopologyEntityList* face_list = primary_topology->FaceList();
        for (int i=0; i<primary_topology->Number_of_Face_Blocks(); ++i) {
          if (!search_data->Is_FaceBlock_FaceFaceSlave(i)) continue;
          int nfaces = face_list->BlockNumEntities(i);
          ContactFace<Real>** Faces = 
            reinterpret_cast<ContactFace<Real>**>(face_list->BlockEntityList(i));
          for (int j=0; j<nfaces; ++j){
            Faces[j]->SetContextBit(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE);
          }
        }
      } //if (do_face_face_search || do_coverage_search)
    
      if (do_elem_elem_search) {
        ContactTopologyEntityList* elem_list = primary_topology->ElemList();
        for (int i=0; i<primary_topology->Number_of_Element_Blocks(); ++i) {
          if (!search_data->Is_ElemBlock_ElemElemSlave(i)) continue;
          int nelems = elem_list->BlockNumEntities(i);
          ContactElement** Elems = 
            reinterpret_cast<ContactElement**>(elem_list->BlockEntityList(i));
          for (int j=0; j<nelems; ++j){
            Elems[j]->SetContextBit(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE);
          }
        }
      } //if (do_elem_elem_search)
      
      if (do_node_node_search) {
        ContactTopologyEntityList* node_list = primary_topology->NodeList();
        for (int i=0; i<primary_topology->Number_of_Node_Blocks(); ++i) {
          if (primary_topology->Node_Block(i)->Type() != POINT) continue;
          if (search_data->Is_NodeBlock_NodeNodeSlave(i)) {
            int nnodes = node_list->BlockNumEntities(i);
            ContactNode<Real>** Nodes = reinterpret_cast<ContactNode<Real>**>(node_list->BlockEntityList(i));
            for (int j=0; j<nnodes; ++j){
              ContactNode<Real>* node = Nodes[j];
              if (node->Ownership()     != ContactTopologyEntity<Real>::OWNED) continue;
              if (node->Physical_Type() == ContactNode<Real>::SHELL_TAB_NODE)  continue;
              if (node->Physical_Type() == ContactNode<Real>::MIXED_TAB_NODE)  continue;
              node->SetContextBit(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE);
            } 
          }
        }
      } //if (do_node_node_search)
      
      if (do_node_face_search) {
        ContactTopologyEntityList* node_list = primary_topology->NodeList();
        ContactTopologyEntityList* face_list = primary_topology->FaceList();
        for (int i=0; i<primary_topology->Number_of_Face_Blocks(); ++i) {
          if (!search_data->Is_FaceBlock_NodeFaceSlave(i)) continue;
          int nfaces = face_list->BlockNumEntities(i);
          ContactFace<Real>** Faces = 
            reinterpret_cast<ContactFace<Real>**>(face_list->BlockEntityList(i));
          for (int j=0; j<nfaces; ++j){
            ContactFace<Real>* face = Faces[j];
            int num_nodes = face->Nodes_Per_Face();
            for (int k=0; k<num_nodes; ++k) {
              ContactNode<Real>* node = face->Node(k);
              if (node->Physical_Type() == ContactNode<Real>::SHELL_TAB_NODE)  continue;
              if (node->Physical_Type() == ContactNode<Real>::MIXED_TAB_NODE)  continue;
              node->SetContextBit(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE);
            }
          }
        }
        for (int i=1; i<primary_topology->Number_of_Node_Blocks(); ++i) {
          if (!search_data->Is_NodeBlock_NodeFaceSlave(i)) continue;
          int nnodes = node_list->BlockNumEntities(i);
          ContactNode<Real>** Nodes = reinterpret_cast<ContactNode<Real>**>(node_list->BlockEntityList(i));
          for (int j=0; j<nnodes; ++j){
            ContactNode<Real>* node = Nodes[j];
            if (node->Physical_Type() == ContactNode<Real>::SHELL_TAB_NODE)  continue;
            if (node->Physical_Type() == ContactNode<Real>::MIXED_TAB_NODE)  continue;
            node->SetContextBit(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE);
          }
        }
      } //if (do_node_face_search)
    } else {
      context_set = true;
      // if no pre-culling is done then the status won't change from 
      // step-to-step so just set the contact once at step zero
      //
      //  No Culling Set all search contexts to true, no need 
      //  to reset the contexts as each are explicity set
      //
      for (int i=0; i<num_primary_nodes; ++i) {
        ContactNode<Real>* node = primary_nodes[i];
        if (node->Physical_Type() == ContactNode<Real>::SHELL_TAB_NODE)  continue;
        if (node->Physical_Type() == ContactNode<Real>::MIXED_TAB_NODE)  continue;
        primary_nodes[i]->SetContextBit(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE);
      }
        
      for (int i=0; i<num_primary_faces; ++i) {
        primary_faces[i]->SetContextBit(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE);
      }
    
      for (int i=0; i<num_primary_elems; ++i) {
        primary_elems[i]->SetContextBit(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE);
      }
    }

    if (tracking_type==LOCAL_TRACKING && tracking_step!=0) {
      for (int i=0; i<num_primary_nodes; ++i) {
        primary_nodes[i]->ClearContextBit(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE);
      }
      ContactTopologyEntityList* face_list = primary_topology->FaceList();
      for (int i=0; i<primary_topology->Number_of_Face_Blocks(); ++i) {
        if (!search_data->Is_FaceBlock_NodeFaceSlave(i)) continue;
        int tracking  = 0x7ffffff;
        int slave_key = primary_topology->Face_Block(i)->Entity_Key();
        int nfaces    = face_list->BlockNumEntities(i);
        ContactFace<Real>** BlockFaces = 
          reinterpret_cast<ContactFace<Real>**>(face_list->BlockEntityList(i));
        for (int j=0; j<primary_topology->Number_of_Face_Blocks(); ++j) {
          int master_key = primary_topology->Face_Block(j)->Entity_Key();
          int interaction_type = (int)search_data->Get_Search_Data(INTERACTION_TYPE, 
                                                                   slave_key,
                                                                   master_key);
          if (interaction_type == NO_INTERACTION) continue;
          tracking = 1;
        } 
        for (int j=0; j<primary_topology->Number_of_Analytic_Surfaces(); ++j) {
          int master_key = primary_topology->Analytic_Surface(j)->Entity_Key();
          int interaction_type = (int)search_data->Get_Search_Data(INTERACTION_TYPE, 
                                                                   slave_key,
                                                                   master_key);
          if (interaction_type == NO_INTERACTION) continue;
          tracking = 1; 
        }
        if (tracking>0) {
          for (int jj=0; jj<nfaces; ++jj){
            ContactFace<Real>* face = BlockFaces[jj];
            int num_nodes = face->Nodes_Per_Face();
            for (int k=0; k<num_nodes; ++k) {
              face->Node(k)->SetContextBit(ContactTopologyEntity<Real>::TRACK_SEARCH_SLAVE);
            }
          }
        }
      }
    }
    
    // have to do a 'swapadd' here because some parallel decompositions can 
    // lead to a node being owned by a processor but it's parent face resides 
    // on another processor
    #ifndef CONTACT_NO_MPI
    if (context_set) {
      contact_swapadd_context(SearchComm, *(primary_topology->Node_Sym_Comm()), *comm_buffer);
    }
    #endif
    
    if (tracking_type==LOCAL_TRACKING && tracking_step!=0 ) {
      for (int i=0; i<num_primary_nodes; ++i) {
        ContactNode<Real>* node = primary_nodes[i];
        if (node->Ownership()!=ContactTopologyEntity<Real>::OWNED) {
          node->ClearContextBit(ContactTopologyEntity<Real>::TRACK_SEARCH_SLAVE );
          node->ClearContextBit(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE);
        } else {
          if (node->CheckContext(ContactTopologyEntity<Real>::TRACK_SEARCH_SLAVE )) {
            node->ClearContextBit(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE);
            if (node->Number_NodeEntity_Interactions(STATE1)==0) {
              node->ClearContextBit(ContactTopologyEntity<Real>::TRACK_SEARCH_SLAVE);
            }
          }
        }
      }
      
      int n_track_nodes  = 0;
      int n_global_nodes = 0;
      num_tracked_interactions = 0;
      for (int i=0; i<num_primary_nodes; ++i) {
        if (primary_nodes[i]->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)) ++n_global_nodes;
        if (primary_nodes[i]->CheckContext(ContactTopologyEntity<Real>::TRACK_SEARCH_SLAVE )) {
          num_tracked_interactions += primary_nodes[i]->Number_NodeEntity_Interactions(STATE1);
          ++n_track_nodes;
        }
      }
          
      int n_track_faces  = 0;
      int n_global_faces = 0;
      for (int i=0; i<num_primary_faces; ++i) {
        if (primary_faces[i]->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)) ++n_global_faces;
        if (primary_faces[i]->CheckContext(ContactTopologyEntity<Real>::TRACK_SEARCH_SLAVE )) ++n_track_faces;
      }
      
      int n_track_elems  = 0;
      int n_global_elems = 0;
      for (int i=0; i<num_primary_elems; ++i) {
        if (primary_elems[i]->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)) ++n_global_elems;
        if (primary_elems[i]->CheckContext(ContactTopologyEntity<Real>::TRACK_SEARCH_SLAVE )) ++n_track_elems;
      }
      
      int num_track_entities  = n_track_nodes+n_track_faces+n_track_elems;
      int num_global_entities = n_global_nodes+n_global_faces+n_global_elems;
      
      num_global_nodes  = num_global_entities;
      num_tracked_nodes = num_track_entities;
      
      #if CONTACT_DEBUG_PRINT_LEVEL>=2
      postream<<"    Tagged "<<n_global_nodes<<" nodes as slaves for the global search\n";
      postream<<"    Tagged "<<n_global_faces<<" faces as slaves for the global search\n";
      postream<<"    Tagged "<<n_global_elems<<" elems as slaves for the global search\n";
      postream<<"    Tagged "<<n_track_nodes<<" nodes as slaves for the track search\n";
      postream<<"    Tagged "<<n_track_faces<<" faces as slaves for the track search\n";
      postream<<"    Tagged "<<n_track_elems<<" elems as slaves for the track search\n";
      #else
      #ifdef CONTACT_HEARTBEAT
      int number_nodes1 = contact_global_sum(n_global_nodes, SearchComm );
      int number_faces1 = contact_global_sum(n_global_faces, SearchComm );
      int number_elems1 = contact_global_sum(n_global_elems, SearchComm );
      int number_nodes2 = contact_global_sum(n_track_nodes, SearchComm );
      int number_faces2 = contact_global_sum(n_track_faces, SearchComm );
      int number_elems2 = contact_global_sum(n_track_elems, SearchComm );
      if (contact_processor_number(SearchComm)==0) {
        std::cout<<"    Tagged "<<number_nodes1<<" nodes as slaves for the global search\n";
        std::cout<<"    Tagged "<<number_faces1<<" faces as slaves for the global search\n";
        std::cout<<"    Tagged "<<number_elems1<<" elems as slaves for the global search\n";
        std::cout<<"    Tagged "<<number_nodes2<<" nodes as slaves for the track search\n";
        std::cout<<"    Tagged "<<number_faces2<<" faces as slaves for the track search\n";
        std::cout<<"    Tagged "<<number_elems2<<" elems as slaves for the track search\n";
      }
      #endif
      #endif
    } else if (context_set) {
      for (int i=0; i<num_primary_nodes; ++i) {
        ContactNode<Real>* node = primary_nodes[i];
        if (node->Ownership()!=ContactTopologyEntity<Real>::OWNED) {
          node->ClearContextBit(ContactTopologyEntity<Real>::TRACK_SEARCH_SLAVE );
          node->ClearContextBit(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE);
        }
      }
      #if CONTACT_DEBUG_PRINT_LEVEL>=2 || defined(CONTACT_HEARTBEAT)
      int n_global_nodes = 0;
      for (int i=0; i<num_primary_nodes; ++i) {
        if (primary_nodes[i]->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)) ++n_global_nodes;
      }
          
      int n_global_faces = 0;
      for (int i=0; i<num_primary_faces; ++i) {
        if (primary_faces[i]->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)) ++n_global_faces;
      }
      
      int n_global_elems = 0;
      for (int i=0; i<num_primary_elems; ++i) {
        if (primary_elems[i]->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)) ++n_global_elems;
      }
      
      num_tracked_nodes = 0;
      
      #if CONTACT_DEBUG_PRINT_LEVEL>=2
      postream<<"    Tagged "<<n_global_nodes<<" nodes as slaves for the global search\n";
      postream<<"    Tagged "<<n_global_faces<<" faces as slaves for the global search\n";
      postream<<"    Tagged "<<n_global_elems<<" elems as slaves for the global search\n";
      #else
      #ifdef CONTACT_HEARTBEAT
      int number_nodes1 = contact_global_sum(n_global_nodes, SearchComm );
      int number_faces1 = contact_global_sum(n_global_faces, SearchComm );
      int number_elems1 = contact_global_sum(n_global_elems, SearchComm );
      if (contact_processor_number(SearchComm)==0) {
        std::cout<<"    Tagged "<<number_nodes1<<" nodes as slaves for the global search\n";
        std::cout<<"    Tagged "<<number_faces1<<" faces as slaves for the global search\n";
        std::cout<<"    Tagged "<<number_elems1<<" elems as slaves for the global search\n";
      }
      #endif
      #endif
      #endif
    }
  }
}

//
// Tag appropriate nodes, faces, and elements for the geometry update.
//
void
ContactSearch::SetContextForGeometryUpdate()
{
  #if CONTACT_DEBUG_PRINT_LEVEL>=2 || defined(CONTACT_HEARTBEAT)
  bool did_update = false;
  #endif
  if (step_number==0 || restart || (tracking_type!=NO_TRACKING && tracking_step<=1)) {
    #if CONTACT_DEBUG_PRINT_LEVEL>=2 || defined(CONTACT_HEARTBEAT)
    did_update = true;
    #endif
    // Set the geometry update context for all nodes
    int num_primary_nodes = primary_topology->Number_of_Nodes();
    ContactNode<Real>** primary_nodes = 
      reinterpret_cast<ContactNode<Real>**>(primary_topology->NodeList()->EntityList());
    for (int i=0; i<num_primary_nodes; ++i) {
      primary_nodes[i]->SetContextBit(ContactTopologyEntity<Real>::GEOMETRY_UPDATE);
    } 
    // Set the geometry update context for all faces
    int num_primary_faces = primary_topology->Number_of_Faces();
    ContactFace<Real>** primary_faces = 
      reinterpret_cast<ContactFace<Real>**>(primary_topology->FaceList()->EntityList());
    for (int i=0; i<num_primary_faces; ++i) {
      primary_faces[i]->SetContextBit(ContactTopologyEntity<Real>::GEOMETRY_UPDATE);
    }
    // Set the geometry update context for all elements
    int num_primary_elems = primary_topology->Number_of_Elements();
    ContactElement** primary_elems = 
      reinterpret_cast<ContactElement**>(primary_topology->ElemList()->EntityList());
    for (int i=0; i<num_primary_elems; ++i) {
      primary_elems[i]->SetContextBit(ContactTopologyEntity<Real>::GEOMETRY_UPDATE);
    }
  } else if (tracking_type==GLOBAL_TRACKING && tracking_step==1) {
    #if CONTACT_DEBUG_PRINT_LEVEL>=2 || defined(CONTACT_HEARTBEAT)
    did_update = true;
    #endif
    // Clear the geometry update context for all nodes
    int num_primary_nodes = primary_topology->Number_of_Nodes();
    ContactNode<Real>** primary_nodes = 
      reinterpret_cast<ContactNode<Real>**>(primary_topology->NodeList()->EntityList());
    for (int i=0; i<num_primary_nodes; ++i) {
      primary_nodes[i]->ClearContextBit(ContactTopologyEntity<Real>::GEOMETRY_UPDATE);
    }
        
    // Clear the geometry update context for all faces
    int num_primary_faces = primary_topology->Number_of_Faces();
    ContactFace<Real>** primary_faces = 
      reinterpret_cast<ContactFace<Real>**>(primary_topology->FaceList()->EntityList());
    for (int i=0; i<num_primary_faces; ++i) {
      primary_faces[i]->ClearContextBit(ContactTopologyEntity<Real>::GEOMETRY_UPDATE);
    }
    
    // Clear the geometry update context for all elements
    int num_primary_elems = primary_topology->Number_of_Elements();
    ContactElement** primary_elems = 
      reinterpret_cast<ContactElement**>(primary_topology->ElemList()->EntityList());
    for (int i=0; i<num_primary_elems; ++i) {
      primary_elems[i]->ClearContextBit(ContactTopologyEntity<Real>::GEOMETRY_UPDATE);
    }
  
    // tag all nodes that are in proximity and all faces that are connected to them
    for(int i=0; i<primary_topology->Number_of_Node_Blocks(); ++i) {
      if (!(primary_topology->Node_Block(i)->Has_Normal_Attributes())) {
        int nnodes = primary_topology->NodeList()->BlockNumEntities(i);
        ContactNode<Real>** nodes  = 
          reinterpret_cast<ContactNode<Real>**>(primary_topology->NodeList()->BlockEntityList(i));
        for (int j=0; j<nnodes; ++j) {
          ContactNode<Real>* node = nodes[j];
          if (node->in_proximity || node->CheckContext(ContactTopologyEntity<Real>::TIED)) {
            if (!node->CheckContext(ContactTopologyEntity<Real>::GHOSTED_FOR_SEARCH)) {
              node->SetContextBit(ContactTopologyEntity<Real>::GEOMETRY_UPDATE);
              int num_faces = node->Number_Face_Connections();
              for (int k=0; k<num_faces; ++k) {
                node->GetFace(k)->SetContextBit(ContactTopologyEntity<Real>::GEOMETRY_UPDATE);
              }
            }
          }
        }
      }
    }
    
    // tag all faces that are in proximity and nodes that are connected to tagged faces
    for (int i=0; i<num_primary_faces; ++i) {
      ContactFace<Real>* face = primary_faces[i];
      if (face->in_proximity || face->CheckContext(ContactTopologyEntity<Real>::TIED)) {
        if (!face->CheckContext(ContactTopologyEntity<Real>::GHOSTED_FOR_SEARCH)) {
          face->SetContextBit(ContactTopologyEntity<Real>::GEOMETRY_UPDATE);
        }
      }
      if (face->CheckContext(ContactTopologyEntity<Real>::GEOMETRY_UPDATE)) {
        for(int j=0 ; j<face->Nodes_Per_Face(); ++j ){
          if (!face->Node(j)->CheckContext(ContactTopologyEntity<Real>::GHOSTED_FOR_SEARCH)) {
            face->Node(j)->SetContextBit(ContactTopologyEntity<Real>::GEOMETRY_UPDATE);
          }
        }
      }
    }
    
    // Set the geometry update context for all elements
    for (int i=0; i<num_primary_elems; ++i) {
      primary_elems[i]->SetContextBit(ContactTopologyEntity<Real>::GEOMETRY_UPDATE);
    }
  }
  
  #if CONTACT_DEBUG_PRINT_LEVEL>=2 || defined(CONTACT_HEARTBEAT)
  if (did_update) {
    int num_tagged_nodes  = 0;
    int num_primary_nodes = primary_topology->Number_of_Nodes();
    ContactNode<Real>** primary_nodes = 
      reinterpret_cast<ContactNode<Real>**>(primary_topology->NodeList()->EntityList());
    for (int i=0; i<num_primary_nodes; ++i) {
      if (primary_nodes[i]->CheckContext(ContactTopologyEntity<Real>::GEOMETRY_UPDATE)) ++num_tagged_nodes;
    }
        
    int num_tagged_faces  = 0;
    int num_primary_faces = primary_topology->Number_of_Faces();
    ContactFace<Real>** primary_faces = 
      reinterpret_cast<ContactFace<Real>**>(primary_topology->FaceList()->EntityList());
    for (int i=0; i<num_primary_faces; ++i) {
      if (primary_faces[i]->CheckContext(ContactTopologyEntity<Real>::GEOMETRY_UPDATE)) ++num_tagged_faces;
    }
    
    int num_tagged_elems  = 0;
    int num_primary_elems = primary_topology->Number_of_Elements();
    ContactElement** primary_elems = 
      reinterpret_cast<ContactElement**>(primary_topology->ElemList()->EntityList());
    for (int i=0; i<num_primary_elems; ++i) {
      if (primary_elems[i]->CheckContext(ContactTopologyEntity<Real>::GEOMETRY_UPDATE)) ++num_tagged_elems;
    }
    #if CONTACT_DEBUG_PRINT_LEVEL>=2
    postream<<"    Tagged "<<num_tagged_nodes<<" of "<<num_primary_nodes<<" nodes for geometry update\n";
    postream<<"    Tagged "<<num_tagged_faces<<" of "<<num_primary_faces<<" faces for geometry update\n";
    postream<<"    Tagged "<<num_tagged_elems<<" of "<<num_primary_elems<<" elems for geometry update\n";
    #else
    #ifdef CONTACT_HEARTBEAT
    int num_nodes_tagged = contact_global_sum(num_tagged_nodes, SearchComm );
    int num_faces_tagged = contact_global_sum(num_tagged_faces, SearchComm );
    int num_elems_tagged = contact_global_sum(num_tagged_elems, SearchComm );
    int num_nodes_total  = contact_global_sum(num_primary_nodes, SearchComm );
    int num_faces_total  = contact_global_sum(num_primary_faces, SearchComm );
    int num_elems_total  = contact_global_sum(num_primary_elems, SearchComm );
    if (contact_processor_number(SearchComm)==0) {
      std::cout<<"    Tagged "<<num_nodes_tagged<<" of "<<num_nodes_total<<" nodes for geometry update\n";
      std::cout<<"    Tagged "<<num_faces_tagged<<" of "<<num_faces_total<<" faces for geometry update\n";
      std::cout<<"    Tagged "<<num_elems_tagged<<" of "<<num_elems_total<<" elems for geometry update\n";
    }
    #endif
    #endif
  }
  #endif
}

//
// Tag all nodes that will be used to determine the ghosting.
// These are all the slave nodes and the nodes whose parent 
// face/element are also slaves.
//
void
ContactSearch::SetContextForGhostingRCB()
{ 
  // if not a parallel run then just return
  if( contact_number_of_processors(SearchComm)==1 ) return;
  
  // if doing any tracking at the time then just return
  if (tracking_type!=NO_TRACKING && tracking_step!=0 ) return;
  
  // if not doing any tracking then only need to do this once
  if (tracking_type==NO_TRACKING && step_number!=0 ) return;
  
  if (tracking_type==GLOBAL_TRACKING) {
    // need to include tied entities since the ghosting is persistent
    int num_primary_nodes = primary_topology->Number_of_Nodes();
    ContactNode<Real>** primary_nodes = 
      reinterpret_cast<ContactNode<Real>**>(primary_topology->NodeList()->EntityList());
    for (int i=0; i<num_primary_nodes; ++i) {
      ContactNode<Real>* node = primary_nodes[i];
      node->ClearContextBit(ContactTopologyEntity<Real>::ACTIVE_FOR_GHOSTING_RCB);
      if (node->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)||
          node->CheckContext(ContactTopologyEntity<Real>::TIED)) node->SetGhostRCB();
    }
    
    int num_primary_faces = primary_topology->Number_of_Faces();
    ContactFace<Real>** primary_faces = 
      reinterpret_cast<ContactFace<Real>**>(primary_topology->FaceList()->EntityList());
    for (int i=0; i<num_primary_faces; ++i) {
      ContactFace<Real>* face = primary_faces[i];
      if (face->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)||
          face->CheckContext(ContactTopologyEntity<Real>::TIED)) {
        int num_nodes = face->Nodes_Per_Face();
        for (int k=0; k<num_nodes; ++k) {
          ContactNode<Real>* node = face->Node(k);
          if (node->Ownership()     != ContactTopologyEntity<Real>::OWNED) continue;
          if (node->Physical_Type() == ContactNode<Real>::SHELL_TAB_NODE)  continue;
          if (node->Physical_Type() == ContactNode<Real>::MIXED_TAB_NODE)  continue;
          node->SetGhostRCB();
        }
      }
    }
  } else {
    // do not need to include tied entities since the ghosting is not persistent
    int num_primary_nodes = primary_topology->Number_of_Nodes();
    ContactNode<Real>** primary_nodes = 
      reinterpret_cast<ContactNode<Real>**>(primary_topology->NodeList()->EntityList());
    for (int i=0; i<num_primary_nodes; ++i) {
      ContactNode<Real>* node = primary_nodes[i];
      node->ClearContextBit(ContactTopologyEntity<Real>::ACTIVE_FOR_GHOSTING_RCB);
      if (node->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)) node->SetGhostRCB();
    }
    
    int num_primary_faces = primary_topology->Number_of_Faces();
    ContactFace<Real>** primary_faces = 
      reinterpret_cast<ContactFace<Real>**>(primary_topology->FaceList()->EntityList());
    for (int i=0; i<num_primary_faces; ++i) {
      ContactFace<Real>* face = primary_faces[i];
      if (face->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)) {
        int num_nodes = face->Nodes_Per_Face();
        for (int k=0; k<num_nodes; ++k) {
          ContactNode<Real>* node = face->Node(k);
          if (node->Ownership()     != ContactTopologyEntity<Real>::OWNED) continue;
          if (node->Physical_Type() == ContactNode<Real>::SHELL_TAB_NODE)  continue;
          if (node->Physical_Type() == ContactNode<Real>::MIXED_TAB_NODE)  continue;
          node->SetGhostRCB();
        }
      }
    }
  }
  
  int num_primary_elems = primary_topology->Number_of_Elements();
  ContactElement** primary_elems = 
    reinterpret_cast<ContactElement**>(primary_topology->ElemList()->EntityList());
  for (int i=0; i<num_primary_elems; ++i) {
    ContactElement* elem = primary_elems[i];
    if (elem->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)) {
      int num_nodes = elem->Nodes_Per_Element();
      for (int k=0; k<num_nodes; ++k) {
        ContactNode<Real>* node = elem->Node(k);
        if (node->Ownership() != ContactTopologyEntity<Real>::OWNED) continue;
        node->SetGhostRCB();
      }
    }
  } 
  
  #if CONTACT_DEBUG_PRINT_LEVEL>=2 || defined(CONTACT_HEARTBEAT)
  int num_rcb_nodes = 0;
  int num_primary_nodes = primary_topology->Number_of_Nodes();
  ContactNode<Real>** primary_nodes = 
    reinterpret_cast<ContactNode<Real>**>(primary_topology->NodeList()->EntityList());
  for (int i=0; i<num_primary_nodes; ++i) {
    if (primary_nodes[i]->IsGhostRCB()) ++num_rcb_nodes;
  }
  #if CONTACT_DEBUG_PRINT_LEVEL>=2 
  postream<<"    Tagged "<<num_rcb_nodes<<" nodes for RCB\n";
  #else
  #ifdef CONTACT_HEARTBEAT
  int n = contact_global_sum(num_rcb_nodes, SearchComm );
  if (contact_processor_number(SearchComm)==0) {
    std::cout<<"    Tagged "<<n<<" nodes for RCB\n";
  }
  #endif
  #endif
  #endif
}

//
// Tag all nodes, faces, and elements if they participate
// in the global search ghosting operation.  Nodes are tagged
// if they are slaves or if their parent face/element is a
// slave. 
//
void
ContactSearch::SetContextForGlobalSearchGhosting()
{ 
  int num_primary_nodes = primary_topology->Number_of_Nodes();
  ContactNode<Real>** primary_nodes = 
    reinterpret_cast<ContactNode<Real>**>(primary_topology->NodeList()->EntityList());
  for (int i=0; i<num_primary_nodes; ++i) {
    primary_nodes[i]->ClearContextBit(ContactTopologyEntity<Real>::ACTIVE_FOR_GLOBAL_GHOSTING);
  }
      
  int num_primary_faces = primary_topology->Number_of_Faces();
  ContactFace<Real>** primary_faces = 
    reinterpret_cast<ContactFace<Real>**>(primary_topology->FaceList()->EntityList());
  for (int i=0; i<num_primary_faces; ++i) {
    primary_faces[i]->ClearContextBit(ContactTopologyEntity<Real>::ACTIVE_FOR_GLOBAL_GHOSTING);
  }
  
  int num_primary_elems = primary_topology->Number_of_Elements();
  ContactElement** primary_elems = 
    reinterpret_cast<ContactElement**>(primary_topology->ElemList()->EntityList());
  for (int i=0; i<num_primary_elems; ++i) {
    primary_elems[i]->ClearContextBit(ContactTopologyEntity<Real>::ACTIVE_FOR_GLOBAL_GHOSTING);
  }
}

//
// Tag all nodes, faces, and elements if they participate
// in the global search ghosting operation 
//
void
ContactSearch::SetContextForTrackSearchGhosting()
{ 
  int num_primary_nodes = primary_topology->Number_of_Nodes();
  ContactNode<Real>** primary_nodes = 
    reinterpret_cast<ContactNode<Real>**>(primary_topology->NodeList()->EntityList());
  for (int i=0; i<num_primary_nodes; ++i) {
    primary_nodes[i]->ClearContextBit(ContactTopologyEntity<Real>::ACTIVE_FOR_TRACK_GHOSTING);
  }
      
  int num_primary_faces = primary_topology->Number_of_Faces();
  ContactFace<Real>** primary_faces = 
    reinterpret_cast<ContactFace<Real>**>(primary_topology->FaceList()->EntityList());
  for (int i=0; i<num_primary_faces; ++i) {
    primary_faces[i]->ClearContextBit(ContactTopologyEntity<Real>::ACTIVE_FOR_TRACK_GHOSTING);
  }
  
  int num_primary_elems = primary_topology->Number_of_Elements();
  ContactElement** primary_elems = 
    reinterpret_cast<ContactElement**>(primary_topology->ElemList()->EntityList());
  for (int i=0; i<num_primary_elems; ++i) {
    primary_elems[i]->ClearContextBit(ContactTopologyEntity<Real>::ACTIVE_FOR_TRACK_GHOSTING);
  }
}

