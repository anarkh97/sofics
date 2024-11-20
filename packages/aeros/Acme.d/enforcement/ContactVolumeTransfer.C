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

#include "ContactErrors.h"
#include "ContactVolumeTransfer.h"
#include "ContactSearch.h"
#include "ContactTopology.h"
#include "ContactElementElementInteraction.h"
#include "ContactElement.h"
#include "Contact_Defines.h"
#include "ContactBoundingBox.h"
#include "ContactScratchManager.h"
#include "cstring"
#ifdef CONTACT_DEBUG_NODE
#include "ContactParOStream.h"
#endif
#include <cmath>
#include <cstdio>

ContactVolumeTransfer::ContactVolumeTransfer( const Real* enf_data, 
					    ContactSearch* Search,
			     ContactSearch::ContactErrorCode& error)
  : ContactEnforcement( error, Search, ContactEnforcement::VolumeTransfer,
			NSIZED, enf_data, false )
{  
  if (error != ContactSearch::NO_ERROR) return;
  // plot variables for regression
  number_global_plot_vars = 0;
  number_nodal_plot_vars = 0;
  number_element_plot_vars = 0;
  SAVE_NVARS = false;
  SAVE_EVARS = false;
  SAVE_CVARS = false;
  final_position = NULL;

#ifndef CONTACT_NO_MPI
  if( number_of_nodes )
    recv_node_var_set = new int[number_of_nodes];
  else
    recv_node_var_set = NULL;
#ifdef CONTACT_TIMINGS
  Register_Enforcement_Timers();
#endif
#endif
  if( contact_processor_number( communicator ) == 0 )
    Display_Enforcement_Data();
  error = ContactSearch::NO_ERROR;
}

ContactVolumeTransfer::ContactVolumeTransfer( ContactSearch* Search,
					    const Real* restart_data,
			     ContactSearch::ContactErrorCode& error)
  : ContactEnforcement( error, Search, ContactEnforcement::VolumeTransfer,
			restart_data )
{  
  if (error != ContactSearch::NO_ERROR) return;
  // plot variables for regression
  number_global_plot_vars = 0;
  number_nodal_plot_vars = 0;
  number_element_plot_vars = 0;
  SAVE_NVARS = false;
  SAVE_EVARS = false;
  SAVE_CVARS = false;
  final_position = NULL;

#ifndef CONTACT_NO_MPI
  if( number_of_nodes )
    recv_node_var_set = new int[number_of_nodes];
    
  else
    recv_node_var_set = NULL;
#ifdef CONTACT_TIMINGS
  Register_Enforcement_Timers();
#endif
#endif

  error = ContactSearch::NO_ERROR;
}

ContactVolumeTransfer::~ContactVolumeTransfer()
{
#ifndef CONTACT_NO_MPI
  delete [] recv_node_var_set;
#endif
}

ContactSearch::ContactErrorCode 
ContactVolumeTransfer::Get_Plot_Variable( Contact_VolTran_Plot_Vars var,
					  Real * buffer )
{
  
  errors->Add_Error_Message( "No CVARS are available" );
  return ContactSearch::UNIMPLEMENTED_FUNCTION;
}

ContactSearch::ContactErrorCode
ContactVolumeTransfer::Extract_General_Restart_Variable( Real* data )
{
  Real* buffer = data;
  Extract_Base_General_Restart_Variable( buffer );
  buffer += Number_Base_General_Restart_Variables();
  return ContactSearch::NO_ERROR;           
}

ContactSearch::ContactErrorCode
ContactVolumeTransfer::Implant_General_Restart_Variable( Real* data )
{
  Real* buffer = data;
  Implant_Base_General_Restart_Variable( buffer );
  buffer += Number_Base_General_Restart_Variables();
  return ContactSearch::NO_ERROR;           
}

ContactSearch::ContactErrorCode
ContactVolumeTransfer::Extract_Nodal_Restart_Variable( int n, Real* data, int* node_ids )
{
  return ContactSearch::NO_ERROR;           
}

ContactSearch::ContactErrorCode
ContactVolumeTransfer::Implant_Nodal_Restart_Variable( int n, Real* data )
{
  return ContactSearch::NO_ERROR;       
}

ContactSearch::ContactErrorCode
ContactVolumeTransfer::Extract_Edge_Restart_Variable( int n, Real* data )
{
  return ContactSearch::NO_ERROR; 
}

ContactSearch::ContactErrorCode
ContactVolumeTransfer::Implant_Edge_Restart_Variable( int n, Real* data )
{
  return ContactSearch::NO_ERROR;   
}

ContactSearch::ContactErrorCode
ContactVolumeTransfer::Extract_Face_Restart_Variable( int n, Real* data )
{
  return ContactSearch::NO_ERROR;        
}

ContactSearch::ContactErrorCode
ContactVolumeTransfer::Implant_Face_Restart_Variable( int n, Real* data )
{
  return ContactSearch::NO_ERROR;     
}

ContactSearch::ContactErrorCode
ContactVolumeTransfer::Extract_Element_Restart_Variable( int n, Real* data )
{
  return ContactSearch::NO_ERROR;        
}

ContactSearch::ContactErrorCode
ContactVolumeTransfer::Implant_Element_Restart_Variable( int n, Real* data )
{
  return ContactSearch::NO_ERROR;     
}

ContactSearch::ContactErrorCode
ContactVolumeTransfer::Extract_Restart_Data( Real* restart_data )
{
  ContactSearch::ContactErrorCode ec;
  // add base class restart data first
  ec = ContactEnforcement::Extract_Restart_Data( restart_data );
  POSTCONDITION( ContactEnforcement::Restart_Size() == Restart_Size() );
  return ec;
}


int
ContactVolumeTransfer::Number_Nodal_Restart_Variables()
{
  int num_vars = 0; // nothing

  return num_vars;
}

int ContactVolumeTransfer::Restart_Size()
{
  int size = 0;
  
  // Get base class restart size first
  size = ContactEnforcement::Restart_Size();
  
  return size;
}

ContactSearch::ContactErrorCode
ContactVolumeTransfer::Compute_Volume_Transfer( int num_node_vars,
						int num_elem_vars,
						const Real* donor_node_vars,
						const Real* donor_elem_vars,
						Real* receiver_node_vars,
						Real* receiver_elem_vars,
						Real* volume_fraction )
{

  ContactSearch::ContactErrorCode error_code = ContactSearch::NO_ERROR;
#ifdef CONTACT_HEARTBEAT 
  if (contact_processor_number(communicator)==0) {
    std::cout << "  Volume Transfer\n";
  }
#endif

  if( num_node_vars == 0 && num_elem_vars == 0 ) {
	  errors->Add_Error_Message("Nothing specified for volume transfer");
	  error_code = ContactSearch::INVALID_DATA;
	  return error_code;
  }

#ifdef CONTACT_DEBUG_NODE
  ContactParOStream& postream = ParOStream();
#endif

  VariableHandle CURRENT_POSITION =
      topology->Variable_Handle( ContactTopology::Current_Position );

  // initialize receiver vectors
  Real* local_receiver_node_vars = new Real[num_node_vars*number_of_nodes];
  Real* local_receiver_elem_vars = new Real[num_elem_vars*number_of_elements];
  Real* local_volume_fraction    = new Real[number_of_elements];
  std::memset( receiver_node_vars, 0, 
	  num_node_vars*number_of_nodes*sizeof(Real) );
  std::memset( receiver_elem_vars, 0, 
	  num_elem_vars*number_of_elements*sizeof(Real) );
  std::memset( volume_fraction, 0, 
	  number_of_elements*sizeof(Real) );
  std::memset( local_receiver_node_vars, 0, 
	  num_node_vars*number_of_nodes*sizeof(Real) );
  std::memset( local_receiver_elem_vars, 0, 
	  num_elem_vars*number_of_elements*sizeof(Real) );
  std::memset( local_volume_fraction, 0, 
	  number_of_elements*sizeof(Real) );

#ifndef CONTACT_NO_MPI
  if( number_of_nodes && num_node_vars )
    std::memset( recv_node_var_set, 0, number_of_nodes*sizeof(int) );
#endif

  // Call base class Set_Up function to prepare for enforcement
  error_code = Set_Up();
  if( error_code )
    return error_code;

  VariableHandle ELEMENT_VOLUME = topology->Variable_Handle( ContactTopology::Element_Volume );
  //
  // Get the scratch memory
  //
  ScratchVariable TOTAL_ELEMENT_DONOR_VARS;
  ScratchVariable TOTAL_NODAL_DONOR_VARS;
  ScratchVariable CLOSEST_ELEMENT_COORDS;

  TOTAL_ELEMENT_DONOR_VARS.Allocate_Scratch(number_of_total_elements, num_elem_vars, ScratchVariable::ZERO_SCRATCH);
  TOTAL_NODAL_DONOR_VARS.  Allocate_Scratch(number_of_total_nodes,    num_node_vars, ScratchVariable::ZERO_SCRATCH);
  CLOSEST_ELEMENT_COORDS.  Allocate_Scratch(number_of_total_nodes,    3,             ScratchVariable::ZERO_SCRATCH);

  // Copy the donor vars into the scratch arrays 
  if( num_node_vars ){
    int inode = 0;
    for( int nn=0 ; nn<topology->Number_of_Node_Blocks() ; nn++ ){
      ContactNode<Real>** BlockNodes = 
        reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->BlockEntityList(nn));
      for (int i=0; i<topology->NodeList()->BlockNumEntities(nn); ++i) {
        int k = BlockNodes[i]->HostGlobalArrayIndex();
        Real* n_data = TOTAL_NODAL_DONOR_VARS.Get_Scratch(BlockNodes[i]);
        for(int j=0 ; j<num_node_vars ; ++j){
          n_data[j] = donor_node_vars[k*num_node_vars+j];
        }
        inode++;
      }
    }
  }
  if( num_elem_vars ){
    int ielem = 0;
    for( int nn=0 ; nn<topology->Number_of_Element_Blocks() ; nn++ ){
      ContactElement** BlockElements = 
        reinterpret_cast<ContactElement**>(topology->ElemList()->BlockEntityList(nn));
      for (int i=0; i<topology->ElemList()->BlockNumEntities(nn); ++i) {
        int k = BlockElements[i]->HostGlobalArrayIndex();
        Real* e_data = TOTAL_ELEMENT_DONOR_VARS.Get_Scratch(BlockElements[i]);
        for(int j=0 ; j<num_elem_vars ; ++j){
          e_data[j] = donor_elem_vars[k*num_elem_vars+j];
        }
        ielem++;
      }
    }
  }

#ifndef CONTACT_NO_MPI
  // Register the imports
  if( num_node_vars ){
    contact_import_scratch_var(communicator, 
                               *Node_AsymComm, 
                               *(search->Get_Comm_Buffer()),
                               TOTAL_NODAL_DONOR_VARS);

  }
  if( num_elem_vars ){
    contact_import_scratch_var(communicator, 
                               *Element_AsymComm, 
                               *(search->Get_Comm_Buffer()),
                               TOTAL_ELEMENT_DONOR_VARS);
  }
#endif

  // perform element transfer and calculate volume fractions of
  // donor to receiver overlap
  number_element_element_interactions =
    topology->Number_ElementElement_Interactions();
  for (int i=0 ; i<number_element_element_interactions ; ++i) {
    ContactElementElementInteraction* ceei = 
      element_element_interaction_list[i];
    ContactElement *donor_element =  ceei->MasterElement();
    POSTCONDITION(donor_element);
    ContactElement *receiver_element = ceei->SlaveElement();
    POSTCONDITION(receiver_element);
    int receiver_index = receiver_element->ProcArrayIndex();

    // accumulate volume fraction
    Real *receiver_volume = 
      receiver_element->Variable(ELEMENT_VOLUME);
    Real Volume_Fraction = ceei->Scalar_Var(ContactElementElementInteraction::VOLUME)/(*receiver_volume);
    local_volume_fraction[receiver_index] += Volume_Fraction;

    // accumulate element vars
    Real* dev = TOTAL_ELEMENT_DONOR_VARS.Get_Scratch(donor_element);
    for(int j=0 ; j<num_elem_vars ; ++j) {
      local_receiver_elem_vars[receiver_index*num_elem_vars+j] += 
    	dev[j]*Volume_Fraction;
    }
  }
  // divide by sum of individual donor volume fractions for final result
  // this yields: 
  //   P_avg   = sum_k{(P_donor_k*volfrac_donor_k)/volfrac_donor_total)}
  //           = sum_k{P_donor_k*(Vol_donor_k_overlap/Vol_donor_total_overlap)}
  //   where:
  //   		 sum_k indicates a sum over all donors k overlapping a receiver
  //             P is the quantity to be volume averaged to the receiver
  //   	         volfrac_donor_k = Vol_donor_k_overlap/Vol_receiver
  //   	         volfrac_donor_total = sum_k{Vol_donor_k_overlap/Vol_receiver}
  //   Note 1: In the following receiver_elem_vars[index] = 
  //   					       sum_k{P_donor_k*volfrac_donor_k}
  for (int i=0 ; i<number_of_elements ; ++i) {
    for(int j=0 ; j<num_elem_vars ; ++j) {
      // make sure we don't have an element from donor mesh since we didn't 
      // calc a volfrac
      if(local_volume_fraction[i]) {
      	local_receiver_elem_vars[i*num_elem_vars+j] =
		local_receiver_elem_vars[i*num_elem_vars+j]/local_volume_fraction[i];
      }
    }
  }

  // perform nodal variable transfer using interpolation transfer
  // donor to receiver overlap	
  Real spatial_tolerance = 1.0e-2; // spatial tolerance for box check
  if( num_node_vars ){
    number_element_element_interactions =
      topology->Number_ElementElement_Interactions();
    for (int i=0 ; i<number_element_element_interactions ; ++i) {
      ContactElementElementInteraction* ceei = 
	element_element_interaction_list[i];
      ContactElement *donor_element =  ceei->MasterElement();
      ContactElement *receiver_element = ceei->SlaveElement();
      Real near_tolerance = Get_Near_Tolerance(ceei);
      
      // loop through receiver element nodes - is it in a donor element?
      int number_donor_node_conn = donor_element->Nodes_Per_Element();
      Real node_scalar[8];
      PRECONDITION( 8 >= number_donor_node_conn );
      int number_receiver_node_conn = receiver_element->Nodes_Per_Element();
      Real local_coords[3];
      Real global_coords[3];
      for(int j=0; j<number_receiver_node_conn; ++j) {
#ifdef CONTACT_DEBUG_NODE
	bool PRINT_THIS_NODE = 
	  topology->Is_a_Debug_Node( receiver_element->Node(j) );
#endif
	int receiver_node_id = receiver_element->Node(j)->ProcArrayIndex();
	Real *position = receiver_element->Node(j)->Variable(CURRENT_POSITION);
	global_coords[0] = position[0];
	global_coords[1] = position[1];
	global_coords[2] = position[2];
#ifdef CONTACT_DEBUG_NODE
	if( PRINT_THIS_NODE ){
	  postream << "Performing Search for Nodal transfer for Node "
		   << receiver_element->Node(j)->Exodus_ID() << "\n";
	  postream << "   Global Coordinates = " << global_coords[0] << " "
		   << global_coords[1] << " " << global_coords[2] << "\n";
	}
#endif
	// check to see if node is in box inscribed by element

        ContactBoundingBox donor_box;
	for (int k=0; k<number_donor_node_conn; ++k) {
	  Real *d_pos = donor_element->Node(k)->Variable(CURRENT_POSITION);
	  donor_box.add_point(d_pos);
	}                                 
        donor_box.add_tolerance(spatial_tolerance);
#ifdef CONTACT_DEBUG_NODE
//	if( PRINT_THIS_NODE )
//	  postream << "   Element Box Coords (min) (max) = ("
//		   << min_x << "," << min_y << "," << min_z << ") ("
//		   << max_x << "," << max_y << "," << max_z << ")\n";
#endif
        if(donor_box.overlap(global_coords)) {
	  donor_element->Compute_Local_Coordinates(CURRENT_POSITION, 
						   global_coords,
						   local_coords);
#ifdef CONTACT_DEBUG_NODE
	  if( PRINT_THIS_NODE )
	    postream << "   Element Local Coords = " << local_coords[0] 
		     << " " << local_coords[1] << "\n"
                     << "   Near Tolerance = " << near_tolerance << "\n";
#endif
	  if( donor_element->
	      Is_Local_Coordinates_Inside_Element(local_coords) ) {
	    CLOSEST_ELEMENT_COORDS.Get_Scratch(receiver_element->Node(j))[0] = 1.0;
	    Real interpolated_value;
	    // interpolate all nodal variables for this node
	    for( int m=0; m<num_node_vars; ++m) {
	      for( int k=0; k< number_donor_node_conn; ++k) {
		Real* dnv = TOTAL_NODAL_DONOR_VARS.Get_Scratch(donor_element->Node(k));
		node_scalar[k] = dnv[m];
	      }
	      receiver_element->Interpolate_Scalar_Value( local_coords, 
							  &(node_scalar[0]),
							  interpolated_value );
#ifdef CONTACT_DEBUG_NODE
	      if( PRINT_THIS_NODE )
		postream << "   Interpolated value = " << interpolated_value 
			 << "\n";
#endif
	      // copy interpolated value into receiver 
#ifndef CONTACT_NO_MPI
	      recv_node_var_set[receiver_node_id] = 1;
#endif
	      local_receiver_node_vars[receiver_node_id*num_node_vars+m] = 
		interpolated_value;
	    }
	  } else if ( CLOSEST_ELEMENT_COORDS.Get_Scratch(receiver_element->Node(j))[0] != 1.0 && 
		      donor_element->Is_Local_Coordinates_Near_Element(&(local_coords[0]), 
								       near_tolerance)
		      == true){
	    bool process = true;
            Real* closest_element = CLOSEST_ELEMENT_COORDS.Get_Scratch(receiver_element->Node(j));
	    if( closest_element[0] == -1.0 ){
	      // determine if this one is closer
	      Real max_outside_old = 0.0;
	      if( closest_element[1] <= -1.0 ) 
		max_outside_old = std::max( max_outside_old, 1.0-closest_element[1] );
	      if( closest_element[1] >= 1.0 )
		max_outside_old = std::max( max_outside_old, closest_element[1]-1.0 );
	      if( closest_element[2] <= -1.0 ) 
		max_outside_old = std::max( max_outside_old, 1.0-closest_element[2] );
	      if( closest_element[2] >= 1.0 )
		max_outside_old = std::max( max_outside_old, closest_element[2]-1.0 );
	      Real max_outside_new = 0.0;
	      if( local_coords[0] <= -1.0 ) 
		max_outside_new = std::max( max_outside_new, 1.0-local_coords[0] );
	      if( local_coords[0] >= 1.0 )
		max_outside_new = std::max( max_outside_new, local_coords[0]-1.0 );
	      if( local_coords[1] <= -1.0 ) 
		max_outside_new = std::max( max_outside_new, 1.0-local_coords[1] );
	      if( local_coords[1] >= 1.0 )
		max_outside_new = std::max( max_outside_new, local_coords[1]-1.0 );
	      if( max_outside_new >= max_outside_old ) process = false;
	    }
	    if( process ){
	      closest_element[0] = -1.0;
	      closest_element[1] = local_coords[0];
	      closest_element[2] = local_coords[1];
	      Real interpolated_value;
	      // interpolate all nodal variables for this node
	      for( int m=0; m<num_node_vars; ++m ) {
		for( int k=0; k< number_donor_node_conn; ++k) {
		  Real* dnv = TOTAL_NODAL_DONOR_VARS.Get_Scratch(donor_element->Node(k));
		  node_scalar[k] = dnv[m];
		}
		receiver_element->Interpolate_Scalar_Value( local_coords, 
							    &(node_scalar[0]),
							    interpolated_value );
#ifdef CONTACT_DEBUG_NODE
		if( PRINT_THIS_NODE )
		  postream << "   Interpolated value = " << interpolated_value 
			   << "\n";
#endif
		// copy interpolated value into receiver 
#ifndef CONTACT_NO_MPI
		recv_node_var_set[receiver_node_id] = 1;
#endif
		local_receiver_node_vars[receiver_node_id*num_node_vars+m] = 
		  interpolated_value;
	      }
	    }
	  }
	}
      }
    }
  }

#ifdef CONTACT_DEBUG_NODE
  postream.flush();
#endif

#ifndef CONTACT_NO_MPI
  if( num_node_vars ){
    swapadd_data_array( local_receiver_node_vars, num_node_vars );
    swapadd_data_array( recv_node_var_set, 1 );
    for(int i=0 ; i<number_of_nodes ; ++i){
      if( recv_node_var_set[i] > 1 ){
	for(int j=0 ; j<num_node_vars ; ++j){
	  local_receiver_node_vars[i*num_node_vars+j] /= recv_node_var_set[i];
	}
      }
    }
  }
#endif

#ifndef CONTACT_NO_EXODUS_OUTPUT
  // Write element variables and volume fraction
  if( num_elem_vars ) {
    Set_EVARS( num_elem_vars, local_receiver_elem_vars );
    Set_EVARS( 1, local_volume_fraction );
  }

  // Write nodal variables
  if( num_node_vars ) {
    Set_NVARS( num_node_vars, local_receiver_node_vars );
  }
#endif

  // Release scratch memory
  if( num_node_vars ){
    for( int nn=0 ; nn<topology->Number_of_Node_Blocks() ; nn++ ){
      ContactNode<Real>** BlockNodes = 
        reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->BlockEntityList(nn));
      for (int i=0; i<topology->NodeList()->BlockNumEntities(nn); ++i) {
        int k1 = BlockNodes[i]->HostGlobalArrayIndex();
        int k2 = BlockNodes[i]->ProcArrayIndex();
        for(int j=0 ; j<num_node_vars ; ++j){
          receiver_node_vars[k1*num_node_vars+j] = 
            local_receiver_node_vars[k2*num_node_vars+j];
        }
      }
    }
  }
  if( num_elem_vars ){
    for( int nn=0 ; nn<topology->Number_of_Element_Blocks() ; nn++ ){
      ContactElement** BlockElements = 
        reinterpret_cast<ContactElement**>(topology->ElemList()->BlockEntityList(nn));
      for (int i=0; i<topology->ElemList()->BlockNumEntities(nn); ++i) {
        int k1 = BlockElements[i]->HostGlobalArrayIndex();
        int k2 = BlockElements[i]->ProcArrayIndex();
        for(int j=0 ; j<num_elem_vars ; ++j){
          receiver_elem_vars[k1*num_elem_vars+j] = 
            local_receiver_elem_vars[k2*num_elem_vars+j];
        } 
        volume_fraction[k1] = local_volume_fraction[k2];
      }
    }
  }
  delete [] local_receiver_node_vars;
  delete [] local_receiver_elem_vars;
  delete [] local_volume_fraction;
  
  // Call the base class clean up function
  Clean_Up();

  return error_code;
}


void ContactVolumeTransfer::Set_EVARS( int num_elem_vars, Real* element_data )
{
  PRECONDITION( num_elem_vars );
  int i, j;

  int old_data_length = 0;
  // if already exists resize and add most recent data to end
  if( element_plot_vars ) {
    old_data_length = number_element_plot_vars*number_of_elements;
    Real *temp_vars = NULL;
    temp_vars = new Real[old_data_length];
    std::memcpy(temp_vars, element_plot_vars, old_data_length*sizeof(Real));
    delete [] element_plot_vars;
    element_plot_vars = NULL;
    number_element_plot_vars += num_elem_vars;
    int new_data_length = number_element_plot_vars*number_of_elements;
    element_plot_vars = new Real[new_data_length];
    std::memcpy(element_plot_vars, temp_vars, old_data_length*sizeof(Real));
    delete [] temp_vars;
  } else {
    old_data_length = 0;
    number_element_plot_vars = num_elem_vars;
    element_plot_vars = new Real[num_elem_vars*number_of_elements];
  }

  // add new data  
  int var_count = old_data_length;
  for( i=0; i<num_elem_vars; ++i) {
    for( j=0; j<number_of_elements; ++j) {
      element_plot_vars[var_count++] = 
	element_data[j*num_elem_vars+i];
    }
  } 
  return;
}

void ContactVolumeTransfer::Set_NVARS( int num_node_vars, Real* nodal_data )
{
  PRECONDITION( num_node_vars );
  int i, j;

  int old_data_length = 0;
  // if already exists resize and add most recent data to end
  if( nodal_plot_vars ) {
    old_data_length = number_nodal_plot_vars*number_of_nodes;
    Real *temp_vars = NULL;
    temp_vars = new Real[old_data_length];
    std::memcpy(temp_vars, nodal_plot_vars, old_data_length*sizeof(Real));
    delete [] nodal_plot_vars;
    nodal_plot_vars = NULL;
    number_nodal_plot_vars += num_node_vars;
    int new_data_length = number_nodal_plot_vars*number_of_nodes;
    nodal_plot_vars = new Real[new_data_length];
    std::memcpy(nodal_plot_vars, temp_vars, old_data_length*sizeof(Real));
    delete [] temp_vars;
  } else {
    old_data_length = 0;
    number_nodal_plot_vars = num_node_vars;
    nodal_plot_vars = new Real[num_node_vars*number_of_nodes];
  }

  // add new data  
  int var_count = old_data_length;
  for( i=0; i<num_node_vars; ++i) {
    for( j=0; j<number_of_nodes; ++j) {
      nodal_plot_vars[var_count++] = 
	nodal_data[j*num_node_vars+i];
    }
  } 
  return;
}


ContactSearch::ContactErrorCode 
ContactVolumeTransfer::Volume_Transfer_Set_Up(  )
{
  ContactSearch::ContactErrorCode error_code = ContactSearch::NO_ERROR;

  // Globally syncronize the error code
  error_code = (ContactSearch::ContactErrorCode) 
    contact_global_error_check(  error_code, communicator );

  return error_code;
}

void ContactVolumeTransfer::Volume_Transfer_Clean_Up()
{
  return ;
}

void ContactVolumeTransfer::Set_CVARS()
{
  PRECONDITION( SAVE_CVARS );

  // Initialize CVARS

  return;

}

ContactSearch::ContactErrorCode
ContactVolumeTransfer::Update_For_Topology_Change( int new_number_of_nodes,
						   int* old_to_new_map )
{
  // We don't need to map final_position to the new topology
  if( final_position ) delete [] final_position;
  if( new_number_of_nodes ) 
    final_position = new Real[3*new_number_of_nodes];
  else
    final_position = NULL;
  number_of_nodes = new_number_of_nodes;
  return ContactSearch::NO_ERROR;
}

void ContactVolumeTransfer::Get_Element_Plot_Variable( int var_num, Real* data )
{
  // variable numbering starts at 0
#ifndef CONTACT_NO_EXODUS_OUTPUT
  PRECONDITION( var_num < number_element_plot_vars );
  std::memcpy( data, element_plot_vars+var_num*number_of_elements, 
	  number_of_elements*sizeof(Real) );
#endif
}

void ContactVolumeTransfer::Get_Nodal_Plot_Variable( int var_num, Real* data )
{
  // variable numbering starts at 0
#ifndef CONTACT_NO_EXODUS_OUTPUT
  PRECONDITION( var_num < number_nodal_plot_vars );
  std::memcpy( data, nodal_plot_vars+var_num*number_of_nodes, 
	  number_of_nodes*sizeof(Real) );
#endif
}

Real ContactVolumeTransfer::Get_Near_Tolerance( ContactElementElementInteraction* ceei )
{
  ContactElement* melem = ceei->MasterElement();
  ContactElement* selem = ceei->SlaveElement();
  int me_key = melem->Entity_Key();
  int se_key = selem->Entity_Key();
  return Enforcement_Data( NEAR_TOLERANCE, me_key, se_key );
}

Real ContactVolumeTransfer::Enforcement_Data( Enforcement_Data_Index EDI,
					      int key_1, int key_2 )
{
  return enforcement_data.Get_Data( EDI, key_1, key_2 );
}

void ContactVolumeTransfer::Display_Enforcement_Data()
{
#ifdef CONTACT_DEBUG_NODE
  std::cout << "Contact Volume Transfer Enforcement Data:" << std::endl;
  std::cout << "  Number of Entity Keys = " << number_entity_keys << std::endl;
  for( int imelems=0; imelems<number_entity_keys; imelems++ ) {
    for( int iselems=0; iselems<number_entity_keys; iselems++ ) {
      std::cout << "    Data for Slave Elements of Entity " << iselems 
	   << " against Master Elements of Entity " << imelems << std::endl;
      std::cout << "      Near Search Tolerance = "
	   << Enforcement_Data( NEAR_TOLERANCE, imelems, iselems ) << std::endl;
    }
  }
#endif
}
