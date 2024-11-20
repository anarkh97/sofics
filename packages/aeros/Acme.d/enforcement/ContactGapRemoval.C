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


#include "ContactErrors.h"
#include "ContactGapRemoval.h"
#include "ContactTopology.h"
#include "ContactNodeFaceInteraction.h"
#include "ContactShellHandler.h"
#include "cstring"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

ContactGapRemoval::ContactGapRemoval( const Real* enf_data, 
					    ContactSearch* Search,
			     ContactSearch::ContactErrorCode& error)
  : ContactEnforcement( error, Search, ContactEnforcement::GapRemoval,
			NSIZED, enf_data ,false )
{  
  if (error != ContactSearch::NO_ERROR) return;
  // plot vars
  number_global_plot_vars = 0;
  number_element_plot_vars = 0;
  number_nodal_plot_vars = 3;

// kinematic partition, friction coefficient
  single_contact_constraints  = 
        new ContactNodeFaceInteraction*[number_of_nodes];  
  max_interactions = search->Max_Interactions_Per_Node();
  multiple_contact_constraints  = 
        new ContactNodeFaceInteraction*[number_of_nodes*max_interactions];
#ifndef CONTACT_NO_EXODUS_OUTPUT
  displ_correction = new Real[number_of_nodes*3];
#else
  displ_correction = NULL;
#endif
  if( contact_processor_number(communicator) == 0 )
    Display_Enforcement_Data();   
  error = ContactSearch::NO_ERROR;
#ifdef CONTACT_TIMINGS
  Register_Enforcement_Timers();
#endif
}

ContactGapRemoval::ContactGapRemoval( ContactSearch* Search,
				      const Real* restart_data,
			     ContactSearch::ContactErrorCode& error)
  : ContactEnforcement( error, Search, ContactEnforcement::GapRemoval,
			restart_data )
{  
  if (error != ContactSearch::NO_ERROR) return;
  // plot vars
  number_global_plot_vars = 0;
  number_element_plot_vars = 0;
  number_nodal_plot_vars = 3;

// kinematic partition, friction coefficient
  single_contact_constraints  = 
        new ContactNodeFaceInteraction*[number_of_nodes];  
  max_interactions = search->Max_Interactions_Per_Node();
  multiple_contact_constraints  = 
        new ContactNodeFaceInteraction*[number_of_nodes*max_interactions];
#ifndef CONTACT_NO_EXODUS_OUTPUT
  displ_correction = new Real[number_of_nodes*3];
#else
  displ_correction = NULL;
#endif
  error = ContactSearch::NO_ERROR;
#ifdef CONTACT_TIMINGS
  Register_Enforcement_Timers();
#endif
}

ContactGapRemoval::~ContactGapRemoval()
{
  delete [] single_contact_constraints;
  delete [] multiple_contact_constraints;
  if( displ_correction ) delete [] displ_correction;
}

void ContactGapRemoval::Get_Nodal_Plot_Variable( int var_num, Real* data )
{
#ifndef CONTACT_NO_EXODUS_OUTPUT
  PRECONDITION( var_num < 3 );
  std::memcpy( data, displ_correction+var_num*number_of_nodes, 
	  number_of_nodes*sizeof(Real) );
#endif
}


ContactSearch::ContactErrorCode
ContactGapRemoval::Extract_General_Restart_Variable( Real* data )
{
  Real* buffer = data;
  Extract_Base_General_Restart_Variable( buffer );
  buffer += Number_Base_General_Restart_Variables();
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactGapRemoval::Implant_General_Restart_Variable( Real* data )
{
  Real* buffer = data;
  Implant_Base_General_Restart_Variable( buffer );
  buffer += Number_Base_General_Restart_Variables();
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactGapRemoval::Extract_Nodal_Restart_Variable( int n, Real* data, int* node_ids )
{
  if( n>Number_Nodal_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactGapRemoval::Implant_Nodal_Restart_Variable( int n, Real* data )
{
  if( n>Number_Nodal_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactGapRemoval::Extract_Edge_Restart_Variable( int n, Real* data )
{
  if( n>Number_Edge_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactGapRemoval::Implant_Edge_Restart_Variable( int n, Real* data )
{
  if( n>Number_Edge_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactGapRemoval::Extract_Face_Restart_Variable( int n, Real* data )
{
  if( n>Number_Face_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactGapRemoval::Implant_Face_Restart_Variable( int n, Real* data )
{
  if( n>Number_Face_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactGapRemoval::Extract_Element_Restart_Variable( int n, Real* data )
{
  if( n>Number_Element_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode
ContactGapRemoval::Implant_Element_Restart_Variable( int n, Real* data )
{
  if( n>Number_Element_Restart_Variables()  || n<=0 ){
    std::sprintf( message, "Variable index %d is unreasonable", n );
    errors->Add_Error_Message(message);
    return ContactSearch::INVALID_ID;
  }
  return ContactSearch::NO_ERROR;
}

ContactSearch::ContactErrorCode 
ContactGapRemoval::Compute_Gap_Removal( int max_iterations, Real trivial_gap,
					Real* Displc  )
{
  ContactSearch::ContactErrorCode error_code = ContactSearch::NO_ERROR;
#ifdef CONTACT_HEARTBEAT 
  if (contact_processor_number(communicator)==0) {
    std::cout << "  Gap Removal\n";
  }
#endif

  int i,j;
  // Call base class Set_Up function to prepare for gap removal
  error_code = Set_Up();
  if( error_code ) return error_code;

  // Get the scratch memory 
  D_COR.         Allocate_Scratch(number_of_total_nodes,          3, ScratchVariable::ZERO_SCRATCH);
  D_CORCOR.      Allocate_Scratch(number_of_total_nodes,          3, ScratchVariable::ZERO_SCRATCH);
  TOTAL_COR.     Allocate_Scratch(number_of_total_nodes,          3, ScratchVariable::ZERO_SCRATCH);
  UPDATED_POS.   Allocate_Scratch(number_of_total_nodes,          3, ScratchVariable::ZERO_SCRATCH);
  GAP_TO_ENFORCE.Allocate_Scratch(number_node_entity_constraints, 1, ScratchVariable::ZERO_SCRATCH);

  std::memset( single_contact_constraints,0,
	  number_of_nodes*sizeof(ContactNodeFaceInteraction*) );
  std::memset( multiple_contact_constraints,0,
	  max_interactions*number_of_nodes*sizeof(ContactNodeFaceInteraction*));
  
  // Build lists grouping the nodes with single and mulitple interactions
  Gapremoval_Set_Up();

#if CONTACT_DEBUG_PRINT_LEVEL>=3
  ParOStream() << "In Compute_Gap_Removal, # single "<< num_w_sing_gap
	       << " # multiple "<< num_w_mult_gap 
	       << " # total " << number_node_entity_constraints
	       << "\n";
  ParOStream().flush();
#endif

  bool gaps_remaining = true;
  int iteration = 1;

  while( gaps_remaining && iteration <= max_iterations ){

    // Zero the corrections starting this step (by nodes)
    // We must do this because of the swapadd that takes place on ghost/phantom
    // nodes.
    for( i=0 ; i<number_of_total_nodes ; ++i){
      Real *d_correction            = D_COR.   Get_Scratch(i);
      Real* d_correction_correction = D_CORCOR.Get_Scratch(i);
      d_correction[0] = 0.0;
      d_correction[1] = 0.0;
      d_correction[2] = 0.0;
      d_correction_correction[0] = 0.0;
      d_correction_correction[1] = 0.0;
      d_correction_correction[2] = 0.0;
    }

    // Compute trial contact displacements
    Compute_Displacement_Correction();

    // SWAPADD the trial contact displacement corrections
    swapadd_node_scratch( D_COR );
    
    // Calculate the contact displacement correction corrections
    Compute_Displacement_Correction_Correction();

    // SWAPADD the slave node acceleration correction
    swapadd_node_scratch( D_CORCOR );
    
    // Compute the displacement (by nodes)
    for( i=0 ; i<number_of_nodes ; ++i){
      Real* d_correction            = D_COR.    Get_Scratch(i);
      Real* d_correction_correction = D_CORCOR. Get_Scratch(i);
      Real* total_cor               = TOTAL_COR.Get_Scratch(i);

      for( j=0 ; j<dimensionality ; ++j){
	total_cor[j] += d_correction[j]+d_correction_correction[j];
      }
    }
    
    Assemble_Corrections();

#ifndef CONTACT_NO_MPI
    contact_import_scratch_var(communicator,
                               *Node_AsymComm,
                               *(search->Get_Comm_Buffer()),
                               TOTAL_COR);
#endif

    Update_Gaps( trivial_gap, gaps_remaining );

#ifndef CONTACT_NO_MPI
    gaps_remaining = contact_global_or( gaps_remaining, communicator );
#endif

#ifdef CONTACT_DEBUG_NODE
    ParOStream() << "completed iteration " << iteration << "\n";
#endif
    iteration++;
  }
  {
    ScratchVariable* scratch_arrays[1];
    scratch_arrays[0] = &TOTAL_COR;
    Copy_Node_Vector_Scratch_to_Host_Arrays( 1, &Displc, scratch_arrays );
  }
#ifndef CONTACT_NO_EXODUS_OUTPUT
  for( i=0 ; i<number_of_nodes; ++i){
    Real* total_cor = TOTAL_COR.Get_Scratch(i);
    for( j=0 ; j<dimensionality ; ++j){
      displ_correction[j*number_of_nodes+i] = total_cor[j];
    }
  }  
#endif

#ifdef CONTACT_DEBUG_NODE
  ContactParOStream& postream = ParOStream();
  for( i=0 ; i<Number_Debug_Nodes() ; ++i){
    ContactNode<Real>* dn = Debug_Node( i );
    if( dn ){
      Real* total_cor = TOTAL_COR.Get_Scratch(dn);
      postream << "Gap Removal Displacement for Node = ("
	       << dn->Exodus_ID() << ") "
	       << total_cor[0] << " " << total_cor[1] << " " << total_cor[2]
	       << "\n";
    }
  }
  postream.flush();
#endif

  // Release the scratch memory
  D_COR.Clear_Scratch();
  D_CORCOR.Clear_Scratch();
  TOTAL_COR.Clear_Scratch();
  UPDATED_POS.Clear_Scratch();
  GAP_TO_ENFORCE.Clear_Scratch();

  // Call the base class clean up function
  Clean_Up();
  
  // NEED TO FIX THIS!!
  // This is just a temporary fix so we can quickly/easily run presto 
  // 2.4beta on the ISL problem with this version of ACME for Presto-PET.  
  // Really needto add step_number to the arg list to each search type 
  // (e.g.static_1_config, etc.) to let the host code manage this value.  
  search->StepNumber(0);

  return (ContactSearch::ContactErrorCode) 
    contact_global_error_check( error_code, communicator );
}

void ContactGapRemoval::Gapremoval_Set_Up()
{
  int j;

  num_w_sing_gap = 0;
  num_w_mult_gap = 0;

  // Pack active contact constraint lists
  int ncc, ncc_new;
  ContactNodeFaceInteraction** temp_interaction_list;
  temp_interaction_list  = new ContactNodeFaceInteraction*[max_interactions];
  int* temp_interaction_mask;
  temp_interaction_mask  = new int[max_interactions];
  ContactNode<Real>* old_node = NULL;
  ContactNode<Real>* cur_node = NULL;
  

  int number_node_face_interactions = topology->Number_NodeFace_Interactions();
  ContactNodeFaceInteraction **node_face_interaction_list = new ContactNodeFaceInteraction*[number_node_face_interactions];

  int i = 0;

  for( ContactNodeFaceInteraction *cnfi = node_entity_list.Node_Face_Iterator_Start();
       cnfi != NULL;
       cnfi = node_entity_list.Node_Face_Iterator_Next()) {

    node_face_interaction_list[i] = cnfi;
    cnfi->EnfArrayIndex(i);
    i++;
  }


  i = 0;
  do {
    ContactNodeFaceInteraction *cnfi=NULL;
    ncc = 0;  // number of valid constraints at node
    if (number_node_face_interactions > 0) {
      cnfi = node_face_interaction_list[i]; // these are ordered by node
      old_node = node_face_interaction_list[i]->Node() ;
      cur_node = old_node ;
    } 
    // collect all valid interactions at node
    while ( cur_node == old_node && cur_node != NULL )  {
      if (cnfi) {
	// valid interaction; penetrating or tied
	temp_interaction_list[ncc++] = cnfi; 
      } 
      if (++i<number_node_face_interactions) {  
         cnfi = node_face_interaction_list[i]; // look at next interaction
         cur_node = cnfi->Node();
      } else {
         cur_node = NULL ;           // no next interaction
      }
    } 
    // add constraints to lists
    if( ncc == 1 ){
        single_contact_constraints[num_w_sing_gap++] = temp_interaction_list[0];
    } else if( ncc > 1 ){
	// Make multiple constraints linearly independent (by nodes w/mult gaps)
	Make_Independent(ncc, temp_interaction_list,temp_interaction_mask );
	ncc_new = 0;
	for( j=0 ; j<ncc ; ++j){
	  cnfi = temp_interaction_list[j];
	  if( cnfi ) {
	    if( temp_interaction_mask[j] ) {
	      multiple_contact_constraints
		[max_interactions*num_w_mult_gap + ncc_new++] = cnfi;
	    }
	  }
	}
	// catch case where muliple constraints become a single constraint
	if (ncc_new == 1) {
	  single_contact_constraints[num_w_sing_gap++] = 
	    multiple_contact_constraints
            [max_interactions*num_w_mult_gap ];
          for( j=0 ; j<ncc ; ++j){
             multiple_contact_constraints
            [max_interactions*num_w_mult_gap + j ] = 0;

          }
	} else if (ncc_new > 1) { 
	  ++num_w_mult_gap;
	}
    }
  } while (i<number_node_face_interactions) ;

  delete [] node_face_interaction_list;
  delete [] temp_interaction_list;
  delete [] temp_interaction_mask;

  // Initialize the gap to enforce

  for( ContactNodeFaceInteraction *cnfi = node_entity_list.Node_Face_Iterator_Start();
       cnfi != NULL;
       cnfi = node_entity_list.Node_Face_Iterator_Next()) {
    // only remove negative gaps (i.e., penetrations)
    *(GAP_TO_ENFORCE.Get_Scratch(cnfi)) = std::min( 0.0, 
      cnfi->Scalar_Var(ContactNodeFaceInteraction::GAP_CUR) +
      cnfi->Scalar_Var(ContactNodeFaceInteraction::GAP_OLD) );
  }
}



void ContactGapRemoval::Compute_Displacement_Correction()
{
  int i,j,k;
  ContactNodeFaceInteraction *cnfi,**mcnfi;

#ifdef CONTACT_DEBUG_NODE
  ContactParOStream& postream = ParOStream();
#endif

  if ((num_w_sing_gap > 0) || (num_w_mult_gap > 0)) {
   // Form unique gap vector
   for( i=0 ; i< num_w_sing_gap ; ++i){
     cnfi               = single_contact_constraints[i];
     Real kin_partition = Kinematic_Partition(cnfi);
     Real gapn_mag = *(GAP_TO_ENFORCE.Get_Scratch(cnfi));
     Real* gapn_dir = 
       cnfi->Vector_Var(ContactNodeFaceInteraction::PUSHBACK_DIR);
     Real* d_correction = D_COR.Get_Scratch(cnfi->Node());
     for( j=0 ; j<dimensionality ; ++j){
       // Use a "-" since gap = -penetration
       d_correction[j] = -kin_partition*gapn_mag*gapn_dir[j];
     }
#ifdef CONTACT_DEBUG_NODE
     if( topology->Is_a_Debug_Node( cnfi->Node() ) ){
       postream << "Computing Displacement Correction for Node "
		<< cnfi->Node()->Exodus_ID() << "\n";
       postream << "    gapn mag   = " << gapn_mag << "\n";
       postream << "    kin part   = " << kin_partition << "\n";
       postream << "    gapn dir   = " << gapn_dir[0] << " " << gapn_dir[1]
		<< " " << gapn_dir[2] << "\n";
       postream << "    correction = " << d_correction[0] << " "
		<< d_correction[1] << " " << d_correction[2] << "\n";
     }
#endif
   }

   int ncc;
   Real mag[3];
   Real gapn[3];
   for( i=0 ; i<num_w_mult_gap ; ++i){
     mcnfi = &(multiple_contact_constraints[i*max_interactions]);
     ncc   =  Number_of_Contact_Constraints (mcnfi);
     for( j=0 ; j<ncc ; ++j){
       cnfi     = multiple_contact_constraints[i*max_interactions+j];
       Real gapn_mag = *(GAP_TO_ENFORCE.Get_Scratch(cnfi));
       Real kin_partition = Kinematic_Partition(cnfi);
       // Use a "-" since gap = -penetration
       mag [j]  = -kin_partition*gapn_mag;
     }

     Unify_gaps( mcnfi, mag, gapn );
     Partition_gap( mcnfi, mag, gapn );

     for( j=0 ; j<ncc ; ++j){
       cnfi               = multiple_contact_constraints[i*max_interactions+j];
       Real* gapn_dir     = cnfi->Vector_Var(ContactNodeFaceInteraction::PUSHBACK_DIR);
       Real* d_correction = D_COR.Get_Scratch(cnfi->Node());
       for (k=0 ; k<dimensionality ; ++k){
	 // Kinematic Partition was accounted for before unify/partition
         if (mag[j] > 0) d_correction[k] += mag[j]*gapn_dir[k];
       }
     }
   }
  }
#ifdef CONTACT_DEBUG_NODE
   postream.flush();
#endif

}



void ContactGapRemoval::Compute_Displacement_Correction_Correction()
{
#ifdef CONTACT_DEBUG_NODE
  ContactParOStream& postream = ParOStream();
#endif

  int i,j,k;
  ContactNodeFaceInteraction *cnfi, **mcnfi;
  if ((num_w_sing_gap > 0) || (num_w_mult_gap > 0)) {
   // correction for surface motion
   for( i=0 ; i<num_w_sing_gap ; ++i){
     cnfi = single_contact_constraints[i];
     Real kin_partition = Kinematic_Partition(cnfi);
     // Use a "-" since gap = -penetration
     Real gapn_mag = -(*(GAP_TO_ENFORCE.Get_Scratch(cnfi)));
     Real* cdispl_dir = 
       cnfi->Vector_Var(ContactNodeFaceInteraction::PUSHBACK_DIR);
     // We add the ms displacement_correction because it already has a minus
     // sign embedded in it due to the fact that pb_dir is for the slave node
     // and is therefore in the opposite direction to the ms correction.
     Real cdispl_mag = gapn_mag - kin_partition*gapn_mag +
                       Displacement_Correction( cnfi );
     Real* d_correction_correction  = D_CORCOR.Get_Scratch(cnfi->Node());
     if( cdispl_mag > 0.0 ){
       for( j=0 ; j<dimensionality ; ++j){
         d_correction_correction[j] += cdispl_mag*cdispl_dir[j];
       }
     }
#ifdef CONTACT_DEBUG_NODE
     if( topology->Is_a_Debug_Node( cnfi->Node() ) ){
       postream << "Computing Displacement Correction Correction for Node "
		<< cnfi->Node()->Exodus_ID() << "\n";
       postream << "    gapn mag   = " << gapn_mag << "\n";
       postream << "    kin part   = " << kin_partition << "\n";
       postream << "    cdis mag   = " << cdispl_mag << "\n";
       postream << "    cdis dir   = " << cdispl_dir[0] << " " << cdispl_dir[1]
		<< " " << cdispl_dir[2] << "\n";
       postream << "    correction = " << d_correction_correction[0] << " "
		<< d_correction_correction[1] << " " 
		<< d_correction_correction[2] << "\n";
     }
#endif
   }
   int ncc;
   Real gapn_mags[3], d_corcor_vec[3];
   for( i=0 ; i<num_w_mult_gap ; ++i){
     mcnfi = &(multiple_contact_constraints[i*max_interactions]);
     ncc   =  Number_of_Contact_Constraints (mcnfi);
     for( j=0 ; j<ncc ; ++j){
       cnfi         = multiple_contact_constraints[i*max_interactions+j];
       gapn_mags[j] = Displacement_Correction( cnfi );
     }

     Unify_gaps (mcnfi, gapn_mags, d_corcor_vec);
     Partition_gap( mcnfi, gapn_mags, d_corcor_vec );

     for( j=0 ; j<ncc ; ++j){
       cnfi         = multiple_contact_constraints[i*max_interactions+j];
       Real kin_partition = Kinematic_Partition(cnfi);
       Real gapn_mag = cnfi->Scalar_Var(ContactNodeFaceInteraction::GAP_CUR) +
	 cnfi->Scalar_Var(ContactNodeFaceInteraction::GAP_OLD);
       Real* cdispl_dir = cnfi->Vector_Var(ContactNodeFaceInteraction::PUSHBACK_DIR);
       Real cdispl_mag = gapn_mag - kin_partition*gapn_mag - gapn_mags[j];

       Real* d_correction_correction  = D_CORCOR.Get_Scratch(cnfi->Node());
       for (k=0 ; k<dimensionality ; ++k){
         d_correction_correction[k] += cdispl_mag*cdispl_dir[k];
       }
     }
   }
  }
#ifdef CONTACT_DEBUG_NODE
  postream.flush();
#endif
}





void ContactGapRemoval::Make_Independent(int ncc,
                                           ContactNodeFaceInteraction** mcnfi,
                                           int *mask)
{
  // if constraints are collinear, choose closest penetrating interaction
  const Real tol = 0.99 ; // tolerance on collinearity of constraints
  Real dot;
  
  if( ncc >= 2 ) {
    mask[0] = 1;
    mask[1] = 1;
    //    Real* pb_1 = mcnfi[0]->Vector_Var(ContactNodeFaceInteraction::PUSHBACK_DIR);
    Real* pb_1 = mcnfi[0]->Vector_Var(ContactNodeFaceInteraction::NORMAL_DIR);
    Real gap_1 = mcnfi[0]->Scalar_Var(ContactNodeFaceInteraction::GAP_CUR) +
      mcnfi[0]->Scalar_Var(ContactNodeFaceInteraction::GAP_OLD);
    //    Real* pb_2 = mcnfi[1]->Vector_Var(ContactNodeFaceInteraction::PUSHBACK_DIR);
    Real* pb_2 = mcnfi[1]->Vector_Var(ContactNodeFaceInteraction::NORMAL_DIR);
    Real gap_2 = mcnfi[1]->Scalar_Var(ContactNodeFaceInteraction::GAP_CUR) +
      mcnfi[1]->Scalar_Var(ContactNodeFaceInteraction::GAP_OLD);
    dot = pb_1[0]*pb_2[0]+pb_1[1]*pb_2[1]+pb_1[2]*pb_2[2];
    if (std::fabs(dot) > tol) {
      if (gap_1 > gap_2) {
        mask[1] = 0;
      } else {
        mask[0] = 0;
      }
    }
    if (ncc == 3) {
      mask[2] = 1;
      //      Real* pb_3=mcnfi[2]->Vector_Var(ContactNodeFaceInteraction::PUSHBACK_DIR);
      Real* pb_3 = mcnfi[2]->Vector_Var(ContactNodeFaceInteraction::NORMAL_DIR);
      Real gap_3=mcnfi[2]->Scalar_Var(ContactNodeFaceInteraction::GAP_CUR) +
	mcnfi[2]->Scalar_Var(ContactNodeFaceInteraction::GAP_OLD);
      dot = pb_3[0]*pb_2[0]+pb_3[1]*pb_2[1]+pb_3[2]*pb_2[2];
      if (std::fabs(dot) > tol) {
	if (gap_3 > gap_2) {
          mask[1] = 0;
	} else {
          mask[2] = 0;
	}
      }
      dot = pb_1[0]*pb_3[0]+pb_1[1]*pb_3[1]+pb_1[2]*pb_3[2];
      if (std::fabs(dot) > tol) {
	if (gap_1 > gap_3) {
          mask[2] = 0;
	} else {
          mask[0] = 0;
	}
      }
    }
  }
}



Real ContactGapRemoval::Kinematic_Partition( ContactNodeFaceInteraction* cnfi )
{
  ContactFace<Real>* face = cnfi->Face();
  int f_key = face->Entity_Key();
  int n_key = (int) 
    cnfi->Scalar_Var(ContactNodeFaceInteraction::NODE_ENTITY_KEY);
  return Enforcement_Data(KINEMATIC_PARTITION,
                          f_key,n_key);
}

int ContactGapRemoval::Number_of_Contact_Constraints ( 
                            ContactNodeFaceInteraction** mcnfi)
{ int ncc= 0;
  if        (mcnfi[2]  != NULL) {
    ncc = 3;
  } else if (mcnfi[1]  != NULL) {
    ncc = 2;
  } else if (mcnfi[0]  != NULL) {
    ncc = 1;
  } 
  return ncc;
}


Real ContactGapRemoval::Enforcement_Data( Enforcement_Data_Index EDI,
					     int key_1, int key_2 )
{
  return enforcement_data.Get_Data( EDI, key_1, key_2 );
}

void ContactGapRemoval::Unify_gaps ( ContactNodeFaceInteraction** mcnfi,
                                       Real* magnitudes, Real* vector)
{ // Solve  {g} [ n_i ] = {g_i}  for g
  Real unify_matrix [3] [3] ;
  Real partition_matrix  [3] [3] ;
  //  Real* pb_1 = mcnfi[0]->Vector_Var(ContactNodeFaceInteraction::PUSHBACK_DIR);
  //  Real* pb_2 = mcnfi[1]->Vector_Var(ContactNodeFaceInteraction::PUSHBACK_DIR);
  Real* pb_1 = mcnfi[0]->Vector_Var(ContactNodeFaceInteraction::NORMAL_DIR);
  Real* pb_2 = mcnfi[1]->Vector_Var(ContactNodeFaceInteraction::NORMAL_DIR);
  Real pb_3 [3] ;
  if ( Number_of_Contact_Constraints(mcnfi) == 2) {
     pb_3[0] = pb_1[1]*pb_2[2] - pb_1[2]*pb_2[1] ;
     pb_3[1] = pb_1[2]*pb_2[0] - pb_1[0]*pb_2[2] ;
     pb_3[2] = pb_1[0]*pb_2[1] - pb_1[1]*pb_2[0] ;
     Real scale = 1.0/std::sqrt(pb_3[0]*pb_3[0]+pb_3[1]*pb_3[1]+pb_3[2]*pb_3[2]);
     pb_3[0] = pb_3[0]*scale ;
     pb_3[1] = pb_3[1]*scale ;
     pb_3[2] = pb_3[2]*scale ;
     magnitudes[2] = 0.0 ;
  } else {
    //     Real* pb3 = mcnfi[2]->Vector_Var(ContactNodeFaceInteraction::PUSHBACK_DIR);
     Real* pb3 = mcnfi[2]->Vector_Var(ContactNodeFaceInteraction::NORMAL_DIR);
     pb_3[0] = pb3[0];
     pb_3[1] = pb3[1];
     pb_3[2] = pb3[2];
  }

  partition_matrix[0] [0] = pb_1[0] ;
  partition_matrix[0] [1] = pb_1[1] ;
  partition_matrix[0] [2] = pb_1[2] ;
  partition_matrix[1] [0] = pb_2[0] ;
  partition_matrix[1] [1] = pb_2[1] ;
  partition_matrix[1] [2] = pb_2[2] ;
  partition_matrix[2] [0] = pb_3[0] ; 
  partition_matrix[2] [1] = pb_3[1] ;
  partition_matrix[2] [2] = pb_3[2] ; 

  Invert_3x3Matrix(partition_matrix,unify_matrix);

  // Unify is the inverse of Partition
  vector[0] = unify_matrix[0] [0]*magnitudes [0] 
            + unify_matrix[0] [1]*magnitudes [1]
            + unify_matrix[0] [2]*magnitudes [2] ;
  vector[1] = unify_matrix[1] [0]*magnitudes [0] 
            + unify_matrix[1] [1]*magnitudes [1]
            + unify_matrix[1] [2]*magnitudes [2] ;
  vector[2] = unify_matrix[2] [0]*magnitudes [0] 
            + unify_matrix[2] [1]*magnitudes [1]
            + unify_matrix[2] [2]*magnitudes [2] ;
}

void ContactGapRemoval::Partition_gap ( ContactNodeFaceInteraction** mcnfi,
                                          Real* magnitudes, Real* vector)
{ // Solve  [n_i]T {g_i} = {g}  for g_i by least squares
  // leastsquares_matrix = [n_i][n_i]T (formed directly)
  Real leastsquares_matrix [3] [3] ;
  Real partition_matrix  [3] [3] ;
  Real* pb_1 = mcnfi[0]->Vector_Var(ContactNodeFaceInteraction::NORMAL_DIR);
  Real* pb_2 = mcnfi[1]->Vector_Var(ContactNodeFaceInteraction::NORMAL_DIR);
  Real  pb_3 [3] ;
  if ( Number_of_Contact_Constraints(mcnfi) == 2) {
     pb_3[0] = pb_1[1]*pb_2[2] - pb_1[2]*pb_2[1] ;
     pb_3[1] = pb_1[2]*pb_2[0] - pb_1[0]*pb_2[2] ;
     pb_3[2] = pb_1[0]*pb_2[1] - pb_1[1]*pb_2[0] ;
     Real scale = 1.0/std::sqrt(pb_3[0]*pb_3[0]+pb_3[1]*pb_3[1]+pb_3[2]*pb_3[2]);
     pb_3[0] = pb_3[0]*scale ;
     pb_3[1] = pb_3[1]*scale ;
     pb_3[2] = pb_3[2]*scale ;
     magnitudes[2] = 0.0 ;
  } else {
     Real* pb3 = mcnfi[2]->Vector_Var(ContactNodeFaceInteraction::NORMAL_DIR);
     pb_3[0] = pb3[0];
     pb_3[1] = pb3[1];
     pb_3[2] = pb3[2];
  }

  Real n1dn2 = 0.0;
  Real n1dn3 = 0.0;
  Real n2dn3 = 0.0;
  for(int j=0 ; j<dimensionality ; ++j){
    n1dn2 += pb_1[j]*pb_2[j];
    n1dn3 += pb_1[j]*pb_3[j];
    n2dn3 += pb_2[j]*pb_3[j];
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

  Real vdn1 = pb_1[0]*vector[0] + pb_1[1]*vector[1] + pb_1[2]*vector[2];
  Real vdn2 = pb_2[0]*vector[0] + pb_2[1]*vector[1] + pb_2[2]*vector[2];
  Real vdn3 = pb_3[0]*vector[0] + pb_3[1]*vector[1] + pb_3[2]*vector[2];

  magnitudes[0] = partition_matrix[0] [0]*vdn1 
                + partition_matrix[0] [1]*vdn2
                + partition_matrix[0] [2]*vdn3 ;
  magnitudes[1] = partition_matrix[1] [0]*vdn1 
                + partition_matrix[1] [1]*vdn2
                + partition_matrix[1] [2]*vdn3 ;
  magnitudes[2] = partition_matrix[2] [0]*vdn1 
                + partition_matrix[2] [1]*vdn2
                + partition_matrix[2] [2]*vdn3 ;
}

Real 
ContactGapRemoval::Displacement_Correction( ContactNodeFaceInteraction* cnfi)
{
  // acceleration at contact point
  int i,j;

  // Get the Node and Face for this interaction
  ContactFace<Real>* face = cnfi->Face();

  // Calculate the face distribution factors
  Real shape_functions[MAX_NODES_PER_FACE];
  PRECONDITION( MAX_NODES_PER_FACE>=face->Nodes_Per_Face() );
  Real* coordinates = cnfi->Vector_Var(ContactNodeFaceInteraction::COORDINATES);
  face->Evaluate_Shape_Functions( coordinates, shape_functions );

  // Compute the normal displacement of the master surface at the contact pt.
  Real* pb_dir = cnfi->Vector_Var(ContactNodeFaceInteraction::PUSHBACK_DIR);
  Real displ_ms_c = 0.0 ;
  for( i=0 ; i<face->Nodes_Per_Face() ; ++i){
    Real* displ_ms = D_COR.Get_Scratch(face->Node(i));
    for( j=0 ; j<dimensionality ; ++j){
      displ_ms_c += pb_dir[j]*displ_ms[j]*shape_functions[i];
    }
  }
  return displ_ms_c;
}


void ContactGapRemoval::Invert_3x3Matrix(Real M [3][3], Real M_inv [3][3])
{ 
   const Real tol = 1.0e-30  ; //tolerance on determinant
   Real cofM [3] [3];
// Cofactors
      cofM[0][0] =   M[1][1]*M[2][2] - M[1][2]*M[2][1] ;
      cofM[1][0] = -(M[0][1]*M[2][2] - M[0][2]*M[2][1]) ;
      cofM[2][0] =   M[0][1]*M[1][2] - M[0][2]*M[1][1] ;
      cofM[0][1] = -(M[1][0]*M[2][2] - M[1][2]*M[2][0]) ;
      cofM[1][1] =   M[0][0]*M[2][2] - M[0][2]*M[2][0] ;
      cofM[2][1] = -(M[0][0]*M[1][2] - M[0][2]*M[1][0]) ;
      cofM[0][2] =   M[1][0]*M[2][1] - M[1][1]*M[2][0] ;
      cofM[1][2] = -(M[0][0]*M[2][1] - M[0][1]*M[2][0]) ;
      cofM[2][2] =   M[0][0]*M[1][1] - M[0][1]*M[1][0] ;
// Determinant
      Real det = M[0][0]*cofM[0][0] + M[0][1]*cofM[0][1] + M[0][2]*cofM[0][2] ;
      if( std::fabs(det) < tol ) {
// ERROR
#if CONTACT_DEBUG_PRINT_LEVEL>=1
	std::cerr << "Determinant = 0 in ContactGapRemoval::Invert_3x3Matrix" 
	     << std::endl;
#endif
         ;
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
}
   

void ContactGapRemoval::Display_Enforcement_Data()
{
#if CONTACT_DEBUG_PRINT_LEVEL>=2
  std::cout << "Gap Removal Enforcement Data:" << std::endl;
  std::cout << "  Number of Entity Keys = " << number_entity_keys << std::endl;
  for( int ifaces=0 ; ifaces<number_entity_keys ; ifaces++ ){
    for( int inodes=0 ; inodes<number_entity_keys ; inodes++ ){
      std::cout << "    Data for Nodes of Entity " << inodes+1
           << " against Faces of Entity " << ifaces+1 << std::endl;
      std::cout << "        Kinematic Partition  = "
           << Enforcement_Data( KINEMATIC_PARTITION,
                               ifaces, inodes ) << std::endl;;
    }
  }
#endif
}                                                                               

ContactSearch::ContactErrorCode
ContactGapRemoval::Update_For_Topology_Change( int new_number_of_nodes,
					       int* old_to_new_map )
{
  // We don't need to map displ_correction to the new topology
  if( displ_correction ) delete [] displ_correction;
  if( new_number_of_nodes ) 
    displ_correction = new Real[3*new_number_of_nodes];
  else
    displ_correction = NULL;
  number_of_nodes = new_number_of_nodes;
  return ContactSearch::NO_ERROR;
}

void ContactGapRemoval::Assemble_Corrections( )
{  
  // If there are no shells, we have removed all the gaps so we can return
  ContactShellHandler* sh = topology->Shell_Handler();
  if( !sh ) return;

#ifdef CONTACT_DEBUG_NODE
  ContactParOStream& postream = ParOStream();
#endif

  // We have shells (at least on some processors), so we need to assemble
  // the displacement corrections to the host nodes and redistribute to the
  // shell nodes.  We will then recompute the gaps and make another pass
  // if necessary.

  int i,j;

  int num_host_code_nodes = sh->Number_Host_Code_Nodes();
  for( i=0 ; i<num_host_code_nodes ; ++i){
    if( sh->Num_Acme_Nodes_for_Host_Node(i) > 1 ){
      Real assembled_correction[3];
      assembled_correction[0] = 0.0;
      assembled_correction[1] = 0.0;
      assembled_correction[2] = 0.0;
#ifdef CONTACT_DEBUG_NODE
      int hi_dbg, lo_dbg;
      sh->Acme_NodeGID_for_Host_Node(i,0,hi_dbg,lo_dbg);
      ContactHostGlobalID gid_dbg(hi_dbg,lo_dbg);
      ContactNode<Real>* node_dbg = static_cast<ContactNode<Real>*>
        (search->Get_Primary_Topology()->NodeList()->Find(gid_dbg));
      bool PRINT_THIS_NODE = topology->Is_a_Debug_Node(node_dbg);
      if( PRINT_THIS_NODE ){
	postream << "Assembling Displacement Corrections for Node "
		 << node_dbg->Exodus_ID() << " (" 
		 << sh->Num_Acme_Nodes_for_Host_Node(i) 
		 << " lofted nodes)\n";
      }			       
#endif
      int num_assembling_nodes = sh->Num_Acme_Nodes_for_Host_Node(i);
      for( j=0 ; j<num_assembling_nodes ; ++j){
        int hi, lo;
        sh->Acme_NodeGID_for_Host_Node(i,j,hi,lo);
        ContactHostGlobalID gid(hi,lo);
        ContactNode<Real>* node = static_cast<ContactNode<Real>*>
          (search->Get_Primary_Topology()->NodeList()->Find(gid));
	Real* displ_cor = TOTAL_COR.Get_Scratch(node);
#ifdef CONTACT_DEBUG_NODE
	if( PRINT_THIS_NODE ){
	  postream << "  Displacement for node " << node->Global_ID() << " = "
		   << displ_cor[0] << " " << displ_cor[1] << " "
		   << displ_cor[2] << "\n";
	}
#endif
	assembled_correction[0] += displ_cor[0];
	assembled_correction[1] += displ_cor[1];
	assembled_correction[2] += displ_cor[2];
      }
      assembled_correction[0] /= num_assembling_nodes;
      assembled_correction[1] /= num_assembling_nodes;
      assembled_correction[2] /= num_assembling_nodes;

#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE ){
	postream << "  Total Correction = " << assembled_correction[0] << " "
		 << assembled_correction[1] << " " 
		 << assembled_correction[2] << "\n";
      }
#endif
      for( j=0 ; j<sh->Num_Acme_Nodes_for_Host_Node(i) ; ++j){
        int hi, lo;
        sh->Acme_NodeGID_for_Host_Node(i,j,hi,lo);
        ContactHostGlobalID gid(hi,lo);
        ContactNode<Real>* node = static_cast<ContactNode<Real>*>
          (search->Get_Primary_Topology()->NodeList()->Find(gid));
	Real* displ_cor = TOTAL_COR.Get_Scratch(node);
	displ_cor[0] = assembled_correction[0];
	displ_cor[1] = assembled_correction[1];
	displ_cor[2] = assembled_correction[2];
      }
    }
  }
#ifdef CONTACT_DEBUG_NODE
  postream.flush();
#endif

}

void ContactGapRemoval::Update_Gaps( double trivial_gap, bool& gaps_remaining )
{
  int i,j,k;

  gaps_remaining = false;
  VariableHandle CURRENT_POSITION =
          topology->Variable_Handle( ContactTopology::Current_Position );

#ifdef CONTACT_DEBUG_NODE
    ContactParOStream& postream = ParOStream();
#endif

  // Now compute the updated positions
  for( i=0 ; i<number_of_total_nodes ; ++i){
    Real* displ_cor   = TOTAL_COR.Get_Scratch(i);
    Real* updated_pos = UPDATED_POS.Get_Scratch(i);

    ContactNode<Real>* node = enforcement_node_list[i];
    Real* cur_cor     = node->Variable( CURRENT_POSITION );
    for( j=0 ; j<dimensionality ; ++j){
      updated_pos[j] = cur_cor[j] + displ_cor[j];
    }
#ifdef CONTACT_DEBUG_NODE
    if( topology->Is_a_Debug_Node( node ) ){
      postream << "Updated position for Node " << node->Exodus_ID()
	       << " = " << updated_pos[0] << " " << updated_pos[1]
	       << " " << updated_pos[2] << "\n";
      postream << "    Orig Pos = " << cur_cor[0] << " " << cur_cor[1] << " "
	       << cur_cor[2] << "\n";
      postream << "    Disp Cor = " << displ_cor[0] << " " << displ_cor[1] 
	       << " " << displ_cor[2] << "\n";
    }
#endif
  }

  for( ContactNodeFaceInteraction *cnfi = node_entity_list.Node_Face_Iterator_Start();
       cnfi != NULL;
       cnfi = node_entity_list.Node_Face_Iterator_Next()) {
    Real* ndir = cnfi->Vector_Var(ContactNodeFaceInteraction::NORMAL_DIR);
    // Compute the total gap in the normal direction
    PRECONDITION( cnfi->Face()->Nodes_Per_Face() <= 8 );
    Real shape_functions[8];
    cnfi->Face()->Evaluate_Shape_Functions( 
	   cnfi->Vector_Var(ContactNodeFaceInteraction::COORDINATES),
	   shape_functions );
    Real x_cp[3] = {0.0,0.0,0.0};
    for( k=0 ; k<cnfi->Face()->Nodes_Per_Face() ; ++k){
      Real* x_mn_K = UPDATED_POS.Get_Scratch(cnfi->Face()->Node(k));
      for( j=0 ; j<dimensionality ; ++j){ 
	x_cp[j] += shape_functions[k]*x_mn_K[j];
      }
    }
    Real* x_sn = UPDATED_POS.Get_Scratch(cnfi->Node());
    Real remaining_gap = 0.0;
    for ( j=0 ; j<dimensionality ; ++j){ 
      remaining_gap += (x_sn[j]-x_cp[j])*ndir[j];
    } 
    if( remaining_gap < 0 ){
      *(GAP_TO_ENFORCE.Get_Scratch(cnfi)) = remaining_gap;
    } else {
      *(GAP_TO_ENFORCE.Get_Scratch(cnfi)) = 0.0;
    }
    if( std::fabs(remaining_gap) > trivial_gap ) gaps_remaining = true;
#ifdef CONTACT_DEBUG_NODE
    if( topology->Is_a_Debug_Node( cnfi->Node() ) ){
      postream << "Gap to Enforce for Node " << cnfi->Node()->Exodus_ID()
	       << " and face " << cnfi->Face()->Global_ID() << " = "
	       << *(GAP_TO_ENFORCE.Get_Scratch(cnfi)) << "\n";
      postream << "    X_n = " << x_sn[0] << " " << x_sn[1] << " " << x_sn[2]
	       << "\n";
      postream << "    X_c = " << x_cp[0] << " " << x_cp[1] << " " << x_cp[2]
	       << "\n";
    }
#endif
  }

#ifdef CONTACT_DEBUG_NODE
  postream.flush();
#endif
}
