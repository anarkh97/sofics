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

  
#include "ContactEnforcement.h"
#include "ContactErrors.h"
#include "ContactEntity.h"
#include "ContactTopology.h"
#include "ContactNodeBlock.h"
#include "ContactFaceBlock.h"
#include "ContactElementBlock.h"
#include "ContactSearch.h"
#include "ContactSearchData.h"
#include "Contact_Communication.h"
#include "ContactSymComm.h"
#include "ContactAsymComm.h"
#include "ContactHostGlobalID.h"
#include "ContactZoltanComm.h"
#include "ContactTopologyEntityHash.h"
#include "ContactNodeFaceInteraction.h"
#include "ContactNodeSurfaceInteraction.h"
#include "ContactShellNode.h"
#include "ContactElementElementInteraction.h"
#include "ContactFaceFaceInteraction.h"
#include "ContactZoltan.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactSequentialAllocator.h"
#ifndef CONTACT_NO_MPI
#include "ContactEnfZoltan.h"
#endif
#include "ContactEnfModel.h"
#include "ContactTDConstantFriction.h"
#include "ContactTDPressureDependent.h"
#include "ContactTDVelocityDependent.h"
#include "ContactTDTied.h"
#include "ContactTDFrictionless.h"
#include "ContactTDSpotWeld.h"
#include "ContactTDSpringWeld.h"
#include "ContactTDAdhesion.h"
#include "ContactTDCohesiveZone.h"
#include "ContactTDJunction.h"
#include "ContactTDThreaded.h"
#include "ContactTDPVDependent.h"
#include "ContactTDAreaWeld.h"
#include "ContactTDUser.h"
#include "ContactTDEnfModel.h"
#include "ContactScratchManager.h"
#ifdef CONTACT_SIMOD_SUPPORT
#include "ContactTDShared.h"
#endif
#include "ContactShellHandler.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>

using namespace std;
using namespace ACME;

#ifndef CONTACT_NO_MPI
//
//  NKC HACK
//  Need to sort the communication lists in some instances for fast removals of duplicates.  Using a merge sort as 
//  This lists provided by Zoltan appear to often be pre sorted (in reverse order) thus quick sort might bomb.  
//  Ideally I'd like to use an STL class to do this, or at least a templatized sort.  These are apparently to be avoided
//  in ACME for the moment to prevent compiler issues.  
//
//  Actually sorting the ContactHostGlobalID objects.  Probally should use a pointer based sort here to prevent moving 
//  potentailly large amount of data around.  Hopefully these objects are small.
//  

void mergeSort(GIDSortClass *numbers, GIDSortClass *temp, int array_size)
{
  m_sort(numbers, temp, 0, array_size - 1);
}


void m_sort(GIDSortClass *numbers, GIDSortClass *temp, int left, int right)
{
  int mid;
  if (right > left) {
    mid = (right + left) / 2;
    m_sort(numbers, temp, left, mid);
    m_sort(numbers, temp, mid+1, right);

    merge(numbers, temp, left, mid+1, right);
  }
}

void merge(GIDSortClass *numbers, GIDSortClass *temp, int left, int mid, int right)
{
  int i, left_end, num_elements, tmp_pos;

  left_end = mid - 1;
  tmp_pos = left;
  num_elements = right - left + 1;

  while ((left <= left_end) && (mid <= right))
  {
    if (numbers[left].key <= numbers[mid].key)
    {
      temp[tmp_pos] = numbers[left];
      tmp_pos = tmp_pos + 1;
      left = left +1;
    }
    else
    {
      temp[tmp_pos] = numbers[mid];
      tmp_pos = tmp_pos + 1;
      mid = mid + 1;
    }
  }

  while (left <= left_end)
  {
    temp[tmp_pos] = numbers[left];
    left = left + 1;
    tmp_pos = tmp_pos + 1;
  }
  while (mid <= right)
  {
    temp[tmp_pos] = numbers[mid];
    mid = mid + 1;
    tmp_pos = tmp_pos + 1;
  }

  for (i=0; i < num_elements; ++i) {
    numbers[right] = temp[right];
    right = right - 1;
  }
}

//
//  Use a binary search to search a list of global ids.  Return the index of the item or -1 if not found.
//
int binary_search( const ContactHostGlobalID &key, GIDSortClass *data, const int &data_length ) {
  int low(0);
  int high(data_length - 1);
  while(low <= high) {
    int middle((low + high) /2);
    if(key == data[middle].key) {
      return middle;
    }  else if(key < data[middle].key) {
      high = middle - 1;
    } else {
      low = middle + 1;
    }
  }
  return -1;
}

#endif

ContactEnforcement::ContactEnforcement( ContactSearch::ContactErrorCode& error,
                                        ContactSearch* Search, 
					ContactEnforcement_Type Type,
					int size_ed_per_pair,
					const Real* Enf_Data,
					bool Use_Enforcement_Models,
					int FM_DATA_INDEX )
  : search(Search), 
    enforcement_data(search->Search_Data()->Num_Search_Entities(), size_ed_per_pair, Enf_Data ), 
    communicator(Search->SearchComm), type(Type), 
    use_enforcement_models(Use_Enforcement_Models), 
    fm_data_index( FM_DATA_INDEX ), initialized( false )
{
  error  = ContactSearch::NO_ERROR;
  errors = new ContactErrors();
  if (search->dimensionality!=3) {
    errors->Add_Error_Message("Invalid Dimesion for enforcement");
    error = ContactSearch::INVALID_DATA;
    return;
  }
  search->Register_Enforcement( this );
  topology = search->primary_topology;
  dimensionality = search->dimensionality;
  number_of_nodes = topology->Number_of_Primary_Nodes();
  ContactShellHandler* shell_handler = topology->Shell_Handler();
  if( shell_handler )
    number_host_code_nodes = shell_handler->Number_Host_Code_Nodes();
  else
    number_host_code_nodes = number_of_nodes;
  number_of_faces = topology->Number_of_Primary_Faces();
  number_of_elements = topology->Number_of_Primary_Elements();
  number_entity_keys = search->Search_Data()->Num_Search_Entities();
#ifdef CONTACT_TD_FACE_FACE_ENF
  number_face_face_interactions = 0;
  allocated_size_face_face_interactions = 0;
  face_face_interaction_list = NULL;
#endif
  number_element_element_interactions = 0;
  allocated_size_element_element_interactions = 0;
  element_element_interaction_list = NULL;
  
  element_plot_vars = NULL;
  nodal_plot_vars = NULL;

#ifndef CONTACT_NO_MPI
  int ierr = 0;
  if (contact_number_of_processors( communicator) != 1)
    zoltan = new ContactEnfZoltan( communicator, 
				   search->zoltan->Get_ZoltanPtr(), 
				   ierr );
  else
    zoltan = NULL;

  number_of_phantom_nodes = 0;
  number_of_phantom_faces = 0;
  number_of_phantom_elements = 0;
  phantom_nodes = NULL;
  phantom_faces = NULL;
  phantom_elements = NULL;
  Node_AsymComm = NULL;
  Node_AsymComm_Full = NULL;
  Face_AsymComm = NULL;
  Element_AsymComm = NULL;
#endif
  scratch_set = false;
  std::memset( num_scratch_vars, 0, NUM_SCRATCH_TYPES*sizeof(int ) );

  num_enforcement_models = 0;
  enforcement_models = NULL;

  enforcement_node_list = NULL;
  enforcement_face_list = NULL;
  enforcement_element_list = NULL;
  //Register_Enforcement_Timers();
}

ContactEnforcement::ContactEnforcement( ContactSearch::ContactErrorCode& error,
                                        ContactSearch* Search,
					ContactEnforcement_Type Type,
					const Real* restart_data,
                                        void** function_pointers )
  : search(Search), enforcement_data(), 
    communicator(Search->SearchComm), 
    type(Type), initialized(true)
{
  error  = ContactSearch::NO_ERROR;
  errors = new ContactErrors();
  if (search->dimensionality!=3) {
    errors->Add_Error_Message("Invalid Dimesion for enforcement");
    error = ContactSearch::INVALID_DATA;
    return;
  }
  // setup some internal data structures
  search->Register_Enforcement( this );
  scratch_set = false;
  topology = search->primary_topology;
  dimensionality = search->dimensionality;
  number_of_nodes = topology->Number_of_Primary_Nodes();
  ContactShellHandler* shell_handler = topology->Shell_Handler();
  if( shell_handler )
    number_host_code_nodes = shell_handler->Number_Host_Code_Nodes();
  else
    number_host_code_nodes = number_of_nodes;
  number_of_faces = topology->Number_of_Primary_Faces();
  number_of_elements = topology->Number_of_Primary_Elements();
  number_entity_keys = search->Search_Data()->Num_Search_Entities();
#ifdef CONTACT_TD_FACE_FACE_ENF
  number_face_face_interactions = 0;
  allocated_size_face_face_interactions = 0;
  face_face_interaction_list = NULL;
#endif
  number_element_element_interactions = 0;
  allocated_size_element_element_interactions = 0;
  element_element_interaction_list = NULL;
   
  element_plot_vars = NULL;
  nodal_plot_vars = NULL;
     
#ifndef CONTACT_NO_MPI
  int ierr = 0;
  if (contact_number_of_processors( communicator) != 1)
    zoltan = new ContactEnfZoltan( communicator, 
				   search->zoltan->Get_ZoltanPtr(), 
				   ierr );
  else
    zoltan = NULL;
  number_of_phantom_nodes = 0;
  number_of_phantom_faces = 0;
  number_of_phantom_elements = 0;
  phantom_nodes = NULL;
  phantom_faces = NULL;
  phantom_elements = NULL;
  Node_AsymComm = NULL;
  Node_AsymComm_Full = NULL;
  Face_AsymComm = NULL;
  Element_AsymComm = NULL;
#endif

  // Now extract the data from restart data
  int offset = 0;
  
  use_enforcement_models = (bool) restart_data[offset++];
  fm_data_index = (int) restart_data[offset++];
  
  offset += enforcement_data.Implant_Restart_Data( &restart_data[offset] );

  enforcement_models = NULL;
  num_enforcement_models = (int) restart_data[offset++];
  if( num_enforcement_models ){
    enforcement_models = new ContactEnfModel*[num_enforcement_models];
    for( int i=0 ; i<num_enforcement_models ; ++i ){
      Enforcement_Model_Types em_type = (Enforcement_Model_Types) 
	restart_data[offset++];
      switch( em_type ){
      case TD_FRICTIONLESS:
	enforcement_models[i] = new ContactTDFrictionless(topology);
	offset += enforcement_models[i]->
	  Implant_Restart_Data( &restart_data[offset] );
	break;
      case TD_CONSTANT_FRICTION:
	enforcement_models[i] = new ContactTDConstantFriction(topology);
	offset += enforcement_models[i]->
	  Implant_Restart_Data( &restart_data[offset] );
	break;
      case TD_PRESSURE_DEPENDENT:
	enforcement_models[i] = new ContactTDPressureDependent(topology);
	offset += enforcement_models[i]->
	  Implant_Restart_Data( &restart_data[offset] );
	break;
      case TD_VELOCITY_DEPENDENT:
	enforcement_models[i] = new ContactTDVelocityDependent(topology);
	offset += enforcement_models[i]->
	  Implant_Restart_Data( &restart_data[offset] );
	break;
      case TD_TIED:
	enforcement_models[i] = new ContactTDTied(topology);
	offset += enforcement_models[i]->
	  Implant_Restart_Data( &restart_data[offset] );
	break;
      case TD_SPOT_WELD:
	enforcement_models[i] = new ContactTDSpotWeld(topology);
	offset += enforcement_models[i]->
	  Implant_Restart_Data( &restart_data[offset] );
	break;
	/*
      case TD_POINT_WELD:
	enforcement_models[i] = new ContactTDSpringWeld(topology);
	offset += enforcement_models[i]->
	  Implant_Restart_Data( &restart_data[offset] );
	break;
	*/
#ifdef CONTACT_SIMOD_SUPPORT
      case TD_SHARED:
	enforcement_models[i] = new ContactTDShared(topology);
	offset += enforcement_models[i]->
	  Implant_Restart_Data( &restart_data[offset] );
//      if (function_pointers == NULL) 
//          return ContactSearch::NULL_FUNCTION_POINTER;
//      void* function_pointer = function_pointers[i];
//      if (function_pointer == NULL) 
//          return ContactSearch::NULL_FUNCTION_POINTER;
        if(function_pointers != NULL) {
          void* function_pointer = function_pointers[i];
          if(function_pointer != NULL) {
            static_cast<ContactTDShared*>
            (enforcement_models[num_enforcement_models-1])
            ->Link_External_Function(function_pointer);
          }
        }
	break;
#endif
      default:
	errors->Add_Error_Message( "Unknown enforcement mode type" );
	POSTCONDITION( 0 );
      }
    }
  }
  scratch_set = false;
  std::memset( num_scratch_vars, 0, NUM_SCRATCH_TYPES*sizeof(int ) );
  
  // Hook up the tables and "failed models"
  for( int i=0 ; i<num_enforcement_models ; ++i ){
    enforcement_models[i]->Initialize_Model( num_enforcement_models,
					     enforcement_models,
					     search->Num_Tables(),
					     search->Tables() );
  }

  //Register_Enforcement_Timers();  
}

ContactEnforcement::~ContactEnforcement() 
{
  if(search != NULL) search->Delete_Enforcement( this );
  delete errors;
  if( element_plot_vars ) delete [] element_plot_vars;
  if( nodal_plot_vars ) delete [] nodal_plot_vars;

#ifndef CONTACT_NO_MPI
  if (zoltan) delete zoltan;
  if( Node_AsymComm ) delete Node_AsymComm;
  if( Node_AsymComm_Full ) delete Node_AsymComm_Full;
  if( Face_AsymComm ) delete Face_AsymComm;
  if( Element_AsymComm ) delete Element_AsymComm;
#endif
  for( int i=0 ; i<num_enforcement_models; ++i ){
    delete enforcement_models[i];
  }
  if( enforcement_models ) delete [] enforcement_models;
}

void ContactEnforcement::Register_Enforcement_Timers() {
#ifdef CONTACT_TIMINGS
  base_class_setup                 = Timer().Register_Timer( CString("    TD Base Class Setup") );
  base_class_init                  = Timer().Register_Timer( CString("      TDBC Initialization") );
  base_class_enf_model_update      = Timer().Register_Timer( CString("      TDBC Enforcement Model Update") );
  if( contact_number_of_processors(communicator) > 1 ){
    base_class_face_import         = Timer().Register_Timer( CString("      TDBC Face Import") );
    base_class_elem_import         = Timer().Register_Timer( CString("      TDBC Elem Import") );
    base_class_node_import         = Timer().Register_Timer( CString("      TDBC Node Import") );
  }
  base_class_enf_node_list         = Timer().Register_Timer( CString("      TDBC Make Node List") );
  base_class_enf_face_list         = Timer().Register_Timer( CString("      TDBC Make Face List") );
  base_class_enf_elem_list         = Timer().Register_Timer( CString("      TDBC Make Elem List") );
  if( contact_number_of_processors(communicator) > 1 ){
    base_class_connections         = Timer().Register_Timer( CString("      TDBC Make Connections") );
    base_class_make_asym_comm      = Timer().Register_Timer( CString("      TDBC Make Asymm Comm") );
    base_class_make_asym_comm1     = Timer().Register_Timer( CString("        TDBC Make Asymm Comm1") );
    base_class_make_asym_comm1a    = Timer().Register_Timer( CString("          TDBC Make Asymm Comm1a") );
    base_class_make_asym_comm1b    = Timer().Register_Timer( CString("          TDBC Make Asymm Comm1b") );
    base_class_make_asym_comm2     = Timer().Register_Timer( CString("        TDBC Make Asymm Comm2") );
    base_class_make_asym_full_comm = Timer().Register_Timer( CString("      TDBC Make ASymm Full Comm") );
  }
  base_class_cleanup             = Timer().Register_Timer( CString("    TDBC Clean Up") );

#endif
}

void ContactEnforcement::Remove_Search_Reference(ContactSearch *search_ref) {
  if(search_ref == search) search = NULL;
}


ContactSearch::ContactErrorCode ContactEnforcement::Initialize()
{
  int i;
  ContactSearch::ContactErrorCode error = ContactSearch::NO_ERROR;
  if( use_enforcement_models ){
    error = Convert_Enf_Model_ID_to_Index();
    if( error ) return error;
    // Initialize all the models
    for( i=0 ; i<num_enforcement_models ; ++i ){
      error = enforcement_models[i]->Initialize_Model( num_enforcement_models,
						       enforcement_models,
						       search->Num_Tables(),
						       search->Tables() );
      if( error ){
	std::sprintf( message, "Initialization of model %d failed", 
		 enforcement_models[i]->ID());
	errors->Add_Error_Message( message );
	return error;
      }
    }
    // Initialized the state data
    for( i=0 ; i<num_enforcement_models ; ++i ){
      error = enforcement_models[i]->Initialize_State_Data();
      if( error ) {
        return error;
      }
    }
  }

  initialized = true;
  return error;
}

int ContactEnforcement::Restart_Size()
{
  int size = 0;

  // Account for use_friction_models, friction model id and extra ghosting dist.
  size += 2;

  // Get the size of the enforcement data
  size += enforcement_data.Restart_Size();

  // Output the enforcement model data
  size += 1; // number of enforcement models
  for( int i=0 ; i < num_enforcement_models ; ++i ){
    size += 1; // model type
    size += enforcement_models[i]->Restart_Size();
  }
  return size;
}

ContactSearch::ContactErrorCode
ContactEnforcement::Extract_Restart_Data( Real* restart_data )
{
  int i;
  int words_added = 0;

  restart_data[words_added++] = use_enforcement_models;
  restart_data[words_added++] = fm_data_index;

  words_added += 
    enforcement_data.Extract_Restart_Data( &restart_data[words_added] );

  restart_data[words_added++] = num_enforcement_models;
  for( i=0 ; i<num_enforcement_models ; ++i ){
    restart_data[words_added++] = enforcement_models[i]->Type();
    words_added += 
      enforcement_models[i]->Extract_Restart_Data( &restart_data[words_added] );
  }

  POSTCONDITION( words_added == ContactEnforcement::Restart_Size() );
  return ContactSearch::NO_ERROR;
}

int ContactEnforcement::Number_of_Errors()
{
  return errors->Number_of_Errors();
}

const char* ContactEnforcement::Error_Message( int i )
{
  // The interface uses a FORTRAN numbering
  return errors->Error_Message( i-1 );
}

ContactSearch::ContactErrorCode ContactEnforcement::Complete_Restart()
{
  // Hook up the tables and "failed models"
  for( int i=0 ; i<num_enforcement_models ; ++i ){
    enforcement_models[i]->Initialize_Model( num_enforcement_models,
					     enforcement_models,
					     search->Num_Tables(),
					     search->Tables() );
  }
  return ContactSearch::NO_ERROR;
}

#ifdef CONTACT_DEBUG_NODE
int ContactEnforcement::Number_Debug_Nodes() 
{
  return search->primary_topology->Number_Debug_Nodes();
}

ContactNode<Real>* ContactEnforcement::Debug_Node(int i)
{
  ContactNode<Real>* node = search->primary_topology->Debug_Node(i); 
  if( node ){
    // Don't return the pointer if we don't own the node
    if( node->Owner() != 
	contact_processor_number(communicator) )
      return (ContactNode<Real>*) NULL;
  }
  return node;
}
#endif

ContactSearch::ContactErrorCode
ContactEnforcement::Add_Enforcement_Model( Enforcement_Model_Types Type,
					   int ID, int* integer_data,
					   Real* real_data,
                                           void* function_pointer )
{
  int i;

  // Check that the ID is valid and unique.
  if( ID <= 0 ){
    errors->Add_Error_Message( "Enforcement Model IDs must be positive" );
    return ContactSearch::INVALID_ID;
  }
  for( i=0 ; i<num_enforcement_models ; ++i ){
    if( ID == enforcement_models[i]->ID() ){
      std::sprintf( message, "Enforcement Model ID %d has already been added", ID );
      errors->Add_Error_Message( message );
      return ContactSearch::INVALID_ID;
    }
  }

  ContactEnfModel** tmp = enforcement_models;
  enforcement_models = new ContactEnfModel*[num_enforcement_models+1];
  if( num_enforcement_models ){
    std::memcpy( enforcement_models, tmp, 
	    num_enforcement_models*sizeof(ContactEnfModel*) );
    delete [] tmp;
  }

  // Create the new model
  switch( Type ){
   case TD_FRICTIONLESS:
    enforcement_models[num_enforcement_models++] = new
      ContactTDFrictionless( ID, integer_data, real_data, topology );
    break;
  case TD_CONSTANT_FRICTION:
    enforcement_models[num_enforcement_models++] = new
      ContactTDConstantFriction( ID, integer_data, real_data, topology );
    break;
  case TD_TIED:
    enforcement_models[num_enforcement_models++] = new
      ContactTDTied( ID, integer_data, real_data, topology );
    break;
  case TD_PRESSURE_DEPENDENT:
    enforcement_models[num_enforcement_models++] = new
      ContactTDPressureDependent( ID, integer_data, real_data, topology );
    break;
  case TD_VELOCITY_DEPENDENT:
    enforcement_models[num_enforcement_models++] = new
      ContactTDVelocityDependent( ID, integer_data, real_data, topology );
    break;
  case TD_SPOT_WELD:
    enforcement_models[num_enforcement_models++] = new
      ContactTDSpotWeld( ID, integer_data, real_data, topology );
    break;
  case TD_SPRING_WELD:
    enforcement_models[num_enforcement_models++] = new
      ContactTDSpringWeld( ID, integer_data, real_data, topology );
    break;
  case TD_ADHESION:
    enforcement_models[num_enforcement_models++] = new
      ContactTDAdhesion( ID, integer_data, real_data, topology );
    break;
  case TD_COHESIVE_ZONE:
    enforcement_models[num_enforcement_models++] = new
      ContactTDCohesiveZone( ID, integer_data, real_data, topology );
    break;
  case TD_JUNCTION:
    enforcement_models[num_enforcement_models++] = new
      ContactTDJunction( ID, integer_data, real_data, topology );
    break;
  case TD_THREADED:
    enforcement_models[num_enforcement_models++] = new
      ContactTDThreaded( ID, integer_data, real_data, topology );
    break;
  case TD_PV_DEPENDENT:
    enforcement_models[num_enforcement_models++] = new
      ContactTDPVDependent( ID, integer_data, real_data, topology );
    break;
  case TD_AREA_WELD:
    enforcement_models[num_enforcement_models++] = new
      ContactTDAreaWeld( ID, integer_data, real_data, topology );
    break;
  case TD_USER:
    enforcement_models[num_enforcement_models++] = new
      ContactTDUser( ID, integer_data, real_data, topology, this );
    break;
#ifdef CONTACT_SIMOD_SUPPORT
  case TD_SHARED:
    enforcement_models[num_enforcement_models++] = new
      ContactTDShared( ID, integer_data, real_data, topology );
    if (function_pointer == NULL) return ContactSearch::NULL_FUNCTION_POINTER;
    static_cast<ContactTDShared*>(enforcement_models[num_enforcement_models-1])
       ->Link_External_Function(function_pointer);
    break;
#endif
  default:
    std::sprintf( message, "Unknown Enforcement Model Type %d", Type );
    errors->Add_Error_Message( message );
    return ContactSearch::INTERNAL_ERROR;
  }
  return ContactSearch::NO_ERROR;
}

void
ContactEnforcement::User_Initialize_Model_Fn(int id, CONTACT_INIT_MODEL_FN * fn)
{
  for( int i=0 ; i<num_enforcement_models ; ++i ){
    if( id == enforcement_models[i]->ID() ){
      ContactTDUser* model = static_cast<ContactTDUser*>
                             (enforcement_models[i]);
      model->Set_Initialize_Model_Fn(fn);
      break;
    }
  }
}

void
ContactEnforcement::User_Initialize_Time_Step_Fn(int id, CONTACT_INIT_TIME_STEP_FN * fn)
{
  for( int i=0 ; i<num_enforcement_models ; ++i ){
    if( id == enforcement_models[i]->ID() ){
      ContactTDUser* model = static_cast<ContactTDUser*>
                             (enforcement_models[i]);
      model->Set_Initialize_Time_Step_Fn(fn);
      break;
    }
  }
}

void
ContactEnforcement::User_Initialize_Node_State_Data_Fn(int id, CONTACT_INIT_NODE_STATE_DATA_FN * fn)
{
  for( int i=0 ; i<num_enforcement_models ; ++i ){
    if( id == enforcement_models[i]->ID() ){
      ContactTDUser* model = static_cast<ContactTDUser*>
                             (enforcement_models[i]);
      model->Set_Initialize_Node_State_Data_Fn(fn);
      break;
    }
  }
}

void
ContactEnforcement::User_Limit_Force_Fn(int id, CONTACT_LIMIT_FORCE_FN * fn)
{
  for( int i=0 ; i<num_enforcement_models ; ++i ){
    if( id == enforcement_models[i]->ID() ){
      ContactTDUser* model = static_cast<ContactTDUser*>
                             (enforcement_models[i]);
      model->Set_Limit_Force_Fn(fn);
      break;
    }
  }
}

void
ContactEnforcement::User_Active_Fn(int id, CONTACT_INTERACTION_ACTIVE_FN * fn)
{
  for( int i=0 ; i<num_enforcement_models ; ++i ){
    if( id == enforcement_models[i]->ID() ){
      ContactTDUser* model = static_cast<ContactTDUser*>
                             (enforcement_models[i]);
      model->Set_Active_Fn(fn);
      break;
    }
  }
}

void
ContactEnforcement::User_Interaction_Type_Fn(int id, CONTACT_INTERACTION_TYPE_FN * fn)
{
  for( int i=0 ; i<num_enforcement_models ; ++i ){
    if( id == enforcement_models[i]->ID() ){
      ContactTDUser* model = static_cast<ContactTDUser*>
                             (enforcement_models[i]);
      model->Set_Interaction_Type_Fn(fn);
      break;
    }
  }
}


ContactSearch::ContactErrorCode
ContactEnforcement::Convert_Enf_Model_ID_to_Index()
{
  ContactSearch::ContactErrorCode error_code = ContactSearch::NO_ERROR;

  int j,k;
  int num_entity_keys = search->Search_Data()->Num_Search_Entities();
  // First reset the IDs to be negative to avoid converting ID 1 to index 2
  // and then setting ID 1 to 2 which sets everything orginially with 1 or
  // 2 to 2
  for( j=0 ; j<num_entity_keys ; ++j ){
    for( k=0 ; k<num_entity_keys ; ++k ){
      int fm_id = (int) enforcement_data.Get_Data( fm_data_index, j, k );
      if( fm_id <= 0 ){
	error_code = ContactSearch::INVALID_ID;
	std::sprintf( message,"Friction Model IDs must be positive. Using %d",
		 fm_id );
	errors->Add_Error_Message( message );
      }	
      enforcement_data.Set_Data( fm_data_index, j, k, -fm_id );
    }
  }
  
  // Convert the IDs to indexes (0-N)
  for( int i=0 ; i<num_enforcement_models ; ++i ){
    int ID = enforcement_models[i]->ID();
    for( j=0 ; j<num_entity_keys ; ++j ){
      for( k=0 ; k<num_entity_keys ; ++k ){
	int fm_id = (int) enforcement_data.Get_Data( fm_data_index, j, k );
	if( fm_id == -ID ) enforcement_data.Set_Data( fm_data_index, j, k, i );
      }
    }
  }

  // Check that we found all models
  for( j=0 ; j<num_entity_keys ; ++j ){
    for( k=0 ; k<num_entity_keys ; ++k ){
      int fm_id = (int) enforcement_data.Get_Data( fm_data_index, j, k );
      if( fm_id < 0 ){
	std::sprintf( message,"Friction Model ID %d not found",-fm_id );
	errors->Add_Error_Message( message );
	error_code = ContactSearch::INVALID_ID;
      }
    }
  }

  return error_code;
}

ContactParOStream& ContactEnforcement::ParOStream() 
{
  return search->postream;
}

ContactCommBuffer* ContactEnforcement::CommBuffer()
{
  return search->comm_buffer;
}

ContactSearch::ContactErrorCode 
ContactEnforcement::Set_Up( bool Update_Enf_Model_State )
{
  //(search->postream).flush();
  //(search->postream)<<"ContactEnforcement::Set_Up()\n";
  
  ContactSearch::ContactErrorCode error_code = ContactSearch::NO_ERROR;
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
  Timer().Start_Timer( base_class_setup );
#endif

  if( !initialized ){
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Start_Timer( base_class_init );
#endif
    error_code = Initialize();
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Stop_Timer( base_class_init );
#endif
    if( error_code ) return error_code;
  }

  if( Update_Enf_Model_State ){
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Start_Timer( base_class_enf_model_update );
#endif
    for( int i=0 ; i<num_enforcement_models ; ++i )
      enforcement_models[i]->Update_State();
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Stop_Timer( base_class_enf_model_update );
#endif
  }

  int my_proc = contact_processor_number( communicator );
  
  //if (contact_processor_number(communicator)==0) {
  //  std::cout<<"  Enforcement Set_Up()\n";
  //  std::cout<<"    HaveGhosting = "<<topology->HaveGhosting()<<"\n";
  //  std::cout<<std::flush;
  //}
  if (topology->HaveGhosting()) {
    ContactSearch::Search_Option_Status status;
    search->Get_Search_Option(ContactSearch::NO_SECONDARY, status, NULL);
    #ifndef CONTACT_NO_MPI
    bool no_secondary = status==ContactSearch::ACTIVE;
    #endif
  
    number_of_nodes    = topology->Number_of_Primary_Nodes();
    number_of_faces    = topology->Number_of_Primary_Faces();
    number_of_elements = topology->Number_of_Primary_Elements();
    ContactNode<Real>**    PrimaryNodes = reinterpret_cast<ContactNode<Real>**>   (topology->PrimaryNodeList()->EntityList());
    ContactFace<Real>**    PrimaryFaces = reinterpret_cast<ContactFace<Real>**>   (topology->PrimaryFaceList()->EntityList());
    ContactElement** PrimaryElems = reinterpret_cast<ContactElement**>(topology->PrimaryElemList()->EntityList());
  
    number_of_total_nodes    = topology->Number_of_Nodes();
    number_of_total_faces    = topology->Number_of_Faces();
    number_of_total_elements = topology->Number_of_Elements();

    int num_nodes = topology->Number_of_Nodes();
    int num_faces = topology->Number_of_Faces();
    int num_elems = topology->Number_of_Elements();
    ContactNode<Real>**    Nodes = reinterpret_cast<ContactNode<Real>**>   (topology->NodeList()->EntityList());
    ContactFace<Real>**    Faces = reinterpret_cast<ContactFace<Real>**>   (topology->FaceList()->EntityList());
    ContactElement** Elems = reinterpret_cast<ContactElement**>(topology->ElemList()->EntityList());
    
    //(search->postream)<<"Have ghosting from search, no_secondary = "<<no_secondary<<"\n";
    //(search->postream)<<"number_of_nodes         = "<<number_of_nodes<<"\n";
    //(search->postream)<<"number_of_faces         = "<<number_of_faces<<"\n";
    //(search->postream)<<"number_of_elems         = "<<number_of_elements<<"\n";
    
    for (int i=0; i<num_nodes; ++i) {
      Nodes[i]->temp_tag = 0;
    }
    for (int i=0; i<num_faces; ++i) {
      Faces[i]->temp_tag = 0;
    }
    for (int i=0; i<num_elems; ++i) {
      Elems[i]->temp_tag = 0;
    }
    for (int i=0; i<topology->Number_of_Node_Blocks(); ++i) {
      ContactBlockEntityList* block_list = topology->Ghosted_Node_Block(i)->NodeList();
      block_list->IteratorStart();
      while (ContactTopologyEntity<Real>* entity=block_list->IteratorForward()) {
        entity->temp_tag |= 1;
      }
    }
    for (int i=0; i<topology->Number_of_Face_Blocks(); ++i) {
      ContactBlockEntityList* block_list = topology->Ghosted_Face_Block(i)->FaceList();
      block_list->IteratorStart();
      while (ContactTopologyEntity<Real>* entity=block_list->IteratorForward()) {
        entity->temp_tag |= 1;
      }
    }
    for (int i=0; i<topology->Number_of_Element_Blocks(); ++i) {
      ContactBlockEntityList* block_list = topology->Ghosted_Element_Block(i)->ElemList();
      block_list->IteratorStart();
      while (ContactTopologyEntity<Real>* entity=block_list->IteratorForward()) {
        entity->temp_tag |= 1;
      }
    }
    //if (status==ContactSearch::INACTIVE) {
    //  for (int i=0; i<topology->Number_of_Node_Blocks(); ++i) {
    //    ContactBlockEntityList* block_list = topology->Ghosted_Node_Block(i)->NodeList();
    //    block_list->IteratorStart();
    //    while (ContactTopologyEntity<Real>* entity=block_list->IteratorForward()) {
    //      entity->temp_tag = 0;
    //    }
    //  }
    //  for (int i=0; i<topology->Number_of_Face_Blocks(); ++i) {
    //    ContactBlockEntityList* block_list = topology->Ghosted_Face_Block(i)->FaceList();
    //    block_list->IteratorStart();
    //    while (ContactTopologyEntity<Real>* entity=block_list->IteratorForward()) {
    //      entity->temp_tag = 0;
    //    }
    //  }
    //  for (int i=0; i<topology->Number_of_Element_Blocks(); ++i) {
    //    ContactBlockEntityList* block_list = topology->Ghosted_Element_Block(i)->ElemList();
    //    block_list->IteratorStart();
    //    while (ContactTopologyEntity<Real>* entity=block_list->IteratorForward()) {
    //      entity->temp_tag = 0;
    //    }
    //  }
    //}
    
    for (int i=0; i<num_nodes; ++i) {
      ContactNode<Real>* node = Nodes[i];
      if (node->Owner()==my_proc) {
        for( int j=0 ; j<node->Number_NodeEntity_Interactions() ; ++j ){
          ContactNodeEntityInteraction *cnei = node->Get_NodeEntity_Interaction( j );
          POSTCONDITION(cnei);
          ContactNodeFaceInteraction* cnfi = dynamic_cast<ContactNodeFaceInteraction*>(cnei);
          if (cnfi) {
            ContactFace<Real>* face = cnfi->Face();
            if (face->Owner()!=my_proc) face->temp_tag |= 2;
            for (int k=0; k<face->Nodes_Per_Face(); ++k) {
              if (face->Node(k)->Owner()!=my_proc) face->Node(k)->temp_tag |= 2;
            }
          }
        }
      }
    }
    if(topology->Number_FaceFace_Interactions()){
      for (int i=0; i<num_faces; ++i) {
        ContactFace<Real>* face = Faces[i];
        ContactInteractionDLL<Real>* interactions = face->Get_FaceFace_Interactions();
        if(interactions != NULL) {
          interactions->IteratorStart();
          while (ContactInteractionEntity<Real>* interaction=interactions->IteratorForward()){
            ContactFaceFaceInteraction<Real>* cffi = static_cast<ContactFaceFaceInteraction<Real>*> (interaction);
            ContactFace<Real>* master_face = cffi->MasterFace();
            if (master_face->Owner()!=my_proc) master_face->temp_tag |= 2;
            for (int k=0; k<master_face->Nodes_Per_Face(); ++k) {
              if (master_face->Node(k)->Owner()!=my_proc) master_face->Node(k)->temp_tag |= 2;
            }
          }
	}
      }
    }
    if(topology->Number_ElementElement_Interactions()){
      for (int i=0; i<num_elems; ++i) {
        ContactElement* element = Elems[i];
        ContactInteractionDLL<Real>* interactions = element->Get_ElementElement_Interactions();
        interactions->IteratorStart();
        while (ContactInteractionEntity<Real>* interaction=interactions->IteratorForward()){
          ContactElementElementInteraction* ceei = static_cast<ContactElementElementInteraction*>(interaction);
          ContactElement* master_elem = ceei->MasterElement();
          if (master_elem->Owner()!=my_proc) master_elem->temp_tag |= 2;
          for (int k=0; k<master_elem->Nodes_Per_Element(); ++k) {
            if (master_elem->Node(k)->Owner()!=my_proc) master_elem->Node(k)->temp_tag |= 2;
          }
        }
      }
    }
    
    int number_of_ghosted_nodes = 0;
    int number_of_ghosted_faces = 0;
    int number_of_ghosted_elems = 0;
    //if (status==ContactSearch::ACTIVE) {
    //  for (int i=0; i<num_nodes; ++i) {
    //    if (Nodes[i]->temp_tag) ++number_of_ghosted_nodes;
    //  }
    //  for (int i=0; i<num_faces; ++i) {
    //    if (Faces[i]->temp_tag) ++number_of_ghosted_faces;
    //  }
    //  for (int i=0; i<num_elems; ++i) {
    //    if (Elems[i]->temp_tag) ++number_of_ghosted_elems;
    //  }
    //} else {
      for (int i=0; i<topology->Number_of_Node_Blocks(); ++i) {
        ContactBlockEntityList* block_list = topology->Ghosted_Node_Block(i)->NodeList();
        block_list->IteratorStart();
        while (ContactTopologyEntity<Real>* entity=block_list->IteratorForward()) {
          if (entity->temp_tag==3) ++number_of_ghosted_nodes;
        }
      }
      for (int i=0; i<topology->Number_of_Face_Blocks(); ++i) {
        ContactBlockEntityList* block_list = topology->Ghosted_Face_Block(i)->FaceList();
        block_list->IteratorStart();
        while (ContactTopologyEntity<Real>* entity=block_list->IteratorForward()) {
          if (entity->temp_tag==3) ++number_of_ghosted_faces;
        }
      }
      for (int i=0; i<topology->Number_of_Element_Blocks(); ++i) {
        ContactBlockEntityList* block_list = topology->Ghosted_Element_Block(i)->ElemList();
        block_list->IteratorStart();
        while (ContactTopologyEntity<Real>* entity=block_list->IteratorForward()) {
          if (entity->temp_tag==3) ++number_of_ghosted_elems;
        }
      }
    //}
    
    //(search->postream)<<"number_of_ghosted_nodes = "<<number_of_ghosted_nodes<<"\n";
    //(search->postream)<<"number_of_ghosted_faces = "<<number_of_ghosted_faces<<"\n";
    //(search->postream)<<"number_of_ghosted_elems = "<<number_of_ghosted_elems<<"\n";
    
    //
    //  Build the node entity interaction list.
    //
    number_node_entity_constraints = topology->Number_NodeEntity_Interactions();
    node_entity_list.Allocate(number_node_entity_constraints, number_of_total_nodes);
    for (int i=0; i<num_nodes; ++i) {
      ContactNode<Real>* node = Nodes[i];
      if (node->Owner()==my_proc) {
        for( int j=0 ; j<node->Number_NodeEntity_Interactions() ; ++j ){
          ContactNodeEntityInteraction *cnei = node->Get_NodeEntity_Interaction( j );
          if (cnei) node_entity_list.Add_Constraint(cnei);
        }
      }
      node_entity_list.Increment_Node();
    }
    
#ifdef CONTACT_TD_FACE_FACE_ENF
    // Build a compact list of the face-face interactions
    number_face_face_interactions = topology->Number_FaceFace_Interactions();
    face_face_interaction_list = new ContactFaceFaceInteraction<Real>* [number_face_face_interactions];
    if(number_face_face_interactions){
      int ffi_index = 0;
      for (int i=0; i<num_faces; ++i) {
        ContactFace<Real>* face = Faces[i];
        ContactInteractionDLL<Real>* interactions = face->Get_FaceFace_Interactions();
        if(interactions != NULL) {
          interactions->IteratorStart();
          while (ContactInteractionEntity<Real>* interaction=interactions->IteratorForward()){
            face_face_interaction_list[ffi_index] = static_cast<ContactFaceFaceInteraction<Real>*> (interaction);
            interaction->ProcIndex(ffi_index);
            ++ffi_index;
          }
	}
      }
    }
#endif    

    // Build a compact list of the element-element interactions
    number_element_element_interactions = topology->Number_ElementElement_Interactions();
    element_element_interaction_list = new ContactElementElementInteraction* [number_element_element_interactions];
    if( topology->Number_ElementElement_Interactions() ){
      int eei_index = 0;
      for (int i=0; i<num_elems; ++i) {
        ContactElement* element = Elems[i];
        ContactInteractionDLL<Real>* interactions = 
          element->Get_ElementElement_Interactions();
        interactions->IteratorStart();
        while (ContactInteractionEntity<Real>* interaction=interactions->IteratorForward()){
          ContactElementElementInteraction* ceei = static_cast<ContactElementElementInteraction*>(interaction);
          element_element_interaction_list[eei_index] = ceei;
          ceei->ProcIndex(eei_index);
          ++eei_index;
        }
      }
    }
  
    // Create the enforcement node list, first part of the array
    // is all the primary topology nodes and the second part of
    // the array is all the ghosted nodes
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Start_Timer( base_class_enf_node_list );
#endif
    //if (status==ContactSearch::ACTIVE) {
    //  number_of_total_nodes = number_of_nodes+number_of_ghosted_nodes;
    //  if( number_of_total_nodes ){
    //    int index = 0;
    //    enforcement_node_list = new ContactNode<Real>*[number_of_total_nodes];
    //    for (int i=0; i<number_of_nodes; ++i) {
    //      enforcement_node_list[index++] = PrimaryNodes[i];
    //    }
    //    for (int i=0; i<num_nodes; ++i) {
    //      if (Nodes[i]->temp_tag) enforcement_node_list[index++] = Nodes[i];
    //    }
    //  } else {
    //    enforcement_node_list = NULL;
    //  } 
    //} else {
      number_of_total_nodes = number_of_nodes+number_of_ghosted_nodes;
      if( number_of_total_nodes ){
        int index = 0;
        enforcement_node_list = new ContactNode<Real>*[number_of_total_nodes];
        for (int i=0; i<number_of_nodes; ++i) {
          enforcement_node_list[index++] = PrimaryNodes[i];
        }
        for (int i=0; i<topology->Number_of_Node_Blocks(); ++i) {
          ContactBlockEntityList* block_list = topology->Ghosted_Node_Block(i)->NodeList();
          block_list->IteratorStart();
          while (ContactTopologyEntity<Real>* entity=block_list->IteratorForward()) {
            if (entity->temp_tag==3) enforcement_node_list[index++] = static_cast<ContactNode<Real>*>(entity);
          }
        }
        POSTCONDITION(number_of_total_nodes==index);
      } else {
        enforcement_node_list = NULL;
      }
    //}
    for(int inode = 0; inode < number_of_total_nodes; ++inode) {
      enforcement_node_list[inode]->EnfArrayIndex(inode);
    }
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Stop_Timer( base_class_enf_node_list );
#endif

    // Create the enforcement face list
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Start_Timer( base_class_enf_face_list );
#endif
    //if (status==ContactSearch::ACTIVE) {
    //  number_of_total_faces = number_of_faces+number_of_ghosted_faces;
    //  if( number_of_total_faces ){
    //    int index = 0;
    //    enforcement_face_list = new ContactFace<Real>*[number_of_total_faces];
    //    for (int i=0; i<number_of_faces; ++i) {
    //      enforcement_face_list[index++] = PrimaryFaces[i];
    //    }
    //    for (int i=0; i<num_faces; ++i) {
    //      if (Faces[i]->temp_tag==3) enforcement_face_list[index++] = Faces[i];
    //    }
    //  } else {
    //    enforcement_face_list = NULL;
    //  }
    //} else {
      number_of_total_faces = number_of_faces+number_of_ghosted_faces;
      if( number_of_total_faces ){
        int index = 0;
        enforcement_face_list = new ContactFace<Real>*[number_of_total_faces];
        for (int i=0; i<number_of_faces; ++i) {
          enforcement_face_list[index++] = PrimaryFaces[i];
        }
        for (int i=0; i<topology->Number_of_Face_Blocks(); ++i) {
          ContactBlockEntityList* block_list = topology->Ghosted_Face_Block(i)->FaceList();
          block_list->IteratorStart();
          while (ContactTopologyEntity<Real>* entity=block_list->IteratorForward()) {
            if (entity->temp_tag==3) enforcement_face_list[index++] = static_cast<ContactFace<Real>*>(entity);
          }
        }
        POSTCONDITION(number_of_total_faces==index);
      } else {
        enforcement_face_list = NULL;
      }
    //}
    for(int iface = 0; iface < number_of_total_faces; ++iface) {
      enforcement_face_list[iface]->EnfArrayIndex(iface);
    }
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Stop_Timer( base_class_enf_face_list );
#endif

    // Create the enforcement element list
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Start_Timer( base_class_enf_elem_list );
#endif
    //if (status==ContactSearch::ACTIVE) {
    //  number_of_total_elements = number_of_elements+number_of_ghosted_elems;
    //  if( number_of_total_elements ){
    //    int index = 0;
    //    enforcement_element_list = new ContactElement*[number_of_total_elements];
    //    for (int i=0; i<number_of_elements; ++i) {
    //      enforcement_element_list[index++] = PrimaryElems[i];
    //    }
    //    for (int i=0; i<num_elems; ++i) {
    //      if (Elems[i]->temp_tag==3) enforcement_element_list[index++] = Elems[i];
    //    }
    //  } else {
    //    enforcement_element_list = NULL;
    //  }
    //} else {
      number_of_total_elements = number_of_elements+number_of_ghosted_elems;
      if( number_of_total_elements ){
        int index = 0;
        enforcement_element_list = new ContactElement*[number_of_total_elements];
        for (int i=0; i<number_of_elements; ++i) {
          enforcement_element_list[index++] = PrimaryElems[i];
        }
        for (int i=0; i<topology->Number_of_Element_Blocks(); ++i) {
          ContactBlockEntityList* block_list = topology->Ghosted_Element_Block(i)->ElemList();
          block_list->IteratorStart();
          while (ContactTopologyEntity<Real>* entity=block_list->IteratorForward()) {
            if (entity->temp_tag==3) enforcement_element_list[index++] = static_cast<ContactElement*>(entity);
          }
        }
        POSTCONDITION(number_of_total_elements==index);
      } else {
        enforcement_element_list = NULL;
      }
    //}
    for(int ielem = 0; ielem < number_of_total_elements; ++ielem) {
      enforcement_element_list[ielem]->EnfArrayIndex(ielem);
    }
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Stop_Timer( base_class_enf_elem_list );
#endif

    //(search->postream)<<"number_of_total_nodes   = "<<number_of_total_nodes<<"\n";
    //(search->postream)<<"number_of_total_faces   = "<<number_of_total_faces<<"\n";
    //(search->postream)<<"number_of_total_elems   = "<<number_of_total_elements<<"\n";
    //(search->postream).flush();

#ifndef CONTACT_NO_MPI
    if (contact_number_of_processors(communicator) > 1) {
  
      #ifdef CONTACT_TIMINGS
      Timer().Start_Timer( base_class_make_asym_comm );
      Timer().Start_Timer( base_class_make_asym_comm1 );
      Timer().Start_Timer( base_class_make_asym_comm1a );
      #endif
      
      ContactZoltanLID zoltanLID;
      ContactZoltanGID zoltanGID;
      LB_ID_TYPE zoltan_lid[ZOLTAN_LID_SIZE];
      LB_ID_TYPE zoltan_gid[ZOLTAN_GID_SIZE];
      int        zoltan_pid;
      ContactZoltanComm GhostNodes_ZoltanComm( ContactZoltanComm::ZOLTAN_IMPORT);
      #ifdef CONTACT_TIMINGS
      Timer().Stop_Timer( base_class_make_asym_comm1a );
      #endif
      if (no_secondary) {
        #ifdef CONTACT_TIMINGS
        Timer().Start_Timer( base_class_make_asym_comm1b );
        #endif
        for (int i=0; i<num_nodes; ++i) {
          if (Nodes[i]->temp_tag==3) {
            zoltan_pid = Nodes[i]->Owner();
            zoltanLID.ZoltanLID(CT_NODE, i, zoltan_lid);
            zoltanGID.ZoltanGID(CT_NODE, &Nodes[i]->Global_ID(), zoltan_gid);
            GhostNodes_ZoltanComm.Add_Import( zoltan_lid, zoltan_gid, zoltan_pid, 0 );
          }
        }
        POSTCONDITION(GhostNodes_ZoltanComm.Num_Import()==number_of_ghosted_nodes);

        int           num_import   = GhostNodes_ZoltanComm.Num_Import();
        ZOLTAN_ID_PTR import_gids  = GhostNodes_ZoltanComm.Import_GIDS();
        ZOLTAN_ID_PTR import_lids  = GhostNodes_ZoltanComm.Import_LIDS();
        int*          import_procs = GhostNodes_ZoltanComm.Import_Procs();
        int           num_export   = -1;
        ZOLTAN_ID_PTR export_gids  = NULL;
        ZOLTAN_ID_PTR export_lids  = NULL;
        int*          export_procs = NULL;
        PRECONDITION(num_import==number_of_ghosted_nodes);
        zoltan->Compute_Destinations(num_import, import_gids, 
                                     import_lids, import_procs, 
                                     &num_export, &export_gids,
                                     &export_lids, &export_procs);
        GhostNodes_ZoltanComm.Set_Export( num_export, export_gids, export_lids, export_procs );
        
        #ifdef CONTACT_TIMINGS
        Timer().Stop_Timer( base_class_make_asym_comm1b );
        Timer().Stop_Timer( base_class_make_asym_comm1 );
        Timer().Start_Timer( base_class_make_asym_comm2 );
        #endif
        Node_AsymComm = new ContactAsymComm( GhostNodes_ZoltanComm,
                                             search->Get_Zoltan()->Get_ZoltanPtr(),
                                             *(topology->NodeList()));
        #ifdef CONTACT_TIMINGS
        Timer().Stop_Timer( base_class_make_asym_comm2 );
        #endif
      } else {
        int index = 0;
        phantom_nodes = new ContactNode<Real>*[number_of_ghosted_nodes];
        for (int i=0; i<topology->Number_of_Node_Blocks(); ++i) {
          ContactBlockEntityList* block_list = topology->Ghosted_Node_Block(i)->NodeList();
          block_list->IteratorStart();
          while (ContactTopologyEntity<Real>* entity=block_list->IteratorForward()) {
            ContactNode<Real>* node = static_cast<ContactNode<Real>*>(entity);
            if (node->temp_tag) {
              phantom_nodes[index++] = node;
              zoltan_pid = node->Owner();
              zoltanLID.ZoltanLID(CT_NODE, i, zoltan_lid);
              zoltanGID.ZoltanGID(CT_NODE, &node->Global_ID(), zoltan_gid);
              GhostNodes_ZoltanComm.Add_Import( zoltan_lid, zoltan_gid, zoltan_pid );
            }
          }
        }
        ContactTopologyEntityHash ghosted_node_hash( number_of_ghosted_nodes, 
                                                    (ContactTopologyEntity<Real>**) phantom_nodes );
        Timer().Stop_Timer( base_class_make_asym_comm1 );
        Timer().Start_Timer( base_class_make_asym_comm2 );
        Node_AsymComm = new ContactAsymComm( GhostNodes_ZoltanComm,
                                             ghosted_node_hash, 
                                             *(topology->NodeList()) );
        #ifdef CONTACT_TIMINGS
        Timer().Stop_Timer( base_class_make_asym_comm2 );
        #endif
      }
      
      
      /*.....
      if( topology->Number_of_Face_Blocks() ){
        ContactTopologyEntityHash phantom_face_hash( number_of_phantom_faces, 
                                                     (ContactTopologyEntity<Real>**) phantom_faces );
        ContactTopologyEntityList *local_faces = topology->FaceList();        


        Face_AsymComm = new ContactAsymComm( *Face_ZoltanComm, phantom_face_hash, *local_faces );
        Face_AsymComm->Set_Index_From_EnfArrayIndex();
      } else {
        Face_AsymComm = NULL;
      }
      if( topology->Number_of_Element_Blocks() > 0 ){

        ContactTopologyEntityHash phantom_element_hash( number_of_phantom_elements, 
                                                       (ContactTopologyEntity<Real>**) phantom_elements );
        ContactTopologyEntityList *local_elements = topology->ElemList();        

        Element_AsymComm = new ContactAsymComm(*Element_ZoltanComm, phantom_element_hash, *local_elements);
        Element_AsymComm->Set_Index_From_EnfArrayIndex();
      } else {
        Element_AsymComm = NULL;
      }
      .....*/
      Face_AsymComm    = NULL;
      Element_AsymComm = NULL;
      #ifdef CONTACT_TIMINGS
      Timer().Stop_Timer( base_class_make_asym_comm );
      Timer().Start_Timer( base_class_make_asym_full_comm );
      #endif
      //
      // Now have root proc communicate total shared nodes to all ghosting
      // and phantoming processors so we can build a ContactSymComm for
      // swap-adds.
      // First we need to count how many phantom nodes we really have
      //
      Node_AsymComm_Full = new ContactAsymComm(*(topology->Node_Asym_Comm()), *Node_AsymComm); 
      //Node_AsymComm->Print("ContactEnforceemnt::Node_AsymComm",search->postream);
      //(search->postream).flush();
      //Node_AsymComm_Full->Print("ContactEnforcement::Node_AsymComm_Full",search->postream);
      //(search->postream).flush(); 
      Node_AsymComm->Set_Index_From_EnfArrayIndex();
      Node_AsymComm_Full->Set_Index_From_EnfArrayIndex();
      #ifdef CONTACT_TIMINGS
      Timer().Stop_Timer( base_class_make_asym_full_comm );
      #endif
    }
#endif
  
  //----------------------------------
  } else {
  
    number_of_total_nodes    = number_of_nodes;
    number_of_total_faces    = number_of_faces;
    number_of_total_elements = number_of_elements;
    
    //(search->postream)<<"Don't have ghosting from search\n";
    //(search->postream)<<"number_nodes            = "<<number_of_nodes<<"\n";
    //(search->postream)<<"number_faces            = "<<number_of_faces<<"\n";
    //(search->postream)<<"number_elems            = "<<number_of_elements<<"\n";
   
    int number_faces_on_proc  = 0;
    int number_faces_off_proc = 0;
    //
    //  Build the node entity interaction list.
    //
    number_node_entity_constraints = topology->Number_NodeEntity_Interactions();
    node_entity_list.Allocate(number_node_entity_constraints, number_of_total_nodes);

    int num_nodes = topology->Number_of_Nodes();
    int num_faces = topology->Number_of_Faces();
    int num_elems = topology->Number_of_Elements();
    ContactNode<Real>**    Nodes = reinterpret_cast<ContactNode<Real>**>   (topology->NodeList()->EntityList());
    ContactFace<Real>**    Faces = reinterpret_cast<ContactFace<Real>**>   (topology->FaceList()->EntityList());
    ContactElement** Elems = reinterpret_cast<ContactElement**>(topology->ElemList()->EntityList());
      
    for (int i=0; i<num_nodes; ++i) {
      ContactNode<Real>* node = Nodes[i];
      for( int j=0 ; j<node->Number_NodeEntity_Interactions() ; ++j ){
        ContactNodeEntityInteraction  *cnei = node->Get_NodeEntity_Interaction( j );
        if(cnei) {
          node_entity_list.Add_Constraint(cnei);
        }
        if(cnei->Get_Type() == ContactNodeEntityInteraction::NODE_FACE_INTERACTION) {
          if( cnei->Get_Entity_Data()->owner == my_proc ){
            ++number_faces_on_proc;
          } else {
            ++number_faces_off_proc;
          }
        }
      }
      node_entity_list.Increment_Node();
    }

#ifdef CONTACT_TD_FACE_FACE_ENF
    // Build a compact list of the face-face interactions
    number_face_face_interactions = topology->Number_FaceFace_Interactions();
    face_face_interaction_list = new ContactFaceFaceInteraction<Real>* [number_face_face_interactions];


    if(number_face_face_interactions){
      int ffi_index = 0;
      for (int i=0; i<num_faces; ++i) {
        ContactFace<Real>* face = Faces[i];
        ContactInteractionDLL<Real>* interactions = face->Get_FaceFace_Interactions();
        if(interactions != NULL) {
          interactions->IteratorStart();
          while (ContactInteractionEntity<Real>* interaction=interactions->IteratorForward()){
            face_face_interaction_list[ffi_index] = static_cast<ContactFaceFaceInteraction<Real>*> (interaction);
            interaction->ProcIndex(ffi_index);
            if( face_face_interaction_list[ffi_index]->MasterFaceEntityData()->owner == my_proc ){
              ++number_faces_on_proc;
            } else {
              ++number_faces_off_proc;
            }
            ++ffi_index;
          }
	}
      }
    }
#endif    

    // Build a compact list of the element-element interactions
    int number_elements_on_proc  = 0;
    int number_elements_off_proc = 0;
    number_element_element_interactions = topology->Number_ElementElement_Interactions();
    element_element_interaction_list = new ContactElementElementInteraction* [number_element_element_interactions];
    if( topology->Number_ElementElement_Interactions() ){
      int eei_index = 0;
      for (int i=0; i<num_elems; ++i) {
        ContactElement* element = Elems[i];
        ContactInteractionDLL<Real>* interactions = 
          element->Get_ElementElement_Interactions();
        interactions->IteratorStart();
        while (ContactInteractionEntity<Real>* interaction=interactions->IteratorForward()){
          ContactElementElementInteraction* ceei = static_cast<ContactElementElementInteraction*>(interaction);
          element_element_interaction_list[eei_index] = ceei;
          ceei->ProcIndex(eei_index);
          if( ceei->MasterElementEntityData()->owner == my_proc )
            ++number_elements_on_proc;
          else
            ++number_elements_off_proc;
          ++eei_index;
        }
      }
    }
    
#ifndef CONTACT_NO_MPI
    ContactZoltanLID zoltanLID;
    ContactZoltanGID zoltanGID;
    ZOLTAN_ID_TYPE zoltan_lid[ZOLTAN_LID_SIZE];
    ZOLTAN_ID_TYPE zoltan_gid[ZOLTAN_GID_SIZE];
    int            num_node_export;
    ZOLTAN_ID_PTR  export_node_gids;
    ZOLTAN_ID_PTR  export_node_lids;
    int*           export_node_procs;
    if( contact_number_of_processors(communicator) > 1 ){
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
      Timer().Start_Timer( base_class_face_import );
#endif
      if( topology->Number_of_Face_Blocks() ){
        // Build a unique list of faces to import
        Face_ZoltanComm = new ContactZoltanComm( ContactZoltanComm::ZOLTAN_IMPORT );
        for( ContactNodeFaceInteraction *cnfi = node_entity_list.Node_Face_Iterator_Start();
             cnfi != NULL;
             cnfi = node_entity_list.Node_Face_Iterator_Next()) {
          int face_proc = cnfi->FaceEntityData()->owner;
          if( face_proc != my_proc ) {
            cnfi->ZoltanFaceLID(zoltan_lid, 1);
            cnfi->ZoltanFaceGID(zoltan_gid);
            Face_ZoltanComm->Add_Import( zoltan_lid, zoltan_gid, face_proc );
          }
        }
#ifdef CONTACT_TD_FACE_FACE_ENF
        for( int i=0 ; i<number_face_face_interactions ; ++i ){
          int face_proc = face_face_interaction_list[i]->MasterFaceEntityData()->owner;
          if( face_proc != my_proc ){
            face_face_interaction_list[i]->ZoltanFaceLID(zoltan_lid, 1);
            face_face_interaction_list[i]->ZoltanFaceGID(zoltan_gid);
            Face_ZoltanComm->Add_Import( zoltan_lid, zoltan_gid, face_proc );
          }
        }
#endif
        number_of_phantom_faces = Face_ZoltanComm->Num_Import();
        
        // Allocate the memory to hold the phantom faces
        if( number_of_phantom_faces )
          phantom_faces = new ContactFace<Real>*[number_of_phantom_faces];
        num_imported_phantom_faces = 0;
        
        // Import the faces
        int           num_face_import   = Face_ZoltanComm->Num_Import();
        ZOLTAN_ID_PTR import_face_gids  = Face_ZoltanComm->Import_GIDS();
        ZOLTAN_ID_PTR import_face_lids  = Face_ZoltanComm->Import_LIDS();
        int*          import_face_procs = Face_ZoltanComm->Import_Procs();
        int           num_face_export   = -1;
        ZOLTAN_ID_PTR export_face_gids  = NULL;
        ZOLTAN_ID_PTR export_face_lids  = NULL;
        int*          export_face_procs = NULL;
        zoltan->Compute_Destinations(num_face_import, import_face_gids, 
                                     import_face_lids, import_face_procs, 
                                     &num_face_export, &export_face_gids,
                                     &export_face_lids, &export_face_procs);
        Face_ZoltanComm->Set_Export(num_face_export, export_face_gids, export_face_lids, export_face_procs);
        zoltan->Set_MultiFaceCallBacks(this);
        zoltan->Help_Migrate( num_face_import, import_face_gids, 
                              import_face_lids, import_face_procs,
                              num_face_export, export_face_gids,
                              export_face_lids, export_face_procs);
        POSTCONDITION( num_imported_phantom_faces == number_of_phantom_faces );
      } else {
        num_imported_phantom_faces = 0;
        number_of_phantom_faces = 0;
        phantom_faces = NULL;
        Face_ZoltanComm = NULL;
      }
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
      Timer().Stop_Timer( base_class_face_import );
#endif

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
      Timer().Start_Timer( base_class_elem_import );
#endif
      if( topology->Number_of_Element_Blocks() > 0 ){
        // Build a unique list of elements to import
        Element_ZoltanComm = new ContactZoltanComm(ContactZoltanComm::ZOLTAN_IMPORT );
        for( int i=0 ; i<number_element_element_interactions ; ++i ){
          int element_proc = element_element_interaction_list[i]->MasterElementEntityData()->owner;
          if( element_proc != my_proc ){
            element_element_interaction_list[i]->ZoltanElementLID(zoltan_lid, 1);
            element_element_interaction_list[i]->ZoltanElementGID(zoltan_gid);
            Element_ZoltanComm->Add_Import( zoltan_lid, zoltan_gid, element_proc );
          }
        }
        number_of_phantom_elements = Element_ZoltanComm->Num_Import();
        
        // Allocate the memory to hold the phantom elements
        if( number_of_phantom_elements )
          phantom_elements = new ContactElement*[number_of_phantom_elements];
        num_imported_phantom_elements = 0;
        
        // Import the elements
        int           num_element_import   = Element_ZoltanComm->Num_Import();
        ZOLTAN_ID_PTR import_element_gids  = Element_ZoltanComm->Import_GIDS();
        ZOLTAN_ID_PTR import_element_lids  = Element_ZoltanComm->Import_LIDS();
        int*          import_element_procs = Element_ZoltanComm->Import_Procs();
        int           num_element_export   = -1;
        ZOLTAN_ID_PTR export_element_gids  = NULL;
        ZOLTAN_ID_PTR export_element_lids  = NULL;
        int*          export_element_procs = NULL;
        zoltan->Compute_Destinations(num_element_import, import_element_gids, 
                                     import_element_lids, import_element_procs, 
                                     &num_element_export, &export_element_gids,
                                     &export_element_lids, &export_element_procs);

     Element_ZoltanComm->Set_Export(num_element_export, export_element_gids, export_element_lids, export_element_procs);
        zoltan->Set_MultiElementCallBacks(this);
        zoltan->Help_Migrate( num_element_import, import_element_gids, 
                              import_element_lids, import_element_procs,
                              num_element_export, export_element_gids,
                              export_element_lids, export_element_procs);
        POSTCONDITION(num_imported_phantom_elements==number_of_phantom_elements);
      } else {
        number_of_phantom_elements = 0;
        num_imported_phantom_elements = 0;
        phantom_elements = NULL;
        Element_ZoltanComm = NULL;
      }
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
      Timer().Stop_Timer( base_class_elem_import );
#endif

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
      Timer().Start_Timer( base_class_node_import );
#endif

      Node_ZoltanComm = new ContactZoltanComm( ContactZoltanComm::ZOLTAN_IMPORT );
      ContactTopologyEntityList* primary_node_list = search->primary_topology->NodeList();
      for (int i=0; i<ZOLTAN_LID_SIZE; ++i) zoltan_lid[i] = 0;
      for( int i=0 ; i<number_of_phantom_faces ; ++i ){
        ContactTopologyEntity<Real>::connection_data* node_info = phantom_faces[i]->NodeInfo();
        for( int j=0 ; j<phantom_faces[i]->Nodes_Per_Face() ; ++j ){
          ContactHostGlobalID GID( node_info[j].host_gid[0], 
                                   node_info[j].host_gid[1] );
          int node_proc = node_info[j].owner;
          if( node_proc != my_proc ) {
            //See if I already have this node (if so, don't add it to be imported)
            if( primary_node_list->Find( GID ) == NULL ){
              if(phantom_faces[i]->FaceType() == ContactSearch::SHELLQUADFACEL4 || phantom_faces[i]->FaceType() == ContactSearch::SHELLTRIFACEL3) {
                zoltanLID.ZoltanLID(CT_SHELL_NODE, node_info[j].owner_proc_array_index, zoltan_lid);
                zoltanGID.ZoltanGID(CT_SHELL_NODE, &GID, zoltan_gid);
              } else {

                zoltanLID.ZoltanLID(CT_NODE, node_info[j].owner_proc_array_index, zoltan_lid);
                zoltanGID.ZoltanGID(CT_NODE, &GID, zoltan_gid);
              }
              Node_ZoltanComm->Add_Import( zoltan_lid, zoltan_gid, node_proc );
            }
          }
        }
      }
      for( int i=0 ; i<number_of_phantom_elements ; ++i ){
        ContactTopologyEntity<Real>::connection_data* node_info = phantom_elements[i]->NodeInfo();
        for( int j=0 ; j<phantom_elements[i]->Nodes_Per_Element() ; ++j ){
          ContactHostGlobalID GID( node_info[j].host_gid[0], 
                                   node_info[j].host_gid[1] );
          int node_proc = node_info[j].owner;
          if( node_proc != my_proc ) {
            //See if I alread have this node (if so, don't add it to be imported) 
            if( primary_node_list->Find( GID ) == NULL ){
              zoltanLID.ZoltanLID(CT_SHELL_NODE, node_info[j].owner_proc_array_index, zoltan_lid);
              zoltanGID.ZoltanGID(CT_SHELL_NODE, &GID, zoltan_gid);
              Node_ZoltanComm->Add_Import( zoltan_lid, zoltan_gid, node_proc );
            }
          }
        }
      }
      number_of_phantom_nodes = Node_ZoltanComm->Num_Import();
      
      // Allocate the memory to hold the phantom nodes
      if( number_of_phantom_nodes ) phantom_nodes = new ContactNode<Real>*[number_of_phantom_nodes];
      num_imported_phantom_nodes = 0;
      
      // Import the nodes
      int           num_node_import   = Node_ZoltanComm->Num_Import();
      ZOLTAN_ID_PTR import_node_gids  = Node_ZoltanComm->Import_GIDS();
      ZOLTAN_ID_PTR import_node_lids  = Node_ZoltanComm->Import_LIDS();
      int*          import_node_procs = Node_ZoltanComm->Import_Procs();
      zoltan->Compute_Destinations(num_node_import, import_node_gids, 
                                   import_node_lids, import_node_procs, 
                                   &num_node_export, &export_node_gids,
                                   &export_node_lids, &export_node_procs);
      Node_ZoltanComm->Set_Export(num_node_export, export_node_gids, export_node_lids, export_node_procs);
      zoltan->Set_MultiNodeCallBacks(this);
      zoltan->Help_Migrate( num_node_import, import_node_gids, 
                            import_node_lids, import_node_procs,
                            num_node_export, export_node_gids,
                            export_node_lids, export_node_procs);
      POSTCONDITION( num_imported_phantom_nodes == number_of_phantom_nodes );
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
      Timer().Stop_Timer( base_class_node_import );
#endif
    } else {
      Node_ZoltanComm = NULL;
      number_of_phantom_nodes = 0;
      number_of_phantom_faces = 0;
      phantom_nodes = NULL;
      phantom_faces = NULL;
    }

#endif

    // Create the enforcement node list
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Start_Timer( base_class_enf_node_list );
#endif
#ifndef CONTACT_NO_MPI
    number_of_total_nodes += number_of_phantom_nodes;
#endif

    if( number_of_total_nodes ){
      int index = 0;
      enforcement_node_list = new ContactNode<Real>*[number_of_total_nodes];
      for (int i=0; i<num_nodes; ++i) {
        ContactNode<Real>* node = Nodes[i];
        enforcement_node_list[index++] = node;
      }
#ifndef CONTACT_NO_MPI
      std::memcpy( enforcement_node_list+number_of_nodes, phantom_nodes,
              number_of_phantom_nodes*sizeof(ContactNode<Real>*) );
#endif
    } else {
      enforcement_node_list = NULL;
    } 

    for(int inode = 0; inode < number_of_total_nodes; ++inode) {
      enforcement_node_list[inode]->EnfArrayIndex(inode);
    }
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Stop_Timer( base_class_enf_node_list );
#endif

    // Create the enforcement face list
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Start_Timer( base_class_enf_face_list );
#endif
#ifndef CONTACT_NO_MPI
    number_of_total_faces += number_of_phantom_faces;
#endif
    if( number_of_total_faces ){
      int index = 0;
      enforcement_face_list = new ContactFace<Real>*[number_of_total_faces];
      for (int i=0; i<num_faces; ++i) {
        enforcement_face_list[index++] = Faces[i];
      }
#ifndef CONTACT_NO_MPI
      std::memcpy( enforcement_face_list+number_of_faces, phantom_faces,
              number_of_phantom_faces*sizeof(ContactFace<Real>*) );
#endif
    } else {
      enforcement_face_list = NULL;
    }

    for(int iface = 0; iface < number_of_total_faces; ++iface) {
      enforcement_face_list[iface]->EnfArrayIndex(iface);
    }

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Stop_Timer( base_class_enf_face_list );
#endif

    // Create the enforcement element list
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Start_Timer( base_class_enf_elem_list );
#endif
#ifndef CONTACT_NO_MPI
    number_of_total_elements += number_of_phantom_elements;
#endif
    if( number_of_total_elements ){
      int index = 0;
      enforcement_element_list = new ContactElement*[number_of_total_elements];
      for (int i=0; i<num_elems; ++i) {
        enforcement_element_list[index++] = Elems[i];
      }
#ifndef CONTACT_NO_MPI
      std::memcpy( enforcement_element_list+number_of_elements, phantom_elements,
              number_of_phantom_elements*sizeof(ContactElement*) );
#endif
    } else {
      enforcement_element_list = NULL;
    }

    for(int ielem = 0; ielem < number_of_total_elements; ++ielem) {
      enforcement_element_list[ielem]->EnfArrayIndex(ielem);
    }
 
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Stop_Timer( base_class_enf_elem_list );
#endif

    //(search->postream)<<"number_of_phantom_nodes = "<<number_of_phantom_nodes<<"\n";
    //(search->postream)<<"number_of_phantom_faces = "<<number_of_phantom_faces<<"\n";
    //(search->postream)<<"number_of_phantom_elems = "<<number_of_phantom_elements<<"\n";
    //(search->postream)<<"number_of_total_nodes   = "<<number_of_total_nodes<<"\n";
    //(search->postream)<<"number_of_total_faces   = "<<number_of_total_faces<<"\n";
    //(search->postream)<<"number_of_total_elems   = "<<number_of_total_elements<<"\n";
    //(search->postream).flush();

#ifndef CONTACT_NO_MPI
    if (contact_number_of_processors(communicator) != 1) {
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
      Timer().Start_Timer( base_class_connections );
#endif
      ContactTopologyEntityHash phantom_node_hash( number_of_phantom_nodes, 
                                                   (ContactTopologyEntity<Real>**) phantom_nodes );
      ContactTopologyEntityList *local_nodes = topology->NodeList();        

      if( number_of_phantom_faces ){
        //
        // Connect the nodes to the phantom faces
        //
        for( int i=0 ; i<number_of_phantom_faces ; ++i ){
          ContactTopologyEntity<Real>::connection_data* node_info = phantom_faces[i]->NodeInfo();
          for( int j=0 ; j<phantom_faces[i]->Nodes_Per_Face() ; ++j ){
            ContactHostGlobalID GID( node_info[j].host_gid[0], 
                                     node_info[j].host_gid[1] );

            ContactNode<Real>* face_node = (ContactNode<Real>*)local_nodes->Find(GID);
            if(face_node == NULL) face_node = (ContactNode<Real>*)phantom_node_hash.find(GID);
            POSTCONDITION( face_node );
            phantom_faces[i]->ConnectNode( j, face_node );
          }
        }
        //
        // Connect the faces to the interactions
        //
        ContactTopologyEntityHash phantom_face_hash( number_of_phantom_faces, 
                           (ContactTopologyEntity<Real>**) phantom_faces );
               

        for( ContactNodeFaceInteraction *cnfi = node_entity_list.Node_Face_Iterator_Start();
             cnfi != NULL;
             cnfi = node_entity_list.Node_Face_Iterator_Next()) {
          if( cnfi->FaceEntityData()->owner != my_proc ){
            ContactHostGlobalID gid( cnfi->FaceEntityData()->host_gid[0], 
                                     cnfi->FaceEntityData()->host_gid[1] );
            cnfi->Connect_Face( phantom_face_hash );
            POSTCONDITION(cnfi->Face());
          } else {
            POSTCONDITION(cnfi->Face());
          }
        }
    
#ifdef CONTACT_TD_FACE_FACE_ENF
        for( int i=0 ; i<number_face_face_interactions ; ++i ){
          if( face_face_interaction_list[i]->MasterFaceEntityData()->owner 
                   != my_proc )
            face_face_interaction_list[i]->Connect_MasterFace(phantom_face_hash);
            POSTCONDITION(face_face_interaction_list[i]->MasterFace());
        }
#endif
      }
    
      if( number_of_phantom_elements ){
        // Connect the nodes to the phantom elements
        for( int i=0 ; i<number_of_phantom_elements ; ++i ){
          ContactTopologyEntity<Real>::connection_data* node_info = phantom_elements[i]->NodeInfo();
          for( int j=0 ; j<phantom_elements[i]->Nodes_Per_Element() ; ++j ){
            ContactHostGlobalID GID( node_info[j].host_gid[0], 
                                     node_info[j].host_gid[1] );


            ContactNode<Real>* element_node = (ContactNode<Real>*)local_nodes->Find(GID);
            if(element_node == NULL) element_node = (ContactNode<Real>*)phantom_node_hash.find(GID);
            POSTCONDITION( element_node );
            phantom_elements[i]->ConnectNode( j, element_node );
          }
        }
        // Connect the phantom elements to the interactions
        ContactTopologyEntityHash phantom_element_hash( number_of_phantom_elements, 
                              (ContactTopologyEntity<Real>**) phantom_elements );
        for( int i=0 ; i<number_element_element_interactions ; ++i ){
          if( element_element_interaction_list[i]->
                 MasterElementEntityData()->owner != my_proc )
            element_element_interaction_list[i]->
              Connect_MasterElement( phantom_element_hash );
            POSTCONDITION(element_element_interaction_list[i]->MasterElement());
        }
      }
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
      Timer().Stop_Timer( base_class_connections );
#endif
    
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
      Timer().Start_Timer( base_class_make_asym_comm );
#endif
      Node_AsymComm = new ContactAsymComm( *Node_ZoltanComm, phantom_node_hash, *local_nodes );
      if( topology->Number_of_Face_Blocks() ){


        ContactTopologyEntityHash phantom_face_hash( number_of_phantom_faces, 
                                                     (ContactTopologyEntity<Real>**) phantom_faces );
        ContactTopologyEntityList *local_faces = topology->FaceList();        


        Face_AsymComm = new ContactAsymComm( *Face_ZoltanComm, phantom_face_hash, *local_faces );
        Face_AsymComm->Set_Index_From_EnfArrayIndex();
      } else {
        Face_AsymComm = NULL;
      }
      if( topology->Number_of_Element_Blocks() > 0 ){

        ContactTopologyEntityHash phantom_element_hash( number_of_phantom_elements, 
                                                       (ContactTopologyEntity<Real>**) phantom_elements );
        ContactTopologyEntityList *local_elements = topology->ElemList();        

        Element_AsymComm = new ContactAsymComm(*Element_ZoltanComm, phantom_element_hash, *local_elements);
        Element_AsymComm->Set_Index_From_EnfArrayIndex();
      } else {
        Element_AsymComm = NULL;
      }
      delete Node_ZoltanComm;
      delete Face_ZoltanComm;
      delete Element_ZoltanComm; 
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
      Timer().Stop_Timer( base_class_make_asym_comm );
#endif
      //
      // Now have root proc communicate total shared nodes to all ghosting
      // and phantoming processors so we can build a ContactSymComm for
      // swap-adds.
      // First we need to count how many phantom nodes we really have
      //
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
      Timer().Start_Timer( base_class_make_asym_full_comm );
#endif
      Node_AsymComm_Full = new ContactAsymComm(*(topology->Node_Asym_Comm()), *Node_AsymComm);  
      Node_AsymComm->Set_Index_From_EnfArrayIndex();
      Node_AsymComm_Full->Set_Index_From_EnfArrayIndex();
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
      Timer().Stop_Timer( base_class_make_asym_full_comm );
#endif
    }
#endif
  }
  
#if !defined (CONTACT_NO_MPI) && CONTACT_DEBUG_PRINT_LEVEL>=3
  if (contact_number_of_processors(communicator) > 1) {
    if (Node_AsymComm) {
      (search->postream).flush();
      Node_AsymComm->Print((char*)"ContactEnforcement::Node_AsymComm",search->postream);
    }
    if (Node_AsymComm_Full) {
      (search->postream).flush();
      Node_AsymComm_Full->Print((char*)"ContactEnforcement::Node_AsymComm_Full",search->postream);
    }
    (search->postream).flush();
  }
#endif

#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
  Timer().Stop_Timer( base_class_setup );
#endif

  return error_code;
}

void ContactEnforcement::Clean_Up()
{
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Start_Timer( base_class_cleanup );
#endif
  if( enforcement_node_list ) delete [] enforcement_node_list;
  enforcement_node_list = NULL;
  if( enforcement_face_list ) delete [] enforcement_face_list;
  enforcement_face_list = NULL;
  if( enforcement_element_list ) delete [] enforcement_element_list;
  enforcement_element_list = NULL;
#ifdef CONTACT_TD_FACE_FACE_ENF
  if( face_face_interaction_list ) delete [] face_face_interaction_list;
  face_face_interaction_list = NULL;
#endif
  if( element_element_interaction_list ) delete [] element_element_interaction_list;
  element_element_interaction_list = NULL;

#ifndef CONTACT_NO_MPI
  int i;
  if( Node_AsymComm ){
    delete Node_AsymComm;
    Node_AsymComm = NULL;
  }
  if( Face_AsymComm ){
    delete Face_AsymComm;
    Face_AsymComm = NULL;
  }
  if( Element_AsymComm ){
    delete Element_AsymComm;
    Element_AsymComm = NULL;
  }
  if( Node_AsymComm_Full) {
    delete Node_AsymComm_Full;
    Node_AsymComm_Full = NULL;
  }

  for( i=0 ; i<number_of_phantom_nodes ; ++i ){
    ContactFixedSizeAllocator* alloc=NULL;
    switch( phantom_nodes[i]->Base_Type() ){
    case CT_NODE:
      alloc = &search->Get_Allocator(ContactSearch::ALLOC_ContactNode);
      phantom_nodes[i]->~ContactNode<Real>();
      alloc->Delete_Frag( phantom_nodes[i] );
      break;
    case CT_SHELL_NODE:
      alloc = &search->Get_Allocator(ContactSearch::ALLOC_ContactShellNode);
      static_cast<ContactShellNode*>(phantom_nodes[i])->~ContactShellNode();
      alloc->Delete_Frag( phantom_nodes[i] );
      break;
    default:
      POSTCONDITION(0);
      break;
    }
    POSTCONDITION(alloc!=NULL);
  }
  if( phantom_nodes ){
    delete [] phantom_nodes;
    phantom_nodes = NULL;
  }

  for( i=0 ; i<number_of_phantom_faces ; ++i ){
    ContactFixedSizeAllocator* alloc = ContactSearch::Get_ContactFaceAlloc(phantom_faces[i]->FaceType(), search->Get_Allocators());
    POSTCONDITION(alloc!=NULL);
    phantom_faces[i]->~ContactFace<Real>();
    alloc->Delete_Frag(phantom_faces[i]);
  }
  if( phantom_faces ){
    delete [] phantom_faces;
    phantom_faces = NULL;
  }

  for( i=0 ; i<number_of_phantom_elements ; ++i ){
    ContactFixedSizeAllocator* alloc=NULL;
    switch (phantom_elements[i]->ElementType()) {
    case ContactSearch::CARTESIANHEXELEMENTL8:
      alloc = &search->Get_Allocator(ContactSearch::ALLOC_ContactCartesianHexElementL8 );
      break;
    case ContactSearch::HEXELEMENTL8:
      alloc = &search->Get_Allocator(ContactSearch::ALLOC_ContactHexElementL8);
      break;
    default:
      POSTCONDITION(0);
      break;
    }
    POSTCONDITION(alloc!=NULL);
    phantom_elements[i]->~ContactElement();
    alloc->Delete_Frag(phantom_elements[i]);
  }

  if( phantom_elements ){
    delete [] phantom_elements;
    phantom_elements = NULL;
  }

  node_entity_list.Clear();

#endif
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
    Timer().Stop_Timer( base_class_cleanup );
#endif
}


void ContactEnforcement::swapadd_node_scratch( ScratchVariable &var )
{
#ifndef CONTACT_NO_MPI
  if( contact_number_of_processors( communicator ) != 1 ){
    PRECONDITION( Node_AsymComm_Full );
    contact_reduce_add_scratch_var( communicator,
				    *Node_AsymComm_Full,
				    *(search->comm_buffer),
				    var);
  }
#endif
}

void ContactEnforcement::swapadd_node_scratch_vars( int num_vars,
						    ScratchVariable** vars)
{
#ifndef CONTACT_NO_MPI
  if( contact_number_of_processors( communicator ) != 1 ){
    PRECONDITION( Node_AsymComm_Full );
    contact_reduce_add_scratch_vars( communicator,
			  	     *Node_AsymComm_Full,
				     *(search->comm_buffer),
				     num_vars,
				     vars);
  }
#endif
}

void ContactEnforcement::swapadd_data_array( Real* data_array, 
					     int size_per_entity )
{
#ifndef CONTACT_NO_MPI
  if( contact_number_of_processors( communicator ) != 1 ){
    contact_swapadd_data_array( communicator,
				*(topology->Node_Sym_Comm()),
				*(search->comm_buffer),
				data_array,
				size_per_entity );
  }
#endif
}

void ContactEnforcement::swapadd_data_array( int* data_array, 
					     int size_per_entity )
{
#ifndef CONTACT_NO_MPI
  if( contact_number_of_processors( communicator ) != 1 ){
    contact_swapadd_data_array( communicator,
				*(topology->Node_Sym_Comm()),
				*(search->comm_buffer),
				data_array,
				size_per_entity );
  }
#endif
}

#ifndef CONTACT_NO_MPI
int ContactEnforcement::Insert_Shared_Node_Import_GID( ContactHostGlobalID& gid,
                                                       GIDSortClass *data_to_add,
                                                       const int &data_to_add_length,
                                                       GIDSortClass *existing_data,
                                                       const int &existing_data_length )
{
  //
  //  NKC PERFORMANCE.  Attempting to yeild identical behavior in this function as before without the N^2 complexity factor.  data_to_add
  //  and existing data are sorted arrays.  This routine should find the first instance of gid in the shared_node_import_gid array and return it's
  //  index.  If it does not exist it should add it to the end of the array.
  //
  //  This is done as follows, 
  //   look for item in existing_data, if it is found return the existing data index
  //   look for item in the data_to_add, 
  //
  //       if it is found and the index >= 0 return that index
  //                                     
  //       if it is found and the index < 0 set the index to the last 
  //       array entry and add that value to the array
  //
  //       if it is not found error out
  //

  int existing_binary_lookup = binary_search(gid, existing_data, existing_data_length);
  if(existing_binary_lookup >= 0) {
    int old_index = existing_data[existing_binary_lookup].old_index;
    return old_index;
  }  else {
    int add_binary_lookup = binary_search(gid, data_to_add, data_to_add_length);
    PRECONDITION(add_binary_lookup >= 0);
    int old_index = data_to_add[add_binary_lookup].old_index;
    if(old_index >= 0) {
      return old_index;
    } else {
      int new_index = count_shared_node_import;
      shared_node_import_gid[count_shared_node_import] = gid;
      ++count_shared_node_import;
      data_to_add[add_binary_lookup].old_index = new_index;
      return new_index;
    }
  }
}
#endif


void ContactEnforcement::Map_Array_to_New_Topology( int new_number_nodes,
						    int old_number_nodes,
						    int* old_to_new_map,
						    Real* old_array, 
						    Real* new_array,
						    int array_size,
						    Local_Array_Ordering order )
{
  PRECONDITION( new_number_nodes == 0 || new_array );
  PRECONDITION( old_number_nodes == 0 || old_array );

  int i,new_index;

  switch( order ){
  case( BY_NODE ):
    for( i=0 ; i<old_number_nodes ; ++i ){
      new_index = old_to_new_map[i];
      if( new_index >= 0 )
	std::memcpy( new_array+new_index*array_size, old_array+i*array_size,
		array_size*sizeof(Real) );
    }
    break;
  case( BY_VARIABLE ):
    for( i=0 ; i<old_number_nodes ; ++i ){
      new_index = old_to_new_map[i];
      if( new_index >= 0 ){
	for( int j=0 ; j<array_size ; ++j )
	  new_array[new_index + j*array_size] = old_array[i + j*array_size];
      }
    }
    break;
  }
}


void ContactEnforcement::
Copy_Host_Scalar_Arrays_to_Node_Scratch( int num_arrays,
					 Real** host_arrays,
					 ScratchVariable** scratch_arrays )
{
  ContactShellHandler* shell_handler = 
    search->Get_Primary_Topology()->Shell_Handler();
  if( shell_handler ){
    int num_host_code_nodes = shell_handler->Number_Host_Code_Nodes();
    for( int i=0 ; i<num_host_code_nodes ; ++i ){
      for (int k = 0; k < shell_handler->Num_Acme_Nodes_for_Host_Node(i); ++k) {
        int hi, lo;
        shell_handler->Acme_NodeGID_for_Host_Node(i,k,hi,lo);
        ContactHostGlobalID gid(hi,lo);
        ContactNode<Real>* node = static_cast<ContactNode<Real>*>
          (search->Get_Primary_Topology()->NodeList()->Find(gid));
	POSTCONDITION(node);
        int node_index = node->EnfArrayIndex();
        for( int j=0 ; j<num_arrays ; ++j ){
          *(scratch_arrays[j]->Get_Scratch(node_index)) = host_arrays[j][i];
	}
      }
    }
  } else {
    for( int i=0 ; i<number_of_nodes ; ++i ){
      int k = enforcement_node_list[i]->HostGlobalArrayIndex();
      for( int j=0 ; j<num_arrays ; ++j ){
        *(scratch_arrays[j]->Get_Scratch(i)) = host_arrays[j][k];
      }
    }
  }
}

void ContactEnforcement::Copy_Variable_to_Scratch(int num_arrays, 
                                                  VariableHandle *handles,
                                                  ScratchVariable **scratch_arrays) {
  for( int i=0 ; i<number_of_total_nodes ; ++i ){
    ContactNode<Real> *node = enforcement_node_list[i];
    for( int j=0 ; j<num_arrays ; ++j ){
      Real* variable = node->Variable(handles[j]);
      Real* scratch  = scratch_arrays[j]->Get_Scratch(i);
      int var_size   = scratch_arrays[j]->Get_Size();
      for(int k = 0; k < var_size; ++k) {
        scratch[k] = variable[k];
      }
    }
  }
}

void  ContactEnforcement::
Copy_Node_Vector_Scratch_to_Host_Arrays( int num_arrays,
				         Real** host_arrays,
					 ScratchVariable** scratch_arrays )
{

  // ASG NOTE
  //  There is a massive assumption in this routine for shell nodes.
  //  It's that the values for all ACME shell nodes created from a single 
  //  host code shell node all have the exact same values for the
  //  variable that is being stored for the host code. We only copy the
  //  values from the first ACME shell node for each host code shell node.
  //  We don't explicitly check this assumption here, but that assumption 
  //  is in force. It is something to confirm if shell problems creep up 
  //  again...
  ContactShellHandler* shell_handler = search->Get_Primary_Topology()->Shell_Handler();
  if( shell_handler ){
    int num_host_code_nodes = shell_handler->Number_Host_Code_Nodes();
    for( int i=0 ; i<num_host_code_nodes ; ++i ){
      int hi, lo;
      shell_handler->Acme_NodeGID_for_Host_Node(i,0,hi,lo);
      ContactHostGlobalID gid(hi,lo);
      ContactNode<Real>* node = static_cast<ContactNode<Real>*>
        (search->Get_Primary_Topology()->NodeList()->Find(gid));
      int node_index = node->EnfArrayIndex();
      POSTCONDITION(node);
      for( int j=0 ; j<num_arrays ; ++j ){
        Real* host = host_arrays[j];
        Real* scr =  scratch_arrays[j]->Get_Scratch(node_index);
        host[3*i  ] = scr[0];
        host[3*i+1] = scr[1];
        host[3*i+2] = scr[2];
      }
    }	
  } else {
    for( int i=0 ; i<number_of_nodes ; ++i ){
      int k = enforcement_node_list[i]->HostGlobalArrayIndex();
      for( int j=0 ; j<num_arrays ; ++j ){
	Real* host = host_arrays[j];
	Real* scr = scratch_arrays[j]->Get_Scratch(i);
	host[3*k  ] = scr[0];
	host[3*k+1] = scr[1];
	host[3*k+2] = scr[2];
      }
    }
  }
}

int  ContactEnforcement::
Get_Node_Host_Array_Index(ContactNode<Real> * node)
{
  ContactShellHandler* shell_handler = 
    search->Get_Primary_Topology()->Shell_Handler();
  if( shell_handler ){
    return shell_handler->Acme_to_Host_Node_Map()[node->fcs_index];
  } else {
    return node->HostGlobalArrayIndex();
  }
}



 
