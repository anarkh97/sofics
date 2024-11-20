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


// This base class has multiple functions.  First of all, it prevents
// ContactSearch from needing to be modified when a new enforcement is added
// because this base class is already a 'friend'.  Second, it manages all of
// the scratch memory and importing/swapadding of data and objects for the 
// derived classes.

#ifndef ContactEnforcement_h_
#define ContactEnforcement_h_

#include "ContactSearch.h"
#include "ContactEnforcementData.h"
#include "ContactTimer.h"
#include "ContactTDUserSubTypes.h"
#include "ContactNodeFaceInteraction.h"
#include "NodeEntityInteractionList.h"

class ContactSearch;
class ContactTopology;
template<typename DataType> class ContactNode;
template<typename DataType> class ContactFace;
class ContactElement;
class ContactAsymComm;
class ContactSymComm;
class ContactZoltanComm;
class ContactNodeFaceInteraction;
class ContactNodeSurfaceInteraction;
template<typename DataType> class ContactFaceFaceInteraction;
class ContactElementElementInteraction;
class ContactEnfModel;
class ContactParOStream;
class ScratchVariable;

#ifndef CONTACT_NO_MPI

class ContactEnfZoltan;
class ContactHostGlobalID;

class GIDSortClass;

extern "C" void mergeSort(GIDSortClass *numbers, GIDSortClass *temp, int array_size);
extern "C" void m_sort   (GIDSortClass *numbers, GIDSortClass *temp, int left, int right);
extern "C" void merge    (GIDSortClass *numbers, GIDSortClass *temp, int left, int mid, int right);
extern "C" int  binary_search(const ContactHostGlobalID &key, GIDSortClass *data, const int &data_length);
#endif

class ContactEnforcement {
  
 public:
 
  enum ContactEnforcement_Type{ TDEnforcement=1, GapRemoval, 
                                TiedKinematics, VolumeTransfer, MPCs, 
				Dummy, TDEnfPenalty
#ifdef CONTACT_TD_FACE_FACE_ENF
				, TDFaceFaceEnf
#endif
                              };

  enum Enforcement_Model_Types { TD_FRICTIONLESS=1,     // 1
				 TD_CONSTANT_FRICTION,  // 2
				 TD_TIED,               // 3
				 TD_SPOT_WELD,          // 4
				 TD_PRESSURE_DEPENDENT, // 5
				 TD_VELOCITY_DEPENDENT, // 6
				 TD_SPRING_WELD,        // 7
				 TD_ADHESION,           // 8
				 TD_COHESIVE_ZONE,      // 9
				 TD_JUNCTION,           // 10
				 TD_THREADED,           // 11
				 TD_PV_DEPENDENT,       // 12
				 TD_AREA_WELD,          // 13
				 TD_USER,               // 14
				 TD_SHARED,             // 15
                                 NUMBER_ENFORCEMENT_MODEL_TYPES };
  

  void Register_Enforcement_Timers();

  virtual int Restart_Size();
  virtual ContactSearch::ContactErrorCode Extract_Restart_Data(Real* restart_data );
  ContactSearch* Search() { return search; };
  ContactTopology* Topology() { return topology; };
  int Num_Models() {return num_enforcement_models;};
  ContactEnfModel** Models() {return enforcement_models;};
  
#ifndef CONTACT_NO_MPI
  MPI_Comm Communicator() { return communicator; };
  int Number_Primary_Nodes() { return number_of_nodes; };
  int& Number_Imported_Phantom_Nodes() { return num_imported_phantom_nodes; };
  int& Number_Imported_Phantom_Faces() { return num_imported_phantom_faces; };
  int& Number_Imported_Phantom_Elements() 
    { return num_imported_phantom_elements; };
  ContactNode<Real>** Phantom_Nodes() { return phantom_nodes; };
  ContactFace<Real>** Phantom_Faces() { return phantom_faces; };
  ContactElement** Phantom_Elements() { return phantom_elements; };

  int * Shared_Node_Export_Info() { return shared_node_export_info; };
  int Max_Num_Shared_Procs() 
    { return max_num_shared_procs; };

  int Insert_Shared_Node_Import_GID(ContactHostGlobalID& gid,
                                    GIDSortClass *data_to_add,
                                    const int &data_to_add_length,
                                    GIDSortClass *existing_array,
                                    const int &existing_array_length );

  ContactHostGlobalID * Shared_Node_Import_GID() 
    { return shared_node_import_gid; };
  int * Shared_Node_Import_Info() { return shared_node_import_info; };
#endif

  // Regression Test Functions
  int Number_of_Global_Plot_Variables() 
    { return number_global_plot_vars; };
  int Number_of_Nodal_Plot_Variables() 
    { return number_nodal_plot_vars; };
  int Number_of_Element_Plot_Variables()
    { return number_element_plot_vars; };
  virtual Real Get_Global_Plot_Variable( int ) = 0;
  virtual void Get_Nodal_Plot_Variable( int, Real* ) = 0;
  virtual void Get_Element_Plot_Variable( int, Real* ) = 0;

  // Restart Functions
  virtual int Number_General_Restart_Variables() = 0;
  virtual int Number_Nodal_Restart_Variables() = 0;
  virtual int Number_Edge_Restart_Variables() = 0;
  virtual int Number_Face_Restart_Variables() = 0;
  virtual int Number_Element_Restart_Variables() = 0;
  virtual ContactSearch::ContactErrorCode
    Extract_General_Restart_Variable( Real* data ) = 0;
  virtual ContactSearch::ContactErrorCode
    Extract_Nodal_Restart_Variable( int n, Real* data, int* node_ids ) = 0;
  virtual ContactSearch::ContactErrorCode
    Extract_Edge_Restart_Variable( int n, Real* data ) = 0;
  virtual ContactSearch::ContactErrorCode
    Extract_Face_Restart_Variable( int n, Real* data ) = 0;
  virtual ContactSearch::ContactErrorCode
    Extract_Element_Restart_Variable( int n, Real* data ) = 0;
  virtual ContactSearch::ContactErrorCode
    Implant_General_Restart_Variable( Real* data ) = 0;
  virtual ContactSearch::ContactErrorCode
    Implant_Nodal_Restart_Variable( int n, Real* data ) = 0;
  virtual ContactSearch::ContactErrorCode
    Implant_Edge_Restart_Variable( int n, Real* data ) = 0;
  virtual ContactSearch::ContactErrorCode
    Implant_Face_Restart_Variable( int n, Real* data ) = 0;
  virtual ContactSearch::ContactErrorCode
    Implant_Element_Restart_Variable( int n, Real* data ) = 0;
  virtual ContactSearch::ContactErrorCode Complete_Restart();

  ContactEnforcement_Type Type(){ return type; };
  
  ContactSearch::ContactErrorCode 
    Add_Enforcement_Model( Enforcement_Model_Types Type, int ID, 
			   int* integer_data, Real* real_data, 
                           void* function_pointer=NULL  );

  virtual ContactSearch::ContactErrorCode
    Update_For_Topology_Change( int new_number_of_nodes,
				int* old_to_new_map ) = 0;
  // Errors
  int Number_of_Errors();
  const char* Error_Message( int i );

  // Timing functions
  ContactTimer& Timer() {return search->timer;};
  
  void User_Initialize_Model_Fn(int id, CONTACT_INIT_MODEL_FN * fn);
  void User_Initialize_Time_Step_Fn(int id, CONTACT_INIT_TIME_STEP_FN *fn);
  void User_Initialize_Node_State_Data_Fn(int id, CONTACT_INIT_NODE_STATE_DATA_FN *fn);
  void User_Limit_Force_Fn(int id, CONTACT_LIMIT_FORCE_FN *fn);
  void User_Active_Fn(int id, CONTACT_INTERACTION_ACTIVE_FN *fn);
  void User_Interaction_Type_Fn(int id, CONTACT_INTERACTION_TYPE_FN *fn);

  void Remove_Search_Reference(ContactSearch *search_ref);

#ifndef CONTACT_NO_MPI
  int Get_Count_Shared_Node_Import() {return count_shared_node_import;};
  ContactHostGlobalID *Get_Shared_Node_Import_Gid() {return shared_node_import_gid;};  
#endif

  bool Initialized() { return initialized;};

 protected:
 
  // Scratch Memory Functions
  enum Scratch_Entity_Type{ NODE_SCRATCH = 0,
			    FACE_SCRATCH,
			    ELEMENT_SCRATCH,
			    NFI_SCRATCH,
			    NSI_SCRATCH,
			    NUM_SCRATCH_TYPES };

  ContactEnforcement( ContactSearch::ContactErrorCode&,
                      ContactSearch*, 
                      ContactEnforcement_Type, 
                      int, 
		      const Real*, 
                      bool , 
		      int FM_DATA_INDEX = -1 );
                      
  // Restart Constructor
  ContactEnforcement( ContactSearch::ContactErrorCode&,
                      ContactSearch*, 
                      ContactEnforcement_Type,
		      const Real*,
                      void** function_pointers = NULL );
                      
  virtual ~ContactEnforcement();
  
  // Parallel Enforcement Functions
  ContactSearch::ContactErrorCode Set_Up( bool Update_Enf_Model_State = false );
  void Clean_Up();
  void swapadd_node_scratch( ScratchVariable &var );
  void swapadd_node_scratch_vars( int num_vars, ScratchVariable** );
  void swapadd_data_array( Real* array, int size_per_entity );
  void swapadd_data_array( int* array, int size_per_entity );

  ContactSearch* search;
  ContactTopology* topology;
  ContactNode<Real>** enforcement_node_list;
  ContactFace<Real>** enforcement_face_list;
  ContactElement** enforcement_element_list;
  ContactEnforcementData enforcement_data;

  int dimensionality;
  int number_of_nodes;
  int number_host_code_nodes;
  int number_of_total_nodes;
  int number_entity_keys;
  int number_of_faces;
  int number_of_total_faces;
  int number_of_elements;
  int number_of_total_elements;

  int number_node_entity_constraints;
  NodeEntityInteractionList node_entity_list;

#ifdef CONTACT_TD_FACE_FACE_ENF
  int number_face_face_interactions;
  ContactFaceFaceInteraction<Real>** face_face_interaction_list;
#endif
  int number_element_element_interactions;
  ContactElementElementInteraction** element_element_interaction_list;

  MPI_Comm communicator;

#ifdef CONTACT_DEBUG_NODE
  // Debug node functions
  int Number_Debug_Nodes();
  ContactNode<Real>* Debug_Node(int i);
#endif
  ContactParOStream& ParOStream();
  ContactCommBuffer* CommBuffer();

  int num_enforcement_models;
  ContactEnfModel** enforcement_models;

  ContactErrors* errors;
  char message[81];

  ContactSearch::ContactErrorCode Initialize();

  enum Local_Array_Ordering { BY_NODE, BY_VARIABLE };
  void Map_Array_to_New_Topology( int new_number_nodes, int old_number_nodes,
				  int* old_to_new_map,
				  Real* old_array, Real* new_array, 
				  int array_size, Local_Array_Ordering order );

  void Copy_Host_Scalar_Arrays_to_Node_Scratch( int num_arrays,
						Real** host_arrays,
						ScratchVariable** scratch_arrays );
  void Copy_Node_Vector_Scratch_to_Host_Arrays( int num_arrays,
					        Real** host_arrays,
					        ScratchVariable** scratch_arrays );

  void Copy_Variable_to_Scratch(int num_arrays, 
                                VariableHandle *handles,
                                ScratchVariable **scratch_arrays);               


  // this returns the host array index for arrays that will be passed
  // back to the host codes.
  int Get_Node_Host_Array_Index(ContactNode<Real> * node);

  // Restart Variables
  int Number_Base_General_Restart_Variables() {return 1;};
  void Extract_Base_General_Restart_Variable( Real* data ) 
    {*data = 0.0;};
  void Implant_Base_General_Restart_Variable( Real* data ) 
    {};

  // regression test
  int number_global_plot_vars;
  int number_element_plot_vars;
  int number_nodal_plot_vars;
  Real* nodal_plot_vars;
  Real* element_plot_vars;
  
  // timing variables
  int base_class_setup;
  int base_class_init;
  int base_class_cleanup;
  int base_class_enf_model_update;
  int base_class_connections;
  int base_class_node_import;
  int base_class_face_import;
  int base_class_elem_import;
  int base_class_enf_node_list;
  int base_class_enf_face_list;
  int base_class_enf_elem_list;
  int base_class_make_asym_comm;
  int base_class_make_asym_comm1;
  int base_class_make_asym_comm1a;
  int base_class_make_asym_comm1b;
  int base_class_make_asym_comm2;
  int base_class_make_asym_full_comm;
  
#ifndef CONTACT_NO_MPI
  // Data needed for creating SymComm for swap-adds
  ContactZoltanComm* Node_ZoltanComm;
  ContactZoltanComm* Face_ZoltanComm;
  ContactZoltanComm* Element_ZoltanComm;
  ContactAsymComm* Node_AsymComm;
  ContactAsymComm* Node_AsymComm_Full;
  ContactAsymComm* Face_AsymComm;
  ContactAsymComm* Element_AsymComm;
#endif

 private:

  ContactEnforcement_Type type;

  int allocated_size_node_surface_interactions;
#ifdef CONTACT_TD_FACE_FACE_ENF
  int allocated_size_face_face_interactions;
#endif
  int allocated_size_element_element_interactions;

  // Scratch memory
  int num_scratch_vars[NUM_SCRATCH_TYPES];
  bool scratch_set;

  bool use_enforcement_models;
  int fm_data_index;
  ContactSearch::ContactErrorCode Convert_Enf_Model_ID_to_Index();
  
  bool initialized;

#ifndef CONTACT_NO_MPI
  // Zoltan Interface
  ContactEnfZoltan* zoltan;

  // Phantom topological entities
  int number_of_phantom_nodes;
  int number_of_phantom_faces;
  int number_of_phantom_elements;
  int num_imported_phantom_faces;
  int num_imported_phantom_nodes;
  int num_imported_phantom_elements;
  ContactNode<Real>** phantom_nodes;
  ContactFace<Real>** phantom_faces;
  ContactElement** phantom_elements;

  int count_shared_node_import;
  ContactHostGlobalID * shared_node_import_gid;
  int * shared_node_import_info;

  int * shared_node_export_info;
  int max_num_shared_procs;

#endif

};

#endif
