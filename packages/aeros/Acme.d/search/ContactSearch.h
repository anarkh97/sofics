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


#ifndef ContactSearch_h_
#define ContactSearch_h_

class ContactLinkList;
class ContactErrors; 
template<typename DataType> class ContactNode;
template<typename DataType> class ContactFace;
class ContactElement;
template<typename DataType> class ContactElem;
class ContactNodeFaceInteraction;
class ContactNodeSurfaceInteraction;
class ContactNodeEntityInteraction;
class ContactShellNode;
class ContactTopology;
class ContactSearchData;
class ContactEnforcement;
class ContactTDEnforcement;
class ContactCommBuffer;
class ContactFixedSizeAllocator;
class ContactSequentialAllocator;
class ContactEnforcement;
class ContactTable;
class ContactBoundingBox;
template<typename DataType> class ContactFaceFaceInteraction;
class ObjectBoundingBox;
class DomainBox;

class ObjectBoundingBoxHierarchy;
class ObjectBoundingBoxHierarchy_Int;


#include "Contact_Defines.h"
#include "ContactParOStream.h"
#include "ContactTimer.h"
#include "ContactTopologyEntity.h"
#include "ContactVector.h"
#include "ContactScratchManager.h"
#include <vector>

#ifndef CONTACT_NO_MPI
class ContactCommBuffer;
class ContactZoltan;
class ContactZoltanCommUtils;
#include "mpi.h"
#define MPI_COMPILE 1
#else
#define MPI_COMPILE 0
#endif
char* ACME_Version();
char* ACME_VersionDate();
int ACME_MPI_Compatibility(int);

void Find_Physical_Face(int num_physical_faces, Real** PF_normals, Real * face_normal, int& physical_face, Real& most_opposed);

class ContactSearch {

 /*!
 This is the documentation for the ContactSearch Base Class
 
 */

 public:

  friend class ContactEnforcement;

  enum ContactNode_Type{ NODE=1, POINT, NNODE_TYPES };
  enum ContactNode_Configuration{ CURRENT_CONFIG=1, PREDICTED_CONFIG };
  enum ContactEdge_Type{ NO_EDGES=0, LINEEDGEL2, LINEEDGEQ3, NEDGE_TYPES };
  enum ContactFace_Type{ QUADFACEL4=1, QUADFACEQ8, QUADFACEQ9, 
                         TRIFACEL3, TRIFACEQ6, 
			 SHELLQUADFACEL4, SHELLTRIFACEL3,
			 LINEFACEL2, LINEFACEQ3,
			 NFACE_TYPES };
  enum ContactElem_Type{ HEXELEML8=1, WEDGEELEML6, NELEM_TYPES };
  enum ContactElement_Type{ CARTESIANHEXELEMENTL8=1, HEXELEMENTL8, 
			    NELEMENT_TYPES };

  enum AnalyticSurface_Type{ PLANE=1, SPHERE, CYLINDER_INSIDE, 
			     CYLINDER_OUTSIDE, NSURF_TYPES };
  enum ContactErrorCode{ NO_ERROR = 0,
			 ID_NOT_FOUND,
			 UNKNOWN_TYPE,
			 INVALID_ID,
			 INVALID_DATA,
                         UNIMPLEMENTED_FUNCTION,
			 ZOLTAN_ERROR,
                         EXODUS_ERROR,
                         INVALID_INTERACTION,
                         NULL_FUNCTION_POINTER,
                         INTERNAL_ERROR };

  enum Search_Data_Index{ INTERACTION_TYPE = 0,
			  SEARCH_NORMAL_TOLERANCE,
			  SEARCH_TANGENTIAL_TOLERANCE,
			  NSIZSD };

  enum Search_Interaction_Type{ NO_INTERACTION = 0,
				SLIDING_INTERACTION,
				TIED_INTERACTION,
                                COVERAGE_INTERACTION,
                                GENERIC_INTERACTION,
                                CPP_INTERACTION,
                                INFINITESIMAL_SLIP_INTERACTION,
                                GLUED_INTERACTION };
                                
  enum Node_Block_Attribute{ RADIUS=1, NORMAL };
  enum Face_Block_Attribute{ SHELL_THICKNESS=1, LOFTING_FACTOR };
  enum Search_Option{ MULTIPLE_INTERACTIONS=0, 
                      NORMAL_SMOOTHING,              // 1
		      COMPUTE_NODE_AREAS,            // 2
                      PARTITION_GAP,                 // 3
                      OLD_DYNAMIC_SEARCH,            // 4
                      ENABLE_TRACKING,               // 5
                      NO_SECONDARY,                  // 6
                      GLOBAL_SEARCH_CULL,            // 7
                      NO_WARPED_VOLUME,              // 8
                      NO_PARALLEL_CONSISTENCY,       // 9
                      AUTO_TOL,                      // 10
                      AGGRESSIVE_TOLERANCES,         // 11
                      SKIP_PHYSICAL_FACES,           // 12
		      EDGE_PHYSICAL_FACES,           // 13
                      SHELL_SIMPLE_LOFTING,          // 14
                      NEW_TIED_ENFORCEMENT,          // 15
                      NUM_OPTS,                      // 16
                      COMPUTE_PARTIALS,              // 17
                      NO_GHOSTING };

  enum Search_Option_Status{ INACTIVE=0, ACTIVE };
  enum Smoothing_Resolution{ USE_NODE_NORMAL=0, USE_EDGE_BASED_NORMAL=1 };
  enum PF_Algorithm{ PF_NONE=0, PF_FACE_BASED=1, PF_EDGE_BASED=2 };
  
  enum Search_Cull{ NO_CULL=0, SLAVE_CULL };
  enum Track_Type{ NO_TRACKING=0, LOCAL_TRACKING, GLOBAL_TRACKING };
  
  enum InteractionState { STATE_NONE = -2,
                          STATE_ALL,
                          STATE_0,
                          STATE_1,
                          NUM_STATES };

  enum Topology { PRIMARY = 1,
                  SECONDARY,
                  SEARCH };
		  
  enum SearchType { STATIC1CONFIG = 1,
                    STATIC1CONFIG_TIED,
                    STATIC2CONFIG,
                    DYNAMIC2CONFIG };
                    
  enum State { NO_STATES=-2, ALL_STATES, STATE0, STATE1, NUMBER_OF_STATES };

  // Generic Constructor
  ContactSearch( int  Dimensionality,
		 int  Number_of_States,
		 int  Number_of_Analytic_Surfaces,
		 int  Number_of_Node_Blocks,
	         const ContactNode_Type* Node_Block_Types,
		 const int* Number_Nodes_in_Blocks,
		 const int* Node_Exodus_IDs,
		 const int* Node_Global_IDs,
		 const Real* coords,
		 int  Number_of_Face_Blocks,
		 const ContactFace_Type* Face_Block_Types,
		 const int* Number_Faces_in_Blocks,
		 const int* Face_Global_IDs,
		 const int* Face_Connectivity,
		 const Real* Face_Lofting_Factors,
		 int Number_of_Elem_Blocks,
		 const ContactElement_Type* Element_Block_Types,
		 const int* Number_Elements_in_Blocks,
		 const int* Element_Global_IDs,
		 const int* Element_Connectivity,
		 int  Number_of_Nodal_Comm_Partners,
		 const int* Nodal_Comm_Proc_IDs,
		 const int* Number_Nodes_to_Partner,
		 const int* Communication_Nodes,
#ifndef CONTACT_NO_MPI
		 MPI_Comm& mpi_communicator,
#else
		 int& mpi_communicator,
#endif
		 ContactErrorCode& error );


  ~ContactSearch();


  static Real ComputeCurvatureFromAngle(const Real angle);

  void UpdateSearch( const int* num_node_deaths_per_block, 
		     const int* node_deaths_global_ids,
                     const int* num_face_deaths_per_block, 
		     const int* face_deaths_global_ids,
                     const int* num_element_deaths_per_block, 
		     const int* element_deaths_global_ids,
                     
                     const int* num_node_births_per_block, 
		     const int* node_births_exodus_ids,
		     const int* node_births_global_ids,
		     const int* number_face_births_per_block, 
		     const int* face_births_global_ids,
		     const int* face_births_connectivity,
		     const int* number_element_births_per_block, 
		     const int* element_births_global_ids,
		     const int* element_births_connectivity,
       
		     const int  num_node_exports,
                     const int* node_export_id_list,
                     const int* node_export_pid,
		     const int  num_face_exports,
                     const int* face_export_id_list,
                     const int* face_export_pid,
		     const int  num_element_exports,
                     const int* element_export_id_list,
                     const int* element_export_pid,
                                 
                     const int num_nodes,
                     const int* node_host_ids,
                     const int num_faces,
                     const int* face_host_ids,
                     const int num_elements,
                     const int* element_host_ids,
                     
		     const int  num_comm_partners,
		     const int* comm_proc_id,
		     const int* number_nodes_to_partner,
		     const int* comm_node,
                     
		     ContactErrorCode& error );

  // Restart functions 
  ContactSearch( Real* restart_data,
	  	 const int* Node_Host_IDs,
	  	 const int* Face_Host_IDs,
	  	 const int* Element_Host_IDs,
#ifndef CONTACT_NO_MPI
		 MPI_Comm& mpi_communicator,
#else
		 int& mpi_communicator,
#endif
		 ContactErrorCode& error );
  int Restart_Size();
  ContactErrorCode Extract_Restart_Data( Real* buffer );

  ContactCommBuffer* Get_Comm_Buffer() const {return comm_buffer;}

  //
  //  Determine if the input face is active for searching in the current time step
  //
  bool Active_Search_Status(ContactFace<Real> *face);



  int Number_General_Restart_Variables();
  int Number_Nodal_Restart_Variables();
  int Number_Edge_Restart_Variables();
  int Number_Face_Restart_Variables();
  int Number_Element_Restart_Variables();
  ContactErrorCode Extract_General_Restart_Variable( Real* data );
  ContactErrorCode Implant_General_Restart_Variable( Real* data );
  ContactErrorCode Extract_Nodal_Restart_Variable( int n, Real* data, int* node_id );
  ContactErrorCode Implant_Nodal_Restart_Variable( int n, Real* data );
  ContactErrorCode Extract_Edge_Restart_Variable( int n, Real* data );
  ContactErrorCode Implant_Edge_Restart_Variable( int n, Real* data );
  ContactErrorCode Extract_Face_Restart_Variable( int n, Real* data );
  ContactErrorCode Implant_Face_Restart_Variable( int n, Real* data );
  ContactErrorCode Extract_Element_Restart_Variable( int n, Real* data );
  ContactErrorCode Implant_Element_Restart_Variable( int n, Real* data );
  ContactErrorCode Complete_Restart();
  //
  //  Remove all tracking data and the like.  This will cause the search to revert to its origional state.
  //  Call this after the various search initializion calls to find overlap or initial tied interactions.
  //
  void Reset_Search();

  void Set_Initial_Gap();

  ContactErrorCode
    Set_Node_Block_Remaining_Gap( int node_blk_id,
				  const Real* gap);

  ContactErrorCode
    Set_Node_Block_Ghosting_Gap( int node_blk_id,
				 const Real* gap);

  ContactErrorCode
    Set_Node_Block_Kinematic_Constraints( int node_blk_id,
					  const int* constraints_per_node,
					  const Real* constraint_vector);

  ContactErrorCode 
    Set_Node_Block_Configuration(ContactNode_Configuration config,
				 int node_blk_id, 
				 const Real* positions );
  
  ContactErrorCode Set_Node_Block_Attributes(Node_Block_Attribute attr,
					     int node_blk_id, 
					     const Real* attributes );
  ContactErrorCode Set_Face_Block_Attributes(Face_Block_Attribute attr,
					     int face_blk_id,
					     const Real* attributes );

  ContactErrorCode Add_Analytic_Surface(AnalyticSurface_Type surface_type, 
					const Real* data);
  ContactErrorCode Set_Analytic_Surface_Configuration( int id,
						       const Real* data );
  
  // Search_Data
  ContactErrorCode Check_Search_Data_Size( int size_data_per_pair,
					   int number_of_entity_keys );

  void Set_Search_Data( const Real* );
  void Set_Search_Data( Search_Data_Index, int, int, Real );

  ContactErrorCode Set_Search_Option( Search_Option, 
                                      Search_Option_Status, 
				      Real* );

  void Get_Search_Option( Search_Option, Search_Option_Status&, Real* );
  
  int Get_Tracking_Status() {return enable_tracking;};
  
  int OffFaceTrackingStatus() {return enable_off_face_tracking;};
  PF_Algorithm PhysicalFaceAlgorithm() {return physical_face_algorithm;};

  ContactSearchData* Search_Data( ) {return search_data;};

  // Search Methods
  void Global_NodeNodeSearch();
  void Global_NodeFaceSearch();
  void Global_FaceFaceSearch();
  void Global_ElemElemSearch();
  
  //!
  //! This is the static_1_configuration search that is used by Calore
  //!
  ContactErrorCode Static_Search_1_Configuration();
  //!
  //! This is the static_2_configuration search that is used by Adagio
  //!
  ContactErrorCode Static_Search_2_Configuration();
  //!
  //! This is the dynamic_2_configuration search that is used by Presto/Alegra
  //!
  ContactErrorCode Dynamic_Search_2_Configuration(Real& dt_old,
                                                  Real& dt );
  ContactErrorCode Dynamic_Tracking_Search_2_Configuration();
  
  ContactErrorCode Static_Search_1_Configuration_Tied();

  void TopLevel_Search( SearchType search_type, Real& dt_old, Real& dt );
  
  void NewGlobalSearch(SearchType search_type, 
                       Real gap_tol, int num_configs,
                       VariableHandle CURRENT_POSITION, 
                       VariableHandle PREDICTED_POSITION,
                       VariableHandle AUGMENTED_POSITION);
  void GlobalSearch(SearchType search_type, 
                    Real gap_tol, int num_configs,
                    VariableHandle CURRENT_POSITION, 
                    VariableHandle PREDICTED_POSITION,
                    VariableHandle AUGMENTED_POSITION);
  void TrackedSearch(SearchType search_type, 
                     Real gap_tol, int num_configs,
                     VariableHandle CURRENT_POSITION, 
                     VariableHandle PREDICTED_POSITION,
                     VariableHandle AUGMENTED_POSITION);
  void Global_NodeNodeSearch(SearchType search_type, int num_configs,
                             VariableHandle CURRENT_POSITION, 
                             VariableHandle PREDICTED_POSITION,
                             VariableHandle AUGMENTED_POSITION);
  void Global_NodeFaceSearch(SearchType search_type, 
                             Real gap_tol, int num_configs,
                             VariableHandle CURRENT_POSITION, 
                             VariableHandle PREDICTED_POSITION,
                             VariableHandle AUGMENTED_POSITION);
  void Global_FaceFaceSearch(SearchType search_type, int num_configs,
                             VariableHandle CURRENT_POSITION, 
                             VariableHandle PREDICTED_POSITION,
                             VariableHandle AUGMENTED_POSITION);
  void Global_ElemElemSearch(SearchType search_type, int num_configs,
                             VariableHandle CURRENT_POSITION, 
                             VariableHandle PREDICTED_POSITION,
                             VariableHandle AUGMENTED_POSITION);
  void Print_Search_Summary();
  void Print_Search_Header(SearchType search_type, Real dt, Real dt_old, ContactTDEnforcement* td_enf);
  void Print_Search_Footer();
  
  ContactErrorCode Delete_All_Interactions(); 

  //
  //  Search helper functions
  //
  void build_nodeface_search_status( int *node_search_status,
                                     int *face_search_status );
  void build_nodeface_search_status_proximity( int *node_search_status,
                                               int *face_search_status );

  void process_node_face_interactions(const int &list_size,
                                      ACME::ContactNode_Vector &list,
                                      ContactFace<Real> *face, 
                                      const VariableHandle &CALCULATION_POSITION,
                                      ContactNode<Real> **Nodes,
                                      const VariableHandle &NODE_NORMAL,
                                      Real &CTOLT1,
                                      Real &CTOLT2,
                                      Real &mag_tol,
                                      int *physical_faces,
                                      ACME::Int_Vector &node_keys,
                                      int interaction_source);


  void Get_Scratch(const int number_of_nodes);
  void Clear_Scratch();

  void Process_Analytic_Surfaces(Real *position,
                                 int *index,
                                 int number_of_nodes,
                                 Real *scratch,
                                 int *rank2,
                                 int *list,
                                 ContactNode<Real>** Nodes, int* node_map,
                                 int interaction_source,
                                 VariableHandle position_var,
                                 VariableHandle position_var_2 = -1);

  void Process_Analytic_Surfaces(
                                 ObjectBoundingBoxHierarchy *node_hierarchy,
                                 ContactNode<Real>** Nodes,
                                 int interaction_source,
                                 VariableHandle position_var,
                                 VariableHandle position_var_2 = -1);


  void Process_Analytic_Surfaces_New(
                                 ObjectBoundingBoxHierarchy **node_hierarchy_ptrs,
                                 ContactNode<Real> **node_block_entities,
                                 int interaction_source,
                                 VariableHandle position_var,
                                 VariableHandle position_var_2 = -1);

  // Errors
  int Number_of_Errors();
  const char* Error_Message( int i );
  
  // Extracting Node-Node Interactions 
  void Size_NodeNode_Interactions( int& num_interactions, int& data_size );
  void Get_NodeNode_Interactions( int* slave_node_block_ids,
				  int* slave_node_indexes_in_block,
				  int* master_node_block_ids,
				  int* master_node_indexex_in_block,
                                  int* master_node_proc, 
				  Real* data );

  // Extracting Node-Face Interactions 
  void Size_NodeFace_Interactions( int& num_interactions, int& data_size );
  void Get_NodeFace_Interactions( int* node_block_ids,
				  int* node_indexes_in_block,
				  int* node_entity_keys,
				  int* face_block_ids,
				  int* face_indexex_in_block,
                                  int* face_proc, 
				  Real* data );

  // Extracting Node-Surface Interactions
  void Size_NodeSurface_Interactions( int& num_interactions, int& data_size );
  void Get_NodeSurface_Interactions( int* node_block_ids,
				     int* node_indexes_in_block,
				     int* analyticsurface_ids,
				     Real* data );

  // Extracting Face-Face Interactions 
  void Size_FaceFace_Interactions( int& num_interactions, int& data_size );
  void Get_FaceFace_Interactions( int* slave_face_block_ids,
				  int* slave_face_indexes_in_block,
                                  int* slave_face_proc,
				  int* master_face_block_ids,
				  int* master_face_indexex_in_block,
                                  int* master_face_proc, 
                                  int* interaction_index,
				  Real* interaction_data );

  // Extracting Face-Coverage Interactions 
  void Size_FaceCoverage_Interactions( int& num_interactions, int& data_size );
  void Get_FaceCoverage_Interactions( int* face_block_ids,
				      int* face_indexes_in_block,
                                      int* interaction_index,
				      Real* interaction_data );

  // Extracting Element-Element Interactions 
  void Size_ElementElement_Interactions( int& num_interactions, int& data_size );
  void Get_ElementElement_Interactions( int* slave_element_block_ids,
				        int* slave_element_indexes_in_block,
				        int* master_element_block_ids,
				        int* master_element_indexex_in_block,
                                        int* master_element_proc, 
				        Real* interaction_data );

  ContactErrorCode Exodus_Output( int Exodus_ID, Real Time );

  inline ContactTopology* Get_Primary_Topology() {return primary_topology;};
  
  inline ContactTopology* Get_Secondary_Topology() {return secondary_topology;};

#ifndef CONTACT_NO_MPI
  inline MPI_Comm Get_Comm() {return SearchComm;};
  inline ContactZoltan* Get_Zoltan() {return zoltan;};
#else
  inline int Get_Comm() {return SearchComm;};
#endif
  
  ContactParOStream& ParOStream() {return postream;};

  ContactErrorCode Add_Debug_Node( int );

  void Display();

  // allocators
  enum AllocatorType { ALLOC_ContactNode = 0,
		       ALLOC_ContactLineEdgeL2,
		       ALLOC_ContactLineEdgeQ3,
		       ALLOC_ContactLineFaceL2,
		       ALLOC_ContactLineFaceQ3,
                       ALLOC_ContactQuadFaceL4,
		       ALLOC_ContactQuadFaceQ8,
		       ALLOC_ContactQuadFaceQ9,
		       ALLOC_ContactTriFaceL3,
		       ALLOC_ContactTriFaceQ6,
		       ALLOC_ContactHexElemL8,
		       ALLOC_ContactWedgeElemL6,
		       ALLOC_ContactNodeNodeInteraction,
		       ALLOC_ContactNodeFaceInteraction,
		       ALLOC_ContactNodeSurfaceInteraction,
		       ALLOC_ContactFaceFaceInteraction,
		       ALLOC_ContactFaceCoverageInteraction,
		       ALLOC_ContactElementElementInteraction,
		       ALLOC_ContactPolyVert,
		       ALLOC_ContactPolyEdge,
		       ALLOC_ContactPoly,
		       ALLOC_ContactFaceFaceGraphNode,
		       ALLOC_ContactCartesianHexElementL8,
		       ALLOC_ContactHexElementL8,
		       ALLOC_ContactShellQuadFaceL4,
		       ALLOC_ContactShellTriFaceL3,
		       ALLOC_ContactShellNode,
                       ALLOC_NUM_ALLOCATED_ENTITIES };
  ContactFixedSizeAllocator* Get_Allocators( ){return allocators;};
  ContactFixedSizeAllocator& Get_Allocator( AllocatorType );
  ContactSequentialAllocator& Get_Scratch_Allocator()
    {return *scratch_allocator;};

  ContactErrorCode Add_Table( int ID, int Num_Points, Real* abscissa,
			      Real* ordinate );
  const ContactTable* Get_Table( int ID );
  int Num_Tables() { return num_tables; };
  ContactTable** Tables() { return tables; };

  static bool Is_a_Shell_Face(ContactFace_Type face_type) {
    if ( face_type == SHELLQUADFACEL4 ||
	 face_type == SHELLTRIFACEL3) return true;
    return false;
  }

  static int                        Number_Nodes_Per_Face       (ContactFace_Type face_type);
  static ContactEdge_Type           Face_Edge_Type              (ContactFace_Type face_type);
  static ContactFace<Real>*         New_ContactFace             (ContactFace_Type type,
                                                                 ContactFixedSizeAllocator *alloc);
  static ContactFixedSizeAllocator* Get_ContactFaceAlloc        (ContactFace_Type type,
                                                                 ContactFixedSizeAllocator *alloc);
  static void                       Delete_ContactTopologyEntity(ContactTopologyEntity<Real>*,
                                                                 ContactFixedSizeAllocator *alloc);

  inline int Max_Interactions_Per_Node() { return 3; };
  
  ContactTimer* Timer() {return &timer;};
  inline int Contact_time() {return contact_time;};
  inline int Owner_help_migrate_size_time() {return owner_help_migrate_size_time;};
  inline int Owner_help_migrate_pack_time() {return owner_help_migrate_pack_time;};
  inline int Owner_help_migrate_unpack_time() {return owner_help_migrate_unpack_time;};
  inline int Interaction_help_migrate_size_time() {return interaction_help_migrate_size_time;};
  inline int Interaction_help_migrate_pack_time() {return interaction_help_migrate_pack_time;};
  inline int Interaction_help_migrate_unpack_time() {return interaction_help_migrate_unpack_time;};

//  void Partition_Between_GapCur_and_GapOld(ContactNodeFaceInteraction* cnfi);
//  void Partition_Between_GapCur_and_GapOld(ContactNodeSurfaceInteraction* cnsi);
  void Partition_Between_GapCur_and_GapOld(ContactNodeEntityInteraction* cnei);

  int Do_NodeNode_Search() { return do_node_node_search;};
  int Do_NodeFace_Search() { return do_node_face_search;};
  int Do_FaceFace_Search() { return do_face_face_search;};
  int Do_Coverage_Search() { return do_coverage_search; };
  int Do_ElemElem_Search() { return do_elem_elem_search;};
  void Create_Processor_Bounding_Boxes(ObjectBoundingBox *proc_box_array, 
                                       int &POSITION, 
                                       const int &total_num_procs);
  int contact_time;
  int search_time;
  int create_secondary_time;
  int cleanup_secondary_time;
  int track_search_time;
  int global_search_time;
  int search_remaining_gap_time;
  int rcb_time;
  int secondary_owner_migration_time;
  int owner_help_migrate_time;
  int owner_help_migrate_size_time;
  int owner_help_migrate_pack_time;
  int owner_help_migrate_unpack_time;
  int secondary_connect_mesh_time;
  int secondary_connect_mesh_phase2_time;
  int secondary_connect_mesh_phase3_time;
  int secondary_connect_mesh_phase4_time;
  int secondary_connect_mesh_phase5_time;

  int interaction_migration_time;
  int interaction_onproc_copy_time;
  int interaction_export_setup_time;
  int interaction_help_migrate_time;
  int interaction_help_migrate_size_time;
  int interaction_help_migrate_pack_time;
  int interaction_help_migrate_unpack_time;
  int baseline_constructors_time;
  int baseline_node_time;
  int baseline_edge_time;
  int baseline_face_time;
  int baseline_size_time;
  int baseline_size_all_time;
  int baseline_pack_time;
  int baseline_pack_all_time;
  int baseline_unpack_time;
  int baseline_unpack_all_time;
  int on_processor_search_time;
  int restore_constraint_continuity_time;
  int augmented_config_time;
  int search_setup_time;
  int search_update_state_time;
  int search_remianing_gap_time;
  int geom_time;
  int context_time;  

  int elem_geom_time;
  int face_charlen_geom_time;
  int face_normal_geom_time;
  int node_normal_geom_time;
  int analytic_surf_geom_time;
  int edge_curve_geom_time;
  int edge_curve_geom_time_serial;
  int edge_smooth_geom_time;
  int shell_loft_geom_time;
  
  int track_ghost_time;
  int track_scratch_time;
  int track_gap_time;
  int track_phys_face_time;
  int track_retr_tracked_time;
  int track_serial_setup_time;
  int track_main_time;
  int track_id_time;
  int track_release_scratch_time;
  int track_cleanup_time;
  int search_scratch_time;
  int search_gap_time;
  int search_phys_face_time;
  int search_retr_tied_time;
  int search_retr_tracked_time;
  int search_serial_setup_time;
  int search_main_time;
  int search_sub_time_1;
  int search_sub_time_2;
  int search_sub_time_3;
  int search_id_time;
  int search_release_scratch_time;
  int delete_ghosting;
  
  void Sort(int, int*);
  
  void StepNumber(int n) {step_number=n;};
  int  StepNumber() {return step_number;};
  void TrackStep(int n) {tracking_step=n;};
  int  TrackStep() {return tracking_step;};

  void Compute_Nodes_in_Box(int& ne,
                            Real* const xmin,
                            Real* const xmax,
                            Real* position_search,
                            int* index,
                            int& number_of_nodes,
                            Real* scratch,
                            int* rank2,
                            int& list_size,
                            int* list);

  void Compute_Remaining_Gap( Real*, Real* ,
                     ContactSearch::Topology, 
                     ContactTopologyEntity<Real>::SearchContext);

  void Compute_Remaining_Gap( Real*, Real* );
  
  Real* MaxNodeMotion()   {return max_node_motion;};
  Real* MaxRemainingGap() {return max_remaining_gaps;};
  
  void Bytes_For_Nodes(Real val) {bytes_nodes=val;};
  void Bytes_For_Edges(Real val) {bytes_edges=val;};
  void Bytes_For_Faces(Real val) {bytes_faces=val;};
  void Bytes_For_Elems(Real val) {bytes_elems=val;};
  
  Real Bytes_For_Nodes() {return bytes_nodes;};
  Real Bytes_For_Edges() {return bytes_edges;};
  Real Bytes_For_Faces() {return bytes_faces;};
  Real Bytes_For_Elems() {return bytes_elems;};

  Real BoxInflation() { return box_inflation; };
  int  AutoTol()      { return (int)auto_tol; };

  void append_invalid_edge(const int n1, const int n2) {invalid_edges.push_back(std::pair<int,int>(n1,n2));};  
  std::vector< std::pair<int, int> >& getInvalidEdges(){return invalid_edges;};

  void SetInitialContext();
  void SetContextForSearchSlaves();
  void ResetContextForSearchSlaves();
  void SetContextForGhostingRCB();
  void SetContextForGlobalSearchGhosting();
  void SetContextForTrackSearchGhosting();
  void SetContextForGeometryUpdate();
  
  Real CaptureMotion() {return capture_motion;};

  int NumConfigs() {return nconfigs;};

  template<typename DataType>
    ContactFaceFaceInteraction<DataType>* Face_Face_Search( ContactFace<DataType>*, ContactFace<DataType>*,
                                                            ContactElem<DataType>*, VariableHandle,
                                                            ContactFixedSizeAllocator*);
  template<typename DataType>
    ContactFaceFaceInteraction<DataType>* Partial_Face_Face_Search( ContactFace<DataType>*, ContactFace<DataType>*,
                                                                    ContactElem<DataType>*, VariableHandle,
                                                                    Real tol, ContactFixedSizeAllocator*);

  template<typename DataType>
    ContactFaceFaceInteraction<DataType>* Second_Partial_Face_Face_Search( ContactFace<DataType>*, ContactFace<DataType>*,
                                                                           ContactElem<DataType>*, VariableHandle,
                                                                           Real tol, ContactFixedSizeAllocator*);

 protected:

  void Register_Enforcement( ContactEnforcement* );
  void Delete_Enforcement( ContactEnforcement* );
  int number_registered_enforcements;
  ContactEnforcement** enforcement;
    
 private:
  
  bool initialized;
  bool initialized_tied;
  bool initialized_context;
  bool initializing_tied;
  
  int step_number;

  int do_node_node_search;
  int do_node_face_search;
  int do_face_face_search;
  int do_coverage_search;
  int do_elem_elem_search;

  
  // Used for timings
  void Register_Timers();

  // A parallel std::ostream object to manage standard out
  ContactParOStream postream;

  // A parallel timer
  ContactTimer timer;

  // A ContactSearch object is not copyable or assignable
  ContactSearch( const ContactSearch& );
  ContactSearch& operator=(const ContactSearch& );
  

  // Secondary Decomposition
  void Create_Search_Topology( int POSITION );
  ContactErrorCode Define_Primary_Interactions();
#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS) && defined(CONTACT_TIMINGS1)
  void BaselineTiming();
#endif

  int dimensionality;
  //int num_entity_keys;
  int num_states;
#ifndef CONTACT_NO_MPI
  MPI_Comm SearchComm;
  ContactZoltan* zoltan;
  void create_zoltan_object(ContactErrorCode&);
#else
  int SearchComm;
#endif

  // Topology Information
  ContactTopology* primary_topology;
  ContactTopology* secondary_topology;
  ContactTopology* search_topology;
  ContactSearchData* search_data;
  ContactErrors* errors;
  ContactCommBuffer* comm_buffer;
  char message[81];

  // Options
  Search_Option_Status multiple_interaction_status;
  Search_Option_Status normal_smoothing_status;
  Search_Option_Status compute_node_areas;
  Search_Option_Status partition_gap_status;
  Search_Option_Status old_dynamic_search;
  Search_Option_Status enable_tracking;
  Search_Option_Status no_secondary;
  Search_Option_Status no_ghosting;
  Search_Option_Status keep_ghosting;
  Search_Option_Status search_cull;
  Search_Option_Status no_warped_volume;
  Search_Option_Status no_parallel_consistency;
  Search_Option_Status auto_tol;
  Search_Option_Status aggressive_tolerances;
  Search_Option_Status skip_physical_faces;
  Search_Option_Status edge_physical_faces;
  Search_Option_Status shell_simple_lofting;
  Search_Option_Status new_tied_enforcement;
  Search_Option_Status compute_partials;

  Real orig_sharp_smooth_angle;
  Real sharp_smooth_curvature;
  Real normal_smoothing_distance;
  Smoothing_Resolution smoothing_resolution;
  Search_Cull global_search_cull;
  Track_Type tracking_type;
  PF_Algorithm physical_face_algorithm;
  int enable_off_face_tracking;
  int computed_partials_order;

  Real max_node_motion[3];
  Real max_node_displacement;
  Real max_remaining_gaps[3];
  Real max_remaining_gap_mag;
  Real reasonable_gap;
  Real box_inflation, gap_inflation;

  ContactErrorCode Interaction_Definition(int num_configs, 
                     ContactSearch::Topology use_topology, 
                     ContactTopologyEntity<Real>::SearchContext status);

  void Retrieve_Glued_Interactions(VariableHandle, ContactSearch::Topology, int flag=0);
//  void Retrieve_Tied_Interactions (VariableHandle, ContactSearch::Topology, ContactTopologyEntity<Real>::SearchContext);
  void Retrieve_Tied_Interactions_From_Primary(VariableHandle);
  void Retrieve_Tracked_Interactions(ContactSearch::SearchType, Real);
  void Initialize_Context();
  void Initialize_Search();
  void Update_Interaction( ContactNode<Real>*, ContactFace<Real>*, int, int, int );
  void Process_Tracking_Face( ContactFace<Real>*, ContactNode<Real>*, int, int, int );
  int Cull_Node_List( int, ContactNode<Real>*, Real*);
  void Process_Interaction( VariableHandle, ContactNodeEntityInteraction*, int& );
  //void Process_Interaction( VariableHandle, ContactNodeEntityInteraction*, int );
  
  void ComputeInteractions(const int list_size, 
                           ACME::ContactNode_Vector &node_list,
                           ContactFace<Real> *face,
                           const VariableHandle &CURRENT_POSITION,
                           const VariableHandle &AUGMENTED_POSITION,
                           const VariableHandle &PREDICTED_POSITION,
                           const VariableHandle &NODE_NORMAL,
                           Real &user_search_tol,
                           Real &CTOLT1,
                           Real &CTOLT2,
                           Real &mag_tol,
                           int *physical_faces,
                           const ACME::Int_Vector &node_keys);

  int CheckForNodeFaceInteraction(ContactNode<Real>*, ContactFace<Real>*,
                                  ContactSearch::SearchType, 
                                  ContactNodeFaceInteraction*, 
                                  int, int, Real);

  int Cull_Node_List_New( ACME::ContactNode_Vector &node_list, int*, ContactFace<Real>*,
                          int, int, int,
                          const ContactBoundingBox &current_face_box,
                          const ContactBoundingBox &predicted_face_box,
                          ACME::Int_Vector &node_keys,
                          const std::vector<Real> &norm_search_tol,
                          const std::vector<Real> &tang_search_tol,
                          const std::vector<bool> &valid_inter,
                          int, int);

  int Cull_Node_List_New1( ACME::ContactNode_Vector &node_list, ContactFace<Real>*,
                           int, int, int,
                           const ContactBoundingBox &current_face_box,
                           const ContactBoundingBox &predicted_face_box,
                           const std::vector<bool> &valid_inter,
                           int, int);

  int Cull_Node_List( int, int*, ContactNode<Real>*, Real*);
  int Dynamic_Process_Method( ContactNode<Real>*, ContactFace<Real>* );
  ContactNodeEntityInteraction* Compete_Interactions(VariableHandle NODE_COORDS,
						   ContactNodeEntityInteraction*,
						   ContactNodeEntityInteraction*);
  ContactNodeFaceInteraction* 
    Compete_Interactions_Connected(ContactNodeFaceInteraction*,
				   ContactNodeFaceInteraction*);
  ContactNodeEntityInteraction* 
    Compete_Interactions_Unconnected(VariableHandle NODE_COORDS,
				     ContactNodeEntityInteraction*,
				     ContactNodeEntityInteraction*);

  void Build_Physical_Face_List_None(ContactSearch::Topology use_topology, 
                                     ContactTopologyEntity<Real>::SearchContext status,
                                     bool use_proximity);
  void Build_Physical_Face_List_FaceWalk(ContactSearch::Topology use_topology, 
                                         ContactTopologyEntity<Real>::SearchContext status,
                                         bool use_proximity);
  void Build_Physical_Face_List_EdgeWalk(ContactSearch::Topology use_topology, 
                                         ContactTopologyEntity<Real>::SearchContext status,
                                         bool use_proximity);
                                         
  void Process_Face_Coverage( void );

  ContactFixedSizeAllocator* allocators;  // array
#if (MAX_FFI_DERIVATIVES > 0)
  ContactFixedSizeAllocator* active_allocators;  // array
#endif
  ContactSequentialAllocator* scratch_allocator;
  void Set_Up_Allocators();
  
  void ResetInteractionHostIDs( );

  ContactErrorCode error_code;

  int num_tables;
  ContactTable** tables;
  
  // Physical Face scratch variable handles
  ScratchVariable NUMBER_PHYSICAL_FACES;
  ScratchVariable PHYSICAL_FACE_LIST_INDEX;
  ScratchVariable PHYSICAL_FACE_NORMAL_1;
  ScratchVariable PHYSICAL_FACE_NORMAL_2;
  ScratchVariable PHYSICAL_FACE_NORMAL_3;

  int num_active_faces;
  int num_active_nodes;
  int num_global_nodes;
  int num_tracked_nodes;
  int num_tracked_interactions;
  
  Real bytes_nodes;
  Real bytes_edges;
  Real bytes_faces;
  Real bytes_elems;
  
  // search sratch memory
  int   ws_size;
  int   max_facets;
  int   data_size;
  Real* ms_coordinates_c;
  Real* ms_coordinates_a;
  Real* ms_coordinates_p;
  Real* ms_normals_c;
  Real* ms_normals_a;
  Real* ms_normals_p;
  Real* ctrl;
  Real* ctrcl;
  Real* ctrcl_facets;
  int*  pushback_dir_flag;
  
  int   global_tracking_interval;
  int   tracking_step;
  bool  restart;
  
  int nconfigs;

  //
  //  Invalid mesh objects to return to host code
  //
  std::vector< std::pair<int, int> > invalid_edges;
  
  Real capture_motion;

 public:
  /**
   * Set_Zoltan_Param passes parameter options (character strings)
   * from the host code to acme and down to zoltan.
   */
  void Set_Zoltan_Param(char *zoltan_parameter, char *zoltan_value);

  /**
   * Set_Zoltan_Method passes method options (character strings)
   * from the host code to acme and down to zoltan.
   */
  void Set_Zoltan_Method(char *zoltan_method);
};

#endif // #ifndef ContactSearch_h_
