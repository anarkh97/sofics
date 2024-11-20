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


#ifndef ContactTopology_h_
#define ContactTopology_h_

#include "ContactEntity.h"
#include "ContactNodeBlock.h"
#include "ContactElementBlock.h"
#include "ContactFace.h"
#include "Contact_Defines.h"
#include "ContactAnalyticSurface.h"
#include "ContactSearch.h"
#include "ContactCommBuffer.h"
#include "ContactParOStream.h"
#ifndef CONTACT_NO_MPI
#include "zoltan.h"
#include "ContactZoltan.h"
#endif
#include "ContactTopologyEntityList.h"
#include "ContactBlockEntityList.h"

class ContactNodeBlock;
class ContactEdgeBlock;
class CString;
class ContactFaceBlock;
class ContactErrors;
class ContactEntityHash;
class ContactParOStream;
#ifndef CONTACT_NO_MPI
class ContactSymComm;
class ContactAsymComm;
class ContactZoltanComm;
class ContactZoltanCommUtils;
class ContactZoltan;
#endif
class ContactEnforcement;
class ContactShellHandler;

class ContactTopology {

 public:

  friend class ContactEnforcement;
  friend class ContactShellHandler;
  
  // Because of the use of this, do NOT remove setting the first entry to -1
  // The enum value, corresponds to the index in the array Var_Handles
  enum VariableHandle_enum{ UNKNOWN_VARIABLE = -1,
#define    NODE_SCALAR_VAR( a,b ) a,
#define    NODE_VECTOR_VAR( a,b ) a,
#define    EDGE_SCALAR_VAR( a,b ) a,
#define    EDGE_VECTOR_VAR( a,b ) a,
#define    FACE_SCALAR_VAR( a,b ) a,
#define    FACE_VECTOR_VAR( a,b ) a,
#define    ELEM_SCALAR_VAR( a,b ) a,
#define    ELEM_VECTOR_VAR( a,b ) a,
#define ELEMENT_SCALAR_VAR( a,b ) a,
#define ELEMENT_VECTOR_VAR( a,b ) a,
#define     NEI_VECTOR_VAR( a,b )
#define     NEI_SCALAR_VAR( a,b )
#define     NNI_SCALAR_VAR( a,b )
#define     EEI_SCALAR_VAR( a,b )
#include "contact_variables.def"
#undef    NODE_SCALAR_VAR
#undef    NODE_VECTOR_VAR
#undef    EDGE_SCALAR_VAR
#undef    EDGE_VECTOR_VAR
#undef    FACE_SCALAR_VAR
#undef    FACE_VECTOR_VAR
#undef    ELEM_SCALAR_VAR
#undef    ELEM_VECTOR_VAR
#undef ELEMENT_SCALAR_VAR
#undef ELEMENT_VECTOR_VAR
#undef     NEI_SCALAR_VAR
#undef     NEI_VECTOR_VAR
#undef     NNI_SCALAR_VAR
#undef     EEI_SCALAR_VAR
			    MAX_VARIABLE_ENUM };

  ContactTopology( ContactErrors*, int dimensionality, 
		   int number_analytic_surfaces, int number_node_blocks, 
		   const ContactSearch::ContactNode_Type*, 
		   const int* nodes_per_block,
		   const int* node_exodus_ids,
		   const int* node_ids, 
		   const double* coords,
		   int num_face_blocks, 
		   const ContactSearch::ContactFace_Type*, 
		   const int* faces_per_block, 
		   const int* face_ids, 
		   const int* face_connectivity,
		   const Real* face_lofting_factors,
		   int num_element_blocks, 
		   const ContactSearch::ContactElement_Type*, 
		   const int* elements_per_block, 
		   const int* element_ids, 
		   const int* element_connectivity,
		   int num_comm_partners,const int*, 
		   const int*, const int*,
		   ContactCommBuffer*,
		   MPI_Comm communicator,
                   ContactSearch*,
		   ContactSearch::ContactErrorCode& error_code);

  ContactTopology( ContactErrors*, int dimensionality,
		   int number_analytic_surfaces, 
		   int num_node_blocks, const ContactSearch::ContactNode_Type*,
		   int num_edge_blocks, const ContactSearch::ContactEdge_Type*,
		   int num_face_blocks, const ContactSearch::ContactFace_Type*,
		   int num_elem_blocks, const ContactSearch::ContactElement_Type*,
		   MPI_Comm communicator,
                   ContactSearch* );

  ~ContactTopology();
  
  void UpdateTopology( ContactErrors* Errors, 
  
		       const int* number_node_deaths_per_block,
		       const int* global_node_death_ids,
		       const int* number_face_deaths_per_block,
		       const int* global_face_death_ids,
		       const int* number_elem_deaths_per_block,
		       const int* global_elem_death_ids,
                       
		       const int* number_node_births_per_block,
		       const int* global_node_birth_ids,
		       const int* exodus_node_birth_ids,
		       const int* number_face_births_per_block,
		       const int* global_face_birth_ids,
		       const int* face_connectivity,
		       const int* number_elem_births_per_block,
		       const int* global_elem_birth_ids,
		       const int* elem_connectivity,
        
		       const int  num_node_exports,
                       const int* node_export_id_list,
                       const int* node_export_pids,
		       const int  num_face_exports,
                       const int* face_export_id_list,
                       const int* face_export_pids,
		       const int  num_element_exports,
                       const int* element_export_id_list,
                       const int* element_export_pids,
                                 
                       const int num_nodes,
                       const int* node_host_ids,
                       const int num_faces,
                       const int* face_host_ids,
                       const int num_elements,
                       const int* element_host_ids,
                       
		       const int number_comm_partners, 
		       const int* comm_proc_ids,
		       const int* number_nodes_to_partner, 
		       const int* comm_nodes,
  
                       ContactParOStream& postream,
		       ContactSearch::ContactErrorCode& error_code);
  
  void TopologyDeath( ContactErrors* Errors, 
  
		      const int* number_node_deaths_per_block,
		      const int* global_node_death_ids,
		      const int* number_face_deaths_per_block,
		      const int* global_face_death_ids,
		      const int* number_elem_deaths_per_block,
		      const int* global_elem_death_ids,
  
                      ContactParOStream& postream,
		      ContactSearch::ContactErrorCode& error_code);
  
  void TopologyBirth( ContactErrors* Errors, 
                       
		      const int* number_node_births_per_block,
		      const int* global_node_birth_ids,
		      const int* exodus_node_birth_ids,
                      ContactType* node_birth_entity_types,
		      const int* number_face_births_per_block,
		      const int* global_face_birth_ids,
		      const int* face_connectivity,
		      const int* number_elem_births_per_block,
		      const int* global_elem_birth_ids,
		      const int* elem_connectivity,
  
                      ContactParOStream& postream,
		      ContactSearch::ContactErrorCode& error_code);
  
  void TopologyDLB( ContactErrors* Errors, 
                       
		    const int  num_node_exports,
                    const int* node_export_id_list,
                    const int* node_export_pids,
		    const int  num_face_exports,
                    const int* face_export_id_list,
                    const int* face_export_pids,
		    const int  num_element_exports,
                    const int* element_export_id_list,
                    const int* element_export_pids,
  
                    ContactParOStream& postream,
		    ContactSearch::ContactErrorCode& error_code);
  
  void UpdateInteractions( ContactErrors* Errors, 
                           ContactParOStream& postream,
		           ContactSearch::ContactErrorCode& error_code);
  
  void Add_Analytic_Surface( ContactAnalyticSurface* );
  ContactSearch::ContactErrorCode Set_Analytic_Surface_Configuration( int,
						  	        const Real* );

  ContactSearch::ContactErrorCode Set_NodeBlk_RemainingGap( int, 
                                                            const Real* );
  ContactSearch::ContactErrorCode Set_NodeBlk_GhostingGap( int, 
                                                           const Real* );
  ContactSearch::ContactErrorCode Set_NodeBlk_KinConstr( int, const int*,
                                                         const Real* );
  ContactSearch::ContactErrorCode Set_NodeBlk_Positions( int, const Real* );
  ContactSearch::ContactErrorCode Set_NodeBlk_Positions_2( int, const Real* );
  ContactSearch::ContactErrorCode Set_NodeBlk_Attributes( 
				      ContactSearch::Node_Block_Attribute,
				      int, const Real* );

  ContactSearch::ContactErrorCode Set_FaceBlk_Attributes(
					 ContactSearch::Face_Block_Attribute,
					 int, const Real* );

  inline int Number_of_Primary_Nodes() {return number_of_primary_nodes;};
  inline int Number_of_Primary_Edges() {return number_of_primary_edges;};
  inline int Number_of_Primary_Faces() {return number_of_primary_faces;};
  inline int Number_of_Primary_Elements() {return number_of_primary_elements;};
  inline int Number_of_Nodes() {return number_of_nodes;};
  inline int Number_of_Edges() {return number_of_edges;};
  inline int Number_of_Faces() {return number_of_faces;};
  inline int Number_of_Elements() {return number_of_elements;};
  inline int Number_of_Node_Blocks() {return number_of_node_blocks;};
  inline int Number_of_Face_Blocks() {return number_of_face_blocks;};
  inline int Number_of_Edge_Blocks() {return number_of_edge_blocks;};
  inline int Number_of_Element_Blocks() {return number_of_element_blocks;};
  inline int Number_of_Analytic_Surfaces() 
    { return number_of_analytic_surfaces; };
  inline int Number_of_Added_Analytic_Surfaces() 
    { return number_of_added_analytic_surfaces; };
  inline bool Have_Shells() { return have_shells; };
  inline Real Min_Characteristic_Length() { return computed_characteristic_length?min_characteristic_length:0.0; };
  inline Real Max_Characteristic_Length() { return computed_characteristic_length?max_characteristic_length:0.0; };
  void Compute_Characteristic_Length();

  void Display(ContactParOStream&);
  void Display_Entities(ContactParOStream&, int flag=0);
  void Display_Ghosted_Entities(ContactParOStream&, int flag=0);
  void Display_NodeNode_Interactions( ContactParOStream&, int state = 0 );
  void Display_NodeEntity_Interactions( ContactParOStream&, int state = 0 );
  void Display_FaceFace_Interactions( ContactParOStream&, int state = 0 );
  void Display_FaceCoverage_Interactions( ContactParOStream&, int state = 0 );
  void Display_ElementElement_Interactions( ContactParOStream&, int state = 0 );
  
  void Display_NodeNode_Interactions_Summary( ContactParOStream&, char* margin=NULL, int state = 0 );
  void Display_NodeEntity_Interactions_Summary( ContactParOStream&, unsigned int context, char* margin=NULL, int state = 0 );
  void Display0_NodeEntity_Interactions_Summary( unsigned int context, char* margin=NULL, int state = 0 );
  void Display_NodeSurface_Interactions_Summary( ContactParOStream&, char* margin=NULL, int state = 0 );
  void Display_FaceFace_Interactions_Summary( ContactParOStream&, char* margin=NULL, int state = 0 );
  void Display_FaceCoverage_Interactions_Summary( ContactParOStream&, char* margin=NULL, int state = 0 );
  void Display_ElementElement_Interactions_Summary( ContactParOStream&, char* margin=NULL, int state = 0 );

  ContactSearch::ContactErrorCode 
    Exodus_Output(int Exodus_ID, Real Time,
		  ContactErrors*,     
		  ContactSearchData*,
		  const ContactSearch::Search_Option_Status&,
		  const ContactSearch::Search_Option_Status&,
		  const Real&,
		  const Real&,
		  const ContactSearch::Smoothing_Resolution&,
		  const ContactSearch::Search_Option_Status&,
		  ContactEnforcement* enf );
  ContactSearch::ContactErrorCode 
    ExodusMesh_Output(int Exodus_ID, 
		  ContactErrors*,     
		  ContactSearchData*
		  );

#ifndef CONTACT_NO_MPI
  void GhostTiedFaces();
  void UpdateTiedFaces();
  void DeleteTiedFaces();
  ContactSearch::ContactErrorCode GhostNFIfaces(int);
  void UpdateGhostedNFIfaces();
  void DeleteGhostedNFIfaces();
  void DoGhosting(VariableHandle, const Real &reasonable_gap);
  void DoGhosting_New_NodeFace(VariableHandle, const Real &reasonable_gap);
  void DoCaptureGhosting_New_NodeFace(VariableHandle, const Real &reasonable_gap);
  void DeleteGhosting();
  void UpdateGhostingSetupForNoSecondary();
  void UpdateGhosting();
  ContactZoltanComm*      GetFaceGhostingComm()  { return GhostFaces_ZoltanComm; };
  ContactZoltanComm*      GetOtherGhostingComm() { return GhostOthers_ZoltanComm; };
  ContactZoltanCommUtils* GetGhostingCommSpec()  { return GhostingCommSpec; };
#endif
 
  inline VariableHandle Variable_Handle( VariableHandle_enum Var) 
    {return( Var_Handles[Var]); };

  // Access functions 
  inline int Dimensionality() { return dimensionality; };
  inline ContactNodeBlock** Node_Blocks() { return node_blocks; };
  inline ContactEdgeBlock** Edge_Blocks() { return edge_blocks; };
  inline ContactFaceBlock** Face_Blocks() { return face_blocks; };
  inline ContactElementBlock** Element_Blocks() { return element_blocks; }
  inline ContactNodeBlock* Node_Block(int i) { return node_blocks[i]; };
  inline ContactEdgeBlock* Edge_Block(int i) { return edge_blocks[i]; };
  inline ContactFaceBlock* Face_Block(int i) { return face_blocks[i]; };
  inline ContactElementBlock* Element_Block(int i) { return element_blocks[i]; };
  inline ContactAnalyticSurface* Analytic_Surface( int i ) 
    { return AnalyticSurfaces[i]; };
  inline ContactAnalyticSurface** Analytic_Surfaces()
    { return AnalyticSurfaces; };
  inline ContactNodeBlock** Ghosted_Node_Blocks() { return ghosted_node_blocks; };
  inline ContactEdgeBlock** Ghosted_Edge_Blocks() { return ghosted_edge_blocks; };
  inline ContactFaceBlock** Ghosted_Face_Blocks() { return ghosted_face_blocks; };
  inline ContactElementBlock** Ghosted_Element_Blocks() { return ghosted_element_blocks; };
  inline ContactNodeBlock*  Ghosted_Node_Block(int i) { return ghosted_node_blocks[i]; };
  inline ContactEdgeBlock*  Ghosted_Edge_Block(int i) { return ghosted_edge_blocks[i]; };
  inline ContactFaceBlock*  Ghosted_Face_Block(int i) { return ghosted_face_blocks[i]; };
  inline ContactElementBlock*  Ghosted_Element_Block(int i) { return ghosted_element_blocks[i]; };
  void Size_NodeNode_Interactions( int&, int&);
  void Get_NodeNode_Interactions( int*, int*, int*, int*, int*, Real* );
  void Size_NodeFace_Interactions( int&, int&);
  void Get_NodeFace_Interactions( int*, int*, int*,int*, int*, int*, Real* );
  void Size_NodeSurface_Interactions( int&, int& );
  void Get_NodeSurface_Interactions( int*, int*, int*, Real* );
  void Size_FaceFace_Interactions( int&, int& );
  void Get_FaceFace_Interactions( int*, int*, int*, int*, int*, int*, int*, Real* );
  void Size_FaceCoverage_Interactions( int&, int& );
  void Get_FaceCoverage_Interactions( int*, int*, int*, Real* );
  void Size_ElementElement_Interactions( int&, int& );
  void Get_ElementElement_Interactions( int*, int*, int*, int*, int*, Real* );
  void Update_State();
  bool Faces_Connected( ContactFace<Real>*, ContactFace<Real>* );
  ContactSearch* Search() { return search; }

  inline void Number_of_Nodes( int n ) { number_of_nodes = n; };
  inline void Number_of_Edges( int n ) { number_of_edges = n; };
  inline void Number_of_Faces( int n ) { number_of_faces = n; };
  inline void Number_of_Elements( int n ) { number_of_elements = n; };
#ifdef CONTACT_DEBUG_NODE
  inline int* Debug_Node_Exodus_Ids() { return debug_node_exodus_ids; };
#endif

  // Functions to build topology
  void Connect_Faces_to_Nodes();
  void Connect_Faces_to_Edges();
  void Construct_and_Connect_Edges( ContactSearch::ContactErrorCode& );
  int Get_Faces_Connected_to_Nodes( ContactNode<Real>*, ContactNode<Real>*,
				    ContactFace<Real>**,ContactFace<Real>**,
				    ContactSearch::ContactErrorCode& );
  void Compute_Owners(ContactSearch::ContactErrorCode&);

  void CleanUp( );
  
#ifndef CONTACT_NO_MPI
  inline ContactAsymComm* Node_Asym_Comm() {return Node_AsymComm; };
  inline ContactSymComm* Node_Sym_Comm() {return Node_SymComm; };
  inline ContactSymComm* Edge_Sym_Comm() {return Edge_SymComm; };
  void Connect_Nodes_to_Faces();
  void Connect_Nodes_to_Elements();

  void Compute_Owners_For_Entity( ContactSymComm*, ContactTopologyEntityList*, int );
  void Compute_Edge_Comm_List(ContactSearch::ContactErrorCode& );
  void Find_Shared_Edge_Candidates( std::vector< std::pair<ContactTopologyEntity<Real>*,int> >* );
  void Complete_Edge_Comm_List( int, std::vector< std::pair<ContactTopologyEntity<Real>*,int> >*, int*);
  int  Compare_Edges( ContactEdge<Real>*, ContactEdge<Real>* );
  void Assign_Secondary_Ownership( ContactZoltan*, VariableHandle );
  inline std::vector<ContactInteractionEntity<Real>::entity_data*>* QueryLinkList()
    { return query_linklist;};
  inline void QueryLinkList( std::vector<ContactInteractionEntity<Real>::entity_data*>* list)
    { query_linklist=list;};
#endif

  inline int* UpdateNodeMap() { return old_to_new_node_map; };

  void Delete_All_Interactions( );

  static void SortEntityList( int, ContactTopologyEntity<Real>** );
  static void SortEntityList1( int, ContactEdge<Real>** );

  //
  //  Check if an edge is less than another edge based on the global IDs
  //  of the edge end nodes.
  //
  static bool EdgeLessThan( ContactEdge<Real>* edge1, ContactEdge<Real>* edge2);

#ifdef CONTACT_DEBUG_NODE
  ContactSearch::ContactErrorCode Add_Debug_Node( int );
  void Display_Debug_Node_IDs( ContactParOStream& );
  void Display_Debug_Nodes( ContactParOStream& );
  bool Is_a_Debug_Node(ContactNode<Real>*);
  int Number_Debug_Nodes() { return number_debug_nodes; };
  ContactNode<Real>* Debug_Node( int i ){
    PRECONDITION( i <= number_debug_nodes && i >= 0);
    return debug_nodes[i];
  };
#endif

  void Get_Scratch();
  void Release_Scratch();

  ContactShellHandler* Shell_Handler() { return shell_handler; };
  
  ContactTopologyEntityList* NodeList() {return node_list;};
  ContactTopologyEntityList* EdgeList() {return edge_list;};
  ContactTopologyEntityList* FaceList() {return face_list;};
  ContactTopologyEntityList* ElemList() {return elem_list;};
  ContactTopologyEntityList* PrimaryNodeList() {return primary_node_list;};
  ContactTopologyEntityList* PrimaryFaceList() {return primary_face_list;};
  ContactTopologyEntityList* PrimaryElemList() {return primary_elem_list;};
  
  int Number_NodeNode_Interactions();
  int Number_NodeFace_Interactions();
  int Number_NodeSurface_Interactions();
  int Number_NodeEntity_Interactions();
  int Number_FaceFace_Interactions();
  int Number_FaceCoverage_Interactions();
  int Number_ElementElement_Interactions();
    

  void Compute_Element_Geometry( int POSITION );
  void Compute_Face_Geometry( int POSITION, bool use_proximity, bool loft_shells );
  void Compute_Surface_Geometry( int POSITION, int num_configs, bool use_proximity );
  double Compute_Curvature( const double*, const double*, const double*, 
			    const double*);

  static Real Compute_Angle( const double*, const double*, const double*, 
			 const double*);

  void Compute_Max_Relative_Node_Motion( Real* );

#ifndef CONTACT_NO_MPI
  void MigrateExportedData();
#endif

  void Set_Search_Options( ContactSearch::Search_Option_Status MI_Status,
			   ContactSearch::Search_Option_Status NS_Status,
			   Real SS_curvature,
			   ContactSearch::Search_Option_Status NPC_Status, 
                           ContactSearch::Search_Option_Status AT_Status){
    multiple_interaction_status = MI_Status;
    normal_smoothing_status = NS_Status;
    sharp_smooth_curvature = SS_curvature;
    no_parallel_consistency = NPC_Status;
    auto_tol = AT_Status;
  }
  void RebuildTopologyEntityList( ContactSearch::Search_Option_Status NPC_Status );
  enum TopologyType{ PRIMARY, SECONDARY };
  TopologyType topology_type;
  
  bool HaveGhosting() {return have_global_ghosting|have_local_ghosting;};

  void Set_Up_Variable_Handles();

 protected:

 private:

  enum NodeNode_Interaction_Data_Enum { NN_DIST = 0, 
					SIZE_NODENODE_INTERACTION_DATA };

  enum NodeFace_Interaction_Data_Enum { NF_COORD1 = 0, NF_COORD2, NF_GAP_CUR,
					/*NF_GAP_OLD, */ NF_PUSHBACK_X,
					NF_PUSHBACK_Y, NF_PUSHBACK_Z,
					NF_NORMAL_X, NF_NORMAL_Y,
					NF_NORMAL_Z, NF_SOURCE, NF_NODE_AREA,
                                        NF_GAP_INT,
					SIZE_NODEFACE_INTERACTION_DATA };

  enum NodeSurface_Interaction_Data_Enum { NS_COORD_X = 0, NS_COORD_Y, 
					   NS_COORD_Z, NS_GAP_CUR, 
					   /*NS_GAP_OLD,*/
					   NS_NORMAL_X, NS_NORMAL_Y, 
					   NS_NORMAL_Z, NS_PFNORMAL_X,
                                           NS_PFNORMAL_Y, NS_PFNORMAL_Z,
                                           NS_ENTITY_KEY, NS_NDAREA,
					   SIZE_NODESURFACE_INTERACTION_DATA };
  
  enum ElementElement_Interaction_Data_Enum { EE_VOLUME = 0, 
					      SIZE_ELEMENTELEMENT_INTERACTION_DATA };
					   
#ifdef CONTACT_DEBUG_NODE
  // Debug node
  int number_debug_nodes;
  ContactHostGlobalID** debug_node_global_ids;
  int* debug_node_exodus_ids;
  ContactNode<Real>** debug_nodes;
#endif

  // Error Handler
  ContactErrors* errors;
  char message[81];
  
  ContactSearch* search;  // owner search object

  // handler for shell objects
  ContactShellHandler* shell_handler;

  // MPI Communicator
  MPI_Comm SearchComm;
  ContactCommBuffer* comm_buffer;

  // Size of the Topology
  int dimensionality;
  int* node_block_ids;
  int number_of_primary_nodes;
  int number_of_primary_edges;
  int number_of_primary_faces;
  int number_of_primary_elements;
  int number_of_nodes;
  int number_of_edges;
  int number_of_faces;
  int number_of_elements;
  int number_of_node_blocks;
  int number_of_edge_blocks;
  int number_of_face_blocks;
  int number_of_element_blocks;
  int number_of_analytic_surfaces;
  int number_of_added_analytic_surfaces;
  int number_of_comm_partners;
  bool have_shells;
  bool computed_characteristic_length;
  Real min_characteristic_length;
  Real max_characteristic_length;

  int* old_to_new_node_map;
  
  // Topological Entities
  ContactNodeBlock** node_blocks;
  ContactEdgeBlock** edge_blocks;
  ContactFaceBlock** face_blocks;
  ContactElementBlock** element_blocks;
  ContactAnalyticSurface** AnalyticSurfaces;
  ContactNodeBlock**    ghosted_node_blocks;
  ContactEdgeBlock**    ghosted_edge_blocks;
  ContactFaceBlock**    ghosted_face_blocks;
  ContactElementBlock** ghosted_element_blocks;

#ifndef CONTACT_NO_MPI
  ContactAsymComm* Node_AsymComm;
  ContactSymComm* Node_SymComm;
  ContactSymComm* Edge_SymComm;
  std::vector<ContactInteractionEntity<Real>::entity_data*> *query_linklist;
  ContactZoltanCommUtils* GhostingCommSpec;
  ContactZoltanComm*      GhostFaces_ZoltanComm;
  ContactZoltanComm*      GhostOthers_ZoltanComm;
  ContactZoltanCommUtils* TiedCommSpec;
  ContactZoltanComm*      TiedFaces_ZoltanComm;
  ContactZoltanComm*      TiedNodes_ZoltanComm;
#endif 
  ContactTopologyEntityList* node_list; 
  ContactTopologyEntityList* edge_list;
  ContactTopologyEntityList* face_list;
  ContactTopologyEntityList* elem_list;
  ContactTopologyEntityList* primary_node_list;
  ContactTopologyEntityList* primary_edge_list;
  ContactTopologyEntityList* primary_face_list;
  ContactTopologyEntityList* primary_elem_list;
  
  bool have_global_ghosting;
  bool have_local_ghosting;
  bool have_tied_ghosting;
  
#ifndef CONTACT_NO_MPI
  int         num_ghosted_import;
  LB_ID_TYPE* ghosted_import_gids;
  LB_ID_TYPE* ghosted_import_lids;
  int*        ghosted_import_pids;
  int         num_ghosted_export;
  LB_ID_TYPE* ghosted_export_gids;
  LB_ID_TYPE* ghosted_export_lids;
  int*        ghosted_export_pids;
  int         num_tied_import;
  LB_ID_TYPE* tied_import_gids;
  LB_ID_TYPE* tied_import_lids;
  int*        tied_import_pids;
  int         num_tied_export;
  LB_ID_TYPE* tied_export_gids;
  LB_ID_TYPE* tied_export_lids;
  int*        tied_export_pids;
#endif

  VariableHandle* Var_Handles;

  // Variables
#define    NODE_SCALAR_VAR(a,b) VariableHandle b;
#define    NODE_VECTOR_VAR(a,b) VariableHandle b;
#define    EDGE_SCALAR_VAR(a,b) VariableHandle b;
#define    EDGE_VECTOR_VAR(a,b) VariableHandle b;
#define    FACE_SCALAR_VAR(a,b) VariableHandle b;
#define    FACE_VECTOR_VAR(a,b) VariableHandle b;
#define    ELEM_SCALAR_VAR(a,b) VariableHandle b;
#define    ELEM_VECTOR_VAR(a,b) VariableHandle b;
#define ELEMENT_SCALAR_VAR(a,b) VariableHandle b;
#define ELEMENT_VECTOR_VAR(a,b) VariableHandle b;
#define     NEI_VECTOR_VAR(a,b)
#define     NEI_SCALAR_VAR(a,b)
#define     NNI_SCALAR_VAR( b,a )
#define     EEI_SCALAR_VAR( b,a )
#include "contact_variables.def"
#undef    NODE_SCALAR_VAR 
#undef    NODE_VECTOR_VAR 
#undef    EDGE_SCALAR_VAR 
#undef    EDGE_VECTOR_VAR 
#undef    FACE_SCALAR_VAR 
#undef    FACE_VECTOR_VAR 
#undef    ELEM_SCALAR_VAR
#undef    ELEM_VECTOR_VAR
#undef ELEMENT_SCALAR_VAR
#undef ELEMENT_VECTOR_VAR
#undef     NEI_SCALAR_VAR  
#undef     NEI_VECTOR_VAR  
#undef     NNI_SCALAR_VAR
#undef     EEI_SCALAR_VAR

  // Search Options (needed for surface geometry)
  ContactSearch::Search_Option_Status multiple_interaction_status;
  ContactSearch::Search_Option_Status no_parallel_consistency;
  ContactSearch::Search_Option_Status normal_smoothing_status;
  ContactSearch::Search_Option_Status auto_tol;
  Real sharp_smooth_curvature;

  // Two configurations specified
  bool two_configurations;

  bool scratch_set;

#ifndef CONTACT_NO_EXODUS_OUTPUT
  ContactSearch::ContactErrorCode 
    Exodus_Output_Results( int, Real, 
			   ContactErrors*,
			   ContactEnforcement* enforcement,
			   int num_proc_nodes,
			   int number_of_elem_blocks, 
			   int Need_Sphere_Element_Block,
			   int number_of_elems,
			   int number_edges_owned,
			   int num_nni, int max_nni,
                           int num_nfi, int max_nfi,
			   int num_nsi, int max_nsi,
			   int num_ffi, int max_ffi,
			   int max_ffi_verts,
			   int num_fci, int max_fci,
			   int max_fci_verts,
			   int num_eei, int max_eei,
			   const ContactSearch::Search_Option_Status&
			   normal_smoothing_status,
			   const ContactSearch::Search_Option_Status&
			   multiple_interaction_status,
			   const Real& sharp_smooth_curvature, 
			   const Real& normal_smoothing_distance, 
			   const ContactSearch::Smoothing_Resolution &
			   smoothing_resolution,
			   const ContactSearch::Search_Option_Status &
			   compute_node_areas );
    
#endif
};

inline bool ContactTopology::EdgeLessThan( ContactEdge<Real>* edge1, ContactEdge<Real>* edge2) {
  PRECONDITION(edge1->Nodes_Per_Edge() == edge2->Nodes_Per_Edge());    
  //
  //  Extract the edge nodes
  //
  ContactHostGlobalID &edge1_node0 = edge1->Node(0)->Global_ID();
  ContactHostGlobalID &edge1_node1 = edge1->Node(1)->Global_ID();
  ContactHostGlobalID &edge2_node0 = edge2->Node(0)->Global_ID();
  ContactHostGlobalID &edge2_node1 = edge2->Node(1)->Global_ID();
  //
  //  Sort the edge end nodes to such that node0 is less than node 1 for both edges
  //
  if(edge1_node0 < edge1_node1) {
    if(edge2_node0 < edge2_node1) {
      if(edge1_node0 < edge2_node0 || (edge1_node0 == edge2_node0 && edge1_node1 < edge2_node1)) {
	return true;
      }
    } else {
      if(edge1_node0 < edge2_node1 || (edge1_node0 == edge2_node1 && edge1_node1 < edge2_node0)) {
	return true;
      }
    }
  } else {
    if(edge2_node0 < edge2_node1) {
      if(edge1_node1 < edge2_node0 || (edge1_node1 == edge2_node0 && edge1_node0 < edge2_node1)) {
	return true;
      }
    } else {
      if(edge1_node1 < edge2_node1 || (edge1_node1 == edge2_node1 && edge1_node0 < edge2_node0)) {
	return true;
      }
    }
  }
  return false;
}


#endif // #ifndef ContactTopology_h_
