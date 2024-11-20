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


#include "Contact_Communication.h"
#include "ContactAnalyticCylinderInside.h"
#include "ContactAnalyticCylinderOutside.h"
#include "ContactAnalyticPlane.h"
#include "ContactAnalyticSphere.h"
#include "ContactSearch.h"
#include "ContactEnforcement.h"
#include "ContactSearchData.h"
#include "ContactTopology.h"
#include "ContactNode.h"
#include "ContactNodeBlock.h"
#include "ContactEdgeBlock.h"
#include "ContactFaceBlock.h"
#include "ContactElementBlock.h"
#include "ContactQuadFaceL4.h"
#include "ContactQuadFaceQ8.h"
#include "ContactQuadFaceQ9.h"
#include "ContactTriFaceL3.h"
#include "ContactTriFaceQ6.h"
#include "ContactShellNode.h"
#include "ContactShellQuadFaceL4.h"
#include "ContactShellTriFaceL3.h"
#include "ContactLineFaceL2.h"
#include "ContactLineFaceQ3.h"
#include "ContactLineEdgeL2.h"
#include "ContactLineEdgeQ3.h"
#include "ContactHexElementL8.h"
#include "ContactNodeNodeInteraction.h"
#include "ContactNodeFaceInteraction.h"
#include "ContactNodeSurfaceInteraction.h"
#include "ContactFaceFaceInteraction.h"
#include "ContactFaceCoverageInteraction.h"
#include "ContactElementElementInteraction.h"
#include "ContactTopologyEntityHash.h"
#include "ContactErrors.h"
#include "ContactParOStream.h"
#include "ContactFixedSizeAllocator.h"
#include "ContactBoundingBox.h"
#include "ContactBoundingBoxHierarchy.h"
#include "ContactBoundingBoxHierarchy_Int.h"
#include "ContactSymComm.h"
#include "ContactCommBuffer.h"

#ifndef CONTACT_NO_MPI
#include "zoltan.h"
#include "ContactZoltan.h"
#include "ContactZoltanComm.h"
#endif

#ifndef CONTACT_NO_EXODUS_OUTPUT
#include "exodusII.h"
#endif

#include <cstring>
#include <cmath>

void
ContactSearch::Create_Search_Topology( VariableHandle POSITION )
{
  if( contact_number_of_processors( SearchComm ) == 1 ){
    search_topology = primary_topology;
    return;
  }
#ifndef CONTACT_NO_MPI
#if !defined (CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
  timer.Start_Timer( create_secondary_time );
#endif
  
  if (no_secondary==INACTIVE) {
  
#if CONTACT_DEBUG_PRINT_LEVEL>=0
    int copy_nodes = 0;
    int copy_faces = 0;
    int copy_elems = 0;
#endif

    int my_proc_id = contact_processor_number(SearchComm);
    ContactZoltanLID zoltanLID;
    ContactZoltanGID zoltanGID;
  
#if CONTACT_DEBUG_PRINT_LEVEL>=2
    postream << "    Secondary Decomposition\n";
#endif

    //=================================================================
    // Use ZOLTAN/RCB to partition the nodes in the secondary topology
    //=================================================================
#ifdef CONTACT_TIMINGS
    timer.Start_Timer( rcb_time );
#endif
    zoltan->Set_GeomCallBacks(this);
    zoltan->Position(POSITION);
    zoltan->Balance();

#if CONTACT_DEBUG_PRINT_LEVEL>=10
    int       num_import   = zoltan->Num_Import();
    LB_ID_PTR import_lids  = zoltan->Import_Local_IDs();
    LB_ID_PTR import_gids  = zoltan->Import_Global_IDs();
    int*      import_procs = zoltan->Import_Procs();
    int       num_export   = zoltan->Num_Export();
    LB_ID_PTR export_lids  = zoltan->Export_Local_IDs();
    LB_ID_PTR export_gids  = zoltan->Export_Global_IDs();
    int*      export_procs = zoltan->Export_Procs();
    if (export_gids!=NULL && export_lids!=NULL && export_procs!=NULL) {
      postream << "RCB EXPORTS:\n";
      for (int i=0; i< num_export; ++i){
	postream << "  " << i << ":  GID = "
		 << "T" << zoltanGID.Type(&export_gids[i*ZOLTAN_GID_SIZE]) << "  "
		 << "H" << zoltanGID.Hi(&export_gids[i*ZOLTAN_GID_SIZE]) << "  "
		 << "L" << zoltanGID.Lo(&export_gids[i*ZOLTAN_GID_SIZE]) << "  "
		 << "    LID = "
		 << "T" << zoltanLID.Type(&export_lids[i*ZOLTAN_LID_SIZE]) << "  "
		 << "I" << zoltanLID.Index(&export_lids[i*ZOLTAN_LID_SIZE]) << "  "
		 << "    PROC = "
		 << export_procs[i]<< "\n";
      }
    }
    if (import_gids!=NULL && import_lids!=NULL && import_procs!=NULL) {
      postream << "RCB IMPORTS:\n";
      for (int i=0; i< num_import; ++i){
	postream << "  " << i << ":  GID = "
		 << "T" << zoltanGID.Type(&import_gids[i*ZOLTAN_GID_SIZE]) << "  "
		 << "H" << zoltanGID.Hi(&import_gids[i*ZOLTAN_GID_SIZE]) << "  "
		 << "L" << zoltanGID.Lo(&import_gids[i*ZOLTAN_GID_SIZE]) << "  "
		 << "    LID = "
		 << "T" << zoltanLID.Type(&import_lids[i*ZOLTAN_LID_SIZE]) << "  "
		 << "H" << zoltanLID.Index(&import_lids[i*ZOLTAN_LID_SIZE]) << "  "
		 << "    PROC = "
		 << import_procs[i]<< "\n";
      }
    }
#endif
#ifdef CONTACT_TIMINGS
    timer.Stop_Timer( rcb_time );
#endif

    /*=========================================================================
      =                                                                       =
      =                T O P _ L E V E L   O B J E C T S                      =
      =                                                                       =
      =========================================================================*/
#ifdef CONTACT_TIMINGS
    timer.Start_Timer( secondary_owner_migration_time );
#endif
    LB_ID_TYPE zoltan_lid[ZOLTAN_LID_SIZE];
    LB_ID_TYPE zoltan_gid[ZOLTAN_GID_SIZE];
    int        zoltan_pid;
  
    int num_primary_faces = primary_topology->Number_of_Faces();
    int num_primary_elems = primary_topology->Number_of_Elements();
    ContactFace<Real>** primary_faces = 
      reinterpret_cast<ContactFace<Real>**>(primary_topology->FaceList()->EntityList());
    ContactElement** primary_elems = 
      reinterpret_cast<ContactElement**>(primary_topology->ElemList()->EntityList());
    
    for (int i=0; i<primary_topology->Number_of_Node_Blocks(); ++i) {
      secondary_topology->Node_Block(i)->NodeList()->UseHash(0);
    }
    for (int i=0; i<primary_topology->Number_of_Edge_Blocks(); ++i) {
      secondary_topology->Edge_Block(i)->EdgeList()->UseHash(0);
    }
    for (int i=0; i<primary_topology->Number_of_Face_Blocks(); ++i) {
      secondary_topology->Face_Block(i)->FaceList()->UseHash(0);
    }
    for (int i=0; i<primary_topology->Number_of_Element_Blocks(); ++i) {
      secondary_topology->Element_Block(i)->ElemList()->UseHash(0);
    }
  
    /*=========================================================================
      =                                                                       =
      =                  A N A L Y T I C   S U R F A C E S                    =
      =                                                                       =
      =========================================================================*/
  
    // Since every processor has all of the analytic surfaces, just create 
    // "copies" of my primary surfaces and give to the secondary topology.
    for(int i=0 ; i<primary_topology->Number_of_Analytic_Surfaces() ; ++i ){
      switch( primary_topology->Analytic_Surface(i)->Surface_Type() ){
      case ContactSearch::PLANE: {
	ContactAnalyticPlane* Plane = reinterpret_cast<ContactAnalyticPlane*>
	  ( primary_topology->Analytic_Surface(i) );
	secondary_topology->Add_Analytic_Surface( 
						 new ContactAnalyticPlane( *Plane ) );
	break;
      }	       
      case ContactSearch::SPHERE:{
	ContactAnalyticSphere* Sphere = reinterpret_cast<ContactAnalyticSphere *>
	  ( primary_topology->Analytic_Surface(i) );
	secondary_topology->Add_Analytic_Surface( 
						 new ContactAnalyticSphere( *Sphere ) );
	break;
      }
      case ContactSearch::CYLINDER_INSIDE:{
	ContactAnalyticCylinderInside* cyl_ins = 
	  reinterpret_cast<ContactAnalyticCylinderInside *>
	  ( primary_topology->Analytic_Surface(i) );
	secondary_topology->Add_Analytic_Surface(
						 new ContactAnalyticCylinderInside( *cyl_ins ) );
	break;
      }
      case ContactSearch::CYLINDER_OUTSIDE:{
	ContactAnalyticCylinderOutside* cyl_out = 
	  reinterpret_cast<ContactAnalyticCylinderOutside *> 
	  ( primary_topology->Analytic_Surface(i) );
	secondary_topology->Add_Analytic_Surface(
						 new ContactAnalyticCylinderOutside( *cyl_out ) );
	break;
      }
      default:
        POSTCONDITION(0);
        break;
      }
    }

    //========================================================================
    // Create asymmetric communication object to store the communication plan
    // for faces.  Note that this is currently allocated to an arbitrary size.
    //  ASG: size must contain both enough to hold all the nodes on the 
    //       faces in the primary and the number of potential interactions
    //       to make sure it is big enough.
    //========================================================================
    ContactZoltanComm Toplevel_ZoltanComm( ContactZoltanComm::ZOLTAN_EXPORT);

    //=================================================
    // In the primary decomposition, assign the 
    // ownership for the secondary decomposition
    //=================================================
    primary_topology->Assign_Secondary_Ownership(zoltan, POSITION);
  
#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
    bytes_nodes = 0.0;
    bytes_faces = 0.0;
    bytes_elems = 0.0;
#endif
  
    //=========================================================================
    //  T O P _ L E V E L   F A C E S
    //=========================================================================
  
    //===================================================
    // Compute tolerance for the drop box to do ghosting 
    //===================================================
    VariableHandle CURRENT_POSITION = primary_topology->Variable_Handle( ContactTopology::Current_Position );
#ifdef CONTACT_DEBUG
    VariableHandle NODE_GHOST_GAP   = primary_topology->Variable_Handle( ContactTopology::Node_Ghost_Gap );
    VariableHandle REMAINING_GAP    = primary_topology->Variable_Handle( ContactTopology::Remaining_Gap );
#endif
    Real user_search_tol = search_data->Max_Search_Tolerance();
    int num_config = 1;
    if (POSITION!=CURRENT_POSITION) num_config = 2;
    //
    //  Determine the actual bounding box of each processors set of elements.  
    //  This may be significantly smaller than box_tol + RCB size.  Bounding 
    //  boxes will be assemebed into a global list which will be reduced and 
    //  copied to all processors.  This will allow culling out ghosts that are 
    //  not strictly nessecary.
    //  
    //  NKC note:  This calculation is not huge, but is non trivial.  Definetly 
    //  worth it for SPH problems or problems where normal smoothing is used and 
    //  the object sizes vary significantly.  There might also may be some work 
    //  reduction for other problem types.  Need to determine if it is significant
    //  enough to warrant peforming the calculation.  Also may need to do similar 
    //  calcs for other portions of tolerance.  For example, huge velocity on one 
    //  node may yeild universally huge ghosting rads, big tolerance on one
    //  interaction could be a problem.  Large adhesion or other enforcement 
    //  tolerance on a small subset of nodes yields very large ghosting radius.
    //
    int total_num_procs = contact_number_of_processors( SearchComm );  
    ObjectBoundingBox *proc_box_array = NULL;
    proc_box_array = new ObjectBoundingBox[total_num_procs];
    Create_Processor_Bounding_Boxes(proc_box_array, POSITION, total_num_procs);
    int hierarchy_size = 2 * total_num_procs - 1;
    ObjectBoundingBoxHierarchy *proc_box_hierarchy = NULL;
    if(hierarchy_size > 0) {
      proc_box_hierarchy = new ObjectBoundingBoxHierarchy[hierarchy_size];
      ObjectBoundingBoxHierarchy::create_hierarchy(proc_box_hierarchy, proc_box_array, total_num_procs);
    }
    int* proc_list = new int [total_num_procs];

    //=======================================================================
    // Ghost all faces to determine where they should be sent.  This section 
    // builds the export part of the asymmetric face communication object.
    //=======================================================================
#if CONTACT_DEBUG_PRINT_LEVEL>=2
    postream<<"Step "<<step_number<<"\n";
    postream<<"  auto_tol        = "<<auto_tol<<"\n";
    postream<<"  box_inflation   = "<<box_inflation<<"\n";
    postream<<"  user_search_tol = "<<user_search_tol<<"\n";
#endif   
    for (int ii=0; ii<num_primary_faces; ++ii) {
      ContactFace<Real>* face = primary_faces[ii];
      //if (!(face->IsMaster())) continue;
      // get bounding box
      ContactBoundingBox face_current_box;
      ContactBoundingBox face_predicted_box;
      ContactBoundingBox object_box;
      face->ComputeBoundingBoxForSearch(num_config, 
				        CURRENT_POSITION,
				        POSITION, 
                                        auto_tol,
                                        box_inflation,
                                        user_search_tol,
	                                face_current_box,
	                                face_predicted_box,
                                        object_box);
      //
      //  NKC note, really want to search each processor bounding box
      //  vs. the current object box.  The processor bounding boxes can overlap.
      //  However, the RCB boxes defined by zoltan will not overlap.
      //  By expanding each object box by the global object size tolerance it ensures
      //  the object box will overlap any potential zoltan RCB bounding box.
      //  May need to not use zoltan and do these calcs by hand later to 
      //  optimize the search on actual processor box size.  Using the global object tolerance
      //  here has the potential to find many processor overlaps that need to by thrown out later
      //  by the proc_box_array overlap calculation
      // 
      //
      int numprocs = 0;
      if(proc_box_hierarchy != NULL) ObjectBoundingBoxHierarchy::search_for_overlap_loop(proc_box_hierarchy, object_box, proc_list, numprocs);
      // loop over the proc_list from the box drop -- if the local 
      // processor is in the list, mark the current face as needing to
      // be copied locally from primary to secondary decomposition (set
      // temp_tag to 1)
      for (int i=0; i<numprocs; ++i) {
	//
	//  Double check that the current bounding box actually overlaps the 
	//  true processor bounding box.  If not, do not ghost this object
	//
	int proc_num = proc_list[i];
	if (proc_num == my_proc_id) {
	  //
	  // on processor - copy from primary to secondary topology
	  //
	  ContactFace<Real>* secondary_face = New_ContactFace(face->FaceType(), allocators);
          POSTCONDITION(secondary_face!=NULL);
#ifdef CONTACT_OLD_XFER
          secondary_face->Copy(face, enable_off_face_tracking);
#else        
          secondary_face->Copy_ForSecondary(face, 
                                            enable_off_face_tracking,
                                            physical_face_algorithm==PF_EDGE_BASED);
#endif
	  ContactFaceBlock* face_block = 
	    secondary_topology->Face_Block(face->BlockID());
	  face_block->Insert_Face(secondary_face);
#if CONTACT_DEBUG_PRINT_LEVEL>=0
	  ++copy_faces;
#endif
        
	} else {
      
	  // off processor - add to communication object
	  zoltan_pid = proc_list[i];
	  face->ZoltanLID(CT_FACE, zoltan_lid);
	  face->ZoltanGID(CT_FACE, zoltan_gid);
	  Toplevel_ZoltanComm.Add_Export(zoltan_lid,zoltan_gid,zoltan_pid,0);
	}
      }
    }
  
    //=========================================================================
    //  T O P _ L E V E L   E L E M E N T S
    //=========================================================================
  

    //=======================================================================
    // Ghost all elements to determine where they should be sent.  This section 
    // builds the export part of the asymmetric element communication object.
    //=======================================================================
    for (int ii=0; ii<num_primary_elems; ++ii) {
      ContactElement* element = primary_elems[ii];
      //if (!(element->IsMaster())) continue;
      // get bounding box
      ContactBoundingBox object_box;
      element->ComputeBoundingBoxForSearch(num_config, 
				           CURRENT_POSITION,
				           POSITION, 
                                           auto_tol,
                                           box_inflation,
                                           user_search_tol,
                                           object_box);
      int numprocs = 0;
      if(proc_box_hierarchy != NULL) ObjectBoundingBoxHierarchy::search_for_overlap_loop(proc_box_hierarchy,object_box, proc_list, numprocs);
      // loop over the proc_list from the box drop -- if the local 
      // processor is in the list, mark the current element as needing to
      // be copied locally from primary to secondary decomposition (set
      // temp_tag to 1)
      for (int i=0; i<numprocs; ++i) {
	//
	//  Double check that the current bounding box actually overlaps the true processor bounding box.  If not, do not
	//  ghost this object
	//
	int proc_num = proc_list[i];
	if (proc_num == my_proc_id) {
	  // on processor - copy from primary to secondary topology
	  ContactElement* secondary_element=NULL;
	  switch(element->ElementType()){
	  case ContactSearch::HEXELEMENTL8 :
	    secondary_element = ContactHexElementL8::new_ContactHexElementL8(
									     allocators );
	    break;
	  case ContactSearch::CARTESIANHEXELEMENTL8 :
	    secondary_element = ContactCartesianHexElementL8::new_ContactCartesianHexElementL8(
											       allocators );
	    break;
	  default:
	    errors->Add_Error_Message("Unknown Element Type in Secondary Decomposition");
	    POSTCONDITION( false );
	  }
          POSTCONDITION(secondary_element!=NULL);	
#ifdef CONTACT_OLD_XFER
          secondary_element->Copy(element);
#else
          secondary_element->Copy_ForSecondary(element);
#endif
	  ContactElementBlock* element_block = 
	    secondary_topology->Element_Block(element->BlockID());
	  element_block->Insert_Element(secondary_element);
#if CONTACT_DEBUG_PRINT_LEVEL>=0
	  ++copy_elems;
#endif
        
	} else {
      
	  // off processor - add to communication object
	  zoltan_pid = proc_list[i];
	  element->ZoltanLID(CT_ELEMENT, zoltan_lid);
	  element->ZoltanGID(CT_ELEMENT, zoltan_gid);
	  Toplevel_ZoltanComm.Add_Export(zoltan_lid,zoltan_gid,zoltan_pid);
	}
      }
    }
  
#if CONTACT_DEBUG_PRINT_LEVEL>=1006
    {
      int nnodes_print = 0; 
      for (int i=0; i<secondary_topology->Number_of_Node_Blocks(); ++i) {
        nnodes_print += secondary_topology->Node_Block(i)->Number_of_Nodes();
      }
      int nfaces_print = 0; 
      for (int i=0; i<secondary_topology->Number_of_Face_Blocks(); ++i) {
        nfaces_print += secondary_topology->Face_Block(i)->Number_of_Faces();
      }
      postream << "Search Topology: After Toplevel Copy\n";
      postream << "  Number of Nodes    = " << nnodes_print << "\n";
      postream << "  Number of Faces    = " << nfaces_print << "\n";
    }
#endif

    //=========================================================================
    //  T O P _ L E V E L   N O D E S   ( B L O C K S   1 - N )
    //=========================================================================
    // In this block, we are only moving "POINTS" not "NODES", 
    // we know they are all ContactNode<Real>'s not ContactShellNodes 
    // so I won't even check
    VariableHandle NODE_RADIUS = primary_topology->Variable_Handle( ContactTopology::Node_Radius );
    int base = 0;
    if (primary_topology->Number_of_Face_Blocks()+
	primary_topology->Number_of_Element_Blocks() > 0) base=1;
    ContactTopologyEntityList* primary_node_list = primary_topology->NodeList();
    for(int i=base; i<primary_topology->Number_of_Node_Blocks(); ++i) {
      ContactNodeBlock* node_block = primary_topology->Node_Block(i);
      if( node_block->Type() == POINT ){
	int nnodes = primary_node_list->BlockNumEntities(i);
	ContactNode<Real>** nodes  = reinterpret_cast<ContactNode<Real>**>(primary_node_list->BlockEntityList(i));
	if (node_block->Has_Radius_Attributes()) {
	  for (int j=0; j<nnodes; ++j) {
	    ContactNode<Real>* node = nodes[j];
	    //
	    //  Skip any nodes not owned by this processor
	    //
            if(node->Owner() != my_proc_id) continue;

	    if (!node->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)) continue;
	    Real* position =  node->Variable(POSITION);
	    Real  radius   = *node->Variable(NODE_RADIUS);

	    ContactBoundingBox object_box;
	    object_box.add_sphere(position, radius);
	    int numprocs = 0;
	    if(proc_box_hierarchy != NULL) ObjectBoundingBoxHierarchy::search_for_overlap_loop(proc_box_hierarchy,object_box, proc_list, numprocs);
	    // loop over the proc_list from the box drop -- if the local 
	    // processor is in the list, mark the current node as needing to
	    // be copied locally from primary to secondary decomposition (set
	    // temp_tag to 1)
	    for (int k=0; k<numprocs; ++k) {
	      //
	      //  Double check that the current bounding box actually overlaps the true processor bounding box.  If not, do not
	      //  ghost this object
	      //
	      if (proc_list[k] == my_proc_id) {
		// on processor - copy from primary to secondary topology
		ContactNode<Real>* new_node = 
		  ContactNode<Real>::new_ContactNode( allocators, NODE );
#ifdef CONTACT_OLD_XFER
                new_node->Copy( node, 1 );
#else
                new_node->Copy_ForSecondary( node, 1 );
#endif
		secondary_topology->Node_Block(i)->Insert_Node(new_node);
#if CONTACT_DEBUG_PRINT_LEVEL>=0
		++copy_nodes;
#endif
	      } else {
		// off processor - add to communication object
		zoltan_pid = proc_list[k];
		node->ZoltanLID(CT_NODE, zoltan_lid);
		node->ZoltanGID(CT_NODE, zoltan_gid);
		Toplevel_ZoltanComm.Add_Export(zoltan_lid,zoltan_gid,zoltan_pid);
	      }   
	    }
	  }
	} else {
	  for (int j=0; j<nnodes; ++j) {
	    ContactNode<Real>* node = nodes[j];

	    //
	    //  Skip any nodes not owned by this processor
	    //
            if(node->Owner() != my_proc_id) continue;
	    if (!node->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)) continue;
	    if (node->Secondary_Owner()==my_proc_id) {
	      // on processor - copy from primary to secondary topology
	      ContactNode<Real>* new_node = 
		ContactNode<Real>::new_ContactNode( allocators, NODE );
#ifdef CONTACT_OLD_XFER
              new_node->Copy( node, 1 );
#else
              new_node->Copy_ForSecondary( node, 1 );
#endif
	      secondary_topology->Node_Block(i)->Insert_Node(new_node);
#if CONTACT_DEBUG_PRINT_LEVEL>=0
	      ++copy_nodes;
#endif
	    } else {
	      // off processor - add to communication object
	      node->ZoltanLID(CT_NODE, zoltan_lid);
	      node->ZoltanGID(CT_NODE, zoltan_gid);
	      zoltan_pid = node->Secondary_Owner();
	      Toplevel_ZoltanComm.Add_Export(zoltan_lid,zoltan_gid,zoltan_pid);
	    }
	  }
	}
      }
    }
    if(proc_box_hierarchy) delete [] proc_box_hierarchy;
    if(proc_box_array) delete [] proc_box_array;
    delete [] proc_list;

    int       num_toplevel_import   = -1;
    LB_ID_PTR import_toplevel_gids  = NULL;
    LB_ID_PTR import_toplevel_lids  = NULL;
    int*      import_toplevel_procs = NULL;
    int       num_toplevel_export   = Toplevel_ZoltanComm.Num_Export();
    LB_ID_PTR export_toplevel_gids  = Toplevel_ZoltanComm.Export_GIDS();
    LB_ID_PTR export_toplevel_lids  = Toplevel_ZoltanComm.Export_LIDS();
    int*      export_toplevel_procs = Toplevel_ZoltanComm.Export_Procs();

    //================================================
    // Migrate all the off processor toplevel objects
    //================================================
#ifdef CONTACT_TIMINGS
    timer.Start_Timer( owner_help_migrate_time );
#endif
    zoltan->Set_ToplevelExportCallBacks(this);
    zoltan->Help_Migrate(num_toplevel_import,  import_toplevel_gids, 
			 import_toplevel_lids, import_toplevel_procs, 
			 num_toplevel_export,  export_toplevel_gids,
			 export_toplevel_lids, export_toplevel_procs);
#ifdef CONTACT_TIMINGS
    timer.Stop_Timer( owner_help_migrate_time );
#endif

    int num_face_import = 0;
    int num_elem_import = 0;
    for (int i=0; i<secondary_topology->Number_of_Face_Blocks(); ++i) {
      num_face_import += secondary_topology->Face_Block(i)->Number_of_Faces();
    }
    for (int i=0; i<secondary_topology->Number_of_Element_Blocks(); ++i) {
      num_elem_import += secondary_topology->Element_Block(i)->Number_of_Elements();
    }
    ContactZoltanComm GhostOthers_ZoltanComm( ContactZoltanComm::ZOLTAN_IMPORT);
  
    int nnodes = primary_topology->NodeList()->NumEntities();
    ContactNode<Real>** Prim_Nodes = reinterpret_cast<ContactNode<Real>**>(primary_topology->NodeList()->EntityList());
    for (int i=0; i<nnodes; ++i) {
      Prim_Nodes[i]->temp_tag = 0;
    }
       
    for(int i=0 ; i<secondary_topology->Number_of_Face_Blocks() ; ++i ){
      ContactTopologyEntity<Real>* entity;
      ContactBlockEntityList* block_face_list = secondary_topology->Face_Block(i)->FaceList();
      block_face_list->IteratorStart();
      while( (entity=block_face_list->IteratorForward()) ){
	ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
	ContactTopologyEntity<Real>::connection_data *node_info = face->NodeInfo();
	POSTCONDITION( node_info );
	for(int j=0 ; j<face->Nodes_Per_Face() ; ++j){
	  if (node_info[j].owner==my_proc_id) {
	    ContactNode<Real>* node = static_cast<ContactNode<Real> *>(primary_topology->NodeList()->Find(&node_info[j]));
	    POSTCONDITION(node);
	    node->temp_tag = 1;
	  } else {
	    zoltan_pid = node_info[j].owner;
	    ContactHostGlobalID GID( node_info[j].host_gid[0], 
				     node_info[j].host_gid[1] );
	    zoltanLID.ZoltanLID(CT_NODE, node_info[j].owner_proc_array_index, zoltan_lid);
	    zoltanGID.ZoltanGID(CT_NODE, &GID, zoltan_gid);
	    GhostOthers_ZoltanComm.Add_Import( zoltan_lid, zoltan_gid, zoltan_pid );
	  }
	}
      }
    }
    for(int i=0 ; i<secondary_topology->Number_of_Element_Blocks() ; ++i ){
      ContactTopologyEntity<Real>* entity;
      ContactBlockEntityList* block_element_list = 
	secondary_topology->Element_Block(i)->ElemList();
      block_element_list->IteratorStart();
      while( (entity=block_element_list->IteratorForward()) ){
	ContactElement* element = static_cast<ContactElement*>(entity);
	ContactTopologyEntity<Real>::connection_data *node_info = element->NodeInfo();
	POSTCONDITION( node_info );
	for(int j=0 ; j<element->Nodes_Per_Element() ; ++j){
	  if (node_info[j].owner==my_proc_id) {
	    ContactNode<Real>* node = static_cast<ContactNode<Real> *>(primary_topology->NodeList()->Find(&node_info[j]));
	    POSTCONDITION(node);
	    node->temp_tag = 1;
	  } else {
	    zoltan_pid = node_info[j].owner;
	    ContactHostGlobalID GID( node_info[j].host_gid[0], 
				     node_info[j].host_gid[1] );
	    zoltanLID.ZoltanLID(CT_NODE, node_info[j].owner_proc_array_index, zoltan_lid);
	    zoltanGID.ZoltanGID(CT_NODE, &GID, zoltan_gid);
	    GhostOthers_ZoltanComm.Add_Import( zoltan_lid, zoltan_gid, zoltan_pid );
	  }
	}
      }
    }
    for (int i=0; i<nnodes; ++i) {
      ContactNode<Real>* node = Prim_Nodes[i];
      if (node->temp_tag) {
	ContactNode<Real>* new_node=NULL;
	switch(node->Base_Type()){
	case CT_NODE :
	  new_node = ContactNode<Real>::new_ContactNode( allocators, NODE );
	  break;
	case CT_SHELL_NODE :
	  new_node = ContactShellNode::new_ContactShellNode( allocators, NODE );
	  break;
	default:
	  errors->Add_Error_Message("Unknown Node Type in Secondary Decomposition");
	  POSTCONDITION( false );
	}
        POSTCONDITION(new_node!=NULL);
	int state = (node->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE) && node->Secondary_Owner()==my_proc_id)?STATE_1:STATE_NONE;
#ifdef CONTACT_OLD_XFER
        new_node->Copy( node, state ); 
#else
        new_node->Copy_ForSecondary( node, state );
#endif
	secondary_topology->Node_Block(node->BlockID())->Insert_Node(new_node);
#if CONTACT_DEBUG_PRINT_LEVEL>=0
	++copy_nodes;
#endif
      }
    }

    num_toplevel_import   = GhostOthers_ZoltanComm.Num_Import();
    import_toplevel_gids  = GhostOthers_ZoltanComm.Import_GIDS();
    import_toplevel_lids  = GhostOthers_ZoltanComm.Import_LIDS();
    import_toplevel_procs = GhostOthers_ZoltanComm.Import_Procs();
    num_toplevel_export   = -1;
    export_toplevel_gids  = NULL;
    export_toplevel_lids  = NULL;
    export_toplevel_procs = NULL;

    //================================================
    // Migrate all the off processor derived objects
    //================================================
#ifdef CONTACT_TIMINGS
    timer.Start_Timer( owner_help_migrate_time );
#endif
    zoltan->Set_ToplevelExportCallBacks(this);
    zoltan->Help_Migrate(num_toplevel_import,  import_toplevel_gids, 
			 import_toplevel_lids, import_toplevel_procs, 
			 num_toplevel_export,  export_toplevel_gids,
			 export_toplevel_lids, export_toplevel_procs);
#ifdef CONTACT_TIMINGS
    timer.Stop_Timer( owner_help_migrate_time );
#endif

#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
    bytes_nodes /= 1024;
    bytes_faces /= 1024;
    bytes_elems /= 1024;
    postream<<"    Data Sent by Zoltan when creating secondary topology...\n";
    postream<<"      nodes:  "<<bytes_nodes<<" Kb\n";
    postream<<"      faces:  "<<bytes_faces<<" Kb\n";
    postream<<"      elems:  "<<bytes_elems<<" Kb\n";
#endif

#ifdef CONTACT_TIMINGS
    timer.Stop_Timer( secondary_owner_migration_time );
#endif
  
#if CONTACT_DEBUG_PRINT_LEVEL>=1006   
    {
      int nnodes_print = 0; 
      for (int i=0; i<secondary_topology->Number_of_Node_Blocks(); ++i) {
        nnodes_print += secondary_topology->Node_Block(i)->Number_of_Nodes();
      }
      int nfaces_print = 0; 
      for (int i=0; i<secondary_topology->Number_of_Face_Blocks(); ++i) {
        nfaces_print += secondary_topology->Face_Block(i)->Number_of_Faces();
      }
      postream << "Search Topology: After Toplevel Migrate\n";
      postream << "  Number of Nodes    = " << nnodes_print << "\n";
      postream << "  Number of Faces    = " << nfaces_print << "\n";
    }
#endif

#if CONTACT_DEBUG_PRINT_LEVEL>=9
    num_toplevel_export   = Toplevel_ZoltanComm.Num_Export();
    export_toplevel_gids  = Toplevel_ZoltanComm.Export_GIDS();
    export_toplevel_lids  = Toplevel_ZoltanComm.Export_LIDS();
    export_toplevel_procs = Toplevel_ZoltanComm.Export_Procs();
    if (num_toplevel_export>0) {
      postream << "TOPLEVEL EXPORTS:\n";
      for (int i=0; i< num_toplevel_export; ++i){
	postream << "  " << i << ":  ";
	switch (zoltanGID.Type(&export_toplevel_gids[i*ZOLTAN_GID_SIZE])) {
	case 1:
	  postream<<"NODE,      \t";
	  break;
	case 2:
	  postream<<"SHELL_NODE,\t";
	  break;
	case 3:
	  postream<<"EDGE,      \t";
	  break;
	case 4:
	  postream<<"FACE,      \t";
	  break;
	case 5:
	  postream<<"ELEM,      \t";
	  break;
	case 6:
	  postream<<"ELEMENT,   \t";
	  break;
	}
	postream<<"GID = ("<<zoltanGID.Hi(&export_toplevel_gids[i*ZOLTAN_GID_SIZE])
		<<", "<<zoltanGID.Lo(&export_toplevel_gids[i*ZOLTAN_GID_SIZE])<<")  \t"
		<<"Index = "<<zoltanLID.Index(&export_toplevel_lids[i*ZOLTAN_LID_SIZE])<<"  \t"
		<<"Proc = "<<export_toplevel_procs[i]<<"\n";
      }
    }
#endif

#ifdef CONTACT_TIMINGS
    timer.Start_Timer( secondary_connect_mesh_time );
    timer.Start_Timer( secondary_connect_mesh_phase2_time );
#endif

    ContactTopologyEntityList* entity_list = secondary_topology->NodeList();
    entity_list->BuildList(secondary_topology->Node_Blocks(),
			   secondary_topology->Number_of_Node_Blocks(),
			   no_parallel_consistency==INACTIVE);
    secondary_topology->Number_of_Nodes(entity_list->NumEntities());
  
    entity_list = secondary_topology->EdgeList();
    entity_list->BuildList(secondary_topology->Edge_Blocks(),
			   secondary_topology->Number_of_Edge_Blocks(),
			   no_parallel_consistency==INACTIVE);
    secondary_topology->Number_of_Edges(entity_list->NumEntities());
  
    entity_list = secondary_topology->FaceList();
    entity_list->BuildList(secondary_topology->Face_Blocks(),
			   secondary_topology->Number_of_Face_Blocks(),
			   no_parallel_consistency==INACTIVE);
    secondary_topology->Number_of_Faces(entity_list->NumEntities());
  
    entity_list = secondary_topology->ElemList();
    entity_list->BuildList(secondary_topology->Element_Blocks(),
			   secondary_topology->Number_of_Element_Blocks(),
			   no_parallel_consistency==INACTIVE);
    secondary_topology->Number_of_Elements(entity_list->NumEntities());
  
#ifdef CONTACT_TIMINGS
    timer.Stop_Timer( secondary_connect_mesh_phase2_time );
    timer.Start_Timer( secondary_connect_mesh_phase3_time );
#endif
  
    //============================================================
    // Setting final enitity ownership in secondary decomposition
    //============================================================
  
    int num_nodes = secondary_topology->Number_of_Nodes();
    ContactNode<Real>** secondary_nodes = 
      reinterpret_cast<ContactNode<Real>**>(secondary_topology->NodeList()->EntityList());
    for (int i=0; i<num_nodes; ++i) {
      if (secondary_nodes[i]->Secondary_Owner()==my_proc_id)  {
	secondary_nodes[i]->Ownership(ContactTopologyEntity<Real>::OWNED);
      } else {
	secondary_nodes[i]->Ownership(ContactTopologyEntity<Real>::NOT_OWNED);
	//PRECONDITION(secondary_nodes[i]->Number_Interactions(1)==0);
      }
    }
  
    int num_faces = secondary_topology->Number_of_Faces();
    ContactFace<Real>** secondary_faces = 
      reinterpret_cast<ContactFace<Real>**>(secondary_topology->FaceList()->EntityList());
    for (int i=0; i<num_faces; ++i) {
      if (secondary_faces[i]->Secondary_Owner()==my_proc_id)  {
	secondary_faces[i]->Ownership(ContactTopologyEntity<Real>::OWNED);
      } else {
	secondary_faces[i]->Ownership(ContactTopologyEntity<Real>::NOT_OWNED);
      }
    }
  
    int num_elems = secondary_topology->Number_of_Elements();
    ContactElement** secondary_elems = 
      reinterpret_cast<ContactElement**>(secondary_topology->ElemList()->EntityList());
    for (int i=0; i<num_elems; ++i) {
      if (secondary_elems[i]->Secondary_Owner()==my_proc_id)  {
	secondary_elems[i]->Ownership(ContactTopologyEntity<Real>::OWNED);
      } else {
	secondary_elems[i]->Ownership(ContactTopologyEntity<Real>::NOT_OWNED);
      }
    }
  
#ifdef CONTACT_TIMINGS
    timer.Stop_Timer( secondary_connect_mesh_phase3_time );
    timer.Start_Timer( secondary_connect_mesh_phase5_time );
#endif

    //==================================================
    // Finish making all the final topology connections
    //================================================== 

    error_code = NO_ERROR;
    if (do_node_face_search) {
      //Make all the connections for the node/entity interactions
      for (int i=0; i<num_nodes; ++i) {
        ContactNode<Real>* node = secondary_nodes[i];
        ContactNodeEntityInteraction** interactions = node->Get_NodeEntity_Interactions(1);


        for (int j=0; j<node->Number_NodeEntity_Interactions(1); ++j) {
          ContactNodeEntityInteraction* cnei = interactions[j];
	  if (cnei->Is_Tied() || cnei->Is_InfSlip() || cnei->Is_Glued()) {
            cnei->Connect_Entity(secondary_topology);
            if (cnei->Entity()==NULL) {
#ifdef CONTACT_DEBUG
              Real* position      = node->Variable(CURRENT_POSITION);
              Real* ghost_gap     = node->Variable(NODE_GHOST_GAP);
              Real* remaining_gap = node->Variable(REMAINING_GAP);
              ContactNodeFaceInteraction* cnfi = dynamic_cast<ContactNodeFaceInteraction*>(cnei);
              postream<<"In ContactSearch::Create_Search_Topology(), step "<<step_number<<"\n"
                      <<"  can't reconnect interaction " << j 
                      <<" with face "<<"("<<cnfi->FaceEntityData()->host_gid[1] <<", "<<cnfi->FaceEntityData()->host_gid[0]<<")\n"
                      <<"    node "<<node->Global_ID()<<" exodus_id:"<< node->Exodus_ID()<<"\n"
                      <<"    position      = ("<<position[0]<<", "<<position[1]<<", "<<position[2]<<")\n"
                      <<"    ghosting_gap  = ("<<ghost_gap[0]<<", "<<ghost_gap[1]<<", "<<ghost_gap[2]<<")\n"
                      <<"    remaining_gap = ("<<remaining_gap[0]<<", "<<remaining_gap[1]<<", "<<remaining_gap[2]<<")\n";
#endif
	      // Can't be found. For now, just delete it. However,
	      // this is a very bad solution which may generate 
	      // different answers in parallel than in serial.
	      // A better solution is to make sure it never happens.
	      // The code is here to make sure it doesn't happen,
	      // but commented out for now
   
	      static int printlimit = 0;
	      if (printlimit < 5) {
		ContactNodeFaceInteraction* cnfi = 
		  dynamic_cast<ContactNodeFaceInteraction*>(cnei);
		std::cout 
		  <<"In ContactSearch::Create_Search_Topology(), step "
		  <<step_number<<"\n"
		  <<"  can't reconnect interaction " << j 
		  <<" with face "<<"("<<cnfi->FaceEntityData()->host_gid[1] 
		  <<", "<<cnfi->FaceEntityData()->host_gid[0]<<")\n"
		  <<"    node "<<node->Global_ID()<<" exodus_id:"
		  << node->Exodus_ID()<<"\n";
		printlimit ++;
	      }
              node->Delete_NodeEntity_Interaction(cnei,1);
	      //   uncomment this to do the right thing, i.e. die
	      //   if we ever reach this condition
              //error_code = INTERNAL_ERROR;
	    }
          }
 	}
      }
    }
    ContactErrorCode err = 
      (ContactErrorCode) contact_global_error_check(error_code, SearchComm);
    if (err != NO_ERROR) {
      errors->Add_Error_Message("Internal Error! Contact a developer.");
      errors->Add_Error_Message("Failed to connect a tied, glued, or infinitesimal");
      errors->Add_Error_Message("slip interaction in the primary, probably due");
      errors->Add_Error_Message("to a missing ghosted entity.");
      error_code = err;
      return;
    }
	

#ifdef CONTACT_TIMINGS
    timer.Stop_Timer( secondary_connect_mesh_phase5_time );
    timer.Start_Timer( secondary_connect_mesh_phase4_time );
#endif
    
    secondary_topology->Connect_Nodes_to_Faces();
    secondary_topology->Connect_Nodes_to_Elements();
    secondary_topology->Connect_Faces_to_Nodes();
  
    for(int i=0; i<primary_topology->Number_of_Node_Blocks(); ++i) {
      int flag = primary_topology->Node_Block(i)->Has_Attributes();
      secondary_topology->Node_Block(i)->Has_Attributes(flag);
    }
  
#if CONTACT_DEBUG_PRINT_LEVEL>=9
    int num_secondary_nfi = 0;
    int num_secondary_nsi = 0;
    int num_secondary_ffi = 0;
    int num_secondary_eei = 0;
    for (int i=0; i<num_nodes; ++i) {
      num_secondary_nfi += secondary_nodes[i]->Number_NodeFace_Interactions(1);
      num_secondary_nsi += secondary_nodes[i]->Number_NodeSurface_Interactions(1);
    }
#endif
  
#ifdef CONTACT_TIMINGS
    timer.Stop_Timer( secondary_connect_mesh_phase4_time );
    timer.Stop_Timer( secondary_connect_mesh_time );
#endif
  
    // DEBUG -- output secondary decomposition mesh. 
    /*...
      int  total_procs = contact_number_of_processors(SearchComm);
      char FileName[1024];
      char temp[81];
      std::strcpy(FileName,"Secondary_Mesh.exo");
      std::strcat(FileName,".");
      std::sprintf(temp,"%d",total_procs);
      std::strcat(FileName,temp);
      std::strcat(FileName,".");
      std::sprintf(temp,"%d",my_proc_id);
      std::strcat(FileName,temp);
      int compws = 8;
      int iows = 8;
      int exodus_id = ex_create( FileName, EX_CLOBBER, &compws, &iows );
      contact_global_sync(SearchComm);
      ...*/
    zoltan->Free_Data();
  
#if CONTACT_DEBUG_PRINT_LEVEL>=9
    secondary_topology->Display_Entities(postream);
#endif
//FIXIT #if CONTACT_DEBUG_PRINT_LEVEL>=2
#if CONTACT_DEBUG_PRINT_LEVEL>=2
    //FIXIT postream << "    Search Topology is Constructed\n";
    postream << "    Step: "<<step_number<<", Search Topology is Constructed\n";
    postream << "      Number of Nodes    = " 
	     << secondary_topology->Number_of_Nodes() 
	     << "   (copied " << copy_nodes
	     << ")\n";
#if CONTACT_DEBUG_PRINT_LEVEL>=9
    for (int i=0; i<num_nodes; ++i) {
      postream << "      ID = " << secondary_nodes[i]->Global_ID() 
	       << "    \tSecondaryOwner = " << secondary_nodes[i]->Secondary_Owner()
	       << "\n";
    }
#endif
    postream << "      Number of Faces    = " 
	     << secondary_topology->Number_of_Faces() 
	     << "   (copied " << copy_faces
	     << ")\n";
#if CONTACT_DEBUG_PRINT_LEVEL>=9
    for (int i=0; i<num_faces; ++i) {
      postream << "      ID = " << secondary_faces[i]->Global_ID() 
	       << "    \tSecondaryOwner = " << secondary_faces[i]->Secondary_Owner()
	       << "\n";
    }
#endif
    postream << "      Number of Elements = " 
	     << secondary_topology->Number_of_Elements() 
	     << "   (copied " << copy_elems
	     << ")\n";
    postream << "      Number of Analytic Surfaces = " 
	     << secondary_topology->Number_of_Analytic_Surfaces() << "\n";
#if CONTACT_DEBUG_PRINT_LEVEL>=9
    postream << "      Number of Node/Face Interactions = " 
	     << num_secondary_nfi << "\n";
    postream << "      Number of Node/Surface Interactions = " 
	     << num_secondary_nsi << "\n";
    postream << "      Number of Face/Face Interactions = " 
	     << num_secondary_ffi << "\n";
    postream << "      Number of Element/Element Interactions = " 
	     << num_secondary_eei << "\n";
#endif
    postream.flush();
#endif

    search_topology = secondary_topology;
  
  } else {
    if(no_ghosting == 0) {
    if(do_node_face_search &&
       !do_face_face_search &&
       !do_node_node_search &&
       !do_coverage_search &&
       !do_elem_elem_search) {
      primary_topology->DoGhosting_New_NodeFace(POSITION, reasonable_gap);
    } else {  
      primary_topology->DoGhosting(POSITION, reasonable_gap);
    }}
    search_topology = primary_topology;  
  }
  
#if !defined(CONTACT_NO_MPI) && defined (CONTACT_TIMINGS)
  timer.Stop_Timer( create_secondary_time );
#endif

#endif
  return;
}


ContactSearch::ContactErrorCode ContactSearch::Define_Primary_Interactions()
{
#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS) && defined(CONTACT_TIMINGS1)
  BaselineTiming();
#endif
  if( contact_number_of_processors(SearchComm)==1 || 
      no_secondary==ACTIVE) {
    return NO_ERROR;
  
#ifndef CONTACT_NO_MPI
  } else {
#ifdef CONTACT_TIMINGS
    timer.Start_Timer( interaction_migration_time );
    timer.Start_Timer( interaction_onproc_copy_time );
#endif

    int my_proc_id = contact_processor_number(SearchComm);
    
    ContactTopologyEntityList* primary_node_list = primary_topology->NodeList();
    ContactTopologyEntityList* primary_face_list = primary_topology->FaceList();
    ContactTopologyEntityList* primary_elem_list = primary_topology->ElemList();
    ContactNode<Real>** PrimaryNodes = 
      reinterpret_cast<ContactNode<Real>**>(primary_node_list->EntityList());
    ContactFace<Real>** PrimaryFaces = 
      reinterpret_cast<ContactFace<Real>**>(primary_face_list->EntityList());
    ContactElement** PrimaryElements = 
      reinterpret_cast<ContactElement**>(primary_elem_list->EntityList());
    
    int num_secondary_nodes    = secondary_topology->Number_of_Nodes();
    int num_secondary_faces    = secondary_topology->Number_of_Faces();
    int num_secondary_elements = secondary_topology->Number_of_Elements();
    ContactTopologyEntityList* secondary_node_list = secondary_topology->NodeList();
    ContactTopologyEntityList* secondary_face_list = secondary_topology->FaceList();
    ContactTopologyEntityList* secondary_elem_list = secondary_topology->ElemList();
    ContactNode<Real>** SecondaryNodes = 
      reinterpret_cast<ContactNode<Real>**>(secondary_node_list->EntityList());
    ContactFace<Real>** SecondaryFaces = 
      reinterpret_cast<ContactFace<Real>**>(secondary_face_list->EntityList());
    ContactElement** SecondaryElements = 
      reinterpret_cast<ContactElement**>(secondary_elem_list->EntityList());

    ContactNodeNodeInteraction*       cnni;
    ContactFaceFaceInteraction<Real>* cffi;
    ContactElementElementInteraction* ceei;
    ContactInteractionEntity<Real>* interaction;

    // copy over any node/node, node/face and node/surface interactions from 
    // the secondary to primary topology that are on the same processor
    int num_node_export = 0;
    for (int n=0; n<num_secondary_nodes; ++n) {
      ContactNode<Real>* secondary_node = SecondaryNodes[n];
      
#ifdef CONTACT_DEBUG_NODE
      bool PRINT_THIS_NODE = primary_topology->Is_a_Debug_Node( secondary_node );
      if( PRINT_THIS_NODE ){
        postream << "  Start of Define_Primary(), Node "
                 << secondary_node->Exodus_ID() 
                 <<" has "<<secondary_node->Number_Interactions()
                 <<" interactions\n";
      }
#endif
      
      //if (secondary_node->StatusFlag()!=ContactTopologyEntity<Real>::GLOBAL_SEARCH) continue;
      if (secondary_node->Number_Interactions()==0) continue;
      if (secondary_node->Ownership() != ContactTopologyEntity<Real>::OWNED) continue;
      if (secondary_node->Owner()==my_proc_id) {
        ContactNode<Real>* primary_node = 
          PrimaryNodes[secondary_node->OwnerProcArrayIndex()];
        POSTCONDITION(primary_node);
        if (secondary_node->Number_Interactions()>0) {
          primary_node->Copy_Interactions(secondary_node, STATE_0);
          ContactInteractionDLL<Real>* nn_interactions = primary_node->Get_NodeNode_Interactions();
          if(nn_interactions != NULL) {
            nn_interactions->IteratorStart();
            while ((interaction=nn_interactions->IteratorForward())){
              cnni = static_cast<ContactNodeNodeInteraction*>(interaction);
              cnni->Connect_MasterNode( *primary_node_list );
            }
	  }
          ContactNodeEntityInteraction** ne_interactions = 
            primary_node->Get_NodeEntity_Interactions();
          for (int i=0; i<primary_node->Number_NodeEntity_Interactions(); ++i) {
            ContactNodeEntityInteraction *cnei = ne_interactions[i];
            cnei->Connect_Entity(primary_topology);
            if(cnei->Get_Entity_Owner()==my_proc_id) {
              POSTCONDITION(cnei->Entity());
            }
            POSTCONDITION(cnei->Node());
          }
        }
        unsigned int context = secondary_node->GetContext() && 
                               (ContactTopologyEntity<Real>::TRACK_SEARCH_SLAVE | ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE);
        primary_node->SetContextBit(context);
      } else {
        ++num_node_export;
      }
    }
    
    // copy over any face/face interactions from the secondary
    // to primary topology that are on the same processor
    int num_face_export = 0;
    for (int n=0; n<num_secondary_faces; ++n) {
      ContactFace<Real>* secondary_face = SecondaryFaces[n];
      if (secondary_face->Number_Interactions()==0) continue;
      if (secondary_face->Owner()==my_proc_id) {
        ContactFace<Real>* primary_face = 
	  PrimaryFaces[secondary_face->OwnerProcArrayIndex()];
        POSTCONDITION(primary_face);
        primary_face->Copy_Interactions(secondary_face);
        ContactInteractionDLL<Real>* interactions = primary_face->Get_FaceFace_Interactions();
        if(interactions == NULL) continue;
        interactions->IteratorStart();
        while ((interaction=interactions->IteratorForward())){
          cffi = static_cast<ContactFaceFaceInteraction<Real>*>(interaction);
          cffi->Connect_MasterFace( *primary_face_list );
        }
      } else {
        ++num_face_export;
      }
    }
    
    // copy over any element/element interactions from the secondary
    // to primary topology that are on the same processor
    int num_elem_export = 0;
    for (int n=0; n<num_secondary_elements; ++n) {
      ContactElement* secondary_element = SecondaryElements[n];
      if (secondary_element->Number_Interactions()==0) continue;
      if (secondary_element->Owner()==my_proc_id) {
        ContactElement* primary_element = 
	  PrimaryElements[secondary_element->OwnerProcArrayIndex()];
        POSTCONDITION(primary_element);
        primary_element->Copy_Interactions(secondary_element);
        ContactInteractionDLL<Real>* interactions = primary_element->Get_ElementElement_Interactions();
        interactions->IteratorStart();
        while ((interaction=interactions->IteratorForward())){
          ceei = static_cast<ContactElementElementInteraction*>(interaction);
          ceei->Connect_MasterElement( *primary_elem_list );
        }
      } else {
        ++num_elem_export;
      }
    }

    int num_exports = num_node_export+num_face_export+num_elem_export;
    int num_total_exports = contact_global_maximum(num_exports, SearchComm);
    
#ifdef CONTACT_TIMINGS
    timer.Stop_Timer( interaction_onproc_copy_time );
#endif
    
    if (num_total_exports>0) {
    
#ifdef CONTACT_TIMINGS
      timer.Start_Timer( interaction_export_setup_time );
#endif
    
      LB_ID_TYPE zoltan_lid[ZOLTAN_LID_SIZE];
      LB_ID_TYPE zoltan_gid[ZOLTAN_GID_SIZE];
      int zoltan_pid;
      ContactZoltanComm InteractionZoltanComm(ContactZoltanComm::ZOLTAN_EXPORT);
    
      // now handle all the node/node, node/face and node/surface
      // interactions that have to go off processor
      for (int n=0; n<num_secondary_nodes; ++n) {
	ContactNode<Real>* secondary_node = SecondaryNodes[n];
	//if (secondary_node->StatusFlag()!=ContactTopologyEntity<Real>::GLOBAL_SEARCH) continue;
	if (secondary_node->Number_Interactions()==0) continue;
	if (secondary_node->Ownership() != ContactTopologyEntity<Real>::OWNED) continue;
	if (secondary_node->Owner()!=my_proc_id) {
	  zoltan_pid = secondary_node->Owner();
	  secondary_node->ZoltanLID(CT_NODE, zoltan_lid);
	  secondary_node->ZoltanGID(CT_NODE, zoltan_gid);
	  InteractionZoltanComm.Add_Export(zoltan_lid,zoltan_gid,zoltan_pid,0);
	}
      }
    
      // now handle all the face/face interactions that have to go off processor
      for (int n=0; n<num_secondary_faces; ++n) {
	ContactFace<Real>* secondary_face = SecondaryFaces[n];
	if (secondary_face->Owner()!=my_proc_id &&
	    secondary_face->Number_Interactions()>0) {
	  zoltan_pid = secondary_face->Owner();
	  secondary_face->ZoltanLID(CT_FACE, zoltan_lid);
	  secondary_face->ZoltanGID(CT_FACE, zoltan_gid);
	  InteractionZoltanComm.Add_Export(zoltan_lid,zoltan_gid,zoltan_pid,0);
	}
      }
    
      // now handle all the element/element interactions that have to go off processor
      for (int n=0; n<num_secondary_elements; ++n) {
	ContactElement* secondary_element = SecondaryElements[n];
	if (secondary_element->Owner()!=my_proc_id &&
	    secondary_element->Number_Interactions()>0) {
	  zoltan_pid = secondary_element->Owner();
	  secondary_element->ZoltanLID(CT_ELEMENT, zoltan_lid);
	  secondary_element->ZoltanGID(CT_ELEMENT, zoltan_gid);
	  InteractionZoltanComm.Add_Export(zoltan_lid,zoltan_gid,zoltan_pid,0);
	}
      }
    
      int       num_import  = -1;
      LB_ID_PTR import_gids  = NULL;
      LB_ID_PTR import_lids  = NULL;
      int*      import_procs = NULL;
      int       num_export   = num_node_export+
	num_face_export+
	num_elem_export;
      LB_ID_PTR export_gids  = InteractionZoltanComm.Export_GIDS();
      LB_ID_PTR export_lids  = InteractionZoltanComm.Export_LIDS();
      int*      export_procs = InteractionZoltanComm.Export_Procs();
      zoltan->Set_InteractionCallBacks(this);
    
#if CONTACT_DEBUG_PRINT_LEVEL>=6
      postream<<"Number of interactions exported = "<<num_export<<"\n";
      postream.flush();
    
      int global_num_export = contact_global_sum(num_export, SearchComm);
      if (contact_processor_number(SearchComm)==0) {
	std::cout<<std::endl;
	std::cout<<"Global number of interactions exported = "<<global_num_export<<std::endl;
	std::cout<<std::endl<<std::flush;
      }
      if (num_export>0) {
	postream << "INTERACTION EXPORTS:\n";
	for (int i=0; i< num_export; ++i){
	  ContactZoltanLID zoltanLID;
	  ContactZoltanGID zoltanGID;
	  postream << "  " << i << ":  ";
	  switch (zoltanGID.Type(&export_gids[i*ZOLTAN_GID_SIZE])) {
	  case 1:
	    postream<<"NODE,      \t";
	    break;
	  case 2:
	    postream<<"SHELL_NODE,\t";
	    break;
	  case 3:
	    postream<<"EDGE,      \t";
	    break;
	  case 4:
	    postream<<"FACE,      \t";
	    break;
	  case 5:
	    postream<<"ELEM,      \t";
	    break;
	  case 6:
	    postream<<"ELEMENT,   \t";
	    break;
	  }
	  postream<<"GID = ("<<zoltanGID.Hi(&export_gids[i*ZOLTAN_GID_SIZE])
		  <<", "<<zoltanGID.Lo(&export_gids[i*ZOLTAN_GID_SIZE])<<")  \t"
		  <<"Index = "<<zoltanLID.Index(&export_lids[i*ZOLTAN_LID_SIZE])<<"  \t"
		  <<"Proc = "<<export_procs[i]<<"\n";
	}
      }
      postream.flush();
#endif

#ifdef CONTACT_TIMINGS
      timer.Stop_Timer( interaction_export_setup_time );
      timer.Start_Timer( interaction_help_migrate_time );
#endif
#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
      bytes_nodes = 0.0;
      bytes_faces = 0.0;
      bytes_elems = 0.0;
#endif
      zoltan->Help_Migrate(num_import,  import_gids, 
			   import_lids, import_procs,
			   num_export,  export_gids,
			   export_lids, export_procs);
#if CONTACT_DEBUG_PRINT_LEVEL>=3 || defined(CONTACT_ANALYZE_DATA_XFER)
      bytes_nodes /= 1024;
      bytes_faces /= 1024;
      bytes_elems /= 1024;
      postream<<"    Data Sent by Zoltan when migrating interactions...\n";
      postream<<"      nodes:  "<<bytes_nodes<<" Kb\n";
      postream<<"      faces:  "<<bytes_faces<<" Kb\n";
      postream<<"      elems:  "<<bytes_elems<<" Kb\n";
#endif
#ifdef CONTACT_TIMINGS
      timer.Stop_Timer( interaction_help_migrate_time );
      timer.Stop_Timer( interaction_migration_time );
#endif
    }
    
#if CONTACT_DEBUG_PRINT_LEVEL>=6
    int num_nni = primary_topology->Number_NodeNode_Interactions();
    int num_nfi = primary_topology->Number_NodeFace_Interactions();
    int num_nsi = primary_topology->Number_NodeSurface_Interactions();
    postream<<"Number of interactions in primary = "<<num_nni+num_nfi+num_nsi<<"\n";
    postream.flush();
    
    int global_num = contact_global_sum(num_nni+num_nfi+num_nsi, SearchComm);
    if (contact_processor_number(SearchComm)==0) {
      std::cout<<std::endl;
      std::cout<<"Global number of interactions in primary = "<<global_num<<std::endl;
      std::cout<<std::endl<<std::flush;
    }
#endif
#ifdef CONTACT_DEBUG_NODE
    int num_primary_nodes = primary_topology->Number_of_Nodes();
    for (int n=0; n<num_primary_nodes; ++n) {
      ContactNode<Real>* primary_node = PrimaryNodes[n];
      bool PRINT_THIS_NODE = primary_topology->Is_a_Debug_Node( primary_node );
      if( PRINT_THIS_NODE ){
        postream << "  End of Define_Primary(), Node "
                 << primary_node->Exodus_ID() 
                 <<" has "<<primary_node->Number_Interactions()
                 <<" interactions\n";
      }
    }
#endif
#endif
  }
  return NO_ERROR;
}

#if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS) && defined(CONTACT_TIMINGS1)
void ContactSearch::BaselineTiming()
{
  int n;
  int size_buf = 0;
  char* buffer0 = new char[10000];
  char* buffer1 = new char[10000];
  char* buffer2 = new char[10000];
  ContactTopologyEntityList* search_node_list = 
    search_topology->NodeList();
  ContactTopologyEntityList* search_edge_list = 
    search_topology->EdgeList();
  ContactTopologyEntityList* search_face_list = 
    search_topology->FaceList();
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(search_node_list->EntityList());
  ContactEdge<Real>** Edges = 
    reinterpret_cast<ContactEdge<Real>**>(search_edge_list->EntityList());
  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(search_face_list->EntityList());
  ContactNode<Real>** node_list = 
    new ContactNode<Real>*[search_node_list->NumEntities()];
  ContactNode<Real>** node_list0 = 
    new ContactNode<Real>*[search_node_list->NumEntities()];
  ContactNode<Real>** node_list1 = 
    new ContactNode<Real>*[search_node_list->NumEntities()];
  ContactEdge<Real>** edge_list = 
    new ContactEdge<Real>*[search_edge_list->NumEntities()];
  ContactEdge<Real>** edge_list0 = 
    new ContactEdge<Real>*[search_edge_list->NumEntities()];
  ContactFace<Real>** face_list = 
    new ContactFace<Real>*[search_face_list->NumEntities()];
  ContactFace<Real>** face_list0 = 
    new ContactFace<Real>*[search_face_list->NumEntities()];
  int nnodes = search_node_list->NumEntities();
  int nedges = search_edge_list->NumEntities();
  int nfaces = search_face_list->NumEntities();
    
    
    
  //=====================================================================
  // Constructor Timing
  //=====================================================================
  
  timer.Start_Timer( baseline_constructors_time );
  
  timer.Start_Timer( baseline_node_time );
  for (int n=0; n<nnodes; ++n) {
    ContactNode<Real>* new_node = ContactNode<Real>::new_ContactNode(allocators, NODE );
    node_list[n] = new_node;
  }
  timer.Stop_Timer( baseline_node_time );
      
  timer.Start_Timer( baseline_edge_time );
  for (int n=0; n<nedges; ++n) {
    ContactEdge<Real>* new_edge;
    switch(Edges[n]->EdgeType()){
    case ContactSearch::LINEEDGEL2 :
      new_edge = ContactLineEdgeL2<Real>::new_ContactLineEdgeL2(
							  allocators[ALLOC_ContactLineEdgeL2] );
      break;
    case ContactSearch::LINEEDGEQ3 :
      new_edge = ContactLineEdgeQ3::new_ContactLineEdgeQ3(
							  allocators[ALLOC_ContactLineEdgeQ3] );
      break;
    }
    edge_list[n] = new_edge;
  }
  timer.Stop_Timer( baseline_edge_time );  
      
  timer.Start_Timer( baseline_face_time );
  for (int n=0; n<nfaces; ++n) {
    ContactFace<Real>* face = Faces[n];
    ContactFace<Real>* new_face = New_ContactFace(face->FaceType());
    face_list[n] = new_face;
  }
  timer.Stop_Timer( baseline_face_time );  
  
  timer.Stop_Timer( baseline_constructors_time );
    
  //=====================================================================
  // Size / Pack / Unpack Timing
  //=====================================================================
  
  int node_size = 0;
  int edge_size = 0;
  int face_size = 0;
  timer.Start_Timer( baseline_size_time );
  for (int n=0; n<nnodes; ++n) {
    int sz     = Nodes[n]->Size();
    node_size += sz;
    size_buf   = std::max(size_buf,sz);
  }
  for (int n=0; n<nedges; ++n) {
    int sz     = Edges[n]->Size();
    edge_size += sz;
    size_buf   = std::max(size_buf,sz);
  }
  for (int n=0; n<nfaces; ++n) {
    int sz     = Faces[n]->Size();
    face_size += sz;
    size_buf   = std::max(size_buf,sz);
  }
  timer.Stop_Timer( baseline_size_time );

  int sums_input[3];
  int sums_output[3];
  sums_input[0] = node_sum;
  sums_input[1] = edge_sum;
  sums_input[2] = face_sum;
  contact_global_sum(sums_input, sums_output, 3, SearchComm);
  int node_sum = sums_output[0];
  int edge_sum = sums_output[1];
  int face_sum = sums_output[2];

  if (contact_processor_number(SearchComm)==0) {
    std::cout<<"Total node buffer size = "<<node_sum<<std::endl;
    std::cout<<"Total edge buffer size = "<<edge_sum<<std::endl;
    std::cout<<"Total face buffer size = "<<face_sum<<std::endl;
    std::cout<<std::flush;
  }
 
  timer.Start_Timer( baseline_pack_time );
  for (int n=0; n<nnodes; ++n) {
    Nodes[n]->Pack( buffer0 );
  }
  for (int n=0; n<nedges; ++n) {
    Edges[n]->Pack( buffer1 );
  }
  for (int n=0; n<nfaces; ++n) {
    Faces[n]->Pack( buffer2 );
  }
  timer.Stop_Timer( baseline_pack_time );
   
  timer.Start_Timer( baseline_unpack_time );
  for (int n=0; n<nnodes; ++n) {
    ContactNode<Real>* new_node = ContactNode<Real>::new_ContactNode(allocators, NODE );
    node_list0[n] = new_node;
    new_node->Unpack( buffer0 );
  }
  for (int n=0; n<nedges; ++n) {
    ContactEdge<Real>* new_edge;
    switch(Edges[n]->EdgeType()){
    case ContactSearch::LINEEDGEL2 :
      new_edge = ContactLineEdgeL2<Real>::new_ContactLineEdgeL2(
							  allocators[ALLOC_ContactLineEdgeL2] );
      break;
    case ContactSearch::LINEEDGEQ3 :
      new_edge = ContactLineEdgeQ3::new_ContactLineEdgeQ3(
							  allocators[ALLOC_ContactLineEdgeQ3] );
      break;
    }
    edge_list0[n] = new_edge;
    new_edge->Unpack( buffer1 );
  }
  for (int n=0; n<nfaces; ++n) {
    ContactFace<Real>* face = Faces[n];
    ContactFace<Real>* new_face = New_ContactFace(face->FaceType());
    face_list0[n] = new_face;
    new_face->Unpack( buffer2 );
  }
  timer.Stop_Timer( baseline_unpack_time );
    
  //=====================================================================
  // Size / Pack / Unpack All Timing
  //=====================================================================
 
  timer.Start_Timer( baseline_size_all_time );
  for (int n=0; n<nnodes; ++n) {
    size_buf = std::max(size_buf,Nodes[n]->Size(0));
  }
  timer.Stop_Timer( baseline_size_all_time );

  timer.Start_Timer( baseline_pack_all_time );
  for (int n=0; n<nnodes; ++n) {
    Nodes[n]->Pack( buffer0, 0 );
  }
  timer.Stop_Timer( baseline_pack_all_time );

  timer.Start_Timer( baseline_unpack_all_time );
  for (int n=0; n<nnodes; ++n) {
    ContactNode<Real>* new_node = ContactNode<Real>::new_ContactNode(allocators, NODE );
    node_list1[n] = new_node;
    new_node->Unpack( buffer0 );
  }
  timer.Stop_Timer( baseline_unpack_all_time );
}
#endif

#ifndef CONTACT_NO_MPI
void ContactTopology::Connect_Nodes_to_Elements( )
{
  ContactElement** Elements = 
    reinterpret_cast<ContactElement**>(elem_list->EntityList());
  for (int i=0; i<number_of_elements; ++i) {
    int num_nodes = Elements[i]->Nodes_Per_Element();
    ContactTopologyEntity<Real>::connection_data *node_info = Elements[i]->NodeInfo();
    POSTCONDITION( node_info );
    for( int j=0 ; j<num_nodes ; ++j ){
      ContactNode<Real>* node = static_cast<ContactNode<Real> *>
	(node_list->Find( &node_info[j] ));
      POSTCONDITION( node );
      Elements[i]->ConnectNode( j,node );
    }
  }
}
#endif

#ifndef CONTACT_NO_MPI
void ContactTopology::Connect_Nodes_to_Faces( )
{
  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(face_list->EntityList());
  for (int i=0; i<number_of_faces; ++i) {
    int num_nodes = Faces[i]->Nodes_Per_Face();
    ContactTopologyEntity<Real>::connection_data *node_info = Faces[i]->NodeInfo();
    POSTCONDITION( node_info );
    for( int j=0 ; j<num_nodes ; ++j ){
      ContactNode<Real>* node = static_cast<ContactNode<Real> *>
	(node_list->Find( &node_info[j] ));
      POSTCONDITION( node );
      Faces[i]->ConnectNode( j,node );
    }
  }
}
#endif

#ifndef CONTACT_NO_MPI

void ContactSearch::Create_Processor_Bounding_Boxes(ObjectBoundingBox *proc_box_array, 
                                                    int &POSITION,
                                                    const int &total_num_procs) {
  //
  //  Save the processor number that each box refers to
  //
  int mwg_cnt = 0;
  for(int iproc = 0; iproc < total_num_procs; ++iproc) {
    proc_box_array[iproc].set_object_number(iproc);
  }
  
  VariableHandle NODE_GHOST_GAP   = primary_topology->Variable_Handle( ContactTopology::Node_Ghost_Gap );
  VariableHandle CURRENT_POSITION = primary_topology->Variable_Handle( ContactTopology::Current_Position );
  int            num_config       = (POSITION!=CURRENT_POSITION)?2:1;
  
  //
  //  Expand the processor bounding boxes to encompase all owned nodes
  //  faces, and elements
  //
  
  if (do_node_face_search || do_node_node_search) {
  VariableHandle NODE_RADIUS = 
    primary_topology->Variable_Handle( ContactTopology::Node_Radius );
  ContactTopologyEntityList* primary_node_list = primary_topology->NodeList();
  for(int i=0; i<primary_topology->Number_of_Node_Blocks(); ++i) {
    ContactNodeBlock* node_block = primary_topology->Node_Block(i);
    //if (node_block->IsSlave()) {
      int nnodes = primary_node_list->BlockNumEntities(i);
      ContactNode<Real>** nodes  = 
        reinterpret_cast<ContactNode<Real>**>(primary_node_list->BlockEntityList(i));
      if (node_block->Type() == POINT) {
        if (node_block->Has_Radius_Attributes()) {
          for (int j=0; j<nnodes; ++j) {
            ContactNode<Real>* node = nodes[j];
            if (node->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)) {
              int secondary_owner = node->Secondary_Owner();
              ObjectBoundingBox *proc_box = &(proc_box_array[secondary_owner]);
              Real* position =  node->Variable(POSITION);
              Real  radius   = *node->Variable(NODE_RADIUS);
              proc_box->add_sphere(position, radius);
              if (POSITION != CURRENT_POSITION) {
                position = node->Variable(CURRENT_POSITION);
                proc_box->add_sphere(position, radius);
                mwg_cnt++;
              }
            }
          }
        } else {
          for (int j=0; j<nnodes; ++j) {
            ContactNode<Real>* node = nodes[j];
            if (node->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)) {
              int secondary_owner = node->Secondary_Owner();
              ObjectBoundingBox *proc_box = &(proc_box_array[secondary_owner]);
              Real* position = node->Variable(POSITION);
              proc_box->add_point(position);
              if (POSITION != CURRENT_POSITION) {
                position = node->Variable(CURRENT_POSITION);
                proc_box->add_point(position);
                mwg_cnt++;
              }
            }
          }
        }
      } else {
        for (int j=0; j<nnodes; ++j) {
          ContactNode<Real>* node = nodes[j];
          if (node->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)) {
            int secondary_owner = node->Secondary_Owner();
            ObjectBoundingBox *proc_box = &(proc_box_array[secondary_owner]);
            ContactBoundingBox node_box;
            node->ComputeBoundingBoxForSearch(num_config,
                                              NODE_GHOST_GAP,
                                              CURRENT_POSITION,
                                              POSITION,
                                              auto_tol,
                                              box_inflation,
                                              node_box);
            proc_box->add_box(node_box);
                mwg_cnt++;
          }
        }
      }
    //}
  }
  }
    
  if (do_face_face_search || do_coverage_search) {
  int num_primary_faces = primary_topology->Number_of_Faces();
  ContactFace<Real>**    primary_faces = reinterpret_cast<ContactFace<Real>**>   (primary_topology->FaceList()->EntityList());
  for (int i=0; i<num_primary_faces; ++i) {
    ContactFace<Real>* face = primary_faces[i];
    if (face->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)) {
      int secondary_owner = face->Secondary_Owner();
      ObjectBoundingBox *proc_box = &(proc_box_array[secondary_owner]);
      ContactBoundingBox object_box;
      int num_nodes = face->Nodes_Per_Face();
      for (int j=0; j<num_nodes; ++j) {
        ContactNode<Real> * node = face->Node(j);
        object_box.add_point(node->Variable(POSITION));
        object_box.add_point(node->Variable(CURRENT_POSITION));
      }
      proc_box->add_box(object_box);
    }
  }
  }
  
  if (do_elem_elem_search) {
  int num_primary_elems = primary_topology->Number_of_Elements();
  ContactElement** primary_elems = reinterpret_cast<ContactElement**>(primary_topology->ElemList()->EntityList());
  for (int i=0; i<num_primary_elems; ++i) {
    ContactElement* elem = primary_elems[i];
    if (elem->CheckContext(ContactTopologyEntity<Real>::GLOBAL_SEARCH_SLAVE)) {
      int secondary_owner = elem->Secondary_Owner();
      ObjectBoundingBox *proc_box = &(proc_box_array[secondary_owner]);

      ContactBoundingBox object_box;
      int num_nodes = elem->Nodes_Per_Element();
      for (int j=0; j<num_nodes; ++j) {
        ContactNode<Real> * node = elem->Node(j);
        object_box.add_point(node->Variable(POSITION));
        object_box.add_point(node->Variable(CURRENT_POSITION));
      }
      proc_box->add_box(object_box);
    }
  }
  }

  //
  //  Perform global reductions to get the true global processor bounding boxes 
  //  for the secondary decomposition
  //
//FIXIT#if CONTACT_DEBUG_PRINT_LEVEL>=3
#if CONTACT_DEBUG_PRINT_LEVEL>=3
  postream<<"    Compute bounding boxes for "<<mwg_cnt<<" nodes\n";
  for (int i=0; i<total_num_procs; ++i) {
    ObjectBoundingBox *proc_box = &(proc_box_array[i]);
    postream<<"    Local Processor "<<i<<" bounding box\n";
    postream<<"      min("<<proc_box->get_x_min()<<", "<<proc_box->get_y_min()<<", "<<proc_box->get_z_min()<<")\n";
    postream<<"      max("<<proc_box->get_x_max()<<", "<<proc_box->get_y_max()<<", "<<proc_box->get_z_max()<<")\n";
  }
#endif
  ObjectBoundingBox::global_box_combine(proc_box_array, total_num_procs, SearchComm);
//FIXIT #if CONTACT_DEBUG_PRINT_LEVEL>=3
#if CONTACT_DEBUG_PRINT_LEVEL>=3
  for (int i=0; i<total_num_procs; ++i) {
    ObjectBoundingBox *proc_box = &(proc_box_array[i]);
    postream<<"    Global Processor "<<i<<" bounding box\n";
    postream<<"      min("<<proc_box->get_x_min()<<", "<<proc_box->get_y_min()<<", "<<proc_box->get_z_min()<<")\n";
    postream<<"      max("<<proc_box->get_x_max()<<", "<<proc_box->get_y_max()<<", "<<proc_box->get_z_max()<<")\n";
  }
#endif
}

#endif

