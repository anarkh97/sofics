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

#include "ContactUtilities.h"
using acme::Normalize;
#include "Contact_Defines.h"
#include "ContactSearch.h"
#include "ContactTopology.h"
#include "ContactNode.h"
#include "ContactEdge.h"
#include "ContactFace.h"
#include "ContactFaceBlock.h"
#include "ContactEdgeBlock.h"
#include "Contact_Communication.h"
#include "ContactParOStream.h"
#include "ContactSymComm.h"
#include "ContactErrors.h"
#include "ContactShellHandler.h"
#include "ContactSearchData.h"
#include "ContactParOStream.h"
#include "ContactTopology.h"

#include <iostream>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>

using namespace std;

void ContactTopology::Compute_Characteristic_Length( )
{
  PRECONDITION( topology_type == PRIMARY );
  // Have each face compute its characteristic length
  if ( !computed_characteristic_length ) {
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    search->Timer()->Start_Timer( search->face_charlen_geom_time );
    #endif
    ContactFace<Real>** Faces = 
      reinterpret_cast<ContactFace<Real>**>(primary_face_list->EntityList());
    min_characteristic_length =  BIGNUM;
    max_characteristic_length = -BIGNUM;
    for (int i=0; i<number_of_primary_faces; ++i) {
      ContactFace<Real>* face = Faces[i];
      face->Compute_CharacteristicLength(Current_Position, CHARACTERISTIC_LENGTH);
      Real* characteristic_length = face->Variable(CHARACTERISTIC_LENGTH);
      min_characteristic_length = std::min(min_characteristic_length,*characteristic_length);
      max_characteristic_length = std::max(max_characteristic_length,*characteristic_length);
    }
    computed_characteristic_length = true;
    Real temp = contact_global_minimum( min_characteristic_length, SearchComm );
    min_characteristic_length = temp;
    temp = contact_global_maximum( max_characteristic_length, SearchComm );
    max_characteristic_length = temp;
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    search->Timer()->Stop_Timer( search->face_charlen_geom_time );
    #endif
  }
}

void ContactTopology::Compute_Element_Geometry( VariableHandle POSITION )
{
  PRECONDITION( topology_type == PRIMARY );

  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  search->Timer()->Start_Timer( search->elem_geom_time );
  #endif
  // compute the element volumes and centroid
  ContactElement** Elements = 
    reinterpret_cast<ContactElement**>(primary_elem_list->EntityList());
  for (int i=0; i<number_of_primary_elements; ++i) {
    Elements[i]->Compute_Volume(POSITION, ELEMENT_VOLUME);
    Elements[i]->Compute_Centroid(POSITION, ELEMENT_CENTROID);
  }
  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  search->Timer()->Stop_Timer( search->elem_geom_time );
  #endif
}

void ContactTopology::Compute_Face_Geometry( VariableHandle POSITION, bool check_context, bool loft_shells )
{
  PRECONDITION( topology_type == PRIMARY );

  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(primary_face_list->EntityList());

  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  search->Timer()->Start_Timer( search->face_normal_geom_time );
  #endif
  
  //
  // Zero the node normals in preperation for assembly
  //
  for(int i=0; i<number_of_node_blocks; ++i) {
    if (!node_blocks[i]->Has_Normal_Attributes()) {
      int nnodes = primary_node_list->BlockNumEntities(i);
      ContactNode<Real>** nodes  = reinterpret_cast<ContactNode<Real>**>(primary_node_list->BlockEntityList(i));
      if(check_context) {
        for (int j=0; j<nnodes; ++j) {
          ContactNode<Real>* node = nodes[j];
          if (node->CheckContext(ContactTopologyEntity<Real>::GEOMETRY_UPDATE)) {
            Real* node_normal = node->Variable(NODE_NORMAL);
            for(int k=0 ; k<dimensionality ; ++k ) {
              node_normal[k] = 0.0;
            }
          }
        }
      } else {
        for (int j=0; j<nnodes; ++j) {
          ContactNode<Real>* node = nodes[j];
          Real* node_normal = node->Variable(NODE_NORMAL);
          for(int k=0 ; k<dimensionality ; ++k ) {
            node_normal[k] = 0.0;
          }
        }
      }
    }
  }

  if(check_context) {
    //
    // Have each face compute its normal and centroid
    //
    bool have_tied = search->Search_Data()->Have_Tied_Interactions();
    for (int i=0; i<number_of_primary_faces; ++i) {
      ContactFace<Real>* face = Faces[i];
      if (face->CheckContext(ContactTopologyEntity<Real>::GEOMETRY_UPDATE)) {
        face->Compute_Centroid(POSITION, CENTROID);
        face->Compute_Normal(POSITION, FACE_NORMAL);
        Real* face_normal = face->Variable(FACE_NORMAL);
        for(int j=0 ; j<face->Nodes_Per_Face(); ++j ){
          Real* node_normal = face->Node(j)->Variable(NODE_NORMAL);
          for(int k=0 ; k<dimensionality ; ++k ) {
            node_normal[k] += face_normal[k];
          }
        }
      } else if (have_tied) {
        face->Compute_Normal(POSITION, FACE_NORMAL);
      }
    }
  } else {
    //
    // Have each face compute its normal and centroid
    //
    for (int i=0; i<number_of_primary_faces; ++i) {
      ContactFace<Real>* face = Faces[i];
      face->Compute_Centroid(POSITION, CENTROID);
      face->Compute_Normal(POSITION, FACE_NORMAL);
      Real* face_normal = face->Variable(FACE_NORMAL);
      for(int j=0 ; j<face->Nodes_Per_Face(); ++j ){
        Real* node_normal = face->Node(j)->Variable(NODE_NORMAL);
        for(int k=0 ; k<dimensionality ; ++k ) {
          node_normal[k] += face_normal[k];
        }
      }
    }
  }  


  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  search->Timer()->Stop_Timer( search->face_normal_geom_time );
  #endif

  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    search->Timer()->Start_Timer( search->node_normal_geom_time );
  #endif

  #ifndef CONTACT_NO_MPI
  // SwapAdd the node normal (i.e., the parallel assembly)
  contact_swapadd_reg_data( SearchComm, *Node_Sym_Comm(), 
                            *comm_buffer, NODE_NORMAL, 3 );
  #endif

  // Normalize the node normal
  for(int i=0; i<number_of_node_blocks; ++i) {
    if (!node_blocks[i]->Has_Normal_Attributes()) {
      int nnodes = primary_node_list->BlockNumEntities(i);
      ContactNode<Real>** nodes  =  reinterpret_cast<ContactNode<Real>**>(primary_node_list->BlockEntityList(i));
      if(check_context) {
        for (int j=0; j<nnodes; ++j) {
          ContactNode<Real>* node = nodes[j];
          if (node->CheckContext(ContactTopologyEntity<Real>::GEOMETRY_UPDATE)) {
            Real* node_normal = node->Variable(NODE_NORMAL);
            Normalize(node_normal);
          }
        }
      } else {
        for (int j=0; j<nnodes; ++j) {
          Real* node_normal = nodes[j]->Variable(NODE_NORMAL);
          Normalize(node_normal);
	}          
      }
    }
  }
  
  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    search->Timer()->Stop_Timer( search->node_normal_geom_time );
  #endif

  // PJSA: Loft the shell nodes to the appropriate location
  if( loft_shells && Have_Shells() ){
    int num_configs = 3;
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      search->Timer()->Start_Timer( search->shell_loft_geom_time );
    #endif
    Shell_Handler()->Loft_Nodes( num_configs,
                                 CURRENT_POSITION,
                                 PREDICTED_POSITION,
                                 AUGMENTED_POSITION,
                                 FACE_NORMAL,
                                 LOFTING_VECTOR );
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      search->Timer()->Stop_Timer( search->shell_loft_geom_time );
    #endif
  }
}

void ContactTopology::Compute_Surface_Geometry( VariableHandle POSITION,
						int  num_configs,
                                                bool check_context )
{
#ifdef CONTACT_DEBUG
  PRECONDITION( topology_type == PRIMARY );
#endif
#if CONTACT_DEBUG_PRINT_LEVEL>=3
  ContactParOStream& postream = search->ParOStream();
#endif

  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(primary_node_list->EntityList());
  ContactEdge<Real>** Edges = 
    reinterpret_cast<ContactEdge<Real>**>(edge_list->EntityList());
  ContactFace<Real>** Faces = 
    reinterpret_cast<ContactFace<Real>**>(primary_face_list->EntityList());

  if (number_of_analytic_surfaces>0) {
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      search->Timer()->Start_Timer( search->analytic_surf_geom_time );
    #endif
    Real box_tol = search->Search_Data()->Max_Search_Tolerance();
    bool doit = false;
    ContactBoundingBox topology_box;
    for (int i=0; i<number_of_analytic_surfaces; ++i) {
      if (AnalyticSurfaces[i]->Surface_Type()==ContactSearch::PLANE) {
        doit = true;
        break;
      }
    }
    if (doit) {
      Real box_inflation = search->BoxInflation();
      for (int i=0; i<number_of_primary_nodes; ++i) {
        ContactNode<Real>* node = Nodes[i];
        ContactBoundingBox object_box;
        node->ComputeBoundingBoxForSearch(num_configs,
                                          REMAINING_GAP,
                                          CURRENT_POSITION,
                                          POSITION,
                                          auto_tol,
                                          box_inflation,
                                          object_box);
        topology_box.add_box(object_box);
      }
      contact_global_boundingbox(topology_box, SearchComm);
      topology_box.add_tolerance(box_tol);
    }
    // Calculate the bounding box for each analytic surface
    #if CONTACT_DEBUG_PRINT_LEVEL>=3
      postream<<"          Topology bounding box\n";
      postream<<"            min("<<topology_box.get_x_min()<<", "
                                  <<topology_box.get_y_min()<<", "
                                  <<topology_box.get_z_min()<<")\n";
      postream<<"            max("<<topology_box.get_x_max()<<", "
                                  <<topology_box.get_y_max()<<", "
                                  <<topology_box.get_z_max()<<")\n";
    #endif
    for (int i=0; i<number_of_analytic_surfaces; ++i) {
      ContactAnalyticSurface* surface = AnalyticSurfaces[i];
      surface->ComputeBoundingBox(&topology_box);
      ContactBoundingBox* bbox = surface->BoundingBox();
      if (surface->Surface_Type()==ContactSearch::PLANE && box_tol==0.0) {
        Real new_tol = 0.0;
        new_tol      = std::min(new_tol, (Real)(bbox->get_x_max()-bbox->get_x_min()));
        new_tol      = std::min(new_tol, (Real)(bbox->get_y_max()-bbox->get_y_min()));
        new_tol      = std::min(new_tol, (Real)(bbox->get_z_max()-bbox->get_z_min()));
        new_tol     *= 0.001;
        bbox->add_tolerance(new_tol);
      } else {
        bbox->add_tolerance(box_tol);
      }
      #if CONTACT_DEBUG_PRINT_LEVEL>=3
        postream<<"          Surface "<<i<<" bounding box\n";
        postream<<"            min("<<bbox->get_x_min()<<", "
                                    <<bbox->get_y_min()<<", "
                                    <<bbox->get_z_min()<<")\n";
        postream<<"            max("<<bbox->get_x_max()<<", "
                                    <<bbox->get_y_max()<<", "
                                    <<bbox->get_z_max()<<")\n";
      #endif
    }
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      search->Timer()->Stop_Timer( search->analytic_surf_geom_time );
    #endif
  }

  Compute_Face_Geometry( POSITION, check_context, false );

  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    search->Timer()->Start_Timer( search->edge_curve_geom_time );
    search->Timer()->Start_Timer( search->edge_curve_geom_time_serial );
  #endif
  //
  // Compute edge curvature (only have them in 3D)
  // Proximity based searches can cull out some of the edges,
  // if no proximity, no need to
  //
  if(check_context) {
    for (int i=0; i<number_of_edges; ++i) {
      ContactEdge<Real> *edge = Edges[i];
      int num_face_conns = edge->Number_Face_Connections();
      PRECONDITION( num_face_conns==1 || num_face_conns==2 );
      if( num_face_conns == 2 ) {
	ContactFace<Real> *face0 = edge->Face(0); 
	ContactFace<Real> *face1 = edge->Face(1); 
	if(face0->CheckContext(ContactTopologyEntity<Real>::GEOMETRY_UPDATE) ||
	   face1->CheckContext(ContactTopologyEntity<Real>::GEOMETRY_UPDATE)) {
	  *edge->Variable(CURVATURE) = 
	    Compute_Curvature( face0->Variable(CENTROID),
			       face0->Variable(FACE_NORMAL),
			       face1->Variable(CENTROID),
			       face1->Variable(FACE_NORMAL));
	  face0->SetEdgeCurvature(CURVATURE, edge);
	  face1->SetEdgeCurvature(CURVATURE, edge);
	}
      } else {
	ContactFace<Real> *face0 = edge->Face(0); 
	if(face0->CheckContext(ContactTopologyEntity<Real>::GEOMETRY_UPDATE)) {
	  *edge->Variable(CURVATURE) = 0.0;
	  face0->SetEdgeCurvature(CURVATURE, edge);
	}
      }
    }
  } else {
    for (int i=0; i<number_of_edges; ++i) {
      ContactEdge<Real> *edge = Edges[i];
      int num_face_conns = edge->Number_Face_Connections();
      PRECONDITION( num_face_conns==1 || num_face_conns==2 );
      if( num_face_conns == 2 ) {
	ContactFace<Real> *face0 = edge->Face(0); 
	ContactFace<Real> *face1 = edge->Face(1); 
	*edge->Variable(CURVATURE) = 
	  Compute_Curvature( face0->Variable(CENTROID),
			     face0->Variable(FACE_NORMAL),
			     face1->Variable(CENTROID),
			     face1->Variable(FACE_NORMAL));
	face0->SetEdgeCurvature(CURVATURE, edge);
	face1->SetEdgeCurvature(CURVATURE, edge);
      } else {
	ContactFace<Real> *face0 = edge->Face(0); 
	*edge->Variable(CURVATURE) = 0.0;
	face0->SetEdgeCurvature(CURVATURE, edge);
      }
    }
  }

  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  search->Timer()->Stop_Timer( search->edge_curve_geom_time_serial );
  #endif
  
  #ifndef CONTACT_NO_MPI
  // Now compute the curvature for shared edges (only needed in parallel)
  if ( contact_number_of_processors(SearchComm)>1 ) {

    // setup for communication
    int mesg = 1005;
    int my_proc = contact_processor_number(SearchComm);
    int num_shared_edges = Edge_SymComm->Size();
    int num_procs = Edge_SymComm->Num_Comm_Partners();
    RequestHandle * recv_handles = new RequestHandle[num_procs];
    RequestHandle * send_handles = new RequestHandle[num_procs];
    char * recv_buff = new char[sizeof(Real)*num_shared_edges*6];
    char * send_buff = new char[sizeof(Real)*num_shared_edges*6];
    Real * recv_edge_data, * send_edge_data, *rb, *sb;
    rb = recv_edge_data = reinterpret_cast<Real*>(recv_buff);
    sb = send_edge_data = reinterpret_cast<Real*>(send_buff);
    
    // register recvs for edges to recv
    for( int i=0 ; i<num_procs ; ++i ){
      int len = Edge_SymComm->Num_to_Proc(i)*6;
      recv_handles[i] = 
        contact_nonblocking_receive( mesg, rb, len, 
                                     Edge_SymComm->Comm_Proc_ID(i),
                                     SearchComm );
      rb+=len;
    }
    
    // sync to make sure all recvs have been posted
    contact_global_sync(SearchComm);
    
    // loop over shared edges and gather up data to send
    Real * ptr_centroids, * ptr_normals;
    ContactTopologyEntity<Real> ** comm_edges = 
      Edge_SymComm->Entity_List(0);
    for(int i=0; i<num_shared_edges; ++i) {
      ContactEdge<Real> * edge = static_cast<ContactEdge<Real> *>(comm_edges[i]);
      ptr_centroids = edge->Face(0)->Variable(CENTROID);
      send_edge_data[i*6]   = ptr_centroids[0];
      send_edge_data[i*6+1] = ptr_centroids[1];
      send_edge_data[i*6+2] = ptr_centroids[2];
      ptr_normals = edge->Face(0)->Variable(FACE_NORMAL);
      send_edge_data[i*6+3] = ptr_normals[0];
      send_edge_data[i*6+4] = ptr_normals[1];
      send_edge_data[i*6+5] = ptr_normals[2];
    }
    
    // post sends
    sb = send_edge_data;
    int proc_index = 0;
    while (proc_index < num_procs && Edge_SymComm->Comm_Proc_ID(proc_index) < my_proc) {
      sb += Edge_SymComm->Num_to_Proc(proc_index)*6;
      proc_index++;
    }
    if (proc_index == num_procs) {
      sb = send_edge_data;
      proc_index = 0;
    }
    for (int i=proc_index, j=0; j<num_procs; ++j) {
      int len = Edge_SymComm->Num_to_Proc(i)*6;
      #ifndef CONTACT_USE_BLOCKING_SEND
      send_handles[i] = contact_nonblocking_send( mesg, sb, len, 
                                                  Edge_SymComm->Comm_Proc_ID(i),
                                                  SearchComm );
      #else
      contact_blocking_send( mesg, sb, len, 
                             Edge_SymComm->Comm_Proc_ID(i),
                             SearchComm );
      #endif
      sb += len;
      if (++i == num_procs) {
        i  = 0;
        sb = send_edge_data;
      }
    }
    // compute concavities for all shared edges
    int index = 0;
    for ( int i = 0; i < num_procs; ++i ) {
      //
      //  Wait until the next message has arrived
      //
      contact_wait_msg_done( recv_handles[i] );
      ContactTopologyEntity<Real> ** edges_for_comm = Edge_SymComm->Entity_List(i);
        for ( int j = 0; j < Edge_SymComm->Num_to_Proc(i); ++j ){
        ContactEdge<Real> * edge = static_cast<ContactEdge<Real> *>(edges_for_comm[j]);
        *edge->Variable(CURVATURE) = 
          Compute_Curvature(&(recv_edge_data[index]),
                            &(recv_edge_data[index+3]),
                            &(send_edge_data[index]),
                            &(send_edge_data[index+3]));
        edge->Face(0)->SetEdgeCurvature(CURVATURE, edge);
        index+=6;
      }
    }

    #ifndef CONTACT_USE_BLOCKING_SEND
    //
    //  Assure that all sends have completed before reusing the buffers
    //
    for( int i=0 ; i<num_procs ; ++i ){
      contact_wait_msg_done( send_handles[i] );
    }
    #endif

    delete[] send_buff;
    delete[] recv_buff;
    delete[] recv_handles;
    delete[] send_handles;

  }
  #endif

  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
  search->Timer()->Stop_Timer( search->edge_curve_geom_time );
  #endif
  
  if (normal_smoothing_status) {
  #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    search->Timer()->Start_Timer( search->edge_smooth_geom_time );
  #endif
  
    for (int i=0; i<number_of_edges; ++i) {
      Edges[i]->Compute_Smoothed_Normal(FACE_NORMAL);
    }
    
    #ifndef CONTACT_NO_MPI
    // Now compute the smoothed normal for shared edges (only needed in parallel)
    if ( contact_number_of_processors(SearchComm) > 1 ) {

      // setup for communication
      int mesg = 1005;
      int my_proc = contact_processor_number(SearchComm);
      int num_shared_edges = Edge_SymComm->Size();
      int num_procs = Edge_SymComm->Num_Comm_Partners();
      RequestHandle * recv_handles = new RequestHandle[num_procs];
      RequestHandle * send_handles = new RequestHandle[num_procs];
      char * recv_buff = new char[sizeof(Real)*num_shared_edges*3];
      char * send_buff = new char[sizeof(Real)*num_shared_edges*3];
      Real * recv_edge_data, * send_edge_data, *rb, *sb;
      rb = recv_edge_data = reinterpret_cast<Real*>(recv_buff);
      sb = send_edge_data = reinterpret_cast<Real*>(send_buff);
      
      // register recvs for edges to recv
      for( int i=0 ; i<num_procs ; ++i ){
        int len = Edge_SymComm->Num_to_Proc(i)*3;
        recv_handles[i] = 
          contact_nonblocking_receive( mesg, rb, len, 
                                       Edge_SymComm->Comm_Proc_ID(i),
                                       SearchComm );
        rb+=len;
      }
      
      // sync to make sure all recvs have been posted
      contact_global_sync(SearchComm);
      
      // loop over shared edges and gather up data to send
      Real* smooth_normal;
      ContactTopologyEntity<Real> ** comm_edges = 
        Edge_SymComm->Entity_List(0);
      for(int i=0; i<num_shared_edges; ++i) {
        ContactEdge<Real> * edge = static_cast<ContactEdge<Real> *>(comm_edges[i]);
        smooth_normal = edge->Variable(SMOOTHED_NORMAL);
        send_edge_data[i*3+0] = smooth_normal[0];
        send_edge_data[i*3+1] = smooth_normal[1];
        send_edge_data[i*3+2] = smooth_normal[2];
      }
      
      // post sends
      sb = send_edge_data;
      int proc_index = 0;
      while (proc_index < num_procs && Edge_SymComm->Comm_Proc_ID(proc_index) < my_proc) {
        sb += Edge_SymComm->Num_to_Proc(proc_index)*3;
        proc_index++;
      }
      if (proc_index == num_procs) {
        sb = send_edge_data;
        proc_index = 0;
      }
      int k = proc_index;
      for (int j=0; j<num_procs; ++j) {
        int len = Edge_SymComm->Num_to_Proc(k)*3;
        #ifndef CONTACT_USE_BLOCKING_SEND
        send_handles[k] = contact_nonblocking_send( mesg, sb, len, 
                                                    Edge_SymComm->Comm_Proc_ID(k),
                                                    SearchComm );
        #else
        contact_blocking_send( mesg, sb, len, 
                               Edge_SymComm->Comm_Proc_ID(k),
                               SearchComm );
        #endif
        sb += len;
        if (++k == num_procs) {
          k  = 0;
          sb = send_edge_data;
        }
      }
      // compute smoothed normals for all shared edges
      int index = 0;
      for ( int i = 0; i < num_procs; ++i ) {
        //
        //  Wait until the next message has arrived
        //
        contact_wait_msg_done( recv_handles[i] );
        ContactTopologyEntity<Real> ** edges_for_comm = Edge_SymComm->Entity_List(i);
        for ( int j = 0; j < Edge_SymComm->Num_to_Proc(i); ++j ){
          ContactEdge<Real> * edge = static_cast<ContactEdge<Real> *>(edges_for_comm[j]);
          smooth_normal     = edge->Variable(SMOOTHED_NORMAL);
          smooth_normal[0] += recv_edge_data[index  ];
          smooth_normal[1] += recv_edge_data[index+1];
          smooth_normal[2] += recv_edge_data[index+2];
          Real mag = std::sqrt( smooth_normal[0]*smooth_normal[0] +
		                smooth_normal[1]*smooth_normal[1] +
		                smooth_normal[2]*smooth_normal[2] );
          if( mag > 0.0 ){
            mag = 1.0/mag;
            smooth_normal[0] *= mag;
            smooth_normal[1] *= mag;
            smooth_normal[2] *= mag;
          }
          index+=3;
        }
      }

      #ifndef CONTACT_USE_BLOCKING_SEND
      //
      //  Assure that all sends have completed before reusing the buffers
      //
      for( int i=0 ; i<num_procs ; ++i ){
        contact_wait_msg_done( send_handles[i] );
      }
      #endif

      delete[] send_buff;
      delete[] recv_buff;
      delete[] recv_handles;
      delete[] send_handles;

    }
    #endif

    for (int i=0; i<number_of_primary_faces; ++i) {
      Faces[i]->SetEdgeSmoothedNormal(SMOOTHED_NORMAL);
    }
    
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
    search->Timer()->Stop_Timer( search->edge_smooth_geom_time );
    #endif

  }

  // Loft the shell nodes to the appropriate location
  if( Have_Shells() ){
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      search->Timer()->Start_Timer( search->shell_loft_geom_time );
    #endif
    Shell_Handler()->Loft_Nodes( num_configs, 
                                 CURRENT_POSITION,
                                 PREDICTED_POSITION,
                                 AUGMENTED_POSITION,
                                 FACE_NORMAL,
                                 LOFTING_VECTOR );
    #if !defined(CONTACT_NO_MPI) && defined(CONTACT_TIMINGS)
      search->Timer()->Stop_Timer( search->shell_loft_geom_time );
    #endif
  }
}

void ContactSearch::Build_Physical_Face_List_None(ContactSearch::Topology use_topology, ContactTopologyEntity<Real>::SearchContext status,
                                                  bool use_proximity)
{
  ContactTopology *topology(0);

  switch(use_topology){
    case ContactSearch::PRIMARY:
      topology = primary_topology;
      break;
    case ContactSearch::SECONDARY :
      topology = secondary_topology;
      break;
    case ContactSearch::SEARCH :
      topology = search_topology;
      break;
    default:
      PRECONDITION(0);
      break;
  }

  int number_of_nodes = topology->Number_of_Nodes();
  ContactNode<Real>** Nodes = reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());

  VariableHandle NODE_NORMAL = topology->Variable_Handle( ContactTopology::Node_Normal );
#ifdef CONTACT_DEBUG_NODE
  VariableHandle FACE_NORMAL = topology->Variable_Handle( ContactTopology::Face_Normal );
#endif

  //
  //  Zero the nodal physical face normals.  
  //
  for (int i=0; i<number_of_nodes; ++i) {
    ContactNode<Real>* node = Nodes[i];
    Real *normals[3];
    *(NUMBER_PHYSICAL_FACES.Get_Scratch(i)) = 1;
    normals[0]     = PHYSICAL_FACE_NORMAL_1.Get_Scratch(i);
    normals[1]     = PHYSICAL_FACE_NORMAL_2.Get_Scratch(i);
    normals[2]     = PHYSICAL_FACE_NORMAL_3.Get_Scratch(i);
    normals[0][0]  = normals[0][1] = normals[0][2] = 0.0;
    normals[1][0]  = normals[1][1] = normals[1][2] = 0.0;
    normals[2][0]  = normals[2][1] = normals[2][2] = 0.0;
    if (!node->CheckContext(status))                           continue;
    if (node->Ownership()     != ContactTopologyEntity<Real>::OWNED) continue;
    if (node->Physical_Type() == ContactNode<Real>::SHELL_TAB_NODE ) continue;
    if (node->Physical_Type() == ContactNode<Real>::MIXED_TAB_NODE ) continue;
    if (node->Number_Face_Connections()==0)                    continue;
    Real* node_normal = node->Variable(NODE_NORMAL);
    normals[0][0]     = node_normal[0]; 
    normals[0][1]     = node_normal[1]; 
    normals[0][2]     = node_normal[2];
#ifdef CONTACT_DEBUG_NODE
    bool PRINT_THIS_NODE = primary_topology->Is_a_Debug_Node( node );
    if( PRINT_THIS_NODE ){
      postream << "Build Physical Face List for Debug Node ("
               << node->Exodus_ID() << ") "<<node->Global_ID()<<"\n";
      postream << "  Nfaces = " << node->Number_Face_Connections() << "\n";
      for( int ii=0 ; ii<node->Number_Face_Connections() ; ++ii ){
        Real * fn = node->GetFace(ii)->Variable(FACE_NORMAL);
        postream << "  Connected to face " 
                 << node->GetFace(ii)->Global_ID()
                 << " with normal " << fn[0] << " " << fn[1] << " " 
                 << fn[2] << "\n";
      }
      int num_pf = (int)*(NUMBER_PHYSICAL_FACES.Get_Scratch(i));
      for( int j=0; j<num_pf; ++j ){
	if( PRINT_THIS_NODE ){
	  postream << "  Physical Face Normal " << j << " has normal "
		   << normals[j][0] << " " << normals[j][1] << " "
		   << normals[j][2] <<"\n";
	}
      }
    }
#endif
  }
#ifdef CONTACT_DEBUG_NODE
  postream.flush();
#endif
  error_code = (ContactErrorCode) contact_global_error_check( error_code, SearchComm );  
}

void ContactSearch::Build_Physical_Face_List_FaceWalk(ContactSearch::Topology use_topology, 
                                                      ContactTopologyEntity<Real>::SearchContext status,
                                                      bool use_proximity)
{
  ContactTopology *topology(0);

  switch(use_topology){
    case ContactSearch::PRIMARY:
      topology = primary_topology;
      break;
    case ContactSearch::SECONDARY :
      topology = secondary_topology;
      break;
    case ContactSearch::SEARCH :
      topology = search_topology;
      break;
    default:
      PRECONDITION(0);
      break;
  }

  int number_of_nodes = topology->Number_of_Nodes();
  ContactNode<Real>** Nodes = reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());
  VariableHandle FACE_NORMAL = topology->Variable_Handle( ContactTopology::Face_Normal );
  VariableHandle NODE_NORMAL = topology->Variable_Handle( ContactTopology::Node_Normal );
  //
  //  Hard code a maximum number of faces that can be connected to a single node.  
  //  Assume an absurdly large value of 256.  The arrays below are required to sort
  //  faces for the multiple interactions case.  There are four face lists.  
  //  Faces for pf1, pf2, pf3, and 'pf4', a set of faces belonging to none of the
  //  other physical faces.
  //
  for (int i=0; i<number_of_nodes; ++i) {
    ContactNode<Real>* node = Nodes[i];
    int num_faces = node->Number_Face_Connections();
    POSTCONDITION(num_faces <= 256);
    Real *normals[3];
    //
    //  Initialize the node's physical face data
    //
    *(NUMBER_PHYSICAL_FACES.Get_Scratch(i)) = 0.0;
    normals[0]     = PHYSICAL_FACE_NORMAL_1.Get_Scratch(i);
    normals[1]     = PHYSICAL_FACE_NORMAL_2.Get_Scratch(i);
    normals[2]     = PHYSICAL_FACE_NORMAL_3.Get_Scratch(i);
    normals[0][0]  = normals[0][1] = normals[0][2] = 0.0;
    normals[1][0]  = normals[1][1] = normals[1][2] = 0.0;
    normals[2][0]  = normals[2][1] = normals[2][2] = 0.0;
    //
    // If we don't own this node we won't process it so don't bother building
    // its physical face data.  If this node doesn't have the status_flag of
    // status  we won't process it.  We also don't need to build the physical
    // face information for a TAB node since it won't be processed either.
    //
    if (use_proximity && !node->in_proximity)                  continue;
    if (!node->CheckContext(status))                           continue;
    if (node->Ownership()     != ContactTopologyEntity<Real>::OWNED) continue;
    if (node->Physical_Type() == ContactNode<Real>::SHELL_TAB_NODE ) continue;
    if (node->Physical_Type() == ContactNode<Real>::MIXED_TAB_NODE ) continue;
    if (num_faces==0)                                          continue;
#ifdef CONTACT_DEBUG_NODE
    bool PRINT_THIS_NODE = primary_topology->Is_a_Debug_Node( node );
    if( PRINT_THIS_NODE ){
      postream << "Build Physical Face List for Debug Node ("
               << node->Exodus_ID() << ") "<<node->Global_ID()<<"\n";
      postream << "  Nfaces = " << num_faces << "\n";
      for( int ii=0 ; ii<node->Number_Face_Connections() ; ++ii ){
        Real * fn = node->GetFace(ii)->Variable(FACE_NORMAL);
        postream << "  Connected to face " 
                 << node->GetFace(ii)->Global_ID()
                 << " with normal " << fn[0] << " " << fn[1] << " " 
                 << fn[2] << "\n";
      }
    }
#endif
    //
    //  The node potentially has multiple phsical faces
    //    normals:    will store the base physical face normal to test against, and at the end the final physical face normal
    //    num_pf:     will store the current number of valid physical faces found, at the end total number of node physical 
    //                faces 
    //    pf_index:   will store the index of the last face added, and at the end the index of the first face of a given
    //                physical face
    //
    //  Optimization notes:  For nearly all problems the vast majority of nodes will have one 
    //                       physical face (interior face nodes), much smaller numbers with 2 (edge),
    //                       3 (corner) or 0/4+(ill defined geometry) physical faces.
    //
    //                       Thus, optimize the algorithm to give best performance for the 1 physical face
    //                       case, and potentially slightly sub-optimal performance for the multi-physical
    //                       face cases
    //
    vector<std::pair<ContactFace<Real>*,int> > &face_list = node->ConnectedFaces();
    //
    //  Extract information for the very first face.
    //  Set the first physical of the node to the first face encountered
    //
    ContactFace<Real> *face0 = face_list[0].first;
    face_list[0].second = 0;
    Real* face0_normal = face0->Variable(FACE_NORMAL);

    Real* face1_normal = NULL;
    Real* face2_normal = NULL;

    int num_pf = 1;
    //
    //  Loop over all remaining faces.  This algorithm does the following:
    //    Sort the faces in up to four groups
    //      [1]  Faces within sharp_smooth_angle of physical face 1
    //      [2]  Faces within sharp_smooth_angle of physical face 2
    //      [3]  Faces within sharp_smooth_angle of physical face 3
    //      [4]  All remaining faces
    //
    //  This first optimized loop assumes that the current node has only one physical face.
    //  If that condition is ever violated break out of the loop and jump down to more general
    //  slightly slower coding.
    //
    int j;
    for (j=1; j<num_faces; ++j) {
      Real *facej_normal = face_list[j].first->Variable(FACE_NORMAL);
      //
      //  Determine if the current face is part of the 1st physical faces (based on angle difference).  
      //  If it is add it to that face.  If not set it is a new physical face and need to switch down
      //  to the more general loop below.
      //  
      //  Note, really want to find out if the current face angle is less than the critical
      //  angle.  
      //
      //  N1.N2 = ||N1|| ||N2|| cos (theta)
      //  Known that ||N1|| == ||N2|| == 1.0;
      //  So theta = acos(N1.N2)
      //
      //  However acos is quite expensive operation!
      //  So, instead compare the cosine of the face angle (N1.N2) to the cosine of the critical
      //  angle (cos_sharp_smooth_angle).  For angles < 180 degrees as long as 
      //  cos_angle_between_faces > cos_sharp_smooth_angle, then it is known that
      //  the true angle between faces is less than the sharp_smooth_angle.  As long as sharp
      //  smooth angle < 180 degrees (which should be a given) this calculation is valid.
      //
      if( acme::ComputeCurvatureFromNormals(face0_normal, facej_normal) < sharp_smooth_curvature ) {
        face_list[j].second = 0;
      } else {
	//
	//  The current node has more than one physical face. Add the current face as the first face of physical face 2.
	//  Prepare other data structures for use in more general loop below.
	//
	normals[1][0] = facej_normal[0];
	normals[1][1] = facej_normal[1];
	normals[1][2] = facej_normal[2];
	for(int k=0; k < j; ++k) {
          Real *facek_normal = face_list[k].first->Variable(FACE_NORMAL);
  	  normals[0][0] += facek_normal[0];
	  normals[0][1] += facek_normal[1];
	  normals[0][2] += facek_normal[2];
	}
	num_pf = 2; 
        face_list[j].second = 1;
        face1_normal = facej_normal;
	break;
      }
    } 
    if(j == num_faces) {
      //
      //  This node has only the single physical face, do simplified list completion, physical face is the same
      //  as the already computed node normal.
      //
      Real* node_normal = node->Variable(NODE_NORMAL);
      normals[0][0] = node_normal[0];
      normals[0][1] = node_normal[1];
      normals[0][2] = node_normal[2];
    } else {
      //
      //  This node has multiple physical faces, need to use more general coding to completely define the remaining
      //  phsyical faces.  Loop over each remaining face.  Add the face to one of the existing physical faces, or
      //  if it matches none, add it as the first face of a new physical face.
      //
      for (int k=j+1; k<num_faces; ++k) {
	Real *facek_normal = face_list[k].first->Variable(FACE_NORMAL);
	//
	//  Determine if the current face is part of one of the existing physical faces.  If it is
	//  add it to that face.  If not set it is a new physical face, if three physical faces
	//  are already defined dump the face into the garbage list.  
	//
        if(acme::ComputeCurvatureFromNormals(facek_normal, face0_normal) < sharp_smooth_curvature) {
          normals[0][0] += facek_normal[0];
	  normals[0][1] += facek_normal[1];
	  normals[0][2] += facek_normal[2];
          face_list[k].second = 0;
	} else if(acme::ComputeCurvatureFromNormals(facek_normal, face1_normal) < sharp_smooth_curvature) {
          normals[1][0] += facek_normal[0];
	  normals[1][1] += facek_normal[1];
	  normals[1][2] += facek_normal[2];
          face_list[k].second = 1;
	} else {
          if(num_pf > 2) {
            if(acme::ComputeCurvatureFromNormals(facek_normal, face2_normal) < sharp_smooth_curvature) {
              normals[2][0] += facek_normal[0];
	      normals[2][1] += facek_normal[1];
	      normals[2][2] += facek_normal[2];
              face_list[k].second = 2;
	    } else {
              face_list[k].second = -1;
	    }
	  } else {
            face2_normal = facek_normal;
            normals[2][0] += facek_normal[0];
            normals[2][1] += facek_normal[1];
	    normals[2][2] += facek_normal[2];
            face_list[k].second = 2;
            num_pf = 3;
	  }
	}
      }
      //
      //  Face lists are all constructed, net physical face normals.
      //
      for( int ipf=0; ipf<num_pf; ++ipf ){
	acme::Normalize(normals[ipf]);
      }
    }
    *(NUMBER_PHYSICAL_FACES.Get_Scratch(i)) = num_pf;
#ifdef CONTACT_DEBUG_NODE
    if( PRINT_THIS_NODE ){
      int npf = (int)*(NUMBER_PHYSICAL_FACES.Get_Scratch(i));
      for( int jj=0; jj<npf; ++jj ){
        if( PRINT_THIS_NODE ){
	  postream << "  Physical Face Normal " << jj << " has normal "
		   << normals[jj][0] << " " << normals[jj][1] << " "
		   << normals[jj][2] <<"\n";
        }
      }
    }
#endif
  }

#ifdef CONTACT_DEBUG_NODE
  postream.flush();
#endif
  error_code = (ContactErrorCode) contact_global_error_check( error_code, SearchComm );  
}

void ContactSearch::Build_Physical_Face_List_EdgeWalk(ContactSearch::Topology use_topology, ContactTopologyEntity<Real>::SearchContext status, 
                                                      bool use_proximity)
{
  //
  // =====================================================
  // =                                                   =
  // =  P H Y S I C A L   F A C E    A L G O R I T H M   =
  // =                                                   =
  // =====================================================
  //
  // The purpose of this routine is to compute the physical faces for a
  // node.  This is an approximate way to get back to a face-face type
  // algorithm without paying the full expense of the face-face imprinting
  // and enforcement.  For each node, the number of physical faces is determined.
  // The number of physical faces determines how many constraints a given node
  // can have (a maximum of three). This version takes into account the tracked
  // search and does not process the physical faces for nodes that are currently
  // being tracked.
  //
  // The steps to compute the physical faces are
  //   0) Initialization of data
  //   1) Find the 3 "sharpest" edges (i.e., those with the largest
  //      angle_between_faces).
  //   2) Classify these 3 "sharpest" edges as non-smooth or smooth based on 
  //      the sharp/non-sharp smoothing angle and count the number of non-smooth
  //      edges.
  //   3) Test the number of non-smooth edges.  If the number is zero, there is
  //      only one physical face so fill the data structure.  If then number is
  //      not zero then keep processing.
  //   4) Reorder Node-Edge back pointers so the edges are visited in
  //      a clockwise order.
  //   5) "Group" the faces into physical faces (including reordering the
  //      Node-Face back pointers) based on the number of non-smooth 
  //      "sharpest" edges (Note: this is a max of 3).  The logic is
  //       # non-smooth           Result
  //           1                  2 or 3 physical faces (depends on the
  //                                               on the pathology)
  //           2                  2 physical faces (seperated by the two
  //                                                non-smooth edges)
  //           3                  3 physical faces (seperated by the three
  //                                                non-smooth edges)
  //

  ContactTopology *topology(0);

  switch(use_topology){
    case ContactSearch::PRIMARY:
      topology = primary_topology;
      break;
    case ContactSearch::SECONDARY :
      topology = secondary_topology;
      break;
    case ContactSearch::SEARCH :
      topology = search_topology;
      break;
    default:
      PRECONDITION(0);
      break;
  }


  int number_of_nodes = topology->Number_of_Nodes();
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(topology->NodeList()->EntityList());

  VariableHandle NODE_NORMAL = 
    topology->Variable_Handle( ContactTopology::Node_Normal );
  VariableHandle FACE_NORMAL = 
    topology->Variable_Handle( ContactTopology::Face_Normal );

  for (int i=0; i<number_of_nodes; i++) {
    ContactNode<Real>* node = Nodes[i];
    Real *normals[3];
    *(NUMBER_PHYSICAL_FACES.Get_Scratch(i)) = 0.0;
    normals[0]     = PHYSICAL_FACE_NORMAL_1.Get_Scratch(i);
    normals[1]     = PHYSICAL_FACE_NORMAL_2.Get_Scratch(i);
    normals[2]     = PHYSICAL_FACE_NORMAL_3.Get_Scratch(i);
    normals[0][0]  = normals[0][1] = normals[0][2] = 0.0;
    normals[1][0]  = normals[1][1] = normals[1][2] = 0.0;
    normals[2][0]  = normals[2][1] = normals[2][2] = 0.0;
    if (!node->CheckContext(status))                           continue;
    if (node->Ownership()     != ContactTopologyEntity<Real>::OWNED) continue;
    if (node->Physical_Type() == ContactNode<Real>::SHELL_TAB_NODE ) continue;
    if (node->Physical_Type() == ContactNode<Real>::MIXED_TAB_NODE ) continue;
    if (node->Number_Face_Connections()==0)                    continue;
#ifdef CONTACT_DEBUG_NODE
    bool PRINT_THIS_NODE = primary_topology->Is_a_Debug_Node( node );
#endif
    if (node->Number_Face_Connections()==1) {
      *(NUMBER_PHYSICAL_FACES.Get_Scratch(i)) = 1;
      Real* node_normal = node->Variable(NODE_NORMAL);
      normals[0][0]     = node_normal[0]; 
      normals[0][1]     = node_normal[1]; 
      normals[0][2]     = node_normal[2];
      node->ConnectedFaces()[0].second = 0;

#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE ){
	postream << "Build Physical Face List New for Debug Node ("
		 << node->Exodus_ID() << ") "<<node->Global_ID()<<"\n";
        postream << "  Nfaces = " << node->Number_Face_Connections() << "\n";
	for( int ii=0 ; ii<node->Number_Face_Connections() ; ii++ ){
	  Real * fn = node->GetFace(ii)->Variable(FACE_NORMAL);
	  postream << "  Connected to face " 
		   << node->GetFace(ii)->Global_ID()
		   << " with normal " << fn[0] << " " << fn[1] << " " 
		   << fn[2] << "\n";
	}
	postream << "  Number of Physical Faces = "<<*(NUMBER_PHYSICAL_FACES.Get_Scratch(i))<<"\n";
      }
#endif

    } else {

#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE ){
	postream << "Build Physical Face List for Debug Node ("
		 << node->Exodus_ID() << ") "<<node->Global_ID()<<"\n";
        postream << "  Nfaces = " << node->Number_Face_Connections() << "\n";
	for( int ii=0 ; ii<node->Number_Face_Connections() ; ii++ ){
	  Real * fn = node->GetFace(ii)->Variable(FACE_NORMAL);
	  postream << "  Connected to face " 
		   << node->GetFace(ii)->Global_ID()
		   << " with normal " << fn[0] << " " << fn[1] << " " 
		   << fn[2] << "\n";
	}
      }
#endif
      
      // ***********************
      // STEP 0: Initialization
      // ***********************
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE ){
        postream << "  Step 0 ------------------------------------------\n";
      }
#endif
      
      int num_physical_face_normals = 0;
      
      int num_faces = node->Number_Face_Connections();
      vector<std::pair<ContactFace<Real>*, int> > &face_list = node->ConnectedFaces();
      int max_num_edge_connections = 2*num_faces;
      typedef struct {
        int id;
        int nfaces;
        ContactFace<Real>* faces[2];
        int face_edge[2];
        ContactNode<Real>* node1;
      } lw_edge;
      
      lw_edge* edge_list = new lw_edge[max_num_edge_connections];
      
      // load the lw_edge data
      int edge_nums[2];
      ContactNode<Real>* edge_nodes[2];
      
      face_list[0].first->GetEdgeInfo(node, edge_nodes, edge_nums);
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE ){
        postream<<"  face_list[0] info:  "<<edge_nodes[0]->Global_ID()<<"  "<<edge_nodes[1]->Global_ID()<<"\n";
      }
#endif
      int num_edges             = 2;
      edge_list[0].id           = 0;
      edge_list[0].nfaces       = 1;
      edge_list[0].faces[0]     = face_list[0].first;
      edge_list[0].node1        = edge_nodes[0];
      edge_list[0].face_edge[0] = edge_nums[0];
      edge_list[1].id           = 1;
      edge_list[1].nfaces       = 1;
      edge_list[1].faces[0]     = face_list[0].first;
      edge_list[1].node1        = edge_nodes[1];
      edge_list[1].face_edge[0] = edge_nums[1];


      for (int ii=1; ii<num_faces; ++ii) {
        bool n1_found = false;
        bool n2_found = false;
        face_list[ii].first->GetEdgeInfo(node, edge_nodes, edge_nums);
#ifdef CONTACT_DEBUG_NODE
        if( PRINT_THIS_NODE ){
          postream<<"  face_list["<<ii<<"] info:  "<<edge_nodes[0]->Global_ID()<<"  "<<edge_nodes[1]->Global_ID()<<"\n";
        }
#endif
        for (int j=0; j<num_edges; ++j) {
          if (edge_list[j].node1 == edge_nodes[0]) {
            n1_found = true;
            edge_list[j].nfaces       = 2;
            edge_list[j].faces[1]     = face_list[ii].first;
            edge_list[j].face_edge[1] = edge_nums[0];
          } else
          if (edge_list[j].node1 == edge_nodes[1]) {
            n2_found = true;
            edge_list[j].nfaces       = 2;
            edge_list[j].faces[1]     = face_list[ii].first;
            edge_list[j].face_edge[1] = edge_nums[1];
          }
        }
#ifdef CONTACT_DEBUG_NODE
        if( PRINT_THIS_NODE ){
          postream<<"  n1_found = "<<n1_found<<"\n";
          postream<<"  n2_found = "<<n2_found<<"\n";
        }
#endif
        if (!n1_found) {
          edge_list[num_edges].id           = num_edges;
          edge_list[num_edges].nfaces       = 1;
          edge_list[num_edges].faces[0]     = face_list[ii].first;
          edge_list[num_edges].node1        = edge_nodes[0];
          edge_list[num_edges].face_edge[0] = edge_nums[0];
          num_edges++;
        }
        if (!n2_found) {
          edge_list[num_edges].id           = num_edges;
          edge_list[num_edges].nfaces       = 1;
          edge_list[num_edges].faces[0]     = face_list[ii].first;
          edge_list[num_edges].node1        = edge_nodes[1];
          edge_list[num_edges].face_edge[0] = edge_nums[1];
          num_edges++;
        }
      }
      // now need to sort edge_list by node1->Global_ID()
      // just use a "straight insertion" sort since N is small
      for (int ii=1; ii<num_edges; ++ii) {
        lw_edge tmp = edge_list[ii];
        int iii = ii-1;
        while (iii>=0 && edge_list[iii].node1->Global_ID()>tmp.node1->Global_ID()) {
          edge_list[iii+1] = edge_list[iii];
          iii--;
        }
        edge_list[iii+1] = tmp;
      }
      lw_edge* old_edge_list = new lw_edge[num_edges];
      for (int ii=0; ii<num_edges; ++ii) old_edge_list[ii] = edge_list[ii];
      
      // *******************************
      // STEP 1: Find "sharpest" edges.
      // *******************************
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE ){
        postream << "  Step 1 ------------------------------------------\n";
      }
#endif
      Real curvature_list[3] = {-1.0, -1.0, -1.0};
      lw_edge sharpest_edges[3];
      for( int ii=0 ; ii<num_edges ; ii++ ){
        Real curvature = fabs(edge_list[ii].faces[0]->GetEdgeCurvature(edge_list[ii].face_edge[0]));
        if( curvature > curvature_list[2] ){
	  curvature_list[2] = curvature;
	  sharpest_edges[2] = edge_list[ii];
	  if( curvature_list[2] > curvature_list[1] ){
	    // switch 1 & 2
	    curvature_list[2] = curvature_list[1];
	    curvature_list[1] = curvature;
	    sharpest_edges[2] = sharpest_edges[1];
	    sharpest_edges[1] = edge_list[ii];
	    if( curvature_list[1] > curvature_list[0] ){
	      // switch 0 & 1
	      curvature_list[1] = curvature_list[0];
	      curvature_list[0] = curvature;
	      sharpest_edges[1] = sharpest_edges[0];
	      sharpest_edges[0] = edge_list[ii];
	    }
	  }
	}
      }
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE ){
        for (int m=0; m<3; m++) {
          postream << "    sharpest_edges["<<m<<"] angle = "<<curvature_list[m]<<"\n";
        }
      }
#endif
      
      // ***************************************
      // STEP 2: Classify the 3 sharpest edges
      // ***************************************
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE ){
        postream << "  Step 2 ------------------------------------------\n";
      }
#endif
      int num_non_smooth_edges = 0;

      for( int ii=0 ; ii<3 ; ii++ ){
	if( curvature_list[ii] > sharp_smooth_curvature )
	  num_non_smooth_edges++;
      }
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE ){
	postream << "    sharp_smooth_curvature = "<<sharp_smooth_curvature<<"\n";
	postream << "    num_non_smooth_edges   = "<<num_non_smooth_edges<<"\n";
      }
#endif
      
      // ********************************************
      // STEP 3: Handle the case of 1 Physical Face
      // ********************************************
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE ){
        postream << "  Step 3 ------------------------------------------\n";
      }
#endif
      bool Need_Normalizing = true;
      if( num_non_smooth_edges == 0 ){
	if( num_edges == 0 ){
	  // A node not connected to any faces
	  num_physical_face_normals = 0;
	} else {
	  // Put all the faces into one physical face
	  num_physical_face_normals = 1;
	  Real* node_normal = node->Variable(NODE_NORMAL);
	  normals[0][0] = node_normal[0];
	  normals[0][1] = node_normal[1];
	  normals[0][2] = node_normal[2];
	  Need_Normalizing = false;
          for(int iface = 0; iface < num_faces; ++iface) {
            face_list[iface].second = 0;
	  }
	}
      } else {
	
	// ***********************
	// STEP 4: Reorder edges
	// ***********************
#ifdef CONTACT_DEBUG_NODE
        if( PRINT_THIS_NODE ){
	  postream << "  Step 4 ------------------------------------------\n";
        }
#endif
	
	int num_edges_with_one_face = 0;
	int edges_with_one_face[2];
	for( int ii=0 ; ii<num_edges ; ii++ ){
	  if( edge_list[ii].nfaces == 1 ){
	    PRECONDITION( num_edges_with_one_face < 2 );
	    edges_with_one_face[num_edges_with_one_face++] = ii;
#ifdef CONTACT_DEBUG_NODE
            if( PRINT_THIS_NODE ){
	      postream << "  edges_with_one_face["<<num_edges_with_one_face-1<<"] = "<<ii
                       <<"   node = "<<edge_list[ii].node1->Global_ID()<<"\n";
            }
#endif
	  }
	}
	
	bool Need_to_Reorder_Edges = true;
	lw_edge starting_edge;
	lw_edge ending_edge;
	ContactFace<Real>* starting_face=NULL;
#ifdef CONTACT_DEBUG_NODE
	if (PRINT_THIS_NODE) {
	  postream << "num_edges_with_one_face " 
	           << num_edges_with_one_face << "\n";
	}
#endif
	switch( num_edges_with_one_face ){
	case 0:{ 
#ifdef CONTACT_DEBUG_NODE
	  if (PRINT_THIS_NODE) {
	    postream << "case 0\n";
	  }
#endif
	  // In this case we have a "circular" mesh around the node
	  // There are two possibilities:
	  //    1) One quadratic edge at a mid-edge node 
	  //    2) Multiple edges at a vertex node
	  if( num_edges == 1 ){
	    Need_to_Reorder_Edges = false;
	  } else {
	    // Use the sharpest edge as the starting edge
	    starting_edge = sharpest_edges[0];
	    ending_edge = starting_edge;
#ifdef CONTACT_DEBUG_NODE
            if( PRINT_THIS_NODE ){
	      postream << "    Starting edge = "<<node->Exodus_ID()<<"  "<<starting_edge.node1->Exodus_ID()<<"\n";
	      postream << "      Face 0: "<<starting_edge.faces[0]->Global_ID()<<"\n";
              for (int mm=0; mm<starting_edge.faces[0]->Nodes_Per_Face(); ++mm) {
                postream<<"        starting_edge.faces[0]->Node("<<mm<<") = "<<starting_edge.faces[0]->Node(mm)->Global_ID()<<"\n";
              }
	      postream << "      Face 1: "<<starting_edge.faces[1]->Global_ID()<<"\n";
              for (int mm=0; mm<starting_edge.faces[1]->Nodes_Per_Face(); ++mm) {
                postream<<"        starting_edge.faces[1]->Node("<<mm<<") = "<<starting_edge.faces[1]->Node(mm)->Global_ID()<<"\n";
              }
            }
#endif
	    // Now determine the starting face.  Only the correct starting face
	    // will have the clockwise edge in the nodes edge list
	    bool Found = false;
	    starting_face = starting_edge.faces[0];
            starting_edge.faces[0]->Clockwise_EdgeNode(starting_edge.face_edge[0], edge_nodes);
            if (edge_nodes[0]->Global_ID() == node->Global_ID()) {
	      for( int ii=0 ; ii<num_edges ; ii++ ){
                if (edge_list[ii].node1->Global_ID()==edge_nodes[1]->Global_ID()) {
		  Found = true;
		  break;
	        }
	      }
	    }
	    if( !Found ) starting_face = starting_edge.faces[1];
	  }
#ifdef CONTACT_DEBUG_NODE
          if( PRINT_THIS_NODE ){
	    postream << "    Starting face = "<<starting_face->Global_ID()<<"\n";
	    postream << "    Starting edge = "<<node->Exodus_ID()<<"  "<<starting_edge.node1->Exodus_ID()<<"\n";
	    postream << "    Ending   edge = "<<node->Exodus_ID()<<"  "<<ending_edge.node1->Exodus_ID()<<"\n";
          }
#endif
	  break;
	}
	case 1:{
#ifdef CONTACT_DEBUG_NODE
	  if (PRINT_THIS_NODE) {
	    postream << "case 1\n";
	  }
#endif
	  // This should only be for a quadratic edge midedge node
	  PRECONDITION( num_edges == 1 );
	  Need_to_Reorder_Edges = false;
	  break;
	}
	case 2:{
#ifdef CONTACT_DEBUG_NODE
	  if (PRINT_THIS_NODE) {
	    postream << "case 2\n";
	  }
#endif
	  starting_edge = edge_list[edges_with_one_face[0]];
	  // Grab the edge clockwise from it and see if its in the list
	  // if it is this is the edge to start with, if not we need to use
	  // the other edge.
	  bool Found_Edge = false;
          PRECONDITION(starting_edge.nfaces==1);
          starting_edge.faces[0]->Clockwise_EdgeNode(starting_edge.face_edge[0], edge_nodes);
#ifdef CONTACT_DEBUG_NODE
          if( PRINT_THIS_NODE ){
            postream<<"edge_num   = "<<starting_edge.face_edge[0]<<"\n";
            postream<<"edge_node0 = "<<edge_nodes[0]->Global_ID()<<"\n";
            postream<<"edge_node1 = "<<edge_nodes[1]->Global_ID()<<"\n";
            for (int mm=0; mm<starting_edge.faces[0]->Nodes_Per_Face(); ++mm) {
              postream<<"starting_edge.faces[0]->Node("<<mm<<") = "<<starting_edge.faces[0]->Node(mm)->Global_ID()<<"\n";
            }
          }
#endif
          if (edge_nodes[0]->Global_ID() == node->Global_ID()) {
          for (int ii=0; ii<num_edges; ++ii) {
            if (edge_list[ii].node1->Global_ID()==edge_nodes[1]->Global_ID()) {
	      Found_Edge = true;
	      ending_edge = edge_list[edges_with_one_face[1]];
	      starting_face = edge_list[edges_with_one_face[0]].faces[0];
#ifdef CONTACT_DEBUG_NODE
              if( PRINT_THIS_NODE ){
	        postream << "    Clockwise edge = "<<node->Exodus_ID()<<"  "<<edge_list[ii].node1->Exodus_ID()<<"\n";
              }
#endif
              break;
            }
          }
          }
#ifdef CONTACT_DEBUG_NODE
          if( PRINT_THIS_NODE ){
	    postream << "    Starting  edge = "<<node->Exodus_ID()<<"  "<<starting_edge.node1->Exodus_ID()<<"\n";
            postream << "    Found_Edge = "<<Found_Edge<<"\n";
	  }
#endif
          if( !Found_Edge ){
	    ending_edge = starting_edge;
	    starting_edge = edge_list[edges_with_one_face[1]];
	    starting_face = edge_list[edges_with_one_face[1]].faces[0];
	  }
#ifdef CONTACT_DEBUG_NODE
          if( PRINT_THIS_NODE ){
	    postream << "    Starting face = "<<starting_face->Global_ID()<<"\n";
	    postream << "    Starting edge = "<<node->Exodus_ID()<<"  "<<starting_edge.node1->Exodus_ID()<<"\n";
	    postream << "    Ending   edge = "<<node->Exodus_ID()<<"  "<<ending_edge.node1->Exodus_ID()<<"\n";
          }
#endif
	  break;
	}
	default:
	  POSTCONDITION(0);
	}
	
	// Now reorder the edges
#ifdef CONTACT_DEBUG_NODE
        if( PRINT_THIS_NODE ){
	  postream << "Need_to_Reorder_Edges = " << Need_to_Reorder_Edges << "\n";
        }
#endif
	if( Need_to_Reorder_Edges ){
	  int ii = 0;
	  bool Finished = false;
	  edge_list[ii++] = starting_edge;
#ifdef CONTACT_DEBUG_NODE
          if( PRINT_THIS_NODE ){
	    postream << "  0) edge_list["<<ii-1<<"] = "<< node->Exodus_ID() << " "
                     << edge_list[ii-1].node1->Exodus_ID() << "\n";
          }
#endif
	  lw_edge next_edge = starting_edge;
	  ContactFace<Real>* next_face = starting_face;
	  while( !Finished ){
	    // FIXIT:  next_edge = next_face->Clockwise_Edge( next_edge );
            int nn = 0;
            if (next_edge.nfaces>1) {
              if (next_face->Global_ID()==next_edge.faces[1]->Global_ID()) nn=1;
            }
            next_face->Clockwise_EdgeNode(next_edge.face_edge[nn], edge_nodes);
#ifdef CONTACT_DEBUG_NODE
            if( PRINT_THIS_NODE ){
              postream<<"edge_num   = "<<next_edge.face_edge[nn]<<"\n";
              postream<<"edge_node0 = "<<edge_nodes[0]->Global_ID()<<"\n";
              postream<<"edge_node1 = "<<edge_nodes[1]->Global_ID()<<"\n";
              for (int mm=0; mm<next_face->Nodes_Per_Face(); ++mm) {
                postream<<"next_face->Node("<<mm<<") = "<<next_face->Node(mm)->Global_ID()<<"\n";
              }
            }
#endif
            PRECONDITION(edge_nodes[0]->Global_ID()==node->Global_ID());
            for (int m=0; m<num_edges; ++m) {
              if (old_edge_list[m].node1->Global_ID()==edge_nodes[1]->Global_ID()) {
                next_edge = old_edge_list[m];
                break;
              }
            }
#ifdef CONTACT_DEBUG_NODE
            if( PRINT_THIS_NODE ){
	      postream << "    next edge = "<<node->Exodus_ID()<<"  "<<next_edge.node1->Exodus_ID()<<"\n";
            }
#endif
	    if( next_edge.id == ending_edge.id ){
	      Finished = true;
	      if( next_edge.id != starting_edge.id ) {
		edge_list[ii++] = next_edge;
#ifdef CONTACT_DEBUG_NODE
                if( PRINT_THIS_NODE ){
	          postream << "  1) edge_list["<<ii-1<<"] = "<< node->Exodus_ID() << " "
                           << edge_list[ii-1].node1->Exodus_ID() << "\n";
                }
#endif
              }
	    } else {
	      edge_list[ii++] = next_edge;
#ifdef CONTACT_DEBUG_NODE
              if( PRINT_THIS_NODE ){
	        postream << "  2) edge_list["<<ii-1<<"] = "<< node->Exodus_ID() << " "
                         << edge_list[ii-1].node1->Exodus_ID() << "\n";
              }
#endif
	      if( next_edge.nfaces != 2 ){
		std::sprintf(message,
			"Problem with mesh connectivity in Build Physical Faces"
			);
		errors->Add_Error_Message(message);
		std::sprintf(message,
			"  Check connectivity at node %d", 
			node->Exodus_ID() );
		errors->Add_Error_Message(message);
		Finished = true;
		error_code = INTERNAL_ERROR;
	      } else {
		if( next_edge.faces[0] != next_face )
		  next_face = next_edge.faces[0];
		else
		  next_face = next_edge.faces[1];
	      }
	    }
	  }
	  if( ii != num_edges ){
	    std::sprintf(message,
		    "Problem with mesh connectivity in Build Physical Faces"
		    );
	    errors->Add_Error_Message(message);
	    std::sprintf(message,
		    "  Check connectivity at node %d", 
		    node->Exodus_ID() );
	    errors->Add_Error_Message(message);
	    Finished = true;
	    error_code = INTERNAL_ERROR;
	  }
	}
	
	if( error_code != NO_ERROR ) break;

	// **********
	// STEP 5:
	// **********
#ifdef CONTACT_DEBUG_NODE
        if( PRINT_THIS_NODE ){
	  postream << "  Step 5\n";
        }
#endif
	switch( num_edges_with_one_face ){
	case 0: {
	  if( num_edges == 1 ){
	    // This should only be the case for a midedge node on a quadratic
	    //  edge which can only have two face connections.
	    PRECONDITION( num_non_smooth_edges == 1);
	    PRECONDITION( num_faces == 2 );
	    PRECONDITION( starting_edge.id = sharpest_edges[0].id );
#ifdef CONTACT_DEBUG_NODE
            if( PRINT_THIS_NODE ){
	      postream << "    midedge node case\n";
            }
#endif
	    // We have two physical faces
	    num_physical_face_normals = 2;
	    Real* face_normal = node->GetFace(0)->Variable(FACE_NORMAL);
	    normals[0][0] = face_normal[0];
	    normals[0][1] = face_normal[1];
	    normals[0][2] = face_normal[2];
	    face_normal = node->GetFace(1)->Variable(FACE_NORMAL);
	    normals[1][0] = face_normal[0];
	    normals[1][1] = face_normal[1];
	    normals[1][2] = face_normal[2];
	    Need_Normalizing = false;

            face_list[0].second = 0;
            for(int i = 1; i < num_faces; ++i) {
              face_list[i].second = 1;
	    }
	  } else {
	    // This is the circular case where we have one, two or three
	    // non-smooth edges.  
	    PRECONDITION( starting_edge.id == ending_edge.id );
	    PRECONDITION( num_non_smooth_edges == 1 ||
			  num_non_smooth_edges == 2 || 
			  num_non_smooth_edges == 3 );
	    PRECONDITION( num_edges = num_faces );
#ifdef CONTACT_DEBUG_NODE
            if( PRINT_THIS_NODE ){
	      postream << "    circular case, num_non_smooth_edges"<<num_non_smooth_edges<<"\n";
            }
#endif
	    switch( num_non_smooth_edges ){
	    case 1:
	      if( curvature_list[1] == curvature_list[2] )
		num_physical_face_normals = 3;
	      else {
		num_physical_face_normals = 2;
		sharpest_edges[2].id = -1;
	      }
	      break;
	    case 2:
	      num_physical_face_normals = 2;
	      sharpest_edges[2].id = -1;
	      break;
	    case 3:
	      num_physical_face_normals = 3;
	      break;
	    }
	    int pf_id = 0;
	    int j = 0;
	    Real* fn;
	    lw_edge next_edge = starting_edge;
	    ContactFace<Real>* next_face = starting_face;
	    for( int ii=0 ; ii<num_edges-1 ; ii++ ){
	      fn = next_face->Variable(FACE_NORMAL);
              face_list[j++].second = pf_id;              
	      normals[pf_id][0] += fn[0];
	      normals[pf_id][1] += fn[1];
	      normals[pf_id][2] += fn[2];
#ifdef CONTACT_DEBUG_NODE
              if( PRINT_THIS_NODE ){
	        postream << "    normal for face "<<next_face->Global_ID()
                         << ":  "<<fn[0]<<"  "<<fn[1]<<"  "<<fn[2]<<"\n";
	        postream << "    normal for physical face "<<pf_id
                         << ":  "<<normals[pf_id][0]<<"  "<<normals[pf_id][1]<<"  "<<normals[pf_id][2]<<"\n";
              }
#endif
	      //FIXIT:  next_edge = next_face->Clockwise_Edge( next_edge );
              int nn = 0;
              if (next_edge.nfaces>1) {
                if (next_face==next_edge.faces[1]) nn=1;
              }
              next_face->Clockwise_EdgeNode(next_edge.face_edge[nn], edge_nodes);
              for (int m=0; m<num_edges; ++m) {
                if (old_edge_list[m].node1==edge_nodes[1]) {
                  next_edge = old_edge_list[m];
                  break;
                }
              }
	      if( next_edge.faces[0] != next_face )
		next_face = next_edge.faces[0];
	      else
		next_face = next_edge.faces[1];
	      if( next_edge.id == sharpest_edges[1].id || 
		  next_edge.id == sharpest_edges[2].id ){
		pf_id++;
		POSTCONDITION( pf_id < num_physical_face_normals );
	      }
	    }
	    fn = next_face->Variable(FACE_NORMAL);
            face_list[j++].second = pf_id;
	    normals[pf_id][0] += fn[0];
	    normals[pf_id][1] += fn[1];
	    normals[pf_id][2] += fn[2];
#ifdef CONTACT_DEBUG_NODE
            if( PRINT_THIS_NODE ){
	      postream << "    normal for face "<<next_face->Global_ID()
                       << ":  "<<fn[0]<<"  "<<fn[1]<<"  "<<fn[2]<<"\n";
	      postream << "    normal for physical face "<<pf_id
                       << ":  "<<normals[pf_id][0]<<"  "<<normals[pf_id][1]<<"  "<<normals[pf_id][2]<<"\n";
            }
#endif
	    POSTCONDITION( j == num_faces );
	  }
	  break;
	}
	case 2: {
	  // This is the non-circular case where we have one or two non-smooth 
	  // edges.  We will ignore the third sharpest angle even its 
	  // non-smooth because we only allow three physical faces. 
	  PRECONDITION( starting_edge.id != ending_edge.id );
#ifdef CONTACT_DEBUG_NODE
          if( PRINT_THIS_NODE ){
	    postream << "    non-circular case\n";
          }
#endif
	  int j = 0;
	  Real* fn;
	  int pf_id = 0;
	  PRECONDITION( edge_list[0].nfaces == 1 );
	  lw_edge next_edge = edge_list[0];
	  ContactFace<Real>* next_face = edge_list[0].faces[0];
	  PRECONDITION( num_edges == num_faces + 1 );
	  // For the case of only one sharp edge, set sharpest_edges[1] = NULL
	  // so that we only get two physical faces
	  if( num_non_smooth_edges == 1 ){
	    sharpest_edges[1].id = -1;
	    num_physical_face_normals = 2;
	  } else
	    num_physical_face_normals = 3;
	  for( int ii=0 ; ii<num_edges-2 ; ii++ ){
	    fn = next_face->Variable(FACE_NORMAL);
            face_list[j++].second = pf_id;
	    normals[pf_id][0] += fn[0];
	    normals[pf_id][1] += fn[1];
	    normals[pf_id][2] += fn[2];
#ifdef CONTACT_DEBUG_NODE
            if( PRINT_THIS_NODE ){
	      postream << "    normal for face "<<next_face->Global_ID()
                       << ":  "<<fn[0]<<"  "<<fn[1]<<"  "<<fn[2]<<"\n";
	      postream << "    normal for physical face "<<pf_id
                       << ":  "<<normals[pf_id][0]<<"  "<<normals[pf_id][1]<<"  "<<normals[pf_id][2]<<"\n";
            }
#endif
	    //FIXIT:  next_edge = next_face->Clockwise_Edge( next_edge );
            int nn = 0;
            if (next_edge.nfaces>1) {
              if (next_face==next_edge.faces[1]) nn=1;
            }
            next_face->Clockwise_EdgeNode(next_edge.face_edge[nn],edge_nodes);
            for (int m=0; m<num_edges; ++m) {
              if (old_edge_list[m].node1==edge_nodes[1]) {
                next_edge = old_edge_list[m];
                break;
              }
            }
	    if( next_edge.faces[0] != next_face )
	      next_face = next_edge.faces[0];
	    else
	      next_face = next_edge.faces[1];
	    if( next_edge.id == sharpest_edges[0].id || 
		next_edge.id == sharpest_edges[1].id ){
	      pf_id++;
	      POSTCONDITION( pf_id < num_physical_face_normals );
	    }
	  }
          face_list[j++].second = pf_id;
	  fn = next_face->Variable(FACE_NORMAL);
	  normals[pf_id][0] += fn[0];
	  normals[pf_id][1] += fn[1];
	  normals[pf_id][2] += fn[2];
#ifdef CONTACT_DEBUG_NODE
          if( PRINT_THIS_NODE ){
	    postream << "    normal for face "<<next_face->Global_ID()
                     << ":  "<<fn[0]<<"  "<<fn[1]<<"  "<<fn[2]<<"\n";
	    postream << "    normal for physical face "<<pf_id
                     << ":  "<<normals[pf_id][0]<<"  "<<normals[pf_id][1]<<"  "<<normals[pf_id][2]<<"\n";
          }
#endif
	  POSTCONDITION( j == num_faces );
	  break;
	}
	default:{
	  POSTCONDITION( 0 );
	  break;
	}
	}
      }
      delete [] edge_list;
      delete [] old_edge_list;

      // all faces have been processed into physical faces... last step is to
      // normalize the physical face normals at node i
      if( Need_Normalizing ){
	for( int j=0; j<num_physical_face_normals; j++ ){
          acme::Normalize(normals[j]);
	}
      } 
#ifdef CONTACT_DEBUG_NODE
      for( int j=0; j<num_physical_face_normals; j++ ){
	if( PRINT_THIS_NODE ){
	  postream << "  Physical Face Normal " << j << " has normal "
		   << normals[j][0] << " " 
                   << normals[j][1] << " "
		   << normals[j][2] << "\n";
	}
      }
#endif
      *(NUMBER_PHYSICAL_FACES.Get_Scratch(i)) = num_physical_face_normals;
#ifdef CONTACT_DEBUG_NODE
      if( PRINT_THIS_NODE ){
	postream << "  Number of Physical Faces = "<<*(NUMBER_PHYSICAL_FACES.Get_Scratch(i))<<"\n";
      }
#endif
    }
  }
#ifdef CONTACT_DEBUG_NODE
  postream.flush();
#endif
  error_code = (ContactErrorCode) contact_global_error_check( error_code,
							      SearchComm );
}

//
//  Convention, Positive curvature is concave, neagative is convex, 0 is flat (or unknown/don't care)
//  Curvature recorded as +/- abs(cos(angle) - 1).
//  Note: cos(0.0) = 1.0, thus curvature varies as:
//    -2.0 = Most Convex
//     0.0 = flat (or unknown/don't care)
//     2.0 = Most Concave
//
//  The curvature is only relavant is if the surface is significantly concave or not
//
Real ContactTopology::Compute_Curvature( const Real* Centroid0,
	 				 const Real* normal0,
			 	   	 const Real* Centroid1,
					 const Real* normal1) {
  //
  // Compute a vector from the centroid of face 0 to face 1
  // Just compare the projections of face normal vectors onto this vector,
  // so no need to normalize the vector
  //
  Real vec_c0_c1[3];
  vec_c0_c1[0] = Centroid1[0] - Centroid0[0];
  vec_c0_c1[1] = Centroid1[1] - Centroid0[1];
  vec_c0_c1[2] = Centroid1[2] - Centroid0[2];

  Real proj0 = acme::Dot(normal0,vec_c0_c1);
  Real proj1 = acme::Dot(normal1,vec_c0_c1);
  if( proj0 < proj1) {
    return -acme::ComputeCurvatureFromNormals(normal0, normal1); //Convex
  } else {
    return  acme::ComputeCurvatureFromNormals(normal0, normal1); //Concave
  }
}

