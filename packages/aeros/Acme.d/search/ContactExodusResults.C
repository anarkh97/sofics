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


#ifndef CONTACT_NO_EXODUS_OUTPUT

#include "ContactTopology.h"
#include "ContactSearch.h"
#include "CString.h"
#include "exodusII.h"
#include "ContactElement.h"
#include "ContactErrors.h"
#include "ContactFace.h"
#include "ContactNode.h"
#include "ContactNodeNodeInteraction.h"
#include "ContactNodeFaceInteraction.h"
#include "ContactNodeSurfaceInteraction.h"
#include "ContactFaceFaceInteraction.h"
#include "ContactFaceCoverageInteraction.h"
#include "ContactElementElementInteraction.h"
#include "ContactFaceBlock.h"
#include "ContactEdgeBlock.h"
#include "ContactSearchData.h"
#include "ContactSymComm.h"
#include "Contact_Communication.h"
#include "ContactEnforcement.h"

// from asg
#include "ContactFixedSizeAllocator.h"

#include <iostream>
#include <cstdio>
#include <cstring>

ContactSearch::ContactErrorCode
ContactTopology::Exodus_Output_Results(int Exodus_ID, Real Time,
				       ContactErrors* error_handler,
				       ContactEnforcement* enforcement,
				       int num_proc_nodes,
				       int number_of_elem_blocks, 
				       int number_of_elems,
				       int number_edges_owned,
				       int Need_Sphere_Element_Block,
				       int num_nni, int max_nni,
                                       int num_nfi, int max_nfi,
				       int num_nsi, int max_nsi,
				       int num_ffi, int max_ffi,
				       int max_ffi_verts,
				       int num_fci, int max_fci,
				       int max_fci_verts,
				       int num_eei, int max_eei,
                                 const ContactSearch::Search_Option_Status&
				       normal_smoothing,
                                 const ContactSearch::Search_Option_Status&
				       multiple_interaction,
                                 const Real& sharp_smooth_ang, 
                                 const Real& normal_smoothing_distance, 
                                 const ContactSearch::Smoothing_Resolution &
				       smoothing_resolution,
                                 const ContactSearch::Search_Option_Status &
				       compute_node_areas )
{
  int i,j,k,l,n,ierr,index;
  ContactNodeFaceInteraction* cnfi;
  ContactFaceFaceInteraction<Real>* cffi;
  ContactTopologyEntity<Real>* entity;
  ContactInteractionEntity<Real>* interaction;
  int ProcArrayIndex;
  ContactNode<Real>** Nodes = 
    reinterpret_cast<ContactNode<Real>**>(primary_node_list->EntityList());
  
  int num_procs = contact_number_of_processors( SearchComm );
  bool PARALLEL = false;
  if( num_procs > 1 ) PARALLEL = true;
  
  //============================================================================
  //  C R E A T E   T I M E S T E P   F O R   R E S U L T S   O U T P U T
  //============================================================================
  // Create the time_step
  ierr = ex_put_time( Exodus_ID, 1, &Time );
  if( ierr!=0 ){
    std::sprintf(message,"Exodus Error %d from ex_put_time",ierr);
    error_handler->Add_Error_Message(message);
    return ContactSearch::EXODUS_ERROR;
  }
  
  //============================================================================
  //  S E T U P   N U M B E R   O F   O U T P U T   V A R I A B L E S
  //============================================================================
  // set up the number of plot variables
  int num_glob_vars = 0;
  int num_node_vars = 0;
  int num_elem_vars = 0;
  num_glob_vars = 11;
  if( enforcement != NULL ) 
    num_glob_vars += enforcement->Number_of_Global_Plot_Variables();
  if( dimensionality == 2){
    num_node_vars = 8;
    num_elem_vars = 2;
  } else {
    if( number_of_face_blocks ) {
      num_elem_vars = 4;
    } else if (number_of_element_blocks) {
      num_elem_vars = 2;
    }
    num_node_vars = 11;
  }
  if (num_nni>0) {
    num_node_vars += 2*max_nni;
  }
  if (num_nfi>0) {
    num_node_vars += 5*max_nfi;
    num_node_vars += 4*dimensionality*max_nfi;
    if( compute_node_areas == ContactSearch::ACTIVE ) num_node_vars += max_nfi;
  }
  if (num_nsi>0) {
    num_node_vars += max_nsi*5;
    num_node_vars += 4*dimensionality*max_nsi;
    if( compute_node_areas == ContactSearch::ACTIVE ) num_node_vars++;
  }
  if (num_ffi>0) {
    num_elem_vars += max_ffi*(2+6*max_ffi_verts);
  }
  if (num_fci>0) {
    num_elem_vars += max_fci*(1+2*max_fci_verts);
  }
  if (num_eei>0) {
    num_elem_vars += 2*max_eei;
  }
#ifndef CONTACT_NO_MPI
  num_node_vars += 4;
  num_elem_vars += 3;
#endif

  int num_nodal_enf_variables;
  int num_element_enf_variables;
  if( enforcement ) {
    num_nodal_enf_variables = enforcement->Number_of_Nodal_Plot_Variables();
    num_element_enf_variables = enforcement->Number_of_Element_Plot_Variables();
  } else {
    num_nodal_enf_variables = 0;
    num_element_enf_variables = 0;
  }

  num_node_vars += num_nodal_enf_variables;
  num_elem_vars += num_element_enf_variables;

  ierr = ex_put_var_param( Exodus_ID, "g", num_glob_vars );
  if( ierr!=0 && num_glob_vars!=0){
    std::sprintf(message,"Exodus Error %d from ex_put_var_param",ierr);
    error_handler->Add_Error_Message(message);
    return ContactSearch::EXODUS_ERROR;
  }
  ierr = ex_put_var_param( Exodus_ID, "n", num_node_vars );
  if( ierr!=0 && num_node_vars!=0 ){
    std::sprintf(message,"Exodus Error %d from ex_put_var_param",ierr);
    error_handler->Add_Error_Message(message);
    return ContactSearch::EXODUS_ERROR;
  }
  ierr = ex_put_var_param( Exodus_ID, "e", num_elem_vars );
  if( ierr!=0 && num_elem_vars!=0 ){
    std::sprintf(message,"Exodus Error %d from ex_put_var_param",ierr);
    error_handler->Add_Error_Message(message);
    return ContactSearch::EXODUS_ERROR;
  }
               
  //============================================================================
  //  C R E A T E   G L O B A L   V A R I A B L E   N A M E S
  //============================================================================
  i=0;
  char** glob_var_names = new char*[num_glob_vars];
  glob_var_names[i++]  = (char*) "num_nn_interactions";
  glob_var_names[i++]  = (char*) "num_nf_interactions";
  glob_var_names[i++]  = (char*) "num_ns_interactions";
  glob_var_names[i++]  = (char*) "num_ff_interactions";
  glob_var_names[i++]  = (char*) "num_fc_interactions";
  glob_var_names[i++]  = (char*) "num_ee_interactions";
  glob_var_names[i++]  = (char*) "mult_interaction_status";
  glob_var_names[i++]  = (char*) "norm_smoothing_status";
  glob_var_names[i++]  = (char*) "smoothing_angle";
  glob_var_names[i++]  = (char*) "smoothing_length";
  glob_var_names[i++]  = (char*) "smoothing_resolution";
  if( enforcement != NULL ){
    if( enforcement->Type() == ContactEnforcement::TDEnforcement ){
      glob_var_names[i++] = (char*) "dt_old";
      glob_var_names[i++] = (char*) "dt";
    }
  }
  POSTCONDITION(i==num_glob_vars);
  ierr = ex_put_var_names( Exodus_ID, "g", num_glob_vars, glob_var_names );
  if( ierr!=0 ){
    std::sprintf(message,"Exodus Error %d from ex_put_var_names(global)",ierr);
    error_handler->Add_Error_Message(message);
    return ContactSearch::EXODUS_ERROR;
  }
  delete [] glob_var_names;

  //============================================================================
  //  C R E A T E   N O D A L   V A R I A B L E   N A M E S
  //============================================================================
  //   Implementation Note: Define an integer variable after assigning the
  //                        name that can be used as the exodus index in
  //                        in the ex_put_nodal_var call.
  char** node_var_names = new char*[num_node_vars];
  index = 0;
  node_var_names[index++] = (char*) "displx";
  int DISPLX = index;
  node_var_names[index++] = (char*) "disply";
  int DISPLY = index;
  if( dimensionality == 3 ) node_var_names[index++]  = (char*) "displz";
  int DISPLZ = index;
  node_var_names[index++] = (char*) "nnormx";
  int NNORMX = index;
  node_var_names[index++] = (char*) "nnormy";
  int NNORMY = index;
  if( dimensionality == 3 ) node_var_names[index++] = (char*) "nnormz";
  int NNORMZ = index;
  node_var_names[index++] = (char*) "numcon";
  int NUMCON = index;
  node_var_names[index++] = (char*) "convecx";
  int CONVECX = index;
  node_var_names[index++] = (char*) "convecy";
  int CONVECY = index;
  if( dimensionality == 3 ) node_var_names[index++] = (char*) "convecz";
  int CONVECZ = index;
  node_var_names[index++] = (char*) "node_ek";
  int NODEEK = index;
  
#ifndef CONTACT_NO_MPI
  node_var_names[index++] = (char*) "Global_ID";
  int NODEGID = index;
  node_var_names[index++] = (char*) "Primary_Owner";
  int GIDPROC = index;
  node_var_names[index++] = (char*) "Primary_Local_ID";
  int GIDINDEX = index;
  node_var_names[index++] = (char*) "Secondary_Owner";
  int SECPROC = index;
#endif
  
  int  var_names_allocated_index = 0;
  int  num_var_names_allocated = 0;
  int* NNI_NODEID = NULL;
  int* NNI_DIST   = NULL;
  if (num_nni>0) {
    NNI_NODEID = new int[max_nni];
    NNI_DIST   = new int[max_nni];
    int  index_save = index;
    int num_to_allocate = 2*max_nni;
    for (i=0; i<num_to_allocate; ++i) {
      node_var_names[index_save+i] = new char[20];
    }
    for (i=0; i<max_nni; ++i) {
      std::sprintf(node_var_names[index],"node_id%d",i+1);
      NNI_NODEID[i] = ++index;
      std::sprintf(node_var_names[index],"distance%d",i+1);
      NNI_DIST[i] = ++index;
    }
    num_var_names_allocated += num_to_allocate;
    var_names_allocated_index = index_save;
  }
  
  int* NFI_FACEID  = NULL;
  int* NFI_ALG     = NULL;
  int* NFI_ND_EK   = NULL;
  int* NFI_GAP_CUR = NULL;
  int* NFI_GAP_OLD = NULL;
  int* NFI_IVECX   = NULL;
  int* NFI_IVECY   = NULL;
  int* NFI_IVECZ   = NULL;
  int* NFI_PBDIRX  = NULL;
  int* NFI_PBDIRY  = NULL;
  int* NFI_PBDIRZ  = NULL;
  int* NFI_NORMX   = NULL;
  int* NFI_NORMY   = NULL;
  int* NFI_NORMZ   = NULL;
  int* NFI_PFNORX  = NULL;
  int* NFI_PFNORY  = NULL;
  int* NFI_PFNORZ  = NULL;
  int* NFI_NDAREA  = NULL;
  if (num_nfi>0) {
    NFI_FACEID = new int[max_nfi];
    NFI_ALG    = new int[max_nfi];
    NFI_ND_EK  = new int[max_nfi];
    NFI_GAP_CUR= new int[max_nfi];
    NFI_GAP_OLD= new int[max_nfi];
    NFI_IVECX  = new int[max_nfi];
    NFI_IVECY  = new int[max_nfi];
    if( dimensionality == 3 )NFI_IVECZ  = new int[max_nfi];
    NFI_PBDIRX = new int[max_nfi];
    NFI_PBDIRY = new int[max_nfi];
    if( dimensionality == 3 )NFI_PBDIRZ = new int[max_nfi];
    NFI_NORMX  = new int[max_nfi];
    NFI_NORMY  = new int[max_nfi];
    if( dimensionality == 3 )NFI_NORMZ  = new int[max_nfi];
    NFI_PFNORX = new int[max_nfi];
    NFI_PFNORY = new int[max_nfi];
    if( dimensionality == 3 )NFI_PFNORZ = new int[max_nfi];
    if( compute_node_areas == ContactSearch::ACTIVE )
      NFI_NDAREA = new int[max_nfi];

    int index_save = index;
    int num_to_allocate = max_nfi*(5+4*dimensionality);
    if( compute_node_areas == ContactSearch::ACTIVE )
      num_to_allocate += max_nfi;
    for (i=0; i<num_to_allocate; ++i) {
      node_var_names[index_save+i] = new char[20];
    }
    for (i=0; i<max_nfi; ++i) {
      std::sprintf(node_var_names[index],"face_id%d",i+1);
      NFI_FACEID[i] = ++index;
      std::sprintf(node_var_names[index],"alg%d",i+1);
      NFI_ALG[i] = ++index;
      std::sprintf(node_var_names[index],"node_ek%d",i+1);
      NFI_ND_EK[i] = ++index;
      std::sprintf(node_var_names[index],"gapcur%d",i+1);
      NFI_GAP_CUR[i] = ++index;
      std::sprintf(node_var_names[index],"gapold%d",i+1);
      NFI_GAP_OLD[i] = ++index;
      std::sprintf(node_var_names[index],"ivec%dx",i+1);
      NFI_IVECX[i] = ++index;
      std::sprintf(node_var_names[index],"ivec%dy",i+1);
      NFI_IVECY[i] = ++index;
      if( dimensionality == 3 ) {
        std::sprintf(node_var_names[index],"ivec%dz",i+1);
        NFI_IVECZ[i] = ++index;
      }
      std::sprintf(node_var_names[index],"pbdir%dx",i+1);
      NFI_PBDIRX[i] = ++index;
      std::sprintf(node_var_names[index],"pbdir%dy",i+1);
      NFI_PBDIRY[i] = ++index;
      if( dimensionality == 3 ) {
        std::sprintf(node_var_names[index],"pbdir%dz",i+1);
        NFI_PBDIRZ[i] = ++index;
      }
      std::sprintf(node_var_names[index],"norm%dx",i+1);
      NFI_NORMX[i] = ++index;
      std::sprintf(node_var_names[index],"norm%dy",i+1);
      NFI_NORMY[i] = ++index;
      if( dimensionality == 3 ) {
        std::sprintf(node_var_names[index],"norm%dz",i+1);
        NFI_NORMZ[i] = ++index;
      }
      std::sprintf(node_var_names[index],"pfnorm%dx",i+1);
      NFI_PFNORX[i] = ++index;
      std::sprintf(node_var_names[index],"pfnorm%dy",i+1);
      NFI_PFNORY[i] = ++index;
      if( dimensionality == 3 ) {
        std::sprintf(node_var_names[index],"pfnorm%dz",i+1);
        NFI_PFNORZ[i] = ++index;
      }
      if( compute_node_areas == ContactSearch::ACTIVE ){
        std::sprintf(node_var_names[index],"ndarea%d",i+1);
        NFI_NDAREA[i] = ++index;
      }
    }
    num_var_names_allocated += num_to_allocate;
    if (var_names_allocated_index==0) var_names_allocated_index = index_save;
  }
  
  int* NSI_SURFACEID = NULL;
  int* NSI_NDEKA     = NULL;
  int* NSI_ALGA      = NULL;
  int* NSI_GAP_CURA  = NULL;
  int* NSI_GAP_OLDA  = NULL;
  int* NSI_IVECAX    = NULL;
  int* NSI_IVECAY    = NULL;
  int* NSI_IVECAZ    = NULL;
  int* NSI_PBDIRAX   = NULL;
  int* NSI_PBDIRAY   = NULL;
  int* NSI_PBDIRAZ   = NULL;
  int* NSI_NORMAX    = NULL;
  int* NSI_NORMAY    = NULL;
  int* NSI_NORMAZ    = NULL;
  int* NSI_PFNORMAX  = NULL;
  int* NSI_PFNORMAY  = NULL;
  int* NSI_PFNORMAZ  = NULL;
  int* NSI_NDAREAA   = NULL;
  if (num_nsi>0) {
    NSI_SURFACEID = new int[max_nsi];
    NSI_NDEKA     = new int[max_nsi];
    NSI_ALGA      = new int[max_nsi];
    NSI_GAP_CURA  = new int[max_nsi];
    NSI_GAP_OLDA  = new int[max_nsi];
    NSI_IVECAX    = new int[max_nsi];
    NSI_IVECAY    = new int[max_nsi];
    if( dimensionality == 3 )NSI_IVECAZ    = new int[max_nsi];
    NSI_PBDIRAX   = new int[max_nsi];
    NSI_PBDIRAY   = new int[max_nsi];
    if( dimensionality == 3 )NSI_PBDIRAZ   = new int[max_nsi];
    NSI_NORMAX    = new int[max_nsi];
    NSI_NORMAY    = new int[max_nsi];
    if( dimensionality == 3 )NSI_NORMAZ    = new int[max_nsi];
    NSI_PFNORMAX  = new int[max_nsi];
    NSI_PFNORMAY  = new int[max_nsi];
    if( dimensionality == 3 )NSI_PFNORMAZ  = new int[max_nsi];
    if( compute_node_areas == ContactSearch::ACTIVE )
      NSI_NDAREAA   = new int[max_nsi];
    
    int index_save = index;
    int num_to_allocate = max_nsi*(5+4*dimensionality);
    if( compute_node_areas == ContactSearch::ACTIVE )
      num_to_allocate += max_nsi;
    for (i=0; i<num_to_allocate; ++i) {
      node_var_names[index_save+i] = new char[20];
    }
    for (i=0; i<max_nsi; ++i) {
      std::sprintf(node_var_names[index], "surface_id%d", i+1);
      NSI_SURFACEID[i] = ++index;
      
      std::sprintf(node_var_names[index], "alga%d", i+1);
      NSI_ALGA[i] = ++index;
      
      std::sprintf(node_var_names[index], "ndeka%d", i+1);
      NSI_NDEKA[i] = ++index;
      
      std::sprintf(node_var_names[index], "gapcura%d", i+1);
      NSI_GAP_CURA[i] = ++index;
      
      std::sprintf(node_var_names[index], "gapolda%d", i+1);
      NSI_GAP_OLDA[i] = ++index;
      
      std::sprintf(node_var_names[index], "iveca%dx", i+1);
      NSI_IVECAX[i] = ++index;
      std::sprintf(node_var_names[index], "iveca%dy", i+1);
      NSI_IVECAY[i] = ++index;
      if( dimensionality == 3) {
        std::sprintf(node_var_names[index], "iveca%dz", i+1);
        NSI_IVECAZ[i] = ++index;
      }
      
      std::sprintf(node_var_names[index], "pbdira%dx", i+1);
      NSI_PBDIRAX[i] = ++index;
      std::sprintf(node_var_names[index], "pbdira%dy", i+1);
      NSI_PBDIRAY[i] = ++index;
      if( dimensionality == 3){
        std::sprintf(node_var_names[index], "pbdira%dz", i+1);
        NSI_PBDIRAZ[i] = ++index;
      }
      
      std::sprintf(node_var_names[index], "norma%dx", i+1);
      NSI_NORMAX[i] = ++index;
      std::sprintf(node_var_names[index], "norma%dy", i+1);
      NSI_NORMAY[i] = ++index;
      if( dimensionality ==3 ){
        std::sprintf(node_var_names[index], "norma%dz", i+1);
        NSI_NORMAZ[i] = ++index;
      }

      std::sprintf(node_var_names[index], "pfnorma%dx", i+1);
      NSI_PFNORMAX[i] = ++index;
      std::sprintf(node_var_names[index], "pfnorma%dy", i+1);
      NSI_PFNORMAY[i] = ++index;
      if( dimensionality == 3) {
        std::sprintf(node_var_names[index], "pfnorma%dz", i+1);
        NSI_PFNORMAZ[i] = ++index;
      }

      if( compute_node_areas == ContactSearch::ACTIVE ){
        std::sprintf(node_var_names[index], "ndareaa%d", i+1);
        NSI_NDAREAA[i] = ++index;
      }
    }
    num_var_names_allocated += num_to_allocate;
    if (var_names_allocated_index==0) var_names_allocated_index = index_save;
  }

  int index_save_node_enf = index;
  // if tdenforcement, use old naming convention else build node variable names
  // this is done to be consistent with old regression tests.
  if( enforcement ) {
    if( enforcement->Type() == ContactEnforcement::TDEnforcement ) {
	node_var_names[index++] = (char*) "Force_x";
	node_var_names[index++] = (char*) "Force_y";
	node_var_names[index++] = (char*) "Force_z";
	node_var_names[index++] = (char*) "Mass";
	node_var_names[index++] = (char*) "Density";
	node_var_names[index++] = (char*) "Wavespeed";
    } else {
      for (i=0; i<num_nodal_enf_variables; ++i) {
	node_var_names[index_save_node_enf+i] = new char[20];
      }
      for( i=0; i<num_nodal_enf_variables; ++i)
	std::sprintf(node_var_names[index++],"NodEnf%d",i+1);
    }
  }
  
  POSTCONDITION( index == num_node_vars );

  ierr = ex_put_var_names( Exodus_ID, "n", num_node_vars, node_var_names );
  if( ierr!=0 ){
    std::sprintf(message,"Exodus Error %d from ex_put_var_names(nodal)",ierr);
    error_handler->Add_Error_Message(message);
    return ContactSearch::EXODUS_ERROR;
  }
  for (i=0; i<num_var_names_allocated; ++i) {
    delete [] node_var_names[var_names_allocated_index+i];
  }

  if( enforcement && enforcement->Type() != ContactEnforcement::TDEnforcement ) {
    for (i=0; i<num_nodal_enf_variables; ++i) {
      delete [] node_var_names[index_save_node_enf+i];
    }
  }

  delete [] node_var_names;
  
  //============================================================================
  //  C R E A T E   E L E M E N T   V A R I A B L E   N A M E S
  //============================================================================
  char** elem_var_names = NULL;
  if (num_elem_vars) elem_var_names = new char*[num_elem_vars];
  int index_save_elem_enf = 0;
  index = 0;
  if( number_of_face_blocks ){
    elem_var_names[index++] = (char*) "fnormx";
    elem_var_names[index++] = (char*) "fnormy";
    if( dimensionality == 3 ){
      elem_var_names[index++] = (char*) "fnormz";
      if( number_of_edge_blocks ){
	elem_var_names[index++] = (char*) "curvature";
      }
    } else {
      if( number_of_edge_blocks ){
	elem_var_names[index++] = (char*) "curvature";
      }
    }
  } else if( number_of_element_blocks ){
    elem_var_names[index++] = (char*) "Volume";
    elem_var_names[index++] = (char*) "NumVolVolOverlaps";
    // add element enforcement variables
    index_save_elem_enf = index;
    for (i=0; i<num_element_enf_variables; ++i) {
      elem_var_names[index_save_elem_enf+i] = new char[20];
    }
    for( i=0; i<num_element_enf_variables; ++i)
      std::sprintf(elem_var_names[index++],"ElmEnf%d",i+1);
  }
  
  int index_save = index;
  int num_extra = max_ffi*(2+6*max_ffi_verts)+
                  max_fci*(1+2*max_fci_verts)+
                  max_eei*2;
  for (i=0; i<num_extra; ++i) {
    elem_var_names[index_save+i] = new char[20];
  }
  
  if (num_ffi>0) {
    for (i=0; i<max_ffi; ++i) {
      std::sprintf(elem_var_names[index++],"FFI%d_FACE_ID",i+1);
      std::sprintf(elem_var_names[index++],"FFI%d_NVERTS",i+1);
      for (j=0; j<max_ffi_verts; ++j) {
        std::sprintf(elem_var_names[index++],"FFI%d_SX%d",i+1,j+1);
        std::sprintf(elem_var_names[index++],"FFI%d_SY%d",i+1,j+1);
        std::sprintf(elem_var_names[index++],"FFI%d_MX%d",i+1,j+1);
        std::sprintf(elem_var_names[index++],"FFI%d_MY%d",i+1,j+1);
        std::sprintf(elem_var_names[index++],"FFI%d_EDGE%d",i+1,j+1);
        std::sprintf(elem_var_names[index++],"FFI%d_FLAG%d",i+1,j+1);
      }
    }
  }
  if (num_fci>0) {
    for (i=0; i<max_fci; ++i) {
      std::sprintf(elem_var_names[index++],"FCI%d_NVERTS",i+1);
      for (j=0; j<max_fci_verts; ++j) {
        std::sprintf(elem_var_names[index++],"FCI%d_X%d",i+1,j+1);
        std::sprintf(elem_var_names[index++],"FCI%d_Y%d",i+1,j+1);
      }
    }
  }
  if (num_eei>0) {
    for (i=0; i<max_eei; ++i) {
      std::sprintf(elem_var_names[index++],"EEI%d_ELEM_ID",i+1);
      std::sprintf(elem_var_names[index++],"EEI%d_VOLUME",i+1);
    }
  }
  
#ifndef CONTACT_NO_MPI
  elem_var_names[num_elem_vars-3] = (char*) "PrimaryOwner";
  elem_var_names[num_elem_vars-2] = (char*) "PrimaryLocalID";
  elem_var_names[num_elem_vars-1] = (char*) "SecondaryOwner";
#endif
  if( num_elem_vars ){
    ierr = ex_put_var_names( Exodus_ID, "e", num_elem_vars, elem_var_names );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_var_names(element)",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    for (i=0; i<num_element_enf_variables; ++i) {
      delete [] elem_var_names[index_save_elem_enf+i];
    }
    for (i=0; i<num_extra; ++i) {
      delete [] elem_var_names[index_save+i];
    }
    delete [] elem_var_names;
  
    //=========================================================================
    //  C R E A T E   E L E M E N T   T R U T H   T A B L E
    //=========================================================================

    // create the element truth table
    
    int* truth = new int[num_elem_vars*number_of_elem_blocks];
    for( i=0 ; i<num_elem_vars*number_of_elem_blocks ; ++i) truth[i] = 0;
    if( dimensionality == 2)
      for( i=0 ; i<num_elem_vars*number_of_elem_blocks ; ++i) truth[i] = 1;
    else{
      for( i=0 ; i<number_of_face_blocks ; ++i){
	j = 5;
	if (face_blocks[i]->Number_of_Faces()>0) {
	  truth[i*num_elem_vars+0] = 1;
	  truth[i*num_elem_vars+1] = 1;
	  truth[i*num_elem_vars+2] = 1;
	  truth[i*num_elem_vars+3] = 0;
	  truth[i*num_elem_vars+4] = 0;
	  for (k=0; k<max_ffi*(2+6*max_ffi_verts); ++k, ++j) {
	    truth[i*num_elem_vars+j] = 1;
	  }
	  for (k=0; k<max_fci*(1+2*max_fci_verts); ++k, ++j) {
	    truth[i*num_elem_vars+j] = 1;
	  }
	} else {
	  truth[i*num_elem_vars+0] = 0;
	  truth[i*num_elem_vars+1] = 0;
	  truth[i*num_elem_vars+2] = 0;
	  truth[i*num_elem_vars+3] = 0;
	  truth[i*num_elem_vars+4] = 0;
	  for (k=0; k<max_ffi*(2+6*max_ffi_verts); ++k, ++j) {
	    truth[i*num_elem_vars+j] = 0;
	  }
	  for (k=0; k<max_fci*(1+2*max_fci_verts); ++k, ++j) {
	    truth[i*num_elem_vars+j] = 0;
	  }
	}
#ifndef CONTACT_NO_MPI
	if (face_blocks[i]->Number_of_Faces()>0) {
	  truth[i*num_elem_vars+num_elem_vars-3] = 1;
	  truth[i*num_elem_vars+num_elem_vars-2] = 1;
	  truth[i*num_elem_vars+num_elem_vars-1] = 1;
	} else {
	  truth[i*num_elem_vars+num_elem_vars-3] = 0;
	  truth[i*num_elem_vars+num_elem_vars-2] = 0;
	  truth[i*num_elem_vars+num_elem_vars-1] = 0;
	}
#endif
      }
      // The next block(s) (i.e., the edges) only have curvature & angle_bf
      for( i=0 ; i<number_of_edge_blocks ; ++i){
	int num_edges_in_block = 0;
#ifndef CONTACT_NO_MPI
	if ( PARALLEL ) {
          ContactEdge<Real>** BlockEdges = 
            reinterpret_cast<ContactEdge<Real>**>(edge_list->BlockEntityList(i));
          for (j=0; j<edge_list->BlockNumEntities(i); ++j) {
	    if ( BlockEdges[j]->Ownership() == ContactTopologyEntity<Real>::OWNED ) 
	      num_edges_in_block++;
          }
	}
	else
#endif
	  num_edges_in_block = edge_blocks[i]->Number_of_Edges();
	j = i+number_of_face_blocks;
	truth[j*num_elem_vars + 0] = 0;
	truth[j*num_elem_vars + 1] = 0;
	truth[j*num_elem_vars + 2] = 0;
	if (num_edges_in_block>0) {
	  truth[j*num_elem_vars + 3] = 1;
	  truth[j*num_elem_vars + 4] = 1;
	} else {
	  truth[j*num_elem_vars + 3] = 0;
	  truth[j*num_elem_vars + 4] = 0;
	}
	for(k=0;k<max_ffi*(2+6*max_ffi_verts)+max_fci*(1+2*max_fci_verts); ++k){
	  truth[j*num_elem_vars+5+k] = 0;
	}
#ifndef CONTACT_NO_MPI
	if (num_edges_in_block>0) {
	  truth[j*num_elem_vars+num_elem_vars-3] = 1;
	  truth[j*num_elem_vars+num_elem_vars-2] = 1;
	  truth[j*num_elem_vars+num_elem_vars-1] = 1;
	} else {
	  truth[j*num_elem_vars+num_elem_vars-3] = 0;
	  truth[j*num_elem_vars+num_elem_vars-2] = 0;
	  truth[j*num_elem_vars+num_elem_vars-1] = 0;
	}
#endif
      }
      if( number_of_node_blocks > 1 ){
        int block_offset = 0;
	for( j=1 ; j<number_of_node_blocks ; ++j){
          if( node_blocks[j]->Type() == ContactSearch::POINT ){
	    int truth_offset = 
	      (number_of_face_blocks+number_of_edge_blocks+(j-1))*num_elem_vars;
	    for (k=0; k<num_elem_vars; ++k) {
	      truth[truth_offset+k] = 0;
	    }
            block_offset++;
          }
	}
      }
      if( number_of_element_blocks ){
	for( j=0 ; j<number_of_element_blocks ; ++j){
	  for( k=0 ; k<num_elem_vars ; ++k){
	    truth[j*num_elem_vars + k ] = 1;
	  }
	}
      }

      if( Need_Sphere_Element_Block ){
	int truth_offset = 
	  (number_of_face_blocks+number_of_edge_blocks+number_of_node_blocks-1)*
	  num_elem_vars;
	for (k=0; k<num_elem_vars; ++k) {
	  truth[truth_offset+k] = 0;
	}
      }
    }
    ierr = ex_put_elem_var_tab( Exodus_ID,number_of_elem_blocks,num_elem_vars,
				truth );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_elem_var_tab",ierr);
      error_handler->Add_Error_Message(message);
      //return ContactSearch::EXODUS_ERROR;
    }
    delete [] truth;
  }
  
  //============================================================================
  //  O U T P U T   G L O B A L   D A T A
  //============================================================================
  Real* gdata = new Real[num_glob_vars];
  for (i=0; i<num_glob_vars; ++i) gdata[i] = 0.0;
  i=0;
  gdata[i++] = num_nni;
  gdata[i++] = num_nfi;
  gdata[i++] = num_nsi;
  gdata[i++] = num_ffi;
  gdata[i++] = num_fci;
  gdata[i++] = num_eei;
  if (multiple_interaction == ContactSearch::ACTIVE) 
    gdata[i++] = 1.0;
  else
    i++;
  if (normal_smoothing == ContactSearch::ACTIVE) 
    gdata[i++] = 1.0;
  else 
    i++;
  if (multiple_interaction_status == ContactSearch::ACTIVE
      || normal_smoothing_status == ContactSearch::ACTIVE ) 
    gdata[i++] = sharp_smooth_ang;
  else
    i++;
  if ( normal_smoothing_status == ContactSearch::ACTIVE ) {
    gdata[i++] = normal_smoothing_distance;
    gdata[i++] = smoothing_resolution;
  } else {
    i+=2;
  }
  if( enforcement != NULL ){
    for( j=0 ; j<enforcement->Number_of_Global_Plot_Variables() ; ++j){
      gdata[i++] = enforcement->Get_Global_Plot_Variable( j );
    }
  }

  ierr = ex_put_glob_vars( Exodus_ID, 1, num_glob_vars, gdata );
  if( ierr!=0 ){
    std::sprintf(message,"Exodus Error %d from ex_put_glob_vars",ierr);
    error_handler->Add_Error_Message(message);
    return ContactSearch::EXODUS_ERROR;
  }
  delete [] gdata;

  //============================================================================
  //  O U T P U T   N O D E   D A T A
  //============================================================================
  Real* node_datax  = new Real[num_proc_nodes];
  Real* node_datay  = new Real[num_proc_nodes];
  Real* node_dataz  = new Real[num_proc_nodes];
  Real* node_datae1 = new Real[num_proc_nodes];
  Real* node_datae2 = new Real[num_proc_nodes];

  // Initialize the memory to zero so we don't have to worry about the
  // SPHERE node memory being uninitialized.
  std::memset( node_datax,  0, num_proc_nodes*sizeof(Real) );
  std::memset( node_datay,  0, num_proc_nodes*sizeof(Real) );
  std::memset( node_dataz,  0, num_proc_nodes*sizeof(Real) );
  std::memset( node_datae1, 0, num_proc_nodes*sizeof(Real) );
  std::memset( node_datae2, 0, num_proc_nodes*sizeof(Real) );

  // Set the displacements
  if( two_configurations ){
    for (n=0; n<number_of_primary_nodes; ++n) {
      Real* pos_1 = Nodes[n]->Variable(CURRENT_POSITION);
      Real* pos_2 = Nodes[n]->Variable(PREDICTED_POSITION);
      node_datax[n] = pos_2[0] - pos_1[0];
      node_datay[n] = pos_2[1] - pos_1[1];
      if( dimensionality == 3 ) node_dataz[n] = pos_2[2] - pos_1[2];
    }
  }

  // Output the displacements
  ierr = ex_put_nodal_var( Exodus_ID, 1, DISPLX, num_proc_nodes, node_datax );
  if( ierr!=0 ){
    std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
    error_handler->Add_Error_Message(message);
    return ContactSearch::EXODUS_ERROR;
  }
  ierr = ex_put_nodal_var( Exodus_ID, 1, DISPLY, num_proc_nodes, node_datay );
  if( ierr!=0 ){
    std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
    error_handler->Add_Error_Message(message);
    return ContactSearch::EXODUS_ERROR;
  }
  if( dimensionality == 3 ){
    ierr = ex_put_nodal_var( Exodus_ID, 1, DISPLZ, num_proc_nodes, 
                             node_dataz );  
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
  }

  // Collect the node normals
  for (n=0; n<number_of_primary_nodes; ++n) {
    Real* node_normal = Nodes[n]->Variable(NODE_NORMAL);
    node_datax[n] = node_normal[0];
    node_datay[n] = node_normal[1];
    if( dimensionality == 3) node_dataz[n] = node_normal[2];
  }

  // Output the node normals
  ierr = ex_put_nodal_var( Exodus_ID, 1, NNORMX, num_proc_nodes, node_datax );
  if( ierr!=0 ){
    std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
    error_handler->Add_Error_Message(message);
    return ContactSearch::EXODUS_ERROR;
  }
  ierr = ex_put_nodal_var( Exodus_ID, 1, NNORMY, num_proc_nodes, node_datay );
  if( ierr!=0 ){
    std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
    error_handler->Add_Error_Message(message);
    return ContactSearch::EXODUS_ERROR;
  }
  if( dimensionality == 3 ){
    ierr = ex_put_nodal_var( Exodus_ID, 1, NNORMZ, num_proc_nodes,
                             node_dataz );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
  }

  // Collect the number of constrained directions at each node
  for (n=0; n<number_of_primary_nodes; ++n) {
    Real num_constr = (Real) *Nodes[n]->Variable(NUM_KIN_CONSTR);
    node_datax[n] = num_constr;
  }
  
  // Output the number of constrained directions at each node
  ierr = ex_put_nodal_var( Exodus_ID, 1, NUMCON, num_proc_nodes, node_datax );
  if( ierr!=0 ){
    std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
    error_handler->Add_Error_Message(message);
    return ContactSearch::EXODUS_ERROR;
  }
  
  // Collect the node constraint vector
  for (n=0; n<number_of_primary_nodes; ++n) {
    Real* node_cons = Nodes[n]->Variable(KIN_CONSTR_VECTOR);
    node_datax[n] = node_cons[0];
    node_datay[n] = node_cons[1];
    if( dimensionality == 3) node_dataz[n] = node_cons[2];
  }

  // Output the node constraint vector
  ierr = ex_put_nodal_var( Exodus_ID, 1, CONVECX, num_proc_nodes, node_datax );
  if( ierr!=0 ){
    std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
    error_handler->Add_Error_Message(message);
    return ContactSearch::EXODUS_ERROR;
  }
  ierr = ex_put_nodal_var( Exodus_ID, 1, CONVECY, num_proc_nodes, node_datay );
  if( ierr!=0 ){
    std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
    error_handler->Add_Error_Message(message);
    return ContactSearch::EXODUS_ERROR;
  }
  if( dimensionality == 3 ){
    ierr = ex_put_nodal_var( Exodus_ID, 1, CONVECZ, num_proc_nodes,
                             node_dataz );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
  }
    
  // Output the node entity key
  for (n=0; n<number_of_primary_nodes; ++n) {
    node_datax[n] = Nodes[n]->Entity_Key();
  }
  ierr = ex_put_nodal_var( Exodus_ID, 1, NODEEK, num_proc_nodes, node_datax );
  
#ifndef CONTACT_NO_MPI  
  // Collect the global id at each node
  for (n=0; n<number_of_primary_nodes; ++n) {
    node_datax[n] = (Real) Nodes[n]->Global_ID().LoInt();
    node_datay[n] = (Real) Nodes[n]->Owner();
    //node_dataz[n] = (Real) Nodes[n]->ProcArrayIndex();
    node_dataz[n] = (Real) Nodes[n]->fcs_index;
  }
  
  // Output the global ids at each node
  ierr = ex_put_nodal_var( Exodus_ID, 1, NODEGID, 
                           num_proc_nodes, node_datax );
  if( ierr!=0 ){
    std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
    error_handler->Add_Error_Message(message);
    return ContactSearch::EXODUS_ERROR;
  } 
  ierr = ex_put_nodal_var( Exodus_ID, 1, GIDPROC, 
                           num_proc_nodes, node_datay );
  if( ierr!=0 ){
    std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
    error_handler->Add_Error_Message(message);
    return ContactSearch::EXODUS_ERROR;
  } 
  ierr = ex_put_nodal_var( Exodus_ID, 1, GIDINDEX, 
                           num_proc_nodes, node_dataz );
  if( ierr!=0 ){
    std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
    error_handler->Add_Error_Message(message);
    return ContactSearch::EXODUS_ERROR;
  } 
    
  // Collect the secondary owner at each node
  for (n=0; n<number_of_primary_nodes; ++n) {
    node_datax[n] = (Real) Nodes[n]->Secondary_Owner();
  }
  
  // Output the secondary at each node
  ierr = ex_put_nodal_var( Exodus_ID, 1, SECPROC, 
                           num_proc_nodes, node_datax );
  if( ierr!=0 ){
    std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
    error_handler->Add_Error_Message(message);
    return ContactSearch::EXODUS_ERROR;
  } 
#endif

  // output the node/node constraint data
  std::memset( node_datax, 0, num_proc_nodes*sizeof(Real) );
  std::memset( node_datay, 0, num_proc_nodes*sizeof(Real) );
  std::memset( node_dataz, 0, num_proc_nodes*sizeof(Real) );
  for (i=0; i<max_nni; ++i) {
    for (n=0; n<number_of_primary_nodes; ++n) {
      //ProcArrayIndex = Nodes[n]->ProcArrayIndex();
      ProcArrayIndex = Nodes[n]->fcs_index;
      ContactNodeNodeInteraction* cnni = 
        Nodes[n]->Get_NodeNode_Interaction( i );
      if( cnni ){
        if( PARALLEL ) {
          node_datax[ProcArrayIndex] = cnni->MasterNodeEntityData()->index_in_owner_proc_array+1;
        } else {
          //node_datax[ProcArrayIndex] = cnni->MasterNode()->ProcArrayIndex()+1;
          node_datax[ProcArrayIndex] = cnni->MasterNode()->fcs_index+1;
        }
        node_datay[ProcArrayIndex] = 
          cnni->Scalar_Var(ContactNodeNodeInteraction::DISTANCE);
      }
    }
#ifndef CONTACT_NO_MPI
    if (PARALLEL) {
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                      node_datay, 1 );
    }
#endif
    ierr = ex_put_nodal_var( Exodus_ID, 1, NNI_NODEID[i],
                             num_proc_nodes, node_datax );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    ierr = ex_put_nodal_var( Exodus_ID, 1, NNI_DIST[i],
                             num_proc_nodes, node_datay );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
  }
  if (NNI_NODEID ) delete [] NNI_NODEID ;
  if (NNI_DIST   ) delete [] NNI_DIST   ;

  // Output the node/face constraint data
  //
  //   If we have shells, the typical way of ensuring parallel consistency
  //   of the ordering of node-face interactions does not work because of
  //   the way that we construct the lofted topology. In this case, we need
  //   to reorder the constraints for each node based on comparing the 
  //   pushback direction.

  vector< vector<ContactNodeFaceInteraction*> >  Reordered_NodeFace_Interactions(number_of_primary_nodes);


  if (Have_Shells()) {
    bool * stored_cnfi = new bool[search->Max_Interactions_Per_Node()];
    for (i =0; i < search->Max_Interactions_Per_Node(); ++i) {
      stored_cnfi[i] = false;
    }
    for (n=0; n<number_of_primary_nodes; ++n) {
      Reordered_NodeFace_Interactions[n].resize(Nodes[n]->Number_NodeEntity_Interactions(), NULL);
      int num_interactions = Nodes[n]->Number_NodeEntity_Interactions();
      if ( num_interactions == 1 ) {
        cnfi = dynamic_cast<ContactNodeFaceInteraction*>(Nodes[n]->Get_NodeEntity_Interaction( 0 ));
        if (cnfi==NULL) continue;
        ContactNodeFaceInteraction* new_cnfi = 
                       ContactNodeFaceInteraction::new_ContactNodeFaceInteraction(
                       search->Get_Allocator(ContactSearch::ALLOC_ContactNodeFaceInteraction),
                       *cnfi);
        Reordered_NodeFace_Interactions[n][0] = new_cnfi;
        new_cnfi->Index(0);
      } else {
        for ( i = 0; i< num_interactions; ++i) {
          ContactHostGlobalID max_id;
          int max_cnfi_id = -1;
          ContactNodeFaceInteraction * max_cnfi = NULL;
          for ( l = 0; l< num_interactions; ++l) {
            if ( stored_cnfi[l] ) continue;
            cnfi = dynamic_cast<ContactNodeFaceInteraction*>(Nodes[n]->Get_NodeEntity_Interaction( l ));
            if (cnfi==NULL) continue;
            if ( max_cnfi_id < 0 ) {
              max_cnfi = cnfi;
              max_cnfi_id = l;
            } else {
              PRECONDITION (max_cnfi);
              Real* pb_dir_max = 
                max_cnfi->Vector_Var(ContactNodeFaceInteraction::PUSHBACK_DIR);
              Real* pb_dir_new = 
                cnfi->Vector_Var(ContactNodeFaceInteraction::PUSHBACK_DIR);
              if ( pb_dir_new[0] < pb_dir_max[0]) continue;
              if ( pb_dir_new[0] == pb_dir_max[0] ) {
                if ( pb_dir_new[1] < pb_dir_max[1]) continue;
                if ( pb_dir_new[1] == pb_dir_max[1] ){
                  if ( pb_dir_new[2] < pb_dir_max[2]) continue;
                }
              }
              max_cnfi = cnfi;
              max_cnfi_id = l;
            }
          }
          POSTCONDITION ( max_cnfi );

          ContactNodeFaceInteraction* new_cnfi = 
                       ContactNodeFaceInteraction::new_ContactNodeFaceInteraction(
                       search->Get_Allocator(ContactSearch::ALLOC_ContactNodeFaceInteraction),
                       *max_cnfi);

          if(Reordered_NodeFace_Interactions[n][i] != NULL) {
            Reordered_NodeFace_Interactions[n][i]->Delete_Frag(search->Get_Allocators());
	  }
          Reordered_NodeFace_Interactions[n][i] = new_cnfi;
          new_cnfi->Index(i);
          stored_cnfi[max_cnfi_id] = true;
        }
        for ( i = 0; i < num_interactions; ++i) {
          stored_cnfi[i] = false;
        }
      }
    } // end loop on nodes
    delete [] stored_cnfi;
  } // end if on shell existence

  // Now do the actual output of constraint data
  for (i=0; i<max_nfi; ++i) {
    for (n=0; n<number_of_primary_nodes; ++n) {
      //ProcArrayIndex = Nodes[n]->ProcArrayIndex();
      ProcArrayIndex = Nodes[n]->fcs_index;
      if (Have_Shells()) {
        if(i >= Reordered_NodeFace_Interactions[n].size()) {
          cnfi = NULL;
	} else {
          cnfi = dynamic_cast<ContactNodeFaceInteraction*>(Reordered_NodeFace_Interactions[n][i]);
	}
      } else {
        cnfi = dynamic_cast<ContactNodeFaceInteraction*>(Nodes[n]->Get_NodeEntity_Interaction( i ));
      }
      if( cnfi ){
        if( PARALLEL ) {
          node_datax[ProcArrayIndex] = cnfi->FaceEntityData()->index_in_owner_proc_array+1;
        } else {
          node_datax[ProcArrayIndex] = cnfi->Face()->fcs_index+1;
        }
        node_datay[ProcArrayIndex] = cnfi->Get_Source();
        node_dataz[ProcArrayIndex] = cnfi->Get_Gap_Cur();
        node_datae1[ProcArrayIndex] = cnfi->Get_Gap_Old();
        // Add 1 to get the fortran indexing
        node_datae2[ProcArrayIndex] = cnfi->Get_Node_Key()+1;
      } else {
        node_datax[ProcArrayIndex]  = 0;
        node_datay[ProcArrayIndex]  = 0;
        node_dataz[ProcArrayIndex]  = 0;
        node_datae1[ProcArrayIndex] = 0;
        node_datae2[ProcArrayIndex] = 0;
      }
    }
#ifndef CONTACT_NO_MPI
    if (PARALLEL) {
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                      node_datay, 1 );
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                      node_dataz, 1 );
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                      node_datae1, 1 );
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                      node_datae2, 1 );
    }
#endif
    ierr = ex_put_nodal_var( Exodus_ID, 1, NFI_FACEID[i],
                             num_proc_nodes, node_datax );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    ierr = ex_put_nodal_var( Exodus_ID, 1, NFI_ALG[i],
                             num_proc_nodes, node_datay );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }    
    ierr = ex_put_nodal_var( Exodus_ID, 1, NFI_GAP_CUR[i], 
                             num_proc_nodes, node_dataz );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    
    ierr = ex_put_nodal_var( Exodus_ID, 1, NFI_GAP_OLD[i], 
                             num_proc_nodes, node_datae1 );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }

    ierr = ex_put_nodal_var( Exodus_ID, 1, NFI_ND_EK[i],
			     num_proc_nodes, node_datae2 );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    
    // Collect and output the pushback direction
    Real pmag0,pmagc;
    for (n=0; n<number_of_primary_nodes; ++n) {
      //ProcArrayIndex = Nodes[n]->ProcArrayIndex();
      ProcArrayIndex = Nodes[n]->fcs_index;
      if (Have_Shells()) {
        if(i >= Reordered_NodeFace_Interactions[n].size()) {
          cnfi = NULL;
	} else {
          cnfi = dynamic_cast<ContactNodeFaceInteraction*>(Reordered_NodeFace_Interactions[n][i]);
	}
      } else {
        cnfi = dynamic_cast<ContactNodeFaceInteraction*>(Nodes[n]->Get_NodeEntity_Interaction( i ));
      }
      if( cnfi ){
        pmag0 = -cnfi->Scalar_Var(ContactNodeFaceInteraction::GAP_CUR);
        Real* pb_dir = 
          cnfi->Vector_Var(ContactNodeFaceInteraction::PUSHBACK_DIR);
        node_datax[ProcArrayIndex] = pb_dir[0];
        node_datay[ProcArrayIndex] = pb_dir[1];
        if( dimensionality == 3 ) node_dataz[ProcArrayIndex] = pb_dir[2];
      } else {
        node_datax[ProcArrayIndex] = 0;
        node_datay[ProcArrayIndex] = 0;
        if( dimensionality == 3 ) node_dataz[ProcArrayIndex] = 0;
      }
    }
#ifndef CONTACT_NO_MPI
    if (PARALLEL) {
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                      node_datax, 1 );
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                      node_datay, 1 );
      if( dimensionality == 3 ) 
        contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                        node_dataz, 1 );
    }
#endif
    ierr = ex_put_nodal_var( Exodus_ID, 1, NFI_PBDIRX[i],
                             num_proc_nodes, node_datax );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    ierr = ex_put_nodal_var( Exodus_ID, 1, NFI_PBDIRY[i],
                             num_proc_nodes, node_datay );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    if( dimensionality == 3 ) {
    ierr = ex_put_nodal_var( Exodus_ID, 1, NFI_PBDIRZ[i],
                             num_proc_nodes, node_dataz );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    }
    
    // Collect and output the interaction vector
    for (n=0; n<number_of_primary_nodes; ++n) {
      //ProcArrayIndex = Nodes[n]->ProcArrayIndex();
      ProcArrayIndex = Nodes[n]->fcs_index;
      if (Have_Shells()) {
        if(i >= Reordered_NodeFace_Interactions[n].size()) {
          cnfi = NULL;
	} else {
          cnfi = dynamic_cast<ContactNodeFaceInteraction*>(Reordered_NodeFace_Interactions[n][i]);
	}
      } else {
        cnfi = dynamic_cast<ContactNodeFaceInteraction*>(Nodes[n]->Get_NodeEntity_Interaction( i ));
      }
      if( cnfi ){
        pmagc = -cnfi->Scalar_Var(ContactNodeFaceInteraction::GAP_CUR);
        pmag0 = -cnfi->Scalar_Var(ContactNodeFaceInteraction::GAP_OLD);
        Real* pb_dir = 
          cnfi->Vector_Var(ContactNodeFaceInteraction::PUSHBACK_DIR);
        node_datax[ProcArrayIndex] = (pmagc+pmag0)*pb_dir[0];
        node_datay[ProcArrayIndex] = (pmagc+pmag0)*pb_dir[1];
        if( dimensionality == 3 ) node_dataz[ProcArrayIndex] = (pmagc+pmag0)*pb_dir[2];
      } else {
        node_datax[ProcArrayIndex] = 0;
        node_datay[ProcArrayIndex] = 0;
        if( dimensionality == 3 ) node_dataz[ProcArrayIndex] = 0;
      }
    }
#ifndef CONTACT_NO_MPI
    if (PARALLEL) {
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                      node_datax, 1 );
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                      node_datay, 1 );
      if( dimensionality == 3 ) 
        contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                        node_dataz, 1 );
    }
#endif
    ierr = ex_put_nodal_var( Exodus_ID, 1, NFI_IVECX[i],
                             num_proc_nodes, node_datax );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    ierr = ex_put_nodal_var( Exodus_ID, 1, NFI_IVECY[i],
                             num_proc_nodes, node_datay );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    if( dimensionality == 3 ) {
    ierr = ex_put_nodal_var( Exodus_ID, 1, NFI_IVECZ[i],
                             num_proc_nodes, node_dataz );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    }
    
    // Collect and output the normal direction
    for (n=0; n<number_of_primary_nodes; ++n) {
      //ProcArrayIndex = Nodes[n]->ProcArrayIndex();
      ProcArrayIndex = Nodes[n]->fcs_index;
      if (Have_Shells()) {
        if(i >= Reordered_NodeFace_Interactions[n].size()) {
          cnfi = NULL;
	} else {
          cnfi = dynamic_cast<ContactNodeFaceInteraction*>(Reordered_NodeFace_Interactions[n][i]);
	}
      } else {
        cnfi = dynamic_cast<ContactNodeFaceInteraction*>(Nodes[n]->Get_NodeEntity_Interaction( i ));
      }
      if( cnfi ){
        Real* normal_dir = 
          cnfi->Vector_Var(ContactNodeFaceInteraction::NORMAL_DIR);
        node_datax[ProcArrayIndex] = normal_dir[0];
        node_datay[ProcArrayIndex] = normal_dir[1];
        if( dimensionality == 3 ) node_dataz[ProcArrayIndex] = normal_dir[2];
      } else {
        node_datax[ProcArrayIndex] = 0;
        node_datay[ProcArrayIndex] = 0;
        if( dimensionality == 3 ) node_dataz[ProcArrayIndex] = 0;
      }
    }
#ifndef CONTACT_NO_MPI
    if (PARALLEL) {
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                      node_datax, 1 );
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                      node_datay, 1 );
      if( dimensionality == 3 ) 
        contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                        node_dataz, 1 );
    }
#endif
    ierr = ex_put_nodal_var( Exodus_ID, 1, NFI_NORMX[i],
                             num_proc_nodes, node_datax );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    ierr = ex_put_nodal_var( Exodus_ID, 1, NFI_NORMY[i],
                             num_proc_nodes, node_datay );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    if( dimensionality == 3 ) {
    ierr = ex_put_nodal_var( Exodus_ID, 1, NFI_NORMZ[i],
                             num_proc_nodes, node_dataz );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    }

    // Collect and output the physical face normal direction
    for (n=0; n<number_of_primary_nodes; ++n) {
      //ProcArrayIndex = Nodes[n]->ProcArrayIndex();
      ProcArrayIndex = Nodes[n]->fcs_index;
      if (Have_Shells()) {
        if(i >= Reordered_NodeFace_Interactions[n].size()) {
          cnfi = NULL;
	} else {
          cnfi = dynamic_cast<ContactNodeFaceInteraction*>(Reordered_NodeFace_Interactions[n][i]);
	}
      } else {
        cnfi = dynamic_cast<ContactNodeFaceInteraction*>(Nodes[n]->Get_NodeEntity_Interaction( i ));
      }
      if( cnfi ){
        Real* normal_dir = cnfi->Get_Physical_Face_Normal();
        node_datax[ProcArrayIndex] = normal_dir[0];
        node_datay[ProcArrayIndex] = normal_dir[1];
        if( dimensionality == 3 ) node_dataz[ProcArrayIndex] = normal_dir[2];
      } else {
        node_datax[ProcArrayIndex] = 0;
        node_datay[ProcArrayIndex] = 0;
        if( dimensionality == 3 ) node_dataz[ProcArrayIndex] = 0;
      }
    }
#ifndef CONTACT_NO_MPI
    if (PARALLEL) {
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                      node_datax, 1 );
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                      node_datay, 1 );
      if( dimensionality == 3 ) 
        contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                        node_dataz, 1 );
    }
#endif
    ierr = ex_put_nodal_var( Exodus_ID, 1, NFI_PFNORX[i],
                             num_proc_nodes, node_datax );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    ierr = ex_put_nodal_var( Exodus_ID, 1, NFI_PFNORY[i],
                             num_proc_nodes, node_datay );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    if( dimensionality == 3 ) {
    ierr = ex_put_nodal_var( Exodus_ID, 1, NFI_PFNORZ[i],
                             num_proc_nodes, node_dataz );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    }
    // Collect and output the node areas (if they were computed)
    if( compute_node_areas == ContactSearch::ACTIVE ){
      for (n=0; n<number_of_primary_nodes; ++n) {
        //ProcArrayIndex = Nodes[n]->ProcArrayIndex();
        ProcArrayIndex = Nodes[n]->fcs_index;
        if (Have_Shells()) {
          if(i >= Reordered_NodeFace_Interactions[n].size()) {
            cnfi = NULL;
	  } else {
            cnfi = dynamic_cast<ContactNodeFaceInteraction*>(Reordered_NodeFace_Interactions[n][i]);
	  }
        } else {
          cnfi = dynamic_cast<ContactNodeFaceInteraction*>(Nodes[n]->Get_NodeEntity_Interaction( i ));
        }
        if( cnfi ){
          Real node_area = cnfi->Get_Node_Area();
          node_datax[ProcArrayIndex] = node_area;
        } else {
          node_datax[ProcArrayIndex] = 0;
        }
      }
#ifndef CONTACT_NO_MPI
      if (PARALLEL) {
        contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                        node_datax, 1 );
      }
#endif
      ierr = ex_put_nodal_var( Exodus_ID, 1, NFI_NDAREA[i],
                               num_proc_nodes, node_datax );
      if( ierr!=0 ){
        std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
        error_handler->Add_Error_Message(message);
        return ContactSearch::EXODUS_ERROR;
      }
    }
  }

  if (NFI_FACEID ) delete [] NFI_FACEID ;
  if (NFI_ALG    ) delete [] NFI_ALG    ;
  if (NFI_ND_EK  ) delete [] NFI_ND_EK  ;
  if (NFI_GAP_CUR) delete [] NFI_GAP_CUR;
  if (NFI_GAP_OLD) delete [] NFI_GAP_OLD;
  if (NFI_IVECX  ) delete [] NFI_IVECX  ;
  if (NFI_IVECY  ) delete [] NFI_IVECY  ;
  if (NFI_IVECZ  ) delete [] NFI_IVECZ  ;
  if (NFI_PBDIRX ) delete [] NFI_PBDIRX ;
  if (NFI_PBDIRY ) delete [] NFI_PBDIRY ;
  if (NFI_PBDIRZ ) delete [] NFI_PBDIRZ ;
  if (NFI_NORMX  ) delete [] NFI_NORMX  ;
  if (NFI_NORMY  ) delete [] NFI_NORMY  ;
  if (NFI_NORMZ  ) delete [] NFI_NORMZ  ;
  if (NFI_PFNORX ) delete [] NFI_PFNORX ;
  if (NFI_PFNORY ) delete [] NFI_PFNORY ;
  if (NFI_PFNORZ ) delete [] NFI_PFNORZ ;
  if (NFI_NDAREA ) delete [] NFI_NDAREA ;
 
  if (Have_Shells()) {
    for (n=0; n<number_of_primary_nodes; ++n) {
      for(i = 0; i < Nodes[n]->Number_NodeEntity_Interactions(); ++i) {
        Reordered_NodeFace_Interactions[n][i]->Delete_Frag(search->Get_Allocators());
      }
    }
  }
 


  // Output the analtyic surface interactions
  for(i = 0; i < max_nsi; ++i){
    
    // Output the node entity key and node area (if computed)
    std::memset(node_datax,  0, sizeof(Real)*num_proc_nodes);
    std::memset(node_datay,  0, sizeof(Real)*num_proc_nodes);
    std::memset(node_dataz,  0, sizeof(Real)*num_proc_nodes);
    std::memset(node_datae1, 0, sizeof(Real)*num_proc_nodes);
    std::memset(node_datae2, 0, sizeof(Real)*num_proc_nodes);
//     for( i=0 ; i<num_proc_nodes /*number_of_primary_nodes*/ ; ++i){
//       node_datax[i] = 0;
//       node_datay[i] = 0;
//       node_dataz[i] = 0;
//     }
    for (n=0; n<number_of_primary_nodes; ++n) {
      ContactNodeSurfaceInteraction* cnsi = dynamic_cast<ContactNodeSurfaceInteraction*>(Nodes[n]->Get_NodeEntity_Interaction(i));
      //index = Nodes[n]->ProcArrayIndex();
      index = Nodes[n]->fcs_index;
      if (cnsi != NULL){
        node_datay[index]  = cnsi->SurfaceID();
        node_dataz[index]  = cnsi->Get_Source();
        node_datax[index]  = cnsi->Get_Node_Key() + 1;
        node_datae1[index] = cnsi->Get_Gap_Cur();
        node_datae2[index] = cnsi->Get_Gap_Old();
      }
    }
#ifndef CONTACT_NO_MPI
    if (PARALLEL) {
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                      node_datax, 1 );
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                      node_datay, 1 );
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                      node_dataz, 1 );
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                      node_datae1, 1 );
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                      node_datae2, 1 );
    }
#endif
    ierr = ex_put_nodal_var( Exodus_ID, 1, NSI_SURFACEID[i],
                             num_proc_nodes, node_datay );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    
    ierr = ex_put_nodal_var( Exodus_ID, 1, NSI_ALGA[i],
                             num_proc_nodes, node_dataz );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }

    ierr = ex_put_nodal_var( Exodus_ID, 1, NSI_NDEKA[i],
                             num_proc_nodes, node_datax );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    
    ierr = ex_put_nodal_var( Exodus_ID, 1, NSI_GAP_CURA[i],
                             num_proc_nodes, node_datae1 );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }

    ierr = ex_put_nodal_var( Exodus_ID, 1, NSI_GAP_OLDA[i],
                             num_proc_nodes, node_datae2 );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }

    // Output the node physical face normal
    std::memset(node_datax,  0, sizeof(Real)*num_proc_nodes);
    std::memset(node_datay,  0, sizeof(Real)*num_proc_nodes);
    std::memset(node_dataz,  0, sizeof(Real)*num_proc_nodes);
//     for( i=0 ; i<number_of_primary_nodes ; ++i){
//       node_datax[i] = 0;
//       node_datay[i] = 0;
//       node_dataz[i] = 0;
//     }
    
    for (n=0; n<number_of_primary_nodes; ++n) {
      ContactNodeSurfaceInteraction* cnsi = dynamic_cast<ContactNodeSurfaceInteraction*>(Nodes[n]->Get_NodeEntity_Interaction(i));
      //index = Nodes[n]->ProcArrayIndex();
      index = Nodes[n]->fcs_index;
      if (cnsi != NULL){
        Real* pf_normala = cnsi->Get_Physical_Face_Normal();
        node_datax[index] = pf_normala[0];
        node_datay[index] = pf_normala[1];
        node_dataz[index] = pf_normala[2];
      }
    }
#ifndef CONTACT_NO_MPI
    if (PARALLEL) {
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                      node_datax, 1 );
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                      node_datay, 1 );
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
                                      node_dataz, 1 );
    }
#endif
    ierr = ex_put_nodal_var( Exodus_ID, 1, NSI_PFNORMAX[i],
                             num_proc_nodes, node_datax );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    ierr = ex_put_nodal_var( Exodus_ID, 1, NSI_PFNORMAY[i],
                             num_proc_nodes, node_datay );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    ierr = ex_put_nodal_var( Exodus_ID, 1, NSI_PFNORMAZ[i],
                             num_proc_nodes, node_dataz );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    } 

    // Output the interaction vector (i.e., product of gap and normal)
    std::memset(node_datax,  0, sizeof(Real)*num_proc_nodes);
    std::memset(node_datay,  0, sizeof(Real)*num_proc_nodes);
    std::memset(node_dataz,  0, sizeof(Real)*num_proc_nodes);
//     for( i=0 ; i<number_of_primary_nodes ; ++i){
//       node_datax[i] = 0;
//       node_datay[i] = 0;
//       node_dataz[i] = 0;
//     }
    
    for (n=0; n<number_of_primary_nodes; ++n) {
      ContactNodeSurfaceInteraction* cnsi = dynamic_cast<ContactNodeSurfaceInteraction*>(Nodes[n]->Get_NodeEntity_Interaction(i));
      //index = Nodes[n]->ProcArrayIndex();
      index = Nodes[n]->fcs_index;
      if (cnsi != NULL){
        Real pmag = -cnsi->Get_Gap_Cur();
        pmag -= cnsi->Get_Gap_Old();
        Real* surf_normal = cnsi->Get_Normal();
        node_datax[index] = surf_normal[0]*pmag;
        node_datay[index] = surf_normal[1]*pmag;
        node_dataz[index] = surf_normal[2]*pmag;
      }
    }
#ifndef CONTACT_NO_MPI
    if (PARALLEL) {
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
				      node_datax, 1 );
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
				      node_datay, 1 );
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
				      node_dataz, 1 );
    }
#endif
    ierr = ex_put_nodal_var( Exodus_ID, 1, NSI_IVECAX[i],
                             num_proc_nodes, node_datax );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    ierr = ex_put_nodal_var( Exodus_ID, 1, NSI_IVECAY[i],
                             num_proc_nodes, node_datay );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    ierr = ex_put_nodal_var( Exodus_ID, 1, NSI_IVECAZ[i],
                             num_proc_nodes, node_dataz );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }

    // Output the push back vector 
    std::memset(node_datax,  0, sizeof(Real)*num_proc_nodes);
    std::memset(node_datay,  0, sizeof(Real)*num_proc_nodes);
    std::memset(node_dataz,  0, sizeof(Real)*num_proc_nodes);
    for(n = 0; n < number_of_primary_nodes; ++n){
      ContactNodeSurfaceInteraction* cnsi = dynamic_cast<ContactNodeSurfaceInteraction*>(Nodes[n]->Get_NodeEntity_Interaction(i));
      //index = Nodes[n]->ProcArrayIndex();
      index = Nodes[n]->fcs_index;
      if(cnsi != NULL){
        Real* pbdir = cnsi->Get_Pushback_Dir();
        node_datax[index] = pbdir[0];
        node_datay[index] = pbdir[1];
        node_dataz[index] = pbdir[2];
      }
    }
    
#ifndef CONTACT_NO_MPI
    if (PARALLEL) {
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
				      node_datax, 1 );
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
				      node_datay, 1 );
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
				      node_dataz, 1 );
    }
#endif
    ierr = ex_put_nodal_var( Exodus_ID, 1, NSI_PBDIRAX[i],
                             num_proc_nodes, node_datax );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    ierr = ex_put_nodal_var( Exodus_ID, 1, NSI_PBDIRAY[i],
                             num_proc_nodes, node_datay );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    ierr = ex_put_nodal_var( Exodus_ID, 1, NSI_PBDIRAZ[i],
                             num_proc_nodes, node_dataz );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }

    // Output the surface normal
    std::memset(node_datax,  0, sizeof(Real)*num_proc_nodes);
    std::memset(node_datay,  0, sizeof(Real)*num_proc_nodes);
    std::memset(node_dataz,  0, sizeof(Real)*num_proc_nodes);
    for(n = 0; n < number_of_primary_nodes; ++n){
      ContactNodeSurfaceInteraction* cnsi = dynamic_cast<ContactNodeSurfaceInteraction*>(Nodes[n]->Get_NodeEntity_Interaction(i));
      //index = Nodes[n]->ProcArrayIndex();
      index = Nodes[n]->fcs_index;
      if(cnsi != NULL){
        Real* surface_normal = cnsi->Get_Normal();
        node_datax[index] = surface_normal[0];
        node_datay[index] = surface_normal[1];
        node_dataz[index] = surface_normal[2];
      }
    }
    
#ifndef CONTACT_NO_MPI
    if (PARALLEL) {
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
				      node_datax, 1 );
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
				      node_datay, 1 );
      contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
				      node_dataz, 1 );
    }
#endif
    ierr = ex_put_nodal_var( Exodus_ID, 1, NSI_NORMAX[i],
                             num_proc_nodes, node_datax );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    ierr = ex_put_nodal_var( Exodus_ID, 1, NSI_NORMAY[i],
                             num_proc_nodes, node_datay );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }
    ierr = ex_put_nodal_var( Exodus_ID, 1, NSI_NORMAZ[i],
                             num_proc_nodes, node_dataz );
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    }

    if( compute_node_areas == ContactSearch::ACTIVE ){
      // Output the node area
      std::memset(node_datax,  0, sizeof(Real)*num_proc_nodes);
      for(n = 0; n < number_of_primary_nodes; ++n){
        ContactNodeSurfaceInteraction* cnsi = dynamic_cast<ContactNodeSurfaceInteraction*>(Nodes[n]->Get_NodeEntity_Interaction(i));
        //index = Nodes[n]->ProcArrayIndex();
        index = Nodes[n]->fcs_index;
        if(cnsi != NULL){
          node_datax[index] = cnsi->Get_Node_Area();
        }
      }

#ifndef CONTACT_NO_MPI
      if (PARALLEL) {
        contact_swapadd_data_array_fcs( SearchComm, *Node_SymComm, *comm_buffer,
				        node_datax, 1 );
      }
#endif
        ierr = ex_put_nodal_var( Exodus_ID, 1, NSI_NDAREAA[i],
			         num_proc_nodes, node_datax );
        if( ierr!=0 ){
	  std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
	  error_handler->Add_Error_Message(message);
	  return ContactSearch::EXODUS_ERROR;
        }
    }

  }
  
  //delete all the node surface arrays
  if(NSI_SURFACEID) delete [] NSI_SURFACEID;
  if(NSI_NDEKA)     delete [] NSI_NDEKA;
  if(NSI_ALGA)      delete [] NSI_ALGA;
  if(NSI_GAP_CURA)  delete [] NSI_GAP_CURA;
  if(NSI_GAP_OLDA)  delete [] NSI_GAP_OLDA;
  if(NSI_IVECAX)    delete [] NSI_IVECAX;
  if(NSI_IVECAY)    delete [] NSI_IVECAY;
  if(NSI_IVECAZ)    delete [] NSI_IVECAZ;
  if(NSI_PBDIRAX)   delete [] NSI_PBDIRAX;
  if(NSI_PBDIRAY)   delete [] NSI_PBDIRAY;
  if(NSI_PBDIRAZ)   delete [] NSI_PBDIRAZ;
  if(NSI_NORMAX)    delete [] NSI_NORMAX;
  if(NSI_NORMAY)    delete [] NSI_NORMAY;
  if(NSI_NORMAZ)    delete [] NSI_NORMAZ;
  if(NSI_PFNORMAX)  delete [] NSI_PFNORMAX;
  if(NSI_PFNORMAY)  delete [] NSI_PFNORMAY;
  if(NSI_PFNORMAZ)  delete [] NSI_PFNORMAZ;
  if(NSI_NDAREAA)   delete [] NSI_NDAREAA;

  // write node variables from enforcement model
  index = index_save_node_enf;
  for( i=0; i<num_nodal_enf_variables; ++i) {
    enforcement->Get_Nodal_Plot_Variable( i, node_datax );
    ierr = ex_put_nodal_var( Exodus_ID, 1, ++index, 
			     num_proc_nodes, node_datax );  
    if( ierr!=0 ){
      std::sprintf(message,"Exodus Error %d from ex_put_nodal_var",ierr);
      error_handler->Add_Error_Message(message);
      return ContactSearch::EXODUS_ERROR;
    } 
  }

  delete [] node_datax;
  delete [] node_datay;
  delete [] node_dataz;
  delete [] node_datae1;
  delete [] node_datae2;

  //============================================================================
  //  O U T P U T   E L E M E N T   D A T A
  //============================================================================
  if( number_of_face_blocks ){
    Real* elem_data1 = new Real[number_of_elems];
    Real* elem_data2 = new Real[number_of_elems];
    Real* elem_data3 = new Real[number_of_elems];
    Real* elem_data4 = new Real[number_of_elems];
    Real* elem_data5 = new Real[number_of_elems];
    Real* elem_data6 = new Real[number_of_elems];
    Real* elem_data7 = new Real[number_of_elems];
    Real* elem_data8 = new Real[number_of_elems];
    
    // Output the face normal
    for( i=0 ; i<number_of_face_blocks ; ++i){
      ContactFace<Real>** BlockFaces = 
        reinterpret_cast<ContactFace<Real>**>(primary_face_list->BlockEntityList(i));
      int num_faces = primary_face_list->BlockNumEntities(i);
      for (j=0; j<num_faces; ++j) {
	Real* fnorm = BlockFaces[j]->Variable(FACE_NORMAL);
	elem_data1[j] = fnorm[0];
	elem_data2[j] = fnorm[1];
	if( dimensionality == 3 ) elem_data3[j] = fnorm[2];
      }
      if (num_faces>0) {
	int ID = i+1;
	ierr = ex_put_elem_var( Exodus_ID, 1, 1, ID, num_faces, elem_data1 );
	if( ierr!=0 ){
	  std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
	  error_handler->Add_Error_Message(message);
                POSTCONDITION(ierr==0);
                delete [] elem_data1;
                delete [] elem_data2;
                delete [] elem_data3;
                delete [] elem_data4;
                delete [] elem_data5;
                delete [] elem_data6;
                delete [] elem_data7;
                delete [] elem_data8;
	  return ContactSearch::EXODUS_ERROR;
	}
	ierr = ex_put_elem_var( Exodus_ID, 1, 2, ID, num_faces, elem_data2 );
	if( ierr!=0 ){
	  std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
	  error_handler->Add_Error_Message(message);
                POSTCONDITION(ierr==0);
                delete [] elem_data1;
                delete [] elem_data2;
                delete [] elem_data3;
                delete [] elem_data4;
                delete [] elem_data5;
                delete [] elem_data6;
                delete [] elem_data7;
                delete [] elem_data8;
	  return ContactSearch::EXODUS_ERROR;
	}
	if( dimensionality == 3 ){
	  ierr = ex_put_elem_var( Exodus_ID, 1, 3, ID, num_faces, elem_data3 );
	  if( ierr!=0 ){
	    std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
	    error_handler->Add_Error_Message(message);
                POSTCONDITION(ierr==0);
                delete [] elem_data1;
                delete [] elem_data2;
                delete [] elem_data3;
                delete [] elem_data4;
                delete [] elem_data5;
                delete [] elem_data6;
                delete [] elem_data7;
                delete [] elem_data8;
	    return ContactSearch::EXODUS_ERROR;
	  }
	}
      }
    }
    
    if( dimensionality == 3 ){
      // Output the edge curvature & angle_bf
      Real* curvature = new Real[number_edges_owned];
      for( i=0 ; i<number_of_edge_blocks ; ++i){
        index = 0;
        ContactEdge<Real>** BlockEdges = 
          reinterpret_cast<ContactEdge<Real>**>(edge_list->BlockEntityList(i));
#ifndef CONTACT_NO_MPI
	if ( PARALLEL) {
          for (j=0; j<edge_list->BlockNumEntities(i); ++j) {
	    if ( BlockEdges[j]->Ownership() == ContactTopologyEntity<Real>::OWNED ) {
	      curvature[index] = *BlockEdges[j]->Variable(CURVATURE);
	      index++;
	    }
	  }
	} else {
#endif
          for (j=0; j<edge_blocks[i]->Number_of_Edges(); ++j) {
	    curvature[j] = *BlockEdges[j]->Variable(CURVATURE);
	  }
	  index = edge_blocks[i]->Number_of_Edges();
#ifndef CONTACT_NO_MPI
	}
#endif
	if (index>0) {
	  int ID = number_of_face_blocks+i+1;
	  ierr = ex_put_elem_var( Exodus_ID, 1, 4, ID, 
				  index, curvature );
	  if( ierr!=0 ){
	    std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
	    error_handler->Add_Error_Message(message);
                POSTCONDITION(ierr==0);
                delete [] elem_data1;
                delete [] elem_data2;
                delete [] elem_data3;
                delete [] elem_data4;
                delete [] elem_data5;
                delete [] elem_data6;
                delete [] elem_data7;
                delete [] elem_data8;
	    return ContactSearch::EXODUS_ERROR;
	  }
	}
      }
      delete [] curvature;
    }
    
    // Output the face-face interactions
    if (num_ffi>0) {
      for( i=0 ; i<number_of_face_blocks ; ++i){
	int num_faces = face_blocks[i]->Number_of_Faces();
	if (num_faces>0) {
	  n = 6;
	  int ID = i+1;
	  for (int k1=0; k1<max_ffi; ++k1) {
            ContactFace<Real>** BlockFaces = 
              reinterpret_cast<ContactFace<Real>**>(primary_face_list->BlockEntityList(i));
            for (j=0; j<primary_face_list->BlockNumEntities(i); ++j) {
	      int jj = 0;
              ContactInteractionDLL<Real>* interactions = BlockFaces[j]->Get_FaceFace_Interactions();
              if(interactions == NULL) continue;
              interactions->IteratorStart();
              while (interaction=interactions->IteratorForward()){
		if (jj==k1) break;
                jj++;
	      }
	      if (interaction) {
		cffi = static_cast<ContactFaceFaceInteraction<Real>*> 
		  (interaction);
		elem_data1[j] = cffi->MasterFaceEntityData()->index_in_owner_proc_array;
		elem_data2[j] = cffi->NumEdges()+1;
	      } else {
		elem_data1[j] = 0;
		elem_data2[j] = 0;
	      }
	    }
	    ierr=ex_put_elem_var(Exodus_ID, 1, n++, ID, num_faces, elem_data1);
	    if( ierr!=0 ){
	      std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
	      error_handler->Add_Error_Message(message);
                POSTCONDITION(ierr==0);
                delete [] elem_data1;
                delete [] elem_data2;
                delete [] elem_data3;
                delete [] elem_data4;
                delete [] elem_data5;
                delete [] elem_data6;
                delete [] elem_data7;
                delete [] elem_data8;
	      return ContactSearch::EXODUS_ERROR;
	    }
	    ierr=ex_put_elem_var(Exodus_ID, 1, n++, ID, num_faces, elem_data2);
	    if( ierr!=0 ){
	      std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
	      error_handler->Add_Error_Message(message);
                POSTCONDITION(ierr==0);
                delete [] elem_data1;
                delete [] elem_data2;
                delete [] elem_data3;
                delete [] elem_data4;
                delete [] elem_data5;
                delete [] elem_data6;
                delete [] elem_data7;
                delete [] elem_data8;
	      return ContactSearch::EXODUS_ERROR;
	    }
	    for (int k2=0; k2<max_ffi_verts; ++k2) {
              j = 0;
              face_blocks[i]->FaceList()->IteratorStart();
              while (entity=face_blocks[i]->FaceList()->IteratorForward()) {
                ContactFace<Real>* face = static_cast<ContactFace<Real>*>(entity);
		int jj = 0;
                ContactInteractionDLL<Real>* interactions = face->Get_FaceFace_Interactions();
                if(interactions != NULL) {
                  interactions->IteratorStart();
                  while (interaction=interactions->IteratorForward()){
		    if (jj==k1) break;
                    jj++;
		  }
		  if (interaction) {
		    cffi = static_cast<ContactFaceFaceInteraction<Real>*> (interaction);
  		    if (k2<cffi->NumEdges()+1) {
		      ContactFaceFaceVertex* vertices = cffi->Get_Vertices();
		      elem_data3[j] = vertices[k2].slave_x;
		      elem_data4[j] = vertices[k2].slave_y;
		      elem_data5[j] = vertices[k2].master_x;
		      elem_data6[j] = vertices[k2].master_y;
		      elem_data7[j] = vertices[k2].slave_edge_id;
		      elem_data8[j] = vertices[k2].master_edge_flag;
		    } else {
		      elem_data3[j] = 0;
		      elem_data4[j] = 0;
		      elem_data5[j] = 0;
		      elem_data6[j] = 0;
		      elem_data7[j] = 0;
		      elem_data8[j] = 0;
		    }
		  } else {
  		    elem_data3[j] = 0;
		    elem_data4[j] = 0;
		    elem_data5[j] = 0;
		    elem_data6[j] = 0;
		    elem_data7[j] = 0;
		    elem_data8[j] = 0;
		  }
		} else {
		  elem_data3[j] = 0;
		  elem_data4[j] = 0;
		  elem_data5[j] = 0;
		  elem_data6[j] = 0;
		  elem_data7[j] = 0;
		  elem_data8[j] = 0;
		}
                j++;
	      }
	      ierr=ex_put_elem_var(Exodus_ID,1,n++,ID,num_faces,elem_data3 );
	      if( ierr!=0 ){
		std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
		error_handler->Add_Error_Message(message);
                POSTCONDITION(ierr==0);
                delete [] elem_data1;
                delete [] elem_data2;
                delete [] elem_data3;
                delete [] elem_data4;
                delete [] elem_data5;
                delete [] elem_data6;
                delete [] elem_data7;
                delete [] elem_data8;
		return ContactSearch::EXODUS_ERROR;
	      }
	      ierr=ex_put_elem_var(Exodus_ID,1,n++,ID,num_faces,elem_data4 );
	      if( ierr!=0 ){
		std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
		error_handler->Add_Error_Message(message);
                POSTCONDITION(ierr==0);
                delete [] elem_data1;
                delete [] elem_data2;
                delete [] elem_data3;
                delete [] elem_data4;
                delete [] elem_data5;
                delete [] elem_data6;
                delete [] elem_data7;
                delete [] elem_data8;
		return ContactSearch::EXODUS_ERROR;
	      }
	      ierr=ex_put_elem_var(Exodus_ID,1,n++,ID,num_faces,elem_data5 );
	      if( ierr!=0 ){
		std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
		error_handler->Add_Error_Message(message);
                POSTCONDITION(ierr==0);
                delete [] elem_data1;
                delete [] elem_data2;
                delete [] elem_data3;
                delete [] elem_data4;
                delete [] elem_data5;
                delete [] elem_data6;
                delete [] elem_data7;
                delete [] elem_data8;
		return ContactSearch::EXODUS_ERROR;
	      }
	      ierr=ex_put_elem_var(Exodus_ID,1,n++,ID,num_faces,elem_data6 );
	      if( ierr!=0 ){
		std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
		error_handler->Add_Error_Message(message);
                POSTCONDITION(ierr==0);
                delete [] elem_data1;
                delete [] elem_data2;
                delete [] elem_data3;
                delete [] elem_data4;
                delete [] elem_data5;
                delete [] elem_data6;
                delete [] elem_data7;
                delete [] elem_data8;
		return ContactSearch::EXODUS_ERROR;
	      }
	      ierr=ex_put_elem_var(Exodus_ID,1,n++,ID,num_faces,elem_data7 );
	      if( ierr!=0 ){
		std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
		error_handler->Add_Error_Message(message);
                POSTCONDITION(ierr==0);
                delete [] elem_data1;
                delete [] elem_data2;
                delete [] elem_data3;
                delete [] elem_data4;
                delete [] elem_data5;
                delete [] elem_data6;
                delete [] elem_data7;
                delete [] elem_data8;
		return ContactSearch::EXODUS_ERROR;
	      }
	      ierr=ex_put_elem_var(Exodus_ID,1,n++,ID,num_faces,elem_data8 );
	      if( ierr!=0 ){
		std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
		error_handler->Add_Error_Message(message);
                POSTCONDITION(ierr==0);
                delete [] elem_data1;
                delete [] elem_data2;
                delete [] elem_data3;
                delete [] elem_data4;
                delete [] elem_data5;
                delete [] elem_data6;
                delete [] elem_data7;
                delete [] elem_data8;
		return ContactSearch::EXODUS_ERROR;
	      }
	    }
	  }
	}
      }
    }
    // Output the face-coverage interactions
    if (num_fci>0) {
      for( i=0 ; i<number_of_face_blocks ; ++i){
	int num_faces = face_blocks[i]->Number_of_Faces();
	if (num_faces>0) {
	  n = 6+max_ffi*(2+6*max_ffi_verts);
	  int ID = i+1;
	  for (int k1=0; k1<max_fci; ++k1) {
            ContactFace<Real>** BlockFaces = 
              reinterpret_cast<ContactFace<Real>**>(primary_face_list->BlockEntityList(i));
            for (j=0; j<primary_face_list->BlockNumEntities(i); ++j) {
	      int jj = 0;
              ContactInteractionDLL<Real>* interactions = BlockFaces[j]->Get_FaceCoverage_Interactions();
              if(interactions != NULL) {
                interactions->IteratorStart();
                while (interaction=interactions->IteratorForward()){
	  	  if (jj==k1) break;
                  jj++;
	        }
	        if (interaction) {
		  ContactFaceCoverageInteraction* cfci =
		    static_cast<ContactFaceCoverageInteraction*> (interaction);
		  elem_data1[j] = cfci->NumVertices();
	        } else {
		  elem_data1[j] = 0;
	        }
	      } else {
                elem_data1[j] = 0;
	      }
	    }
	    ierr = ex_put_elem_var( Exodus_ID, 1, n++, ID, num_faces, elem_data1 );
	    if( ierr!=0 ){
	      std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
	      error_handler->Add_Error_Message(message);
                POSTCONDITION(ierr==0);
                delete [] elem_data1;
                delete [] elem_data2;
                delete [] elem_data3;
                delete [] elem_data4;
                delete [] elem_data5;
                delete [] elem_data6;
                delete [] elem_data7;
                delete [] elem_data8;
	      return ContactSearch::EXODUS_ERROR;
	    }
	    for (int k2=0; k2<max_fci_verts; ++k2) {
              for (j=0; j<primary_face_list->BlockNumEntities(i); ++j) {
		int jj = 0;
                ContactInteractionDLL<Real>* interactions = BlockFaces[j]->Get_FaceCoverage_Interactions();
                if(interactions != NULL) {
		  interactions->IteratorStart();
		  while (interaction=interactions->IteratorForward()){
		    if (jj==k1) break;
		    jj++;
		  }
		  if (interaction) {
		    ContactFaceCoverageInteraction* cfci =
		      static_cast<ContactFaceCoverageInteraction*> (interaction);
		    if (k2<cfci->NumVertices()) {
		      int kk = 0;
		      ContactFaceCoverageVertex* ll_node;
		      for( ll_node=cfci->Head(); ll_node; 
			   ll_node=ll_node->next, kk++ ){
			if (kk==k2) break;
		      }
		      elem_data2[j] = ll_node->slave_x;
		      elem_data3[j] = ll_node->slave_y;
		    } else {
		      elem_data2[j] = 0;
		      elem_data3[j] = 0;
		    }
		  } else {
		    elem_data2[j] = 0;
		    elem_data3[j] = 0;
		  }
		} else {
                  elem_data2[j] = 0;
		  elem_data3[j] = 0;
		}
	      }
	      ierr = ex_put_elem_var( Exodus_ID, 1, n++, ID, num_faces, elem_data2 );
	      if( ierr!=0 ){
		std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
		error_handler->Add_Error_Message(message);
                POSTCONDITION(ierr==0);
                delete [] elem_data1;
                delete [] elem_data2;
                delete [] elem_data3;
                delete [] elem_data4;
                delete [] elem_data5;
                delete [] elem_data6;
                delete [] elem_data7;
                delete [] elem_data8;
		return ContactSearch::EXODUS_ERROR;
	      }
	      ierr = ex_put_elem_var( Exodus_ID, 1, n++, ID, num_faces, elem_data3 );
	      if( ierr!=0 ){
		std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
		error_handler->Add_Error_Message(message);
                POSTCONDITION(ierr==0);
                delete [] elem_data1;
                delete [] elem_data2;
                delete [] elem_data3;
                delete [] elem_data4;
                delete [] elem_data5;
                delete [] elem_data6;
                delete [] elem_data7;
                delete [] elem_data8;
		return ContactSearch::EXODUS_ERROR;
	      }
	    }
	  }
	}
      }
    }
    delete [] elem_data1;
    delete [] elem_data2;
    delete [] elem_data3;
    delete [] elem_data4;
    delete [] elem_data5;
    delete [] elem_data6;
    delete [] elem_data7;
    delete [] elem_data8;

#ifndef CONTACT_NO_MPI
    int max_number = std::max(number_of_elems,number_edges_owned);
    Real* elem_data = new Real[max_number];
    // Output the face/edge primary owner
    for( i=0 ; i<number_of_face_blocks ; ++i){
      int num_faces = primary_face_list->BlockNumEntities(i);
      ContactFace<Real>** BlockFaces = 
        reinterpret_cast<ContactFace<Real>**>(primary_face_list->BlockEntityList(i));
      for (j=0; j<num_faces; ++j) {
	elem_data[j] = BlockFaces[j]->Owner();
      }
      if (num_faces>0) {
	int ID = i+1;
	ierr = ex_put_elem_var( Exodus_ID, 1, num_elem_vars-2, 
				ID, num_faces, elem_data );
	if( ierr!=0 ){
	  std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
	  error_handler->Add_Error_Message(message);
	  return ContactSearch::EXODUS_ERROR;
	}
      }
    }
    for( i=0 ; i<number_of_edge_blocks ; ++i){
      index = 0;
      ContactEdge<Real>** BlockEdges = 
        reinterpret_cast<ContactEdge<Real>**>(edge_list->BlockEntityList(i));
      for (j=0; j<edge_list->BlockNumEntities(i); ++j) {
	if ( BlockEdges[j]->Ownership() == ContactTopologyEntity<Real>::OWNED ) {
	  elem_data[index] = BlockEdges[j]->Owner();
	  index++;
	}
      }
      if (index>0) {
	int ID = number_of_face_blocks+i+1;
	ierr = ex_put_elem_var( Exodus_ID, 1, num_elem_vars-2, 
				ID, index, elem_data );
	if( ierr!=0 ){
	  std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
	  error_handler->Add_Error_Message(message);
	  return ContactSearch::EXODUS_ERROR;
	}
      }
    }
    
    // Output the face/edge primary local_id
    for( i=0 ; i<number_of_face_blocks ; ++i){
      int num_faces = primary_face_list->BlockNumEntities(i);
      ContactFace<Real>** BlockFaces = 
        reinterpret_cast<ContactFace<Real>**>(primary_face_list->BlockEntityList(i));
      for (j=0; j<num_faces; ++j) {
	//elem_data[j] = BlockFaces[j]->ProcArrayIndex();
	elem_data[j] = BlockFaces[j]->fcs_index;
      }
      if (num_faces>0) {
	int ID = i+1;
	ierr = ex_put_elem_var( Exodus_ID, 1, num_elem_vars-1, 
				ID, num_faces, elem_data );
	if( ierr!=0 ){
	  std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
	  error_handler->Add_Error_Message(message);
	  return ContactSearch::EXODUS_ERROR;
	}
      }
    }
    for( i=0 ; i<number_of_edge_blocks ; ++i){
      index = 0;
      ContactEdge<Real>** BlockEdges = 
        reinterpret_cast<ContactEdge<Real>**>(edge_list->BlockEntityList(i));
      for (j=0; j<edge_list->BlockNumEntities(i); ++j) {
	if ( BlockEdges[j]->Ownership() == ContactTopologyEntity<Real>::OWNED ) {
	  //elem_data[index] = BlockEdges[j]->ProcArrayIndex();
	  elem_data[index] = BlockEdges[j]->fcs_index;
	  index++;
	}
      }
      if (index>0) {
	int ID = number_of_face_blocks+i+1;
	ierr = ex_put_elem_var( Exodus_ID, 1, num_elem_vars-1, 
				ID, index, elem_data );
	if( ierr!=0 ){
	  std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
	  error_handler->Add_Error_Message(message);
	  return ContactSearch::EXODUS_ERROR;
	}
      }
    }
    
    // Output the face/edge secondary owner
    for( i=0 ; i<number_of_face_blocks ; ++i){
      int num_faces = primary_face_list->BlockNumEntities(i);
      ContactFace<Real>** BlockFaces = 
        reinterpret_cast<ContactFace<Real>**>(primary_face_list->BlockEntityList(i));
      for (j=0; j<num_faces; ++j) {
	elem_data[j] = BlockFaces[j]->Secondary_Owner();
      }
      if (num_faces>0) {
	int ID = i+1;
	ierr = ex_put_elem_var( Exodus_ID, 1, num_elem_vars, 
				ID, num_faces, elem_data );
	if( ierr!=0 ){
	  std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
	  error_handler->Add_Error_Message(message);
	  return ContactSearch::EXODUS_ERROR;
	}
      }
    }
    for( i=0 ; i<number_of_edge_blocks ; ++i){
      index = 0;
      ContactEdge<Real>** BlockEdges = 
        reinterpret_cast<ContactEdge<Real>**>(edge_list->BlockEntityList(i));
      for (j=0; j<edge_list->BlockNumEntities(i); ++j) {
	if ( BlockEdges[j]->Ownership() == ContactTopologyEntity<Real>::OWNED ) {
	  elem_data[index] = BlockEdges[j]->Secondary_Owner();
	  index++;
	}
      }
      if (index>0) {
	int ID = number_of_face_blocks+i+1;
	ierr = ex_put_elem_var( Exodus_ID, 1, num_elem_vars, 
				ID, index, elem_data );
	if( ierr!=0 ){
	  std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
	  error_handler->Add_Error_Message(message);
	  return ContactSearch::EXODUS_ERROR;
	}
      }
    }
    delete [] elem_data;
#endif
  } else if( number_of_element_blocks ){
    Real* elem_data1 = new Real[number_of_elems];
    Real* elem_data2 = new Real[number_of_elems];
    for( i=0 ; i<number_of_element_blocks ; ++i){
      int num_elements = primary_elem_list->BlockNumEntities(i);
      ContactElement** BlockElements = 
        reinterpret_cast<ContactElement**>(primary_elem_list->BlockEntityList(i));
      for (j=0; j<num_elements; ++j) {
	Real volume = *(BlockElements[j]->Variable( ELEMENT_VOLUME ));
	elem_data1[j] = volume;
      }
      if (num_elements>0) {
	int ID = i+1;
	ierr = ex_put_elem_var( Exodus_ID, 1, 1, 
				ID, num_elements, elem_data1 );
	if( ierr!=0 ){
	  std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
	  error_handler->Add_Error_Message(message);
	  return ContactSearch::EXODUS_ERROR;
	}
      }
    }
    
    for( i=0 ; i<number_of_element_blocks ; ++i){
      int num_elements = primary_elem_list->BlockNumEntities(i);
      ContactElement** BlockElements = 
        reinterpret_cast<ContactElement**>(primary_elem_list->BlockEntityList(i));
      for (j=0; j<num_elements; ++j) {
        int jj = BlockElements[j]->Get_ElementElement_Interactions()->NumEntities();
	elem_data1[j] = jj;
      }
      if (num_elements>0) {
	int ID = i+1;
	ierr = ex_put_elem_var( Exodus_ID, 1, 2, 
				ID, num_elements, elem_data1 );
	if( ierr!=0 ){
	  std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
	  error_handler->Add_Error_Message(message);
	  return ContactSearch::EXODUS_ERROR;
	}
      }
    }

    for( i=0; i<num_element_enf_variables; ++i) {
      enforcement->Get_Element_Plot_Variable( i, elem_data1 );
      int element_count = 0;
      for( j=0 ; j<number_of_element_blocks ; ++j){
	int num_elements = element_blocks[j]->Number_of_Elements();
	if (num_elements>0) {
	  int ID = j+1;
	  ierr = ex_put_elem_var( Exodus_ID, 1, 3+i, 
				  ID, num_elements, &elem_data1[element_count] );
	  if( ierr!=0 ){
	    std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
	    error_handler->Add_Error_Message(message);
	    return ContactSearch::EXODUS_ERROR;
	  }
	}
	element_count += num_elements;
      }
    }      
    
    // Output the element-element interactions
    if (num_eei>0) {
      for( i=0 ; i<number_of_element_blocks ; ++i){
	int num_elements = primary_elem_list->BlockNumEntities(i);
	if (num_elements>0) {
	  n = 3+num_element_enf_variables;
	  int ID = i+1;
	  for (int k1=0; k1<max_eei; ++k1) {
            ContactElement** BlockElements = 
              reinterpret_cast<ContactElement**>(primary_elem_list->BlockEntityList(i));
            for (j=0; j<num_elements; ++j) {
	      int jj = 0;
              ContactInteractionDLL<Real>* interactions = BlockElements[j]->Get_ElementElement_Interactions();
              interactions->IteratorStart();
              while (interaction=interactions->IteratorForward()){
		if (jj==k1) break;
                jj++;
	      }
	      if (interaction) {
		ContactElementElementInteraction* ceei = 
                    static_cast<ContactElementElementInteraction*>(interaction);
		elem_data1[j] = ceei->MasterElementEntityData()->index_in_owner_proc_array;
		elem_data2[j] = ceei->Scalar_Var(ContactElementElementInteraction::VOLUME);
	      } else {
		elem_data1[j] = 0;
		elem_data2[j] = 0;
	      }
	    }
	    ierr=ex_put_elem_var(Exodus_ID, 1, n++, ID, num_elements, 
				 elem_data1);
	    if( ierr!=0 ){
	      std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
	      error_handler->Add_Error_Message(message);
	      return ContactSearch::EXODUS_ERROR;
	    }
	    ierr=ex_put_elem_var(Exodus_ID, 1, n++, ID, num_elements, 
				 elem_data2);
	    if( ierr!=0 ){
	      std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
	      error_handler->Add_Error_Message(message);
	      return ContactSearch::EXODUS_ERROR;
	    }
	  }
	}
      }
    }
    
#ifndef CONTACT_NO_MPI
    
    // Output the element primary owner
    for( i=0 ; i<number_of_element_blocks ; ++i){
      int num_elements = primary_elem_list->BlockNumEntities(i);
      ContactElement** BlockElements = 
        reinterpret_cast<ContactElement**>(primary_elem_list->BlockEntityList(i));
      for (j=0; j<num_elements; ++j) {
	elem_data1[j] = BlockElements[j]->Owner();
      }
      if (num_elements>0) {
	int ID = i+1;
	ierr = ex_put_elem_var( Exodus_ID, 1, num_elem_vars-2, 
				ID, num_elements, elem_data1 );
	if( ierr!=0 ){
	  std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
	  error_handler->Add_Error_Message(message);
	  return ContactSearch::EXODUS_ERROR;
	}
      }
    }
    
    // Output the element primary local_id
    for( i=0 ; i<number_of_element_blocks ; ++i){
      int num_elements = primary_elem_list->BlockNumEntities(i);
      ContactElement** BlockElements = 
        reinterpret_cast<ContactElement**>(primary_elem_list->BlockEntityList(i));
      for (j=0; j<num_elements; ++j) {
	//elem_data1[j] = BlockElements[j]->ProcArrayIndex();
	elem_data1[j] = BlockElements[j]->fcs_index;
      }
      if (num_elements>0) {
	int ID = i+1;
	ierr = ex_put_elem_var( Exodus_ID, 1, num_elem_vars-1, 
				ID, num_elements, elem_data1 );
	if( ierr!=0 ){
	  std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
	  error_handler->Add_Error_Message(message);
	  return ContactSearch::EXODUS_ERROR;
	}
      }
    }
    
    // Output the element secondary owner
    for( i=0 ; i<number_of_element_blocks ; ++i){
      int num_elements = primary_elem_list->BlockNumEntities(i);
      ContactElement** BlockElements = 
        reinterpret_cast<ContactElement**>(primary_elem_list->BlockEntityList(i));
      for (j=0; j<num_elements; ++j) {
	elem_data1[j] = BlockElements[j]->Secondary_Owner();
      }
      if (num_elements>0) {
	int ID = i+1;
	ierr = ex_put_elem_var( Exodus_ID, 1, num_elem_vars, 
				ID, num_elements, elem_data1 );
	if( ierr!=0 ){
	  std::sprintf(message,"Exodus Error %d from ex_put_elem_var",ierr);
	  error_handler->Add_Error_Message(message);
	  return ContactSearch::EXODUS_ERROR;
	}
      }
    }
#endif
    delete [] elem_data1;
    delete [] elem_data2;
  }
  return ContactSearch::NO_ERROR;
}

#endif
