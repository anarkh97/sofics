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


#include "ContactTopology.h"
#include "ContactSearchData.h"
#include "ContactBoundingBox.h"
#include "ContactParOStream.h"

#include <cstddef>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>

ContactSearchData::ContactSearchData(ContactTopology* topology)
{
  num_elem_blocks       = topology->Number_of_Element_Blocks();
  num_face_blocks       = topology->Number_of_Face_Blocks();
  if (num_elem_blocks+num_face_blocks) {
    num_node_blocks     = topology->Number_of_Node_Blocks()-1;
  } else {
    num_node_blocks     = topology->Number_of_Node_Blocks();
  }
  num_analytic_surfaces = topology->Number_of_Analytic_Surfaces();
  num_search_entities   = num_elem_blocks + 
                          num_face_blocks + 
                          num_node_blocks + 
                          num_analytic_surfaces;
  size_search_data = ContactSearch::NSIZSD;

  data_array = 
    new Real[size_search_data*num_search_entities*num_search_entities];
  interaction = 
    new int[num_search_entities*num_search_entities];
  intersection = 
    new ContactBoundingBox[num_search_entities*num_search_entities];
    
  max_search_tolerance              = 0.0;
  max_search_normal_tolerance       = 0.0;
  max_search_tangential_tolerance   = 0.0;
  search_for_interactions           = false;
  have_tied_interactions            = false;
  have_only_tied_interactions       = false;
  have_node_node_interactions       = false;
  have_node_face_interactions       = false;
  have_face_face_interactions       = false;
  have_face_coverage_interactions   = false;
  have_element_element_interactions = false;
}

ContactSearchData::~ContactSearchData(){
  if( data_array )   delete [] data_array;
  if( interaction )  delete [] interaction;
  if( intersection ) delete [] intersection;
}

void ContactSearchData::Set_Search_Data( const Real* data )
{

  std::memcpy( data_array,data,size_search_data*num_search_entities*
	       num_search_entities*sizeof(Real) );
  Compute_Max_Search_Tolerances();
  Set_Tied_Interaction_Status();
  Set_Search_Types();
  if (have_node_node_interactions      ||
      have_node_face_interactions      ||
      have_face_face_interactions      ||
      have_face_coverage_interactions  ||
      have_element_element_interactions ) search_for_interactions = true;
}

void ContactSearchData::Compute_Max_Search_Tolerances()
{
  max_search_tolerance = 0;
  max_search_normal_tolerance = 0;
  max_search_tangential_tolerance = 0;
  for( int i=0 ; i<num_search_entities ; ++i){
    for( int j=0 ; j<num_search_entities ; ++j){
      Real norm_tol = Get_Search_Data( ContactSearch::SEARCH_NORMAL_TOLERANCE,i,j );
      max_search_normal_tolerance = std::max(max_search_normal_tolerance,norm_tol);
      Real tan_tol = Get_Search_Data( ContactSearch::SEARCH_TANGENTIAL_TOLERANCE,i,j );
      max_search_tangential_tolerance = std::max(max_search_tangential_tolerance,tan_tol);
    }
  }
  max_search_tolerance = std::max(max_search_normal_tolerance, 
                             max_search_tangential_tolerance);
}

void ContactSearchData::Set_Tied_Interaction_Status()
{
  have_only_tied_interactions = true;
  for( int i=0 ; i<num_search_entities ; ++i){
    for( int j=0 ; j<num_search_entities ; ++j){
      int type = (int)Get_Search_Data( ContactSearch::INTERACTION_TYPE, i, j );
      interaction[j*num_search_entities+i] = type;
      if( !(type == ContactSearch::TIED_INTERACTION ||
            type == ContactSearch::INFINITESIMAL_SLIP_INTERACTION ||
            type == ContactSearch::NO_INTERACTION)) {
        have_only_tied_interactions = false;
      }
      if( type == ContactSearch::TIED_INTERACTION ||
          type == ContactSearch::INFINITESIMAL_SLIP_INTERACTION ) {
        have_tied_interactions = true;
      }
    }
  }
}

void ContactSearchData::Set_Search_Types()
{
  int i, j;
  int slave;
  int master;
  have_node_node_interactions = false;
  have_node_face_interactions = false;
  have_face_face_interactions = false;
  have_face_coverage_interactions = false;
  have_element_element_interactions = false;
  if (num_elem_blocks>0) {
    slave = 0;
    for( i=0 ; i<num_elem_blocks ; i++, slave++ ){
      master = 0;
      for( j=0 ; j<num_elem_blocks ; j++, master++ ){
   	if( Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) == 
   	    ContactSearch::GENERIC_INTERACTION ){
   	  have_element_element_interactions = true;
   	  break;
   	}
      }
      if (have_element_element_interactions) break;
    }
  }
  if (num_face_blocks>0) {
    slave = num_elem_blocks;
    for( i=0 ; i<num_face_blocks ; i++, slave++ ){
      master = num_elem_blocks;
      for( j=0 ; j<num_face_blocks ; j++, master++ ){
   	if( Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) == 
   	    ContactSearch::COVERAGE_INTERACTION ){
   	  have_face_face_interactions = true;
   	  have_face_coverage_interactions = true;
 	  break;
   	}
   	if( Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) == 
   	    ContactSearch::GENERIC_INTERACTION ){
   	  have_face_face_interactions = true;
   	}
      }
      if (have_face_coverage_interactions) break;
    }
  }
  if (num_face_blocks) {
    slave = num_elem_blocks;
    for( i=0 ; i<num_face_blocks+num_node_blocks ; i++, slave++ ){
      master = num_elem_blocks;
      for( j=0 ; j<num_face_blocks+num_node_blocks ; j++, master++ ){
        if( Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
              ContactSearch::TIED_INTERACTION ||
            Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
              ContactSearch::INFINITESIMAL_SLIP_INTERACTION ||
            Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
              ContactSearch::GLUED_INTERACTION ||
            Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
              ContactSearch::SLIDING_INTERACTION ||
            Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
              ContactSearch::CPP_INTERACTION ){
          have_node_face_interactions = true;
          break;
        }
      }
      if (have_node_face_interactions) break;
      master = num_elem_blocks+num_face_blocks+num_node_blocks;
      for( j=0 ; j<num_analytic_surfaces ; j++, master++ ){
        if( Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
              ContactSearch::TIED_INTERACTION ||
            Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
              ContactSearch::INFINITESIMAL_SLIP_INTERACTION ||
            Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
              ContactSearch::GLUED_INTERACTION ||
            Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
              ContactSearch::SLIDING_INTERACTION ||
            Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
              ContactSearch::CPP_INTERACTION ){
          have_node_face_interactions = true;
          break;
        }
      }
      if (have_node_face_interactions) break;
    }
  }
  if (num_node_blocks) {
    slave = num_elem_blocks+num_face_blocks;
    for( i=0 ; i<num_node_blocks ; i++, slave++ ){
      master = num_elem_blocks+num_face_blocks;
      for( j=0 ; j<num_node_blocks ; j++, master++ ){
        if( Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) == 
            ContactSearch::GENERIC_INTERACTION ){
          have_node_node_interactions = true;
          break;
        }
      }
      if (have_node_node_interactions) break;
    }
  }
}

bool ContactSearchData::Master_Has_Only_Tied_Interactions(int face_key)
{
  bool status = true;
  for( int i=0 ; i<num_search_entities ; ++i){
    int type = (int)
      Get_Search_Data( ContactSearch::INTERACTION_TYPE, i, face_key );
    if( type == ContactSearch::NO_INTERACTION ||
	type == ContactSearch::TIED_INTERACTION ||
	type == ContactSearch::INFINITESIMAL_SLIP_INTERACTION) continue;
    status = false;
    break;
  }
  return status;
}

bool ContactSearchData::Slave_Has_Only_Tied_Interactions(int node_key)
{
  bool status = true;
  for( int i=0 ; i<num_search_entities ; ++i){
    int type = (int)
      Get_Search_Data( ContactSearch::INTERACTION_TYPE, node_key, i );
    if( type == ContactSearch::NO_INTERACTION ||
	type == ContactSearch::TIED_INTERACTION ||
	type == ContactSearch::INFINITESIMAL_SLIP_INTERACTION) continue;
    status = false;
    break;
  }
  return status;
}

bool ContactSearchData::Is_NodeBlock_Slave(int node_block)
{
  return (Is_NodeBlock_NodeNodeSlave(node_block) || Is_NodeBlock_NodeFaceSlave(node_block));
}

bool ContactSearchData::Is_FaceBlock_Slave(int face_block)
{
  return (Is_FaceBlock_NodeFaceSlave(face_block) || Is_FaceBlock_FaceFaceSlave(face_block));
}

bool ContactSearchData::Is_ElemBlock_Slave(int elem_block)
{
  return (Is_ElemBlock_ElemElemSlave(elem_block));
}

bool ContactSearchData::Is_NodeBlock_Master(int node_block)
{
  return (Is_NodeBlock_NodeNodeMaster(node_block));
}

bool ContactSearchData::Is_FaceBlock_Master(int face_block)
{
  return (Is_FaceBlock_NodeFaceMaster(face_block) || Is_FaceBlock_FaceFaceMaster(face_block));
}

bool ContactSearchData::Is_ElemBlock_Master(int elem_block)
{
  return (Is_ElemBlock_ElemElemMaster(elem_block));
}

bool ContactSearchData::Is_NodeBlock_NodeNodeSlave(int node_block)
{
  if (num_node_blocks>0) {
    int j;
    int offset = (num_elem_blocks+num_face_blocks)>0?1:0;
    int slave = num_elem_blocks+num_face_blocks+node_block-offset;
    for( j=0 ; j<num_node_blocks ; ++j){
      int master = num_elem_blocks+num_face_blocks+j;
      if( Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::GENERIC_INTERACTION ){
        return true;
      }
    }
  }
  return false;
}

bool ContactSearchData::Is_NodeBlock_NodeFaceSlave(int node_block)
{
  if (num_node_blocks>0 && node_block>0) {
    int j;
    int offset = (num_elem_blocks+num_face_blocks)>0?1:0;
    int slave  = num_elem_blocks+num_face_blocks+node_block-offset;
    for( j=0 ; j<num_face_blocks ; ++j){
      int master = num_elem_blocks+j;
      if( Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::TIED_INTERACTION ||
          Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::INFINITESIMAL_SLIP_INTERACTION ||
          Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::GLUED_INTERACTION ||
          Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::SLIDING_INTERACTION ||
          Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::CPP_INTERACTION ){
        return true;
      }
    }
    for( j=0 ; j<num_analytic_surfaces ; ++j){
      int master = num_elem_blocks+num_face_blocks+num_node_blocks+j;
      if( Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::TIED_INTERACTION ||
          Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::INFINITESIMAL_SLIP_INTERACTION ||
          Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::GLUED_INTERACTION ||
          Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::SLIDING_INTERACTION ||
          Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::CPP_INTERACTION ){
        return true;
      }
    }
  }
  return false;
}

bool ContactSearchData::Is_FaceBlock_NodeFaceSlave(int face_block)
{
  if (num_face_blocks) {
    int j;
    int slave = num_elem_blocks+face_block;
    for( j=0 ; j<num_face_blocks+num_node_blocks ; ++j){
      int master = num_elem_blocks+j;
      if( Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::TIED_INTERACTION ||
          Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::INFINITESIMAL_SLIP_INTERACTION ||
          Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::GLUED_INTERACTION ||
          Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::SLIDING_INTERACTION ||
          Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::CPP_INTERACTION ){
        return true;
      }
    }
    for( j=0 ; j<num_analytic_surfaces ; ++j){
      int master = num_elem_blocks+num_face_blocks+num_node_blocks+j;
      if( Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::TIED_INTERACTION ||
          Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::INFINITESIMAL_SLIP_INTERACTION ||
          Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::GLUED_INTERACTION ||
          Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::SLIDING_INTERACTION ||
          Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::CPP_INTERACTION ){
        return true;
      }
    }
  }
  return false;
}

bool ContactSearchData::Is_FaceBlock_FaceFaceSlave(int face_block)
{
  if (num_face_blocks) {
    int slave = num_elem_blocks+face_block;
    for( int j=0 ; j<num_face_blocks ; ++j){
      int master = num_elem_blocks+j;
      if( Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
          ContactSearch::GENERIC_INTERACTION ) return true;
    }
  }
  return false;
}

bool ContactSearchData::Is_ElemBlock_ElemElemSlave(int elem_block)
{
  if (num_elem_blocks) {
    int slave = elem_block;
    for( int j=0 ; j<num_elem_blocks ; ++j){
      int master = j;
      if( Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
          ContactSearch::GENERIC_INTERACTION ) return true;
    }
  }
  return false;
}

bool ContactSearchData::Is_NodeBlock_NodeNodeMaster(int node_block)
{
  if (num_node_blocks>0) {
    int j;
    int offset = (num_elem_blocks+num_face_blocks)>0?1:0;
    int master = num_elem_blocks+num_face_blocks+node_block-offset;
    for( j=0 ; j<num_node_blocks ; ++j){
      int slave = num_elem_blocks+num_face_blocks+j;
      if( Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::GENERIC_INTERACTION ){
        return true;
      }
    }
  }
  return false;
}

bool ContactSearchData::Is_FaceBlock_NodeFaceMaster(int face_block)
{
  if (num_face_blocks) {
    int j;
    int master = num_elem_blocks+face_block;
    for( j=0 ; j<num_face_blocks+num_node_blocks ; ++j){
      int slave = num_elem_blocks+j;
      if( Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::TIED_INTERACTION ||
          Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::INFINITESIMAL_SLIP_INTERACTION ||
          Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::GLUED_INTERACTION ||
          Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::SLIDING_INTERACTION ||
          Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
            ContactSearch::CPP_INTERACTION ){
        return true;
      }
    }
  }
  return false;
}

bool ContactSearchData::Is_FaceBlock_FaceFaceMaster(int face_block)
{
  if (num_face_blocks) {
    int master = num_elem_blocks+face_block;
    for( int j=0 ; j<num_face_blocks ; ++j){
      int slave = num_elem_blocks+j;
      if( Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
          ContactSearch::GENERIC_INTERACTION ) return true;
    }
  }
  return false;
}

bool ContactSearchData::Is_ElemBlock_ElemElemMaster(int master)
{
  if (num_elem_blocks) {
    for( int slave=0 ; slave<num_elem_blocks ; slave++ ){
      if( Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master ) ==
          ContactSearch::GENERIC_INTERACTION ) return true;
    }
  }
  return false;
}

void ContactSearchData::SetAllButTied()
{
  search_for_interactions = false;
  for( int i=0 ; i<num_search_entities ; ++i){
    for( int j=0 ; j<num_search_entities ; ++j){
      int index = j*num_search_entities+i;
      if( interaction[index] == ContactSearch::TIED_INTERACTION ||
          interaction[index] == ContactSearch::INFINITESIMAL_SLIP_INTERACTION) {  
        Set_Search_Data( ContactSearch::INTERACTION_TYPE, i, j, ContactSearch::NO_INTERACTION );
      } else {
        Set_Search_Data( ContactSearch::INTERACTION_TYPE, i, j, interaction[index] );
        if (interaction[index] != ContactSearch::NO_INTERACTION) search_for_interactions = true;
      }
    }
  }
}

void ContactSearchData::SetOnlyTied()
{
  search_for_interactions = false;
  for( int i=0 ; i<num_search_entities ; ++i){
    for( int j=0 ; j<num_search_entities ; ++j){
      int index = j*num_search_entities+i;
      if( interaction[index] == ContactSearch::TIED_INTERACTION ||
          interaction[index] == ContactSearch::INFINITESIMAL_SLIP_INTERACTION) {     
        Set_Search_Data( ContactSearch::INTERACTION_TYPE, i, j, interaction[index] );
        search_for_interactions = true;
      } else {
        Set_Search_Data( ContactSearch::INTERACTION_TYPE, i, j, ContactSearch::NO_INTERACTION );
      }
    }
  }
}

void ContactSearchData::SetAll()
{
  search_for_interactions = false;
  for( int i=0 ; i<num_search_entities ; ++i){
    for( int j=0 ; j<num_search_entities ; ++j){
      int index = j*num_search_entities+i;
      Set_Search_Data( ContactSearch::INTERACTION_TYPE, i, j, interaction[index] );
      if (interaction[index] != ContactSearch::NO_INTERACTION) search_for_interactions = true;
    }
  }
}

void ContactSearchData::Display(ContactParOStream& postream)
{
#if CONTACT_DEBUG_PRINT_LEVEL>=2
  const char* type[8]  = { "NO_INTERACTION",
		     "SLIDING_INTERACTION",
		     "TIED_INTERACTION",
                     "COVERAGE_INTERACTION",
                     "GENERIC_INTERACTION",
                     "CPP_INTERACTION",
                     "INFINITESIMAL_SLIP_INTERACTION",
                     "GLUED_INTERACTION" };
  postream << "Search Data:" << "\n";
  postream << "  Number of Entity Keys = " << num_search_entities << "\n";
  postream << "    Number of Element Blocks    = " << num_elem_blocks << "\n";
  postream << "    Number of Face Blocks       = " << num_face_blocks << "\n";
  postream << "    Number of Node Blocks       = " << num_node_blocks << "\n";
  postream << "    Number of Analytic Surfaces = " << num_analytic_surfaces << "\n";
  for( int slave=0 ; slave<num_search_entities ; slave++ ){
    for( int master=0 ; master<num_search_entities ; master++ ){
      postream << "    Data for Slave Entity " << slave+1 
	       << " against Master Entity " << master+1 << "\n";
      postream << "        Interaction Type     = " 
	       << type[(int)Get_Search_Data( ContactSearch::INTERACTION_TYPE, slave, master )] << "\n";
      postream << "        Normal Tolerance     = " 
	       << Get_Search_Data( ContactSearch::SEARCH_NORMAL_TOLERANCE, slave, master ) << "\n";
      postream << "        Tangential Tolerance = " 
	       << Get_Search_Data( ContactSearch::SEARCH_TANGENTIAL_TOLERANCE, slave, master ) << "\n";
    }
  }
  postream<<"    have_only_tied_interactions       = "<<have_only_tied_interactions<<"\n";
  postream<<"    have_node_node_interactions       = "<<have_node_node_interactions<<"\n";
  postream<<"    have_node_face_interactions       = "<<have_node_face_interactions<<"\n";
  postream<<"    have_face_face_interactions       = "<<have_face_face_interactions<<"\n";
  postream<<"    have_face_coverage_interactions   = "<<have_face_coverage_interactions<<"\n";
  postream<<"    have_element_element_interactions = "<<have_element_element_interactions<<"\n";
#endif
}


void ContactSearchData::Display0()
{
  const char* type[8]  = { "NO_INTERACTION",
		     "SLIDING_INTERACTION",
		     "TIED_INTERACTION",
                     "COVERAGE_INTERACTION",
                     "GENERIC_INTERACTION",
                     "CPP_INTERACTION",
                     "INFINITESIMAL_SLIP_INTERACTION",
                     "GLUED_INTERACTION" };
  std::cout << "Search Data:" << std::endl;
  std::cout << "  Number of Entity Keys = " << num_search_entities << std::endl;
  std::cout << "    Number of Element Blocks    = " << num_elem_blocks << std::endl;
  std::cout << "    Number of Face Blocks       = " << num_face_blocks << std::endl;
  std::cout << "    Number of Node Blocks       = " << num_node_blocks << std::endl;
  std::cout << "    Number of Analytic Surfaces = " << num_analytic_surfaces << std::endl;
  for( int slave=0 ; slave<num_search_entities ; slave++ ){
    for( int master=0 ; master<num_search_entities ; master++ ){
      std::cout << "    Data for Slave Entity " << slave+1 
	        << " against Master Entity " << master+1 << std::endl;
      std::cout << "        Interaction Type     = " 
	        << type[(int)Get_Search_Data( ContactSearch::INTERACTION_TYPE, 
			       slave, master )] << std::endl;;
      std::cout << "        Normal Tolerance     = " 
	        << Get_Search_Data( ContactSearch::SEARCH_NORMAL_TOLERANCE, 
			       slave, master ) << std::endl;
      std::cout << "        Tangential Tolerance = " 
	        << Get_Search_Data( ContactSearch::SEARCH_TANGENTIAL_TOLERANCE,
			       slave, master ) << std::endl;
    }
  }
  std::cout<<"    have_only_tied_interactions       = "<<have_only_tied_interactions<<std::endl;
  std::cout<<"    have_node_node_interactions       = "<<have_node_node_interactions<<std::endl;
  std::cout<<"    have_node_face_interactions       = "<<have_node_face_interactions<<std::endl;
  std::cout<<"    have_face_face_interactions       = "<<have_face_face_interactions<<std::endl;
  std::cout<<"    have_face_coverage_interactions   = "<<have_face_coverage_interactions<<std::endl;
  std::cout<<"    have_element_element_interactions = "<<have_element_element_interactions<<std::endl;
}
