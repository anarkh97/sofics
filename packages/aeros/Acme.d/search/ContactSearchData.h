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


#ifndef ContactSearchData_h_
#define ContactSearchData_h_

#include "Contact_Defines.h"
#include "ContactSearch.h"

class ContactParOStream;
class ContactTopology;
class ContactBoundingBox;

class ContactSearchData {

 public:
  ContactSearchData(ContactTopology*);
  ~ContactSearchData();

  //Functions to set and access search data
  void Set_Search_Data( const Real* );
  inline void Set_Search_Data( ContactSearch::Search_Data_Index index,
                               int slave_key, int master_key, Real value )
    {
      // master_key == row
      // slave_key  == column
      PRECONDITION (slave_key  >= 0);
      PRECONDITION (master_key >= 0);
      PRECONDITION (slave_key  < num_search_entities);
      PRECONDITION (master_key < num_search_entities);
      PRECONDITION (index    <  size_search_data );
      int i = (master_key*num_search_entities+slave_key)*size_search_data+index;
      data_array[i] = value;
    };                
  inline Real Get_Search_Data( ContactSearch::Search_Data_Index index,
			       int slave_key, int master_key)
    {
      // master_key == row
      // slave_key  == column
      PRECONDITION (slave_key  >= 0);
      PRECONDITION (master_key >= 0);
      PRECONDITION (slave_key  < num_search_entities);
      PRECONDITION (master_key < num_search_entities);
      PRECONDITION (index    <  size_search_data );
      int i = (master_key*num_search_entities+slave_key)*size_search_data+index;
      return data_array[i];
    };
  inline ContactBoundingBox* Get_Intersection( int slave_key, int master_key)
    {
      // master_key == row
      // slave_key  == column
      PRECONDITION (slave_key  >= 0);
      PRECONDITION (master_key >= 0);
      PRECONDITION (slave_key  < num_search_entities);
      PRECONDITION (master_key < num_search_entities);
      int i = master_key*num_search_entities+slave_key;
      return &intersection[i];
    };
  inline void Set_Intersection( int slave_key, int master_key,
                                const ContactBoundingBox& box)
    {
      // master_key == row
      // slave_key  == column
      PRECONDITION (slave_key  >= 0);
      PRECONDITION (master_key >= 0);
      PRECONDITION (slave_key  < num_search_entities);
      PRECONDITION (master_key < num_search_entities);
      int i = master_key*num_search_entities+slave_key;
      intersection[i].set_box(box);
    };

  bool Master_Has_Only_Tied_Interactions(int master_key);
  bool Slave_Has_Only_Tied_Interactions(int slave_key);
  inline Real Max_Search_Tolerance()
      {return max_search_tolerance;};
  inline Real Max_Search_Normal_Tolerance()
      {return max_search_normal_tolerance;};
  inline Real Max_Search_Tangential_Tolerance()
      {return max_search_tangential_tolerance;};

  inline int Num_Search_Entities() { return num_search_entities; };

  inline Real* Search_Data()                      { return data_array; };
  inline bool Have_Tied_Interactions()            { return have_tied_interactions; };
  inline bool Have_Only_Tied_Interactions()       { return have_only_tied_interactions; };
  inline bool Have_Node_Node_Interactions()       { return have_node_node_interactions; };
  inline bool Have_Node_Face_Interactions()       { return have_node_face_interactions; };
  inline bool Have_Face_Face_Interactions()       { return have_face_face_interactions; };
  inline bool Have_Face_Coverage_Interactions()   { return have_face_coverage_interactions; };
  inline bool Have_Element_Element_Interactions() { return have_element_element_interactions; };
  inline bool SearchForInteractions()             { return search_for_interactions; };
  
  void Compute_Max_Search_Tolerances();
  void Set_Tied_Interaction_Status();
  void Set_Search_Types();
  
  void SetAllButTied();
  void SetOnlyTied();
  void SetAll();
  
  bool Is_NodeBlock_Slave(int);
  bool Is_FaceBlock_Slave(int);
  bool Is_ElemBlock_Slave(int);
  bool Is_NodeBlock_Master(int);
  bool Is_FaceBlock_Master(int);
  bool Is_ElemBlock_Master(int);
  
  bool Is_NodeBlock_NodeNodeSlave(int);
  bool Is_NodeBlock_NodeFaceSlave(int);
  bool Is_FaceBlock_NodeFaceSlave(int);
  bool Is_FaceBlock_FaceFaceSlave(int);
  bool Is_ElemBlock_ElemElemSlave(int);
  
  bool Is_NodeBlock_NodeNodeMaster(int);
  bool Is_FaceBlock_NodeFaceMaster(int);
  bool Is_FaceBlock_FaceFaceMaster(int);
  bool Is_ElemBlock_ElemElemMaster(int);
  
  inline void No_Self_Contact( )
    {
      for (int i=0; i<num_search_entities; ++i) {
        int slave_key  = i;
        int master_key = i;
        int j          = (master_key*num_search_entities+slave_key)*
                          size_search_data+ContactSearch::INTERACTION_TYPE;
        data_array[j]  = ContactSearch::NO_INTERACTION;
      }
    };                
  
  void Display(ContactParOStream&);
  void Display0();

 private:
  
  ContactSearchData(ContactSearchData&);
  ContactSearchData& operator=(ContactSearchData&);

  int num_search_entities;
  int num_elem_blocks;
  int num_face_blocks;
  int num_node_blocks;
  int num_analytic_surfaces;
  int size_search_data;
  Real* data_array;
  int*  interaction;
  ContactBoundingBox* intersection;
  Real max_search_tolerance;
  Real max_search_normal_tolerance;
  Real max_search_tangential_tolerance;
  bool have_tied_interactions;
  bool have_only_tied_interactions;
  bool have_node_node_interactions;
  bool have_node_face_interactions;
  bool have_face_face_interactions;
  bool have_face_coverage_interactions;
  bool have_element_element_interactions;
  bool search_for_interactions;
  
};

#endif // #ifndef ContactSearchData_h_
