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

#ifndef NodeEntityInteractionList_h_
#define NodeEntityInteractionList_h_

#include <Contact_Defines.h>
#include <ContactNodeEntityInteraction.h>
#include <ContactNodeFaceInteraction.h>
#include <ContactNodeSurfaceInteraction.h>
#include <vector>

class ContactSearch;
class ContactTopology;
template<typename DataType> class ContactNode;
template<typename DataType> class ContactFace;
class ContactElement;
class ContactAsymComm;
class ContactSymComm;
class ContactZoltanComm;
template<typename DataType> class ContactFaceFaceInteraction;
class ContactElementElementInteraction;
class ContactEnfModel;
class ContactParOStream;
class ScratchVariable;

class Node_Constraint_Group {
 public:
  Node_Constraint_Group(const int global_entity_list_index_,
                        const int num_cnei_,
                        ContactNodeEntityInteraction** cnei_list_,
                        int *constrained_node_index_,
                        int *num_face_nodes_per_interaction_,
                        int *first_face_node_index_,
                        std::vector<int> &face_node_indexes_,
                        std::vector<Real> &face_shape_functions_) :

    global_entity_list_index(global_entity_list_index_),
    num_cnei(num_cnei_),
    cnei_list(cnei_list_),
    constrained_node_index(constrained_node_index_),
    num_face_nodes_per_interaction(num_face_nodes_per_interaction_),
    first_face_node_index(first_face_node_index_),
    face_node_indexes(&face_node_indexes_),
    face_shape_functions(&face_shape_functions_) {}

  inline bool Valid_Group() const {
    if(cnei_list) {
      return true;
    } else {
      return false;
    }
  }

  inline int Num_Interactions() const {
    return num_cnei;
  }

  inline int Get_Num_Face_Nodes(const int i) const {
    return num_face_nodes_per_interaction[i];
  }

  inline int Get_Face_Node_Index(const int interaction_num, const int node_num) const {
    return (*face_node_indexes)[first_face_node_index[interaction_num] + node_num];
  }

  inline Real Get_Face_Shape_Function(const int interaction_num, const int node_num) const {
    return (*face_shape_functions)[first_face_node_index[interaction_num] + node_num];
  }

  inline ContactNodeEntityInteraction *Get_Interaction(const int i) const {
    return cnei_list[i];
  }

  inline void Set_Interaction(const int i, ContactNodeEntityInteraction *new_inter) const {
    cnei_list[i] = new_inter;
  }


  inline int Get_Index(const int i) const {
    return global_entity_list_index + i;
  }

  inline void Set_Constrained_Node_Index(const int new_index) {
    *constrained_node_index = new_index;
  }

  inline int Get_Constrained_Node_Index() const { return *constrained_node_index;}

 private:
  int global_entity_list_index;
  int num_cnei;
  ContactNodeEntityInteraction** cnei_list;
  int *constrained_node_index;
  int *num_face_nodes_per_interaction;
  int *first_face_node_index;
  std::vector<int> *face_node_indexes;
  std::vector<Real> *face_shape_functions;
};



//
//  Holds all of the node entity interactions
//
class NodeEntityInteractionList {
 public:
  NodeEntityInteractionList();
 
  ~NodeEntityInteractionList();
  void Clear();
  void Allocate(const int list_size_, const int num_nodes_);

  inline int Size() const { return list_size;}

  inline ContactNodeFaceInteraction *Node_Face_Iterator_Start() {
    current_index = 0;
    while(current_index < list_size) {
      ContactNodeEntityInteraction *cnei = list[current_index];
      if (cnei->Get_Type()==ContactNodeEntityInteraction::NODE_FACE_INTERACTION) {
        ContactNodeFaceInteraction *cnfi = reinterpret_cast<ContactNodeFaceInteraction*>(cnei);
        return cnfi;
      }
      current_index++;    
    }
    return NULL;
  }

  inline ContactNodeFaceInteraction *Node_Face_Iterator_Next() {
    current_index++;
    while(current_index < list_size) {
      ContactNodeEntityInteraction *cnei = list[current_index];
      if (cnei->Get_Type()==ContactNodeEntityInteraction::NODE_FACE_INTERACTION) {
        ContactNodeFaceInteraction *cnfi = reinterpret_cast<ContactNodeFaceInteraction*>(cnei);
        return cnfi;
      }
      current_index++;    
    }
    return NULL;
  }

  inline ContactNodeEntityInteraction* Iterator_Start(){
    current_index = 0;
    while(current_index < list_size) {
      if(list[current_index] != NULL) return list[current_index];
      current_index++;
    }
    return NULL;
  }

  inline ContactNodeEntityInteraction* Iterator_Next(){
    current_index++;
    while(current_index < list_size) {
      if(list[current_index] != NULL) return list[current_index];
      current_index++;
    }
    return NULL;
  }

  inline int Iterator_Index() {
    return current_index;
  }

  inline Node_Constraint_Group Node_Group_Start() {
    current_index = 0;
    current_node = 0;

    while(current_node < num_nodes) {
      if(current_index >= list_size) {
        Node_Constraint_Group calc(-1, -1, NULL, NULL, NULL, NULL, face_node_indexes, face_shape_functions);
        return calc;
      }
      //
      //  Count num valid constraints
      //
      int num_valid = 0;
      for(int i = 0; i < num_constraints_per_node[current_node]; ++i) {
        if(list[current_index + i]) num_valid++;
      }
      if(num_valid > 0) {
        Node_Constraint_Group calc(current_index, 
                                   num_valid,
                                   list + current_index,
                                   constrained_node_indexes + current_node,
                                   num_face_nodes_per_interaction + current_index,
                                   first_face_node_index + current_index,
                                   face_node_indexes,
                                   face_shape_functions);
        return calc;
      }
      current_index += num_constraints_per_node[current_node];
      current_node++;
    }
    Node_Constraint_Group calc(-1, -1, NULL, NULL, NULL, NULL, face_node_indexes, face_shape_functions);
    return calc;
  }

  inline Node_Constraint_Group Node_Group_Next() {

    while(current_node < num_nodes) {
      current_index += num_constraints_per_node[current_node];
      current_node++;
      if(current_index >= list_size) {
        Node_Constraint_Group calc(-1, -1, NULL, NULL, NULL, NULL, face_node_indexes, face_shape_functions);
        return calc;
      }

      //
      //  Count num valid constraints
      //
      int num_valid = 0;
      for(int i = 0; i < num_constraints_per_node[current_node]; ++i) {
        if(list[current_index + i]) num_valid++;
      }
      if(num_valid > 0) {
        Node_Constraint_Group calc(current_index, 
                                   num_valid,
                                   list + current_index,
                                   constrained_node_indexes + current_node,
                                   num_face_nodes_per_interaction + current_index,
                                   first_face_node_index + current_index,
                                   face_node_indexes,
                                   face_shape_functions);
        return calc;
      }
    }
    Node_Constraint_Group calc(-1, -1, NULL, NULL, NULL, NULL, face_node_indexes, face_shape_functions);
    return calc;
  }

  inline void Increment_Node() {
    current_node++;
  }

  inline int Current_Index() const {
    return current_index;
  }

  inline void Add_Constraint(ContactNodeEntityInteraction *cnei) {
    list[current_index++] = cnei;
    list_size = current_index;
    num_constraints_per_node[current_node]++;
  }

  //
  //  Allocate and define the fast lookup arrays for this list.  The fast lookup arrays allow quick access
  //  to the scratch arrays for the various parts of this interaction.
  //
  void Finalize();

 private:
  int num_nodes;
  int list_size;
  ContactNodeEntityInteraction **list;
  int *num_constraints_per_node;
  int *constrained_node_indexes;
  int current_index;
  int current_node;
  int *num_face_nodes_per_interaction;
  int *first_face_node_index;
  std::vector<int> face_node_indexes;
  std::vector<Real> face_shape_functions;
};

#endif
