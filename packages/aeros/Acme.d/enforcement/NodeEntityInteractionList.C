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



#include "NodeEntityInteractionList.h"
#include "ContactNodeFaceInteraction.h"

using namespace std;

NodeEntityInteractionList::NodeEntityInteractionList() : 
  num_nodes(0),
  list_size(0),
  list(NULL),
  num_constraints_per_node(NULL),
  constrained_node_indexes(NULL),
  current_index(-1),
  current_node(-1),
  num_face_nodes_per_interaction(NULL),
  first_face_node_index(NULL)
{}

NodeEntityInteractionList::~NodeEntityInteractionList() {
  Clear();
}

//
//  Allocate and define the fast lookup arrays for this list.  The fast lookup arrays allow quick access
//  to the scratch arrays for the various parts of this interaction.
//
void NodeEntityInteractionList::Finalize() {
  int total_face_nodes = 0;
  num_face_nodes_per_interaction = new int[list_size];
  first_face_node_index = new int[list_size];
  for(Node_Constraint_Group cnei_list = this->Node_Group_Start();
      cnei_list.Valid_Group();
      cnei_list = this->Node_Group_Next()) {
    cnei_list.Set_Constrained_Node_Index(cnei_list.Get_Interaction(0)->Node()->EnfArrayIndex());
    for(int i = 0; i < cnei_list.Num_Interactions(); ++i) {
      ContactNodeEntityInteraction *cnei = cnei_list.Get_Interaction(i);
      int cnei_index = cnei_list.Get_Index(i);
      if (cnei->Get_Type()==ContactNodeEntityInteraction::NODE_FACE_INTERACTION) {
        ContactNodeFaceInteraction *cnfi = static_cast<ContactNodeFaceInteraction*>(cnei);
	ContactFace<Real>* face = cnfi->Face();
	const int num_face_nodes = face->Nodes_Per_Face();
	num_face_nodes_per_interaction[cnei_index] = num_face_nodes;
	first_face_node_index[cnei_index] = total_face_nodes;
        Real shape_functions[MAX_NODES_PER_FACE];
        PRECONDITION( MAX_NODES_PER_FACE>=face->Nodes_Per_Face() );
        Real* coordinates = cnfi->Vector_Var(ContactNodeFaceInteraction::COORDINATES);
        face->Evaluate_Shape_Functions(coordinates, shape_functions);
	for(int inode = 0; inode < num_face_nodes; ++inode) {
	  face_node_indexes.push_back(face->Node(inode)->EnfArrayIndex());
          face_shape_functions.push_back(shape_functions[inode]);
	}
	total_face_nodes += num_face_nodes;
      } else {
	num_face_nodes_per_interaction[cnei_index] = 0;
	first_face_node_index[cnei_index] = total_face_nodes;
      }
    }
  }
}

void NodeEntityInteractionList::Clear() {
  if(list) delete [] list;
  list = NULL;
  list_size = 0;
  if(num_constraints_per_node) delete [] num_constraints_per_node;
  num_constraints_per_node = NULL;
  if(constrained_node_indexes) delete [] constrained_node_indexes;
  constrained_node_indexes = NULL;
  if(num_face_nodes_per_interaction) delete [] num_face_nodes_per_interaction;
  num_face_nodes_per_interaction = NULL;
  if(first_face_node_index) delete [] first_face_node_index;
  first_face_node_index = NULL;
  face_node_indexes.clear();
  face_shape_functions.clear();
}

void NodeEntityInteractionList::Allocate(const int list_size_, const int num_nodes_) {
  PRECONDITION(list_size_ >= 0);
  PRECONDITION(num_nodes_ >= 0);

  num_nodes = num_nodes_;
  list = new ContactNodeEntityInteraction*[list_size_];
  list_size = 0;
  num_constraints_per_node = new int[num_nodes];
  std::memset(num_constraints_per_node, 0, num_nodes * sizeof(int));
  constrained_node_indexes = new int[num_nodes];
  current_index = 0;
  current_node = 0;
}

