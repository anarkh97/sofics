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


#ifndef ContactWedgeElementL6_h_
#define ContactWedgeElementL6_h_

#include "ContactElement.h"

template<typename DataType> class ContactNode;
template<typename DataType> class ContactEdge;
template<typename DataType> class ContactFace;
class ContactFixedSizeAllocator;


template<typename DataType>
class ContactWedgeElemL6 : public ContactElem<DataType> {

 public:
  ContactWedgeElemL6(int blk_indx=-1, int host_indx_in_block=-1, int key=-1);
  static ContactWedgeElemL6<DataType>* new_ContactWedgeElemL6(ContactFixedSizeAllocator&,
                    int blk_indx=-1, int host_indx_in_block=-1, int key=-1);
  ~ContactWedgeElemL6();
  void BuildTopology(int, int, int, ContactFixedSizeAllocator*);
  void DeleteTopology(ContactFixedSizeAllocator*);
  void UpdateTopology(ContactFace<DataType>*, VariableHandle, VariableHandle,
                      VariableHandle, Real, bool use_node_normals=false);
  void Compute_Partial_Face_Normal( int, VariableHandle, VariableHandle, DataType (*)[3], Real, DataType (*)[3], DataType * );
  void Compute_Second_Partial_Face_Normal( int, VariableHandle, VariableHandle, DataType (*)[3], DataType (*)[3], Real,
                                           DataType (*)[3], DataType * );
  int Nodes_Per_Element() { return 6; };
  int Edges_Per_Element() { return 9; };
  int Faces_Per_Element() { return 5; };
  void Evaluate_Shape_Functions( DataType*, DataType* );
  void Compute_Local_Coordinates( DataType, VariableHandle, VariableHandle,
				  VariableHandle, DataType*, DataType* );
  bool Compute_Local_Coordinates( VariableHandle, DataType*, DataType* );
  void Compute_Global_Coordinates( VariableHandle, DataType*, DataType* );
  bool Is_Local_Coordinates_Inside_Element( DataType* );
  bool Is_Local_Coordinates_Near_Element( DataType*, DataType );
  ContactSearch::ContactNode_Type Node_Type() 
    {return ContactSearch::NODE;};
  ContactSearch::ContactEdge_Type Edge_Type() 
    {return ContactSearch::LINEEDGEL2;};
  ContactSearch::ContactFace_Type Face_Type(int i) 
    {return faces[i]->FaceType();};

  static void Compute_Shape_Functions( DataType local_coords[4], 
                                       DataType Shape_Funcs[6] );
  
  static void Compute_Shape_Derivatives( DataType local_coords[4],
                                         DataType shape_derivatives[3][6] );
  
  static bool Compute_Local_Coords( DataType node_positions[6][3],
				    DataType global_coords[3],
				    DataType local_coords[4] );
  
  static void Compute_Global_Coords( DataType node_positions[6][3],
				     DataType local_coords[4],
				     DataType global_coords[3] );

  static void Interpolate_Scalar( DataType  local_coords[4],
				  DataType  node_scalars[6],
				  DataType& interpolated_scalar );

  static void Interpolate_Vector( DataType local_coords[4],
				  DataType node_vectors[6][3],
				  DataType interpolated_vector[3] );


  inline ContactNode<DataType>** Nodes() {return nodes;};
  inline ContactEdge<DataType>** Edges() {return edges;};
  inline ContactFace<DataType>** Faces() {return faces;};
  inline int* Node_Ids() {return node_ids;};
  inline int* Edge_Ids() {return edge_ids;};
  inline int* Face_Ids() {return face_ids;};

  using ContactElem<DataType>::Node;
  using ContactElem<DataType>::Edge;
  using ContactElem<DataType>::Face;

 protected:
 private:
  ContactNode<DataType>* nodes[6];
  ContactEdge<DataType>* edges[9];
  ContactFace<DataType>* faces[5];
  int node_ids[6];
  int edge_ids[9];
  int face_ids[5];

};


#endif // ifdef ContactWedgeElementL6_h_

