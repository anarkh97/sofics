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


#ifndef ContactHexElementL8_h_
#define ContactHexElementL8_h_

#include "ContactElement.h"

template<typename DataType> class ContactNode;
template<typename DataType> class ContactEdge;
template<typename DataType> class ContactFace;
class ContactFixedSizeAllocator;


template<typename DataType>
class ContactHexElemL8 : public ContactElem<DataType> {

 public:
  ContactHexElemL8(int blk_indx=-1, int indx_in_block=-1, int key=-1);
  static ContactHexElemL8<DataType>* new_ContactHexElemL8(ContactFixedSizeAllocator&,
                    int blk_indx=-1, int indx_in_block=-1, int key=-1);
  ~ContactHexElemL8();
  void BuildTopology(int, int, int, ContactFixedSizeAllocator*);
  void DeleteTopology(ContactFixedSizeAllocator*);
  void UpdateTopology(ContactFace<DataType>*, VariableHandle, VariableHandle,
                      VariableHandle, Real, bool use_node_normals=false);
  int Nodes_Per_Element() { return 8; };
  int Edges_Per_Element() { return 12; };
  int Faces_Per_Element() { return 6; };
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
                      
  static void Compute_Shape_Functions( DataType local_coords[3], 
                                       DataType Shape_Funcs[8] );
  
  static void Compute_Shape_Derivatives( DataType local_coords[3],
                                         DataType shape_derivatives[3][8] );
  
  static bool Compute_Local_Coords( DataType node_positions[8][3],
				    DataType global_coords[3],
				    DataType local_coords[3] );
  
  static void Compute_Global_Coords( DataType node_positions[8][3],
				     DataType local_coords[3],
				     DataType global_coords[3] );

  static void Interpolate_Scalar( DataType  local_coords[3],
				  DataType  node_scalars[8],
				  DataType& interpolated_scalar );

  static void Interpolate_Vector( DataType local_coords[3],
				  DataType node_vectors[8][3],
				  DataType interpolated_vector[3] );

  inline ContactNode<DataType> **Nodes() {return nodes;};
  inline ContactEdge<DataType> **Edges() {return edges;};
  inline ContactFace<DataType> **Faces() {return faces;};
  inline int *Node_Ids() {return node_ids;};
  inline int *Edge_Ids() {return edge_ids;};
  inline int *Face_Ids() {return face_ids;};

  using ContactElem<DataType>::Node;
  using ContactElem<DataType>::Edge;
  using ContactElem<DataType>::Face;

 private:
  ContactNode<DataType>* nodes[8];
  ContactEdge<DataType>* edges[12];
  ContactFace<DataType>* faces[6];
  int node_ids[16];
  int edge_ids[24];
  int face_ids[12];

};


class ContactCartesianHexElementL8 : public ContactElement {

 public:
  ContactCartesianHexElementL8(ContactFixedSizeAllocator*, 
                               int blk_indx=-1, 
			       int indx_in_block=-1, int key=-1 );
  static ContactCartesianHexElementL8* new_ContactCartesianHexElementL8(
	ContactFixedSizeAllocator*,
        int blk_indx=-1, 
	int indx_in_block=-1, int key=-1 );
  ~ContactCartesianHexElementL8();
  ContactCartesianHexElementL8();
  
  virtual void TetDice(int &ntets, Real thex[][4][3], VariableHandle);
  virtual void Compute_Volume( VariableHandle, VariableHandle );
  virtual void Compute_Centroid( VariableHandle, VariableHandle );
  virtual int Nodes_Per_Element() { return 8; };

  void Evaluate_Shape_Functions( Real*, Real* );
  void Compute_Local_Coordinates( Real, VariableHandle, VariableHandle,
				  VariableHandle, Real*, Real* );
  void Compute_Local_Coordinates( VariableHandle, Real*, Real* );
  void Compute_Global_Coordinates( VariableHandle, Real*, Real* );
  bool Is_Local_Coordinates_Inside_Element( Real* );
  bool Is_Local_Coordinates_Near_Element( Real*, Real );
  void Interpolate_Scalar_Value( Real*, Real*, Real& );
                      
  static void Compute_Shape_Functions( Real local_coords[3], 
                                       Real Shape_Funcs[8] );
  
  static void Compute_Local_Coords( Real node_positions[8][3],
				    Real global_coords[3],
				    Real local_coords[3] );
  
  static void Compute_Global_Coords( Real node_positions[8][3],
				     Real local_coords[3],
				     Real global_coords[3] );

  static void Interpolate_Scalar( Real  local_coords[3],
				  Real  node_scalars[8],
				  Real& interpolated_scalar );

  inline ContactNode<Real>** Nodes() {return nodes;};
  inline connection_data* NodeInfo() {return Node_Info;};

 private:
  ContactNode<Real>* nodes[8];
  connection_data Node_Info[16];

};

class ContactHexElementL8 : public ContactElement {

 public:
  ContactHexElementL8(ContactFixedSizeAllocator*, 
                      int blk_indx=-1,
		      int indx_in_block=-1, int key=-1 );
  static ContactHexElementL8* new_ContactHexElementL8(
	ContactFixedSizeAllocator*, 
        int blk_indx=-1, 
	int indx_in_block=-1, int key=-1 );
  ~ContactHexElementL8();
  ContactHexElementL8();
  
  virtual void TetDice(int &ntets, Real thex[5][4][3], VariableHandle);
  virtual void Compute_Volume( VariableHandle, VariableHandle );  
  virtual void Compute_Centroid( VariableHandle, VariableHandle );  
  virtual int Nodes_Per_Element() { return 8; };

  void Evaluate_Shape_Functions( Real*, Real* );
  void Compute_Local_Coordinates( Real, VariableHandle, VariableHandle,
				  VariableHandle, Real*, Real* );
  void Compute_Local_Coordinates( VariableHandle, Real*, Real* );
  void Compute_Global_Coordinates( VariableHandle, Real*, Real* );
  bool Is_Local_Coordinates_Inside_Element( Real* );
  bool Is_Local_Coordinates_Near_Element( Real*, Real );
  void Interpolate_Scalar_Value( Real *,
				 Real *,
				 Real& );
                      
  static void Compute_Shape_Functions( Real local_coords[3], 
                                       Real Shape_Funcs[8] );
  
  static void Compute_Shape_Derivatives( Real local_coords[3],
                                         Real shape_derivatives[3][8] );
  
  static void Compute_Local_Coords( Real node_positions[8][3],
				    Real global_coords[3],
				    Real local_coords[3] );
  
  static void Compute_Global_Coords( Real node_positions[8][3],
				     Real local_coords[3],
				     Real global_coords[3] );

  static void Interpolate_Scalar( Real  local_coords[3],
				  Real  node_scalars[8],
				  Real& interpolated_scalar );

  inline ContactNode<Real>** Nodes() {return nodes;};
  inline connection_data* NodeInfo() {return Node_Info;};

 private:
  ContactNode<Real>* nodes[8];
  connection_data Node_Info[8];

};

#endif // ifdef ContactHexElementL8_h_

