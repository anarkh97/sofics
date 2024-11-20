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


#ifndef ContactTriFaceL3_h_
#define ContactTriFaceL3_h_

#include "ContactFace.h"
#include "ContactEdge.h"

class ContactFixedSizeAllocator;
template<typename DataType> class ContactNode;
template<typename DataType> class ContactEdge;


/* This class represents the linear three node triangle face with the
   following fortran numbering convention.

                              2
                              /\
                             /  \
                            /    \
                        E2 /      \ E1
                          /        \
                         /          \
                        0------------1
                              E0
*/

template<typename DataType>
class ContactTriFaceL3 : public ContactFace<DataType> {
 public:
  ContactTriFaceL3( ContactFixedSizeAllocator*,int blk_indx=-1, 
                    int indx_in_block=-1, int key=-1 );
  static ContactTriFaceL3<DataType>* new_ContactTriFaceL3(ContactFixedSizeAllocator*,
                                                int blk_indx=-1, 
                                                int indx_in_block=-1, int key=-1);
  ~ContactTriFaceL3( );
  ContactSearch::ContactEdge_Type Edge_Type() 
    {return ContactSearch::LINEEDGEL2;};
  void Get_Edge_Nodes( int, ContactNode<DataType>**);
  int Get_Edge_Number( ContactNode<DataType>** );
  int Get_Edge_Number( DataType* );

  void Compute_Normal(VariableHandle, VariableHandle );
  void Compute_Partial_Face_Normal(VariableHandle, DataType (*)[3] );
  void Compute_Second_Partial_Face_Normal(VariableHandle, DataType (*)[3] );
  void Compute_Normal(VariableHandle, DataType*, DataType* );
  void Compute_Normal(DataType**, DataType*, DataType* );
  void Compute_CharacteristicLength(VariableHandle, VariableHandle );
  void Compute_Centroid(VariableHandle, VariableHandle );
  void Compute_Edge_Normal( VariableHandle, VariableHandle,
					int , DataType*);
  void Compute_Local_Coordinates( DataType, VariableHandle, VariableHandle,
				  VariableHandle, DataType*, DataType* );
  void Compute_Local_Coordinates( VariableHandle, DataType*, DataType* );
  void Compute_Partial_Local_Coordinates_1( VariableHandle, DataType*, DataType (*)[2] );
  void Compute_Partial_Local_Coordinates_2( VariableHandle, DataType*, DataType[2], DataType[2], DataType[2] );
  void Compute_Second_Partial_Local_Coordinates_1( VariableHandle, DataType*, DataType (*)[2] );
  void Compute_Second_Partial_Local_Coordinates_2( VariableHandle, DataType*, DataType[2], DataType[2], DataType[2],
                                                   DataType[2], DataType[2], DataType[2] );
  void Compute_Second_Partial_Local_Coordinates_12(VariableHandle, DataType*, DataType (*)[2], DataType (*)[2], DataType (*)[2] );
  void Compute_Global_Coordinates( VariableHandle, DataType*, DataType* );
  void Evaluate_Shape_Functions( DataType* local_coords, DataType* shape_funcs );
  bool Is_Inside_Face( DataType* local_coords );
  inline bool IsPlanar(VariableHandle) { return true; };
  ContactFace<DataType>* Neighbor( DataType* local_coords );
  void Get_Close_Edges( DataType*, int&, int&, int& );
  void FacetDecomposition(int &, 
                          DataType*, DataType*, VariableHandle,
                          DataType*, DataType*, VariableHandle,
                          DataType*, DataType*, VariableHandle);
  void FacetStaticRestriction(int, DataType*, DataType*, DataType*, DataType*);
  void FacetDynamicRestriction(int, DataType*, DataType*);

  void Smooth_Normal( VariableHandle, VariableHandle, VariableHandle, 
		      VariableHandle, ContactSearch::Smoothing_Resolution,
		      DataType, DataType*, DataType*, DataType );
  
  void Compute_Node_Areas( VariableHandle, VariableHandle, DataType* );
                      
  int FaceEdge_Intersection(VariableHandle, ContactEdge<DataType>*, DataType*);
  
  static void Compute_Shape_Functions( DataType local_coords[3], 
                                       DataType shape_funcs[3] );
  
  static void Compute_Shape_Derivatives( DataType local_coords[3],
                                         DataType shape_derivatives[2][3] );
                                   
  void Compute_Local_Coords( DataType node_positions[MAX_NODES_PER_FACE][3],
				    DataType global_coords[3],
				    DataType local_coords[3] );
  
  static void Compute_Global_Coords( DataType node_positions[3][3],
				     DataType local_coords[3],
				     DataType global_coords[3] );

  static void Interpolate_Scalar( DataType  local_coords[3],
				  DataType  node_scalars[3],
				  DataType& interpolated_scalar );

  static void Interpolate_Vector( DataType local_coords[3],
				  DataType node_vectors[3][3],
				  DataType interpolated_vector[3] );

  using ContactTopologyEntity<DataType>::Variable;
  using ContactFace<DataType>::Node;
  using ContactFace<DataType>::Nodes_Per_Face;
  using ContactFace<DataType>::GetEdgeCurvature;
  using ContactFace<DataType>::LENGTH;

 protected:
 private:
  ContactNode<DataType>* nodes[3];
  ContactEdge<DataType>* edges[3];
  typename ContactTopologyEntity<DataType>::connection_data Node_Info[3];
  typename ContactTopologyEntity<DataType>::connection_data Edge_Info[3];
};

#endif // ifndef ContactTriFaceL3_h_
