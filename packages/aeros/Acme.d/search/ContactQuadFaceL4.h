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


#ifndef ContactQuadFaceL4_h_
#define ContactQuadFaceL4_h_

#include "ContactFace.h"
#include "ContactEdge.h"

template<typename DataType> class ContactNode;
template<typename DataType> class ContactEdge;
class ContactFixedSizeAllocator;

/* This class represents the bilinear four node quadrilateral face with the
   following fortran numbering convention.

                             E2
                       3-------------2
                       |             |
                       |             |
                    E3 |             | E1
                       |             |
                       0-------------1
                             E0
*/

template<typename DataType>
class ContactQuadFaceL4 : public ContactFace<DataType> {

 public:
  ContactQuadFaceL4(ContactFixedSizeAllocator*, int blk_indx=-1, 
                    int indx_in_block=-1, int key=-1);
  static ContactQuadFaceL4<DataType>* new_ContactQuadFaceL4(ContactFixedSizeAllocator*,
                                                  int blk_indx=-1, 
                                                  int indx_in_block=-1, 
                                                  int key=-1);
  ~ContactQuadFaceL4();
  ContactSearch::ContactEdge_Type Edge_Type() 
    {return ContactSearch::LINEEDGEL2;};
  void Get_Edge_Nodes( int, ContactNode<DataType>** );
  int Get_Edge_Number( ContactNode<DataType>** );
  int Get_Edge_Number( DataType* );

  void Compute_Normal(VariableHandle, VariableHandle );
  void Compute_Partial_Face_Normal(VariableHandle, DataType (*)[3] );
  void Compute_Second_Partial_Face_Normal(VariableHandle, DataType (*)[3] );
  void Compute_Normal(VariableHandle, DataType*, DataType* );
  void Compute_Normal(DataType**, DataType*, DataType* );
  void Compute_CharacteristicLength(VariableHandle, VariableHandle);
  void Compute_Centroid(VariableHandle, VariableHandle );
  void Compute_Edge_Normal( VariableHandle, VariableHandle,
					int, DataType*);
  void Compute_Local_Coordinates( DataType, VariableHandle, VariableHandle,
				  VariableHandle, DataType*, DataType* );
  void Compute_Local_Coordinates( VariableHandle, DataType*, DataType* );
  void Compute_Global_Coordinates( VariableHandle, DataType*, DataType* );
  void Evaluate_Shape_Functions( DataType* Local_Coords, DataType* Shape_Funcs );
  bool Is_Inside_Face( DataType* local_coords );
  bool IsPlanar(VariableHandle);
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
                                       DataType Shape_Funcs[4] );
  
  static void Compute_Shape_Derivatives( DataType local_coords[3],
                                         DataType shape_derivatives[2][4] );
  
  void Compute_Local_Coords( DataType node_positions[MAX_NODES_PER_FACE][3],
			     DataType global_coords[3],
			     DataType local_coords[3] );
  
  static void Compute_Quad_Local_Coords(DataType node_positions[MAX_NODES_PER_FACE][3],
				    DataType global_coords[3],
				    DataType local_coords[3] );


  static void Compute_Global_Coords( DataType node_positions[4][3],
				     DataType local_coords[3],
				     DataType global_coords[3] );

  static void Interpolate_Scalar( DataType  local_coords[3],
				  DataType  node_scalars[4],
				  DataType& interpolated_scalar );

  static void Interpolate_Vector( DataType local_coords[3],
				  DataType node_vectors[4][3],
				  DataType interpolated_vector[3] );

  using ContactTopologyEntity<DataType>::Variable;
  using ContactFace<DataType>::Node;
  using ContactFace<DataType>::Nodes_Per_Face;
  using ContactFace<DataType>::GetEdgeCurvature;
  using ContactFace<DataType>::LENGTH;
  using ContactFace<DataType>::MSPARAM;
  using ContactFace<DataType>::IPENMAG;
  using ContactFace<DataType>::ILOCATION;
  using ContactFace<DataType>::ICTIMC;

 protected:
 private:
  ContactNode<DataType>* nodes[4];
  ContactEdge<DataType>* edges[4];
  typename ContactTopologyEntity<DataType>::connection_data Node_Info[4];
  typename ContactTopologyEntity<DataType>::connection_data Edge_Info[4];

};


#endif // ifdef ContactQuadFaceL4_h_
