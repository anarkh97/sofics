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


#ifndef ContactLineFaceL2_h_
#define ContactLineFaceL2_h_

#include "ContactFace.h"
#include "ContactEdge.h"

template<typename DataType> class ContactNode;
template<typename DataType> class ContactEdge;
class ContactFixedSizeAllocator;

/* This class represents the line two node face with the following fortran
   numbering convention.


       1-------------2

*/

class ContactLineFaceL2 : public ContactFace<Real> {
  
 public:
  ContactLineFaceL2(ContactFixedSizeAllocator*, int blk_indx=-1, 
                    int indx_in_block=-1, int key=-1);
  static ContactLineFaceL2* new_ContactLineFaceL2(ContactFixedSizeAllocator*,
                                                  int blk_indx=-1, 
                                                  int indx_in_block=-1, 
                                                  int key=-1);
  ~ContactLineFaceL2();

  ContactSearch::ContactEdge_Type Edge_Type() 
    {return ContactSearch::NO_EDGES;};
  void Get_Edge_Nodes(int, ContactNode<Real>**) 
    {PRECONDITION(0);};
  int Get_Edge_Number(ContactNode<Real>** ) 
    {PRECONDITION(0); return -1;};
  int Get_Edge_Number( Real* )
    {PRECONDITION(0); return -1;};
  void Compute_Normal(VariableHandle,VariableHandle);
  void Compute_Normal(VariableHandle, Real*, Real* );
  void Compute_Normal(Real**, Real*, Real* );
  void Compute_CharacteristicLength(VariableHandle, VariableHandle);
  void Compute_Centroid(VariableHandle,VariableHandle);
  void Compute_Edge_Normal( VariableHandle, VariableHandle,
					int , Real*);
  void Compute_Local_Coordinates( Real, VariableHandle, VariableHandle,
			 	  VariableHandle, Real*, Real* );
  void Compute_Local_Coordinates( VariableHandle, Real*, Real* );
  void Compute_Global_Coordinates( VariableHandle, Real*, Real* );
  void Evaluate_Shape_Functions( Real* local_coords, Real* shape_funcs );
  bool Is_Inside_Face( Real* local_coords );
  inline bool IsPlanar(VariableHandle) {return true;};
  ContactFace<Real>* Neighbor( Real* local_coords );

  void Get_Close_Edges( Real*, int&, int&, int& );
  void FacetDecomposition(int &, 
                          Real*, Real*, VariableHandle,
                          Real*, Real*, VariableHandle,
                          Real*, Real*, VariableHandle);
  void FacetStaticRestriction(int, Real*, Real*, Real*, Real*);
  void FacetDynamicRestriction(int, Real*, Real*);

  void Smooth_Normal( VariableHandle,VariableHandle, VariableHandle, 
		      VariableHandle, 
		      ContactSearch::Smoothing_Resolution,
		      Real,  Real*, Real*, Real );
  int FaceEdge_Intersection(VariableHandle, ContactEdge<Real>*, Real*);
  
  void Compute_Node_Areas(VariableHandle, VariableHandle, Real*);
                      
  static void Compute_Shape_Functions( Real local_coords[3], 
                                       Real Shape_Funcs[2] );
  
  static void Compute_Shape_Derivatives( Real local_coords[3],
                                         Real shape_derivatives[2][2] );
  
  void Compute_Local_Coords( Real node_positions[MAX_NODES_PER_FACE][3],
				    Real global_coords[3],
				    Real local_coords[3] );
  
  static void Compute_Global_Coords( Real node_positions[2][3],
				     Real local_coords[3],
				     Real global_coords[3] );

  static void Interpolate_Scalar( Real  local_coords[3],
				  Real  node_scalars[2],
				  Real& interpolated_scalar );

  static void Interpolate_Vector( Real local_coords[3],
				  Real node_vectors[2][3],
				  Real interpolated_vector[3] );

 protected:
 private:
  ContactNode<Real>* nodes[2];
  ContactEdge<Real>* edges[2];
  connection_data Node_Info[2];
  connection_data Edge_Info[2];
};

#endif // ContactLineFaceL2_h_
