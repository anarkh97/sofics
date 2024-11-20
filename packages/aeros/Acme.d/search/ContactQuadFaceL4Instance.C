// $Id$

#include "ContactQuadFaceL4.C"

template
ContactQuadFaceL4<Real>::ContactQuadFaceL4( ContactFixedSizeAllocator* alloc,
                                      int Block_Index, 
				      int Index_in_Block, int key );

template
ContactQuadFaceL4<Real>* ContactQuadFaceL4<Real>::new_ContactQuadFaceL4(
                        ContactFixedSizeAllocator* alloc,
                        int Block_Index, int Index_in_Block, int key);

template
void ContactQuadFaceL4_SizeAllocator<Real>(ContactFixedSizeAllocator& alloc);

template
ContactQuadFaceL4<Real>::~ContactQuadFaceL4();

template
void ContactQuadFaceL4<Real>::Compute_Normal(VariableHandle CURRENT_POSITION,
				       VariableHandle FACE_NORMAL );

template
void ContactQuadFaceL4<Real>::Compute_Partial_Face_Normal(VariableHandle CURRENT_POSITION, Real (*dface_normal)[3]);

template
void ContactQuadFaceL4<Real>::Compute_Second_Partial_Face_Normal(VariableHandle CURRENT_POSITION, Real (*dface_normal)[3] );

template
void ContactQuadFaceL4<Real>::Compute_Normal(VariableHandle POSITION,
                                       Real* normal, Real* local_coords );

template
void ContactQuadFaceL4<Real>::Compute_Normal(Real** nodal_positions,
                                       Real* local_coords, Real* normal );

template
void ContactQuadFaceL4<Real>::Compute_CharacteristicLength(VariableHandle CURRENT_POSITION,
				                     VariableHandle CHARACTERISTIC_LENGTH );

template
void ContactQuadFaceL4<Real>::Compute_Centroid( VariableHandle CURRENT_POSITION,
					  VariableHandle CENTROID );

template
void ContactQuadFaceL4<Real>::Get_Edge_Nodes( int i, ContactNode<Real>** node );

template
int ContactQuadFaceL4<Real>::Get_Edge_Number( ContactNode<Real>** edge_nodes );

template
int ContactQuadFaceL4<Real>::Get_Edge_Number( Real* local_coords );

template
void 
ContactQuadFaceL4<Real>::Compute_Edge_Normal( VariableHandle POSITION,
					VariableHandle FACE_NORMAL,
					int Edge,
					Real* edge_normal);

template
bool ContactQuadFaceL4<Real>::Is_Inside_Face( Real* local_coords );

template
ContactFace<Real>* ContactQuadFaceL4<Real>::Neighbor( Real* local_coords );

template
void ContactQuadFaceL4<Real>::Get_Close_Edges( Real* local_coords, int& number,
					 int& edge_1_id, int& edge_2_id );

template
bool ContactQuadFaceL4<Real>::IsPlanar(VariableHandle POSITION);

template
void ContactQuadFaceL4<Real>::FacetDecomposition(int& nfacets,
          Real* coordinates0, Real* normals0, VariableHandle POSITION0,
          Real* coordinates1, Real* normals1, VariableHandle POSITION1,
          Real* coordinates2, Real* normals2, VariableHandle POSITION2);

template
void ContactQuadFaceL4<Real>::FacetStaticRestriction(int nfacets, Real* coordinates, 
                                         Real* normals, Real* ctrcl_facets, 
                                         Real* ctrcl);

template
void ContactQuadFaceL4<Real>::FacetDynamicRestriction(int nfacets, 
                                                Real* ctrcl_facets, 
                                                Real* ctrcl);

template
void ContactQuadFaceL4<Real>::Smooth_Normal( VariableHandle CURRENT_POSITION,
				       VariableHandle NODE_NORMAL,
				       VariableHandle FACE_NORMAL,
				       VariableHandle CURVATURE,
			ContactSearch::Smoothing_Resolution resolution,
				       Real percentage,
				       Real* coordinates,
				       Real* smooth_normal,
                                       Real critical_curvature );

template
void ContactQuadFaceL4<Real>::Compute_Node_Areas( VariableHandle POSITION, 
                                             VariableHandle FACE_NORMAL,
                                             Real* node_areas );

template
int ContactQuadFaceL4<Real>::FaceEdge_Intersection(VariableHandle POSITION,
                                             ContactEdge<Real>* edge, Real* coords);

template
void ContactQuadFaceL4<Real>::Evaluate_Shape_Functions( Real* local_coords,
						  Real* shape_functions );

template
void ContactQuadFaceL4<Real>::Compute_Global_Coordinates( VariableHandle POSITION,
						    Real* local_coords,
						    Real* global_coords );

template
void ContactQuadFaceL4<Real>::Compute_Local_Coordinates( Real Config_Param,
					           VariableHandle POSITION0, 
					           VariableHandle POSITION1, 
						   VariableHandle FACE_NORMAL,
						   Real* global_coords,
						   Real* local_coords );

template
void ContactQuadFaceL4<Real>::Compute_Local_Coordinates( VariableHandle POSITION,
						   Real* global_coords,
						   Real* local_coords );

template
void ContactQuadFaceL4<Real>::Compute_Shape_Functions( Real* local_coords,
						 Real* shape_functions );

template
void ContactQuadFaceL4<Real>::Compute_Shape_Derivatives( Real* local_coords,
						   Real shape_derivs[2][4] );

template
void 
ContactQuadFaceL4<Real>::Compute_Local_Coords( Real node_positions[MAX_NODES_PER_FACE][3], 
					 Real global_coords[3],
					 Real local_coords[3] );

template
void 
ContactQuadFaceL4<Real>::Compute_Quad_Local_Coords( Real node_positions[MAX_NODES_PER_FACE][3], 
					 Real global_coords[3],
					 Real local_coords[3] );

template
void ContactQuadFaceL4<Real>::Compute_Global_Coords( Real node_positions[4][3],
					       Real local_coords[2],
					       Real global_coords[3] );

template
void  ContactQuadFaceL4<Real>::Interpolate_Scalar( Real  local_coords[2],
					     Real  node_scalars[4],
					     Real& interpolated_scalar );

template
void  ContactQuadFaceL4<Real>::Interpolate_Vector( Real local_coords[2],
					     Real node_vectors[4][3],
					     Real interpolated_vector[3] );

#if (MAX_FFI_DERIVATIVES > 0)
template
ContactQuadFaceL4<ActiveScalar>::ContactQuadFaceL4( ContactFixedSizeAllocator* alloc,
                                      int Block_Index, 
				      int Index_in_Block, int key );

template
ContactQuadFaceL4<ActiveScalar>* ContactQuadFaceL4<ActiveScalar>::new_ContactQuadFaceL4(
                        ContactFixedSizeAllocator* alloc,
                        int Block_Index, int Index_in_Block, int key);

template
void ContactQuadFaceL4_SizeAllocator<ActiveScalar>(ContactFixedSizeAllocator& alloc);

template
ContactQuadFaceL4<ActiveScalar>::~ContactQuadFaceL4();

template
void ContactQuadFaceL4<ActiveScalar>::Compute_Normal(VariableHandle CURRENT_POSITION,
				       VariableHandle FACE_NORMAL );

template
void ContactQuadFaceL4<ActiveScalar>::Compute_Normal(VariableHandle POSITION,
                                       ActiveScalar* normal, ActiveScalar* local_coords );

template
void ContactQuadFaceL4<ActiveScalar>::Compute_Normal(ActiveScalar** nodal_positions,
                                       ActiveScalar* local_coords, ActiveScalar* normal );

template
void ContactQuadFaceL4<ActiveScalar>::Compute_CharacteristicLength(VariableHandle CURRENT_POSITION,
				                     VariableHandle CHARACTERISTIC_LENGTH );

template
void ContactQuadFaceL4<ActiveScalar>::Compute_Centroid( VariableHandle CURRENT_POSITION,
					  VariableHandle CENTROID );

template
void ContactQuadFaceL4<ActiveScalar>::Get_Edge_Nodes( int i, ContactNode<ActiveScalar>** node );

template
int ContactQuadFaceL4<ActiveScalar>::Get_Edge_Number( ContactNode<ActiveScalar>** edge_nodes );

template
int ContactQuadFaceL4<ActiveScalar>::Get_Edge_Number( ActiveScalar* local_coords );

template
void 
ContactQuadFaceL4<ActiveScalar>::Compute_Edge_Normal( VariableHandle POSITION,
					VariableHandle FACE_NORMAL,
					int Edge,
					ActiveScalar* edge_normal);

template
bool ContactQuadFaceL4<ActiveScalar>::Is_Inside_Face( ActiveScalar* local_coords );

template
ContactFace<ActiveScalar>* ContactQuadFaceL4<ActiveScalar>::Neighbor( ActiveScalar* local_coords );

template
void ContactQuadFaceL4<ActiveScalar>::Get_Close_Edges( ActiveScalar* local_coords, int& number,
					 int& edge_1_id, int& edge_2_id );

template
bool ContactQuadFaceL4<ActiveScalar>::IsPlanar(VariableHandle POSITION);

template
void ContactQuadFaceL4<ActiveScalar>::FacetDecomposition(int& nfacets,
          ActiveScalar* coordinates0, ActiveScalar* normals0, VariableHandle POSITION0,
          ActiveScalar* coordinates1, ActiveScalar* normals1, VariableHandle POSITION1,
          ActiveScalar* coordinates2, ActiveScalar* normals2, VariableHandle POSITION2);

template
void ContactQuadFaceL4<ActiveScalar>::FacetStaticRestriction(int nfacets, ActiveScalar* coordinates, 
                                         ActiveScalar* normals, ActiveScalar* ctrcl_facets, 
                                         ActiveScalar* ctrcl);

template
void ContactQuadFaceL4<ActiveScalar>::FacetDynamicRestriction(int nfacets, 
                                                ActiveScalar* ctrcl_facets, 
                                                ActiveScalar* ctrcl);

template
void ContactQuadFaceL4<ActiveScalar>::Smooth_Normal( VariableHandle CURRENT_POSITION,
				       VariableHandle NODE_NORMAL,
				       VariableHandle FACE_NORMAL,
				       VariableHandle CURVATURE,
			ContactSearch::Smoothing_Resolution resolution,
				       ActiveScalar percentage,
				       ActiveScalar* coordinates,
				       ActiveScalar* smooth_normal,
                                       ActiveScalar critical_curvature );

template
void ContactQuadFaceL4<ActiveScalar>:: Compute_Node_Areas( VariableHandle POSITION, 
                                             VariableHandle FACE_NORMAL,
                                             ActiveScalar* node_areas );

template
int ContactQuadFaceL4<ActiveScalar>::FaceEdge_Intersection(VariableHandle POSITION,
                                             ContactEdge<ActiveScalar>* edge, ActiveScalar* coords);

template
void ContactQuadFaceL4<ActiveScalar>::Evaluate_Shape_Functions( ActiveScalar* local_coords,
						  ActiveScalar* shape_functions );

template
void ContactQuadFaceL4<ActiveScalar>::Compute_Global_Coordinates( VariableHandle POSITION,
						    ActiveScalar* local_coords,
						    ActiveScalar* global_coords );

template
void ContactQuadFaceL4<ActiveScalar>::Compute_Local_Coordinates( ActiveScalar Config_Param,
					           VariableHandle POSITION0, 
					           VariableHandle POSITION1, 
						   VariableHandle FACE_NORMAL,
						   ActiveScalar* global_coords,
						   ActiveScalar* local_coords );

template
void ContactQuadFaceL4<ActiveScalar>::Compute_Local_Coordinates( VariableHandle POSITION,
						   ActiveScalar* global_coords,
						   ActiveScalar* local_coords );

template
void ContactQuadFaceL4<ActiveScalar>::Compute_Shape_Functions( ActiveScalar* local_coords,
						 ActiveScalar* shape_functions );

template
void ContactQuadFaceL4<ActiveScalar>::Compute_Shape_Derivatives( ActiveScalar* local_coords,
						   ActiveScalar shape_derivs[2][4] );

template
void 
ContactQuadFaceL4<ActiveScalar>::Compute_Local_Coords( ActiveScalar node_positions[MAX_NODES_PER_FACE][3], 
					 ActiveScalar global_coords[3],
					 ActiveScalar local_coords[3] );

template
void 
ContactQuadFaceL4<ActiveScalar>::Compute_Quad_Local_Coords( ActiveScalar node_positions[MAX_NODES_PER_FACE][3], 
					 ActiveScalar global_coords[3],
					 ActiveScalar local_coords[3] );

template
void ContactQuadFaceL4<ActiveScalar>::Compute_Global_Coords( ActiveScalar node_positions[4][3],
					       ActiveScalar local_coords[2],
					       ActiveScalar global_coords[3] );

template
void  ContactQuadFaceL4<ActiveScalar>::Interpolate_Scalar( ActiveScalar  local_coords[2],
					     ActiveScalar  node_scalars[4],
					     ActiveScalar& interpolated_scalar );

template
void  ContactQuadFaceL4<ActiveScalar>::Interpolate_Vector( ActiveScalar local_coords[2],
					     ActiveScalar node_vectors[4][3],
					     ActiveScalar interpolated_vector[3] );
#endif
