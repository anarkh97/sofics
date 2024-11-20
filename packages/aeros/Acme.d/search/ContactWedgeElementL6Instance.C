// $Id$

#include "ContactWedgeElementL6.C"

template
ContactWedgeElemL6<Real>::ContactWedgeElemL6( int Block_Index, 
				        int Host_Index_in_Block, int key ); 

template
ContactWedgeElemL6<Real>* ContactWedgeElemL6<Real>::new_ContactWedgeElemL6(
                        ContactFixedSizeAllocator& alloc,
                        int Block_Index, int Host_Index_in_Block, int key);

template
void ContactWedgeElemL6_SizeAllocator<Real>(ContactFixedSizeAllocator& alloc);

template
ContactWedgeElemL6<Real>::~ContactWedgeElemL6();

template
void ContactWedgeElemL6<Real>::BuildTopology(int nID, int eID, int fID,
				       ContactFixedSizeAllocator* allocators);

template
void ContactWedgeElemL6<Real>::DeleteTopology(ContactFixedSizeAllocator* allocators);

template
void ContactWedgeElemL6<Real>::UpdateTopology(ContactFace<Real>* face, 
					VariableHandle POSITION,
					VariableHandle FACE_NORMAL,
					VariableHandle NODE_NORMAL,
					Real tol, bool use_node_normals);

template
void ContactWedgeElemL6<Real>::Compute_Partial_Face_Normal(int i, VariableHandle CURRENT_POSITION,
                                                           VariableHandle FACE_NORMAL,
                                                           Real (*face_dface_normal)[3], Real tol,
                                                           Real (*elem_facei_dface_normal)[3], Real *dd );

template
void ContactWedgeElemL6<Real>::Compute_Second_Partial_Face_Normal( int i, VariableHandle CURRENT_POSITION,
                                                                   VariableHandle FACE_NORMAL,
                                                                   Real (*face_dface_normal)[3],
                                                                   Real (*face_d2face_normal)[3], Real tol,
                                                                   Real (*elem_d2face_normal)[3], Real *d2d );

template
bool ContactWedgeElemL6<Real>::Is_Local_Coordinates_Inside_Element( Real* local_coords );

template
bool ContactWedgeElemL6<Real>::Is_Local_Coordinates_Near_Element( Real* local_coords, Real tolerance );

template
void ContactWedgeElemL6<Real>::Evaluate_Shape_Functions( Real* local_coords,
						   Real* shape_functions );

template
void ContactWedgeElemL6<Real>::Compute_Global_Coordinates( VariableHandle POSITION,
						     Real* local_coords,
						     Real* global_coords );

template
void ContactWedgeElemL6<Real>::Compute_Local_Coordinates( Real Config_Param,
						    VariableHandle POSITION0, 
						    VariableHandle POSITION1, 
						    VariableHandle FACE_NORMAL,
						    Real* global_coords,
						    Real* local_coords );

template
bool ContactWedgeElemL6<Real>::Compute_Local_Coordinates( VariableHandle POSITION,
						    Real* global_coords,
						    Real* local_coords );

template
void ContactWedgeElemL6<Real>::Compute_Shape_Functions( Real* local_coords,
						  Real* shape_functions );

template
void ContactWedgeElemL6<Real>::Compute_Shape_Derivatives( Real* local_coords,
						    Real  shape_derivs[3][6] );

template
bool ContactWedgeElemL6<Real>::Compute_Local_Coords( Real node_positions[8][3], 
					       Real* global_coords,
					       Real* local_coords );

template
void ContactWedgeElemL6<Real>::Compute_Global_Coords( Real node_positions[6][3],
						Real local_coords[4],
						Real global_coords[3] );

template
void ContactWedgeElemL6<Real>::Interpolate_Scalar( Real  local_coords[4],
					     Real  node_scalars[6],
					     Real& interpolated_scalar );

template
void ContactWedgeElemL6<Real>::Interpolate_Vector( Real local_coords[4],
					     Real node_vectors[6][3],
					     Real interpolated_vector[3] );

#if (MAX_FFI_DERIVATIVES > 0)
template
ContactWedgeElemL6<ActiveScalar>::ContactWedgeElemL6( int Block_Index, 
				        int Host_Index_in_Block, int key ); 

template
ContactWedgeElemL6<ActiveScalar>* ContactWedgeElemL6<ActiveScalar>::new_ContactWedgeElemL6(
                        ContactFixedSizeAllocator& alloc,
                        int Block_Index, int Host_Index_in_Block, int key);

template
void ContactWedgeElemL6_SizeAllocator<ActiveScalar>(ContactFixedSizeAllocator& alloc);

template
ContactWedgeElemL6<ActiveScalar>::~ContactWedgeElemL6();

template
void ContactWedgeElemL6<ActiveScalar>::BuildTopology(int nID, int eID, int fID,
				       ContactFixedSizeAllocator* allocators);

template
void ContactWedgeElemL6<ActiveScalar>::DeleteTopology(ContactFixedSizeAllocator* allocators);

template
void ContactWedgeElemL6<ActiveScalar>::UpdateTopology(ContactFace<ActiveScalar>* face, 
					VariableHandle POSITION,
					VariableHandle FACE_NORMAL,
					VariableHandle NODE_NORMAL,
					Real tol, bool use_node_normals);

template
bool ContactWedgeElemL6<ActiveScalar>::Is_Local_Coordinates_Inside_Element( ActiveScalar* local_coords );

template
bool ContactWedgeElemL6<ActiveScalar>::Is_Local_Coordinates_Near_Element( ActiveScalar* local_coords, ActiveScalar tolerance );

template
void ContactWedgeElemL6<ActiveScalar>::Evaluate_Shape_Functions( ActiveScalar* local_coords,
						   ActiveScalar* shape_functions );

template
void ContactWedgeElemL6<ActiveScalar>::Compute_Global_Coordinates( VariableHandle POSITION,
						     ActiveScalar* local_coords,
						     ActiveScalar* global_coords );

template
void ContactWedgeElemL6<ActiveScalar>::Compute_Local_Coordinates( ActiveScalar Config_Param,
						    VariableHandle POSITION0, 
						    VariableHandle POSITION1, 
						    VariableHandle FACE_NORMAL,
						    ActiveScalar* global_coords,
						    ActiveScalar* local_coords );

template
bool ContactWedgeElemL6<ActiveScalar>::Compute_Local_Coordinates( VariableHandle POSITION,
						    ActiveScalar* global_coords,
						    ActiveScalar* local_coords );

template
void ContactWedgeElemL6<ActiveScalar>::Compute_Shape_Functions( ActiveScalar* local_coords,
						  ActiveScalar* shape_functions );

template
void ContactWedgeElemL6<ActiveScalar>::Compute_Shape_Derivatives( ActiveScalar* local_coords,
						    ActiveScalar  shape_derivs[3][6] );

template
bool ContactWedgeElemL6<ActiveScalar>::Compute_Local_Coords( ActiveScalar node_positions[8][3], 
					       ActiveScalar* global_coords,
					       ActiveScalar* local_coords );

template
void ContactWedgeElemL6<ActiveScalar>::Compute_Global_Coords( ActiveScalar node_positions[6][3],
						ActiveScalar local_coords[4],
						ActiveScalar global_coords[3] );

template
void ContactWedgeElemL6<ActiveScalar>::Interpolate_Scalar( ActiveScalar  local_coords[4],
					     ActiveScalar  node_scalars[6],
					     ActiveScalar& interpolated_scalar );

template
void ContactWedgeElemL6<ActiveScalar>::Interpolate_Vector( ActiveScalar local_coords[4],
					     ActiveScalar node_vectors[6][3],
					     ActiveScalar interpolated_vector[3] );
#endif
