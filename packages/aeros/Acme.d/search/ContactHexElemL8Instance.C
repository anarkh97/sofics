// $Id$

#include "ContactHexElemL8.C"

template
ContactHexElemL8<Real>::ContactHexElemL8( int Block_Index, 
				    int Host_Index_in_Block, int key );

template
ContactHexElemL8<Real>* ContactHexElemL8<Real>::new_ContactHexElemL8(
                        ContactFixedSizeAllocator& alloc,
                        int Block_Index, int Host_Index_in_Block, int key);

template
void ContactHexElemL8_SizeAllocator<Real>(ContactFixedSizeAllocator& alloc);

template
ContactHexElemL8<Real>::~ContactHexElemL8();

template
void ContactHexElemL8<Real>::BuildTopology(int nID, int eID, int fID,
				     ContactFixedSizeAllocator* allocators);
template
void ContactHexElemL8<Real>::DeleteTopology(ContactFixedSizeAllocator* allocators);

template
void ContactHexElemL8<Real>::UpdateTopology(ContactFace<Real>* face, 
				      VariableHandle POSITION,
				      VariableHandle FACE_NORMAL,
				      VariableHandle NODE_NORMAL,
				      Real tol, bool use_node_normals);

template
bool
ContactHexElemL8<Real>::Is_Local_Coordinates_Inside_Element( Real* local_coords );

template
bool
ContactHexElemL8<Real>::Is_Local_Coordinates_Near_Element( Real* local_coords, Real tolerance );

template
void ContactHexElemL8<Real>::Evaluate_Shape_Functions( Real* local_coords,
						  Real* shape_functions );

template
void ContactHexElemL8<Real>::Compute_Global_Coordinates( VariableHandle POSITION,
						   Real* local_coords,
						   Real* global_coords );

template
void ContactHexElemL8<Real>::Compute_Local_Coordinates( Real Config_Param,
						  VariableHandle POSITION0, 
						  VariableHandle POSITION1, 
						  VariableHandle FACE_NORMAL,
						  Real* global_coords,
						  Real* local_coords );

template
bool ContactHexElemL8<Real>::Compute_Local_Coordinates( VariableHandle POSITION,
						  Real* global_coords,
						  Real* local_coords );

template
void ContactHexElemL8<Real>::Compute_Shape_Functions( Real* local_coords,
						Real* shape_functions );

template
void ContactHexElemL8<Real>::Compute_Shape_Derivatives( Real* local_coords,
						  Real shape_derivs[3][8] );

template
bool 
ContactHexElemL8<Real>::Compute_Local_Coords( Real node_positions[8][3], 
					Real global_coords[3],
					Real local_coords[3] );

template
void ContactHexElemL8<Real>::Compute_Global_Coords( Real node_positions[8][3],
					      Real local_coords[3],
					      Real global_coords[3] );

template
void  ContactHexElemL8<Real>::Interpolate_Scalar( Real  local_coords[3],
					    Real  node_scalars[8],
					    Real& interpolated_scalar );

template
void  ContactHexElemL8<Real>::Interpolate_Vector( Real local_coords[3],
					    Real node_vectors[8][3],
					    Real interpolated_vector[3] );

#if (MAX_FFI_DERIVATIVES > 0)
template
ContactHexElemL8<ActiveScalar>::ContactHexElemL8( int Block_Index, 
				    int Host_Index_in_Block, int key );

template
ContactHexElemL8<ActiveScalar>* ContactHexElemL8<ActiveScalar>::new_ContactHexElemL8(
                        ContactFixedSizeAllocator& alloc,
                        int Block_Index, int Host_Index_in_Block, int key);

template
void ContactHexElemL8_SizeAllocator<ActiveScalar>(ContactFixedSizeAllocator& alloc);

template
ContactHexElemL8<ActiveScalar>::~ContactHexElemL8();

template
void ContactHexElemL8<ActiveScalar>::BuildTopology(int nID, int eID, int fID,
				     ContactFixedSizeAllocator* allocators);
template
void ContactHexElemL8<ActiveScalar>::DeleteTopology(ContactFixedSizeAllocator* allocators);

template
void ContactHexElemL8<ActiveScalar>::UpdateTopology(ContactFace<ActiveScalar>* face, 
				      VariableHandle POSITION,
				      VariableHandle FACE_NORMAL,
				      VariableHandle NODE_NORMAL,
				      Real tol, bool use_node_normals);

template
bool
ContactHexElemL8<ActiveScalar>::Is_Local_Coordinates_Inside_Element( ActiveScalar* local_coords );

template
bool
ContactHexElemL8<ActiveScalar>::Is_Local_Coordinates_Near_Element( ActiveScalar* local_coords, ActiveScalar tolerance );

template
void ContactHexElemL8<ActiveScalar>::Evaluate_Shape_Functions( ActiveScalar* local_coords,
						  ActiveScalar* shape_functions );

template
void ContactHexElemL8<ActiveScalar>::Compute_Global_Coordinates( VariableHandle POSITION,
						   ActiveScalar* local_coords,
						   ActiveScalar* global_coords );

template
void ContactHexElemL8<ActiveScalar>::Compute_Local_Coordinates( ActiveScalar Config_Param,
						  VariableHandle POSITION0, 
						  VariableHandle POSITION1, 
						  VariableHandle FACE_NORMAL,
						  ActiveScalar* global_coords,
						  ActiveScalar* local_coords );

template
bool ContactHexElemL8<ActiveScalar>::Compute_Local_Coordinates( VariableHandle POSITION,
						  ActiveScalar* global_coords,
						  ActiveScalar* local_coords );

template
void ContactHexElemL8<ActiveScalar>::Compute_Shape_Functions( ActiveScalar* local_coords,
						ActiveScalar* shape_functions );

template
void ContactHexElemL8<ActiveScalar>::Compute_Shape_Derivatives( ActiveScalar* local_coords,
						  ActiveScalar shape_derivs[3][8] );

template
bool 
ContactHexElemL8<ActiveScalar>::Compute_Local_Coords( ActiveScalar node_positions[8][3], 
					ActiveScalar global_coords[3],
					ActiveScalar local_coords[3] );

template
void ContactHexElemL8<ActiveScalar>::Compute_Global_Coords( ActiveScalar node_positions[8][3],
					      ActiveScalar local_coords[3],
					      ActiveScalar global_coords[3] );

template
void  ContactHexElemL8<ActiveScalar>::Interpolate_Scalar( ActiveScalar  local_coords[3],
					    ActiveScalar  node_scalars[8],
					    ActiveScalar& interpolated_scalar );

template
void  ContactHexElemL8<ActiveScalar>::Interpolate_Vector( ActiveScalar local_coords[3],
					    ActiveScalar node_vectors[8][3],
					    ActiveScalar interpolated_vector[3] );
#endif
