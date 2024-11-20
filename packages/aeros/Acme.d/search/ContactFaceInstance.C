// $Id$

#include "ContactFace.C"

template
ContactFace<Real>::ContactFace(ContactFixedSizeAllocator* alloc,
                         ContactSearch::ContactFace_Type Type, 
			 int Block_ID, int Host_Index_in_Block, int key,
                         ContactNode<Real> **node_list_,
                         ContactEdge<Real> **edge_list_,
                         connection_data *node_info_list_,
                         connection_data *edge_info_list_);

template
ContactFace<Real>::~ContactFace();

template
void ContactFace<Real>::Initialize_Lookup_Arrays();

template
void ContactFace<Real>::SetNeighborFacesInfo();

template
void ContactFace<Real>::ConnectNode(const int num, ContactNode<Real>* node );

template
void ContactFace<Real>::ConnectEdge(const int num, ContactEdge<Real>* edge );

#ifndef CONTACT_NO_MPI
template
int ContactFace<Real>::Size_Interactions(int state);

template
void ContactFace<Real>::Pack_Interactions( char* buffer, int state );

template
void ContactFace<Real>::Unpack_Interactions( char* buffer, int state );

template
void ContactFace<Real>::Copy_Interactions( ContactFace<Real>* src, int state );

template
int ContactFace<Real>::Size_Interactions_ForSecondary(int state);

template
void ContactFace<Real>::Pack_Interactions_ForSecondary( char* buffer, int state );

template
void ContactFace<Real>::Unpack_Interactions_ForSecondary( char* buffer, int state );

template
void ContactFace<Real>::Copy_Interactions_ForSecondary( ContactFace<Real>* src, int state );

template
Real ContactFace<Real>::MaxSize( VariableHandle POSITION );
#endif

template
ContactFaceFaceInteraction<Real>* 
ContactFace<Real>::Get_FaceFace_Interaction(int interaction_number, int state );

template
void 
ContactFace<Real>::Store_FaceFace_Interaction( 
             ContactFaceFaceInteraction<Real>* ffi, int state );

template
void 
ContactFace<Real>::Delete_FaceFace_Interaction( 
             ContactFaceFaceInteraction<Real>* ffi, int state );

template
void 
ContactFace<Real>::Store_FaceCoverage_Interaction( 
             ContactFaceCoverageInteraction* fci, int state );

template
void 
ContactFace<Real>::Display_FaceFace_Interactions( ContactParOStream& postream, int state );

template
void 
ContactFace<Real>::Display_FaceCoverage_Interactions( ContactParOStream& postream, int state );

template
void
ContactFace<Real>::Update_Interactions( );

template
void
ContactFace<Real>::SetEdgeCurvature(VariableHandle CURVATURE);

template
void
ContactFace<Real>::SetEdgeCurvature(VariableHandle &CURVATURE,
                              ContactEdge<Real> *edge);

template
Real
ContactFace<Real>::GetEdgeCurvature(int i);

template
void 
ContactFace<Real>::GetEdgeInfo(ContactNode<Real>* node, 
                         ContactNode<Real>** edge_nodes, 
                         int* edge_nums);

template
void
ContactFace<Real>::SetEdgeSmoothedNormal(VariableHandle SMOOTHED_NORMAL);

template
void
ContactFace<Real>::GetEdgeSmoothedNormal(int i, Real* smoothed_normal);

template
void 
ContactFace<Real>::ComputeBoundingBoxForSearch(const int num_configs,
                                         const VariableHandle &NODE_COORD_START,
                                         const VariableHandle &NODE_COORD_END,
                                         const int  auto_tol,
                                         const Real box_inflation,
                                         const Real user_tol,
                                         ContactBoundingBox &box_c,
                                         ContactBoundingBox &box_p,
                                         ContactBoundingBox &box_s);

template
void 
ContactFace<Real>::ComputeBoundingBoxForSearch(const int num_configs,
                                         const VariableHandle &NODE_COORD_START,
                                         const VariableHandle &NODE_COORD_END,
                                         const int  auto_tol,
                                         const Real box_inflation,
					 const Real* max_node_motion,
                                         const Real max_remaining_gap_mag,
                                         const Real user_search_tol,
                                         ContactBoundingBox &box_c,
                                         ContactBoundingBox &box_p,
                                         ContactBoundingBox &box_s);

#if (MAX_FFI_DERIVATIVES > 0)
template
ContactFace<ActiveScalar>::ContactFace(ContactFixedSizeAllocator* alloc,
                         ContactSearch::ContactFace_Type Type,
                         int Block_ID, int Host_Index_in_Block, int key,
                         ContactNode<ActiveScalar> **node_list_,
                         ContactEdge<ActiveScalar> **edge_list_,
                         connection_data *node_info_list_,
                         connection_data *edge_info_list_);

template
ContactFace<ActiveScalar>::~ContactFace();

template
void ContactFace<ActiveScalar>::Initialize_Lookup_Arrays();

template
void ContactFace<ActiveScalar>::SetNeighborFacesInfo();

template
void ContactFace<ActiveScalar>::ConnectNode(const int num, ContactNode<ActiveScalar>* node );

template
void ContactFace<ActiveScalar>::ConnectEdge(const int num, ContactEdge<ActiveScalar>* edge );

template
ActiveScalar
ContactFace<ActiveScalar>::GetEdgeCurvature(int i);

template
void
ContactFace<ActiveScalar>::GetEdgeSmoothedNormal(int i, ActiveScalar* smoothed_normal);

template
ContactFaceFaceInteraction<ActiveScalar>*
ContactFace<ActiveScalar>::Get_FaceFace_Interaction(int interaction_number, int state );

template
void
ContactFace<ActiveScalar>::Store_FaceFace_Interaction(
             ContactFaceFaceInteraction<ActiveScalar>* ffi, int state );

#endif
