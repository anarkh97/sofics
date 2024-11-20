// $Id$

#include "ContactEdge.C"

template
ContactEdge<Real>::ContactEdge( ContactSearch::ContactEdge_Type Type, 
                          int Block_Index, int Host_Index_in_Block,
                          ContactNode<Real> **node_list_);

template
ContactEdge<Real>::~ContactEdge();

template
void ContactEdge<Real>::ConnectFace( ContactFace<Real>* face );

template
void ContactEdge<Real>::Initialize_Lookup_Arrays();

template
void ContactEdge<Real>::Smooth_Normal( VariableHandle FACE_NORMAL, 
				 Real* coords, Real* smooth_normal );

template
void ContactEdge<Real>::Compute_Smoothed_Normal( VariableHandle FACE_NORMAL );

template
void ContactEdge<Real>::Pack( char* buffer );

template
void ContactEdge<Real>::Copy( ContactEdge* src );

#if (MAX_FFI_DERIVATIVES > 0)
template
ContactEdge<ActiveScalar>::ContactEdge( ContactSearch::ContactEdge_Type Type,
                          int Block_Index, int Host_Index_in_Block,
                          ContactNode<ActiveScalar> **node_list_);

template
ContactEdge<ActiveScalar>::~ContactEdge();

template
void ContactEdge<ActiveScalar>::ConnectFace( ContactFace<ActiveScalar>* face );

template
void ContactEdge<ActiveScalar>::Initialize_Lookup_Arrays();
#endif
