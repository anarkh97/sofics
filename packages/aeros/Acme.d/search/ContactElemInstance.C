// $Id$

#include "ContactElem.C"

template
ContactElem<Real>::ContactElem(ContactSearch::ContactElem_Type Type, 
			 int Block_Index, int Host_Index_in_Block, int key);

template
ContactElem<Real>::~ContactElem();

template
void ContactElem<Real>::ConnectNode( int num, ContactNode<Real>* node );

template
void ContactElem<Real>::ConnectEdge( int num, ContactEdge<Real>* edge );

template
void ContactElem<Real>::ConnectFace( int num, ContactFace<Real>* face );

template
int ContactElem<Real>::Size();

template
void ContactElem<Real>::Pack( char* buffer );

template
void ContactElem<Real>::Unpack( char* buffer );

#ifndef CONTACT_NO_MPI
template
Real ContactElem<Real>::MaxSize( VariableHandle POSITION );
#endif

#if (MAX_FFI_DERIVATIVES > 0)
template
ContactElem<ActiveScalar>::ContactElem(ContactSearch::ContactElem_Type Type,
                         int Block_Index, int Host_Index_in_Block, int key);

template
ContactElem<ActiveScalar>::~ContactElem();

template
void ContactElem<ActiveScalar>::ConnectNode( int num, ContactNode<ActiveScalar>* node );

template
void ContactElem<ActiveScalar>::ConnectEdge( int num, ContactEdge<ActiveScalar>* edge );

template
void ContactElem<ActiveScalar>::ConnectFace( int num, ContactFace<ActiveScalar>* face );

template
int ContactElem<ActiveScalar>::Size();
#endif
