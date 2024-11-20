// $Id$

#include "ContactLineEdgeL2.C"

template
ContactLineEdgeL2<Real>::ContactLineEdgeL2( int Blk_Index, int Host_Index_in_Blk );

template
ContactLineEdgeL2<Real>* ContactLineEdgeL2<Real>::new_ContactLineEdgeL2(
		    ContactFixedSizeAllocator& alloc,
		    int Block_Index, int Host_Index_in_Block );

template
void ContactLineEdgeL2_SizeAllocator<Real>(ContactFixedSizeAllocator& alloc);

template
ContactLineEdgeL2<Real>::~ContactLineEdgeL2();

#if (MAX_FFI_DERIVATIVES > 0)
template
ContactLineEdgeL2<ActiveScalar>::ContactLineEdgeL2( int Blk_Index, int Host_Index_in_Blk );

template
ContactLineEdgeL2<ActiveScalar>* ContactLineEdgeL2<ActiveScalar>::new_ContactLineEdgeL2(
                    ContactFixedSizeAllocator& alloc,
                    int Block_Index, int Host_Index_in_Block );

template
void ContactLineEdgeL2_SizeAllocator<ActiveScalar>(ContactFixedSizeAllocator& alloc);

template
ContactLineEdgeL2<ActiveScalar>::~ContactLineEdgeL2();
#endif
