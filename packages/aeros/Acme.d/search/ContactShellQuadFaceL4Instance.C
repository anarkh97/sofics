// $Id$

#include "ContactShellQuadFaceL4.C"

template
ContactShellQuadFaceL4<Real>::ContactShellQuadFaceL4( ContactFixedSizeAllocator* alloc,
                                    int Block_Index,
                                    int Index_in_Block,int key );

template
ContactShellQuadFaceL4<Real>* ContactShellQuadFaceL4<Real>::new_ContactShellQuadFaceL4(
                        ContactFixedSizeAllocator* alloc,
                        int Block_Index, int Index_in_Block, int key);

template
void ContactShellQuadFaceL4_SizeAllocator<Real>(ContactFixedSizeAllocator& alloc);

template
ContactShellQuadFaceL4<Real>::~ContactShellQuadFaceL4();

#if (MAX_FFI_DERIVATIVES > 0)
template
ContactShellQuadFaceL4<ActiveScalar>::ContactShellQuadFaceL4( ContactFixedSizeAllocator* alloc,
                                    int Block_Index,
                                    int Index_in_Block,int key );

template
ContactShellQuadFaceL4<ActiveScalar>* ContactShellQuadFaceL4<ActiveScalar>::new_ContactShellQuadFaceL4(
                        ContactFixedSizeAllocator* alloc,
                        int Block_Index, int Index_in_Block, int key);

template
void ContactShellQuadFaceL4_SizeAllocator<ActiveScalar>(ContactFixedSizeAllocator& alloc);

template
ContactShellQuadFaceL4<ActiveScalar>::~ContactShellQuadFaceL4();
#endif
