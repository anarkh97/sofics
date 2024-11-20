// $Id$

#include "ContactShellTriFaceL3.C"

template
ContactShellTriFaceL3<Real>::ContactShellTriFaceL3( ContactFixedSizeAllocator* alloc,
                                    int Block_Index,
                                    int Index_in_Block,int key );

template
ContactShellTriFaceL3<Real>* ContactShellTriFaceL3<Real>::new_ContactShellTriFaceL3(
                        ContactFixedSizeAllocator* alloc,
                        int Block_Index, int Index_in_Block, int key);

template
void ContactShellTriFaceL3_SizeAllocator<Real>(ContactFixedSizeAllocator& alloc);

template
ContactShellTriFaceL3<Real>::~ContactShellTriFaceL3();

#if (MAX_FFI_DERIVATIVES > 0)
template
ContactShellTriFaceL3<ActiveScalar>::ContactShellTriFaceL3( ContactFixedSizeAllocator* alloc,
                                    int Block_Index,
                                    int Index_in_Block,int key );

template
ContactShellTriFaceL3<ActiveScalar>* ContactShellTriFaceL3<ActiveScalar>::new_ContactShellTriFaceL3(
                        ContactFixedSizeAllocator* alloc,
                        int Block_Index, int Index_in_Block, int key);

template
void ContactShellTriFaceL3_SizeAllocator<ActiveScalar>(ContactFixedSizeAllocator& alloc);

template
ContactShellTriFaceL3<ActiveScalar>::~ContactShellTriFaceL3();
#endif
