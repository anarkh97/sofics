// $Id$

#include "ContactFaceFaceSearch.C"

#define INSTANTIATION_HELPER(T) \
template ContactFaceFaceInteraction<T>* ContactSearch::Face_Face_Search<T>(ContactFace<T>* slave_face, \
                                                 ContactFace<T>* master_face, \
                                                 ContactElem<T>* element, \
                                                 VariableHandle POSITION, \
                                                 ContactFixedSizeAllocator*);

INSTANTIATION_HELPER(Real);
#if (MAX_FFI_DERIVATIVES > 0)
INSTANTIATION_HELPER(ActiveScalar);
#endif
