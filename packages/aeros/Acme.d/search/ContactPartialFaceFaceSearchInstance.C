// $Id$

#include "ContactPartialFaceFaceSearch.C"

#define INSTANTIATION_HELPER(T) \
template ContactFaceFaceInteraction<T>* ContactSearch::Partial_Face_Face_Search<T>(ContactFace<T>* slave_face, \
                                                 ContactFace<T>* master_face, \
                                                 ContactElem<T>* element, \
                                                 VariableHandle POSITION, \
                                                 Real tol, \
                                                 ContactFixedSizeAllocator*);

INSTANTIATION_HELPER(Real);
