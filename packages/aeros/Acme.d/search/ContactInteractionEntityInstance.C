// $Id$

#include "ContactInteractionEntity.C"

template
ContactInteractionEntity<Real>::ContactInteractionEntity(Real *data_array_, ContactType base_type_);

template
ContactInteractionEntity<Real>::~ContactInteractionEntity();

template
void 
ContactInteractionEntity<Real>::Display(ContactParOStream& postream);

#if (MAX_FFI_DERIVATIVES > 0)
template
ContactInteractionEntity<ActiveScalar>::ContactInteractionEntity(ActiveScalar *data_array_, ContactType base_type_);

template
ContactInteractionEntity<ActiveScalar>::~ContactInteractionEntity();

template
void
ContactInteractionEntity<ActiveScalar>::Display(ContactParOStream& postream);
#endif
