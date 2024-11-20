// $Id$

#include "ContactTopologyEntity.C"

template
ContactTopologyEntity<Real>::ContactTopologyEntity(Real *data_array_, 
                                             const ContactType base_type_);

template
ContactTopologyEntity<Real>::ContactTopologyEntity( int Block_ID, 
			                      int host_index_in_block,
                                              Real *data_array_, 
                                              const ContactType base_type_);

template
ContactTopologyEntity<Real>::~ContactTopologyEntity();

template
void 
ContactTopologyEntity<Real>::Display(ContactParOStream& postream);

#if (MAX_FFI_DERIVATIVES > 0)
template
ContactTopologyEntity<ActiveScalar>::ContactTopologyEntity( int Block_ID,
                                              int host_index_in_block,
                                              ActiveScalar *data_array_,
                                              const ContactType base_type_);

template
ContactTopologyEntity<ActiveScalar>::~ContactTopologyEntity();
#endif

