#ifndef PITA_NLDRIVERIMPL_H
#define PITA_NLDRIVERIMPL_H

#include "NlDriver.h"
#include "Fwk.h"
#include "DynamState.h"

namespace Pita {

class NlDriverImpl : public NlDriver {
public:
  EXPORT_PTRINTERFACE_TYPES(NlDriverImpl);

protected:
  explicit NlDriverImpl(PitaNonLinDynamic * problemDescriptor) :
    NlDriver(problemDescriptor)
  {}

  DynamState initialState() const;
};

} // end namespace Pita

#endif /* PITA_NLDRIVERIMPL_H */
