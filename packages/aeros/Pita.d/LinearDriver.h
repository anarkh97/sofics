#ifndef PITA_LINEARDRIVER_H
#define PITA_LINEARDRIVER_H

#include "Fwk.d/Ptr.h"
#include "Fwk.d/PtrInterface.h"
#include "Fwk.d/Macros.h"

class SingleDomainDynamic;

namespace Pita {

// Abstract base class to isolate main from the Pita module
class LinearDriver : public Fwk::PtrInterface<LinearDriver> {
public:
  EXPORT_PTRINTERFACE_TYPES(LinearDriver);

  SingleDomainDynamic * probDesc() const { return probDesc_; }
  
  // Main routine
  virtual void solve() = 0;

protected:
  explicit LinearDriver(SingleDomainDynamic * pbDesc) :
    probDesc_(pbDesc) {}

private:
  SingleDomainDynamic * probDesc_;

  DISALLOW_COPY_AND_ASSIGN(LinearDriver);
};

} // end namespace Pita

extern Pita::LinearDriver::Ptr linearPitaDriverNew(SingleDomainDynamic * pbDesc);
extern Pita::LinearDriver::Ptr linearReversiblePitaDriverNew(SingleDomainDynamic * pbDesc);

#endif /* PITA_LINEARDRIVER_H */
