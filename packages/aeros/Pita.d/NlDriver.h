#ifndef PITA_NLDRIVER_H
#define PITA_NLDRIVER_H

#include "Fwk.d/Ptr.h"
#include "Fwk.d/PtrInterface.h"
#include "Fwk.d/Macros.h"

namespace Pita {

class PitaNonLinDynamic;
  
class NlDriver : public Fwk::PtrInterface<NlDriver> {
public:
  EXPORT_PTRINTERFACE_TYPES(NlDriver);

  PitaNonLinDynamic * probDesc() const { return probDesc_; }

  // Main routine
  virtual void solve() = 0;

protected:
  explicit NlDriver(PitaNonLinDynamic * problemDescriptor) :
    probDesc_(problemDescriptor) {}

private:
  PitaNonLinDynamic * probDesc_;

  DISALLOW_COPY_AND_ASSIGN(NlDriver);
};
  
} // end namespace Pita

extern Pita::NlDriver::Ptr nlReversiblePitaDriverNew(Pita::PitaNonLinDynamic * problemDescriptor);

#endif /* PITA_NLDRIVER_H */
