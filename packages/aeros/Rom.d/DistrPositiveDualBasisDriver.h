#ifndef ROM_DISTRPOSITIVEDUALBASISDRIVER_H
#define ROM_DISTRPOSITIVEDUALBASISDRIVER_H

#include "DriverInterface.h"

#if defined (USE_SCALAPACK) && defined (USE_EIGEN3)
class Communicator;
class DistrInfo;

namespace Rom {

class DistrPositiveDualBasisDriver : public MultiDomainDynam, public DriverInterface {
public:
  virtual void solve();
  
  DistrPositiveDualBasisDriver(Domain *, Communicator *);

private:
  Domain *domain_;
  Communicator *comm_;
};

} /* end namespace Rom */
#endif

Rom::DriverInterface *distrPositiveDualBasisDriverNew(Domain *);

#endif /* ROM_DISTRPOSITIVEDUALBASISDRIVER_H */
