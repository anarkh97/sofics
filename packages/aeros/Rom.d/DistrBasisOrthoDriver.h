#ifndef ROM_DISTRBASISORTHODRIVER_H
#define ROM_DISTRBASISORTHODRIVER_H

#include "DriverInterface.h"

class Communicator;
class DistrInfo;

namespace Rom {

class DistrBasisOrthoDriver : public MultiDomainDynam, public DriverInterface {
public:
  virtual void solve();
  
  DistrBasisOrthoDriver(Domain *, Communicator *);
  const DistrInfo& vectorSize() const;

private:
  Domain *domain_;
  Communicator *comm_;
};

} /* end namespace Rom */

Rom::DriverInterface *distrBasisOrthoDriverNew(Domain *);

#endif /* ROM_DISTRBASISORTHODRIVER_H */
