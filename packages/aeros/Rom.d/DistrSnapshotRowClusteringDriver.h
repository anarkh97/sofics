#ifndef ROM_DISTRSNAPSHOTROWCLUSTERINGDRIVER_H
#define ROM_DISTRSNAPSHOTROWCLUSTERINGDRIVER_H

#include "DriverInterface.h"
#include "DistrSnapshotRowClusteringSolver.h"

#if defined (USE_SCALAPACK) && defined (USE_EIGEN3)
class Communicator;
class DistrInfo;

namespace Rom {

class DistrSnapshotRowClusteringDriver : public MultiDomainDynam, public DriverInterface {
public:
  virtual void solve();
  
  DistrSnapshotRowClusteringDriver(Domain *, Communicator *);

private:
  Domain *domain_;
  Communicator *comm_;

};

} /* end namespace Rom */
#endif

Rom::DriverInterface *distrSnapshotRowClusteringDriverNew(Domain *);

#endif /* ROM_DISTRSNAPSHOTROWCLUSTERINGDRIVER_H */
