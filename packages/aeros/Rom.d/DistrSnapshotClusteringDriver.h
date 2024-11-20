#ifndef ROM_DISTRSNAPSHOTCLUSTERINGDRIVER_H
#define ROM_DISTRSNAPSHOTCLUSTERINGDRIVER_H

#include "DriverInterface.h"
#include "DistrSnapshotClusteringSolver.h"

#if defined (USE_SCALAPACK) && defined (USE_EIGEN3)
class Communicator;
class DistrInfo;

namespace Rom {

class DistrSnapshotClusteringDriver : public MultiDomainDynam, public DriverInterface {
public:
  virtual void solve();
  
  DistrSnapshotClusteringDriver(Domain *, Communicator *);

private:
  Domain *domain_;
  Communicator *comm_;

};

} /* end namespace Rom */
#endif

Rom::DriverInterface *distrSnapshotClusteringDriverNew(Domain *);

#endif /* ROM_DISTRSNAPSHOTCLUSTERINGDRIVER_H */
