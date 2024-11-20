#ifndef PITA_LINEARDRIVERIMPL_H
#define PITA_LINEARDRIVERIMPL_H

#include "LinearDriver.h"

#include "Fwk.h"
#include "DynamState.h"

class GeoSource;
class Domain;
class SolverInfo;
class Communicator;

namespace Pita {

class LinearDriverImpl : public LinearDriver {
public:
  EXPORT_PTRINTERFACE_TYPES(LinearDriverImpl);
  
  // Independent from global state
  GeoSource * geoSource() const { return geoSource_; }
  Domain * domain() const { return domain_; }
  SolverInfo * solverInfo() const { return solverInfo_; }
  Communicator * baseComm() const { return baseComm_; }

protected:
  LinearDriverImpl(SingleDomainDynamic *, GeoSource *, Domain *, SolverInfo *, Communicator *);

  DynamState initialSeed() const;

private:
  GeoSource * geoSource_;
  Domain * domain_;
  SolverInfo * solverInfo_;
  Communicator * baseComm_;

  DISALLOW_COPY_AND_ASSIGN(LinearDriverImpl);
};

} // end namespace Pita

#endif /* PITA_LINEARDRIVERIMPL_H */
