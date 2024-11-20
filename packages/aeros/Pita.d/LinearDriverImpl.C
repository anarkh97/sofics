#include "LinearDriverImpl.h"

#include <Problems.d/DynamDescr.h>
#include <Driver.d/SysState.h>

namespace Pita {

LinearDriverImpl::LinearDriverImpl(SingleDomainDynamic * pbDesc,
                                   GeoSource * geoSource,
                                   Domain * domain,
                                   SolverInfo * solverInfo,
                                   Communicator * baseComm) :
  LinearDriver(pbDesc),
  geoSource_(geoSource),
  domain_(domain),
  solverInfo_(solverInfo),
  baseComm_(baseComm)
{}

DynamState
LinearDriverImpl::initialSeed() const {
  int vectorSize = probDesc()->solVecInfo();
  DynamState result = DynamState(vectorSize);
  Vector & init_d = result.displacement();
  Vector & init_v = result.velocity();
  Vector init_a(vectorSize);
  Vector init_vp(vectorSize);
  
  SysState<Vector> dummyState(init_d, init_v, init_a, init_vp); 
  probDesc()->getInitState(dummyState);

  return result;
}

} /* end namespace Pita */
