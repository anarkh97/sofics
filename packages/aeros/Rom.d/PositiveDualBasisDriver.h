#ifndef ROM_POSITIVE_DUAL_BASIS_DRIVER_H
#define ROM_POSITIVE_DUAL_BASIS_DRIVER_H

#include "DriverInterface.h"

#include "Problems.d/DynamDescr.h"
#include "VecBasisOps.h"
#include "FileNameInfo.h"
#include "VecBasisFile.h"

namespace Rom {

class PositiveDualBasisDriver : public SingleDomainDynamic, public DriverInterface {
public:
  virtual void solve();
  
  explicit PositiveDualBasisDriver(Domain *);

private:
  void preProcess();

};

} /* end namespace Rom */

Rom::DriverInterface *positiveDualBasisDriverNew(Domain *);

#endif /* ROM_POSITIVE_DUAL_BASIS_DRIVER_H */
