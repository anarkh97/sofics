#ifndef ROM_BASIS_ORTHODRIVER_H
#define ROM_BASIS_ORTHODRIVER_H

#include "DriverInterface.h"

#include "Problems.d/DynamDescr.h"
#include "VecBasisOps.h"
#include "FileNameInfo.h"
#include "VecBasisFile.h"

namespace Rom {

class BasisOrthoDriver : public SingleDomainDynamic, public DriverInterface {
public:
  virtual void solve();
  
  explicit BasisOrthoDriver(Domain *);

private:
  void preProcess();

};

} /* end namespace Rom */

Rom::DriverInterface *basisOrthoDriverNew(Domain *);

#endif /* ROM_BASIS_ORTHODRIVER_H */
