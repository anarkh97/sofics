#ifndef ROM_DISTRELEMENTSAMPLINGDRIVER_H
#define ROM_DISTRELEMENTSAMPLINGDRIVER_H

#include "DriverInterface.h"
#include "ParallelSparseNonNegativeLeastSquaresSolver.h"
#include <memory>
#include "BasisId.h"

class Communicator;
class DistrInfo;

namespace Rom {

class DistrElementSamplingDriver : public MultiDomainDynam, public DriverInterface {
public:
  virtual void solve();
  void computeSolution(Vector *solutions, double relativeTolerance, bool verboseFlag = true);
  
  DistrElementSamplingDriver(Domain *, Communicator *);
  const DistrInfo& vectorSize() const;

private:
  Communicator *comm_;
  ParallelSparseNonNegativeLeastSquaresSolver *solver_;

  void buildDomainCdsa();
  void subMakeMass(int isub, SparseMatrix **subM);

  void addContactElems(std::vector<int> &sampleElemIds, std::map<int, double> &weights);
};

} /* end namespace Rom */

Rom::DriverInterface *distrElementSamplingDriverNew(Domain *);

#endif /* ROM_DISTRELEMENTSAMPLINGDRIVER_H */
