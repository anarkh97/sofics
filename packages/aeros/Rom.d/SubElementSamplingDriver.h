#ifndef ROM_SUBELEMENTSAMPLINGDRIVER_H
#define ROM_SUBELEMENTSAMPLINGDRIVER_H

#include "ElementSamplingDriver.h"
#include <vector>

namespace Rom {

class SubElementSamplingDriver : public ElementSamplingDriver<std::vector<double>,size_t> {
public:
  explicit SubElementSamplingDriver(Domain *);
  void getGlobalWeights(Vector &solution, std::vector<double> &lweights, std::vector<int> &lelemids, bool verboseFlag = true);
  void preProcess();
};

} // end namespace Rom

Rom::DriverInterface *subElementSamplingDriverNew(Domain *);

#endif /* ROM_SUBELEMENTSAMPLINGDRIVER_H */
