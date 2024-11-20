#ifndef PITA_NAMEDTASK_H
#define PITA_NAMEDTASK_H

#include "Fwk.h"
#include "Types.h"

namespace Pita {

class NamedTask : public NamedInterface {
public:
  EXPORT_PTRINTERFACE_TYPES(NamedTask);

  IterationRank iteration() const { return iteration_; }
  virtual void iterationIs(IterationRank i) = 0;

protected:
  explicit NamedTask(const String & name) :
    NamedInterface(name),
    iteration_(0)
  {}

  void setIteration(IterationRank i) { iteration_ = i; }

private:
  IterationRank iteration_;
};

} /* end namespace Pita */

#endif /* PITA_NAMEDTASK_H */
