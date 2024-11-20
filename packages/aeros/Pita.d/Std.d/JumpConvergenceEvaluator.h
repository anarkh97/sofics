#ifndef PITA_STD_JUMPCONVERGENCEEVALUATOR_H
#define PITA_STD_JUMPCONVERGENCEEVALUATOR_H

#include "Fwk.h"

#include "../NamedTask.h"
#include "../Seed.h"
#include "../DynamOps.h"
#include "SliceMapping.h"

#include <Comm.d/Communicator.h>

#include "../SimpleBuffer.h"

#include <vector>
#include <map>

namespace Pita { namespace Std {

class JumpConvergenceEvaluator : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(JumpConvergenceEvaluator);

  virtual void iterationIs(IterationRank iter) = 0; // overriden

  void localJumpIs(SliceRank rank, Seed * jump);
  Seed * localJump(SliceRank rank) const;

protected:  
  JumpConvergenceEvaluator(const String & name, SliceMapping * mapping);

  SliceMapping::Ptr mapping_;

  typedef std::map<SliceRank, Seed::Ptr> JumpMap;
  JumpMap localJump_;
};

class TrivialConvergenceEvaluator : public JumpConvergenceEvaluator {
public:
  EXPORT_PTRINTERFACE_TYPES(TrivialConvergenceEvaluator);

  virtual void iterationIs(IterationRank iter); // overriden

  static Ptr New(SliceMapping * mapping) {
    return new TrivialConvergenceEvaluator(mapping);
  }

protected:
  explicit TrivialConvergenceEvaluator(SliceMapping * mapping);
};

class AccumulatedJumpConvergenceEvaluator : public JumpConvergenceEvaluator {
public:
  EXPORT_PTRINTERFACE_TYPES(AccumulatedJumpConvergenceEvaluator);

  virtual void iterationIs(IterationRank iter); // overriden

  static Ptr New(double targetRatio, const DynamOps * metric,
                 SliceMapping * mapping, Communicator * timeComm) {
      return new AccumulatedJumpConvergenceEvaluator(targetRatio, metric,
                                                     mapping, timeComm);
  }

protected: 
  AccumulatedJumpConvergenceEvaluator(double targetRatio, const DynamOps * metric,
                                      SliceMapping * mapping, Communicator * timeComm);

private:
  double targetRatio_;
  DynamOps::PtrConst metric_;
  Communicator * timeCommunicator_;

  std::vector<double> targetEstimate_;
  std::vector<double> currentEstimate_;

  SimpleBuffer<double> buffer_;
};

} /* end namespace Std */ } /* end namespace Pita */

#endif /* PITA_STD_JUMPCONVERGENCEEVALUATOR_H */
