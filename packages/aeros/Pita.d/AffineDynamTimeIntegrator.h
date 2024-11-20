#ifndef PITA_AFFINEDYNAMTIMEINTEGRATOR_H
#define PITA_AFFINEDYNAMTIMEINTEGRATOR_H

#include "DynamTimeIntegrator.h"

namespace Pita {

class AffineDynamTimeIntegrator : public DynamTimeIntegrator {
public:
  EXPORT_PTRINTERFACE_TYPES(AffineDynamTimeIntegrator);

  enum ExternalForceStatus {
    NONHOMOGENEOUS = 0,
    HOMOGENEOUS
  };

  ExternalForceStatus externalForceStatus() const { return externalForceStatus_; }
  virtual void externalForceStatusIs(ExternalForceStatus efs) = 0;

protected:
  explicit AffineDynamTimeIntegrator(size_t vectorSize, ExternalForceStatus forceStatus) :
    externalForceStatus_(forceStatus),
    DynamTimeIntegrator(vectorSize)
  {}

  void setExternalForceStatus(ExternalForceStatus efs) { externalForceStatus_ = efs; }

private:
  ExternalForceStatus externalForceStatus_;
};

} /* end namespace Pita */

#endif /* PITA_AFFINEDYNAMTIMEINTEGRATOR_H */
