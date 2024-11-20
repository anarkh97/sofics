#ifndef PITA_FINEINTEGRATORMANAGER_H
#define PITA_FINEINTEGRATORMANAGER_H

#include "Fwk.h"
#include "HalfSliceId.h"

#include "../DynamTimeIntegrator.h"

namespace Pita { namespace Hts {

template <typename IntegratorType>
class GenFineIntegratorManager : public Fwk::PtrInterface<GenFineIntegratorManager<IntegratorType> > {
public:
  EXPORT_PTRINTERFACE_TYPES(GenFineIntegratorManager);
  
  Seconds fineTimeStepSize() const { return fineTimeStepSize_; }
  IntegratorType * fineIntegrator(Direction direction) const;

protected:
  explicit GenFineIntegratorManager(Seconds fineTimeStep);

  virtual IntegratorType * createFineIntegrator(Direction direction) const = 0;

private:
  Seconds fineTimeStepSize_;

  mutable Fwk::Ptr<IntegratorType> forwardFineIntegrator_;
  mutable Fwk::Ptr<IntegratorType> backwardFineIntegrator_;
  
  DISALLOW_COPY_AND_ASSIGN(GenFineIntegratorManager);
};

template <typename IntegratorType>
GenFineIntegratorManager<IntegratorType>::GenFineIntegratorManager(Seconds fineTimeStepSize) :
  fineTimeStepSize_(fineTimeStepSize),
  forwardFineIntegrator_(NULL),
  backwardFineIntegrator_(NULL)
{}

template <typename IntegratorType>
IntegratorType *
GenFineIntegratorManager<IntegratorType>::fineIntegrator(Direction direction) const {
  IntegratorType * result;
  
  switch (direction) {
    case FORWARD:
      if (!forwardFineIntegrator_) {
        forwardFineIntegrator_ = createFineIntegrator(direction);
      }
      result = forwardFineIntegrator_.ptr();
      break;

    case BACKWARD:
      if (!backwardFineIntegrator_) {
        backwardFineIntegrator_ = createFineIntegrator(direction);
      }
      result = backwardFineIntegrator_.ptr();
      break;

    default:
      result = NULL;
  }

  return result;
}

} /* end namespace Hts */ } /* end namespace Hts */

#endif /* PITA_FINEINTEGRATORMANAGER_H */
