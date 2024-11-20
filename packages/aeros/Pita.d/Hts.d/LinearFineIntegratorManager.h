#ifndef PITA_HTS_LINEARFINEINTEGRATORMANAGER_H
#define PITA_HTS_LINEARFINEINTEGRATORMANAGER_H

#include "FineIntegratorManager.h"

#include "../LinearDynamOps.h"

namespace Pita { namespace Hts {

template <typename GenAlphaIntegratorType>
class LinearFineIntegratorManager : public GenFineIntegratorManager<GenAlphaIntegratorType> {
public:
  EXPORT_PTRINTERFACE_TYPES(LinearFineIntegratorManager);

  LinearDynamOps::Manager * dynamOpsManager() const { return dynamOpsManager_.ptr(); }
  const GeneralizedAlphaParameter & parameter(Direction direction) const;

  static Ptr New(LinearDynamOps::Manager * dynOpsMgr, const GeneralizedAlphaParameter & forwardParam) {
    return new LinearFineIntegratorManager(dynOpsMgr, forwardParam);
  }

protected:
  LinearFineIntegratorManager(LinearDynamOps::Manager * dom, const GeneralizedAlphaParameter & fp);

  virtual GenAlphaIntegratorType * createFineIntegrator(Direction direction) const; // Overriden

private:
  LinearDynamOps::Manager::Ptr dynamOpsManager_;
  GeneralizedAlphaParameter forwardParameter_;
  GeneralizedAlphaParameter backwardParameter_;
};

template <typename GenAlphaIntegratorType>
LinearFineIntegratorManager<GenAlphaIntegratorType>::LinearFineIntegratorManager(LinearDynamOps::Manager * dom, const GeneralizedAlphaParameter & fp) :
  GenFineIntegratorManager<GenAlphaIntegratorType>(fp.timeStepSize()),
  dynamOpsManager_(dom),
  forwardParameter_(fp),
  backwardParameter_(GeneralizedAlphaParameter(-fp.timeStepSize(), fp.rhoInfinity()))
{}

template <typename GenAlphaIntegratorType>
const GeneralizedAlphaParameter &
LinearFineIntegratorManager<GenAlphaIntegratorType>::parameter(Direction direction) const {
  switch (direction) {
    case NO_DIRECTION: // Fall through
    case FORWARD:
      return forwardParameter_;
      break;
    case BACKWARD:
      return backwardParameter_; 
      break;
  }

  throw Fwk::InternalException("In LinearFineIntegratorManager::parameter");
}

template <typename GenAlphaIntegratorType>
GenAlphaIntegratorType *
LinearFineIntegratorManager<GenAlphaIntegratorType>::createFineIntegrator(Direction direction) const {
  if (direction == NO_DIRECTION) {
    return NULL;
  }

  return new GenAlphaIntegratorType(this->dynamOpsManager(), parameter(direction));
}


} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_HTS_LINEARFINEINTEGRATORMANAGER_H */
