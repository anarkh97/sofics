#ifndef PITA_CORRECTIONPROPAGATOR_H
#define PITA_CORRECTIONPROPAGATOR_H

#include "Fwk.h"
#include "NamedTask.h"

namespace Pita {

template <typename S> class SharedState;

template <typename IS, typename FS = IS>
class CorrectionPropagator : public NamedTask {
public:
  EXPORT_PTRINTERFACE_TYPES(CorrectionPropagator);
  typedef Fwk::GenManagerInterface<CorrectionPropagator<IS, FS> *, String> Manager;
  
  // Inputs  
  const SharedState<IS> * jump() const { return jump_.ptr(); }
  const SharedState<IS> * correction() const { return correction_.ptr(); }
 
  virtual void jumpIs(const SharedState<IS> * as) { setJump(as); }
  virtual void correctionIs(const SharedState<IS> * sc) { setCorrection(sc); }

  // Outputs
  SharedState<FS> * nextCorrection() const { return nextCorrection_.ptr(); }
 
  virtual void nextCorrectionIs(SharedState<FS> * nc) { setNextCorrection(nc); }
  
protected:
  explicit CorrectionPropagator(const String & name) :
    NamedTask(name)
  {}
  
  void setCorrection(const SharedState<IS> * c) { correction_ = c; }
  void setJump(const SharedState<IS> * j) { jump_ = j; }
  void setNextCorrection(SharedState<FS> * nc) { nextCorrection_ = nc; }
  
private:
  Fwk::Ptr<const SharedState<IS> > jump_;
  Fwk::Ptr<const SharedState<IS> > correction_;
  Fwk::Ptr<SharedState<FS> > nextCorrection_;

  DISALLOW_COPY_AND_ASSIGN(CorrectionPropagator);
};
  
} /* end namespace Pita */

#endif /* PITA_CORRECTIONPROPAGATOR_H */
