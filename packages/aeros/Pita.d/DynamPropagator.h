#ifndef PITA_DYNAMPROPAGATOR_H
#define PITA_DYNAMPROPAGATOR_H

#include "Fwk.h"
#include "Types.h"
#include "DynamState.h"

namespace Pita {

class DynamPropagator : public Fwk::PtrInterface<DynamPropagator> {
public:
  EXPORT_PTRINTERFACE_TYPES(DynamPropagator);

  size_t vectorSize() const { return vectorSize_; }
  const DynamState & initialState() const { return initialState_; }
  const DynamState & finalState() const { return finalState_; }
  
  virtual void initialStateIs(const DynamState & is) = 0;
  
  class Notifiee : public Fwk::BaseMultiNotifiee<const DynamPropagator> {
  public:
    EXPORT_PTRINTERFACE_TYPES(Notifiee);

    virtual void onInitialState() {}
    virtual void onFinalState()   {}

  protected:
    explicit Notifiee(const DynamPropagator * notifier) :
      Fwk::BaseMultiNotifiee<const DynamPropagator>(notifier)
    {}
  };

  // Notification interface
  virtual void lastNotifieeIs(Notifiee * notifiee) const { this->notifierDelegate().lastNotifieeIs(notifiee); }
  virtual void notifieeDel(Notifiee * notifiee) const { this->notifierDelegate().notifieeDel(notifiee); }

protected:
  explicit DynamPropagator(size_t vectorSize) :
    vectorSize_(vectorSize),
    initialState_(vectorSize_),
    finalState_(vectorSize_)
  {}

  void setVectorSize(size_t v) { vectorSize_ = v; }
  void setInitialState(DynamState is) { initialState_ = is; }
  void setFinalState(const DynamState & fs) { finalState_ = fs; }

  void initialStateNotify() { this->notifierDelegate().lastNotificationIs(&Notifiee::onInitialState); }
  void finalStateNotify() { this->notifierDelegate().lastNotificationIs(&Notifiee::onFinalState); } 

  GenNotifierDelegate<Notifiee> & notifierDelegate() const { return const_cast<DynamPropagator *>(this)->notifierDelegate_; }

private:
  size_t vectorSize_;
  DynamState initialState_;
  DynamState finalState_;
  GenNotifierDelegate<Notifiee> notifierDelegate_;

  DISALLOW_COPY_AND_ASSIGN(DynamPropagator);
};

} // end namespace Pita

#endif /* PITA_DYNAMPROPAGATOR_H */
