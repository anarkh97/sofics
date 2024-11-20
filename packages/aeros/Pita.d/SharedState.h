#ifndef PITA_SHAREDSTATE_H
#define PITA_SHAREDSTATE_H

#include "Fwk.h"
#include "Types.h"

namespace Pita {

/* SharedStateRoot */

class SharedStateRoot : public Fwk::NamedInterface {
public:
  EXPORT_PTRINTERFACE_TYPES(SharedStateRoot);

  enum Status {
    INACTIVE = 0,
    ACTIVE,
    CONVERGED,
    SPECIAL
  };

  Status status() const { return status_; }
  IterationRank iteration() const { return iteration_; }
  
  virtual void statusIs(Status s) = 0;
  virtual void iterationIs(IterationRank i) = 0;

protected:
  explicit SharedStateRoot(const Fwk::String & name) :
    Fwk::NamedInterface(name),
    status_(INACTIVE),
    iteration_(0)
  {}

  void setStatus(Status s) { status_ = s; }
  void setIteration(IterationRank i) { iteration_ = i; }

private:
  Status status_;
  IterationRank iteration_;
  
  DISALLOW_COPY_AND_ASSIGN(SharedStateRoot);
};

inline
OStream & operator<<(OStream & os, SharedStateRoot::Status s) {
  String out;
  switch(s) {
    case SharedStateRoot::INACTIVE:  out = "Inactive";  break;
    case SharedStateRoot::ACTIVE:    out = "Active";    break;
    case SharedStateRoot::CONVERGED: out = "Converged"; break;
    case SharedStateRoot::SPECIAL:   out = "Special";   break;
    default:                         out = "Unknown";
  }
  return os << out; 
}

/* SharedState<StateType> */

template <typename S>
class SharedState : public SharedStateRoot {
public:
  EXPORT_PTRINTERFACE_TYPES(SharedState);

  typedef S StateType;
  class Manager;
  
  const StateType & state() const { return state_; }

  virtual void stateIs(const StateType & s); 
  virtual void statusIs(Status s);
  virtual void iterationIs(IterationRank i);

  class NotifieeConst : public Fwk::BaseMultiNotifiee<const SharedState, NotifieeConst> {
  public:
    EXPORT_PTRINTERFACE_TYPES(NotifieeConst);

    virtual void onState() {}
    virtual void onStatus() {}
    virtual void onIteration() {}

  protected:
    explicit NotifieeConst(const SharedState * notifier = NULL) :
      Fwk::BaseMultiNotifiee<const SharedState, NotifieeConst>(notifier)
    {}
  };

  // Notifier implementation
  void lastNotifieeIs(NotifieeConst * notifiee) const { notifierDelegate().lastNotifieeIs(notifiee); }
  void notifieeDel(NotifieeConst * notifiee) const { notifierDelegate().notifieeDel(notifiee); }

protected:
  explicit SharedState(const Fwk::String & name);

  void setState(const StateType & s) { state_ = s; }
  
  GenNotifierDelegate<NotifieeConst> & notifierDelegate() const { return const_cast<SharedState *>(this)->notifierDelegate_; }
  
  friend class Manager;

private:
  StateType state_;
  
  // Notifier implementation
  GenNotifierDelegate<NotifieeConst> notifierDelegate_;
};

template <typename S>
class SharedState<S>::Manager : public Fwk::GenNamedInterfaceManager<SharedState<S> > {
public:
  EXPORT_PTRINTERFACE_TYPES(typename SharedState<S>::Manager);

  static typename SharedState<S>::Manager::Ptr New() {
    return new Manager();
  }

protected:
  Manager();

  virtual SharedState<S> * createNewInstance(const String & name); // Overriden
};


template <typename S>
SharedState<S>::SharedState(const Fwk::String & name) :
  SharedStateRoot(name),
  state_()
{}

template <typename S>
void
SharedState<S>::stateIs(const SharedState<S>::StateType & s) {
  setState(s);
  notifierDelegate_.lastNotificationIs(&NotifieeConst::onState);
}

template <typename S>
void
SharedState<S>::iterationIs(IterationRank s) {
  setIteration(s);
  notifierDelegate_.lastNotificationIs(&NotifieeConst::onIteration);
}

template <typename S>
void
SharedState<S>::statusIs(SharedState<S>::Status s) {
  setStatus(s);
  notifierDelegate_.lastNotificationIs(&NotifieeConst::onStatus);
}

template <typename S>
SharedState<S>::Manager::Manager() {
  // Nothing to do
}

template <typename S>
SharedState<S> *
SharedState<S>::Manager::createNewInstance(const String & name) {
  return new SharedState<S>(name); 
}

} // end namespace Pita

#endif /* PITA_SHAREDSTATE_H */
