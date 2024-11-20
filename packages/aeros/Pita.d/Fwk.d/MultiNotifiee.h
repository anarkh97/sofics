#ifndef FWK_MULTINOTIFIEE_H
#define FWK_MULTINOTIFIEE_H

#include "PtrInterface.h"
#include "Ptr.h"
#include "Notifiee.h"

#include <set>
#include <algorithm>
#include <functional>

namespace Fwk {

template<typename Notifier, typename ActualNotifiee = typename Notifier::Notifiee>
class BaseMultiNotifiee : public RootNotifiee {
public:
  explicit BaseMultiNotifiee(Notifier * n = NULL) : notifier_(n) {
    if (n != NULL) {
      n->lastNotifieeIs(static_cast<ActualNotifiee*>(this));
    }
  }
	
  virtual ~BaseMultiNotifiee() {
    if (notifier_) {
      notifier_->notifieeDel(static_cast<ActualNotifiee*>(this));
    }
  }
	
  Notifier * notifier() const {
    return notifier_.ptr();
  }
	
  virtual void notifierIs(Notifier * n) {
    if (notifier_ != n) {
      if (notifier_) {
        notifier_->notifieeDel(static_cast<ActualNotifiee*>(this));
      }
      notifier_ = n;
      if (n != NULL) {
        n->lastNotifieeIs(static_cast<ActualNotifiee*>(this));
      }
    }
  }
	
private:
	Ptr<Notifier> notifier_;
};


template <typename Notifiee>
class GenNotifierDelegate {
public:
  typedef void (Notifiee::*NotificationType)(void);
  // NotificationType refers to a member (void -> void) function of the class Notifiee

  void lastNotifieeIs(Notifiee * notifiee) {
    notifiee_.insert(notifiee);
  }

  void notifieeDel(Notifiee * notifiee) {
    notifiee_.erase(notifiee);
  }

  void lastNotificationIs(NotificationType f) {
    std::for_each(notifiee_.begin(), notifiee_.end(), std::mem_fun(f));
  }

private:
  std::set<Notifiee *> notifiee_;
};

 
} // end namespace Fwk   

#endif /* FWK_MULTINOTIFIEE_H */
