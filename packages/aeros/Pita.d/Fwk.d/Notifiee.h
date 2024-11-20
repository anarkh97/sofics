#ifndef FWK_NOTIFIEE_H
#define FWK_NOTIFIEE_H

#include "PtrInterface.h"
#include "Ptr.h"

namespace Fwk {

class RootNotifiee : public PtrInterface<RootNotifiee> {
  /* Deliberately empty */
};

template<typename Notifier, typename ActualNotifiee = typename Notifier::Notifiee>
class BaseNotifiee : public RootNotifiee {
public:
  BaseNotifiee(Notifier* n = NULL) : notifier_(n) {
    if (n != NULL) {
      if (n->lastNotifiee()) {
        n->lastNotifiee()->notifierIs(NULL);
      }
      n->lastNotifieeIs(static_cast<ActualNotifiee*>(this));
    }
  }
	
  ~BaseNotifiee() {
    if (notifier_ != NULL) {
      notifier_->lastNotifieeIs(0);
    }
  }
	
  Ptr<Notifier> notifier() const {
    return notifier_;
  }
	
  void notifierIs(Ptr<Notifier> n) {
    if (notifier_ != n) {
      if (notifier_ != NULL) {
        notifier_->lastNotifieeIs(0);
      }
      notifier_ = n;
      if (n != NULL) {
        if (n->lastNotifiee()) {
          n->lastNotifiee()->notifierIs(NULL);
        }
        n->lastNotifieeIs(static_cast<ActualNotifiee*>(this));
      }
    }
  }
	
private:
	Ptr<Notifier> notifier_;
};
 
} //end namespace Fwk   

#endif /* FWK_NOTIFIEE_H */
