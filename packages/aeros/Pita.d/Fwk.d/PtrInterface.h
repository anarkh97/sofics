#ifndef FWK_PTRINTERFACE_H
#define FWK_PTRINTERFACE_H

namespace Fwk {

template <class T>
class PtrInterface {
public:
  PtrInterface() : ref_(0) {}
  unsigned long references() const { return ref_; }
  inline const PtrInterface * newRef() const { ++ref_; return this; }
  inline void deleteRef() const { if (--ref_ == 0) onZeroReferences(); }
protected:
  inline void refDec() { --ref_; }
  virtual ~PtrInterface() {}
  virtual void onZeroReferences() const { delete this; }
private:
  mutable long unsigned ref_;
};

}

#endif /* FWK_PTRINTERFACE_H */
