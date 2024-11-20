#ifndef FWK_NAMEDINTERFACE_H
#define FWK_NAMEDINTERFACE_H

#include "PtrInterface.h"
#include "Types.h"

namespace Fwk {

class NamedInterface : public PtrInterface<NamedInterface> {
public:
  const String & name() const { return name_; }
protected:
	explicit NamedInterface(const String & name) : name_(name) { }
private:
  String name_;
};

} // end namespace Pita

#endif
