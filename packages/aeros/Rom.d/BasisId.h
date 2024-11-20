#ifndef ROM_BASISID_H
#define ROM_BASISID_H

#include <string>

namespace Rom {

class BasisId {
public:
  enum Type  { STATE, RESIDUAL, JACOBIAN, FORCE, ACCELERATION, VELOCITY, INTERNALSTATE, DUALSTATE, CONSTRAINT, MUSTATE };
  enum Level { SNAPSHOTS, POD, ROB};

  Type  type()  const { return type_; }
  Level level() const { return level_; }

  BasisId(Type type, Level level) :
    type_(type), level_(level)
  {}

private:
  Type type_;
  Level level_;
};

std::string toString(BasisId::Type);
std::string toString(BasisId::Level);

} /* end namespace Rom */

#endif /* ROM_BASISID_H */
