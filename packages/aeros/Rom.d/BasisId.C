#include "BasisId.h"

namespace Rom {

std::string
toString(BasisId::Type t) {
  static const std::string str[]  = { "state", "res", "jac", "for", "accel", "veloc", "internal_state", "constraint" };
  return str[t];
}

std::string
toString(BasisId::Level l) {
  static const std::string str[] = { "snap", "pod", "rob"};
  return str[l];
}

} /* end namespace Rom */
