#ifndef FWK_TYPES_H
#define FWK_TYPES_H

#include <string>
#include <sstream>
#include <fstream>
#include <ostream>
#include <istream>

namespace Fwk {

typedef std::ostream OStream;
typedef std::istream IStream;
typedef std::ostringstream OStringStream;
typedef std::istringstream IStringStream;
typedef std::stringstream StringStream;
typedef std::string String;

template <typename T>
inline
String toString(const T & v) {
  OStringStream s;
  s << v;
  return s.str();
}

} // end namespace Fwk

#endif /* FWK_TYPES_H */
