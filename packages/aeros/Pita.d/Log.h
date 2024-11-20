#ifndef PITALOG_H
#define PITALOG_H

#include "Fwk.d/Types.h"

namespace Pita {

class Log {
public:
  template <typename T> Log & operator<<(const T & t) {
    Fwk::OStream & os = outStream();
    os << t;
#ifndef NDEBUG
    os.flush();
#endif /* NDEBUG */
    return *this;
  }
  
  virtual ~Log() {}

protected:
  virtual Fwk::OStream & outStream() = 0;
};

Log & log();
  
} // end namespace Pita

#endif /* PITALOG_H */
