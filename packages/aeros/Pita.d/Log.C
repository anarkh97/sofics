#include "Log.h"
#include <Comm.d/Communicator.h>
#include <sstream>
#include <fstream>

extern Communicator * structCom;

namespace Pita {

class LogImpl : public Log {
public:
  LogImpl() {
    std::stringstream s;
    s << std::string("debug.") << structCom->myID();
    outStream_.open(s.str().c_str());
    outStream_.precision(12);
  }

  virtual ~LogImpl() { outStream_.close(); }

protected:
  virtual Fwk::OStream & outStream() { return outStream_; }

private:
  std::ofstream outStream_;
};

Log &
log() {
  static LogImpl log_;
  return log_;
}
  
} // end namespace Pita
