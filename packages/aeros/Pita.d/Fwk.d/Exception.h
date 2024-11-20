#ifndef FWK_EXCEPTION_H
#define FWK_EXCEPTION_H

#include <cstring>

namespace Fwk {

class Exception {
public:
  const char * what() const { return what_; }

  explicit Exception(const char * message = "") {
    buildMessage("Exception ", message);
  }

  virtual ~Exception() {}

protected:
  Exception(const char * base, const char * message) {
    buildMessage(base, message);
  }

private:  
  void buildMessage(const char * base, const char * message) {
    std::strcpy(what_, base);
    std::strcat(what_, message);
  }

  char what_[256];
};


class InternalException : public Exception {
public:
  explicit InternalException(const char * message = "") :
    Exception("Internal Exception ", message)
  {}

protected:
  InternalException(const char * base, const char * message) : 
    Exception(base, message)
  {}
};


class RangeException : public Exception {
public:
  explicit RangeException(const char * message = "") :
    Exception("Range Exception ", message)
  {}

protected:
  RangeException(const char * base, const char * message) :
    Exception(base, message)
  {}
};


class NameInUseException : public Exception {
public:
  explicit NameInUseException(const char * message = "") :
    Exception("Name In Use Exception ", message)
  {}

protected:
  NameInUseException(const char * base, const char * message) :
    Exception(base, message)
  {}
};


class NetworkException : public Exception {
public:
  explicit NetworkException(const char * message = "") :
    Exception("Network Exception ", message)
  {}

protected:
  NetworkException(const char * base, const char * message) :
    Exception(base, message)
  {}
};


class MemoryException : public Exception {
public:
  explicit MemoryException(const char * message = "") :
    Exception("Memory Exception ", message)
  {}

protected:
  MemoryException(const char * base, const char * message) :
    Exception(base, message)
  {}
};

} // end namespace Fwk

#endif
