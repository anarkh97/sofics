#ifndef ROM_PTRPTRITERADAPTERUTILS_H
#define ROM_PTRPTRITERADAPTERUTILS_H

#include <cstddef>

namespace Rom {

// Handles T** as a regular iterator
template <typename T>
class PtrPtrIterAdapter {
public:
  explicit PtrPtrIterAdapter(T **ptr) :
    ptr_(ptr)
  {}

  bool operator==(const PtrPtrIterAdapter &other) const {
    return this->ptr_ == other.ptr_;
  }
  
  bool operator!=(const PtrPtrIterAdapter &other) const {
    return this->ptr_ != other.ptr_;
  }

  const T &operator*() const { return **ptr_; }
  T &operator*() { return **ptr_; }
 
  const T *operator->() const { return *ptr_; }
  T *operator->() { return *ptr_; }
  
  PtrPtrIterAdapter &operator++() { ++ptr_; return *this; }
  PtrPtrIterAdapter operator++(int) { PtrPtrIterAdapter temp(*this); ++ptr_; return temp; }

  const T &operator[](std::size_t i) const { return *(ptr_[i]); }
  T &operator[](std::size_t i) { return *(ptr_[i]); }

private:
  T **ptr_;
};

} /* end namespace Rom */

#endif /* ROM_PTRPTRITERADAPTERUTILS_H */
