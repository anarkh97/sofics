#ifndef UTILS_VECTORTOARRAY_H
#define UTILS_VECTORTOARRAY_H

#include <vector>
#include <cstddef>

// Avoid undefined behavior when converting an empty vector to a C-style array

template <typename T>
inline
const T *
toArray(const std::vector<T> &v) {
  return !v.empty() ? &v[0] : NULL;
}

template <typename T>
inline
T *
toArray(std::vector<T> &v) {
  return !v.empty() ? &v[0] : NULL;
}

#endif /* UTILS_VECTORTOARRAY_H */
