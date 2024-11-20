#ifndef ROM_CONNECTIVITYUTILS_H
#define ROM_CONNECTIVITYUTILS_H

#include <Utils.d/Connectivity.h>

#include <set>
#include <algorithm>
#include <cassert>

namespace Rom {

template <typename InputIterator, typename OutputIterator>
OutputIterator
neighborhood(const Connectivity &connectivity,
             InputIterator first, InputIterator last,
             OutputIterator result) {
  Connectivity & conn_fix = const_cast<Connectivity &>(connectivity);
  std::set<int> found;

  for (InputIterator it = first; it != last; ++it) {
    const int key = *it;
    assert(key < conn_fix.csize()); 
    found.insert(key);

    auto span = conn_fix[key];
    found.insert(span.begin(), span.end());
  }

  return std::copy(found.begin(), found.end(), result);
}

template <typename InputIterator, typename OutputIterator>
OutputIterator
connections(const Connectivity &connectivity,
             InputIterator first, InputIterator last,
             OutputIterator result) {
  Connectivity & conn_fix = const_cast<Connectivity &>(connectivity);
  std::set<int> found;

  for (InputIterator it = first; it != last; ++it) {
    const int key = *it;
    assert(key < conn_fix.csize());

    auto span = conn_fix[key];
    found.insert(span.begin(), span.end());
  }

  return std::copy(found.begin(), found.end(), result);
}

} /* end namespace Rom */

#endif /* ROM_CONNECTIVITYUTILS_H */
