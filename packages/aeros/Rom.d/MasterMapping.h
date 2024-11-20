#ifndef ROM_MASTERMAPPING_H
#define ROM_MASTERMAPPING_H

#include <vector>

namespace Rom {

struct IndexPair {
  int local, global;

  IndexPair(int l, int g) :
    local(l), global(g)
  {}
};

class MasterMapping {
public:
  template <typename IdxFwdIt, typename BoolFwdIt>
  MasterMapping(IdxFwdIt globalFirst, IdxFwdIt globalLast, BoolFwdIt masterFirst);

  typedef std::vector<IndexPair>::const_iterator IndexPairIterator;

  IndexPairIterator begin() const { return mapping_.begin(); }
  IndexPairIterator end()   const { return mapping_.end();   }

private:
  std::vector<IndexPair> mapping_;
};

template <typename IdxFwdIt, typename BoolFwdIt>
MasterMapping::MasterMapping(IdxFwdIt globalFirst, IdxFwdIt globalLast, BoolFwdIt masterFirst) {
  BoolFwdIt masterIt = masterFirst;
  int localIdx = 0;
  for (IdxFwdIt globalIt = globalFirst; globalIt != globalLast; ++globalIt) {
    if (*masterIt++) {
      mapping_.push_back(IndexPair(localIdx, *globalIt));
    }
    ++localIdx;
  }
}

} // end namespace Rom

#endif /* ROM_MASTERMAPPING_H */
