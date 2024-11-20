#ifndef PITA_HALFSLICEID_H
#define PITA_HALFSLICEID_H

#include "Fwk.h"
#include "Types.h"

namespace Pita { namespace Hts {

class HalfSliceId {
public:
  HalfSliceId(HalfSliceRank r = HalfSliceRank(-1),
              Direction d = NO_DIRECTION) :
    rank_(r),
    direction_(d)
  {}

  HalfSliceRank rank() const { return rank_; }
  Direction direction() const { return direction_; }

  bool operator==(const HalfSliceId & other) const {
    return (rank() == other.rank()) && (direction() == other.direction());
  }

  bool operator<(const HalfSliceId & other) const {
    return (rank() < other.rank()) ? true : (rank() == other.rank() && direction() < other.direction());
  }

private:
  HalfSliceRank rank_;
  Direction direction_;
};

inline
OStream &
operator<<(OStream & out, Direction d) {
  out << (d == FORWARD ? 'F' : 'B');
  return out;
}

inline
OStream &
operator<<(OStream & out, const HalfSliceId & id) {
  out << id.rank() << id.direction();
  return out; 
}

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HALFSLICEID_H */ 
