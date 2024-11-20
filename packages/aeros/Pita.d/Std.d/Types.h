#ifndef PITA_STD_TYPES_H
#define PITA_STD_TYPES_H

#include "../Types.h"

namespace Pita { namespace Std {

enum SeedType {
  UNDEFINED_SEED = 0,
  MAIN_SEED,
  PROPAGATED_SEED,
  SEED_JUMP,
  SEED_CORRECTION
};

inline
Fwk::OStream &
operator<<(Fwk::OStream & out, SeedType t) {
  static const char * const table = "UMPJC?";
  out << table[t <= SEED_CORRECTION ? t : SEED_CORRECTION + 1];
  return out;
}


template <typename T>
class GenId {
public:
  T type() const { return type_; }
  SliceRank rank() const { return rank_;}

  GenId(T t, SliceRank r) :
    type_(t), rank_(r)
  {}

  bool operator==(const GenId<T> & other) const {
    return (type() == other.type()) && (rank() == other.rank());
  }

  bool operator<(const GenId<T> & other) const {
  return (rank() < other.rank()) ? true : ((rank() == other.rank()) ? (type() < other.type()) : false);
}

private:
  T type_;
  SliceRank rank_;
};

typedef GenId<SeedType> SeedId;

inline
Fwk::OStream &
operator<<(Fwk::OStream & out, const SeedId & id) {
  out << id.type() << id.rank();
  return out;
}

class CommId {
public:
  SeedId seed() const { return seed_; }
  CpuRank cpu() const { return cpu_; }

  CommId(const SeedId & seed, CpuRank cpu) :
    seed_(seed),
    cpu_(cpu)
  {}

  bool operator==(const CommId & other) const {
    return (seed() == other.seed()) && (cpu() == other.cpu()); 
  }

private:
  SeedId seed_;
  CpuRank cpu_;
};

inline
bool
operator<(const CommId & a, const CommId & b) {
  return (a.seed() == b.seed()) ? (a.cpu() < b.cpu()) : (a.seed() < b.seed());
}

inline
Fwk::OStream &
operator<<(Fwk::OStream & out, const CommId & c) {
  out << c.seed().type() << c.seed().rank() << "->" << c.cpu();
  return out;
}

} /* end namespace Std */ } /* end namespace Pita */

#endif /* PITA_STD_TYPES_H */
