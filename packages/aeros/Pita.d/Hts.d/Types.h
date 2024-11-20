#ifndef PITA_HTS_TYPES_H
#define PITA_HTS_TYPES_H

#include "../Types.h"

namespace Pita { namespace Hts {

class HalfTimeSlice;
typedef Fwk::Numeric<HalfTimeSlice, int> HalfSliceCount;

class HalfSliceRank : public Fwk::Interval<HalfTimeSlice, int, HalfSliceCount> {
public:
  explicit HalfSliceRank(int v = 0) : Fwk::Interval<HalfTimeSlice, int, HalfSliceCount>(v) {}

  HalfSliceRank next() const { return HalfSliceRank(value() + 1); }
  HalfSliceRank previous() const { return HalfSliceRank(value() - 1); }
};

inline
HalfSliceRank
operator+(const HalfSliceRank & a, const HalfSliceCount & d) {
  return HalfSliceRank(a.value() + d.value());
}

inline 
HalfSliceRank
operator-(const HalfSliceRank & a, const HalfSliceCount & d) {
  return HalfSliceRank(a.value() - d.value());
}

enum Direction {
  NO_DIRECTION = 0,
  FORWARD,
  BACKWARD
};

enum SeedType {
  UNDEFINED_SEED = 0,
  MAIN_SEED,
  LEFT_SEED,
  RIGHT_SEED,
  SEED_JUMP,
  SEED_CORRECTION
};

inline
Fwk::OStream &
operator<<(Fwk::OStream & out, SeedType t) {
  static const char * const table = "UMLRJC?";
  out << table[t <= SEED_CORRECTION ? t : SEED_CORRECTION + 1];
  return out;
}

template <typename T>
class GenId {
public:
  T type() const { return type_; }
  HalfSliceRank rank() const { return rank_;}

  GenId(T t, HalfSliceRank r) :
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
  HalfSliceRank rank_;
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

class FullTimeSlice;
class FullSliceCount : public Fwk::Numeric<FullTimeSlice, int> {
public:
  explicit FullSliceCount(int v = 0) : Fwk::Numeric<FullTimeSlice, int>(v) {}

  operator HalfSliceCount() { return HalfSliceCount(value() * 2); }
};

} /* end namespace Hts */ } /* end namespace Pita */

#endif /* PITA_HTS_TYPES_H */
