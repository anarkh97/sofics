#ifndef PITA_OLD_NLTIMESLICE_H
#define PITA_OLD_NLTIMESLICE_H

#include "DynamState.h"
#include "DynamStateSet.h"
#include "RankDeficientProjection.h"

#include <Math.d/Vector.h>

namespace Pita { namespace Old {

class PitaNonLinDynamic;

struct NLTimeSlice
{
  typedef GenVector<double> VecType;
  typedef DynamState<double> State;
  typedef DynamStateSet<double> StateSet;

  NLTimeSlice() : sliceRank(-1), initialTime(0.0), finalTime(0.0), converged(false), active(false), projector(0, defaultTolerance) {}
  NLTimeSlice(const PitaNonLinDynamic & probDesc, int rank);
  void initialize(const PitaNonLinDynamic & probDesc, int rank);
  void clearAllBases();

  int sliceRank;
  double initialTime;
  double finalTime;
  bool converged;
  bool active;
  
  State seedState;
  State propState;
  State jumpState;
  State nextSeedState;
  State oldSeedState;
  RankDeficientProjection projector;
  StateSet propBase;
  StateSet localBase;
  
  // Data necessary to rebuild K when used as metric for orthogonalization
  VecType locDispOG;
  double locTimeOG;

  static const double defaultTolerance;
};

} /* end namespace Old */ } /* end namespace Pita */

#endif
