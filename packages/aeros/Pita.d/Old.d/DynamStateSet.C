#include "DynamStateSet.h"

#include <Math.d/SparseMatrix.h>
#include <Solvers.d/Solver.h>

#include <cstdio>

using std::fprintf;

namespace Pita { namespace Old {

// Memory management

template <typename Scalar>
void DynamStateSet<Scalar>::init_(int vectorSize, int maxNumStates)
{
  stateSet_.reserve((maxNumStates > 0) ? maxNumStates : 0);
  vectorSize_ = (vectorSize > 0) ? vectorSize : 0;
}

template <typename Scalar>
void DynamStateSet<Scalar>::merge(const DynamStateSet<Scalar> &dss)
{
  if (dss.vectorSize_ != vectorSize_)
  {
    fprintf(stderr, "Warning -- in DynamStateSet::merge : vectorSize mismatch -- function aborted\n");
    return;
  }

  stateSet_.reserve(stateSet_.size() + dss.stateSet_.size());
  typename std::vector<State>::iterator it;
  for (it = dss.stateSet_.begin(); it != dss.stateSet_.end(); ++it)
  {
    stateSet_.push_back(*it);
  }
}

// Get raw data to contiguous buffer
template <typename Scalar>
void DynamStateSet<Scalar>::getRaw(Scalar * buffer) const
{
  int shift = 2 * vectorSize_;
  int numStates = stateSet_.size();
  for (int i = 0; i < numStates; ++i)
  {
    stateSet_[i].getRaw(buffer + i * shift);
  }
}

} /* end namespace Old */ } /* end namespace Pita */
