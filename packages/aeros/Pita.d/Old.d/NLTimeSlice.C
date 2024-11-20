#include "NLTimeSlice.h"
#include "PitaNonLinDynam.h"

namespace Pita { namespace Old {

const double NLTimeSlice::defaultTolerance = 1.0e-6;

NLTimeSlice::NLTimeSlice(const PitaNonLinDynamic & probDesc, int rank)
  : projector(const_cast<PitaNonLinDynamic &>(probDesc).solVecInfo(), defaultTolerance)
{
  initialize(probDesc, rank);
}

void NLTimeSlice::initialize(const PitaNonLinDynamic & probDesc, int rank)
{
  double sliceSpan = probDesc.getCoarseDt();
  int vectorSize = const_cast<PitaNonLinDynamic &>(probDesc).solVecInfo();
  
  sliceRank = rank;
  initialTime = rank * sliceSpan;
  finalTime = initialTime + sliceSpan;
  converged = false;
  active = false;

  projector.relativeToleranceIs(probDesc.getProjectionTolerance());
  projector.vectorSizeIs(vectorSize);

  propBase.reset(vectorSize, 0);
  localBase.reset(vectorSize, probDesc.getJratio() + 1);
  seedState.reset(vectorSize);
  propState.reset(vectorSize);
  jumpState.reset(vectorSize);
  nextSeedState.reset(vectorSize, 0.0);
  oldSeedState.reset(vectorSize);
}

void NLTimeSlice::clearAllBases()
{
  projector.stateDel();
  propBase.clear();
  localBase.clear();
}

} /* end namespace Old */ } /* end namespace Pita */
