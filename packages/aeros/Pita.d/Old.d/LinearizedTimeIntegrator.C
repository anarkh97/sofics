#include "LinearizedTimeIntegrator.h"
#include "PitaNonLinDynam.h"
#include "DynamStateSet.h"

namespace Pita { namespace Old {

  LinearizedTimeIntegrator::LinearizedTimeIntegrator(PitaNonLinDynamic & pbDesc) :
  probDesc(pbDesc),
  delta(pbDesc.getDelta()),
  oneOverDelta(1.0 / pbDesc.getDelta()),
  midTimeDisp(pbDesc.solVecInfo()),
  temp(pbDesc.solVecInfo()),
  newState(pbDesc.solVecInfo())
{}

void
LinearizedTimeIntegrator::stepLinearizedIntegrate(DynamStateSet<double> & base) const
{
  int imax = base.numStates();
  for (int i = 0; i < imax; ++i)
  {
    probDesc.pitaTimers.newIteration();
    temp.linC(1.0, base[i].disp(), delta, base[i].vel());
    probDesc.zeroRotDofs(temp);
    const_cast<GenSparseMatrix<VecType::DataType>*>(probDesc.getMassMatrix())->mult(temp, midTimeDisp);
    probDesc.getSolver()->reSolve(midTimeDisp);
    temp.linC(-1.0, base[i].disp(), 2.0, midTimeDisp);
    newState.disp() = temp;
    temp -= base[i].disp();
    temp *= oneOverDelta;
    temp -= base[i].vel();
    probDesc.zeroRotDofs(temp);
    newState.vel() = temp;
    base[i] = newState;
  }
}

} /* end namespace Old */ } /* end namespace Pita */
