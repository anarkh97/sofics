#ifndef PITA_OLD_LINEARIZED_TIME_INTEGRATOR_H
#define PITA_OLD_LINEARIZED_TIME_INTEGRATOR_H

#include "DynamState.h"

#include <Math.d/Vector.h>

namespace Pita { namespace Old {

class PitaNonLinDynamic;
template <typename Scalar> class DynamStateSet;

class LinearizedTimeIntegrator
{
public:
  typedef double Scalar;
  typedef GenVector<Scalar> VecType;
  explicit LinearizedTimeIntegrator(PitaNonLinDynamic &);
  void stepLinearizedIntegrate(DynamStateSet<double> &) const;

private:
  PitaNonLinDynamic & probDesc;
  double delta;
  double oneOverDelta;
  mutable VecType midTimeDisp;
  mutable VecType temp;
  mutable DynamState<Scalar> newState;
};

} /* end namespace Old */ } /* end namespace Pita */

#endif
