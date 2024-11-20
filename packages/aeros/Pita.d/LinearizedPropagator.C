#include "LinearizedPropagator.h"
#include "NlDynamTimeIntegrator.h" 
#include "PitaNonLinDynam.h"
#include <Solvers.d/Solver.h>
#include <Math.d/SparseMatrix.h>

namespace Pita {

LinearizedPropagator::LinearizedPropagator(const NlDynamTimeIntegrator * baseIntegrator) :
  baseIntegrator_(baseIntegrator),
  temp1_(baseIntegrator->vectorSize()),
  temp2_(baseIntegrator->vectorSize())
{}

size_t
LinearizedPropagator::vectorSize() const {
  return baseIntegrator_->vectorSize();
}

const PitaNonLinDynamic *
LinearizedPropagator::probDesc() const {
  return baseIntegrator_->probDesc();
}

DynamState &
LinearizedPropagator::finalState(DynamState & state) const {
  double dt = baseIntegrator_->timeStepSize().value();
  Vector & disp = state.displacement();
  Vector & velo = state.velocity();

  // u_{n+1/2} = (M + (dt^2/4) K)^{-1} * M * Z * (u_n + (dt/2) v_n)
  temp1_.linC(disp, dt / 2.0, velo);
  probDesc()->zeroRotDofs(temp1_);
  const_cast<GenSparseMatrix<double>*>(const_cast<PitaNonLinDynamic*>(probDesc())->getMassMatrix())->mult(temp1_, temp2_);
  const_cast<PitaNonLinDynamic*>(probDesc())->getSolver()->reSolve(temp2_);

  // Half displacement increment: (1/2) * (u_{n+1} - u_{n}) = u_{n+1/2} - u_n
  temp2_ -= disp;

  // u_{n+1} = 2 * u_{n+1/2} - u_n = u_n + 2 * (u_{n+1/2} - u_n)
  disp.linAdd(2.0, temp2_);
  
  // v_{n+1} = Z * ((4/dt) * (u_{n+1/2} - u_n) - v_n)
  temp2_ *= (4.0 / dt);
  velo.diff(temp2_, velo);
  probDesc()->zeroRotDofs(velo);
  
  return state;
}

} /* end namespace Pita */
