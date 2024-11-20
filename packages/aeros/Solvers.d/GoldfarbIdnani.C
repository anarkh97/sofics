#ifdef USE_EIGEN3
#include <iostream>
#include <Solvers.d/GoldfarbIdnani.h>

template<>
void
GoldfarbIdnaniQpSolver<WrapEiSparseMat<double>,double>::solve(const double* _rhs, double* _sol)
{
  Eigen::Map<const VectorXd> rhs(_rhs,n+p+m);
  Eigen::Map<VectorXd> sol(_sol,n+p+m);

  VectorXd g0(n), ce0(p), ci0(m);
  for(int i = 0; i < neqs(); ++i) {
    switch(doftype[i]) {
      case 0 : g0(dofmap[i])  = -rhs(i); break;
      case 1 : ce0(dofmap[i]) = -rhs(i); break;
      case 2 : ci0(dofmap[i]) =  rhs(i); break;
    }
  }

  VectorXd x(n), lambda(p), mu(m);
  double traceG = diagG.sum();
  double f = Eigen::solve_quadprog2(this->solver, traceG, g0, CE, ce0, CI, ci0, x, &lambda, &mu,
                                    tol, &this->solver.permutationPinv());

  for(int i = 0; i < neqs(); ++i) {
    switch(doftype[i]) {
      case 0 : sol(i) = x(dofmap[i]); break;
      case 1 : sol(i) = -lambda(dofmap[i]); break;
      case 2 : sol(i) = mu(dofmap[i]); break;
    }
  }
}

template<>
void
GoldfarbIdnaniQpSolver<WrapEiSparseMat<complex<double> >,complex<double> >::solve(const complex<double>*, complex<double>*)
{
  std::cerr << "GoldfarbIdnaniQpSolver not implemented for complex\n";
}

#endif
