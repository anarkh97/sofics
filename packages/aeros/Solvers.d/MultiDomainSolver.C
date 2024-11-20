#include <complex>
#include <Driver.d/Communicator.h>
#include <Driver.d/SubDomain.h>
#include <Solvers.d/MultiDomainSolver.h>

template<class Scalar>
void
MultiDomainSolver<Scalar>::reSolve(GenDistrVector<Scalar> &solution)
{
  GenDistrVector<Scalar> rhs(solution);
  timespec tS1,tS2;
#ifdef MDS_TIMING
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &tS1);
#endif
  //------------------------------------------
  solve(rhs, solution);
  //------------------------------------------
#ifdef MDS_TIMING
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &tS2);
  std::cout << "              **** total solve time = " << (tS2.tv_nsec-tS1.tv_nsec)/1e9 << std::endl;
#endif
}

template<class Scalar>
void
MultiDomainSolver<Scalar>::solve(const GenDistrVector<Scalar> &rhs, GenDistrVector<Scalar> &solution)
{
  Scalar *f_g = new Scalar[neq];
  Scalar *u_g = new Scalar[neq];

  // 1. f_g = L^T*f
  multLT(rhs, f_g);

  // 2. u = arg min f(v) = 1/2 v^T K v - f^T v subj to B*u = 0 (for example)
  solve(f_g, u_g);

  // 3. u = L*u_g
  multL(u_g, solution);

  delete [] f_g;
  delete [] u_g;
}

template<class Scalar>
double
MultiDomainSolver<Scalar>::getFNormSq(GenDistrVector<Scalar> &f)
{
  Scalar *f_g = new Scalar[neq];
  multLT(f, f_g);
  Scalar ret = 0;
  for(int i = 0; i < neq; ++i) ret += f_g[i]*f_g[i];
  delete [] f_g;
  return ScalarTypes::norm(ret);
}

template<class Scalar>
void
MultiDomainSolver<Scalar>::multLTinv(Scalar *f_g, GenDistrVector<Scalar> &f) const
{
  // f = L^{-T}*f_g = L*(L^T*L)^{-1}*f_g
  for(int i = 0; i < nsub; ++i) sd[i]->multLTinv(f_g, f.subData(i));
}

template<class Scalar>
void
MultiDomainSolver<Scalar>::multL(Scalar *u_g, GenDistrVector<Scalar> &u) const
{
  // u = L*u_g
  for(int i = 0; i < nsub; ++i) sd[i]->multL(u_g, u.subData(i));
}

template<class Scalar>
void
MultiDomainSolver<Scalar>::multLinv(GenDistrVector<Scalar> &u, Scalar *u_g) const
{
  // u_g = L^{-1}*u = (L^T*L)^{-1}*L^T*u
  for(int i = 0; i < neq; ++i) u_g[i] = 0;
  for(int i = 0; i < nsub; ++i) sd[i]->multAddLinv(u.subData(i), u_g);
#ifdef DISTRIBUTED
  if(com) com->globalSum(neq, u_g);
#endif
}

template<class Scalar>
void
MultiDomainSolver<Scalar>::multLT(const GenDistrVector<Scalar> &f, Scalar *f_g) const
{
  // f_g = L^T*f
  for(int i = 0; i < neq; ++i) f_g[i] = 0;
  for(int i = 0; i < nsub; ++i) sd[i]->multAddLT(f.subData(i), f_g);
#ifdef DISTRIBUTED
  if(com) com->globalSum(neq, f_g);
#endif
}

template class MultiDomainSolver<double>;
template class MultiDomainSolver<std::complex<double>>;