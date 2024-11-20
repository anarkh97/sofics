#ifdef USE_EIGEN3
#include <Eigen/Core>
#include <Timers.d/GetTime.h>
#include <algorithm>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <list>
#include <utility>
#include <vector>
#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

struct double_int {
  double val;
  int rank;
};

struct long_int {
  long index;
  int sub;
};

bool operator== (const long_int& lhs, const long_int& rhs);

// Parallel implementation of Non-negative Conjugate Gradient Pursuit using either MPI, OpenMP or both.
// This is a non-negative constrained variant of Conjugate Gradient Pursuit, and generates identical iterates to
// Lawson & Hanson's NNLS (non-negative least squares) algorithm.
// References:
// 1. Blumensath, Thomas, and Michael E. Davies. "Gradient pursuits." Signal Processing, IEEE Transactions on 56.6 (2008): 2370-2382.
// 2. Lawson, C. L., & Hanson, R. J. (1974). Solving least squares problems (Vol. 161). Englewood Cliffs, NJ: Prentice-hall.

Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1>
pomp(const std::vector<Eigen::Map<Eigen::MatrixXd> >&A, const Eigen::Ref<const Eigen::VectorXd> &b, double& rnorm, const long int n,
       long int &info, double maxsze, double maxite, double reltol, bool verbose, bool scaling, bool positive, double &dtime)
{
  // each A[i] is the columnwise block of the global A matrix assigned to a subdomain on this mpi process
  // each x[i] of the return value x is the corresponding row-wise block of the global solution vector
  // n is the number of columns in the global A matrix (note this is only used to define the stopping criteria)
  using namespace Eigen;
  int myrank;
#ifdef USE_MPI
  const MPI_Comm mpicomm = MPI_COMM_WORLD;
  MPI_Comm_rank(mpicomm, &myrank);
#else
  myrank = 0;
#endif
  struct double_int s;
  struct long_int p;

  const int nsub = A.size(); // number of subdomains assigned to this mpi process
  const long int m = b.rows();
  const long int maxvec = std::min(m, (long int)(maxsze*n));
  const long int maxit = maxite*n;
  std::vector<long int> maxlocvec(nsub); for(int i=0; i<nsub; ++i) maxlocvec[i] = std::min(A[i].cols(),maxvec);
  double bnorm = b.norm();
  double abstol = reltol*bnorm;

  Array<VectorXd,Dynamic,1> x_(nsub), y(nsub), g(nsub), h(nsub), g_(nsub), S(nsub), t(nsub);
  for(int i=0; i<nsub; ++i) {
    x_[i] = VectorXd::Zero(maxlocvec[i]);
    y[i].resize(maxlocvec[i]);
    g_[i].resize(maxlocvec[i]);
    S[i].resize(A[i].cols());
    if(scaling) for(int j=0; j<A[i].cols(); ++j) { double s = A[i].col(j).norm(); S[i][j] = (s != 0) ? 1/s : 0; }
    else S[i].setOnes();
    t[i].resize(maxlocvec[i]);
  }
  VectorXd r(m), DtGDinv(maxvec), z(maxvec), a(maxvec);
  r = b;
  rnorm = bnorm;
  info = (n < 0) ? 2 : 1;
  a.setZero();
  Array<MatrixXd,Dynamic,1> B(nsub), D(nsub), GD(nsub);
  for(int i=0; i<nsub; ++i) {
    B[i]  = MatrixXd::Zero(m,maxlocvec[i]);
    D[i]  = MatrixXd::Zero(maxlocvec[i],maxvec);
    GD[i] = MatrixXd::Zero(maxlocvec[i],maxvec);
  }
  Matrix<double,Dynamic,Dynamic,ColMajor> BD(m,maxvec);
  MatrixXd z_loc(maxvec,nsub), c_loc(m,nsub);
  long int k = 0; // k is the dimension of the (global) set of selected indices
  std::vector<long int> l(nsub); // l[i] is the dimension of the subset of selected indices local to a subdomain.
  std::list<std::pair<int,long_int> > gindices; // global indices
  std::list<std::pair<int,long_int> > nld_indices;
  std::vector<std::vector<long int> > indices(nsub); // local indices

  Array<long int,Dynamic,1> jk(nsub);
  ArrayXd gmax(nsub), ymin(nsub), alpha(nsub);

  long int iter   = 0; // number of iterations
  long int downIt = 0;
  while(true) {

    if(myrank == 0 && verbose) {
      std::cout << "Iteration = " << std::setw(9) << iter << "    "
                << "Downdate = " << std::setw(9) << downIt << "    "
                << "Active set size = " << std::setw(9) << k << "    "
                << "Residual norm = " << std::setw(13) << std::scientific << std::uppercase << rnorm << "    "
                << "Target = " << std::setw(13) << std::scientific << std::uppercase << abstol << std::endl;
      std::cout.unsetf(std::ios::scientific);
      std::cout.unsetf(std::ios::uppercase);
    }

    if(iter >= maxit) { info = 3; break; }

#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) {
      g[i] = S[i].asDiagonal()*(A[i].transpose()*r);
      // make sure the index has not already been selected
      h[i] = g[i]; for(long int j=0; j<l[i]; ++j) h[i][indices[i][j]] = -std::numeric_limits<double>::max();
      // also make sure near linear dependent indices are not selected
      for(std::list<std::pair<int,long_int> >::iterator it = nld_indices.begin(); it != nld_indices.end(); ++it)
        if(it->first == myrank && it->second.sub == i) h[i][it->second.index] = -std::numeric_limits<double>::max();
      gmax[i] = (A[i].cols() > 0) ? h[i].maxCoeff(&jk[i]) : -std::numeric_limits<double>::max();
    }
    int ik; // subdomain which has the max coeff.
    s.val  = (nsub > 0) ? gmax.maxCoeff(&ik) : -std::numeric_limits<double>::max();
    s.rank = myrank;
#ifdef USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &s, 1, MPI_DOUBLE_INT, MPI_MAXLOC, mpicomm);
#endif
    if(s.rank == myrank) {
      p.index = jk[ik];
      p.sub   = ik;
    }
#ifdef USE_MPI
    MPI_Bcast(&p, 1, MPI_LONG_INT, s.rank, mpicomm);
#endif
    if(s.val <= 0) break;
    if(s.rank == myrank) {
      B[ik].col(l[ik]) = S[ik][jk[ik]]*A[ik].col(jk[ik]);
      // update GD due to extra column added to B (note: B.col(i)*D.row(i).head(k) = 0, so BD does not need to be updated)
      GD[ik].row(l[ik]).head(k) = B[ik].col(l[ik]).transpose()*BD.leftCols(k);
      indices[ik].push_back(jk[ik]);
      l[ik]++;
    }
    gindices.push_back(std::pair<int,long_int>(s.rank,p));

#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) {
      for(long int j=0; j<l[i]; ++j) g_[i][j] = g[i][indices[i][j]];
      z_loc.col(i).head(k) = GD[i].topLeftCorner(l[i],k).transpose()*g_[i].head(l[i]);
    }
    z.head(k) = z_loc.topRows(k).rowwise().sum();
#ifdef USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, z.data(), k, MPI_DOUBLE, MPI_SUM, mpicomm);
#endif

    z.head(k) = (DtGDinv.head(k).asDiagonal()*z.head(k)).eval();
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) {
      Block<MatrixXd,Dynamic,1,true> d = D[i].col(k);
      d.head(l[i]) = g_[i].head(l[i]) - D[i].topLeftCorner(l[i],k)*z.head(k);
      c_loc.col(i) = B[i].leftCols(l[i])*d.head(l[i]);
    }
    Block<MatrixXd,Dynamic,1,true> c = BD.col(k) = c_loc.rowwise().sum();
#ifdef USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, c.data(), m, MPI_DOUBLE, MPI_SUM, mpicomm);
#endif
    DtGDinv[k] = 1/c.squaredNorm();
    a[k] = r.dot(c)*DtGDinv[k];
    if(a[k] < 0) { // check for near linear dependence
      nld_indices.push_back(std::pair<int,long_int>(s.rank,p));
      gindices.pop_back();
      if(s.rank == myrank) {
        indices[ik].pop_back();
        l[ik]--;
      }
      continue;
    } else nld_indices.clear();
    r -= a[k]*c;
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) {
      Block<MatrixXd,Dynamic,1,true> d = D[i].col(k);
      y[i].head(l[i]) = x_[i].head(l[i]) + a[k]*d.head(l[i]);
      GD[i].col(k).head(l[i]) = B[i].leftCols(l[i]).transpose()*c;
    }
    k++;
    iter++;
    if((rnorm <= abstol || k+nld_indices.size() == maxvec) && positive) {
      while(true) {
        iter++;
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
        for(int i=0; i<nsub; ++i) {
          ymin[i] = (l[i] > 0) ? y[i].head(l[i]).minCoeff() : std::numeric_limits<double>::max();
        }
        double minCoeff = (nsub > 0) ? ymin.minCoeff() : std::numeric_limits<double>::max();
#ifdef USE_MPI
        MPI_Allreduce(MPI_IN_PLACE, &minCoeff, 1, MPI_DOUBLE, MPI_MIN, mpicomm);
#endif
        if(minCoeff < 0) {
          dtime -= getTime();
          downIt++;
          if(myrank == 0 && verbose) {
          std::cout << "Iteration = " << std::setw(9) << iter << "    "
                    << "Downdate = " << std::setw(9) << downIt << "    "
                    << "Active set size = " << std::setw(9) << k << "    "
                    << "Residual norm = " << std::setw(13) << std::scientific << std::uppercase << rnorm << "    "
                    << "Target = " << std::setw(13) << std::scientific << std::uppercase << abstol << std::endl;
          std::cout.unsetf(std::ios::scientific);
          std::cout.unsetf(std::ios::uppercase);
          }
          // compute maximum feasible step length in the direction (y-x_) and corresponding index in active set jk[i] for each subdomain
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
          for(int i=0; i<nsub; ++i) {
            for(long int j=0; j<l[i]; ++j) t[i][j] = (y[i][j] >= 0) ? std::numeric_limits<double>::max() : -x_[i][j]/(y[i][j]-x_[i][j]);
            alpha[i] = (l[i] > 0) ? t[i].head(l[i]).minCoeff(&jk[i]) : std::numeric_limits<double>::max();
          }

          int ik; // subdomain which has the maximum feasible step length.
          s.val = (nsub > 0) ? alpha.minCoeff(&ik) : std::numeric_limits<double>::max();
          s.rank = myrank;
#ifdef USE_MPI
          MPI_Allreduce(MPI_IN_PLACE, &s, 1, MPI_DOUBLE_INT, MPI_MINLOC, mpicomm);
#endif
          if(s.rank == myrank) {
            p.index = indices[ik][jk[ik]];
            p.sub = ik;
            indices[ik].erase(indices[ik].begin() + jk[ik]); // remove index jk[ik] from the local active set of ik-th subdomain
            for(int j=jk[ik]; j<l[ik]-1; ++j) { x_[ik][j] = x_[ik][j+1]; y[ik][j] = y[ik][j+1]; } // erase jk[ik]-th element from x_[ik] and y[ik]
            x_[ik][l[ik]-1] = 0; y[ik][l[ik]-1] = 0;
            l[ik]--;
            //std::cout << "removing index " << p.index << " from subdomain " << ik << " on process with rank " << myrank << std::endl;
          }

#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
          // update x_ (note: this is used only when there are two or more consecutive downdate iterations)
          for(int i=0; i<nsub; ++i) {
            x_[i].head(l[i]) += s.val*(y[i].head(l[i])-x_[i].head(l[i]));
          }

#ifdef USE_MPI
          MPI_Bcast(&p, 1, MPI_LONG_INT, s.rank, mpicomm);
#endif
          // remove index from the global active set
          std::list<std::pair<int,long_int> >::iterator pos = std::find(gindices.begin(), gindices.end(), std::pair<int,long_int>(s.rank,p));
          std::list<std::pair<int,long_int> >::iterator fol = gindices.erase(pos);

          // Note: it is necessary to re-G-orthogonalize the basis D now, project the solution y onto the new basis and compute the corresponding residual r.
          // This is done here by starting from the column of D pointed to by fol (because the ones before this are already G-orthogonal), and then
          // following what is the essentially same procedure that is used above to construct the original basis, with a few optimizations when possible.
          k = std::distance(gindices.begin(), fol);
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
          for(int i=0; i<nsub; ++i) {
            y[i].head(l[i]).setZero();
            l[i] = 0; for(std::list<std::pair<int,long_int> >::iterator it = gindices.begin(); it != fol; ++it) { if(it->first == myrank && it->second.sub == i) l[i]++; }
            D[i].bottomRightCorner(maxlocvec[i]-l[i],maxvec-k).setZero();
            y[i].head(l[i]) = D[i].topLeftCorner(l[i],k)*a.head(k);
          }
          r = b - BD.leftCols(k)*a.head(k);
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
          for(int i=0; i<nsub; ++i) {
            g_[i].head(l[i]) = B[i].leftCols(l[i]).transpose()*r;
          }
          for(std::list<std::pair<int,long_int> >::iterator it = fol; it != gindices.end(); ++it) {
            if(it->first == myrank) {
              int i = it->second.sub;
              B[i].col(l[i]) = S[i][it->second.index]*A[i].col(it->second.index);
              g_[i][l[i]] = B[i].col(l[i]).transpose()*r;
              // update GD due to extra column added to B (note: B.col(i)*D.row(i).head(k) = 0, so BD does not need to be updated)
              GD[i].row(l[i]).head(k) = B[i].col(l[i]).transpose()*BD.leftCols(k);
              l[i]++;
            }

#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
            for(int i=0; i<nsub; ++i) {
              z_loc.col(i).head(k) = GD[i].topLeftCorner(l[i],k).transpose()*g_[i].head(l[i]);
            }
            z.head(k) = z_loc.topRows(k).rowwise().sum();
#ifdef USE_MPI
            MPI_Allreduce(MPI_IN_PLACE, z.data(), k, MPI_DOUBLE, MPI_SUM, mpicomm);
#endif

            z.head(k) = (DtGDinv.head(k).asDiagonal()*z.head(k)).eval();
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
            for(int i=0; i<nsub; ++i) {
              Block<MatrixXd,Dynamic,1,true> d = D[i].col(k);
              d.head(l[i]) = g_[i].head(l[i]) - D[i].topLeftCorner(l[i],k)*z.head(k);
              c_loc.col(i) = B[i].leftCols(l[i])*d.head(l[i]);
            }
            Block<MatrixXd,Dynamic,1,true> c = BD.col(k) = c_loc.rowwise().sum();
#ifdef USE_MPI
            MPI_Allreduce(MPI_IN_PLACE, c.data(), m, MPI_DOUBLE, MPI_SUM, mpicomm);
#endif

            DtGDinv[k] = 1/c.squaredNorm();
            a[k] = r.dot(c)*DtGDinv[k];
            r -= a[k]*c;
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
            for(int i=0; i<nsub; ++i) {
              Block<MatrixXd,Dynamic,1,true> d = D[i].col(k);
              y[i].head(l[i]) += a[k]*d.head(l[i]);
              GD[i].col(k).head(l[i]) = B[i].leftCols(l[i]).transpose()*c;
              g_[i].head(l[i]) -= a[k]*GD[i].col(k).head(l[i]);
            }
            k++;
          }
          dtime += getTime();
        }
        else {
          rnorm = r.norm();
          break;
        }
      }
    }

#if defined(_OPENMP)
#pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) {
      x_[i].head(l[i]) = y[i].head(l[i]);
    }

    if(rnorm <= abstol || k+nld_indices.size() == maxvec) break;
    rnorm = r.norm();
  }

  dtime /= 1000.0;
  if(myrank == 0 && verbose) std::cout.flush();

  Array<VectorXd,Dynamic,1> x(nsub);
  for(int i=0; i<nsub; ++i) {
    x[i] = VectorXd::Zero(A[i].cols());
    for(long int j=0; j<l[i]; ++j) x[i][indices[i][j]] = S[i][indices[i][j]]*x_[i][j];
  }
  return x;
}

#endif
