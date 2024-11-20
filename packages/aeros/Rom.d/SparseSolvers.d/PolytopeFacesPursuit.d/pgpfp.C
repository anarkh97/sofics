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
// This is a non-negative constrained variant of Conjugate Gradient Pursuit, and generates identical iterates (up to round-off error) to
// Lawson & Hanson's NNLS (non-negative least squares) algorithm.
// References:
// 1. Blumensath, Thomas, and Michael E. Davies. "Gradient pursuits." Signal Processing, IEEE Transactions on 56.6 (2008): 2370-2382.
// 2. Lawson, C. L., & Hanson, R. J. (1974). Solving least squares problems (Vol. 161). Englewood Cliffs, NJ: Prentice-hall.

Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1>
pgpfp(std::vector<Eigen::Map<Eigen::MatrixXd> >&A, Eigen::Ref<Eigen::VectorXd> b, double& rnorm, const long int n,
       long int &info, double maxsze, double maxite, double reltol, bool verbose, bool scaling, bool center, bool positive, double &dtime)
{
  // each A[i] is the columnwise block of the global A matrix assigned to a subdomain on this mpi process
  // each x[i] of the return value x is the corresponding row-wise block of the global solution vector
  // n is the number of columns in the global A matrix (note this is only used to define the stopping criteria)
  using namespace Eigen;
  const int nsub = A.size(); // number of subdomains assigned to this mpi process
  int myrank, numproc;
#ifdef USE_MPI
  MPI_Comm mpicomm;
  MPI_Comm_split(MPI_COMM_WORLD, int(nsub > 0), 0, &mpicomm);
  if(nsub == 0) return Array<VectorXd,Dynamic,1>(0,1);
  MPI_Comm_rank(mpicomm, &myrank);
  MPI_Status status;
  MPI_Comm_size(mpicomm, &numproc);
#else
  myrank = 0;
  numproc = 1;
#endif
  struct double_int s;
  struct long_int p, q;

  const long int m = b.rows();
  const long int maxvec = std::min(m, (long int)(maxsze*n));
  const long int maxit = maxite*n;
  std::vector<long int> maxlocvec(nsub); for(int i=0; i<nsub; ++i) maxlocvec[i] = std::min(A[i].cols(),maxvec);
  double bnorm = b.norm();
  double abstol = reltol*bnorm;

  int *counts = new int[numproc];
  int *displs = new int[numproc];
  for(int i=0; i<numproc; ++i) { counts[i] = m/numproc + ((i < m%numproc) ? 1 : 0); displs[i] = (i==0) ? 0 : displs[i-1]+counts[i-1]; }
#ifdef USE_MPI
  VectorXd c(m), buf1(m), buf2(maxvec), buf3(counts[myrank]); // temporary buffers
#endif
  Array<VectorXd,Dynamic,1> x_(nsub), y(nsub), g(nsub), h(nsub), g_(nsub), S(nsub), t(nsub), g1(nsub), g2(nsub);
  for(int i=0; i<nsub; ++i) {
    x_[i] = VectorXd::Zero(maxlocvec[i]);
    y[i].resize(maxlocvec[i]);
    g_[i].resize(maxlocvec[i]);
    g1[i].resize(A[i].cols());
    S[i].resize(A[i].cols());
    if(center) {
      for(long int j=0; j<A[i].cols(); ++j) { A[i].col(j).array() -= A[i].col(j).mean(); }
      b.array() -= b.mean();
    }
    if(scaling) for(long int j=0; j<A[i].cols(); ++j) { double s = A[i].col(j).norm(); S[i][j] = (s != 0) ? 1/s : 0; }
    else S[i].setOnes();
    t[i].resize(maxlocvec[i]);
  }
  VectorXd r(m), vertex(m), DtGDinv(maxvec), z(maxvec), a(maxvec), lambda(maxvec);
  r = b;
  vertex.setZero();
  rnorm = bnorm;
  info = (n < 0) ? 2 : 1;
  a.setZero();
  lambda.setZero();
  Array<MatrixXd,Dynamic,1> B(nsub), D(nsub);
  Array<Matrix<double,Dynamic,Dynamic,RowMajor>,Dynamic,1> GD(nsub);
  for(int i=0; i<nsub; ++i) {
    B[i].resize(m,maxlocvec[i]);
    D[i].resize(maxlocvec[i],maxvec);
    GD[i].resize(maxlocvec[i],maxvec);
  }

  Matrix<double,Dynamic,Dynamic,ColMajor> BD(counts[myrank],maxvec);
  MatrixXd z_loc(maxvec,nsub), c_loc(m,nsub);
  long int k = 0; // k is the dimension of the (global) set of selected indices
  std::vector<long int> l(nsub); // l[i] is the dimension of the subset of selected indices local to a subdomain.
  std::list<std::pair<int,long_int> > gindices; // global indices
  std::list<std::pair<int,long_int> > nld_indices;
  std::vector<int> nld_setKey;
  std::vector<std::vector<long int> > indices(nsub); // local indices
  std::vector<std::vector<int> > setKey(nsub);
  std::vector<std::vector<std::pair<int,long_int> > > perm(nsub); // column permutation
  for(int i=0; i<nsub; ++i) { // initialize the column permutation
    perm[i].reserve(A[i].cols());
    for(long int j=0; j<A[i].cols(); ++j) {
      p.sub = i;
      p.index = j;
      perm[i][j] = std::pair<int,long_int>(myrank,p);
    }
  }

  Array<long int,Dynamic,1> jk(nsub);
  Array<int,Dynamic,1> setMax(nsub);
  Array<double,Dynamic,1> lamMax(nsub);
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

    if(rnorm <= abstol || k+nld_indices.size() == maxvec) break;
    if(iter >= maxit) { info = 3; break; }

#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) {
      g[i]  = S[i].asDiagonal()*(A[i].transpose()*r);
      g2[i] = S[i].asDiagonal()*(A[i].transpose()*vertex);

      for(long int col = 0; col != A[i].cols(); col++) {
        double num = g[i][col];
        double den = g2[i][col];
        if(num > 0. && den != 1.0) {
          g1[i][col] = num/(1.0-den);
          if(!positive)
            g2[i][col] = -std::numeric_limits<double>::max();
        } else if (num < 0. && den != -1.0) {
          g1[i][col] = -std::numeric_limits<double>::max();
          if(!positive)
            g2[i][col] = (-1.0)*num/(1.0+den);
        } else {
          g1[i][col] = -std::numeric_limits<double>::max();
          if(!positive)
            g2[i][col] = -std::numeric_limits<double>::max();
        }
      }
      // make sure that element has not already been selected
      for(long int j=0; j<l[i]; ++j) {
        if(setKey[i][j] == 1) {
          g1[i][indices[i][j]] = -std::numeric_limits<double>::max();
        } else if(!positive) {
          g2[i][indices[i][j]] = -std::numeric_limits<double>::max();
        }
      }
      // also make sure near linear dependent indices are not selected
      std::vector<int>::iterator keyIt = nld_setKey.begin();
      for(std::list<std::pair<int,long_int> >::iterator it = nld_indices.begin(); it != nld_indices.end(); ++it, ++keyIt) {
        if(it->first == myrank && it->second.sub == i) {
          if(*keyIt == 1) {
            g1[i][it->second.index] = -std::numeric_limits<double>::max();
          } else if(!positive) {
            g2[i][it->second.index] = -std::numeric_limits<double>::max();
          }
        }
      }

      int whichSet = 1;
      long int position1 = 0;
      long int position2 = 0;

      double lam  = 0.0;
      double lam1 = g1[i].maxCoeff(&position1);
      double lam2 = (!positive) ? g2[i].maxCoeff(&position2) : 0;

      if(A[i].cols() > 0) {
        if(lam1 > lam2 || positive) {
          lam = lam1;
          jk[i] = position1;
        } else {
          lam = lam2;
          jk[i] = position2;
          whichSet = -1;
        }
      } else {
        lam = -std::numeric_limits<double>::max();
      }

      lamMax[i] = lam;
      setMax[i] = whichSet;
    }
    int ik; // subdomain which has the max coeff.
    int Set; // correct sign of selected set
    s.val  = (nsub > 0) ? lamMax.maxCoeff(&ik) : -std::numeric_limits<double>::max();
    s.rank = myrank;
#ifdef USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &s, 1, MPI_DOUBLE_INT, MPI_MAXLOC, mpicomm);
#endif
    if(s.val <= 0) break;
    lambda[k] = 1.0/s.val;

    q.index = 0; q.sub = myrank;
    for(int i=0; i<nsub; ++i) q.index += l[i];
#ifdef USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &q, 1, MPI_LONG_INT, MPI_MINLOC, mpicomm); // find mpi process with minimum inactive set
#endif
    //if(myrank == 0) std::cerr << "rank of process with fewest elements (" << q.index << ") is " << q.sub << std::endl;

#ifdef USE_MPI
    if(q.sub != s.rank) {
      if(s.rank == myrank) { // exchange column of A with process with fewest elements
        MPI_Sendrecv_replace(const_cast<double*>(A[ik].col(jk[ik]).data()), m, MPI_DOUBLE,
                             q.sub, 0, q.sub, 0, mpicomm, &status);                                            // send column of A
        if(scaling) MPI_Sendrecv_replace(&S[ik][jk[ik]], 1, MPI_DOUBLE, q.sub, 1, q.sub, 1, mpicomm, &status); // send scaling factor
        MPI_Sendrecv_replace(&perm[ik][jk[ik]].first, 1, MPI_INT, q.sub, 2, q.sub, 2, mpicomm, &status);       // send permutation 1
        MPI_Sendrecv_replace(&perm[ik][jk[ik]].second, 1, MPI_LONG_INT, q.sub, 3, q.sub, 3, mpicomm, &status); // send permutation 2
      }
      if(q.sub == myrank) {
        MPI_Sendrecv_replace(const_cast<double*>(A[ik].col(jk[ik]).data()), m, MPI_DOUBLE, 
                             s.rank, 0, s.rank, 0, mpicomm, &status);                                            // recieve column A 
        if(scaling) MPI_Sendrecv_replace(&S[ik][jk[ik]], 1, MPI_DOUBLE, s.rank, 1, s.rank, 1, mpicomm, &status); // recieve scaling factor
        MPI_Sendrecv_replace(&perm[ik][jk[ik]].first, 1, MPI_INT, s.rank, 2, s.rank, 2, mpicomm, &status);       // recieve permuation 1
        MPI_Sendrecv_replace(&perm[ik][jk[ik]].second, 1, MPI_LONG_INT, s.rank, 3, s.rank, 3, mpicomm, &status); // recieve permutation 2
        p.index = jk[ik];
        p.sub   = ik;
      }
      s.rank = q.sub;
    }
#endif

    if(s.rank == myrank) {
      Set = setMax[ik];
      B[ik].col(l[ik]) = double(Set)*S[ik][jk[ik]]*A[ik].col(jk[ik]);
      // update GD due to extra column added to B (note: B.col(i)*D.row(i).head(k) = 0, so BD does not need to be updated)
#ifdef USE_MPI
      MPI_Scatterv(B[ik].col(l[ik]).data(), counts, displs, MPI_DOUBLE, buf3.data(), counts[myrank], MPI_DOUBLE, s.rank, mpicomm);
      buf2.head(k) = buf3.transpose()*BD.leftCols(k);
      MPI_Reduce(buf2.data(), GD[ik].row(l[ik]).data(), int(k), MPI_DOUBLE, MPI_SUM, s.rank, mpicomm);
#else
      GD[ik].row(l[ik]).head(k) = B[ik].col(l[ik]).transpose()*BD.leftCols(k);
#endif
      indices[ik].push_back(jk[ik]);
      setKey[ik].push_back(Set);
      l[ik]++;
      p.index = jk[ik];
      p.sub   = ik;
    }
#ifdef USE_MPI
    else {
      MPI_Scatterv(0, counts, displs, MPI_DOUBLE, buf3.data(), counts[myrank], MPI_DOUBLE, s.rank, mpicomm);
      buf2.head(k) = buf3.transpose()*BD.leftCols(k);
      MPI_Reduce(buf2.data(), 0, int(k), MPI_DOUBLE, MPI_SUM, s.rank, mpicomm);
    }

    MPI_Bcast(&Set, 1, MPI_INT, s.rank, mpicomm);
    MPI_Bcast(&p, 1, MPI_LONG_INT, s.rank, mpicomm);
#endif
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
      d.setZero();
      d.head(l[i]) = g_[i].head(l[i]) - D[i].topLeftCorner(l[i],k).triangularView<Upper>()*z.head(k);
      c_loc.col(i) = B[i].leftCols(l[i])*d.head(l[i]);
    }
#ifdef USE_MPI
    c = c_loc.rowwise().sum();
    MPI_Allreduce(MPI_IN_PLACE, c.data(), m, MPI_DOUBLE, MPI_SUM, mpicomm);
    BD.col(k) = c.segment(displs[myrank],counts[myrank]);
#else
    Block<MatrixXd,Dynamic,1,true> c = BD.col(k) = c_loc.rowwise().sum();
#endif
    DtGDinv[k] = 1/c.squaredNorm();
    a[k] = r.dot(c)*DtGDinv[k];
    if(a[k] < 0) { // check for near linear dependence
      nld_indices.push_back(std::pair<int,long_int>(s.rank,p));
      nld_setKey.push_back(Set);
      gindices.pop_back();
      if(s.rank == myrank) {
        indices[ik].pop_back();
        l[ik]--;
      }
      continue;
    } else nld_indices.clear();
    r -= a[k]*c;
    vertex += lambda[k]*a[k]*c;
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) {
      Block<MatrixXd,Dynamic,1,true> d = D[i].col(k);
      y[i].head(l[i]) = x_[i].head(l[i]) + a[k]*d.head(l[i]);
      GD[i].col(k).head(l[i]) = B[i].leftCols(l[i]).transpose()*c;
    }
    k++;

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
        // compute maximum feasible step length in the direction (y-x_) and corresponding index in active set jk[i] for each subdomain
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
        for(int i=0; i<nsub; ++i) {
          for(long int j=0; j<l[i]; ++j) t[i][j] = (y[i][j] >= 0) ? std::numeric_limits<double>::max() : -x_[i][j]/(y[i][j]-x_[i][j]);
          alpha[i] = (l[i] > 0) ? t[i].head(l[i]).minCoeff(&jk[i]) : std::numeric_limits<double>::max();
        }

        int ik; // subdomain which has the smallest maximum feasible step length.
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
        long int j = k;
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
        for(int i=0; i<nsub; ++i) {
          y[i].head(l[i]).setZero();
          l[i] = 0; for(std::list<std::pair<int,long_int> >::iterator it = gindices.begin(); it != fol; ++it) { if(it->first == myrank && it->second.sub == i) l[i]++; }
          y[i].head(l[i]) = D[i].topLeftCorner(l[i],k).triangularView<Upper>()*a.head(k);
        }
#ifdef USE_MPI
        if(k < gindices.size()/2){ //if removing an early column, recompute residual
          buf1.segment(displs[myrank],counts[myrank]) = BD.leftCols(k)*a.head(k);
          MPI_Allgatherv(MPI_IN_PLACE, counts[myrank], MPI_DOUBLE, buf1.data(), counts, displs, MPI_DOUBLE, mpicomm);
          r = b - buf1;
        } else { // if removing a late column, update residual instead. 
          buf1.segment(displs[myrank],counts[myrank]) = BD.block(0,k,counts[myrank],gindices.size()-k+1)*a.segment(k,gindices.size()-k+1);
          MPI_Allgatherv(MPI_IN_PLACE, counts[myrank], MPI_DOUBLE, buf1.data(), counts, displs, MPI_DOUBLE, mpicomm);
          r += buf1;
        }
        buf1.segment(displs[myrank],counts[myrank]) = BD.leftCols(k)*lambda.head(k).asDiagonal()*a.head(k);
        MPI_Allgatherv(MPI_IN_PLACE, counts[myrank], MPI_DOUBLE, buf1.data(), counts, displs, MPI_DOUBLE, mpicomm);
        vertex = buf1; 
#else
        if(k < gindices.size()/2) 
          r = b - BD.leftCols(k)*a.head(k);
        else 
          r += BD.block(0,k,m,gindices.size()-k+1)*a.segment(k,gindices.size()-k+1);

        vertex = BD.leftCols(k)*lambda.head(k).asDiagonal()*a.head(k);
#endif
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
        for(int i=0; i<nsub; ++i) {
          g_[i].head(l[i]) = B[i].leftCols(l[i]).transpose()*r;
        }
        for(std::list<std::pair<int,long_int> >::iterator it = fol; it != gindices.end(); ++it) {
          if(it->first == myrank) {
            int i = it->second.sub;
            int slot = std::distance(indices[i].begin(), std::find(indices[i].begin(), indices[i].end(), it->second.index));
            B[i].col(l[i]) = double(setKey[i][slot])*S[i][it->second.index]*A[i].col(it->second.index);
            double num = B[i].col(l[i]).transpose()*r;
            double den = B[i].col(l[i]).dot(vertex);
            g_[i][l[i]] = num;
            lambda[k] = (1.0-den)/num;
            // update GD due to extra column added to B (note: B.col(i)*D.row(i).head(k) = 0, so BD does not need to be updated)
            if(s.rank == myrank && p.sub == i) GD[i].row(l[i]).head(j) = GD[i].row(l[i]+1).head(j);
#ifdef USE_MPI
            MPI_Scatterv(B[i].col(l[i]).data(), counts, displs, MPI_DOUBLE, buf3.data(), counts[myrank], MPI_DOUBLE, it->first, mpicomm);
            buf2.head(k-j) = buf3.transpose()*BD.block(0,j,counts[myrank],k-j);
            MPI_Reduce(buf2.data(), GD[i].row(l[i]).data()+j, int(k-j), MPI_DOUBLE, MPI_SUM, it->first, mpicomm);
#else
            GD[i].row(l[i]).segment(j,k-j) = B[i].col(l[i]).transpose()*BD.block(0,j,counts[myrank],k-j);
#endif
            l[i]++;
          }
#ifdef USE_MPI
          else {
            MPI_Scatterv(0, counts, displs, MPI_DOUBLE, buf3.data(), counts[myrank], MPI_DOUBLE, it->first, mpicomm);
            buf2.head(k-j) = buf3.transpose()*BD.block(0,j,counts[myrank],k-j);
            MPI_Reduce(buf2.data(), 0, int(k-j), MPI_DOUBLE, MPI_SUM, it->first, mpicomm);
          }
#endif

#ifdef USE_MPI
          MPI_Bcast(&lambda[k], 1, MPI_DOUBLE, it->first, mpicomm);
#endif

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
            d.setZero();
            d.head(l[i]) = g_[i].head(l[i]) - D[i].topLeftCorner(l[i],k).triangularView<Upper>()*z.head(k);
            c_loc.col(i) = B[i].leftCols(l[i])*d.head(l[i]);
          }
#ifdef USE_MPI
          c = c_loc.rowwise().sum();
          MPI_Allreduce(MPI_IN_PLACE, c.data(), m, MPI_DOUBLE, MPI_SUM, mpicomm);
          BD.col(k) = c.segment(displs[myrank],counts[myrank]);
#else
          Block<MatrixXd,Dynamic,1,true> c = BD.col(k) = c_loc.rowwise().sum();
#endif
          DtGDinv[k] = 1/c.squaredNorm();
          a[k] = r.dot(c)*DtGDinv[k];
          r -= a[k]*c;
          vertex += lambda[k]*a[k]*c;
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
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
        for(int i=0; i<nsub; ++i) {
          x_[i].head(l[i]) = y[i].head(l[i]);
        }
        break;
      }
    }

    rnorm = r.norm();
  }

  delete [] counts;
  delete [] displs;

  dtime /= 1000.0;
  if(myrank == 0 && verbose) std::cout.flush();

  Array<VectorXd,Dynamic,1> x(nsub);
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
  for(int i=0; i<nsub; ++i) {
    x[i] = VectorXd::Zero(A[i].cols());
    for(long int j=0; j<l[i]; ++j) x[i][indices[i][j]] = S[i][indices[i][j]]*x_[i][j];
  }

#ifdef USE_MPI
  Array<VectorXd,Dynamic,1> Px(nsub);
  int *toSend = new int[numproc];
  int maxindex = 0;
  for(int i=0; i<numproc; ++i) toSend[i] = 0;
  for(int i=0; i<nsub; ++i) {
    Px[i] = VectorXd::Zero(A[i].cols());
    for(long int j=0; j<l[i]; ++j) {
      const int &dest = perm[i][indices[i][j]].first;
      if(dest != myrank) toSend[dest] += 1;
    }
    maxindex = std::max(int(A[i].cols()),maxindex);
  }
  MPI_Allreduce(MPI_IN_PLACE, toSend, numproc, MPI_INT, MPI_SUM, mpicomm);
  MPI_Allreduce(MPI_IN_PLACE, &maxindex, 1, MPI_INT, MPI_MAX, mpicomm);
  for(int i=0; i<nsub; ++i) {
    for(long int j=0; j<l[i]; ++j) {
      const long int &k = indices[i][j];
      const int &dest = perm[i][k].first;
      if(dest != myrank) {
        double_int toto;
        toto.val = x[i][k];
        toto.rank = perm[i][k].second.sub*maxindex + int(perm[i][k].second.index);
        MPI_Send(&toto, 1, MPI_DOUBLE_INT, dest, 0, mpicomm);
      }
      else {
        Px[perm[i][k].second.sub][perm[i][k].second.index] = x[i][k];
      }
    }
  }
  for(int i=0; i<toSend[myrank]; ++i) {
    double_int toto;
    MPI_Recv(&toto, 1, MPI_DOUBLE_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, mpicomm, &status);
    Px[toto.rank/maxindex][toto.rank%maxindex] = toto.val;
  }
  delete [] toSend;

  return Px;
#else
  return x;
#endif
}

#endif
