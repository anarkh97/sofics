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

Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1>
pnncgp(std::vector<Eigen::Map<Eigen::MatrixXd> >&A, Eigen::Ref<Eigen::VectorXd> b, double& rnorm, const long int n,
       long int &info, double maxsze, double maxite, double reltol, bool verbose, bool scaling, bool center, double &dtime, 
       std::list<std::pair<int,long_int> > &hotIndices);

Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1>
splh(const std::vector<Eigen::Map<Eigen::MatrixXd> >&A, const Eigen::Ref<const Eigen::VectorXd> &b, double& rnorm, const long int n,
       long int &info, double maxsze, double maxite, double reltol, bool verbose, bool scaling, bool positive, double &dtime,
       int npMax, int scpkMB, int scpkNB, int scpkMP, int scpkNP, int option);

bool operator== (const long_int& lhs, const long_int& rhs);

// Parallel implementation of Non-negative LASSO based on Conjugate Gradient Pursuit using either MPI, OpenMP or both.
// References:
// 1. Blumensath, Thomas, and Michael E. Davies. "Gradient pursuits." Signal Processing, IEEE Transactions on 56.6 (2008): 2370-2382.
// 2. Efron, Bradley, Hastie, Trevor, Johnstone, Iain, and Tibshirani, Robert. "Least Angle Regression" The Annals of Statistics

Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,1>
pcglars(std::vector<Eigen::Map<Eigen::MatrixXd> >&A, Eigen::Ref<Eigen::VectorXd> b, double& rnorm, const long int n,
       long int &info, double maxsze, double maxite, double reltol, bool verbose, bool scaling, bool center, bool project, double &dtime)
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

  // set up relavent problem sizes
  const long int m = b.rows();
  const long int maxvec = std::min(m, (long int)(maxsze*n));
  const long int maxit = maxite*n;
  std::vector<long int> maxlocvec(nsub); for(int i=0; i<nsub; ++i) maxlocvec[i] = std::min(A[i].cols(),maxvec);

  // stoping criteria
  double bnorm = b.norm();
  double abstol = reltol*bnorm;

  int *counts = new int[numproc]; // buffer for row-wise partition sizes, static quantity
  int *displs = new int[numproc]; // buffer for starting point in memory for vectors that are distributed row-wise 
  for(int i=0; i<numproc; ++i) { counts[i] = m/numproc + ((i < m%numproc) ? 1 : 0); displs[i] = (i==0) ? 0 : displs[i-1]+counts[i-1]; }

  VectorXd c(m), buf1(m), buf2(maxvec), buf3(counts[myrank]); // temporary buffers
   
  // allocate space for working vectors and containers
  Array<VectorXd,Dynamic,1> wA(nsub), y(nsub), g(nsub), h(nsub), minbuffer(nsub), S(nsub), t(nsub);
  VectorXd r(m), a(maxvec), update(m);
  Array<MatrixXd,Dynamic,1> B(nsub), D(nsub);
  Matrix<double,Dynamic,Dynamic,ColMajor> BD(counts[myrank],maxvec);
  MatrixXd c_loc(m,nsub);

  std::vector<long int> l(nsub);                                  // l[i] is the cardinality of the set of selected indices local to a subdomain.
  std::list<std::pair<int,long_int> > gindices;                   // global indices
  std::list<std::pair<int,long_int> > nld_indices;                // linearly dependent vectors
  std::vector<std::vector<long int> > indices(nsub);              // local indices
  std::vector<std::vector<std::pair<int,long_int> > > perm(nsub); // column permutation

  Array<long int,Dynamic,1> jk(nsub);
  ArrayXd gmax(nsub), ymin(nsub);

  double grad, DtGDinv, C, oneNwA;
  int ik; 

  // resizing and initialization
  for(int i=0; i<nsub; ++i) {
    wA[i] = VectorXd::Zero(maxlocvec[i]); // solution of system inversion
    y[i] = VectorXd::Zero(maxlocvec[i]);  // lasso solution
    S[i].resize(A[i].cols());             // column scaling
    if(center) { // center predictors and target
      if(myrank == 0) std::cout << " ... Centering Covariates ... " << std::endl;
      for(long int j=0; j<A[i].cols(); ++j) { A[i].col(j).array() -= A[i].col(j).mean(); }
      b.array() -= b.mean();
    }
    // scale predictors
    if(scaling) for(long int j=0; j<A[i].cols(); ++j) { double s = A[i].col(j).norm(); S[i][j] = (s != 0) ? 1/s : 0; }
    else S[i].setOnes();
    t[i].resize(maxlocvec[i]); // minimum angle at which lasso solution changes sign
    t[i].setZero();
  }

  r = b; // target
  rnorm = bnorm;
  info = (n < 0) ? 2 : 1;
  a.setZero();
  for(int i=0; i<nsub; ++i) {
    B[i].resize(m,maxlocvec[i]);
    B[i].setZero();
    D[i].resize(maxlocvec[i],maxvec);
    D[i].setZero();
  }

  for(int i=0; i<nsub; ++i) { // initialize the column permutation
    perm[i].reserve(A[i].cols());
    for(long int j=0; j<A[i].cols(); ++j) {
      p.sub = i;
      p.index = j;
      perm[i][j] = std::pair<int,long_int>(myrank,p);
    }
  }

  long int k      = 0; // k is the dimension of the (global) set of selected indices
  long int iter   = 0; // number of iterations
  long int downIt = 0; // number of downdates

  // ensure that an element is not selected immediately after being dropped
  bool     dropId  = false;
  std::pair<int,long_int> blockId; 

  // initialize starting set with maximally correlated element
  {
    // form correlation vector
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) {
      g[i] = S[i].asDiagonal()*(A[i].transpose()*r);
      gmax[i] = (A[i].cols() > 0) ? g[i].maxCoeff(&jk[i]) : -std::numeric_limits<double>::max();
    }
    // subdomain which has the max coeff.
    s.val  = (nsub > 0) ? gmax.maxCoeff(&ik) : -std::numeric_limits<double>::max();
    s.rank = myrank;
#ifdef USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &s, 1, MPI_DOUBLE_INT, MPI_MAXLOC, mpicomm);
#endif

    C = s.val;

    if(s.rank == myrank){
      B[ik].col(l[ik]) = S[ik][jk[ik]]*A[ik].col(jk[ik]);        // add local vector
      indices[ik].push_back(jk[ik]);                             // update local indices
      p.index = jk[ik];
      p.sub = ik;                                                
      grad = 1;                                                  // set last element in gradient
      Block<MatrixXd,Dynamic,1,true> d = D[ik].col(k);           
      d.setZero();
      d[l[ik]] = grad;                                       // set first conjugate directioin
      buf1 = B[ik].col(l[ik])*d.head(l[ik]+1); 
      DtGDinv = 1/buf1.squaredNorm();
      buf1 *= DtGDinv;
      a[k] = grad*grad*DtGDinv;
#ifdef USE_MPI
      MPI_Scatterv(buf1.data(), counts, displs, MPI_DOUBLE, buf3.data(), counts[myrank], MPI_DOUBLE, s.rank, mpicomm);
#else
      buf3 = buf1; 
#endif 
      l[ik]++;
    }
#ifdef USE_MPI
    else {
      MPI_Scatterv(0, counts, displs, MPI_DOUBLE, buf3.data(), counts[myrank], MPI_DOUBLE, s.rank, mpicomm);
    }
#endif
    
#ifdef USE_MPI
    MPI_Bcast(&a[k], 1, MPI_DOUBLE, s.rank, mpicomm);
#endif
     
    BD.col(k) = buf3; 
#ifdef USE_MPI
    MPI_Bcast(&p, 1, MPI_LONG_INT, s.rank, mpicomm);
#endif 
    gindices.push_back(std::pair<int,long_int>(s.rank,p)); // update global indices

  }

  while(true) {

    rnorm = r.norm();

    if(myrank == 0 && verbose) {
      std::cout << "Iteration = " << std::setw(9) << iter << "    "
                << "Downdate = " << std::setw(9) << downIt << "    "
                << "Active set size = " << std::setw(9) << k << "    "
                << "Max Correlation = " << std::setw(13) << std::scientific << std::uppercase << C << "    "
                << "Residual norm = " << std::setw(13) << std::scientific << std::uppercase << rnorm << "    "
                << "Target = " << std::setw(13) << std::scientific << std::uppercase << abstol << std::endl;
      std::cout.unsetf(std::ios::scientific);
      std::cout.unsetf(std::ios::uppercase);
    }

    if(rnorm <= abstol || k+nld_indices.size() == maxvec) break;
    if(iter >= maxit) { info = 3; break; }

    if(!dropId){ // solve square system if we havent just recomputed it
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
      for(int i=0; i<nsub; ++i) {
         wA[i].head(l[i]) += a[k]*D[i].col(k).head(l[i]);
      }
    }

#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) {// compute scalar
      gmax[i] = (l[i] > 0) ? wA[i].head(l[i]).sum() : 0;
    }
    oneNwA = gmax.sum();
#if USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &oneNwA, 1, MPI_DOUBLE, MPI_SUM, mpicomm);
#endif
    oneNwA = 1.0/sqrt(oneNwA); // all processors now have oneNwA

#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) { // compute update vector
       c_loc.col(i) = B[i].leftCols(l[i])*wA[i].head(l[i]);
    }
    update = c_loc.rowwise().sum();
#if USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, update.data(), m, MPI_DOUBLE, MPI_SUM, mpicomm); // all processor now have update, it will be reused later
#endif
    update *= oneNwA;

#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) {
      h[i] = S[i].asDiagonal()*(A[i].transpose()*update);
    }
    
    // compute LARS angle
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) {
      minbuffer[i] = (C - g[i].array())/(oneNwA - h[i].array());
      for(std::vector<long int>::iterator it = indices[i].begin(); it != indices[i].end(); ++it) minbuffer[i][*it] = std::numeric_limits<double>::max(); // max out inactive set
      for(long int row = 0; row < minbuffer[i].rows(); ++row) minbuffer[i][row] = (minbuffer[i][row] <= 0) ? std::numeric_limits<double>::max() : minbuffer[i][row];// max out non-positive members  
    }
    if(dropId && (blockId.first == myrank)) minbuffer[blockId.second.sub][blockId.second.index] = std::numeric_limits<double>::max();// max out previously rejected element

    // compute lasso angle
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) {
       t[i].head(l[i]) = -1.0*y[i].head(l[i]).array()/(oneNwA*wA[i].head(l[i]).array());
       for(long int ele = 0; ele < l[i]; ++ ele) t[i][ele] = (t[i][ele] <= 0) ? std::numeric_limits<double>::max() : t[i][ele];
       gmax[i] = (l[i] > 0) ? t[i].head(l[i]).minCoeff(&jk[i]) : std::numeric_limits<double>::max(); 
    }

    s.val = (nsub > 0) ? gmax.minCoeff(&ik) : std::numeric_limits<double>::max();
    s.rank = myrank;
#ifdef USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &s, 1, MPI_DOUBLE_INT, MPI_MINLOC, mpicomm); // find minimum among all procs
#endif
    double gamma_tilde = s.val; // all processors now have lasso angle
    blockId.first = s.rank;     // block on this processor if lasso angle is smallest
    if(s.rank == myrank) { 
      blockId.second.sub = ik; 
      blockId.second.index = jk[ik];
    }
    double gamma1; 
    while(true) { // loop to ensure linear independence
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
      for(int i=0; i<nsub; ++i) {
         for(std::list<std::pair<int,long_int> >::iterator it = nld_indices.begin(); it != nld_indices.end(); ++it)    // max out linear dependency
           if(it->first == myrank && it->second.sub == i) minbuffer[i][it->second.index] = std::numeric_limits<double>::max();
         gmax[i] = minbuffer[i].minCoeff(&jk[i]); 
      }

      s.val = (nsub > 0) ? gmax.minCoeff(&ik) : std::numeric_limits<double>::max();
      s.rank = myrank;
#ifdef USE_MPI
      MPI_Allreduce(MPI_IN_PLACE, &s, 1, MPI_DOUBLE_INT, MPI_MINLOC, mpicomm); // find minimum among all procs
#endif
      gamma1 = (s.val < gamma_tilde) ? s.val : gamma_tilde;
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
      for(int i=0; i<nsub; ++i) { // update solution
        y[i].head(l[i]).array() += gamma1*oneNwA*wA[i].head(l[i]).array();
      }

      if(s.val > gamma_tilde){ // reject an element if lasso angle less than lars angle
        dropId = true;
        gamma1 = gamma_tilde;
        break;
      } else {                 // add an element if lars angle less than lasso angle
        // exchange columns
        // find mpi process with minimum inactive set 
        q.index = 0; q.sub = myrank;
        for(int i=0; i<nsub; ++i) q.index += l[i];
#ifdef USE_MPI
        MPI_Allreduce(MPI_IN_PLACE, &q, 1, MPI_LONG_INT, MPI_MINLOC, mpicomm); 

        if(q.sub != s.rank) {
          if(s.rank == myrank) { // exchange column of A with process with fewest elements
            MPI_Sendrecv_replace(const_cast<double*>(A[ik].col(jk[ik]).data()), m, MPI_DOUBLE,
                                 q.sub, 0, q.sub, 0, mpicomm, &status);                                              // send column of A
            if(scaling) MPI_Sendrecv_replace(&S[ik][jk[ik]], 1, MPI_DOUBLE, q.sub, 1, q.sub, 1, mpicomm, &status);   // send scaling factor
            MPI_Sendrecv_replace(&perm[ik][jk[ik]].first, 1, MPI_INT, q.sub, 2, q.sub, 2, mpicomm, &status);         // send permutation 1
            MPI_Sendrecv_replace(&perm[ik][jk[ik]].second, 1, MPI_LONG_INT, q.sub, 3, q.sub, 3, mpicomm, &status);   // send permutation 2
          }
          if(q.sub == myrank) {
            MPI_Sendrecv_replace(const_cast<double*>(A[ik].col(jk[ik]).data()), m, MPI_DOUBLE,
                                 s.rank, 0, s.rank, 0, mpicomm, &status);                                            // recieve column A 
            if(scaling) MPI_Sendrecv_replace(&S[ik][jk[ik]], 1, MPI_DOUBLE, s.rank, 1, s.rank, 1, mpicomm, &status); // recieve scaling factor
            MPI_Sendrecv_replace(&perm[ik][jk[ik]].first, 1, MPI_INT, s.rank, 2, s.rank, 2, mpicomm, &status);       // recieve permuation 1
            MPI_Sendrecv_replace(&perm[ik][jk[ik]].second, 1, MPI_LONG_INT, s.rank, 3, s.rank, 3, mpicomm, &status); // recieve permutation 2
            p.index = jk[ik]; //broadcast new location
            p.sub   = ik;
          }
          s.rank = q.sub;
        }
#endif     
          
        // add new column and update gradient vector
        if(s.rank == myrank) {
          indices[ik].push_back(jk[ik]);
          B[ik].col(l[ik]) = S[ik][jk[ik]]*A[ik].col(jk[ik]);
#ifdef USE_MPI
          MPI_Scatterv(B[ik].col(l[ik]).data(), counts, displs, MPI_DOUBLE, buf3.data(), counts[myrank], MPI_DOUBLE, s.rank, mpicomm);
#else
          buf3 = B[ik].col(l[ik]);
#endif
          p.index = jk[ik];
          p.sub   = ik;
        }
#ifdef USE_MPI
        else {
          MPI_Scatterv(0, counts, displs, MPI_DOUBLE, buf3.data(), counts[myrank], MPI_DOUBLE, s.rank, mpicomm); // vector is now scattered
        }
        MPI_Bcast(&p, 1, MPI_LONG_INT, s.rank, mpicomm);
#endif

        gindices.push_back(std::pair<int,long_int>(s.rank,p)); // udpate global indices set
        grad = buf3.dot(update.segment(displs[myrank],counts[myrank])); 
#ifdef USE_MPI
        MPI_Allreduce(MPI_IN_PLACE, &grad, 1, MPI_DOUBLE, MPI_SUM, mpicomm);
#endif        
        grad = 1.0 - grad/oneNwA;                                // all processors now have grad
        buf3 *= grad*-1; 
        
        buf2.head(k+1) = BD.leftCols(k+1).transpose()*buf3;
  
#ifdef USE_MPI
        MPI_Allreduce(MPI_IN_PLACE, buf2.data(), maxvec, MPI_DOUBLE, MPI_SUM, mpicomm); // all processor now have buf2 to compute DtGDinv, it will be reused later
#endif

#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
        for(int i=0; i<nsub; ++i) { // compute new conjugate direction
          D[i].col(k+1).setZero();
          D[i].col(k+1).head(l[i]) = D[i].topLeftCorner(l[i],k+1).triangularView<Upper>()*buf2.head(k+1);
        }

        if(s.rank == myrank) {
          D[ik].col(k+1)[l[ik]] = grad; // set last element of new vector
          l[ik]++;
        }

#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
        for(int i=0; i<nsub; ++i) {
          c_loc.col(i) = B[i].leftCols(l[i])*D[i].col(k+1).head(l[i]);
        }
        buf1 = c_loc.rowwise().sum();
#ifdef USE_MPI
        MPI_Allreduce(MPI_IN_PLACE, buf1.data(), m, MPI_DOUBLE, MPI_SUM, mpicomm); // all processor now have buf1 to compute DtGDinv, it will be reused later
#endif
        DtGDinv = 1/buf1.squaredNorm();                                            // all processors now have DtGDinv
        BD.col(k+1) = DtGDinv*buf1.segment(displs[myrank],counts[myrank]);         // BD is now updated on every processor

        a[k+1] = grad*grad*DtGDinv;                                                // all processors now have new step length

        // check for linear dependence
        if(a[k+1] != a[k+1] /* check for nan*/ || a[k+1] <= std::numeric_limits<double>::min()  /* or too close to current columns*/ || grad <= 0.0) {
          nld_indices.push_back(std::pair<int,long_int>(s.rank,p));
          gindices.pop_back();
          if(s.rank == myrank) {
            indices[ik].pop_back();
            l[ik]--; 
          }
          continue;
        } else nld_indices.clear();

        dropId = false;

        break; // break from linear dependence loop 
      } 
    } // end linear independence loop
    // update solution estimates
    C -= gamma1*oneNwA; // all processes have C 
    r -= gamma1*update; // all processes have r (the residual)
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) { // compute new correlation vector 
      g[i] = S[i].asDiagonal()*(A[i].transpose()*r);
    }
 
    k++;
    iter++;

    if(dropId) { // if gamma_tilde is selected, remove that index and downdate
      dtime -= getTime(); 
      downIt++;

      // remove selected index from global and local indice sets
        if(blockId.first == myrank) {
          p.index = indices[blockId.second.sub][blockId.second.index];
          p.sub = blockId.second.sub;
          indices[blockId.second.sub].erase(indices[blockId.second.sub].begin() + blockId.second.index); // remove index jk[ik] from the local active set of ik-th subdomain
          for(int j=blockId.second.index; j<l[blockId.second.sub]-1; ++j) { y[blockId.second.sub][j] = y[blockId.second.sub][j+1]; } // erase jk[ik]-th element from y[ik]
          y[blockId.second.sub][l[blockId.second.sub]-1] = 0;
          l[blockId.second.sub]--;
          blockId.second.index = p.index; // since the index was removed from set, set block index now
        }      

#ifdef USE_MPI
        MPI_Bcast(&p, 1, MPI_LONG_INT, blockId.first, mpicomm);
#endif
        // remove index from the global active set
        std::list<std::pair<int,long_int> >::iterator pos = std::find(gindices.begin(), gindices.end(), std::pair<int,long_int>(blockId.first,p));
        std::list<std::pair<int,long_int> >::iterator fol = gindices.erase(pos);

        // Note: it is necessary to re-G-orthogonalize the basis D now, project the solution y onto the new basis and compute the corresponding residual r.
        // This is done here by starting from the column of D pointed to by fol (because the ones before this are already G-orthogonal), and then
        // following what is the essentially same procedure that is used above to construct the original basis, with a few optimizations when possible.
        k = std::distance(gindices.begin(), fol);

        // update wA accross all processors
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
        for(int i=0; i<nsub; ++i) {
           wA[i].setZero();
           l[i] = 0; for(std::list<std::pair<int,long_int> >::iterator it = gindices.begin(); it != fol; ++it) { if(it->first == myrank && it->second.sub == i) l[i]++; }
           wA[i].head(l[i]) = D[i].topLeftCorner(l[i],k).triangularView<Upper>()*a.head(k);
        } 

        // reconjugate all vectors after dropped element
        for(std::list<std::pair<int,long_int> >::iterator it = fol; it != gindices.end(); ++it,k++) {

#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
          for(int i=0; i<nsub; ++i) {// compute update vector
             c_loc.col(i) = B[i].leftCols(l[i])*wA[i].head(l[i]);
          }
          update = c_loc.rowwise().sum();
#ifdef USE_MPI
          MPI_Allreduce(MPI_IN_PLACE, update.data(), m, MPI_DOUBLE, MPI_SUM, mpicomm); // all processor now have update, it will be reused later
#endif

          // add new column and update gradient vector
          if(it->first == myrank) {
            int i = it->second.sub;
            B[i].col(l[i]) = S[i][it->second.index]*A[i].col(it->second.index);
#ifdef USE_MPI
            MPI_Scatterv(B[i].col(l[i]).data(), counts, displs, MPI_DOUBLE, buf3.data(), counts[myrank], MPI_DOUBLE, it->first, mpicomm);
#else
            buf3 = B[i].col(l[i]);
#endif
          }
#ifdef USE_MPI
          else {
             MPI_Scatterv(0, counts, displs, MPI_DOUBLE, buf3.data(), counts[myrank], MPI_DOUBLE, it->first, mpicomm); // vector is now scattered
          }
#endif

          grad = buf3.dot(update.segment(displs[myrank],counts[myrank])); 

#ifdef USE_MPI
          MPI_Allreduce(MPI_IN_PLACE, &grad, 1, MPI_DOUBLE, MPI_SUM, mpicomm);
#endif
          
          grad = 1.0 - grad; // all processors now have grad 
          buf3 *= grad*-1;    

          buf2.head(k) = BD.leftCols(k).transpose()*buf3;

#ifdef USE_MPI
          MPI_Allreduce(MPI_IN_PLACE, buf2.data(), maxvec, MPI_DOUBLE, MPI_SUM, mpicomm); // all processor now have buf2 to compute DtGDinv, it will be reused later
#endif

#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
          for(int i=0; i<nsub; ++i) {
            D[i].col(k).setZero();
            D[i].col(k).head(l[i]) = D[i].topLeftCorner(l[i],k).triangularView<Upper>()*buf2.head(k);
          }

          if(it->first == myrank) {
            int i = it->second.sub;
            D[i].col(k)[l[i]] = grad; // set last element of new vector
            l[i]++;
          }

#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
          for(int i=0; i<nsub; ++i) {
            c_loc.col(i) = B[i].leftCols(l[i])*D[i].col(k).head(l[i]);
          }
          buf1 = c_loc.rowwise().sum();
#if USE_MPI
          MPI_Allreduce(MPI_IN_PLACE, buf1.data(), m, MPI_DOUBLE, MPI_SUM, mpicomm); // all processor now have buf1 to compute DtGDinv, it will be reused later
#endif
          DtGDinv = 1/buf1.squaredNorm(); // all processors now have DtGDinv 
          BD.col(k) = DtGDinv*buf1.segment(displs[myrank],counts[myrank]); // BD is now updated on every processor

          a[k] = grad*grad*DtGDinv; // all processors now have the new step length

#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
          for(int i=0; i<nsub; ++i) {
             wA[i].head(l[i]) += a[k]*D[i].col(k).head(l[i]);
          }

      }
      k--;
      dtime += getTime();
    }// end if 

  }// end while

  delete [] counts;
  delete [] displs;

  dtime /= 1000.0;
  if(myrank == 0 && verbose) std::cout.flush();

  Array<VectorXd,Dynamic,1> x(nsub);
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
  for(int i=0; i<nsub; ++i)
    x[i] = VectorXd::Zero(A[i].cols());

  if(project){

    if(myrank == 0)
      std::cout << "Performing Non-negative least squares projection" << std::endl;

    Array<VectorXd,Dynamic,1> dummyx(nsub);
    std::vector<Eigen::Map<Eigen::MatrixXd> > subsetA(nsub, Eigen::Map<Eigen::MatrixXd>(NULL,0,0));
    for(int i=0; i<nsub; ++i) {// pass only the selected subset of elements to the parallel NNLS solver
      new (&subsetA[i]) Eigen::Map<Eigen::MatrixXd>(B[i].data(),B[i].rows(),l[i]);
    }
#if defined(USE_MPI) && defined(USE_SCALAPACK)
    dummyx = splh(subsetA, b, rnorm, n, info, maxsze, maxit, reltol, true, scaling, false, dtime, 0, 0, 0, 0, 0, 0);
#else
    std::list<std::pair<int,long_int> > hotIndices;
    dummyx = pnncgp(subsetA, b, rnorm, n, info, maxsze, maxite, reltol, true, scaling, center, dtime, hotIndices);
#endif
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i)
      for(long int j=0; j<l[i]; ++j) x[i][indices[i][j]] = S[i][indices[i][j]]*dummyx[i][j];

  } else {

#if defined(_OPENMP)
  #pragma omp parallel for schedule(static,1)
#endif
    for(int i=0; i<nsub; ++i) 
      for(long int j=0; j<l[i]; ++j) x[i][indices[i][j]] = S[i][indices[i][j]]*y[i][j];

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
