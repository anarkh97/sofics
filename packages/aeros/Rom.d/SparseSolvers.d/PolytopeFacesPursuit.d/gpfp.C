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

// This function solves the Basis Pursuit problem
//
// min ||x||_0 s.t. A*x = b
//
// which is equivalent to 
//
// min ||s||_1 s.t. [A,-A]*s = b & s >= 0 
//
// and its dual problem
//
// min b^T*c s.t. [A,-A]^T*c <= 1
//
// where s_i = x_i > 0 for i = 1:N
// and   s_i = x_i < 0 for i = N+1:2*N
//
// The solution framework is known as Basis Pursuit and the algorithm is called Gradient Polytope 
// Faces Pursuit. This is a particular implementation of PFP which uses the Conjugate Gradient method to 
// solve a system of equations at each iteration whose dimensionality changes from one iteration to the 
// next, thus only matrix vector products are required and the scheme is easily parallelizable.
// This algorithm produces identical residuals to that of the standard Polytope Faces Pursuit
//
// References:
// 1. Plumbley, M & Gretsistas, I. "Gradient Polytope Faces Pursuit for Large Sparse Recovery Problems" ICASSP
// 2. Blumensath, T & Davies, M. "Gradient Pursuits" IEEE
//
// ARGUMENTS:
// A        = where A*x= b
// b        = target vector 
// rnorm    = residual norm 
// maxsze   = maximum allowable sparsity
// maxite   = maximum allowable iterations
// reltol   = stopping criteria
// verbose  = verbose flag
// scaling  = flag to turn on/off unit normalization of columns of A
// positive = turn positivity constraint on/off

Eigen::VectorXd
nncgp(Eigen::Ref<Eigen::MatrixXd> A, Eigen::Ref<Eigen::VectorXd> b, double& rnorm,
      long int &info, double maxsze, int &maxEle, double maxite, double reltol, bool verbose, bool scaling, bool center, bool reverse, double &dtime, std::vector<long int> &indices);

Eigen::VectorXd
gpfp(Eigen::Ref<Eigen::MatrixXd> A, Eigen::Ref<Eigen::VectorXd> b, double& rnorm,
     long int &info, double maxsze, int &maxEle, double maxite, double reltol, bool verbose, bool scaling, bool center, bool project, bool positive, double &dtime)
{
  using namespace Eigen;

  const long int m = b.rows();
  const long int maxvec = (maxEle > 0) ? std::min(m,(long int)maxEle) : std::min(m, (long int)(maxsze*A.cols()));
  const long int maxit = maxite*A.cols();
  double bnorm = b.norm();
  double abstol = reltol*bnorm;

  VectorXd x_(maxvec), y(maxvec), r(A.rows()), vertex(A.rows()), g1(A.cols()), g2(A.cols()), h(A.cols()), DtGDinv(maxvec), g_(maxvec), lambda(maxvec), a(maxvec), S(A.cols()), t(maxvec);
  x_.setZero();
  y.setZero();
  vertex.setZero();
  r = b;
  rnorm = bnorm;
  info = 1;
  MatrixXd B(A.rows(),maxvec), D(maxvec,maxvec), GD(maxvec,maxvec);
  Matrix<double,Dynamic,Dynamic,ColMajor> BD(A.rows(),maxvec);
  B.setZero(); D.setZero(); BD.setZero(); GD.setZero();
  long int k = 0;
  std::vector<long int> indices;
  std::vector<int>      setKey;
  std::vector<long int> nld_indices;
  std::vector<long int> nld_setKey;

  if(center){
   std::cout << "Centering Covariates" << std::endl;
   for(int i=0; i<A.cols(); ++i) { A.col(i).array() -= A.col(i).mean();}
   b.array() -= b.mean();
  }

  if(scaling) for(int i=0; i<A.cols(); ++i) S[i] = 1/A.col(i).norm();
  else S.setOnes();

  dtime = 0;
  int iter   = 0; // number of iterations
  int downIt = 0; // number of downdates
  while(true) {

    if(verbose) {
      std::cout << "Iteration = " << std::setw(9) << iter << "    "
                << "Downdate = " << std::setw(9) << downIt << "    "
                << "Active set size = " << std::setw(9) << k << "    "
                << "Residual norm = " << std::setw(13) << std::scientific << std::uppercase << rnorm << "    "
                << "Target = " << std::setw(13) << std::scientific << std::uppercase << abstol << std::endl;
      std::cout.unsetf(std::ios::scientific);
      std::cout.unsetf(std::ios::uppercase);
    }

    if((rnorm <= abstol && maxEle == 0) || k+nld_indices.size() == maxvec) {maxEle = k; break;}
    if(iter >= maxit) { info = 3; break; }

    g1.setZero();
    g2.setZero();
    int Set = 1;
    long int position  = 0;
    long int position1 = 0;
    long int position2 = 0;

    // compute step length along current face
    g1 = S.asDiagonal()*(A.transpose()*r);
    g2 = S.asDiagonal()*(A.transpose()*vertex);

    for(long int col = 0; col != A.cols(); col++) {
      double num = g1(col);
      double den = g2(col);
      if(num > 0. && den != 1.0) {
        h(col) = num/(1.0-den);
        if(!positive)
          g2(col) = -std::numeric_limits<double>::max();
      } else if(num < 0. && den != -1.0) {
        h(col) = -std::numeric_limits<double>::max();
        if(!positive)
          g2(col) = (-1.0)*num/(1.0+den);
      } else {
        h(col) = -std::numeric_limits<double>::max();
        if(!positive)
          g2(col) = -std::numeric_limits<double>::max();
      }
    }
    // make sure the index has not already been selected
    for(long int j=0; j<k; ++j) {
      if(setKey[j] == 1) {
     h(indices[j]) = -std::numeric_limits<double>::max();
      } else if(!positive){
        g2(indices[j]) = -std::numeric_limits<double>::max();
      }
    }

    // also make sure near linear dependent indices are not selected
    for(long int j=0; j<nld_indices.size(); ++j) {
      if(nld_setKey[j] == 1)
        h[nld_indices[j]] = -std::numeric_limits<double>::max();
      else if(!positive)
        g2[nld_indices[j]] = -std::numeric_limits<double>::max();
    }

    double lam  = 0.0;
    double lam1 = h.maxCoeff(&position1);
    double lam2 = 0.;
    if(!positive)
      lam2 = g2.maxCoeff(&position2);

    if ((lam1 > lam2) || positive) {
      lam = lam1;
      position = position1;
    } else {
      lam = lam2;
      position = position2;
      Set = -1;
    }

    // add step length for vertex estimate update to array
    lambda[k] = 1./lam;

    // remove maximum element from the correct active set and place
    // in correct dual set, also store appropriate column of [A,-A]
    setKey.push_back(Set);
    B.col(k) = double(Set)*S[position]*A.col(position);
    indices.push_back(position);

    // update GD due to extra column added to B (note: B.col(i)*D.row(i).head(k) = 0, so BD does not need to be updated)
    GD.row(k).head(k) = B.col(k).transpose()*BD.leftCols(k);

    for(long int j=0; j<k+1; ++j) g_[j] = double(setKey[j])*g1[indices[j]];
    Block<MatrixXd,Dynamic,1,true> d = D.col(k), c = BD.col(k);
    d.head(k+1) = g_.head(k+1) - D.topLeftCorner(k+1,k).triangularView<Upper>()*(DtGDinv.head(k).asDiagonal()*(GD.topLeftCorner(k+1,k).transpose()*g_.head(k+1))); // direction

    c = B.leftCols(k+1)*d.head(k+1);
    GD.col(k).head(k+1) = B.leftCols(k+1).transpose()*c;
    DtGDinv[k] = 1./c.squaredNorm();
    a[k] = r.dot(c)*DtGDinv[k]; // step length
    // check for near linear dependence
    if(a[k] < 0) { nld_indices.push_back(position); nld_setKey.push_back(Set); indices.pop_back(); setKey.pop_back(); continue; } else { nld_indices.clear(); nld_setKey.clear(); }
    y.head(k+1) = x_.head(k+1) + a[k]*d.head(k+1); // candidate solution
    r -= a[k]*c; // residual
    vertex += lambda[k]*a[k]*c;
    k++;
    while(true) {
      iter++;
      int i = 0;
      // find location of minimum coefficient and check for positivity
      double minCoeff = y.head(k).minCoeff(&i); 
      if(minCoeff < 0.) {
        dtime -= getTime();
        downIt++;

        // compute maximum feasible step length in the direction (y-x_) and corresponding index in the active set, i
        for(long int j=0; j<k; ++j) t[j] = (y[j] >= 0) ? std::numeric_limits<double>::max() : -x_[j]/(y[j]-x_[j]);
        double alpha = t.head(k).minCoeff(&i);
 
        // remove index i from the active set
        std::vector<long int>::iterator fol = indices.erase(indices.begin()+i);
        setKey.erase(setKey.begin()+i);

        // update x_ (note: this is used only when there are two or more consecutive downdate iterations)
        for(int j=0; j<i; ++j) x_[j] += alpha*(y[j]-x_[j]);
        for(int j=i; j<k-1; ++j) x_[j] = x_[j+1] + alpha*(y[j+1]-x_[j+1]); 
        x_[k-1] = 0;

        // Note: it is necessary to re-G-orthogonalize the basis D now, project the solution x_ onto the new basis and compute the corresponding residual r.
        // This is done here by starting from the column of D pointed to by fol (because the ones before this are already G-orthogonal), and then
        // following what is the essentially same procedure that is used above to construct the original basis, with a few optimizations when possible.
        y.segment(i,k-i).setZero();
        k = i;
        y.head(k) = D.topLeftCorner(k,k).triangularView<Upper>()*a.head(k);
        r = b - BD.leftCols(k)*a.head(k);
        vertex = BD.leftCols(k)*lambda.head(k).asDiagonal()*a.head(k);
        g_.head(k) = B.leftCols(k).transpose()*r;
        for(std::vector<long int>::iterator it = fol; it != indices.end(); ++it) {
          B.col(k) = double(setKey[k])*S[*it]*A.col(*it);
          g_[k] = B.col(k).transpose()*r;
          lambda[k] = (1.0-B.col(k).dot(vertex))/(g_[k]);
          GD.row(k).head(i) = GD.row(k+1).head(i);
          GD.row(k).segment(i,k-i) = B.col(k).transpose()*BD.block(0,i,m,k-i);
          Block<MatrixXd,Dynamic,1,true> d = D.col(k), c = BD.col(k); 
          d.head(k+1) = g_.head(k+1) - D.topLeftCorner(k+1,k).triangularView<Upper>()*(DtGDinv.head(k).asDiagonal()*(GD.topLeftCorner(k+1,k).transpose()*g_.head(k+1)));
          c = B.leftCols(k+1)*d.head(k+1);
          GD.col(k).head(k+1) = B.leftCols(k+1).transpose()*c;
          DtGDinv[k] = 1/c.squaredNorm();
          a[k] = r.dot(c)*DtGDinv[k];
          y.head(k+1) += a[k]*d.head(k+1);
          r -= a[k]*c;
          vertex += lambda[k]*a[k]*c;
          k++;
          g_.head(k) -= a[k-1]*GD.col(k-1).head(k);
        }
        dtime += getTime();
      }
      else {
        x_.head(k) = y.head(k);
        break;
      }
    }

    rnorm = r.norm();
  }

  if(project) { // can do a non-negative least squares projection onto the non-negative lasso basis
    std::cout << "*** PROJECTING SOLUTION ON TO SELECTED BASIS ***" << std::endl;
    maxEle = 0;
    std::vector<long int> newIndices;
    y.head(k) = nncgp(B.leftCols(k), b, rnorm, info, maxsze, maxEle, maxite, reltol, verbose, scaling, center, false, dtime, newIndices);
  }

  r = b - B.leftCols(k)*y.head(k);
  std::cout << "Projected Residual = " << r.norm() <<  std::endl;

  dtime /= 1000.0;
  if(verbose) std::cout.flush();

  VectorXd x = VectorXd::Zero(A.cols());
  long int element = 0;
  for(std::vector<int>::iterator it = setKey.begin(); it != setKey.end(); it++) {
    x[indices[element]] = S[indices[element]]*double(*it)*y[element];
    element++;
  }
  return x;
}

#endif
