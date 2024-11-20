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

// Sequential implementation of Non-negative Conjugate Gradient Pursuit.
// This is a non-negative constrained variant of Conjugate Gradient Pursuit, and generates identical iterates to
// Lawson & Hanson's NNLS (non-negative least squares) algorithm.
// References:
// 1. Blumensath, Thomas, and Michael E. Davies. "Gradient pursuits." Signal Processing, IEEE Transactions on 56.6 (2008): 2370-2382.
// 2. Lawson, C. L., & Hanson, R. J. (1974). Solving least squares problems (Vol. 161). Englewood Cliffs, NJ: Prentice-hall.

Eigen::VectorXd
nncgp(Eigen::Ref<Eigen::MatrixXd> A, Eigen::Ref<Eigen::VectorXd> b, double& rnorm,
      long int &info, double maxsze, int &maxEle, double maxite, double reltol, bool verbose, bool scaling, bool center, bool reverse, double &dtime, std::vector<long int> &hotIndices)
{
  using namespace Eigen;

  std::vector<long int> indices;

  // if doing a hot start, pass in indice set of pre-selected elements and compute weights for new target
  bool hotStart = false;

  const long int m = b.rows();
  const long int maxvec = (maxEle > 0) ? std::min(m,(long int)maxEle) : std::min(m, (long int)(maxsze*A.cols()));
  const long int maxit = maxite*A.cols();
  double bnorm = b.norm();
  double abstol = reltol*bnorm;

  VectorXd x_(maxvec), y(maxvec), r(A.rows()), g(A.cols()), h(A.cols()), DtGDinv(maxvec), g_(maxvec), a(maxvec), S(A.cols()), t(maxvec);
  x_.setZero();
  r = b;
  rnorm = bnorm;
  info = 1;
  MatrixXd B(A.rows(),maxvec), D(maxvec,maxvec), GD(maxvec,maxvec);
  Matrix<double,Dynamic,Dynamic,ColMajor> BD(A.rows(),maxvec);
  B.setZero(); D.setZero(); BD.setZero(); GD.setZero();
  long int k = 0;
  int sizeHS = 0;
  if(hotIndices.size() > 0) {
    hotStart = true;
    sizeHS = hotIndices.size();
  }
    
  std::vector<long int> nld_indices;

  //center covariates and target so that the mean is 0
  if(center){
   std::cout << "Centering Covariates" << std::endl;
   for(int i=0; i<A.cols(); ++i) { A.col(i).array() -= A.col(i).mean();}
   b.array() -= b.mean();
  }

  if(scaling) for(int i=0; i<A.cols(); ++i) { double s = A.col(i).norm(); S[i] = (s != 0) ? 1/s : 0; }
  else S.setOnes();

  dtime      = 0; // time spent in downdates
  int iter   = 0; // number of iterations
  int downIt = 0; // number of downdates

  if(hotStart) {
    std::cout << " ... Computing Hot Start Quantities ... " << std::endl;
  }

  while(true) {   // primary loop
    long int i; 
    if(verbose) {
      std::cout << "Iteration = " << std::setw(9) << iter << "    "
                << "Downdate = " << std::setw(9) << downIt << "    "
                << "Active set size = " << std::setw(9) << k << "    "
                << "Residual norm = " << std::setw(13) << std::scientific << std::uppercase << rnorm << "    "
                << "Target = " << std::setw(13) << std::scientific << std::uppercase << abstol << std::endl;
      std::cout.unsetf(std::ios::scientific);
      std::cout.unsetf(std::ios::uppercase);
    }

    if((rnorm <= abstol && maxEle == 0) || k+nld_indices.size() == maxvec ) {maxEle = k; break;}
    if(iter >= maxit) { maxEle = k; info = 3; break; }

    if(hotStart) { // if starting with an initial element set, onl multiply appropriate column
      g.setZero();
      for(int prek = 0; prek < k+1; ++prek){
        g(hotIndices[prek]) = S(hotIndices[prek])*A.col(hotIndices[prek]).transpose()*r;
      }
      i  = hotIndices[k]; 
      if(k >= sizeHS){
        hotStart = false;
        std::cout << " ... Hot Start Complete ... " << std::endl;
      }
    } else {
      g = S.asDiagonal()*(A.transpose()*r); // gradient
      h = g; for(long int j=0; j<k; ++j) h[indices[j]] = -std::numeric_limits<double>::max(); // make sure the index has not already been selected
      for(long int j=0; j<nld_indices.size(); ++j) h[nld_indices[j]] = -std::numeric_limits<double>::max(); // also make sure near linear dependent indices are not selected
      double gi = h.maxCoeff(&i);
      if(gi <= 0) {maxEle = k; break;}
    }
    B.col(k) = S[i]*A.col(i);
    indices.push_back(i);
    // update GD due to extra column added to B (note: B.col(i)*D.row(i).head(k) = 0, so BD does not need to be updated)
    GD.row(k).head(k) = B.col(k).transpose()*BD.leftCols(k);
    for(long int j=0; j<k+1; ++j) g_[j] = g[indices[j]];
    Block<MatrixXd,Dynamic,1,true> d = D.col(k), c = BD.col(k);
    d.head(k+1) = g_.head(k+1) - D.topLeftCorner(k+1,k).triangularView<Upper>()*(DtGDinv.head(k).asDiagonal()*(GD.topLeftCorner(k+1,k).transpose()*g_.head(k+1))); // direction

    c = B.leftCols(k+1)*d.head(k+1);
    GD.col(k).head(k+1) = B.leftCols(k+1).transpose()*c;
    DtGDinv[k] = 1/c.squaredNorm();
    a[k] = r.dot(c)*DtGDinv[k]; // step length
    if(a[k] < 0) { nld_indices.push_back(i); indices.pop_back(); continue; } else nld_indices.clear(); // check for near linear dependence
    y.head(k+1) = x_.head(k+1) + a[k]*d.head(k+1); // candidate solution
    r -= a[k]*c; // residual
    k++;
    while(true) { // secondary loop
      iter++;
      if(y.head(k).minCoeff() < 0) { // check if constraints are violated
        std::vector<long int>::iterator fol;
        dtime -= getTime();  
        downIt++;
        // compute maximum feasible step length in the direction (y-x_) and corresponding index in the active set, i
        for(long int j=0; j<k; ++j) t[j] = (y[j] >= 0) ? std::numeric_limits<double>::max() : -x_[j]/(y[j]-x_[j]);
        double alpha = t.head(k).minCoeff(&i);

        // remove index i from the active set
        fol = indices.erase(indices.begin()+i);

        if(hotStart) {
          sizeHS--;
          hotIndices.erase(hotIndices.begin()+i);
        }

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
        if(k < indices.size()/2)
          r = b - BD.leftCols(k)*a.head(k);
        else 
          r = r + BD.block(0,i,m,indices.size()-i+1)*a.segment(i,indices.size()-i+1);
 
        g_.head(k) = B.leftCols(k).transpose()*r;

        std::vector<long int> indicesMap;
        MatrixXd GDcopy;
        if(reverse){//experimental: shuffle columns to try and get better performance
          int row = 0; 
          // generate random reordering
          for(std::vector<long int>::iterator it = fol; it != indices.end(); ++it, ++row){ 
             indicesMap.push_back(row); 
          }          

          // make copy of GD so we don't have to recompute
          GDcopy.resize(indicesMap.size(),i);
          // reorder indices vector
          std::vector<long int> indcopy(fol,indices.end());
          std::reverse(indicesMap.begin(),indicesMap.end());
          int dummy = 0;
          for(std::vector<long int>::iterator it = indicesMap.begin(); it != indicesMap.end(); ++it, ++dummy) {
            indices[i+dummy] = indcopy[*it];
          }
          //copy GD block into container
          GDcopy = GD.block(k+1,0,indicesMap.size(),i);
        }

        std::vector<long int>::iterator itMap = indicesMap.begin();
        for(std::vector<long int>::iterator it = fol; it != indices.end(); ++it, ++itMap) {
          B.col(k) = S[*it]*A.col(*it);
          g_[k] = B.col(k).transpose()*r;
          if(reverse){
            GD.row(k).head(i) = GDcopy.row(*itMap).head(i);
          } else {
            GD.row(k).head(i) = GD.row(k+1).head(i);
          }
          GD.row(k).segment(i,k-i) = B.col(k).transpose()*BD.block(0,i,m,k-i);
          Block<MatrixXd,Dynamic,1,true> d = D.col(k), c = BD.col(k); 
          d.head(k+1) = g_.head(k+1) - D.topLeftCorner(k+1,k).triangularView<Upper>()*(DtGDinv.head(k).asDiagonal()*(GD.topLeftCorner(k+1,k).transpose()*g_.head(k+1)));
          c = B.leftCols(k+1)*d.head(k+1);
          GD.col(k).head(k+1) = B.leftCols(k+1).transpose()*c;
          DtGDinv[k] = 1/c.squaredNorm();
          a[k] = r.dot(c)*DtGDinv[k];
          y.head(k+1) += a[k]*d.head(k+1);
          r -= a[k]*c;
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

  dtime /= 1000.0;
  if(verbose) std::cout.flush();

  hotIndices.clear();
  hotIndices = indices;

  VectorXd x = VectorXd::Zero(A.cols());
  for(long int j=0; j<k; ++j) x[indices[j]] = S[indices[j]]*x_[j];
  return x;
}

#endif
