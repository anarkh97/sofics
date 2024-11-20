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

Eigen::VectorXd
nncgp(Eigen::Ref<Eigen::MatrixXd> A, Eigen::Ref<Eigen::VectorXd> b, double& rnorm,
      long int &info, double maxsze, int &maxEle, double maxite, double reltol, bool verbose, bool scaling, bool center, bool reverse, double &dtime, std::vector<long int> &indices);

Eigen::VectorXd
cglars(Eigen::Ref< Eigen::MatrixXd> A, Eigen::Ref< Eigen::VectorXd> b, double& rnorm,
      long int &info, double maxsze, int &maxEle, double maxite, double reltol, bool verbose, bool scaling, bool center, bool project, bool positive, double &dtime)
{
  using namespace Eigen;

  // problem dimensions
  const long int m = b.rows();  
  const long int maxvec = (maxEle > 0) ? std::min(m,(long int)maxEle)+1 : std::min(m, (long int)(maxsze*A.cols()))+1;
  const long int maxit = maxite*A.cols();

  // allocate 
  VectorXd wA(maxvec), ylar(maxvec), t(maxvec);
  VectorXd crlt(A.cols()), h(A.cols()), S(A.cols()), minBuffer(A.cols());
  VectorXd residual(A.rows()), update(A.rows()), a(maxvec);

  MatrixXd B(A.rows(),maxvec), D(maxvec,maxvec);
  Matrix<double,Dynamic,Dynamic,ColMajor> BD(A.rows(),maxvec);
  double grad, DtGDinv;
  std::vector<long int> indices;     // inactive set
  std::vector<long int> nld_indices; // linearly dependent set
 
  // initialize
  wA.setZero(); ylar.setZero(); t.setZero();
  minBuffer.setZero();
  B.setZero(); D.setZero(); BD.setZero(); 

  info = 1;

  // center covariates and target so that the mean is 0
  if(center){
   std::cout << "Centering Covariates" << std::endl;
   for(int i=0; i<A.cols(); ++i) { A.col(i).array() -= A.col(i).mean();}
   b.array() -= b.mean();
  }

  // scale columns to norm 1
  if(scaling) for(int i=0; i<A.cols(); ++i) { double s = A.col(i).norm(); S[i] = (s != 0) ? 1/s : 0; }
  else S.setOnes();

  dtime = 0;
  residual = b;
  double bnorm  = b.norm();
  double abstol = reltol*bnorm;
  double C;  

  // ensure that an element is not selected immediately after being dropped
  bool     dropId  = false;
  long int blockId = 0;

  long int k       = 0; // current column index <-> Active set size -1
  long int iter    = 0; // number of iterations
  long int downIt  = 0; // number of downdates


  // initialize starting set with maximaly corellated element
  {
    crlt = S.asDiagonal()*(A.transpose()*b);
    long int i;
    C = crlt.maxCoeff(&i);
    indices.push_back(i);

    // initialize CG
    B.col(k)       = S[indices[0]]*(A.col(indices[0]));
    Block<MatrixXd,Dynamic,1,true> d = D.col(k), c = BD.col(k); 
    grad           = 1.0;
    d[k]           = grad;
    c              = B.leftCols(k+1)*d.head(k+1);
    DtGDinv        = 1/c.squaredNorm();
    c             *= DtGDinv;
    a[k]           = grad*grad*DtGDinv;
  }

  while(true) {

    rnorm = residual.norm();

    if(verbose) {
      std::cout << "Iteration = " << std::setw(9) << iter << "    "
                << "Downdate = " << std::setw(9) << downIt << "    "
                << "Active set size = " << std::setw(9) << k << "    "
                << "Max Correlation = " << std::setw(13) << std::scientific << std::uppercase << C << "    "
                << "Residual norm = " << std::setw(13) << std::scientific << std::uppercase << rnorm << "    "
                << "Target = " << std::setw(13) << std::scientific << std::uppercase << abstol << std::endl;
      std::cout.unsetf(std::ios::scientific);
      std::cout.unsetf(std::ios::uppercase);
    }

    if((rnorm <= abstol && maxEle == 0) || k+nld_indices.size()+1 == maxvec) {maxEle = k; break;}
    if(iter >= maxit) { info = 3; break; }

    long int i = 0; // dummy integer for stuff

    if(!dropId){// don't update solution if we just dropped an inactive element
      // update the solution to the sytem G^T*G*w = 1
      Block<MatrixXd,Dynamic,1,true> d = D.col(k);
      wA.head(k+1) = wA.head(k+1) + a[k]*d.head(k+1);
    }

    if(!dropId && wA[k] <= 0) {
       std::cout << "wA = " << wA.head(k+1).transpose() << std::endl;
       std::cout << "C = " << std::endl;
       std::cout << "DtGDinv = " << DtGDinv << std::endl;
       std::cout << "a = " << a[k] << std::endl;
       for(std::vector<long int>::iterator it = indices.begin(); it != indices.end(); ++it) std::cout << crlt[*it] <<  " ";
       std::cout << "\n *** Roundoff Error *** " << std::endl;
       exit(-1);
       break;
    }
    double oneNwA = 1.0/sqrt(wA.head(k+1).sum());

    // compute smallest angle at which a new covariant becomes dominant
    update = oneNwA*B.leftCols(k+1)*wA.head(k+1);
    h      = S.asDiagonal()*(A.transpose()*update);

    double gamma1, gamma_tilde;

    //compute LARS angle
    minBuffer = (C - crlt.array())/(oneNwA - h.array());
    for(std::vector<long int>::iterator it = indices.begin();     it != indices.end();     ++it) minBuffer[*it] = std::numeric_limits<double>::max();  // max out inactive set
    for(long int row = 0; row < minBuffer.rows(); ++row) minBuffer[row] = (minBuffer[row] <= 0) ? std::numeric_limits<double>::max() : minBuffer[row]; // max out non-positive members
    if(dropId) minBuffer[blockId] = std::numeric_limits<double>::max();                                                                                // max out previously rejected member 

    // compute smallest angle at which ylars changes sign
    long int j;
    t.setZero();
    t.head(k+1) = -1.0*ylar.head(k+1).array()/(oneNwA*wA.head(k+1).array());
    for(long int ele = 0; ele < k+1; ++ele) t[ele] = (t[ele] <= 0.) ? std::numeric_limits<double>::max() : t[ele];
    gamma_tilde = (k > 0) ? t.head(k+1).minCoeff(&j) : std::numeric_limits<double>::max(); // skip this step on first iteration

    while(true) { // loop to ensure linear independence
      for(std::vector<long int>::iterator it = nld_indices.begin(); it != nld_indices.end(); ++it) minBuffer[*it] = std::numeric_limits<double>::max();  // max out linearly dependent members
      gamma1 = minBuffer.minCoeff(&i); // compute step length 

      if(gamma_tilde < gamma1){
        dropId = true;
        gamma1 = gamma_tilde;
        i = j; // drop index if gamma_tilde selected

        blockId = indices[i];
        break; // break from linear dependence loop
      } else {
        indices.push_back(i); // add index if gamma1 selected

        B.col(k+1) = S[indices[k+1]]*(A.col(indices[k+1]));

        // update BD due to extra column added to B (note: B.col(i)*D.row(i).head(k) = 0, so BD does not need to be updated)
        grad = oneNwA - (B.col(k+1).transpose()*update);
        grad = grad/oneNwA;
         
        Block<MatrixXd,Dynamic,1,true> d_next = D.col(k+1), c = BD.col(k+1); 
        d_next.head(k+1) = D.topLeftCorner(k+1,k+1).triangularView<Upper>()*(BD.leftCols(k+1).transpose()*B.col(k+1)*grad*-1);
        d_next[k+1] = grad; 
        c = B.leftCols(k+2)*d_next.head(k+2);
        DtGDinv = 1/c.squaredNorm();
        c *= DtGDinv; 
        a[k+1] = grad*grad*DtGDinv;

        //if diagonal element is too small, then column is near linearly dependent
        if(a[k+1] != a[k+1] /* check for nan*/ || a[k+1] <= std::numeric_limits<double>::min()  /* or too close to current columns*/ || grad <= 0. /* or already in subspace */) {
          nld_indices.push_back(i); indices.pop_back();
          std::cout << "*** Rejecting selected covariant [" << i << "] ***" << std::endl;
          continue;
        } else {
          nld_indices.clear();
        }

        dropId = false;
        break; // break from linear dependence loop
      }
    }
    C        -= gamma1*oneNwA;
    residual -= gamma1*update;
    // update solution and estimate
    ylar.head(k+1) = ylar.head(k+1).array() + gamma1*(oneNwA*wA.head(k+1).array());
    crlt = S.asDiagonal()*(A.transpose()*residual); // correlations

    k++;
    iter++;
    if(dropId) {// if gamma_tilde is selected, remove that index and downdate 
      dtime -= getTime();
      downIt++;

      k = i;

      //remove selected index
      std::vector<long int>::iterator fol = indices.erase(indices.begin()+k);

      //zero out row
      ylar[k] = 0;
      wA.setZero();
      wA.head(k) = D.topLeftCorner(k,k).triangularView<Upper>()*a.head(k);

      //reconjugate vectors
      for(std::vector<long int>::iterator it = fol; it != indices.end(); ++it,k++){
        B.col(k) = B.col(k+1);
        // set gradient
        grad = 1 - B.col(k).transpose()*B.leftCols(k)*wA.head(k);
        // reconjugate vector
        Block<MatrixXd,Dynamic,1,true> d = D.col(k), c = BD.col(k);
        d.head(k) = D.topLeftCorner(k,k).triangularView<Upper>()*(BD.leftCols(k).transpose()*B.col(k)*grad*-1);
        d[k] = grad;
        c = B.leftCols(k+1)*d.head(k+1);
        DtGDinv = 1/c.squaredNorm();
        c *= DtGDinv;
        // update step length
        a[k] = grad*grad*DtGDinv;
        wA.head(k+1) = wA.head(k+1) + a[k]*d.head(k+1);
        // transfer solution 
        ylar[k] = ylar[k+1];
      }
      B.col(k).setZero();
      ylar[k] = 0;
      k--;
      dtime += getTime();
    }

  }// END LASSO loop

  if(project) { // can do a non-negative least squares projection onto the non-negative lasso basis
    std::cout << "*** PROJECTING SOLUTION ON TO SELECTED BASIS ***" << std::endl;
    maxEle = 0;
    std::vector<long int> newIndices; 
    ylar.head(k) = nncgp(B.leftCols(k), b, rnorm, info, maxsze, maxEle, maxite, reltol, verbose, scaling, center, false, dtime, newIndices);
  }


  residual = b - B.leftCols(k)*ylar.head(k);
  std::cout << "Projected Residual = " << residual.norm() <<  std::endl;

  dtime /= 1000.0;
  if(verbose) std::cout.flush();

  VectorXd x = VectorXd::Zero(A.cols());
  for(long int j=0; j<k; ++j) {
     x[indices[j]] = S[indices[j]]*ylar[j];}
  return x;

}

#endif
