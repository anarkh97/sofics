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

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
};

Eigen::VectorXd
nncgp(Eigen::Ref<Eigen::MatrixXd> A, Eigen::Ref<Eigen::VectorXd> b, double& rnorm,
      long int &info, double maxsze, int &maxEle, double maxite, double reltol, bool verbose, bool scaling, bool center, bool reverse, double &dtime, std::vector<long int> &indices);

Eigen::VectorXd
lars(Eigen::Ref< Eigen::MatrixXd> A, Eigen::Ref< Eigen::VectorXd> b, double& rnorm,
      long int &info, double maxsze, int &maxEle, double maxite, double reltol, bool verbose, bool scaling, bool center, bool project, bool positive, double &dtime)
{
  using namespace Eigen;
  
  const long int m = b.rows();  
  const long int maxvec = (maxEle > 0) ? std::min(m,(long int)maxEle)+1 : std::min(m, (long int)(maxsze*A.cols()))+1;
  const long int maxit = maxite*A.cols();

  // allocate
  VectorXd wA(maxvec), rhs(maxvec), vk(maxvec), ylar(maxvec), t(maxvec);
  VectorXd crlt(A.cols()), h(A.cols()), S(A.cols()), minBuffer(A.cols());
  VectorXd colK(A.rows()), residual(A.rows()), update(A.rows());
 
  // initialize
  wA.setZero(); rhs.setOnes(); vk.setZero(); ylar.setZero(); t.setZero(); 
  minBuffer.setZero();
  colK.setZero();  

  info = 1;

  MatrixXd B(A.rows(),maxvec), R(maxvec,maxvec);
  B.setZero(); R.setZero();//B is storage for inactive columns and R is cholesky factorization of Gramm matrix

  std::vector<long int> indices;
  std::vector<long int> nld_indices;

  //center covariates and target so that the mean is 0
  if(center){
     std::cout << "Centering Covariates" << std::endl;
     for(int i=0; i<A.cols(); ++i) { A.col(i).array() -= A.col(i).mean();}
     b.array() -= b.mean();
  }

  //scale covariates so that the 2norm is 1
  if(scaling) for(int i=0; i<A.cols(); ++i) { double s = A.col(i).norm(); S[i] = (s != 0) ? 1/s : 0; }
  else S.setOnes();

  dtime = 0;
  residual = b; 
  double bnorm  = b.norm();
  double abstol = reltol*bnorm;
  double C;

  //initialize starting set with maximaly corellated element
  {
    crlt = S.asDiagonal()*(A.transpose()*b);
    long int i;
    C = crlt.maxCoeff(&i);
    indices.push_back(i);
    
    //intitialize cholesky factorization
    B.col(0)     = S[indices[0]]*(A.col(indices[0]));
    double diagK = B.col(0).squaredNorm();
    double r     = sqrt(diagK);
    R(0,0)       = r;
  }
  
  //ensure that an element is not selected immediately after being dropped
  bool     dropId  = false;
  long int blockId = 0; 

  long int k       = 0; // current column index <-> Active set size -1
  long int iter    = 0; // number of iterations
  long int downIt  = 0; // number of downdates
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

    crlt = S.asDiagonal()*(A.transpose()*residual); // correlations

    long int i = 0; // dummy integer for stuff

    // solve the symetric system of equations
    wA.head(k+1)  = R.topLeftCorner(k+1,k+1).triangularView<Upper>().transpose().solve(rhs.head(k+1));
    wA.head(k+1)  = R.topLeftCorner(k+1,k+1).triangularView<Upper>().solve(wA.head(k+1));
    if(!dropId && wA[k] <= 0) {
       std::cout << "wA = " << wA.head(k+1).transpose() << std::endl;
       std::cout << "R.diag = " << R.diagonal().head(k+1).transpose() << std::endl;
       std::cout << "C = " << std::endl;
       for(std::vector<long int>::iterator it = indices.begin(); it != indices.end(); ++it) std::cout << crlt[*it] <<  " ";
       std::cout << "\n *** Roundoff Error *** " << std::endl;
       exit(-1);
       break;
    }
    double oneNwA = 1.0/sqrt(wA.head(k+1).sum());
    wA.head(k+1) *= oneNwA;
 
    // compute smallest angle at which a new covariant becomes dominant
    update = B.leftCols(k+1)*wA.head(k+1);
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
    t.head(k+1) = -1.0*ylar.head(k+1).array()/(wA.head(k+1).array());
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

        //update cholesky factorization R'*R = B'*B where R is upper triangular
        double diagK   = B.col(k+1).squaredNorm();
        colK.head(k+1) = B.leftCols(k+1).transpose()*(B.col(k+1));
        vk.head(k+1)   = R.topLeftCorner(k+1,k+1).triangularView<Upper>().transpose().solve(colK.head(k+1));
        double r       = sqrt(diagK - vk.head(k+1).squaredNorm());

        Block<MatrixXd,Dynamic,1,true> colR = R.col(k+1);
        colR.head(k+1) = vk.head(k+1);
        colR[k+1]      = r;

        //if diagonal element is too small, then column is near linearly dependent
        if(r != r /* check for nan*/ || r <= std::numeric_limits<double>::min() /* or too close to current columns*/ ) { 
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
    C -= gamma1*oneNwA;
    // update solution and estimate
    ylar.head(k+1) = ylar.head(k+1).array() + gamma1*(wA.head(k+1).array());
    residual -= gamma1*update;

    k++;
    iter++;
    if(dropId) {// if gamma_tilde is selected, remove that index and downdate Cholesky factorization
      dtime -= getTime();
      downIt++;

      long int leftSide = k;
      k = i;

      //remove selected index
      std::vector<long int>::iterator fol = indices.erase(indices.begin()+k);
       
      //zero out column
      R.col(k).head(k+1).setZero();
      ylar[k] = 0;

      //now zero out diagonal
      for(std::vector<long int>::iterator it = fol; it != indices.end(); ++it,k++){
        //construct Givens rotation
        double l = R.block<2,1>(k,k+1).norm();
        double c = R(k,k+1)/l;
        double s = -1.0*R(k+1,k+1)/l;
        Matrix2d G; G << c, -s, s, c;
        //apply to block of R to zero out old diagonal element
        R.block(k,k+1,2,leftSide-k-1) = G*R.block(k,k+1,2,leftSide-k-1);
        R.col(k) = R.col(k+1);
        B.col(k) = B.col(k+1);
        ylar[k] = ylar[k+1]; 
      }
      R.col(k).head(k+1).setZero();
      R.row(k).head(k+1).setZero();
      B.col(k).setZero();
      ylar[k] = 0;
      k--;
      dtime += getTime();
    }

  }

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
