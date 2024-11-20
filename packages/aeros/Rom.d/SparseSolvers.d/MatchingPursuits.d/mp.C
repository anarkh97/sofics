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
#include <set>

Eigen::VectorXd
nncgp(Eigen::Ref<Eigen::MatrixXd> A, Eigen::Ref<Eigen::VectorXd> b, double& rnorm,
      long int &info, double maxsze, int &maxEle, double maxite, double reltol, bool verbose, bool scaling, bool center, bool reverse, double &dtime, std::vector<long int> &indices);

Eigen::VectorXd
mp(Eigen::Ref<Eigen::MatrixXd> A, Eigen::Ref<Eigen::VectorXd> b, double& rnorm,
      long int &info, double maxsze, double maxite, double reltol, bool verbose, bool scaling, bool center, bool positive)
{
  using namespace Eigen;
  
  const long int m = b.rows();
  const long int maxvec = std::min(m, (long int)(maxsze*A.cols()));
  const long int maxit = maxite*A.cols();

  // allocate
  VectorXd ymp(maxvec);
  VectorXd crlt(A.cols()), S(A.cols());
  VectorXd residual(A.rows());
 
  MatrixXd B(A.rows(),maxvec);

  // intitialize
  ymp.setZero();  
  B.setZero();

  info = 1;

  std::vector<long int> indices;
  std::set<long int> nld_indices;

  for(int i=0; i<A.cols(); ++i) { double s = A.col(i).norm(); S[i] = (s != 0) ? 1/s : 0; }

  double  bnorm = b.norm();
  double abstol = reltol*bnorm;
  rnorm = bnorm;
  residual = b;
  double C = 0;

  long int i       = 0; // max coefficient index
  long int k       = 0; // current column index <-> Active set size -1
  long int iter    = 0; // number of iterations
  while(true) {

    if(verbose) {
      std::cout << "Iteration = " << std::setw(9) << iter << "    "
                << "Active set size = " << std::setw(9) << nld_indices.size() << "    "
                << "Max Correlation = " << std::setw(13) << std::scientific << std::uppercase << C << "    "
                << "Residual norm = " << std::setw(13) << std::scientific << std::uppercase << rnorm << "    "
                << "Target = " << std::setw(13) << std::scientific << std::uppercase << abstol << std::endl;
      std::cout.unsetf(std::ios::scientific);
      std::cout.unsetf(std::ios::uppercase);
    }

    if(rnorm <= abstol || k == maxvec) break;
    if(iter >= maxit) { info = 3; break; }

    crlt = S.asDiagonal()*(A.transpose()*residual); // correlations

    if(positive){
      for(std::set<long int>::iterator it = nld_indices.begin(); it!= nld_indices.end(); ++it) crlt[*it] = -std::numeric_limits<double>::max(); 
      C = crlt.maxCoeff(&i);
      B.col(k) = A.col(i);
      if(C <= 0) break;
    } else {
      for(std::set<long int>::iterator it = nld_indices.begin(); it!= nld_indices.end(); ++it) crlt[*it] = -std::numeric_limits<double>::min();
      crlt.cwiseAbs().maxCoeff(&i);
      C = crlt[i];
    }

    indices.push_back(i);  
    nld_indices.insert(i);

    ymp[k] = C; 

    residual -= S[i]*A.col(i)*C; 

    rnorm = residual.norm();

    k++;
    iter++;
  }

  double dtime = 0;
  int dummy = 0;
  std::vector<long int> newIndices; 
  ymp.head(k) = nncgp(B.leftCols(k), b, rnorm, info, maxsze, dummy, maxite, reltol, verbose, scaling, center, false, dtime, newIndices);  

  if(verbose) std::cout.flush();

  VectorXd x = VectorXd::Zero(A.cols());
  for(long int j=0; j<k; ++j) x[indices[j]] = S[indices[j]]*ymp[j];
  return x;
}

#endif
