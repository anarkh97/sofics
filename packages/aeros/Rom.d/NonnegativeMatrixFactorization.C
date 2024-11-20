#include "NonnegativeMatrixFactorization.h"

#include <Utils.d/linkfc.h>

#include <cmath>
#include <iostream>
#include <vector>
#include <cstdio>

#ifdef USE_EIGEN3
#include <Eigen/Core>
#include <Eigen/Dense>
#endif

extern "C" {
  void _FORTRAN(spnnls)(double *a, const long int *mda, const long int *m, const long int *n,
                        const double *b, double *x, const double *reltol, double *rnorm, double *w,
                        double *zz, double *zz2, long int *index, long int *mode, long int *prtflg,
                        long int *sclflg, const double *maxsze, const double *maxite, double *dtime);
}

namespace Rom {

NonnegativeMatrixFactorization::NonnegativeMatrixFactorization(int maxBasisDimension, int method) :
  rowCount_(0),
  colCount_(0),
  basisDimension_(maxBasisDimension),
  maxBasisDimension_(maxBasisDimension),
  numRandInit_(1),
  method_(method),
  matrixBuffer_(0),
  robBuffer_(0),
  maxIter_(100),
  nmfcAlpha(0.0),
  nmfcBeta(0.0),
  nmfcGamma(0.0)
{
  bFile = fopen("LambdaBasisFile.m","w+");
}

void
NonnegativeMatrixFactorization::matrixSizeIs(int rows, int cols) {
  matrixBuffer_.sizeIs(rows * cols);

  rowCount_ = rows;
  colCount_ = cols;
}

void
NonnegativeMatrixFactorization::robSizeIs(int rows, int cols) {
  robBuffer_.sizeIs(rows * cols);
  maxBasisDimensionIs(cols);
}

void
NonnegativeMatrixFactorization::maxIterIs(int maxIter) {
  maxIter_ = maxIter;
}

void
NonnegativeMatrixFactorization::toleranceIs(double tol) {
  tol_ = tol;
}

void
NonnegativeMatrixFactorization::nmfcAlphaIs(double alp) {
  nmfcAlpha = alp;
}

void
NonnegativeMatrixFactorization::nmfcBetaIs(double bet) {
  nmfcBeta = bet;
}

void
NonnegativeMatrixFactorization::nmfcGammaIs(double gam) {
  nmfcGamma = gam;
}

void
NonnegativeMatrixFactorization::maxBasisDimensionIs(int maxBasisDimension) {
  maxBasisDimension_ = maxBasisDimension;
}

void
NonnegativeMatrixFactorization::basisDimensionIs(int basisDimension) {
  basisDimension_ = basisDimension;
}

void 
NonnegativeMatrixFactorization::numRandInitIs(int numRandInit) {
  numRandInit_ = numRandInit;
}

void
NonnegativeMatrixFactorization::solve(int basisDimensionRestart) {
#ifdef USE_EIGEN3
  Eigen::Map<Eigen::MatrixXd> X(matrixBuffer_.array(), rowCount_, colCount_);

  std::vector<int> cols;
  for(int i=0; i<colCount_; ++i) if(!(X.col(i).array() == 0).all()) cols.push_back(i);
  std::cerr << "X has " << X.cols() << " columns of which " << cols.size() << " are non-zero\n";

  std::vector<int> rows;
  for(int i=0; i<rowCount_; ++i) if(!(X.row(i).array() == 0).all()) rows.push_back(i);
  std::cerr << "X has " << X.rows() << " rows of which " << rows.size() << " are non-zero\n";

  Eigen::Map<Eigen::MatrixXd> ROB(robBuffer_.array(), rowCount_, maxBasisDimension_);
  int m = rows.size();
  int n = cols.size();
  const int &k = basisDimension_;
  std::cerr << "Factoring data into " << k << " nonnegative columns " << std::endl;
  Eigen::MatrixXd A(m,n);
  FILE * pFile;
  pFile = fopen("LambdaSnapFile.m","w+");
  fprintf(pFile,"A = [\n");
  for(int i=0; i<m; ++i){
    for(int j=0; j<n; ++j){
      double elem = X(rows[i],cols[j]);
      double ebuf = elem*elem;
      A(i,j) = sqrt (ebuf); // note -ve is due to sign convention (Lagrange multipliers are negative in Aero-S)
      fprintf(pFile," %7.6e",A(i,j));
    }
    fprintf(pFile,"\n");
  }
  fprintf(pFile,"];");
  fclose(pFile); 
  Eigen::MatrixXd W(m,k), H(k,n);

  double res, index;
  Eigen::MatrixXd Err(A);
  switch(method_) { 
    default: case 1 : { // NMF ROB
      Eigen::MatrixXd W_min(W);
      double res_min;
      double nrmA = A.norm();
      for (int iRandom=0; iRandom<numRandInit_; ++iRandom) {
        // Initialization
        if (basisDimensionRestart>0) { 
          std::cout << "Initial guess based on smaller factorization" << std::endl;
          for (int i=0; i<m; ++i)
            for (int j=0; j<basisDimensionRestart; ++j)
              W(i,j) = ROB(rows[i],j);
          W.rightCols(k-basisDimensionRestart) = 0.5*(Eigen::MatrixXd::Random(m,k-basisDimensionRestart)+Eigen::MatrixXd::Ones(m,k-basisDimensionRestart));
        } else {
          W = 0.5*(Eigen::MatrixXd::Random(m,k) + Eigen::MatrixXd::Ones(m,k));
        }

        Eigen::MatrixXd W_copy(W);
        for(int i=0; i<maxIter_; ++i) {
          H = solveNNLS_MRHS(W, A);
          W.transpose() = solveNNLS_MRHS(H.transpose(),A.transpose());
          Eigen::MatrixXd Z = W*H; 
          Err = A-Z;
          index = findColumnWithLargestMagnitude(Err);
          res = Err.norm()/nrmA;

          double maskErr = 0.0;
          // enforce equality of non-zero components
          for(int col = 0; col < n; ++col){
            for(int row = 0; row < m; ++row) {
              if(A(row,col) > 1e-16){
                maskErr += pow(A(row,col) - Z(row,col),2.0);
              }
            }
          }

          double inc = (W-W_copy).norm();
          fprintf(stderr,"\r Iteration: %d, increment: %1.4e, Norm(P(W*H - X))/norm(X): %1.4e, rel. residual %1.4e, Stopping tolerance: %1.4e", i+1, inc, sqrt(maskErr)/nrmA, res, tol_);
          if(inc < tol_) break;
          W_copy = W;
        }
        fprintf(stderr,"\n");
        if (numRandInit_>1) {
          if (iRandom==0) {
            res_min = res;
            W_min = W;
          }
          else {
            if (res_min>res) {
              res_min = res;
              W_min = W;
            }
            if (iRandom==numRandInit_-1)
              W = W_min;
              std::cout << "NNMF basis retained with rel. residual = " << res_min << std::endl;
          }
        }
      }
    
    } break;
    case 2 : { // Greedy ROB
      W.setZero();
      H.setZero();
      // first vector is vector with largest magnitude
      index = findColumnWithLargestMagnitude(A);
      for (int i=1; i<=k; ++i) {
        W.col(i-1) = A.col(index);
        H.topRows(i) = solveNNLS_MRHS(W.leftCols(i), A);
        Err = A-W.leftCols(i)*H.topRows(i); 
        res = Err.norm()/A.norm();
        index = findColumnWithLargestMagnitude(Err);
        std::cout << "greedy iteration = " << i << ", rel. residual = " << res << ", maximum error = " << Err.col(index).norm() << std::endl;
      }
      // normalize vectors in W
      for (int i=0; i<k; ++i) W.col(i).normalize();
    } break;
    case 4 : { // Alternating direction algorithm for Nonnegative matrix completion,  Xu, Yin, Wen, Zhang
      // initialize working arrays
      std::cout << "... Nonnegative Matrix Completion ..." << std::endl;
      W.setZero(); H = Eigen::MatrixXd::Random(k,n);
      Eigen::MatrixXd U(m,k), V(k,n); U.setZero(); V.setZero();
      Eigen::MatrixXd Lambda(m,k), Pi(k,n); Lambda.setZero(); Pi.setZero();
      Eigen::MatrixXd Z = A; 
     
      // set solver coefficients 
      double alpha = (nmfcAlpha > 1e-16) ? nmfcAlpha : Z.lpNorm<Eigen::Infinity>();
      double beta  = (nmfcBeta  > 1e-16) ? nmfcBeta  : Z.lpNorm<Eigen::Infinity>();
      double gamma = (nmfcGamma > 1e-16) ? nmfcGamma : 1.618;
   
      fprintf(stderr,"... α =% 3.2e, β =%3.2e, γ= %3.2e ...\n", alpha, beta, gamma);

      Eigen::MatrixXd rhs1; rhs1.resize(m,k); rhs1.setZero();
      Eigen::MatrixXd rhs2; rhs2.resize(k,n); rhs2.setZero();
      Eigen::MatrixXd lhs;  lhs.resize(k,k);  lhs.setZero();

      Eigen::MatrixXd W_copy(W);

      double nrmX = A.norm();
      //begin nonnegative matrix completion
      for (int i = 0; i <= maxIter_; ++i) {

        W_copy = W;  
        // fist update W
        rhs1 = Z*H.transpose() + alpha*U - Lambda;
        lhs  = H*H.transpose() + alpha*Eigen::MatrixXd::Identity(k,k);
  
        W.transpose() = lhs.selfadjointView<Eigen::Lower>().llt().solve(rhs1.transpose());

        // then update H
        rhs2 = W.transpose()*Z + beta*V - Pi;
        lhs  = W.transpose()*W + beta*Eigen::MatrixXd::Identity(k,k);

        H = lhs.selfadjointView<Eigen::Lower>().llt().solve(rhs2);

        // compute current approximation
        Z = W*H;

        // update constraints
        U = W + (1/alpha)*Lambda;
        V = H + (1/beta)*Pi;

        double stpcrt = 0.0;
        // enforce equality of non-zero components
        for(int col = 0; col < n; ++col){
          for(int row = 0; row < m; ++row) {
            if(A(row,col) > 1e-16){
              stpcrt += pow(A(row,col) - Z(row,col),2.0);
              Z(row,col) = A(row,col);
            }
          }
        }

        // threshold constraints 
        for(int col = 0; col < k; ++col){
          for(int row = 0; row < m; ++row) {
            if(U(row,col) < 0.)
              U(row,col) = 0.;
          }
        } 

        for(int col = 0; col < n; ++col){
          for(int row = 0; row < k; ++row) {
            if(V(row,col) < 0.)
              V(row,col) = 0.;
          }
        }

        // update Lagrange multipliers
        Lambda += gamma*alpha*(W - U);
        Pi     += gamma*beta*(H - V);

        double inc = (W-W_copy).norm();
        fprintf(stderr,"\r Iteration: %d, increment: %1.4e, Norm(P(W*H - X))/norm(X) %1.4e, Stopping tolerance: %1.4e", i, inc, sqrt(stpcrt)/nrmX, tol_);

        if((sqrt(stpcrt) < tol_*nrmX) || (i > maxIter_) || (inc < tol_)){
          W = U; 
          if(sqrt(stpcrt) <tol_*nrmX)
            fprintf(stderr,"\nStopping criteria reached\n");
          if(i > maxIter_)
            fprintf(stderr,"\nMaximum iteration exceeded\n");
          break;
        }

      }
      fprintf(stderr,"\n");
      W = U;  
    } break;
  }

  fprintf(bFile,"W%d = [\n",k);
  // copy W into buffer
  ROB.setZero();
  for(int i=0; i<m; ++i){
    for(int j=0; j<k; ++j){
      ROB(rows[i],j) = W(i,j);
      fprintf(bFile," %7.6e", W(i,j));
    }
    fprintf(bFile,"\n");
  }
  fprintf(bFile,"];\n");
#endif
}

#ifdef USE_EIGEN3
Eigen::MatrixXd
NonnegativeMatrixFactorization::solveNNLS_MRHS(const Eigen::Ref<const Eigen::MatrixXd> &A, const Eigen::Ref<const Eigen::MatrixXd> &B)
{
  Eigen::MatrixXd X(A.cols(),B.cols());  
  for(int j=0; j<B.cols(); ++j) {
    X.col(j) = solveNNLS(A, B.col(j));
  }
  return X;
}

int
NonnegativeMatrixFactorization::findColumnWithLargestMagnitude(const Eigen::Ref<const Eigen::MatrixXd> &X)
{
  int index = 0;
  X.colwise().norm().maxCoeff(&index);
  return index;
}

Eigen::VectorXd
NonnegativeMatrixFactorization::solveNNLS(const Eigen::Ref<const Eigen::MatrixXd> &_A, const Eigen::Ref<const Eigen::VectorXd> &_b)
{
  Eigen::MatrixXd A(_A);
  const long int m = A.rows();
  const long int n = A.cols();
  Eigen::VectorXd b(_b);
  Eigen::VectorXd x(n);
  const double reltol = 1e-16;
  double rnorm;
  Eigen::Array<double,Eigen::Dynamic,1> w(n);
  Eigen::Array<double,Eigen::Dynamic,1> zz(m);
  Eigen::Array<double,Eigen::Dynamic,1> zz2(n);
  Eigen::Array<long int,Eigen::Dynamic,1> index(n);
  long int mode;
  long int prtflg = 0;
  long int scaflg = 1;
  const double maxsze = 1.0;
  const double maxite = 3.0;
  double dtime;

  _FORTRAN(spnnls)(A.data(), &m, &m, &n, b.data(), x.data(), &reltol, &rnorm, w.data(),
                   zz.data(), zz2.data(), index.data(), &mode, &prtflg, &scaflg, &maxsze, &maxite, &dtime);

  if(mode != 1) std::cerr << "Error: spnnls unsuccessful, mode = " << mode << std::endl;

  return x;
}
#endif

} /* end namespace Rom */
