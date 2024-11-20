#include <Solvers.d/GmresOrthoSet.h>
#include <iostream>
#include <iomanip>

template<class Scalar, class AnyVector, class AnyOperator, class LeftPreconditioner, class RightPreconditioner>
GmresSolver<Scalar, AnyVector, AnyOperator, LeftPreconditioner, RightPreconditioner>::GmresSolver(
  int _maxit, double _tol, AnyOperator *_op, void (AnyOperator::*_matvec)(AnyVector &, AnyVector &),
  LeftPreconditioner *_leftprec, void (LeftPreconditioner::*_applyLeft)(const AnyVector &, AnyVector &),
  RightPreconditioner *_rightprec, void (RightPreconditioner::*_applyRight)(const AnyVector &, AnyVector &),
  FSCommunicator *_com)
 : maxit(_maxit), tol(_tol), op(_op), matvec(_matvec), leftprec(_leftprec), applyLeft(_applyLeft), 
   rightprec(_rightprec), applyRight(_applyRight), oSetGMRES(NULL),
   maxortho(1000), printNumber(1), verbose(1)
{
  com = _com;
  rank = (_com) ? _com->cpuNum() : 0;
}

template<class Scalar, class AnyVector, class AnyOperator, class LeftPreconditioner, class RightPreconditioner>
GmresSolver<Scalar, AnyVector, AnyOperator, LeftPreconditioner, RightPreconditioner>::~GmresSolver()
{
  if(oSetGMRES) delete oSetGMRES;
}

template<class Scalar, class AnyVector, class AnyOperator, class LeftPreconditioner, class RightPreconditioner>
void
GmresSolver<Scalar, AnyVector, AnyOperator, LeftPreconditioner, RightPreconditioner>::reSolve(AnyVector &x)
{
  AnyVector b(x), x0(x), v(x), w(x), z(x);
  if(!oSetGMRES) oSetGMRES = new GmresOrthoSet<Scalar>(x.size(), maxortho, com);
  else oSetGMRES->reset();

  int &iter = m_info[0] = 0;

  double alpha = 0.0, beta; // parameters used to define normwise backward error
/*
  v = Matrix<Scalar, Dynamic, 1>::Ones(v.rows(),1); // v.setRandom();
  for(int i=0; i<20; ++i) { // power iterations to compute norm of M^{-1}*A
    (op->*matvec)(v, w); // for non-hermition A this should be A^H*v
    (leftprec->*applyLeft)(w, w); // for non-hermition M^{-1} this should be M{-1}^H*w
    (leftprec->*applyLeft)(w, w);
    (op->*matvec)(w, z);
    double n = (ip->*norm)(z);
    v = z; v /= n;
    alpha = sqrt(abs((ip->*dot)(v,z))); // for no left prec and symmetric op, alpha = abs((ip->*dot)(v,z));
  }
*/
  if(leftprec) (leftprec->*applyLeft)(b, v); else v = b;
  beta = v.norm();
//  std::cerr  << "alpha = " << alpha << ", beta = " << beta << std::endl;
 
  // x^0 = 0 (arbitrary)
  x.zero(); //x.setRandom();

  if(verbose && rank == 0 && printNumber > 0)
    std::cerr << " Iteration   Backward Error\n";

  bool stop = false;
  while(!stop) {
    x0 = x;

    // compute the residual v = M^{-1}(b-Ax^0)
    if(rightprec) (rightprec->*applyRight)(x0, z); else z = x0;
    (op->*matvec)(z,w); // w = A*x
    v = b - w; 
    if(leftprec) (leftprec->*applyLeft)(v, v);  // v = M^{-1}*r

    double resGMRES = oSetGMRES->init(v.data()); // v /= ||v||

    // Arnoldi iteration (Algorithm see Saad SISC) 
    for (int j = 0; j < maxortho; j++, iter++) {

      // check stopping criteria
      double error = resGMRES/(alpha*x.norm() + beta); // gmres normwise backward error
      //if(error <= iterInfo->cntl[iter::tol]) { (op->*matvec)(x, Ax); error = (ip->*norm)(M^{-1}*(b-Ax))/(alpha*(ip->*norm)(x) + beta); } // normwise backward error
      stop = (iter == maxit || error <= tol);

      if(verbose && rank == 0 && (stop || ((printNumber > 0) && (iter % printNumber == 0))))
        std::cerr << " " << std::setw(5) << iter << "     " << std::scientific << std::setprecision(8) << "  " << std::setw(14) << error << "  " << std::endl;

      if(stop) break;

      if(rightprec) (rightprec->*applyRight)(v, z); else z = v;
      (op->*matvec)(z, w); // w = A*v
      if(leftprec) (leftprec->*applyLeft)(w, w); // w = M^-1*w

      // Do Arnoldi step 
      resGMRES = oSetGMRES->ModorthoAdd(w.data(), v.data());

      // update solution (now at every iteration for backward error)
      oSetGMRES->solution(w.data()); // w = -V*y
      x = x0 - w;
    }

    if(!stop && verbose && rank == 0) std::cerr << "restarting GMRES\n";
  }
  if(rightprec) (rightprec->*applyRight)(x, x);
}

