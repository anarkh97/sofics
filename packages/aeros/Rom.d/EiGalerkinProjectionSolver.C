#ifndef ROM_EIGALERKINPROJECTIONSOLVER_C
#define ROM_EIGALERKINPROJECTIONSOLVER_C

#ifdef USE_EIGEN3
#include "EiGalerkinProjectionSolver.h"
#include "DistrVecBasisOps.h"
#include "PodProjectionSolver.h"
#include <stdexcept>
#include <vector>
#include <Solvers.d/eiquadprog.hpp>
#include <Solvers.d/ParallelSolver.h>
#include <Driver.d/Mpc.h>
#include <Paral.d/SubDOp.h>

namespace Rom {

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::GenEiSparseGalerkinProjectionSolver(
    const Connectivity *cn,
    const DofSetArray *dsa, const ConstrainedDSA *c_dsa, bool selfadjoint, double tol):
  GenEiSparseMatrix<Scalar>(cn, dsa, c_dsa, selfadjoint), 
  spMat(NULL),
  K(NULL),
  cdsa_(c_dsa),
  basisSize_(0), 
  dualBasisSize_(0),
  projectionBasis_(NULL), 
  dualProjectionBasis_(NULL),
  Empirical(false),
  selfadjoint_(selfadjoint),
  contact_(false),
  tol_(tol),
  startCol_(0),
  blockCols_(0),
  startDualCol_(0),
  dualBlockCols_(0),
  startq(0),
  endq(0),
  qsize(0),
  myID(0),
  grpSize(1),
  useX0(false)
{
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::GenEiSparseGalerkinProjectionSolver(
	const Connectivity *cn,
	const DofSetArray *dsa, const ConstrainedDSA *c_dsa, int numSub_,
	GenSparseMatrix<Scalar> **spMat_, bool selfadjoint, double tol,
	int grpSize_):
  GenEiSparseMatrix<Scalar>(cn, dsa, c_dsa, selfadjoint),
  spMat(spMat_),
  cdsa_(c_dsa),
  basisSize_(0),
  numSub(numSub_),
  dualBasisSize_(0),
  projectionBasis_(NULL), 
  dualProjectionBasis_(NULL),
  Empirical(false),
  selfadjoint_(selfadjoint),
  contact_(false),
  tol_(tol),
  startCol_(0),
  blockCols_(0),
  startDualCol_(0),
  dualBlockCols_(0),
  startq(0),
  endq(0),
  qsize(0),
  myID(0),
  grpSize(grpSize_),
  useX0(false)
{
  K = new GenSubDOp<Scalar>(numSub,spMat);

  // if running with mpi, split the communicator and allocate working arrays
#ifdef USE_MPI
  if(structCom && grpSize > 1){
    int world_rank = structCom->myID();
    int world_size = structCom->numCPUs();
    int color      = world_rank / grpSize;
    int remainder  = world_size % grpSize;
    // check their are no groups of size 1
    if(remainder == 1) {
      // if this color group is of size 1, subtract one to put in previous group
      if( color*grpSize == world_size ) 
        color -= 1; 
    }
    MPI_Comm_split(MPI_COMM_WORLD, color, world_rank, &JacobiComm);
    MPI_Comm_rank(JacobiComm, &myID);
    recvcounts.resize(grpSize);
    displs.resize(grpSize);
    filePrint(stderr," ... Reduced System Solver: %d procs ... \n", grpSize);
  }
#endif
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::setLocalBasis(int startCol, int blockCols) 
{
  startCol_  = startCol;
  blockCols_ = blockCols;
  reducedMatrix_.resize(blockCols_,blockCols_);
  reducedMatrix_.setZero();
#ifdef USE_MPI
  if(structCom && grpSize > 1){
    x0.resize(blockCols_); x0.setZero();

    qsize  = ceil(blockCols_/grpSize); //first compute ceiling of qsize
    startq = myID*qsize;

    if(myID == grpSize - 1)
      endq = blockCols_ - 1;
    else
      endq = (myID+1)*qsize - 1;

    qsize = endq - startq + 1; // then truncate to get actual qsize for this processor

    recvcounts[myID] = qsize;
    displs[myID]     = startq;

    MPI_Allreduce(MPI_IN_PLACE,  &recvcounts.data()[0], recvcounts.size(), MPI_INT, MPI_SUM, JacobiComm);
    MPI_Allreduce(MPI_IN_PLACE,  &displs.data()[0],     displs.size(),     MPI_INT, MPI_SUM, JacobiComm);
  } 
#endif
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::setLocalDualBasis(int startDualCol, int dualBlockCols)
{
  startDualCol_  = startDualCol;  // set which columns are to be used in the reduced constraint matrix
  dualBlockCols_ = dualBlockCols;
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::zeroAll()
{
  GenEiSparseMatrix<Scalar>::zeroAll();
  reducedMatrix_.setZero();
  if(K) K->zeroAll();
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::storeReducedMass(Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> &VtMV_)
{
  VtMV = VtMV_;
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::addReducedMass(double Mcoef)
{
  if(VtMV.rows() > 0)
    reducedMatrix_ += Mcoef*VtMV.block(startCol_,startCol_,blockCols_,blockCols_);
  else
    reducedMatrix_ += Mcoef*Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic>::Identity(blockCols_, blockCols_);
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::addToReducedMatrix(const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> &ContributionMat, 
                                                                                      double Coef)
{
  reducedMatrix_ += Coef*ContributionMat.block(startCol_,startCol_,blockCols_,blockCols_);
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::addMPCs(Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> &ContributionMat, Eigen::Matrix<Scalar,Eigen::Dynamic,1> &WtRhs, double Kcoef)
{
  if(!selfadjoint_) { std::cerr << "Error: unsymmetric solver is not supported for contact ROM with Lagrange Multipliers\n"; exit(-1); }

  reducedConstraintMatrix_ = Kcoef*ContributionMat; 
  reducedConstraintRhs0_ = reducedConstraintRhs_ = Kcoef*WtRhs; 
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::addLMPCs(int numLMPC, LMPCons **lmpc, double Kcoef)
{
  if(numLMPC > 0 && !selfadjoint_) { std::cerr << "Error: unsymmetric solver is not supported for contact ROM with Lagrange Multipliers\n"; exit(-1); }

  // loop through mpcs and multiply by dual basis
  // take (C^T*W) \in R^(N x m) and multiply by V^T 

  std::vector<Eigen::Triplet<Scalar> > tripletList;

  Eigen::SparseMatrix<Scalar> C(numLMPC, cdsa_->size());
  Eigen::Matrix<Scalar,Eigen::Dynamic,1> g(numLMPC);

  for(int i=0; i<numLMPC; ++i) {
    for(int j=0; j<lmpc[i]->nterms; ++j) {
      int cdof = cdsa_->locate(lmpc[i]->terms[j].nnum, 1 << lmpc[i]->terms[j].dofnum);
      if(cdof > -1) {
        tripletList.push_back(Eigen::Triplet<Scalar>(i, cdof, Scalar(lmpc[i]->terms[j].coef.r_value)));
      }
    }
    g[i] = lmpc[i]->rhs.r_value;
  }

  C.setFromTriplets(tripletList.begin(), tripletList.end());

  // this is only called once, so set W and V to all columns, 
  // if using local basis for either W or V, 
  // then select the correct block diagonal element
  const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> &V = projectionBasis_->basis();
  const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> &W = dualProjectionBasis_->basis();
  reducedConstraintMatrix_ = Kcoef*W.transpose()*C*V; 
  reducedConstraintRhs0_ = reducedConstraintRhs_ = Kcoef*W.transpose()*g; 
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::addModalLMPCs(double Kcoef, int Wcols, std::vector<double>::const_iterator it, std::vector<double>::const_iterator it_end)
{
  filePrint(stderr," ... Using Modal LMPCs              ...\n");
  dualBasisSize_  = dualBlockCols_ = Wcols;
  int counter = 0; int column = 0; int row = 0;
  reducedConstraintMatrix_.setZero(dualBasisSize_,basisSize_); //allocate enough space for all local bases
  reducedConstraintRhs0_.setZero(dualBasisSize_);         
  // set reduced Constraint Matrix
  while(counter < dualBasisSize_*basisSize_) {
    reducedConstraintMatrix_(row,column) = *it; row++; counter++; it++;
    if(row == dualBasisSize_) {
      row = 0;
      column++;
    }
  }
  reducedConstraintMatrix_ *= Kcoef;
  // set reduced Constraint RHS
  row = 0;
  while(it != it_end) {
    reducedConstraintRhs0_(row) = *it; row++; it++;
  }

  reducedConstraintRhs0_ *= Kcoef;
  reducedConstraintRhs_   = reducedConstraintRhs0_;

  reducedConstraintForce_.setZero(basisSize_);
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::updateLMPCs(GenVecType<Scalar> &_q)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > q(_q.data(), _q.size()); 
  reducedConstraintRhs_ = reducedConstraintRhs0_ - reducedConstraintMatrix_*q;
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::projectionBasisIs(GenVecBasis<Scalar,GenVecType> &reducedBasis)
{
/*  if (reducedBasis.globalVectorSize() != GenEiSparseMatrix<Scalar>::neqs()) {
    fprintf(stderr,"Basis Vector length: %d, Number of Unconstrained DOFs: %d \n",reducedBasis.globalVectorSize(), GenEiSparseMatrix<Scalar>::neqs());
    throw std::domain_error("Vectors of the reduced basis have the wrong size");
  }*/

  projectionBasis_ = &reducedBasis;
  basisSize_       = reducedBasis.vectorCount();

  reducedMatrix_.setZero(basisSize_, basisSize_);

  // local bases: solver uses all bases unless setLocalBasis is called
  startCol_  = 0;
  blockCols_ = basisSize_;

#ifdef USE_MPI
  if(structCom && grpSize > 1){
    x0.resize(blockCols_); x0.setZero();

    qsize  = ceil(blockCols_/grpSize); //first compute ceiling of qsize
    startq = myID*qsize;

    if(myID == grpSize - 1)
      endq = blockCols_ - 1;
    else
      endq = (myID+1)*qsize - 1;

    qsize = endq - startq + 1; // then truncate to get actual qsize for this processor

    recvcounts[myID] = qsize;
    displs[myID]     = startq;

    MPI_Allreduce(MPI_IN_PLACE,  &recvcounts.data()[0], recvcounts.size(), MPI_INT, MPI_SUM, JacobiComm);
    MPI_Allreduce(MPI_IN_PLACE,  &displs.data()[0],     displs.size(),     MPI_INT, MPI_SUM, JacobiComm);
  }
#endif
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::dualProjectionBasisIs(GenVecBasis<Scalar,GenVecType> &dualReducedBasis)
{
  dualProjectionBasis_ = &dualReducedBasis;
  dualBasisSize_       = dualReducedBasis.vectorCount();
  reducedConstraintMatrix_.setZero(dualBasisSize_, basisSize_);
  reducedConstraintForce_.setZero(basisSize_);

  // local basis: dual solver uses all columns unless setLocalDualBasis is called
  startDualCol_  = 0; 
  dualBlockCols_ = dualBasisSize_;
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::dualProjectionBasisIs(std::vector<std::map<int, double> > &dualReducedBasis)
{
  //dualProjectionBasis_ = &dualReducedBasis;
  dualBasisSize_       = dualReducedBasis.size();
  reducedConstraintMatrix_.setZero(dualBasisSize_, basisSize_);
  reducedConstraintForce_.setZero(basisSize_);

  // local basis: dual solver uses all columns unless setLocalDualBasis is called
  startDualCol_  = 0;
  dualBlockCols_ = dualBasisSize_;
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::EmpiricalSolver()
{
  std::cout << " ...    Empirical Solver Selected   ..." << std::endl;
  Empirical = true;
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::factor()
{
  LocalBasisType V = projectionBasis_->basis().block(0,startCol_,projectionBasis_->size(),blockCols_);

  if(selfadjoint_ && !Empirical) {
    reducedMatrix_.template triangularView<Eigen::Lower>()
    += V.transpose()*(this->M.template selfadjointView<Eigen::Upper>()*V);
    if(dualBlockCols_ > 0) c1_ = reducedMatrix_.trace();
    llt_.compute(reducedMatrix_);
  }
  else {
    if(Empirical) {
      lu_.compute(reducedMatrix_);
    } else {
      reducedMatrix_ += V.transpose()*(this->M*V);
      lu_.compute(reducedMatrix_);
    }
  }
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::refactor()
{
  // do nothing
}

template <>
void
GenEiSparseGalerkinProjectionSolver<double,GenDistrVector,GenParallelSolver<double> >::refactor() // only called in parallel code branch
{ // for parallel implicit ROM, V^T*M*V is computed in problem descriptor and passed through
  // addReducedMass and addToReducedMatrix
  if(selfadjoint_ && !Empirical) {
   
    GenFullSquareMatrix<double> K_reduced; // local data structure
    calculateReducedStiffness(*K, *projectionBasis_, K_reduced, selfadjoint_); // parallel multiplication

    // upper half of K_reduced is filled, but data structure is row major, Krmap is Column major
    Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > Krmap(K_reduced.data(), blockCols_, blockCols_);
    reducedMatrix_.triangularView<Eigen::Lower>() 
    += Krmap; 

    //std::cout << "Kr = \n" << Krmap << std::endl;

    if(dualBlockCols_ > 0) c1_ = reducedMatrix_.trace();
    if (structCom && grpSize > 1) // if using block Jacobi, only factorize diagonal sub-matrix
      llt_.compute(reducedMatrix_.block(startq,startq,qsize,qsize));
    else
      llt_.compute(reducedMatrix_);
  }
  else {
    if(Empirical) {
      lu_.compute(reducedMatrix_); 
    } else { 

      GenFullSquareMatrix<double> K_reduced;
      calculateReducedStiffness(*K, *projectionBasis_, K_reduced);

      Eigen::Map< Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > Krmap(K_reduced.data(), blockCols_, blockCols_);
      reducedMatrix_ += Krmap;

      //std::cout << "Kr = \n" << Krmap << std::endl;

      if(structCom && grpSize > 1)
        lu_.compute(reducedMatrix_.block(startq,startq,qsize,qsize));
      else
        lu_.compute(reducedMatrix_);
    }
  }
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
Eigen::Matrix<Scalar, Eigen::Dynamic, 1>
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::blockJacobi(Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > &b) const {
#ifdef USE_MPI
  // x is the right hand side. Initial guess is the zero vector 
  //compute stopping citeria
  double target = 1e-6*b.segment(startq,qsize).squaredNorm();
  MPI_Allreduce(MPI_IN_PLACE,  &target, 1, MPI_DOUBLE, MPI_SUM, JacobiComm);

  // allocate and initialize working arrays
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> dx(blockCols_);
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1>  x(blockCols_);
  // initialize the local residual and solution increment
  if(useX0) {
    if(selfadjoint_) {
      dx.segment(startq,qsize) = -reducedMatrix_.block(  startq,     0,                qsize,startq)*x0.head(startq)
                                - reducedMatrix_.block(endq + 1,startq,blockCols_ - endq - 1, qsize).transpose()*x0.tail(blockCols_ - endq - 1) 
                                + b.segment(startq,qsize);
    } else {
      dx.segment(startq,qsize) = -reducedMatrix_.block(startq,       0,qsize,               startq)*x0.head(startq)
                                - reducedMatrix_.block(startq,endq + 1,qsize,blockCols_ - endq - 1)*x0.tail(blockCols_ - endq - 1) 
                                + b.segment(startq,qsize);
    }
  } else {
    dx.segment(startq,qsize) = b.segment(startq,qsize);
  }

  if(selfadjoint_)
    x.segment(startq,qsize) = dx.segment(startq,qsize) = (llt_.solve(dx.segment(startq,qsize))).eval();
  else
    x.segment(startq,qsize) = dx.segment(startq,qsize) = (lu_.solve(dx.segment(startq,qsize))).eval();

  if(useX0)
    dx.segment(startq,qsize) -= x0.segment(startq,qsize);

  MPI_Allgatherv(MPI_IN_PLACE, recvcounts[myID], MPI_DOUBLE, dx.data(), &recvcounts.data()[0], &displs.data()[0], MPI_DOUBLE, JacobiComm);

  int maxIt = 1000; // make it big cause yolo
  for(int i = 0; i < maxIt; ++i){

     // get increment norm
     double dxnorm = dx.segment(startq,qsize).squaredNorm();
     MPI_Allreduce(MPI_IN_PLACE,  &dxnorm, 1, MPI_DOUBLE, MPI_SUM, JacobiComm); 
     if(dxnorm <= target){
       // filePrint(stderr,"Converged after %d iterations (target: %3.2e, ||dx||: %3.2e)\n",i, target, dxnorm);
       break;
     }

     // update local segment of state increment
     if(selfadjoint_){// K is lower triangular, use symmetry for multiplication
       dx.segment(startq,qsize) = -reducedMatrix_.block(  startq,     0,                qsize,startq)*dx.head(startq)
                                 - reducedMatrix_.block(endq + 1,startq,blockCols_ - endq - 1, qsize).transpose()*dx.tail(blockCols_ - endq - 1);
     } else {
       dx.segment(startq,qsize) = -reducedMatrix_.block(startq,       0,qsize,               startq)*dx.head(startq)
                                 - reducedMatrix_.block(startq,endq + 1,qsize,blockCols_ - endq - 1)*dx.tail(blockCols_ - endq - 1);
     }

     // solve small system
     if(selfadjoint_)
       dx.segment(startq,qsize) = (llt_.solve(dx.segment(startq,qsize))).eval();
     else
       dx.segment(startq,qsize) =  (lu_.solve(dx.segment(startq,qsize))).eval();

     // update state estimate
     x.segment(startq,qsize) += dx.segment(startq,qsize);

     // communicate increment to other processors. 
     MPI_Allgatherv(MPI_IN_PLACE, recvcounts[myID], MPI_DOUBLE, dx.data(), &recvcounts.data()[0], &displs.data()[0], MPI_DOUBLE, JacobiComm);

     if(i == maxIt - 1)
       filePrint(stderr,"Warning: Parallel ROM Solver did not Converge after %d iterations (target: %3.2e, ||dx||: %3.2e )\n",i, target, dxnorm);

  }
  // after convergence, communicate state estimate
  MPI_Allgatherv(MPI_IN_PLACE, recvcounts[myID], MPI_DOUBLE, x.data(), &recvcounts.data()[0], &displs.data()[0], MPI_DOUBLE, JacobiComm);

  return x;
#endif
  throw "Unimplemented blockJacobi was called.";
}



template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::reSolve(GenVecType<Scalar> &rhs)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > x(rhs.data()+startCol_, blockCols_);

  //std::cout << "rhs = \n" << x.transpose() << std::endl; 

  // llt_  : precomputed cholesky decomposition
  // c1_   : estimate for condition number of dynamic stiffness matrix
  // g0    : residual  
  // CE    : Equality constraints
  // ce0   : rhs for equality constraints
  // CI    : Inequality constraints
  // ci0   : rhs for inequality constraints
  // _x    : primal solution
  // Lambda: multipliers associated with equality constraints
  // Mu    : multipliers associated with inequality constraints

  //double dummyTime = -1.0*getTime();
  if(dualBlockCols_ > 0 && contact_) {
    Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> CE(0,0), CI = -reducedConstraintMatrix_.block(startDualCol_,startCol_,dualBlockCols_,blockCols_).transpose();
    Eigen::Matrix<Scalar,Eigen::Dynamic,1> g0 = -x, ce0(0,1), _x(blockCols_), Lambda(0,1), Mu(dualBlockCols_);
    solve_quadprog2(llt_, c1_, g0, CE, ce0, CI, reducedConstraintRhs_.segment(startDualCol_,dualBlockCols_), _x, &Lambda, &Mu, tol_);
    x = _x;
    reducedConstraintForce_.setZero();
    reducedConstraintForce_.segment(startCol_,blockCols_) = reducedConstraintMatrix_.block(startDualCol_,startCol_,dualBlockCols_,blockCols_).transpose()*Mu;
  }
  else if(selfadjoint_ && !Empirical) {
    if(structCom && grpSize > 1)
      x0 = x = blockJacobi(x);
    else
      llt_.solveInPlace(x); 
  } else  {
    if(structCom && grpSize > 1)
      x0 = x = blockJacobi(x);
    else
      x = (lu_.solve(x)).eval();
  }

  //std::cout << "x = \n" << x.transpose() << std::endl;

  // zero out unused parts of rhs
  for(int i=0; i<startCol_; ++i) rhs[i] = 0;
  for(int i=startCol_+blockCols_; i<rhs.size(); ++i) rhs[i] = 0;
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
void
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::solve(const GenVecType<Scalar> &rhs, GenVecType<Scalar> &sol)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > b(rhs.data()+startCol_, blockCols_), x(sol.data()+startCol_, blockCols_);
  sol.zero();

  //double dummyTime = -1.0*getTime();
  if(dualBlockCols_ > 0 && contact_) {
    Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> CE(0,0), CI = -reducedConstraintMatrix_.block(startDualCol_,startCol_,dualBlockCols_,blockCols_).transpose();
    Eigen::Matrix<Scalar,Eigen::Dynamic,1> g0 = -b, ce0(0,1), _x(blockCols_), Lambda(0,1), Mu(dualBlockCols_);
    solve_quadprog2(llt_, c1_, g0, CE, ce0, CI, reducedConstraintRhs_.segment(startDualCol_,dualBlockCols_), _x, &Lambda, &Mu, tol_);
    x = _x;
    reducedConstraintForce_.setZero();
    reducedConstraintForce_.segment(startCol_,blockCols_) = reducedConstraintMatrix_.block(startDualCol_,startCol_,dualBlockCols_,blockCols_).transpose()*Mu;
  }
  else if(selfadjoint_ && !Empirical){
     if(structCom && grpSize > 1){
       useX0 = true;
       x0 = x = blockJacobi(b);
       useX0 = true;
     } else {
       x = llt_.solve(b);
     }
  } else {
     if(structCom && grpSize > 1) {
       useX0 = true;
       x0 = x = blockJacobi(b);
       useX0 = true;
     } else {
       x = lu_.solve(b);
     }
  }

  //dummyTime += getTime();
  //reducedSolveTime += dummyTime/1000.0;
  //nSolve++;
  //filePrint(stderr,"Average Solve Time: %3.2e \n",reducedSolveTime/double(nSolve));
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
double
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::getResidualNorm(const GenVecType<Scalar> &v)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > vj(v.data()+startCol_, blockCols_);
  return vj.norm();
}

template <typename Scalar, template<typename> class GenVecType, class BaseSolver>
double
GenEiSparseGalerkinProjectionSolver<Scalar,GenVecType,BaseSolver>::getFNormSq(GenVecType<Scalar> &v)
{
  Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > vj(v.data()+startCol_, blockCols_);
  double dummy = vj.norm();
  return dummy*dummy;
}

} /* end namespace Rom */

template class Rom::GenEiSparseGalerkinProjectionSolver<double, GenDistrVector, GenParallelSolver<double> >;
template class Rom::GenEiSparseGalerkinProjectionSolver<std::complex<double>, GenDistrVector, GenParallelSolver<std::complex<double> > >;

#endif /* USE_EIGEN3 */

#endif /* ROM_EIGALERKINPROJECTIONSOLVER_C */
