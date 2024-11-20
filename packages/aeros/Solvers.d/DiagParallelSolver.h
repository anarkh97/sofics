#ifndef DIAG_PARALLEL_SOLVER
#define DIAG_PARALLEL_SOLVER

#include <Driver.d/Communicator.h>
#include <Solvers.d/ParallelSolver.h>
#include <iostream>

template<class Scalar> class GenDistrVector;
template<class Scalar> class GenDistrVectorSet;
template<class Scalar> class GenSolver;

class Timings;

template <class Scalar> 
class DiagParallelSolver : public GenParallelSolver<Scalar>
{
    int nSub;
    GenSubDomain<Scalar> **sd;
    GenSolver<Scalar> **solver;
    mutable FSCommPattern<Scalar> *vPat;
    mutable Timings times;
    
    void dispatchInterface(int iSub);
    void assembleDiag(int iSub);
    void dispatchForce(int iSub, const GenDistrVector<Scalar> &) const;
    void assembleForce(int iSub, GenDistrVector<Scalar> &) const;
    void subReSolve(int iSub, GenDistrVector<Scalar>&) const;
    void subSquareRootMult(int iSub, GenDistrVector<Scalar> &v);
    void subInverseSquareRootMult(int iSub, GenDistrVector<Scalar> &v);
  public:
    DiagParallelSolver(int _nSub, GenSubDomain<Scalar> **_sd,
		       GenSolver<Scalar> **_sol, Connectivity *cpuToSub, FSCommunicator *);
    virtual ~DiagParallelSolver() { delete vPat; for(int i=0; i<nSub; ++i) delete solver[i]; delete [] solver; }
    
    void reSolve(GenDistrVector<Scalar> &) override;
    void squareRootMult(GenDistrVector<Scalar> &) override;
    void inverseSquareRootMult(GenDistrVector<Scalar> &) override;
    double getSolutionTime() const override { return times.solve; }
    void solve(const GenDistrVector<Scalar> &, GenDistrVector<Scalar> &) override;
    Timings& getTimers() override { return times; }
};

template<class Scalar>
DiagParallelSolver<Scalar>::DiagParallelSolver(int _nSub, GenSubDomain<Scalar> **_sd,
      GenSolver<Scalar> **_sol, Connectivity *cpuToSub, FSCommunicator *com)
{
 nSub = _nSub;
 sd = _sd;
 solver = _sol;
 
 int iSub;
 // need a FSCommPattern as the DiagParallelSolver is build assembled
 vPat = new FSCommPattern<Scalar>(com, cpuToSub, com->cpuNum(), FSCommPattern<Scalar>::CopyOnSend);
 for(iSub=0; iSub<this->nSub; ++iSub) this->sd[iSub]->setDofPlusCommSize(this->vPat);
 this->vPat->finalize();
 
 execParal(nSub, this, &DiagParallelSolver<Scalar>::dispatchInterface);
 this->vPat->exchange();
 execParal(nSub, this, &DiagParallelSolver<Scalar>::assembleDiag);

 times.solve = 0;
}

template<class Scalar>
void
DiagParallelSolver<Scalar>::reSolve(GenDistrVector<Scalar> &rhsSol)
{
  times.solve -= getTime();
//  execParal(nSub, this, &DiagParallelSolver<Scalar>::dispatchForce, rhsSol);
    threadManager->callParal(nSub, [this, &rhsSol](int s) { dispatchForce(s, rhsSol); });
  this->vPat->exchange();
  execParal(nSub, this, &DiagParallelSolver<Scalar>::assembleForce, rhsSol);

  execParal(nSub, this, &DiagParallelSolver<Scalar>::subReSolve, rhsSol);
  times.solve += getTime();
}

template<class Scalar>
void
DiagParallelSolver<Scalar>::squareRootMult(GenDistrVector<Scalar> &rhsSol)
{
  execParal(nSub, this, &DiagParallelSolver<Scalar>::subSquareRootMult, rhsSol);
}

template<class Scalar>
void
DiagParallelSolver<Scalar>::inverseSquareRootMult(GenDistrVector<Scalar> &rhsSol)
{
  execParal(nSub, this, &DiagParallelSolver<Scalar>::subInverseSquareRootMult, rhsSol);
}

template<class Scalar>
void
DiagParallelSolver<Scalar>::solve(const GenDistrVector<Scalar> &rhs, GenDistrVector<Scalar> &solution)
{
  times.solve -= getTime();
  solution=rhs;
  execParal(nSub, this, &DiagParallelSolver<Scalar>::dispatchForce, solution);
  this->vPat->exchange();
  execParal(nSub, this, &DiagParallelSolver<Scalar>::assembleForce, solution);

  execParal(nSub, this, &DiagParallelSolver<Scalar>::subReSolve, solution);
  times.solve += getTime();
}


template<class Scalar>
void
DiagParallelSolver<Scalar>::dispatchInterface(int iSub)
{
  Scalar *diagVals = new Scalar[solver[iSub]->neqs()];
  GenSparseMatrix<Scalar> *sps = dynamic_cast<GenSparseMatrix<Scalar> * >(solver[iSub]);
  if(sps == 0)
     throw "Using DiagParallelSolver with an inadequate local solver\n";
  for(int i = 0; i < solver[iSub]->neqs(); ++i)
    diagVals[i] = sps->diag(i);
  sd[iSub]->extractAndSendInterf(diagVals, vPat);
  delete [] diagVals;
}

template<class Scalar>
void
DiagParallelSolver<Scalar>::assembleDiag(int iSub)
{
  Scalar *diagVals = new Scalar[solver[iSub]->neqs()];
  for(int i = 0; i < solver[iSub]->neqs(); ++i)
    diagVals[i] = 0;
  sd[iSub]->assembleInterf(diagVals, vPat);
  GenSparseMatrix<Scalar> *sps = dynamic_cast<GenSparseMatrix<Scalar> * >(solver[iSub]);
  int count = 0;
  double small = 1e-5;
  for(int i = 0; i < solver[iSub]->neqs(); ++i) {
    sps->diag(i) += diagVals[i];
    if(sps->diag(i) == 0.0) {
      sps->diag(i) = small;
      count++;
    }
  }
  if(count > 0) std::cerr << " *** WARNING: " << count << " zero diagonal/s detected in subdomain mass matrix set to " << small << std::endl;
  delete [] diagVals;
}

template<class Scalar>
void
DiagParallelSolver<Scalar>::dispatchForce(int iSub, const GenDistrVector<Scalar> &v) const
{
  const Scalar *subV = v.subData(iSub);
  sd[iSub]->extractAndSendInterf(subV, vPat);
}

template<class Scalar>
void
DiagParallelSolver<Scalar>::assembleForce(int iSub, GenDistrVector<Scalar> &v) const
{
  Scalar *subV = v.subData(iSub);
  sd[iSub]->assembleInterf(subV, vPat);
}

template<class Scalar>
void
DiagParallelSolver<Scalar>::subReSolve(int iSub, GenDistrVector<Scalar> &v) const
{
  solver[iSub]->reSolve(v.subData(iSub));
}

template<class Scalar>
void
DiagParallelSolver<Scalar>::subSquareRootMult(int iSub, GenDistrVector<Scalar> &v)
{
  dynamic_cast<GenSparseMatrix<Scalar>*>(solver[iSub])->squareRootMult(v.subData(iSub));  
}

template<class Scalar>
void
DiagParallelSolver<Scalar>::subInverseSquareRootMult(int iSub, GenDistrVector<Scalar> &v)
{
  dynamic_cast<GenSparseMatrix<Scalar>*>(solver[iSub])->inverseSquareRootMult(v.subData(iSub));
}

#endif
