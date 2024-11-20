#ifndef __SOLVER_FACTORY_H__
#define __SOLVER_FACTORY_H__

#include <Solvers.d/Solver.h>
#include <Solvers.d/SolverCntl.h>
#include <map>
#include <memory>
#include <string>

class EqNumberer;
class DofSetArray;
class ConstrainedDSA;
template <class Scalar> class GenSparseMatrix;
class FSCommunicator;

template <typename Scalar>
struct SolverAndMatrix {
	GenSolver<Scalar> *solver;
	std::unique_ptr<GenSparseMatrix<Scalar>> sparseMatrix;
};

template<class Scalar>
class GenSolverFactory
{
 public:
  virtual ~GenSolverFactory() {}

  // used for feti coarse solvers: sequential or parallel direct solver only
  virtual GenSolver<Scalar>* createSolver(const Connectivity *con, const EqNumberer *eqnum, const SolverCntl &cntl,
                                          GenSparseMatrix<Scalar> *&sparse, int ngrbm = 0,
                                          FSCommunicator *com = 0, std::string name = "") const;
  virtual GenSolver<Scalar>* createDistSolver(const Connectivity *con, const EqNumberer *eqnum, const SolverCntl &cntl,
                                              GenSparseMatrix<Scalar> *&sparse,
                                              FSCommunicator *com = 0, std::string name = "") const;
  // used for global solvers and feti local solvers: direct or iterative solver
  virtual GenSolver<Scalar>* createSolver(const Connectivity *con, const DofSetArray *dsa, const ConstrainedDSA *cdsa,
                                          const SolverCntl &cntl, GenSparseMatrix<Scalar> *&sparse, Rbm *rbm,
                                          GenSparseMatrix<Scalar> *&spp, GenSolver<Scalar> *&prec,
                                          FSCommunicator *com = 0, std::string name = "") const;
	// used for feti Kii solver: sequential direct solver only
	virtual SolverAndMatrix<Scalar> createSolver(const Connectivity *con, const DofSetArray *dsa, int *map,
	                                             const SolverCntl &cntl,
	                                             std::string name = "") const;

  static GenSolverFactory* getFactory();

};

extern std::unique_ptr<GenSolverFactory<double> >   solverFactory;
extern std::unique_ptr<GenSolverFactory<DComplex> > solverFactoryC;

#endif

