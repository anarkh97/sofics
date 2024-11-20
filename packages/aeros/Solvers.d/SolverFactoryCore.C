#include <Solvers.d/SolverFactory.h>

template<>
GenSolverFactory<double>* GenSolverFactory<double>::getFactory() { return solverFactory.get(); }

template<>
GenSolverFactory<DComplex>* GenSolverFactory<DComplex>::getFactory() { return solverFactoryC.get(); }

