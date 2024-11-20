#ifdef USE_EIGEN3
#include "EiGalerkinProjectionSolver.C"

namespace Rom {

template
GenEiSparseGalerkinProjectionSolver<double>
::GenEiSparseGalerkinProjectionSolver(const Connectivity*, const DofSetArray*, const ConstrainedDSA*, bool, double);

template
GenEiSparseGalerkinProjectionSolver<complex<double> >
::GenEiSparseGalerkinProjectionSolver(const Connectivity*, const DofSetArray*, const ConstrainedDSA*, bool, double);

template
void
GenEiSparseGalerkinProjectionSolver<double>
::zeroAll();

template
void
GenEiSparseGalerkinProjectionSolver<complex<double> >
::zeroAll();

template
void
GenEiSparseGalerkinProjectionSolver<double>
::addReducedMass(double);

template
void
GenEiSparseGalerkinProjectionSolver<complex<double> >
::addReducedMass(double);

template
void
GenEiSparseGalerkinProjectionSolver<double>
::addToReducedMatrix(const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &, double);

template
void
GenEiSparseGalerkinProjectionSolver<complex<double> >
::addToReducedMatrix(const Eigen::Matrix<complex<double>,Eigen::Dynamic,Eigen::Dynamic> &, double);

template
void
GenEiSparseGalerkinProjectionSolver<double>
::addLMPCs(int, LMPCons**, double);

template 
void
GenEiSparseGalerkinProjectionSolver<complex<double> >
::addLMPCs(int, LMPCons**, double);

template
void 
GenEiSparseGalerkinProjectionSolver<double>
::updateLMPCs(GenVector<double> &);

template 
void
GenEiSparseGalerkinProjectionSolver<complex<double> >
::updateLMPCs(GenVector<complex<double> > &);

template
void
GenEiSparseGalerkinProjectionSolver<double>
::projectionBasisIs(GenVecBasis<double>&);

template
void
GenEiSparseGalerkinProjectionSolver<complex<double> >
::projectionBasisIs(GenVecBasis<complex<double> >&);

template
void
GenEiSparseGalerkinProjectionSolver<double>
::EmpiricalSolver();

template
void
GenEiSparseGalerkinProjectionSolver<complex<double> >
::EmpiricalSolver();

template
void
GenEiSparseGalerkinProjectionSolver<double>
::factor();

template
void
GenEiSparseGalerkinProjectionSolver<complex<double> >
::factor();
   
template
void
GenEiSparseGalerkinProjectionSolver<double>
::reSolve(GenVector<double>&);

template
void
GenEiSparseGalerkinProjectionSolver<complex<double> >
::reSolve(GenVector<complex<double> >&);

template
void
GenEiSparseGalerkinProjectionSolver<double>
::solve(const GenVector<double>&, GenVector<double>&);

template
void
GenEiSparseGalerkinProjectionSolver<complex<double> >
::solve(const GenVector<complex<double> >&, GenVector<complex<double> >&);

template
double
GenEiSparseGalerkinProjectionSolver<double>
::getResidualNorm(const GenVector<double> &v);

template 
double
GenEiSparseGalerkinProjectionSolver<complex<double> >
::getResidualNorm(const GenVector<complex<double> > &v);

} /* end namespace Rom */
#endif
