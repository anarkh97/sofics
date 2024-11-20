#include <Driver.d/OpMake.C>

#define OPMAKE_INSTANTIATION_HELPER(Scalar) \
template \
void \
Domain::makeSparseOps<Scalar>(AllOps<Scalar>&, double, double, double, GenSparseMatrix<Scalar>*,\
                              FullSquareMatrix*, FullSquareMatrix*, FullSquareMatrix*);\
template \
GenDBSparseMatrix<Scalar> * \
Domain::constructDBSparseMatrix<Scalar>(DofSetArray*, Connectivity*);\
\
template \
GenNBSparseMatrix<Scalar> * \
Domain::constructNBSparseMatrix<Scalar>();\
\
template \
GenCuCSparse<Scalar> * \
Domain::constructCuCSparse<Scalar>(DofSetArray*);\
\
template \
GenCuCSparse<Scalar> * \
Domain::constructCCSparse<Scalar>(DofSetArray*);\
\
template \
GenBLKSparseMatrix<Scalar> * \
Domain::constructBLKSparseMatrix<Scalar>(DofSetArray*, Rbm*);\
\
template \
void \
Domain::buildPreSensitivities<Scalar>(AllSensitivities<Scalar>&, Scalar*);\
\
template \
void \
Domain::buildPostSensitivities<Scalar>(GenSolver<Scalar>*, GenSparseMatrix<Scalar>*,\
                                       GenSparseMatrix<Scalar>*, AllSensitivities<Scalar>&,\
                                       GenVector<Scalar>*, Scalar*, bool, GeomState*, GeomState*, Corotator **, bool);\
\
template \
void \
Domain::buildNLPostSensitivities<Scalar>(GenSolver<Scalar>*, AllSensitivities<Scalar>&,\
                                         GeomState*, GeomState*, Corotator **, bool);\
\
template \
void \
Domain::buildOps<Scalar>(AllOps<Scalar>&, double, double, double,\
                         Rbm*, FullSquareMatrix*, FullSquareMatrix*,\
                         FullSquareMatrix*, bool);\
\
template \
void \
Domain::rebuildOps<Scalar>(AllOps<Scalar>&, double, double, double,\
                           Rbm*, FullSquareMatrix*, FullSquareMatrix*,\
                           FullSquareMatrix*, bool);\
\
template \
void \
Domain::getSolverAndKuc<Scalar>(AllOps<Scalar>&, FullSquareMatrix*, Rbm*, bool);\
\
template \
void \
Domain::makeStaticOpsAndSolver<Scalar>(AllOps<Scalar>&, double, double,\
                                       double, GenSolver<Scalar>*&, GenSparseMatrix<Scalar>*&,\
                                       Rbm*, FullSquareMatrix*, FullSquareMatrix*, FullSquareMatrix*);\
\
template \
void \
Domain::makeDynamicOpsAndSolver<Scalar>(AllOps<Scalar>&, double, double,\
                                        double, GenSolver<Scalar>*&, GenSparseMatrix<Scalar>*&,\
                                        Rbm*, FullSquareMatrix*, FullSquareMatrix*, FullSquareMatrix*);\
\
template \
void \
Domain::addGravityForce<Scalar>(GenVector<Scalar>&);\
\
template \
void \
Domain::addGravityForceSensitivity<Scalar>(GenVector<Scalar>&);\
\
template \
void \
Domain::addPressureForce<Scalar>(GenVector<Scalar>&, int, double);\
\
template \
void \
Domain::addAtddnbForce<Scalar>(GenVector<Scalar>&, int, double);\
\
template \
void \
Domain::addAtdrobForce<Scalar>(GenVector<Scalar>&, int, double);\
\
template \
void \
Domain::addThermalForce<Scalar>(GenVector<Scalar>&);\
\
template \
void \
Domain::addMpcRhs<Scalar>(GenVector<Scalar>&, double);\
\
template \
void \
Domain::scaleDisp<Scalar>(Scalar*);\
\
template \
void \
Domain::scaleInvDisp<Scalar>(Scalar*);\
\
template \
void \
Domain::scaleDisp<Scalar>(Scalar*, double);\
\
template \
int Domain::mergeDistributedDisp<Scalar>(Scalar (*)[11], Scalar*, Scalar*, Scalar (*)[11]);\
\
template \
void Domain::forceDistributedContinuity<Scalar>(Scalar*, Scalar (*)[11]);\
\
template \
void \
Domain::buildRHSForce<Scalar>(GenVector<Scalar>&, GenSparseMatrix<Scalar>*);\
\
template \
void \
Domain::computeReactionForce<Scalar>(GenVector<Scalar>&, GenVector<Scalar>&,\
                                     GenSparseMatrix<Scalar>*, GenSparseMatrix<Scalar>*);\
\
template \
void \
Domain::buildRHSForce<Scalar>(GenVector<Scalar>&, GenVector<Scalar>&,\
                              GenSparseMatrix<Scalar>*, GenSparseMatrix<Scalar>*, \
                              GenSparseMatrix<Scalar>**, \
                              GenSparseMatrix<Scalar>**, \
                              GenSparseMatrix<Scalar>**, \
                              GenSparseMatrix<Scalar>**, \
                              double, double, GeomState*);\
\
template \
void \
Domain::buildFreqSweepRHSForce<Scalar>(GenVector<Scalar>&, GenSparseMatrix<Scalar>*,\
                                       GenSparseMatrix<Scalar>**, \
                                       GenSparseMatrix<Scalar>**, \
                                       int, double);\
\
template \
void \
Domain::buildDeltaK(double w0, double w, GenSparseMatrix<Scalar> *deltaK, \
                                         GenSparseMatrix<Scalar> *deltaKuc); \
\
template \
void \
Domain::assembleATDROB<Scalar>(GenSparseMatrix<Scalar>*, AllOps<Scalar>*, double);\
\
template \
void \
Domain::assembleSommer<Scalar>(GenSparseMatrix<Scalar>*, AllOps<Scalar>*);\
\
template \
void \
Domain::computeSommerDerivatives<Scalar>(double, double, int, int*, FullSquareMatrix&,\
                                         DComplex**, double, double, AllOps<Scalar>*);\
\
template \
void \
Domain::updateMatrices<Scalar>(AllOps<Scalar>*, GenSparseMatrix<Scalar>*, int*, int*,\
                               FullSquareMatrix*, FullSquareMatrix*, double);\
\
template \
void \
Domain::updateDampingMatrices<Scalar>(AllOps<Scalar>*, int*, FullSquareMatrix*,\
                                      FullSquareMatrix*, double, int);\
\
template \
int \
Domain::processDispTypeOutputs<Scalar>(OutputInfo&, Scalar (*)[11], int,\
                                       int, double, double, int);\
\
template \
int \
Domain::processOutput<Scalar>(OutputInfo::Type&, GenVector<Scalar>&, Scalar*, int,\
                              double, double, int);\
\
template \
void \
Domain::postProcessing<Scalar>(GenVector<Scalar>&, Scalar*, GenVector<Scalar>&,\
                               int, int, double, double,\
                               GenSparseMatrix<Scalar>*, GenSparseMatrix<Scalar>*);\
\
template \
void \
Domain::addConstantForceSensitivity<Scalar>(GenVector<Scalar>&, GenSparseMatrix<Scalar>*);\
\
template \
void \
Domain::computeConstantForce<Scalar>(GenVector<Scalar>&, GenSparseMatrix<Scalar>*);\
\
template \
void \
Domain::computeExtForce<Scalar>(GenVector<Scalar>&, double, GenSparseMatrix<Scalar>*,\
                                ControlInterface*, GenSparseMatrix<Scalar>*, \
                                double, GenSparseMatrix<Scalar>*);\
\
template \
void \
Domain::computeExtForce4<Scalar>(GenVector<Scalar>&, const GenVector<Scalar>&,\
                                 double, GenSparseMatrix<Scalar>*, ControlInterface*,\
                                 GenSparseMatrix<Scalar>*, double, GenSparseMatrix<Scalar>*);\

OPMAKE_INSTANTIATION_HELPER(double)
OPMAKE_INSTANTIATION_HELPER(std::complex<double>)

#ifdef USE_EIGEN3
template
GenEiSparseMatrix<double,Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Upper> > *
Domain::constructEiSparse<double>(DofSetArray*, Connectivity*, bool);

template
GenEiSparseMatrix<std::complex<double>,Eigen::SimplicialLLT<Eigen::SparseMatrix<std::complex<double> >,Eigen::Upper> > *
Domain::constructEiSparse<std::complex<double> >(DofSetArray*, Connectivity*, bool);

template
GenEiSparseMatrix<double,Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Upper> > *
Domain::constructEiSparseMatrix<double,Eigen::SimplicialLLT<Eigen::SparseMatrix<double>,Eigen::Upper> >(DofSetArray*, Connectivity*, bool);

template 
GenEiSparseMatrix<std::complex<double>,Eigen::SimplicialLLT<Eigen::SparseMatrix<std::complex<double> >,Eigen::Upper> > * 
Domain::constructEiSparseMatrix<std::complex<double>,Eigen::SimplicialLLT<Eigen::SparseMatrix<std::complex<double> >,Eigen::Upper> >(DofSetArray*, Connectivity*, bool);

template 
GenEiSparseMatrix<double,Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>,Eigen::Upper> > * 
Domain::constructEiSparseMatrix<double,Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>,Eigen::Upper> >(DofSetArray*, Connectivity*, bool);

template
GenEiSparseMatrix<std::complex<double>,Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double> >,Eigen::Upper> > *
Domain::constructEiSparseMatrix<std::complex<double>,Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double> >,Eigen::Upper> >(DofSetArray*, Connectivity*, bool);

template 
GenEiSparseMatrix<double,Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int> > > * 
Domain::constructEiSparseMatrix<double,Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int> > >(DofSetArray*, Connectivity*, bool);

template
GenEiSparseMatrix<std::complex<double>,Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double> >,Eigen::COLAMDOrdering<int> > > *
Domain::constructEiSparseMatrix<std::complex<double>,Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double> >,Eigen::COLAMDOrdering<int> > >(DofSetArray*, Connectivity*, bool);

template 
GenEiSparseMatrix<double,Eigen::SparseQR<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int> > > * 
Domain::constructEiSparseMatrix<double,Eigen::SparseQR<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int> > >(DofSetArray*, Connectivity*, bool);

template
GenEiSparseMatrix<std::complex<double>,Eigen::SparseQR<Eigen::SparseMatrix<std::complex<double> >,Eigen::COLAMDOrdering<int> > > *
Domain::constructEiSparseMatrix<std::complex<double>,Eigen::SparseQR<Eigen::SparseMatrix<std::complex<double> >,Eigen::COLAMDOrdering<int> > >(DofSetArray*, Connectivity*, bool);

#ifdef EIGEN_CHOLMOD_SUPPORT
template 
GenEiSparseMatrix<double,Eigen::CholmodDecomposition<Eigen::SparseMatrix<double>,Eigen::Upper> > * 
Domain::constructEiSparseMatrix<double,Eigen::CholmodDecomposition<Eigen::SparseMatrix<double>,Eigen::Upper> >(DofSetArray*, Connectivity*, bool);

template
GenEiSparseMatrix<std::complex<double>,Eigen::CholmodDecomposition<Eigen::SparseMatrix<std::complex<double> >,Eigen::Upper> > *
Domain::constructEiSparseMatrix<std::complex<double>,Eigen::CholmodDecomposition<Eigen::SparseMatrix<std::complex<double> >,Eigen::Upper> >(DofSetArray*, Connectivity*, bool);
#endif

#ifdef EIGEN_UMFPACK_SUPPORT
template 
GenEiSparseMatrix<double,Eigen::UmfPackLU<Eigen::SparseMatrix<double> > > * 
Domain::constructEiSparseMatrix<double,Eigen::UmfPackLU<Eigen::SparseMatrix<double> > >(DofSetArray*, Connectivity*, bool);

template
GenEiSparseMatrix<std::complex<double>,Eigen::UmfPackLU<Eigen::SparseMatrix<std::complex<double> > > > *
Domain::constructEiSparseMatrix<std::complex<double>,Eigen::UmfPackLU<Eigen::SparseMatrix<std::complex<double> > > >(DofSetArray*, Connectivity*, bool);
#endif

#ifdef EIGEN_SPQR_SUPPORT
template 
GenEiSparseMatrix<double,Eigen::SPQR<Eigen::SparseMatrix<double> > > * 
Domain::constructEiSparseMatrix<double,Eigen::SPQR<Eigen::SparseMatrix<double> > >(DofSetArray*, Connectivity*, bool);

template
GenEiSparseMatrix<std::complex<double>,Eigen::SPQR<Eigen::SparseMatrix<std::complex<double> > > > *
Domain::constructEiSparseMatrix<std::complex<double>,Eigen::SPQR<Eigen::SparseMatrix<std::complex<double> > > >(DofSetArray*, Connectivity*, bool);
#endif

#ifdef EIGEN_SUPERLU_SUPPORT
template 
GenEiSparseMatrix<double,Eigen::SuperLU<Eigen::SparseMatrix<double> > > * 
Domain::constructEiSparseMatrix<double,Eigen::SuperLU<Eigen::SparseMatrix<double> > >(DofSetArray*, Connectivity*, bool);

template
GenEiSparseMatrix<std::complex<double>,Eigen::SuperLU<Eigen::SparseMatrix<std::complex<double> > > > *
Domain::constructEiSparseMatrix<std::complex<double>,Eigen::SuperLU<Eigen::SparseMatrix<std::complex<double> > > >(DofSetArray*, Connectivity*, bool);
#endif
#endif
