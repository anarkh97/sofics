#ifndef _MD_EIGEN_DESCR_H_
#define _MD_EIGEN_DESCR_H_
#include <Driver.d/DecDomain.h>

template <class Scalar> class GenDecDomain;
typedef GenDecDomain<double> DecDomain;
template <class Scalar> class GenDynamMat;
typedef GenDynamMat<double> DynamMat;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
template <class Scalar> class GenVectorSet;
typedef GenVectorSet<double> VectorSet;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
class StaticTimers;
template <class Scalar> class GenParallelSolver;
typedef GenParallelSolver<double> ParallelSolver;
template <class Scalar> class GenSolver;
typedef GenSolver<double> Solver;
template <class Scalar> class GenParallelSolver;

//CBM
template <class Scalar> class GenDistrVectorSet;
typedef GenDistrVectorSet<double> DistrVectorSet;
template <class Scalar> class GenDistrVector;
typedef GenDistrVector<double> DistrVector;
class DistrGeomState;
class DistrInfo;


// Multiple Domain Eigenvalue Problem Post Processor
template <class Scalar>
class GenMultiDomainEigenPostProcessor
{
    GenDecDomain<Scalar> *decDomain;
    GenParallelSolver<Scalar> *solver;
    StaticTimers *times;

  public:
    GenMultiDomainEigenPostProcessor(GenDecDomain<Scalar> *d, GenParallelSolver<Scalar> *s,
                                StaticTimers* _times=0)
       { decDomain = d; solver = s; times = _times; }

    void eigenOutput(GenVector<Scalar>& eigValues, GenDistrVectorSet<Scalar>& eigVectors, int convEig = 0);
};

// Multiple Domain Eigenvale Problem Descriptor

template <class Scalar>
class GenMultiDomainEigen 
{
    Domain *domain;
    GenDecDomain<Scalar> *decDomain;
    GenParallelSolver<Scalar> *solver;
    StaticTimers *times;

 public:

    // Constructor
    GenMultiDomainEigen(Domain *d);
//    ~GenMultiDomainEigen() { delete decDomain; delete times; } 
    ~GenMultiDomainEigen() { delete times; }

    int solVecSize();
    DistrInfo &solVecInfo();
    void preProcess();
    bool getFilter() { return domain->solInfo().readmodeCalled; }

    GenMultiDomainEigenPostProcessor<Scalar> *getPostProcessor(); 

    void buildEigOps(MDDynamMat &dMat);
    void reBuild(MDDynamMat &dMat);
    void error(int subspacesize, int numrbm); 
    void printTimers(GenParallelSolver<double> *solver);
    int  getNumEigen();
    int  getEigenSolverType();
    void getSubSpaceInfo(int& subspacesize, int& maxIter,
                         double& tolEig, double& tolJac, bool &explicitK);
    void getMemoryK(int iSub, long *memory);
    void getMemoryPrec(int iSub, long *memory);
    void convertModeDataToVecSet(DistrVectorSet& vModeData) { filePrint(stderr," ... GenMultiDomainEigen::convertModeDataToVecSet has not been implemented.\n"); }
};

typedef GenMultiDomainEigen<double> MultiDomainEigen;
typedef GenMultiDomainEigenPostProcessor<double> MultiDomainEigenPostProcessor;

#ifdef _TEMPLATE_FIX_
#include <Paral.d/MDEigen.C>
#endif

#endif
