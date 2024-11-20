#ifndef _EIGENDESCR_H_
#define _EIGENDESCR_H_
#include <Driver.d/Domain.h>

//class Domain;
template <class Scalar> class GenDynamMat;
typedef GenDynamMat<double> DynamMat;
template <class Scalar> class GenVector;
typedef GenVector<double> Vector;
template <class Scalar> class GenVectorSet;
typedef GenVectorSet<double> VectorSet;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
template <class Scalar> class GenDBSparseMatrix;
typedef GenDBSparseMatrix<double> DBSparseMatrix;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
class StaticTimers;
class Corotator;
class GeomState;


// Single Domain Eigenvalue Problem Post Processor

class SDEigenPostProcessor {
    Domain *domain;
    StaticTimers *times;
    double *bcx;
  public:
    // Constructor
    SDEigenPostProcessor(Domain *_domain,StaticTimers *_times, double *_bcx=0) 
    { domain = _domain; times = _times; bcx = _bcx;}

    void eigenOutput(Vector& eigvalues, VectorSet& eigVectors, int convEig = 0);
#ifdef USE_EIGEN3
    void eigenQROutput(Eigen::MatrixXd& X, Eigen::MatrixXd& Q, Eigen::MatrixXd& R);
#endif
};


// Single Domain Eigenvale Problem Descriptor

class SingleDomainEigen {

    Domain       *domain;
    StaticTimers *times;

    // members for nonlinear eigen problem
    FullSquareMatrix *kelArray;
    FullSquareMatrix *melArray;
    FullSquareMatrix *geomKelArray; // For Buckling analysis

    Corotator **allCorot;
    GeomState *geomState;

    double *bcx;

 public:

    // Constructor
    SingleDomainEigen(Domain *d) { domain = d; kelArray = 0; melArray = 0; 
                                   geomState = 0; allCorot = 0; geomKelArray = 0; }

    SDEigenPostProcessor *getPostProcessor();

    FullSquareMatrix * getkelArray() { return kelArray; }

    int  solVecSize();
    int  solVecInfo();
    int  getNumEigen();
    int  getEigenSolverType();
    int  getEigenSolverSubType();
    bool getFilter() { return domain->solInfo().readmodeCalled; }

    void preProcess();
    void getEigenInfo(int& subspacesize, int& numEig, double& tolEig, 
                      double& tolJac);
    void buildEigOps( DynamMat &eM );
    void reBuild(DynamMat &eM);//CBM
    void error(int subspacesize, int numrbm); 
    void initQ(Vector* Z, int nsub);
    void printTimers(Solver *solver);

    void getSubSpaceInfo(int& subspacesize, int& maxIter, 
                         double& tolEig, double& tolJac, bool &explicitK);
    void convertModeDataToVecSet(VectorSet& vModeData);
};

#endif
