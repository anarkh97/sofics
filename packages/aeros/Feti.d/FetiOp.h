#ifndef _FETIOP_H_
#define _FETIOP_H_

#include <Threads.d/Paral.h>

template <class Scalar> class GenFetiOpControler;
template <class Scalar> class GenFetiSolver;
template <class Scalar> class FetiBaseClass;
class IntFullM;
template <class Scalar> class GenSubDomain;
template <class Scalar> class GenSparseMatrix;

template<class Scalar>
class GenFetiOp : public TaskDescr 
{
   protected:
    GenSubDomain<Scalar> *sd = nullptr;
    GenSolver<Scalar> *solver = nullptr;
    GenSolver<Scalar> *K = nullptr;
    GenSparseMatrix<Scalar> *KasSparse = nullptr;
    Rbm *rbm = nullptr; // For dynamics only and nonlinear
    int numNeighb;
    int halfOffset;
    std::vector<int> alphaOffset;
	std::vector<int> betaOffset;
    int isFeti2, isDynamic;

    int numRBM = 0;            // total number of RBMs
    std::vector<double> locRBMs ;       // local RBMs for the whole domain
    std::vector<double> locInterfRBMs ; // local RBMs on the interface in local linear form
                           // The modes are consecutive in the array

    std::vector<int> neighbNumRBMs;
//    double **neighbRBMs = nullptr;

    int crnDofSize;        // number of additional lagrange multipliers
    IntFullM *BClocal = nullptr;

    GenFetiOpControler<Scalar> *control = nullptr;

    int QGisLocal;  // Wether or not QG is local

    Scalar *interfBuff = nullptr;
    FSCommPattern<Scalar> *vPat = nullptr;

   public:
    GenFetiOp() {}
    GenFetiOp(GenSubDomain<Scalar> *, GenFetiOpControler<Scalar> *, 
              int, int, FSCommPattern<Scalar> *, Rbm * =0);
    virtual ~GenFetiOp();
    void setSysMatrix(GenSolver<Scalar> *k, GenSparseMatrix<Scalar> *ks) { K = k; KasSparse = ks; }
    void run() override;
    void runFor(int) override { throw "Illegal operation called on GenFetiOp"; }
    void clean_up();

    void localSolve();
    void sendInterfRBM(FSCommPattern<Scalar> *rbmPat);
    void sendNumNeighbRBM(FSCommPattern<int> *sPat);
    void getNumNeighbRBM(FSCommPattern<int> *sPat);
    void getNeighbQGs(FSCommPattern<Scalar> *rbmPat);
    void getGtMult();
    void getGtQMult();
    void getGtFMult();
    void subAlphaG();
    void subAlphaG1();
    void subAlphaG2();
    void subAlphaGQ();
    void subNuC();
    void initializeCRNs(FSCommPattern<int> *sPat);
    void assembleGtCs();
    void getNumNeighbCRNs(FSCommPattern<int> *sPat);
    void assembleGtQGs();
    void makeCoarseSet();
    void computeFiBC();
    void getNeighbFC();
    void assembleGtFCs();
    void assembleCtFCs();
    void computeNeighborFGs();
    void getCtFMult();
    void getCtMult();
    void subNuFC();
    void subAlphaFG();
    void setHalfOffset(int a) { halfOffset = a; }
    void reSetAlphaOffsets(int *v);
    void setAlphaOffsets(int *v);
    void setBetaOffsets(int *v);

    int  getNumRBM()     { return numRBM;      }
    int  getcrnDofSize() { return crnDofSize;  }
    int  numEdge()       { return numNeighb;   } 

    void setglobalSum(void *);
    double res;
	friend class GenFetiSolver<Scalar>;
	friend class FetiBaseClass<Scalar>;

};

typedef GenFetiOp<double> FetiOp;

#endif


