#ifndef _EIGENPROBTYPE_H_
#define _EIGENPROBTYPE_H_

//HB: I have noticed that the eigen vectors are store twice which can cause early memory 
//    limitation for large model and/or large number of eigen values/vectors requested.
//    For instance, the eigen/Ritz vectors are first computed and store in the "Z vector set and
//    then copied to their final destination vector set eigVec. We may be able to use/re-use the 
//    "Z vector" set for storing the eigen vectors -> reduction by 2 of the amount of memory required 

template <class Scalar> class GenVector; 
typedef class GenVector<double> Vector;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;

template <class EigOps,
          class VecType,
          class VecSet,
          class PostProcessor,
          class ProblemDescriptor>
class EigenSolver {

  protected:

    int numEig;
    int nrmod;
    int totalEig;
    int origSubSize;
 
    Vector  *eigVal;
    VecSet  *eigVec;

    EigOps  *eM;

    PostProcessor     *postProcessor;

    ProblemDescriptor *probDesc;

  public:

    EigenSolver()                      { probDesc=0; postProcessor = 0; eigVal = 0; eigVec = 0; eM = 0; }
    EigenSolver(ProblemDescriptor *p)  { probDesc=p; postProcessor = 0; eigVal = 0; eigVec = 0; eM = 0; }

//    static EigenSolver * buildEigenSolver(ProblemDescriptor *p);

    Vector  * getpeigval() { return eigVal ; }
    VecSet  & getpeigvec() { return *eigVec; }
    EigOps  & geteM()      { return *eM ;    }

    void getJacobi(double *kappa, double * mu, FullSquareMatrix &xx,
                   double *eigVal, int nsmax, int subSpaceSize, double tolJac);
    void ortho(VecType *v1, VecType *vr, int nsub, int nrbm);
    void ortho(VecSet& v1, VecSet& vr, int nsub, int nrbm);
    void absoluteInnerproductNormalized(const VecType& v1, const VecType& v2, double &result);
    void pickMostCorrelatedModes(Vector &, VecSet &);
    void setUp();
//    void cleanup();

    virtual void initialize()=0;
    virtual void solve()=0;
    void performQR(Vector *,VecSet *, int);
};

template <class EigOps,
          class VecType,
          class VecSet,
          class PostProcessor,
          class ProblemDescriptor>
class SubSpaceSolver: public EigenSolver <EigOps,
                                          VecType,
                                          VecSet,
                                          PostProcessor,
                                          ProblemDescriptor>
{
  private:
    int subSpaceSize;
    int nsmax;
    int nsub;

    double tolEig;
    double tolJac;
    bool explicitK;

    VecSet *Q;
    VecSet *Z;

    Vector *subVal;
    Vector *subOld;

  public:
    SubSpaceSolver()                      {this->probDesc=0; Q = 0; Z = 0; subVal = 0; subOld = 0;}
    SubSpaceSolver(ProblemDescriptor *p)  {this->probDesc=p; Q = 0; Z = 0; subVal = 0; subOld = 0;}

    void initialize();
    void solve();

};

template <class EigOps,
          class VecType,
          class VecSet,
          class PostProcessor,
          class ProblemDescriptor>
class SymArpackSolver : public EigenSolver <EigOps,
                                           VecType,
                                           VecSet,
                                           PostProcessor,
                                           ProblemDescriptor>
{
  private:
    int subSpaceSize;
    int nsmax;
    int nsub;

    double tolEig;
    double tolJac;
    bool explicitK;

    VecSet* Q; // Lanczos vectors
    VecSet* Z; // Ritz vectors

  public:
    SymArpackSolver()                      {this->probDesc=0; Q = 0; Z = 0; }
    SymArpackSolver(ProblemDescriptor *p)  {this->probDesc=p; Q = 0; Z = 0; }

    void initialize() {};
    void solve();
    void rebuildSolver(double frequency);
};

template <class EigOps,
          class VecType,
          class VecSet,
          class PostProcessor,
          class ProblemDescriptor>
class LOBPCGSolver: public EigenSolver <EigOps,
                                       VecType,
                                       VecSet,
                                       PostProcessor,
                                       ProblemDescriptor>
{
  private:
    int subSpaceSize;
    int nsmax;
    int nsub;

    double tolEig;
    double tolJac;
    bool explicitK;

    VecSet* Q; 
    VecSet* Z; 

    Vector *subVal;

  public:
    LOBPCGSolver()                      {this->probDesc=0; Q = 0; Z = 0; subVal = 0; }
    LOBPCGSolver(ProblemDescriptor *p)  {this->probDesc=p; Q = 0; Z = 0; subVal = 0; }

    void initialize();
    void solve();
};


#ifdef _TEMPLATE_FIX_
#include <Driver.d/EigenProbType.C>
#endif

#endif
