#ifndef _MD_NL_STATIC_H_
#define _MD_NL_STATIC_H_

template <class Scalar> class GenDecDomain;
typedef GenDecDomain<double> DecDomain;
class Domain;
template <class Scalar> class GenParallelSolver;
typedef GenParallelSolver<double> ParallelSolver;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
class Corotator;
class StaticTimers;
class DistrInfo;
template <class Scalar> class GenDistrVector;
typedef GenDistrVector<double> DistrVector;
class DistrGeomState;
template <class Scalar> class GenMultiDomainPostProcessor;
typedef GenMultiDomainPostProcessor<double> MultiDomainPostProcessor;
template <class Scalar> class GenFetiDPSolver;

class MDNLStatic 
{
    Domain     *domain;
    DecDomain  *decDomain;
    ParallelSolver *solver;
    GenFetiDPSolver<double> *fetiSolver;

    FullSquareMatrix **kelArray;
    Corotator ***allCorot;

    DistrVector *reactions;

    double firstDv;
    double firstRes;
    double tolerance;
    StaticTimers *times;
    int numSystems;
    double deltaLambda;

    std::map<std::pair<int,int>, double> *mu; // lagrange multipliers for the contact surfaces
    std::vector<double> *lambda; // lagrange multipliers for all the other constraints

    bool myDecDomain;
    bool updateCS;

 public:
    // Constructor
    MDNLStatic(Domain *d, DecDomain *dd=0);
    virtual ~MDNLStatic();

    DistrInfo& solVecInfo();
    DistrInfo& sysVecInfo();
    DistrInfo& elemVecInfo();
    int checkConvergence(int iter, double normDv, double normRes);
    void updateContactSurfaces(DistrGeomState& geomState, DistrGeomState *refState);
    void updateStates(DistrGeomState *refState, DistrGeomState& geomState, double lambda);
    int  getMaxit();
    double getScaleFactor();
    double getDeltaLambda0();
    double getMaxLambda();
    virtual void getRHS(DistrVector &);
    ParallelSolver *getSolver();

    void printTimers();

    void staticOutput(DistrGeomState *geomState, double lambda, 
                      DistrVector &force, DistrVector &glRes, DistrGeomState *refState);

    MultiDomainPostProcessor *getPostProcessor();

    void preProcess(bool factor = true);

    int reBuild(int iter, int step, DistrGeomState& geomState);

    virtual DistrGeomState* createGeomState();

    void updatePrescribedDisplacement(DistrGeomState *geomState, double l=1.0);

    void initializeParameters(int step, DistrGeomState *geomState);
    void updateParameters(DistrGeomState *geomState);
    bool checkConstraintViolation(double &err, DistrGeomState *geomState);

    double getStiffAndForce(DistrGeomState& geomState, DistrVector& residual, 
                            DistrVector& elementInternalForce, DistrVector& gRes,
                            double lambda = 1.0, DistrGeomState *refState = NULL, bool forceOnly = false);

    double getTolerance();

    LinesearchInfo& linesearch();
    double getEnergy(double lambda, DistrVector& force, DistrGeomState* geomState)
      { std::cerr << "MDNLStatic::getEnergy is not implemented\n"; return 0; }

    double getResidualNorm(DistrVector &vec, DistrGeomState &geomState);
    bool getResizeFlag();

  private:
    void getSubStiffAndForce(int isub, DistrGeomState &geomState,
                             DistrVector &res, DistrVector &elemIntForce, double lambda,
                             DistrGeomState *refState, bool forceOnly);

    void makeSubCorotators(int isub);
    void deleteSubCorotators(int isub);
    void makeSubKelArrays(int isub);
    void deleteSubKelArrays(int isub);
    void makeSubDofs(int isub);
    void updatePrescribedDisp(int isub, DistrGeomState& geomState);
    void subGetRHS(int isub, DistrVector& rhs);
    void addConstraintForces(int isub, DistrVector &vec);
    void getConstraintMultipliers(int isub);
    void subUpdateStates(int isub, DistrGeomState *refState, DistrGeomState *geomState, double time);
    void subInitializeMultipliers(int isub, DistrGeomState& geomState);
    void subInitializeParameters(int isub);
    void subUpdateMultipliers(int isub, DistrGeomState& geomState);
    void subUpdateParameters(int isub);
    void clean();
};

#endif
