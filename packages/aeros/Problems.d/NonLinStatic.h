#ifndef _NON_LIN_STATIC_H_
#define _NON_LIN_STATIC_H_

class Domain;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
template <class Scalar> class GenSolver;
typedef GenSolver<double> Solver;
template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;
class StaticTimers;
class GeomState;
class Corotator;
template<class T, class VectorType, class SolverType> class SingleDomainPostProcessor;

class NonLinStatic {
    Domain *domain;
    double *bcx;
    Solver *solver;
    SparseMatrix *spm;
    Solver *prec;
    SparseMatrix *spp;
    FullSquareMatrix *kelArray;
    Corotator **allCorot;
    double firstRes;
    double firstDv;
    double tolerance;
    StaticTimers *times;
    Vector *reactions;
    bool updateCS;

 public:
    // Constructor
    NonLinStatic(Domain *d);
    ~NonLinStatic();

    int  solVecInfo();
    int  sysVecInfo();
    int  elemVecInfo();
    int  getMaxit();
    double getScaleFactor();  // only nlstatic
    double getDeltaLambda0(); // only nlstatic
    double getMaxLambda();    // only maxlambda
    virtual void getRHS(Vector &rhs); 
    void preProcess(bool factor = true);
    Solver *getSolver();
    SingleDomainPostProcessor<double,Vector,Solver> *getPostProcessor();

    int reBuild(int iter, int step, GeomState& geomState);
    virtual GeomState* createGeomState();

    virtual void staticOutput(GeomState *geomState, double lambda, Vector& force, Vector &, GeomState *refState);
    int checkConvergence(int iter, double normDv, double residualNorm);

    void updateContactSurfaces(GeomState& geomState);
    void updateStates(GeomState *refState, GeomState& geomState, double lambda);

    double getStiffAndForce(GeomState& geomState, Vector& residual, 
                            Vector& elementInternalForce, Vector &,
                            double lambda = 1, GeomState *refState = NULL, bool forceOnly = false);

    void updatePrescribedDisplacement(GeomState *geomState, double lambda = 1);
    void initializeParameters(int step, GeomState *geomState);
    void updateParameters(GeomState *geomState);
    bool checkConstraintViolation(double &err, GeomState *geomState);

    void printTimers();

    double getTolerance(); 

    LinesearchInfo& linesearch(); 
    double getEnergy(double lambda, Vector& force, GeomState* geomState);

    double getResidualNorm(Vector &res, GeomState &geomState);
    bool getResizeFlag();

  private:
    void clean();
};

#endif
