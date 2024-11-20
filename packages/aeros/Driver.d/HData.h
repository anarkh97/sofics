#ifndef _HDATA_H_
#define _HDATA_H_

#include <memory>
#include <Element.d/Sommerfeld.d/SommerElement.h>
#include <Math.d/ComplexD.h>
#include <list>
#include <iostream>

template <class Scalar> class GenCuCSparse;
typedef GenCuCSparse<DComplex> CuCComplexSparse;
template <class Scalar> class GenSolver;
typedef GenSolver<DComplex> ComplexSolver;
template <class Scalar> class GenDBSparseMatrix;
typedef GenDBSparseMatrix<double> DBSparseMatrix;
typedef GenDBSparseMatrix<DComplex> DBComplexSparseMatrix;
template <class Scalar> class GenSkyMatrix;
typedef GenSkyMatrix<double> SkyMatrix;
template <class Scalar> class GenBLKSparseMatrix;
typedef GenBLKSparseMatrix<double> BLKSparseMatrix;
template <class Scalar> class GenSparseMatrix;
typedef GenSparseMatrix<DComplex> ComplexSparseMatrix;
typedef GenSparseMatrix<double> SparseMatrix;
template <class Scalar> class GenVector;
typedef GenVector<DComplex> ComplexVector;
typedef GenVector<double> Vector;
template <class T> class ResizeArray;

class SommerElement;
class Domain;
class Connectivity;
class CoordSet;
struct ComplexBCond;


// This class is designed to handle data particular to Helmholtz

class HData
{
public:
    static double coupledScaling, cscale_factor, cscale_factor2;
    int nffp;
    int limitffp;
    double *waveDirections;
    int numWaveDirections;
    int numFFPDirections;
    double *ffpDirections;
    int numKirchhoffLocations;
    double *kirchhoffLocations;
    int numFrequencies;
    int iWaveDir;
    int iRHS;
    int implicitFlag;
    int PMLFlag;
    int pointSourceFlag;
    // Bayliss-Turkel condition
    // 0 == "zero order"; 1 == 1st ord. Bayliss-Turkel; 2 == 2nd. ord. Bayliss-Turkel
    int sommerfeldType;
    int curvatureFlag; // Flag to check for explicit curvatures
    double curvatureConst1; // value of the explicit curvatures
    double curvatureConst2; // value of the explicit curvatures
    ComplexSolver *csolver;    // this domains complex solver
    int numComplexDirichlet;   // number of complex dirichlet bc
    ComplexBCond* cdbc;        // set of those bc
    int numComplexNeuman;      // number of complex Neuman bc
    ComplexBCond* cnbc;        // set of complex Neuman bc
    int numSommer;             // number of Sommerfeld bc
    int numScatter;          // number of element on scatterer
    int numNeum;
    int numWet;
    ResizeArray<SommerElement *> sommer; // set of Sommerfeld bc's
    ResizeArray<SommerElement *> scatter; // set of scatterer elements for ffp
    ResizeArray<SommerElement *> neum; // set of Neumann bc's
    ResizeArray<SommerElement *> wet; // set of matching wet interface eles
    int numSBoundNodes;
    ResizeArray<int> sBoundNodes;
    double (*normalScatter)[3];
    double *lenScatter;
    ComplexD *dudnScatter;
    CuCComplexSparse *KucC;
    DBComplexSparseMatrix *KccC;
    std::unique_ptr<DBSparseMatrix> Kss;
    BLKSparseMatrix *Mff;
    //SkyMatrix *Mff;
    int *sBoundMap;
    int numComplexLMPC;        // number of complex Linear Multi-Point Constraints
    bool sommerChecked;

    HData();
    virtual ~HData();

    void make_bc(Domain *dom, int *bc, ComplexD *bcx);
    void makeKss(Domain*);
    void ffp(Domain*, int ndir, ComplexD *ffp, double (*dir)[3], ComplexD *u, bool direction);
    template<class Scalar>
    void wError(Domain*, double *l2err, double *h1err, double *l2, double *h1, Scalar *u);
    void addSommerElem(int en, int type, double cc, int nn, int *nodes);
    void addWetElem(int en, int type, double cc, int nn, int *nodes);
    void addScatterElem(int en, int type, double cc, int nn, int *nodes);
    void addNeumElem(int en, int type, double cc, int nn, int *nodes, PressureBCond* = 0);
    void addhScatter(int n, double *v);
    void inithScatter(int m, int n);
    int  setComplexDirichlet(int,ComplexBCond *);
    int  setComplexNeuman(int, ComplexBCond * );
    void addSommer(SommerElement *ele);
    void addScatter(SommerElement *ele);
    void addNeum(SommerElement *ele);
    void addWet(SommerElement *ele);
    void addSBoundNodes();
    void addSBoundNode(int);
    void checkSommerTypeBC(Domain *domain, const Connectivity *_elemToNode = 0, const Connectivity *_nodeToElem = 0);
    int  nCDirichlet() const { return numComplexDirichlet; }
    int  numSommerfeld() const { return numSommer; }
    int  numSSN() const { return numSommer + numScatter + numNeum + numWet; }
    int  getNumWet() { return numWet; }
    void makeSommerConnectivities();
    // returns the curvature of the boundary for the Bayliss-Turkel bc
    void getCurvatures(Domain *dom);
    void getCurvatures3D(Domain *dom);
    void getCurvatures3Daccurate(Domain *dom);
    std::unique_ptr<Connectivity> somElemToNode;
    std::unique_ptr<Connectivity> somNodeToElem;
    std::unique_ptr<Connectivity> somNodeToNode;
    void getTau(double [3], double [3], double [3]);

    double *curvatures;
    double *curvaturesH;
    double *curvaturesLapH;
    double *curvaturesK;
    double *curvatures_e;
    double *curvatures_f;
    double *curvatures_g;
    double (*curvatures_normal)[3];
    double (*curvatures_tau1)[3];
    double (*curvatures_tau2)[3];
    int    *nodeToSommerNodeMap;

//     void   setWaveNumber(double _kappa) { kappa = _kappa; }
//     double getWaveNumber() { return kappa; }
    void   setFFP(int);
    void   setFFP(int,int);
    void   setFFPDirections(double d1, double d2, double d3);
    void   setKirchhoffLocations(double x, double y, double z);
    void   addFrequencies1(double w0, double deltaw, int n_deltaw, bool minus_flag = true);
    void   addFrequencies2(double w0, double w1, int n_deltaw, bool first_flag = true);
    void   addFrequencies(double w_start, double w_end, int n_coarse, int n_fine);
    void   addFrequencies(double w_end, int n_fine);
    void   addFrequency(double w);
    void   addCoarseFrequency(double w);
    void   initFreqSweep(double w0);

    std::map<int,std::list<double> *> set_of_frequencies;
    std::map<int,std::list<double> *> set_of_coarse_frequencies;
    void setFrequencySet(int i) {
        if (set_of_frequencies.find(i)==set_of_frequencies.end())
            set_of_frequencies[i] = new std::list<double>();
        if (set_of_coarse_frequencies.find(i)==set_of_coarse_frequencies.end())
            set_of_coarse_frequencies[i] = new std::list<double>();
        frequencies = set_of_frequencies[i];
        coarse_frequencies = set_of_coarse_frequencies[i];
        initFreqSweep(0.0);
    }
    std::list<double> *frequencies;
    std::list<double> *coarse_frequencies;
    bool isCoarseGridSolve;

    // Shape Opt
    int ihScatter;
    int numhScatter;
    int lenhScatter;
    double (*hScatter)[3];
    Connectivity *scaElemToNode;
    Connectivity *scaNodeToElem;
    int *scaToEl;
    int *subScaToSca;

    // Implicit Dirichlet boundary conditions routine
    void setWaveDirections(int , double d1, double d2);
    void setWaveDirections(int, double d1, double d2, double d3);

    // Functions used by the implicit Dirichlet/Neumann boundary condition
    double *getWaveDirection() { return waveDirections + 3*iWaveDir;}
    int getNumWaveDirections() const { return numWaveDirections;}
    int getImplicitFlag() const { return implicitFlag; }

    // virtual void updateMatrices(SparseMatrix*,int*,FullSquareMatrix*,FullSquareMatrix*);
    virtual void addFSRHS(ComplexVector &force);
    virtual void getEigen(SparseMatrix*) { ; }

    // output
    void outputFFP(ComplexVector& sol, int iInfo);
    void outputFFP(Vector& sol, int iInfo)
    { std::cerr << " *** WARNING: HData::outputFFP(Vector& sol, int iInfo) is not implemented \n"; }
};

#endif
