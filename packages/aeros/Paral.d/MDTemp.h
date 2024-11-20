#ifndef _MD_TEMP_H_
#define _MD_TEMP_H_

#include <Paral.d/MDDynam.h>
#include <Driver.d/TempProbType.h>
#include <Driver.d/Dynam.h>

class MDTempDynamPostProcessor {
    DecDomain *decDomain;
    DistrGeomState *geomState;
    Corotator ***allCorot;
  public:
    MDTempDynamPostProcessor(DecDomain *d, DistrGeomState *g, Corotator ***c)
     : decDomain(d), geomState(g), allCorot(c) {}
    void tempdynamOutput(int, MDDynamMat&, DistrVector&, TempState<DistrVector>&);
};

class MultiDomainTemp {
    Domain *domain;
    DecDomain *decDomain;
    FullSquareMatrix **kelArray;
    Corotator ***allCorot;
    DistrGeomState *geomState;
    MDDynamMat *dynMat;

  public:
    MultiDomainTemp(Domain *d);
    ~MultiDomainTemp();

    MDDynamMat buildOps(double, double, double);

    MDTempDynamPostProcessor *getPostProcessor();
    DistrInfo & solVecInfo();
    int getTimeIntegration();
    int getAeroheatFlag();
    int getThermohFlag();
    int getHzemFlag();
    int getZEMFlag();

    void temptrProject(DistrVector &f);
    void tempProject(DistrVector &v);
    void getTempTimes(double &dtemp, double &tmax);
    void getTempNeum(double &epsiln);
    void tempInitState(TempState<DistrVector> &);
    void computeExtForce(DistrVector &, double t, int tIndex, DistrVector &);
    void preProcess();
    void aeroHeatPreProcess(DistrVector& d_n, DistrVector& v_n, DistrVector& v_p);
    void thermohPreProcess(DistrVector& d_n, DistrVector& v_n, DistrVector& v_p);
    void getSteadyStateParam(int &steadyFlag, int &steadyMin, int &steadMax,
                             double &steadyTol);
    void getQuasiStaticParameters(double &maxVel, double &qsbeta);
    void getInternalForce(DistrVector&, DistrVector&);
    void getInitialTime(int &tIndex, double &initialTime);

    int getModeDecompFlag();
    void modeDecomp(double t, int tIndex, DistrVector& d_n);

    int cmdComHeat(int cmdFlag);

  private:
    void subTempInitState(int isub, TempState<DistrVector> &inState);
    void subComputeExtForce(int isub, DistrVector &ext_f, double t, int tIndex, DistrVector &prev_f);
    void makeAllDOFs(int isub);
    void initSubPrescribedTemperature(int isub);
    void makeSubCorotators(int isub);
    void makeSubElementArrays(int isub);
    void subGetInternalForce(int isub, DistrVector& d, DistrVector& f);
};

#endif
