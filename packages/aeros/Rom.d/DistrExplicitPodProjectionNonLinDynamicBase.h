#ifndef ROM_DISTREXPLICITPODPROJECTIONNONLINDYNAMICBASE_H
#define ROM_DISTREXPLICITPODPROJECTIONNONLINDYNAMICBASE_H

#include <Paral.d/MDDynam.h>
#include <Driver.d/GeoSource.h>
#include <cstdio>
#include "DistrVecBasis.h"

namespace Rom {

class MultiDomDynPodPostProcessor : public MultiDomDynPostProcessor {

public:
 
  MultiDomDynPodPostProcessor(DecDomain *, StaticTimers*, DistrGeomState *, Corotator ***);
  ~MultiDomDynPodPostProcessor();

  void dynamOutput(int, double, MDDynamMat &, DistrVector &, DistrVector *aeroF, SysState<DistrVector> &);
  void printPODSize(int);
  void makeSensorBasis(DistrVecBasis *);
  
private:

  OutputInfo *oinfo;
  int numOutInfo;
  int podSize;

  void subBuildSensorNodeVector(int);
  void subPrintSensorValues(int, GenDistrVector<double> &, OutputInfo *, double *, int);
  double subGetPrescribedSensorValue(int, int, int, int);

  DistrVecBasis *SensorBasis;
  GenDistrVector<double> *DispSensorValues;
  GenDistrVector<double> *AccSensorValues;
  GenDistrVector<double> *VelSensorValues;

  DofSetArray **all_cdsa;

  bool DispSensor;
  bool AccSensor;
  bool VelSensor;

  std::vector<std::vector<int> > nodeVector;
};

class DistrExplicitPodProjectionNonLinDynamicBase : public MultiDomainDynam {
public:
  explicit DistrExplicitPodProjectionNonLinDynamicBase(Domain *);

  // Overriding via hiding
  void preProcess(); // Additional pre-processing
  const DistrInfo &solVecInfo();
  DistrInfo &reducedVecInfo();
  void getInitState(SysState<DistrVector> &);
  void updateState(double, DistrVector&, DistrVector&);
  void getContactForce(DistrVector &d_n, DistrVector &dinc,
                       DistrVector &ctc_f, double t_n_p,
                       double dt, double dt_old);
  void computeExtForce2(SysState<DistrVector> &distState,
                        DistrVector &f, DistrVector &cnst_f, int tIndex,
                        double t, DistrVector *aero_f=0,
                        double gamma=0.5, double alphaf=0.0);
  void computeStabilityTimeStep(double&, MDDynamMat&);
  void getConstForce(DistrVector& v);
  void getInternalForce(DistrVector &d, DistrVector &f, double t, int tIndex);
  void printFullNorm(DistrVector &);
  MultiDomDynPodPostProcessor *getPostProcessor();
  MDDynamMat * buildOps(double, double, double); // Reduced-order matrix solver

  //local basis functions
  void initLocalBasis(DistrVector &q);
  int  selectLocalBasis(DistrVector &q);
  void setLocalBasis(DistrVector &q, DistrVector &qd);
  void projectLocalBases(int i, int j, DistrVector &q);
  virtual void setLocalReducedMesh(int j) {}

protected:
  Domain        *domain;
  DistrInfo      reducedInfo;
  DistrVecBasis  normalizedBasis_;
  DistrVecBasis  centroids;
  //dummy Variables to fascillitate computation on Reduced Coordinates
  StaticTimers * times;
  DistrVector  * fExt;
  DistrVector  * fInt;
  DistrVector  * cnst_fBig;
  DistrVector  * aero_fBig;
  DistrVector  * d_n;
  DistrVector  * v_n;
  DistrVector  * a_n;
  DistrVector  * v_p;
  DistrVector  * tempVec;
  DistrVector  * ctc_f;
  SysState<DistrVector> *dummyState;
  bool haveRot;
  int stableCount;
  GenParallelSolver<double> * fullMassSolver;
  void reduceDisp(DistrVector &d, DistrVector &dr) const;

#ifdef USE_EIGEN3
  Eigen::Array<Eigen::MatrixXd,Eigen::Dynamic,Eigen::Dynamic> VtV;
  Eigen::Array<double,Eigen::Dynamic,Eigen::Dynamic> d;
  Eigen::Array<Eigen::VectorXd,Eigen::Dynamic,Eigen::Dynamic> w;
#endif

private:
  MultiDomDynPodPostProcessor *mddPostPro;
  // Disallow copy and assignment
  DistrExplicitPodProjectionNonLinDynamicBase(const DistrExplicitPodProjectionNonLinDynamicBase &);
  DistrExplicitPodProjectionNonLinDynamicBase &operator=(const DistrExplicitPodProjectionNonLinDynamicBase &);
  void subInitStiff(int isub);

  int localBasisId; 
  std::vector<int> locBasisVec;
};

} // end namespace Rom

#endif /* ROM_DISTREXPLICITPODPROJECTIONNONLINDYNAMICBASE_H */
