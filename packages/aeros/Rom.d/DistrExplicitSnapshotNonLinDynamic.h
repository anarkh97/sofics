#ifndef ROM_DISTREXPLICITSNAPSHOTNONLINDYNAMIC_H
#define ROM_DISTREXPLICITSNAPSHOTNONLINDYNAMIC_H

#include <Paral.d/MDDynam.h>

#include <Problems.d/DynamProbTraits.h>

namespace Rom {

class DistrExplicitSnapshotNonLinDynamic : public MultiDomainDynam {
public:
  explicit DistrExplicitSnapshotNonLinDynamic(Domain *);

  // Overriding via hiding
  void preProcess(); // Required additional pre-processing
  // Intercept passing information (2 overloads)
  void computeExtForce2(SysState<DistrVector> &distState,
                        DistrVector &f, DistrVector &cnst_f,
                        int tIndex, double t,
                        DistrVector *aero_f,
                        double gamma, double alphaf);
  void computeExtForce2(SysState<DistrVector> &distState,
                        DistrVector &f, DistrVector &cnst_f,
                        int tIndex, double t,
                        DistrVector *aero_f);

  // Added functionality
  void currentTimeIs(double t);
  void stateSnapshotAdd(const DistrVector &state);
  void accelerationSnapshotAdd(const DistrVector &accel);
  void velocitySnapshotAdd(const DistrVector &veloc);
  void forceSnapshotAdd(const DistrVector &f);

  ~DistrExplicitSnapshotNonLinDynamic();

  bool collectState;
  bool collectAccel;
  bool collectVeloc;
  bool collectForce;

private:

  Domain *domain_;

  class SnapshotHandler;
  friend class SnapshotHandler;
  SnapshotHandler * snapshotHandler_;

  // Disallow copy and assignment
  DistrExplicitSnapshotNonLinDynamic(const DistrExplicitSnapshotNonLinDynamic &);
  DistrExplicitSnapshotNonLinDynamic &operator=(const DistrExplicitSnapshotNonLinDynamic &);
};

inline
void
DistrExplicitSnapshotNonLinDynamic::computeExtForce2(SysState<DistrVector> &distState,
                                                     DistrVector &f, DistrVector &cnst_f,
                                                     int tIndex, double t,
                                                     DistrVector *aero_f,
                                                     double gamma, double alphaf) {
  currentTimeIs(t);
  MultiDomainDynam::computeExtForce2(distState, f, cnst_f, tIndex, t, aero_f, gamma, alphaf);
}

inline
void
DistrExplicitSnapshotNonLinDynamic::computeExtForce2(SysState<DistrVector> &distState,
                                                     DistrVector &f, DistrVector &cnst_f,
                                                     int tIndex, double t,
                                                     DistrVector *aero_f) {
  currentTimeIs(t);
  MultiDomainDynam::computeExtForce2(distState, f, cnst_f, tIndex, t, aero_f);
}

} // end namespace Rom

inline
void
handleDisplacement(Rom::DistrExplicitSnapshotNonLinDynamic &probDesc, DistrVector &d) {
  if(probDesc.collectState)
    probDesc.stateSnapshotAdd(d);
}

inline
void
handleAcceleration(Rom::DistrExplicitSnapshotNonLinDynamic &probDesc, DistrVector &d) {
  if(probDesc.collectAccel)
    probDesc.accelerationSnapshotAdd(d);
}

inline
void
handleVelocity(Rom::DistrExplicitSnapshotNonLinDynamic &probDesc, DistrVector &d) {
  if(probDesc.collectVeloc)
    probDesc.velocitySnapshotAdd(d);
}

inline
void
handleForce(Rom::DistrExplicitSnapshotNonLinDynamic &probDesc, DistrVector &d) {
  if(probDesc.collectForce)
    probDesc.forceSnapshotAdd(d);
}

#endif /* ROM_DISTREXPLICITSNAPSHOTNONLINDYNAMIC_H */
