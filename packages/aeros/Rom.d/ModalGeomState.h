#ifndef _MODAL_GEOM_STATE_H
#define _MODAL_GEOM_STATE_H

#include <Math.d/Vector.h>
#include <Feti.d/DistrVector.h>

namespace Rom {
class PodProjectionNonLinDynamic;
class LumpedPodProjectionNonLinDynamic;

class DistrPodProjectionNonLinDynamic;
class DistrLumpedPodProjectionNonLinDynamic;
}

template<class VecType>
class GenModalGeomState {
/* class used by PodProjectionNonLinDynamic for storing the geometric state of the structure
   passed as a template parameter GeomType to class NLDynamSolver defined in
   Driver.d/NLDynamProbType.[Ch]
*/
protected:

  VecType q;          // modal deformation coefficients
  VecType vel;
  VecType acc;

  int numFlex;

public:

  GenModalGeomState(int numModes);
  GenModalGeomState(DistrInfo &dinfo);
  GenModalGeomState(const GenModalGeomState& mgs);

  void update(const VecType &, int SO3param = 0);
  void get_inc_displacement(VecType &incDsp, GenModalGeomState &stepState, bool)
    {
      incDsp = q - stepState.q;
    }
  void get_tot_displacement(VecType &totVec, bool rescaled = true) {}
  void midpoint_step_update(VecType &veloc_n, VecType &accel_n, double delta, GenModalGeomState &ss,
                            double beta, double gamma, double alphaf, double alpham,
                            bool zeroRot);

  void printState(const char* = "");
  void setVelocity(const VecType &_vel, int SO3param = 0) { vel = _vel; }
  void setAcceleration(const VecType &_acc, int SO3param = 0) { acc = _acc; } 
  void setVelocityAndAcceleration(VecType &_vel, VecType &_acc) { vel = _vel; acc = _acc; }
  void push_forward(VecType &) {}
  void pull_back(VecType &) {}
  void print() {}

  VecType *getModalq() {return &q; }

  friend class Rom::PodProjectionNonLinDynamic;
  friend class Rom::LumpedPodProjectionNonLinDynamic;
  
  friend class Rom::DistrPodProjectionNonLinDynamic;
  friend class Rom::DistrLumpedPodProjectionNonLinDynamic;
};

typedef GenModalGeomState<Vector> ModalGeomState;
typedef GenModalGeomState<DistrVector> DistrModalGeomState;

#endif
