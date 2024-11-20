#ifndef _CONTROL_INTERFACE_H_
#define _CONTROL_INTERFACE_H_

#include "Math.d/Vector.h"
template <class Vector> class SysState;
class SingleDomainDynamic;

class ControlInterface {

  protected:
    double dt;	       // time step size
    
  public:
  
    // set time step
    void setDt(double h) { dt = h; }

    // Control force initialization routine
    virtual void init(double *displacement, double *velocity, double *acc,
                      SingleDomainDynamic * probDesc=0) = 0;

    // Control force routine
    virtual void ctrl(double *dis, double *vel, double *acc, double *f,
                      double time=0, 
		      SysState<Vector> *state=0, Vector *ext_f=0) = 0;

    // User defined force routine
    virtual void usd_forc(double time, double *userDefineForc) = 0;

    // User defined displacement routine
    virtual void usd_disp(double time, double *userDefineDisp,
                          double *userDefineVel, double *userDefineAcc) = 0;
};

#endif
