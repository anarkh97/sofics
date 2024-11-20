#include <iostream>
#include <cmath>

#include "Control.d/ControlInterface.h"

class MyControl : public ControlInterface {

  public:
  
    void init(double *displacement, double *velocity, double *acceleration,
              SingleDomainDynamic *probDesc=0);

    void ctrl(double *displacement, double *velocity, double *acceleration, 
              double *force, double time=0, SysState<Vector> *state=0, Vector *ext_f=0);
	      
    void usd_disp(double time, double *userDefineDisp, double *userDefineVel,
                  double *userDefineAccel);

    void usd_forc(double time, double *userDefineForce);
    
};


ControlInterface *controlObj = new MyControl();


void
MyControl::ctrl(double *disp, double *vel, double *accel, double *force,
                double time, SysState<Vector> *state, Vector *ext_f)
{
}

void
MyControl::init(double *disp, double *vel, double *accel, SingleDomainDynamic *probDesc)
{
}

void
MyControl::usd_disp(double time, double *userDefineDisp, double *userDefineVel,
                    double *userDefineAccel)
{
 std::cout << "Applying time-dependent prescribed displacement at time = " << time << std::endl;

 userDefineDisp[0] = 4.0*sin(time);
 userDefineDisp[1] = 4.0*sin(time);
 userDefineDisp[2] = 4.0*sin(time);
 userDefineVel[0] = 0.0;
 userDefineVel[1] = 0.0;
 userDefineVel[2] = 0.0;
 userDefineAccel[0] = 0.0;
 userDefineAccel[1] = 0.0;
 userDefineAccel[2] = 0.0;
}

void
MyControl::usd_forc(double time, double *usdForce)
{
}
