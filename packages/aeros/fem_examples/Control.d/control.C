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

const double F0 = 1.0;
const double omega = 1.0;

void
MyControl::ctrl(double *disp, double *vel, double *accel, double *force,
                double time, SysState<Vector> *state, Vector *ext_f)
{
 // disp  = sensor displacement
 // vel   = sensor velocity
 // accel = sensor acceleration
 std::cerr << "Applying Control Force at time = " << time << std::endl;
 force[0] = F0*std::cos(omega*disp[0]);
 force[1] = F0*std::cos(omega*disp[1]);
 force[2] = F0*std::cos(omega*disp[2]);
}

void
MyControl::init(double *disp, double *vel, double *accel, SingleDomainDynamic *probDesc)
{
}

void
MyControl::usd_disp(double time, double *userDefineDisp, double *userDefineVel,
                    double *userDefineAccel)
{
}

void
MyControl::usd_forc(double time, double *usdForce)
{
}
