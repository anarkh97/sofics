#include <iostream>
#include <cmath>
#include <Control.d/ControlInterface.h>

class MyControl : public ControlInterface {

  public:
 
    // initialization routine
    void init(double *displacement, double *velocity, double *acceleration,
              SingleDomainDynamic * probDesc=0);

    // actuator force
    void ctrl(double *displacement, double *velocity, double *acceleration, 
              double *force, double time=0, SysState<Vector> *state=0, 
	      Vector *ext_f=0);
	      
    void usd_disp(double time, double *userDefineDisplacement, 
                  double *userDefineVelocity);

    void usd_forc(double time, double *userDefineForce);

};

ControlInterface *controlObj = new MyControl();

void
MyControl::init(double *displacement, double *velocity, double *acceleration,
                SingleDomainDynamic * probDesc)
{
 // blank intentionally
}

void
MyControl::ctrl(double *displacement, double *velocity, double *acceleration, 
                double *force, double time, SysState<Vector> *state, 
                Vector *ext_f)
{
  // implement this function for ACTUATOR controls, i.e. forces which are a function of the state at the previous converged timestep
}

void
MyControl::usd_disp(double time, double *userDefineDisp, double *userDefineVel)
{
  // implement this function for USDF controls, i.e. prescribed displacements which are a function of time
}

void
MyControl::usd_forc(double time, double *usdForce)
{
 // implement this function for USDF constrols, i.e. forces which are a function of time
 // don't forget to include: LOAD "./control.so" in the aeros input file
 // type make to compile control.C into control.so
 double w = 10;
 double f = 1.0;
 usdForce[0] = f*sin(w*time);
}

