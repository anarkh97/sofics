#include <cstdio>
#include <cmath>

#include "ControlInterface.h"

class MyControl : public ControlInterface {

  public:
  
    // initialization routine
    void init(double *displacement, double *velocity, double *acceleration,
              SingleDomainDynamic * probDesc=0);

    // control force
    void ctrl(double *displacement, double *velocity, double *acceleration, 
              double *force, double time=0, SysState<Vector> *state=0, 
	      Vector *ext_f=0);
	      
    void usd_disp(double time, double *userDefineDisplacement, 
                  double *userDefineVelocity, double *userDefineAcceleration);

    void usd_forc(double time, double *userDefineForce);
    
};


// Define new control Object (i.e. control force function, user defined
// force function, user defined displacement function)

ControlInterface *controlObj = new MyControl();

// disp  = displacement
// vel   = velocity
// accel = acceleration

void
MyControl::ctrl(double *displacement, double *velocity, double *acceleration, 
              double *force, double time, SysState<Vector> *state, 
	      Vector *ext_f)
{

}


// init = initialization routine

void
MyControl::init(double *displacement, double *velocity, double *acceleration,
                SingleDomainDynamic * probDesc)
{

}

// usd_disp = user defined displacements
//
// This function allows the user to define a time-dependent displacement
// function for selected degrees of freedom in a structural model.
// The user selects the degrees of freedom from the input file using the
// USDD input command and then works with a local numbering system inside
// of the usd_disp() function. So if the user selects 3 degrees of freedom,
// they will be numbered sequentially beginning with 0 in the same order
// they were selected in the input file.

void
MyControl::usd_disp(double time, double *userDefineDisp, double *userDefineVel,
                    double *userDefineAcc)
{
 // blank intentionally
}

// usd_forc = user defined forces
//
// This function allows the user to define a time-dependent forcing 
// function for selected degrees of freedom in a structural model.
// The user selects the degrees of freedom from the input file with the
// USDF input command and works with a local numbering system as defined
// above in the description of the usd_disp() command.

void
MyControl::usd_forc(double time, double *usdForce)
{

 double F0    = 1.0;
 double omega = 2.0;

 //usdForce[0] = F0*cos(omega*time);
 //usdForce[1] = F0*cos(omega*time);
 //usdForce[2] = F0*cos(omega*time);
 //usdForce[3] = F0*cos(omega*time);

 usdForce[0] = F0*time;
 usdForce[1] = F0*time;
 usdForce[2] = F0*time;
 usdForce[3] = F0*time;

 fprintf(stderr,"Time %e TimeStep %e Applied Force %e\n",time,dt,usdForce[0]);

}
