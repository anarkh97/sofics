#ifndef _CONTROLLAWINFO_H_
#define _CONTROLLAWINFO_H_

class BCond;

// Control Law Information class, see Control.d/control.C
// to implement user defined control forces, user defined forces or
// user defined displacements

struct ControlLawInfo
{
  char *fileName;
  char *routineName;
  int numSensor;        // number of sensors
  BCond *sensor;
  int numActuator;      // number of actuators
  BCond *actuator;
  int numUserDisp;      // number of user defined displacements
  BCond *userDisp;
  int numUserForce;     // number of user defined forces
  BCond *userForce;

  ControlLawInfo();
  ~ControlLawInfo();
  void print();
  void makeGlobalClaw(ControlLawInfo *subClaw);
};

#endif
