#include <Driver.d/ControlLawInfo.h>
#include <Element.d/Element.h>
#include <cstdio>

ControlLawInfo::ControlLawInfo()
{
  fileName = routineName = 0;
  numSensor = numActuator = numUserDisp = numUserForce = 0;
  sensor = 0; actuator = 0; userDisp = 0; userForce = 0;
}

ControlLawInfo::~ControlLawInfo()
{
  if(sensor) delete [] sensor;
/*
  if(actuator) delete [] actuator;
  if(userDisp) delete [] userDisp;
  if(userForce) delete [] userForce;
*/
}

void ControlLawInfo::print()
{
  fprintf(stderr, " Number of Sensors: %d\n", numSensor);
  fprintf(stderr, " Number of Actuators: %d\n", numActuator);
  fprintf(stderr, " Number of UserForces: %d\n", numUserForce);
  fprintf(stderr, " Number of UserDisps: %d\n", numUserDisp);
  fprintf(stderr, " filename: %s\n", fileName);
  fprintf(stderr, " routine: %s\n", routineName);

  for (int j = 0; j < numUserDisp; j++)
    fprintf(stderr, " usdd[%d][%d] = %f\n", userDisp[j].nnum, userDisp[j].dofnum, userDisp[j].val);
}

void ControlLawInfo::makeGlobalClaw(ControlLawInfo *subClaw)
{
  //only the numbers are augmented for sower input
  numSensor += subClaw->numSensor;
  numActuator += subClaw->numActuator;
  numUserForce += subClaw->numUserForce;
  numUserDisp += subClaw->numUserDisp;
}
