#include <Driver.d/SysState.h>
#include <Driver.d/Domain.h>
#include <Corotational.d/utilities.h>
#include "PodProjectionNonLinDynamic.h"

extern GeoSource * geoSource;

namespace Rom {

SDDynamPodPostProcessor::SDDynamPodPostProcessor(Domain *d, double *bcx, double *vcx, double *acx,
                                                 StaticTimers* times, GeomState *geomState, Corotator **corot) :
    SDDynamPostProcessor(d, bcx, vcx, acx, times, geomState, corot),
    DispSensorValues(NULL),
    AccSensorValues(NULL),
    VelSensorValues(NULL)
{
  oinfo = geoSource->getOutputInfo();
  numOutInfo = geoSource->getNumOutInfo();

  DispSensor = false;
  AccSensor  = false;
  VelSensor  = false;

  for (int iOut = 0; iOut < numOutInfo; iOut++) {
    switch(oinfo[iOut].type) {
      case OutputInfo::Accel6 : case OutputInfo::Acceleration :
        if(oinfo[iOut].nodeNumber != -1) {
          AccSensor = true;
          if(oinfo[iOut].type == OutputInfo::Accel6 && (oinfo[iOut].angularouttype != OutputInfo::total || oinfo[iOut].rescaling)) {
            // in this case rotation vector and angular velocity are required for angular acceleration transformation
            DispSensor = true;
            VelSensor = true;
          }
        } else {
          oinfo[iOut].filptr = fopen(oinfo[iOut].filename, "wb");
          if(!oinfo[iOut].filptr) {
            fprintf(stderr," *** ERROR: Cannot open %s, exiting...\n", oinfo[iOut].filename);
            exit(-1);
          }
          filePrint(oinfo[iOut].filptr, "0\n"); 
        }
        break;
      case OutputInfo::Disp6DOF : case OutputInfo::Displacement : case OutputInfo::Temperature :
        if(oinfo[iOut].nodeNumber != -1) {
          DispSensor = true;
        } else {
          oinfo[iOut].filptr = fopen(oinfo[iOut].filename, "wb");
          if(!oinfo[iOut].filptr) {
            fprintf(stderr," *** ERROR: Cannot open %s, exiting...\n", oinfo[iOut].filename);
            exit(-1);
          }
          filePrint(oinfo[iOut].filptr, "1\n");
        }
        break;
      case OutputInfo::Velocity6 : case OutputInfo::Velocity : case OutputInfo::TemperatureFirstTimeDerivative :
        if(oinfo[iOut].nodeNumber != -1) {
          VelSensor = true;
          if(oinfo[iOut].type == OutputInfo::Velocity6 && (oinfo[iOut].angularouttype != OutputInfo::total || oinfo[iOut].rescaling)) {
            // in this case rotation vector is required for angular velocity transformation
            DispSensor = true;
          }
        } else {
          oinfo[iOut].filptr = fopen(oinfo[iOut].filename, "wb");
          if(!oinfo[iOut].filptr) {
            fprintf(stderr," *** ERROR: Cannot open %s, exiting...\n", oinfo[iOut].filename);
            exit(-1);
          }
          filePrint(oinfo[iOut].filptr, "2\n");
        }
        break; 
      default:
        filePrint(stderr, " *** WARNING: Online ROM output only supports GDISPLAC, GVELOCIT, and GACCELER.\n");
        filePrint(stderr, "     output type selected is %d \n", oinfo[iOut].type);
    }
  }

  if(DispSensor || VelSensor || AccSensor) {
    geoSource->openSensorOutputFiles();
  }

  buildSensorNodeVector();
}

SDDynamPodPostProcessor::~SDDynamPodPostProcessor() {

  if(DispSensorValues) delete DispSensorValues;
  if(AccSensorValues) delete AccSensorValues;
  if(VelSensorValues) delete VelSensorValues;
}

void
SDDynamPodPostProcessor::printPODSize(int PODsize) {

  podSize = PODsize;

  for (int iOut = 0; iOut < numOutInfo; iOut++) {
    switch(oinfo[iOut].type) {
      case OutputInfo::Accel6 : case OutputInfo::Acceleration :
        if(oinfo[iOut].nodeNumber == -1)
          filePrint(oinfo[iOut].filptr, "%d\n", PODsize);
        break;
      case OutputInfo::Disp6DOF : case OutputInfo::Displacement :
        if(oinfo[iOut].nodeNumber == -1)
          filePrint(oinfo[iOut].filptr, "%d\n", PODsize);
        break;
      case OutputInfo::Velocity6 : case OutputInfo::Velocity :
        if(oinfo[iOut].nodeNumber == -1)
          filePrint(oinfo[iOut].filptr, "%d\n", PODsize);
        break;
      default:
        break;
    }
  }
}

void
SDDynamPodPostProcessor::makeSensorBasis(VecBasis *fullBasis) {

  SensorBasis = fullBasis;
  SensorBasis->makeSparseBasis2(nodeVector, domain->getCDSA());

  // allocate space for containers to hold sensor values
  if(DispSensor) {
    DispSensorValues = new GenVector<double>(SensorBasis->vectorInfo());
    DispSensorValues->zero();
  }
  if(AccSensor) {
    AccSensorValues = new GenVector<double>(SensorBasis->vectorInfo());
    AccSensorValues->zero();
  }
  if(VelSensor) {
    VelSensorValues = new GenVector<double>(SensorBasis->vectorInfo());
    VelSensorValues->zero();
  }
}
  
void
SDDynamPodPostProcessor::buildSensorNodeVector() {

  //load vector of sensor nodes converted to local numbering
  for (int iOut = 0; iOut < numOutInfo; iOut++) {
    if(oinfo[iOut].nodeNumber != -1) {
      nodeVector.push_back(oinfo[iOut].nodeNumber);
    }
  }

  //if multiple outputs are requested for a single node, cull redundancies from node vector
  std::sort(nodeVector.begin(), nodeVector.end());
  std::vector<int>::iterator packedNodeIt = std::unique(nodeVector.begin(), nodeVector.end());
  nodeVector.resize(packedNodeIt-nodeVector.begin());
}

double
SDDynamPodPostProcessor::GetPrescribedSensorValue(int locNode, int j, double *bcdata)
{
  int dof = -1;
  switch(j) {
    case 0 : dof = domain->getDSA()->locate(locNode, DofSet::Xdisp); break;
    case 1 : dof = domain->getDSA()->locate(locNode, DofSet::Ydisp); break;
    case 2 : dof = domain->getDSA()->locate(locNode, DofSet::Zdisp); break;
    case 3 : dof = domain->getDSA()->locate(locNode, DofSet::Xrot); break;
    case 4 : dof = domain->getDSA()->locate(locNode, DofSet::Yrot); break;
    case 5 : dof = domain->getDSA()->locate(locNode, DofSet::Zrot); break;
  }

  return (dof != -1) ? bcdata[dof] : 0.0;
}

void
SDDynamPodPostProcessor::printSensorValues(GenVector<double> &SensorData, OutputInfo *OINFO, double *time, double *bcdata) {

  // XXX prescribed values in vcx and acx are convected.
  int locNode = OINFO->nodeNumber;
  int ndofs;
  int dofs[6];
  if(OINFO->type == OutputInfo::Disp6DOF || OINFO->type == OutputInfo::Velocity6 || OINFO->type == OutputInfo::Accel6) {
    ndofs = domain->getCDSA()->number(locNode, DofSet::XYZdisp | DofSet::XYZrot, dofs);
  }
  else if(OINFO->type == OutputInfo::Temperature || OINFO->type == OutputInfo::TemperatureFirstTimeDerivative) {
    ndofs = domain->getCDSA()->number(locNode, DofSet::Temp, dofs);
  }
  else {
    ndofs = domain->getCDSA()->number(locNode, DofSet::XYZdisp, dofs);
  }
  double *data = new double[ndofs];
  for(int j = 0; j < ndofs; ++j) {
    data[j] = (dofs[j] != -1) ? SensorData[dofs[j]] : GetPrescribedSensorValue(locNode,j,bcdata);
  }
  // transform rotation vector, if necessary
  if(OINFO->type == OutputInfo::Disp6DOF && (OINFO->rotvecouttype != OutputInfo::Euler || OINFO->rescaling)) {
    double rten[3][3], psi[3] = { data[3], data[4], data[5] };
    if(OINFO->rescaling) vec_to_mat(psi, rten);
    tran_rvec(rten, psi, OINFO->rescaling, OINFO->rotvecouttype, data+3);
  }
  // transform angular velocity, if necessary
  else if(OINFO->type == OutputInfo::Velocity6 && (OINFO->angularouttype != OutputInfo::total || OINFO->rescaling)) {
    double rten[3][3], psi[3], psidot[3] = { data[3], data[4], data[5] };
    for(int j = 0; j < 3; ++j) psi[j] = (dofs[3+j] != -1) ? (*DispSensorValues)[dofs[3+j]] : GetPrescribedSensorValue(locNode,3+j,bcx);
    if(OINFO->rescaling) vec_to_mat(psi, rten);
    tran_veloc(rten, psi, psidot, 2, OINFO->angularouttype, false, OINFO->rescaling, data+3);
  }
  // transform angular acceleration, if necessary
  else if(OINFO->type == OutputInfo::Accel6 && (OINFO->angularouttype != OutputInfo::total || OINFO->rescaling)) {
    double rten[3][3], psi[3], psidot[3], psiddot[3] = { data[3], data[4], data[5] };
    for(int j = 0; j < 3; ++j) {
      psi[j] = (dofs[3+j] != -1) ? (*DispSensorValues)[dofs[3+j]] : GetPrescribedSensorValue(locNode,3+j,bcx);
      psidot[j] = (dofs[3+j] != -1) ? (*VelSensorValues)[dofs[3+j]] : GetPrescribedSensorValue(locNode,3+j,vcx);
    }
    if(OINFO->rescaling) vec_to_mat(psi, rten);
    tran_accel(rten, psi, psidot, psiddot, 2, OINFO->angularouttype, false, OINFO->rescaling, data+3);
  }
  int w = OINFO->width;
  int p = OINFO->precision; 
  fprintf(OINFO->filptr, "  % *.*E  ", w, p, *time);
  for(int j = 0; j < ndofs; ++j) {
    fprintf(OINFO->filptr, " % *.*E", w, p, data[j]);
  }
  fprintf(OINFO->filptr, "\n");
  fflush(OINFO->filptr);
  delete [] data;
}

void
SDDynamPodPostProcessor::dynamOutput(int tIndex, double t, DynamMat &dynOps, Vector &externalForce,
                                     Vector *AeroF, SysState<Vector>& systemState) {

  //all MPI processes have a full copy of reduced coordinates, only master processes needs to print
  int p = std::numeric_limits<double>::digits10+1;

  bool DispProjected = false, AccProjected = false, VelProjected = false;

  for(int iOut = 0; iOut < numOutInfo; iOut++) {

    if(tIndex % oinfo[iOut].interval == 0) {

      switch(oinfo[iOut].type) {
         case OutputInfo::Accel6 : case OutputInfo::Acceleration :
           {
             if(oinfo[iOut].nodeNumber == -1) {
               filePrint(oinfo[iOut].filptr, "   %.*e\n", p, t); // print timestamp
               for(int i = 0; i < podSize; i++) {
                 filePrint(oinfo[iOut].filptr, "%.*e ", p, systemState.getAccel()[i]);
               }
               filePrint(oinfo[iOut].filptr, "\n");
             }
             else {
               if(!AccProjected) {
                 SensorBasis->expand2(systemState.getAccel(), *AccSensorValues);
                 AccProjected = true;
               }
               if(oinfo[iOut].type == OutputInfo::Accel6 && (oinfo[iOut].angularouttype != OutputInfo::total || oinfo[iOut].rescaling)) {
                 // in this case rotation vector and angular velocity are required for angular acceleration transformation
                 if(!DispProjected) {
                   SensorBasis->expand2(systemState.getDisp(), *DispSensorValues);
                   DispProjected = true;
                 }
                 if(!VelProjected) {
                   SensorBasis->expand2(systemState.getVeloc(), *VelSensorValues);
                   VelProjected = true;
                 }
               }
               printSensorValues(*AccSensorValues, &oinfo[iOut], &t, acx);
             }
           }
           break;
         case OutputInfo::Disp6DOF : case OutputInfo::Displacement : case OutputInfo::Temperature :
           {
             if(oinfo[iOut].nodeNumber == -1) {
               filePrint(oinfo[iOut].filptr, "   %.*e\n", p, t); // print timestamp
               for(int i = 0; i < podSize; i++) {
                 filePrint(oinfo[iOut].filptr, "%.*e ", p, systemState.getDisp()[i]);
               }
               filePrint(oinfo[iOut].filptr, "\n");
             }
             else {
               if(!DispProjected) {
                 SensorBasis->expand2(systemState.getDisp(), *DispSensorValues);
                 DispProjected = true;
               }
               printSensorValues(*DispSensorValues, &oinfo[iOut], &t, bcx);
             }
           }
           break;
         case OutputInfo::Velocity6 : case OutputInfo::Velocity : case OutputInfo::TemperatureFirstTimeDerivative :
           {
             if(oinfo[iOut].nodeNumber == -1) {
               filePrint(oinfo[iOut].filptr, "   %.*e\n", p, t); // print timestamp
               for(int i = 0; i < podSize; i++) {
                 filePrint(oinfo[iOut].filptr, "%.*e ", p, systemState.getVeloc()[i]);
               }
               filePrint(oinfo[iOut].filptr, "\n");
             }
             else {
               if(!VelProjected) {
                 SensorBasis->expand2(systemState.getVeloc(), *VelSensorValues);
                 VelProjected = true;
               }
               if(oinfo[iOut].type == OutputInfo::Velocity6 && (oinfo[iOut].angularouttype != OutputInfo::total || oinfo[iOut].rescaling)) {
                 // in this case rotation vector is required for angular velocity transformation
                 if(!DispProjected) {
                   SensorBasis->expand2(systemState.getDisp(), *DispSensorValues);
                   DispProjected = true;
                 }
               }
               printSensorValues(*VelSensorValues, &oinfo[iOut], &t, vcx);
             }
           }
           break;
         default:
           filePrint(stderr, " ... ROM output only supports Acceleration, Displacement, and Velocity ... \n");
      }
    }
  }
}

} /* end namespace Rom */
