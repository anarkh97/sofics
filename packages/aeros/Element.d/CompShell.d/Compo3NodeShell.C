#include <cstdio>
#include <cstdlib>
#include <Utils.d/dbg_alloca.h>
#include <cmath>

#include <Element.d/CompShell.d/Compo3NodeShell.h>
#include <Element.d/State.h>
#include <Corotational.d/Shell3Corotator.h>
#include <Corotational.d/utilities.h>
#include <Corotational.d/GeomState.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/matrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/pstress.h>
#include <Hetero.d/InterpPoint.h>
#include <Utils.d/Conwep.d/BlastLoading.h>
#include <Utils.d/dbg_alloca.h>

extern "C" {
void _FORTRAN(compms)(double*, double*, double*, double*, 
                      const double&, double*, const int&, const int&, 
                      const int&, const int&, int*, double*,    
                      double*, const int&, const int&, const int&, 
                      const int&, double*, double*, const int&, 
                      double&, double&, const int&);

void _FORTRAN(compvms)(const int&, const int&, const int&, const int&,
                       const int&, double&, double&, double*,
                       double*, double*, double*, double*,
                       double*, double*, const int&, const int&,
                       const int&, double*, int*, double*, 
                       double*, const int&, int*, const int&,
                       const int&, const int&, const int&, const int&,
                       const int&, const int&, const int&, double&, double&);

void _FORTRAN(compst)(const double&, const int&, double*, double*,
                      const int&, const double&, double*, double*, 
                      double*, const int&, const int&, const int&, 
                      double*, int*, double*, double*, 
                      const int&, const int&, const int&, const int&,
                      int&);

void _FORTRAN(intes3)(const int&, const double&, const double&,
                      const int&, double*, double*);

void _FORTRAN(compthmfr)(double*, double*, double*, double&, double*, double*,
                         double*, double*, double*, double*, double*, double*, int&, 
			 double*, double *, const int&, int&);
}

extern int quietFlag;

bool Compo3NodeShell::Wzero_density = true;
bool Compo3NodeShell::Wthermal_force = true;

Compo3NodeShell::Compo3NodeShell(int* nodenums)
{
  nn[0] = nodenums[0];
  nn[1] = nodenums[1];
  nn[2] = nodenums[2];
  pbc = 0;
  type = 0;
  cFrame = 0;
  numLayers = 0;
}

Element *
Compo3NodeShell::clone()
{
  return new Compo3NodeShell(*this);
}

void
Compo3NodeShell::renum(const int *table)
{
  if(table[nn[0]] < 0 || table[nn[1]] < 0 || table[nn[2]] < 0)  {
    fprintf(stderr,"no mapping for these nodes %d %d %d \n",nn[0],nn[1],nn[2]);
    return;
  }
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
}

void
Compo3NodeShell::renum(EleRenumMap& table)
{
  if(table[nn[0]] < 0 || table[nn[1]] < 0 || table[nn[2]] < 0)  {
    fprintf(stderr,"no mapping for these nodes %d %d %d \n",nn[0],nn[1],nn[2]);
    return;
  }
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
}

void
Compo3NodeShell::getVonMises(Vector &stress, Vector &weight, CoordSet &cs,
                             Vector &elDisp, int strInd, int surface,
                             double *ndTemps, double ylayer, double zlayer, int avgnum)

{
  weight = 1.0;
  if(strInd == -1) return;

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);

  double x[3], y[3], z[3], h[3];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  h[0] = h[1] = h[2] = prop->eh;

  const int maxstr =  7;
  const int maxgus =  3;

  double elStress[3][7];
  double *laysigm = NULL;// new double[numLayers*maxstr*maxgus];

  int *laysgid = (int*)dbg_alloca(sizeof(int)*5*numLayers); 

  int i;
  for(i=0; i<5*numLayers; ++i)
    laysgid[i] = 0;

  int strainFlg = 0;
  if(strInd > 6 ) strainFlg = 1;

  double* disp = elDisp.data();

  int cfrm = (cFrame) ? 1 : 0;

  double thermalStrain1 = 0.0;
  double thermalStrain2 = 0.0; 
  if(ndTemps && type != 1) {
    if(type == 2 || type == 3) {
      // Compute the average thermal strain. Note that various approximations
      // are made in compvms so this additional approximation is not entirely
      // inappropriate IMHO. If accurate stress computation is required, use
      // element 15 instead of element 20.
      int numLayCoeff = 12;
      for(i=0; i<numLayers; ++i) {
        double dt = (ndTemps[0]+ndTemps[1]+ndTemps[2])/3 - prop->Ta; 
        thermalStrain1 += layData[numLayCoeff*i + 9]*dt;
        thermalStrain2 += layData[numLayCoeff*i + 10]*dt;
      }
      thermalStrain1 /= numLayers;
      thermalStrain2 /= numLayers;
    }
    else { // isotropic material
      double dt = (ndTemps[0]+ndTemps[1]+ndTemps[2])/3 - prop->Ta;
      thermalStrain1 = thermalStrain2 = (prop->W)*dt;
    }
  }

  _FORTRAN(compvms)(1, 1, 1, maxstr, maxgus, prop->E, prop->nu, h,
                    x, y, z, disp, (double*)elStress, (double*)laysigm,
                    1, numLayers, 1, (double*)cCoefs, (int*)idlay,
                    layData, cFrame, 0, laysgid, 1, type, 1, cfrm, numLayers,
                    1, strainFlg, surface, thermalStrain1, thermalStrain2);

  if(strInd < 7) {
    stress[0] = elStress[0][strInd];
    stress[1] = elStress[1][strInd];
    stress[2] = elStress[2][strInd];
  }
  else {
    stress[0] = elStress[0][strInd-7];
    stress[1] = elStress[1][strInd-7];
    stress[2] = elStress[2][strInd-7];
  }
}

void
Compo3NodeShell::getAllStress(FullM &stress, Vector &weight, CoordSet &cs,
                              Vector &elDisp, int strInd, int surface,
                              double *ndTemps)
{
  weight = 1.0;

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);

  double x[3], y[3], z[3], h[3];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  h[0] = h[1] = h[2] = prop->eh;

  const int maxstr =  7;
  const int maxgus =  3;

  double elStress[3][7];
  double *laysigm = NULL;// new double[numLayers*maxstr*maxgus];

  int *laysgid = (int*)dbg_alloca(sizeof(int)*5*numLayers);

  int i;
  for(i=0; i<5*numLayers; ++i)
    laysgid[i] = 0;

  double* disp = elDisp.data();

  int cfrm = (cFrame) ? 1 : 0;

  double thermalStrain1 = 0.0;
  double thermalStrain2 = 0.0;
  if(ndTemps && type != 1) {
    if(type == 2 || type == 3) { 
      // Compute the average thermal strain. Note that various approximations
      // are made in compvms so this additional approximation is not entirely
      // inappropriate IMHO. If accurate stress computation is required, use
      // element 15 instead of element 20.
      int numLayCoeff = 12;
      for(i=0; i<numLayers; ++i) {
        double dt = (ndTemps[0]+ndTemps[1]+ndTemps[2])/3 - prop->Ta;                   
        thermalStrain1 += layData[numLayCoeff*i + 9]*dt;
        thermalStrain2 += layData[numLayCoeff*i + 10]*dt;
      }
      thermalStrain1 /= numLayers;
      thermalStrain2 /= numLayers;
    }
    else { // isotropic material
      double dt = (ndTemps[0]+ndTemps[1]+ndTemps[2])/3 - prop->Ta;
      thermalStrain1 = thermalStrain2 = (prop->W)*dt;
    }
  }

  _FORTRAN(compvms)(1, 1, 1, maxstr, maxgus, prop->E, prop->nu, h,
                    x, y, z, disp, (double*)elStress, (double*)laysigm,
                    1, numLayers, 1, (double*)cCoefs, (int*)idlay, layData,
                    cFrame, 0, laysgid, 1, type, 1, cfrm, numLayers,
                    1, strInd, surface, thermalStrain1, thermalStrain2);

  // Store all Stress or all Strain as defined by strInd
  int j;
  for(i=0; i<3; ++i) {
    for(j=0; j<6; ++j) {
      stress[i][j] = elStress[i][j];
    }
  }

  // Get Element Principals
  double svec[6], pvec[3];
  for(j=0; j<6; ++j) {
    svec[j] = stress[0][j];
  }
  // Convert Engineering to Tensor Strains
  if(strInd != 0) {
    svec[3] /= 2;
    svec[4] /= 2;
    svec[5] /= 2;
  }

  pstress(svec,pvec);
  for(i=0; i<3; ++i) {
    for(j=0; j<3; ++j) {
      stress[i][j+6] = pvec[j];
    }
  }
}

void
Compo3NodeShell::getGravityForce(CoordSet& cs, double *gravityAcceleration, 
                                 Vector& gravityForce, int gravflg, GeomState *geomState)
{
  if (prop == NULL) {
    gravityForce.zero();
    return;
  }

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);

  double x[3], y[3], z[3], h[3], ElementMassMatrix[18][18];

  double grvfor[3];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  h[0] = h[1] = h[2] = prop->eh;

  int cfrm = (cFrame) ? 1 : 0;
  int grvflg = 1, masflg = 0;

  double totmas = 0;
  double area = 0;

  _FORTRAN(compms)(x, y, z, h, prop->rho, (double *)ElementMassMatrix,
                   18,numLayers, 1, 1, (int *)idlay, layData, cFrame,
                   1, type, 1, cfrm, gravityAcceleration, grvfor, grvflg,
                   totmas, area, masflg);

  // scale gravity force by number of nodes
  grvfor[0] /= 3.0;
  grvfor[1] /= 3.0;
  grvfor[2] /= 3.0;

  double mx[3],my[3],mz[3];
  int i;
  for(i=0; i<3; ++i) {
    mx[i]=0.0;
    my[i]=0.0;
    mz[i]=0.0;
  }

  // Lumped
  if(gravflg == 0) {

  }
  // Consistent or lumped with fixed end moments.  Compute treating shell as 3 beams.
  else {

    double T1[3],T2[3],T3[3];
    // Vector 1 from Node 1->2
    T1[0] = x[1] - x[0];
    T1[1] = y[1] - y[0];
    T1[2] = z[1] - z[0];
    normalize( T1 );
    // Vector 2 from Node 1->3
    T2[0] = x[2] - x[0];
    T2[1] = y[2] - y[0];
    T2[2] = z[2] - z[0];
    normalize( T2 );
    // Local Z-axis as cross between V1 and V2
    crossprod( T1, T2, T3 );
    normalize( T3);

    int beam, beamnode[3][2];
    beamnode[0][0] = 0;
    beamnode[0][1] = 1;
    beamnode[1][0] = 0;
    beamnode[1][1] = 2;
    beamnode[2][0] = 1;
    beamnode[2][1] = 2;

    for(beam=0; beam<3; ++beam) {
      double length, dx, dy, dz, localg[3];
      int n1, n2;
      n1 = beamnode[beam][0];
      n2 = beamnode[beam][1];
      dx = x[n2] - x[n1];
      dy = y[n2] - y[n1];
      dz = z[n2] - z[n1];
      length = sqrt(dx*dx + dy*dy + dz*dz);
      // Local X-axis from Node 1->2
      T1[0] = x[n2] - x[n1];
      T1[1] = y[n2] - y[n1];
      T1[2] = z[n2] - z[n1];
      normalize( T1 );
      // Local Y-axis as cross between Z and X
      crossprod( T3, T1, T2 );
      normalize( T2);

      for(i=0; i<3; ++i)
        localg[i] = 0.0;
      for(i=0; i<3; ++i) {
        localg[0] += T1[i]*grvfor[i];
        localg[1] += T2[i]*grvfor[i];
        localg[2] += T3[i]*grvfor[i];
      }
      double lmy,lmz;
      if (gravflg == 2) { // consistent
        lmy = -localg[2]*length/12.0;
        lmz = localg[1]*length/12.0;
      }
      else { // lumped with fixed-end moments
        lmy = -localg[2]*length/16.0;
        lmz = localg[1]*length/16.0;
      }
      mx[n1] += ((T2[0]*lmy) + (T3[0]*lmz));
      my[n1] += ((T2[1]*lmy) + (T3[1]*lmz));
      mz[n1] += ((T2[2]*lmy) + (T3[2]*lmz));
      mx[n2] -= ((T2[0]*lmy) + (T3[0]*lmz));
      my[n2] -= ((T2[1]*lmy) + (T3[1]*lmz));
      mz[n2] -= ((T2[2]*lmy) + (T3[2]*lmz));
    }
  }

  // set gravity force
  gravityForce[0]  = grvfor[0];
  gravityForce[1]  = grvfor[1];
  gravityForce[2]  = grvfor[2];
  gravityForce[3]  = mx[0];
  gravityForce[4]  = my[0];
  gravityForce[5]  = mz[0];
  gravityForce[6]  = grvfor[0];
  gravityForce[7]  = grvfor[1];
  gravityForce[8]  = grvfor[2];
  gravityForce[9]  = mx[1];
  gravityForce[10] = my[1];
  gravityForce[11] = mz[1];
  gravityForce[12] = grvfor[0];
  gravityForce[13] = grvfor[1];
  gravityForce[14] = grvfor[2];
  gravityForce[15] = mx[2];
  gravityForce[16] = my[2];
  gravityForce[17] = mz[2];
}

double
Compo3NodeShell::getMass(const CoordSet &cs) const
{ 
  if(prop == NULL) return 0.0;

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);

  double x[3], y[3], z[3], h[3], ElementMassMatrix[18][18];
  double *gravityAcceleration=NULL, *grvfor=NULL;

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  h[0] = h[1] = h[2] = prop->eh;

  int cfrm = (cFrame) ? 1 : 0;
  int grvflg = 0, masflg = 1;

  double totmas = 0.0;
  double area = 0.0;

  _FORTRAN(compms)(x, y, z, h, prop->rho, (double *)ElementMassMatrix,
                   18,numLayers, 1, 1, (int *)idlay, layData, cFrame,
                   1, type, 1, cfrm, gravityAcceleration, grvfor, grvflg,
                   totmas, area, masflg);

  return totmas;
}

double
Compo3NodeShell::getMassThicknessSensitivity(CoordSet &cs)
{ 
  if(prop == NULL) return 0.0;

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);

  double x[3], y[3], z[3], h[3], ElementMassMatrix[18][18];
  double *gravityAcceleration=NULL, *grvfor=NULL;

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  h[0] = h[1] = h[2] = prop->eh;

  int cfrm = (cFrame) ? 1 : 0;
  int grvflg = 0, masflg = 1;

  double totmas = 0.0;
  double area = 0.0;

  _FORTRAN(compms)(x, y, z, h, prop->rho, (double *)ElementMassMatrix,
                   18,numLayers, 1, 1, (int *)idlay, layData, cFrame,
                   1, type, 1, cfrm, gravityAcceleration, grvfor, grvflg,
                   totmas, area, masflg);

  return totmas/prop->eh;
}

FullSquareMatrix
Compo3NodeShell::massMatrix(const CoordSet &cs, double *mel, int cmflg) const
{
  if(prop == NULL) {
    FullSquareMatrix ret(18,mel);
    ret.zero();
    return ret;
  }

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);

  double x[3], y[3], z[3], h[3];

  double *gravityAcceleration=NULL, *grvfor=NULL;

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  h[0] = h[1] = h[2] = prop->eh;

  int cfrm = (cFrame) ? 1 : 0;
  int grvflg = 0, masflg = 0;
  double totmas = 0;
  double area = 0;

  // check if the density is negative or zero
  if(prop->rho <= 0.0 && (type == 0 || type == 1) && quietFlag == 0 && Wzero_density) {
    fprintf(stderr," *** WARNING: Composite shell element # %d has zero or negative density.\n"
                   "              Use command-line option -q to suppress this warning.\n", getGlNum()+1);
    Wzero_density = false;
  }

  _FORTRAN(compms)(x, y, z, h, prop->rho, (double *)mel,
                   18,numLayers, 1, 1, (int *)idlay, layData, cFrame,
                   1, type, 1, cfrm, gravityAcceleration, grvfor, grvflg,
                   totmas, area, masflg);

  FullSquareMatrix ret(18,mel);

  return ret;
}

void
Compo3NodeShell::setCompositeData(int _type, int nlays, double *lData,
                                  double *coefs, double *frame)
{
 type = _type;
 numLayers = nlays;
 if(nlays > 0) {
  idlay = new int[nlays][5];
  int il;
  for(il = 0; il < nlays; ++il) {
    idlay[il][0] = 1;
    idlay[il][1] = nlays;
    idlay[il][2] = il+1;
  }
  layData = lData;
 } else
   idlay = NULL;
 cFrame = frame;
 cCoefs = coefs;
}

#define PI 3.14159265358979

double *
Compo3NodeShell::setCompositeData2(int _type, int nlays, double *lData,
                                   double *coefs, CoordSet &cs, double theta)
{
 // variant where cframe is not pre-defined but calculated from nodal coordinates and angle theta
 // theta is the angle in degrees between node1-node2 and the material x axis
 type = _type;
 numLayers = nlays;
 if(nlays > 0) {
  idlay = new int[nlays][5];
  int il;
  for(il = 0; il < nlays; ++il) {
    idlay[il][0] = 1;
    idlay[il][1] = nlays;
    idlay[il][2] = il+1;
  }
  layData = lData;
 } else
   idlay = NULL;
 cCoefs = coefs;

 // compute cFrame
 cFrame = new double[9];
                                                                                                        
 auto &nd1 = cs.getNode(nn[0]);
 auto &nd2 = cs.getNode(nn[1]);
 auto &nd3 = cs.getNode(nn[2]);
 
 double ab[3], ac[3], x[3], y[3], z[3];
 ab[0] = nd2.x - nd1.x; ab[1] = nd2.y - nd1.y; ab[2] = nd2.z - nd1.z;
 ac[0] = nd3.x - nd1.x; ac[1] = nd3.y - nd1.y; ac[2] = nd3.z - nd1.z;
 // x = AB
 x[0] = ab[0]; x[1] = ab[1]; x[2] = ab[2];
 // z = AB cross AC
 z[0] = ab[1]*ac[2]-ab[2]*ac[1];
 z[1] = ab[2]*ac[0]-ab[0]*ac[2];
 z[2] = ab[0]*ac[1]-ab[1]*ac[0];
 // y = z cross x
 y[0] = z[1]*x[2]-z[2]*x[1];
 y[1] = z[2]*x[0]-z[0]*x[2];
 y[2] = z[0]*x[1]-z[1]*x[0];
 // rotate x and y about z
 theta *= PI/180.; // convert to radians
 double c = cos(theta), s = sin(theta);

 // use Rodrigues' Rotation Formula to rotate x and y about z by an angle theta
 double R[3][3];
 normalize(x); normalize(y); normalize(z); double wx = z[0], wy = z[1], wz = z[2];
 R[0][0] = c + wx*wx*(1-c);
 R[0][1] = wx*wy*(1-c)-wz*s;
 R[0][2] = wy*s+wx*wz*(1-c);
 R[1][0] = wz*s+wx*wy*(1-c);
 R[1][1] = c+wy*wy*(1-c);
 R[1][2] = -wx*s+wy*wz*(1-c);
 R[2][0] = -wy*s+wx*wz*(1-c);
 R[2][1] = wx*s+wy*wz*(1-c);
 R[2][2] = c+wz*wz*(1-c);

 cFrame[0] = R[0][0]*x[0] + R[0][1]*x[1] + R[0][2]*x[2];
 cFrame[1] = R[1][0]*x[0] + R[1][1]*x[1] + R[1][2]*x[2];
 cFrame[2] = R[2][0]*x[0] + R[2][1]*x[1] + R[2][2]*x[2];
 cFrame[3] = R[0][0]*y[0] + R[0][1]*y[1] + R[0][2]*y[2];
 cFrame[4] = R[1][0]*y[0] + R[1][1]*y[1] + R[1][2]*y[2];
 cFrame[5] = R[2][0]*y[0] + R[2][1]*y[1] + R[2][2]*y[2];
 cFrame[6] = z[0];
 cFrame[7] = z[1];
 cFrame[8] = z[2];
 
 return cFrame;
}

void
Compo3NodeShell::getCFrame(CoordSet &cs, double cFrame[3][3]) const
{
  if(Compo3NodeShell::cFrame) {
    fprintf(stderr," *** WARNING: Compo3NodeShell::getCFrame is not implemented\n");
  }
  cFrame[0][0] = cFrame[1][1] = cFrame[2][2] = 1.;
  cFrame[0][1] = cFrame[0][2] = cFrame[1][0] = cFrame[1][2] = cFrame[2][0] = cFrame[2][1] = 0.;
}

FullSquareMatrix
Compo3NodeShell::stiffness(const CoordSet &cs, double *d, int flg) const
{
  if(prop == NULL) {
    FullSquareMatrix ret(18,d);
    ret.zero();
    return ret;
  }

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);

  double x[3], y[3], z[3], h[3];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
  h[0] = h[1] = h[2] = prop->eh;

  // check if the thickness is negative or zero
  if(h[0] <= 0.0 && type == 0) {
    fprintf(stderr," *** ERROR: Composite shell element # %d has zero or negative thickness. Exiting...\n", getGlNum()+1);
    exit(-1);
  }

  int cfrm = (cFrame) ? 1 : 0;

  _FORTRAN(compst)(prop->E, 1, h, (double *)d, 18,
                   prop->nu, x, y, z, 1, numLayers, 1, (double*) cCoefs,
                   (int *)idlay,
                   layData, cFrame, 1, type, 1, cfrm, flg);

  FullSquareMatrix ret(18,d);

  return ret;
}

int
Compo3NodeShell::numNodes() const
{
  return 3;
}

int*
Compo3NodeShell::nodes(int *p) const
{
  if(p == 0) p = new int[3];
  p[0] = nn[0];
  p[1] = nn[1];
  p[2] = nn[2];
  return p;
}

int
Compo3NodeShell::numDofs() const
{
  return 18;
}

int*
Compo3NodeShell::dofs(DofSetArray &dsa, int *p) const
{
  if(p == 0) p = new int[18];

  dsa.number(nn[0],DofSet::XYZdisp | DofSet::XYZrot, p  );
  dsa.number(nn[1],DofSet::XYZdisp | DofSet::XYZrot, p+6);
  dsa.number(nn[2],DofSet::XYZdisp | DofSet::XYZrot, p+12);

  return p;
}

void
Compo3NodeShell::markDofs(DofSetArray &dsa) const
{
  dsa.mark(nn, 3, DofSet::XYZdisp | DofSet::XYZrot);
}

void
Compo3NodeShell::
computeDisp(CoordSet&, State &state, const InterpPoint &ip,
            double *res, GeomState *gs)
{
  const double *gp = ip.xy;
  double xyz[3][6];
  state.getDV(nn[0], xyz[0], xyz[0]+3);
  state.getDV(nn[1], xyz[1], xyz[1]+3);
  state.getDV(nn[2], xyz[2], xyz[2]+3);

  int j;
  for(j = 0; j < 6; ++j)
    res[j] = (1-gp[0]-gp[1]) * xyz[0][j] + gp[0]*xyz[1][j] + gp[1]*xyz[2][j];
}

void
Compo3NodeShell::getFlLoad(CoordSet &, const InterpPoint &ip, double *flF, 
                           double *resF, GeomState *gs)
{
  const double *gp = ip.xy;
  int i;
  for(i = 0; i < 3; ++i) {
    resF[i+3] = resF[i+9] = resF[i+15] = 0;
    resF[i] = (1-gp[0]-gp[1]) * flF[i];
    resF[6+i] = gp[0] * flF[i];
    resF[12+i] = gp[1] * flF[i];
  }
}

Corotator *
Compo3NodeShell::getCorotator(CoordSet &cs, double *kel, int fitAlgShell, int)
{
  int flag = 0;
  FullSquareMatrix myStiff = stiffness(cs, kel, flag);
  return new Shell3Corotator(nn[0],nn[1],nn[2],myStiff, fitAlgShell);
}

int
Compo3NodeShell::getTopNumber() const
{
  return 120;
}

void
Compo3NodeShell::computePressureForce(CoordSet& cs, Vector& elPressureForce,
                                      GeomState *geomState, int cflg, double time) {
     double pressure = pbc->val;
     // Check if Conwep is being used. If so, add the pressure from the blast loading function.
     if (pbc->conwep && pbc->conwepswitch) {
       double* CurrentElementNodePositions = (double*) dbg_alloca(sizeof(double)*3*4);
       int Offset;
       for(int i = 0; i < 4; ++i) {
         Offset = i*3;
         if (i==3) {
           CurrentElementNodePositions[Offset+0] = cs[nn[2]]->x;
           CurrentElementNodePositions[Offset+1] = cs[nn[2]]->y;
           CurrentElementNodePositions[Offset+2] = cs[nn[2]]->z;
         }
         else {
           CurrentElementNodePositions[Offset+0] = cs[nn[i]]->x;
           CurrentElementNodePositions[Offset+1] = cs[nn[i]]->y;
           CurrentElementNodePositions[Offset+2] = cs[nn[i]]->z;
         }
       }
       pressure += BlastLoading::ComputeShellPressureLoad(CurrentElementNodePositions, time, *(pbc->conwep));
     }
     double px = 0.0;
     double py = 0.0;
     double pz = 0.0;

     double mx[3],my[3],mz[3];
     int i;
     for(i=0; i<3; ++i) {
       mx[i]=0.0;
       my[i]=0.0;
       mz[i]=0.0;
     }

     // Compute area of shell
     auto &nd1 = cs.getNode(nn[0]);
     auto &nd2 = cs.getNode(nn[1]);
     auto &nd3 = cs.getNode(nn[2]);

     double r1[3], r2[3], r3[3], v1[3], v2[3], normal[3];

     r1[0] = nd1.x; r1[1] = nd1.y; r1[2] = nd1.z;
     r2[0] = nd2.x; r2[1] = nd2.y; r2[2] = nd2.z;
     r3[0] = nd3.x; r3[1] = nd3.y; r3[2] = nd3.z;

     v1[0] = r3[0] - r1[0];
     v1[1] = r3[1] - r1[1];
     v1[2] = r3[2] - r1[2];

     v2[0] = r2[0] - r1[0];
     v2[1] = r2[1] - r1[1];
     v2[2] = r2[2] - r1[2];

     // Compute normal to shell using vector cross product rule
     crossprod(v2, v1, normal);

     double magnitude = sqrt(normal[0]*normal[0] + normal[1]*normal[1]
                                                  + normal[2]*normal[2]);
     double area = 0.5*magnitude;

     // compute pressure force per node
     double pressureForce = pressure * area / 3.0;

     double xn[3][3];
       for(i=0; i<3; ++i) {
         xn[0][i] = r1[i];
         xn[1][i] = r2[i];
         xn[2][i] = r3[i];
       }

     double mcoef = 8; // 12 ???

     if(!geomState) {

       // compute unit normal to shell surface

       normal[0] /= magnitude;
       normal[1] /= magnitude;
       normal[2] /= magnitude;

       px = pressureForce*normal[0];
       py = pressureForce*normal[1];
       pz = pressureForce*normal[2];

       if (cflg == 1) {

         int beam, beamnode[3][2];
         beamnode[0][0] = 0;
         beamnode[0][1] = 1;
         beamnode[1][0] = 0;
         beamnode[1][1] = 2;
         beamnode[2][0] = 1;
         beamnode[2][1] = 2;

         for(beam=0; beam<3; ++beam) {
           double length, dx, dy, dz;
           int n1, n2;
           n1 = beamnode[beam][0];
           n2 = beamnode[beam][1];
           dx = xn[n2][0] - xn[n1][0];
           dy = xn[n2][1] - xn[n1][1];
           dz = xn[n2][2] - xn[n1][2];
           length = sqrt(dx*dx + dy*dy + dz*dz);
           // Local X-axis from Node 1->2
           for(i=0; i<3; i++ ) v1[i] = xn[n2][i] - xn[n1][i];
           normalize( v1 );
           // Local Y-axis as cross between Z and X
           crossprod( normal, v1, v2 );
           normalize( v2 );

           double lmy = -pressureForce*length/mcoef; // mcoef = 8 seems to work better for cylinder shell model

           mx[n1] += (v2[0]*lmy);
           my[n1] += (v2[1]*lmy);
           mz[n1] += (v2[2]*lmy);
           mx[n2] -= (v2[0]*lmy);
           my[n2] -= (v2[1]*lmy);
           mz[n2] -= (v2[2]*lmy);
         }
       }

       elPressureForce[0]  = px;
       elPressureForce[1]  = py;
       elPressureForce[2]  = pz;
       elPressureForce[3]  = mx[0];
       elPressureForce[4]  = my[0];
       elPressureForce[5]  = mz[0];

       elPressureForce[6]  = px;
       elPressureForce[7]  = py;
       elPressureForce[8]  = pz;
       elPressureForce[9]  = mx[1];
       elPressureForce[10] = my[1];
       elPressureForce[11] = mz[1];

       elPressureForce[12] = px;
       elPressureForce[13] = py;
       elPressureForce[14] = pz;
       elPressureForce[15] = mx[2];
       elPressureForce[16] = my[2];
       elPressureForce[17] = mz[2];
     }

     //if local coordinates are needed for nonlinear analysis
     if (geomState) {

        //compute centroid
        double xc0[3];
        double t0[3][3],xl0[3][3];
        int nod;
        for (i=0;i<3;i++)
        xc0[i] = ( xn[0][i] + xn[1][i] + xn[2][i] )/3.0;

        // Compute t0 transformation matrix with x axis along side 1-2 
        for( i=0; i<3; i++ ) t0[0][i] = xn[1][i] - xn[0][i];
        normalize( t0[0] );

        // local y axis
        for( i=0; i<3; i++ ) t0[1][i] = xn[2][i] - xn[0][i];
        crossprod( t0[0], t0[1], t0[2] );
        normalize( t0[2] );

        // local z axis
        crossprod( t0[2], t0[0], t0[1] );
        normalize( t0[1] );

        // Compute local coordinates of undeformed element
        for( nod=0; nod<3; nod++) {
             for( i=0; i<3; i++ ) {
                 xl0[nod][i] = t0[i][0]*(xn[nod][0] - xc0[0])
                              +t0[i][1]*(xn[nod][1] - xc0[1])
                              +t0[i][2]*(xn[nod][2] - xc0[2]);
                }
             }

        elPressureForce[0]  = 0;
        elPressureForce[1]  = 0;
        elPressureForce[2]  = pressureForce;
        elPressureForce[3]  = static_cast<double>(cflg)*pressureForce*( xl0[2][1] - xl0[0][1]
                                                  +xl0[1][1] - xl0[0][1])/mcoef;
        elPressureForce[4]  = static_cast<double>(cflg)*pressureForce*( xl0[0][0] - xl0[2][0]
                                                  +xl0[0][0] - xl0[1][0])/mcoef;
        elPressureForce[5]  = 0;

        elPressureForce[6]  = 0;
        elPressureForce[7]  = 0;
        elPressureForce[8]  = pressureForce;
        elPressureForce[9]  = static_cast<double>(cflg)*pressureForce*( xl0[0][1] - xl0[1][1]
                                                  +xl0[2][1] - xl0[1][1])/mcoef;
        elPressureForce[10] = static_cast<double>(cflg)*pressureForce*( xl0[1][0] - xl0[0][0]
                                                  +xl0[1][0] - xl0[2][0])/mcoef;
        elPressureForce[11] = 0;

        elPressureForce[12] = 0;
        elPressureForce[13] = 0;
        elPressureForce[14] = pressureForce;
        elPressureForce[15] = static_cast<double>(cflg)*pressureForce*( xl0[1][1] - xl0[2][1]
                                                  +xl0[0][1] - xl0[2][1])/mcoef;
        elPressureForce[16] = static_cast<double>(cflg)*pressureForce*( xl0[2][0] - xl0[1][0]
                                                  +xl0[2][0] - xl0[0][0])/mcoef;
        elPressureForce[17] = 0;

     }

}

double *
Compo3NodeShell::getMidPoint(CoordSet &cs)
{ 
  double * midPoint = new double[3];

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);

  midPoint[0] = ( nd1.x + nd2.x + nd3.x ) / 3.0; 
  midPoint[1] = ( nd1.y + nd2.y + nd3.y ) / 3.0; 
  midPoint[2] = ( nd1.z + nd2.z + nd3.z ) / 3.0; 

  return midPoint;
}

//--------------------------------------------------------------------------

void
Compo3NodeShell::getThermalForce(CoordSet& cs, Vector& ndTemps,
				Vector &elThermalForce, int glflag, 
				GeomState *)
{  
  //check to see that the coefficent of thermal expansions will exist 
  if(prop == NULL || type == 1) {
    if(type == 1 && quietFlag == 0 && Wthermal_force) {
      fprintf(stderr," *** WARNING: Thermal forces are not computed for composite shell element if\n" 
                     "              the constitutive matrix is inputted using the COEF sub-command.\n"
                     "              Element type 20 (2020) should be changed to type 15 (1515).\n"
                     "              Use command-line option -q to suppress this warning.\n");
      Wthermal_force = false;
    }
    elThermalForce.zero();
    return;
  }

  int numlay = 1;
  if(type != 0) numlay = numLayers;
   
  int i,j;
  double x[3],y[3],z[3];

  //get nodal temperatures and relate to reference temperature
   
  double meant = 0.0; //determine the average nodal temperature
   
  for(i=0; i<3; i++) meant += ndTemps[i]/3;

  double* thk   = (double*) alloca(sizeof(double)*numlay);
  double* emk1  = (double*) alloca(sizeof(double)*numlay);
  double* emk2  = (double*) alloca(sizeof(double)*numlay);
  double* nuk   = (double*) alloca(sizeof(double)*numlay);
  double* ctek1 = (double*) alloca(sizeof(double)*numlay);
  double* ctek2 = (double*) alloca(sizeof(double)*numlay);
  double* phik  = (double*) alloca(sizeof(double)*numlay);
  double* tak   = (double*) alloca(sizeof(double)*numlay);

  //create a vector which contains the thicknesses, Young's modulus, 
  //Poisson's ratio and the coefficent of thermal expansion of each layer 
  //the sign convention here for zvector is opposite that of the reference 
  //Note: the reference temperature is the same for each layer.
  int numLayCoeff = 12;
  if (type != 0) {
    for(i=0; i<numLayers; ++i) {
      emk1[i]  = layData[numLayCoeff*i    ];
      emk2[i]  = layData[numLayCoeff*i + 1];
      nuk[i]   = layData[numLayCoeff*i + 2];
      thk[i]   = layData[numLayCoeff*i + 7];
      ctek1[i] = layData[numLayCoeff*i + 9];
      ctek2[i] = layData[numLayCoeff*i +10];
      phik[i]  = layData[numLayCoeff*i + 8];
      tak[i]   = prop->Ta;
    }
  }  
  else {
    //isotropic material
    thk[0] = prop->eh;
    emk1[0] = prop->E;
    emk2[0] = prop->E;
    nuk[0] = prop->nu;
    ctek1[0] = prop->W;
    ctek2[0] = prop->W;
    phik[0] = 0.0;
    tak[0] = prop->Ta;
  }
 
  // Compute Node's coordinates
  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);
  
  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  // Now compute the elemental load vectors 
  // Note: Unlike element 15, this element assumes that the shear thermal expansion coefficient w12 is zero.
 
  elThermalForce.zero();

  int cfrm = (cFrame) ? 1 : 0;
     
  _FORTRAN(compthmfr)(x, y, z, meant, ctek1, ctek2, emk1, emk2, nuk, thk, phik, tak, 
                      numlay, cFrame, (double *)elThermalForce.data(), cfrm, glflag);
}
