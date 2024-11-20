#ifdef USE_EIGEN3
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <iomanip>
#include <iostream>

#include <Corotational.d/GeomState.h>
#include <Corotational.d/Shell3Corotator.h>
#include <Corotational.d/utilities.h>
#include <Element.d/FelippaShell.d/FelippaShell.h>
#include <Element.d/NonLinearity.d/ExpMat.h>
#include <Element.d/NonLinearity.d/MaterialWrapper.h>
#include <Element.d/NonLinearity.d/NLMembrane.h>
#include <Element.d/Utils.d/SolidElemUtils.h>
#include <Element.d/State.h>
#include <Hetero.d/InterpPoint.h>
#include <Material.d/IsotropicLinearElasticJ2PlasticPlaneStressMaterial.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/matrix.h>
#include <Math.d/BLAS.h>
#include <Utils.d/Conwep.d/BlastLoading.h>
#include <Utils.d/dbg_alloca.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/pstress.h>
#include <Utils.d/SolverInfo.h>

#include <Element.d/FelippaShell.d/ShellElementTemplate.hpp>
#include <Element.d/FelippaShell.d/EffMembraneTriangle.hpp>
#include <Element.d/FelippaShell.d/AndesBendingTriangle.hpp>

typedef ShellElementTemplate<double,EffMembraneTriangle,AndesBendingTriangle> Impl;

extern int verboseFlag;
extern SolverInfo &solInfo;

int FelippaShell::sflg = 1;
int FelippaShell::tflg = 1;

FelippaShell::FelippaShell(int* nodenums)
{
  nn[0] = nodenums[0];
  nn[1] = nodenums[1];
  nn[2] = nodenums[2];
  type = 0;
  cFrame = 0;
  pbc = 0;
  gpmat = 0;
  nmat = 0;
}

FelippaShell::~FelippaShell()
{
  if(nmat && nmat != gpmat) delete nmat;
  if(gpmat) delete gpmat;
}

Element *
FelippaShell::clone()
{
  return new FelippaShell(*this);
}

void
FelippaShell::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
}

void
FelippaShell::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
}

void
FelippaShell::getVonMises(Vector &stress, Vector &weight, CoordSet &cs,
                          Vector &elDisp, int strInd, int surface,
                          double *ndTemps, double ylayer, double zlayer,
                          int avgnum)

{
  getVonMisesImpl(stress, weight, cs, elDisp, strInd, surface, ndTemps, ylayer, zlayer, 1);
}

void
FelippaShell::getNLVonMises(Vector& stress, Vector& weight, GeomState &curState,
                            GeomState *refState, CoordSet& c0, int strIndex, int surface,
                            double ylayer, double zlayer, int avgnum, int measure)
{
  Vector elDisp(numDofs(), 0.);

  // Get Nodes original coordinates (C0 configuration)
  Node &node1 = c0.getNode( n1 );
  Node &node2 = c0.getNode( n2 );
  Node &node3 = c0.getNode( n3 );

  // Get Nodes current coordinates (C0n configuration)
  NodeState &ns1 = curState[ n1 ];
  NodeState &ns2 = curState[ n2 ];
  NodeState &ns3 = curState[ n3 ];

  double xl0[3][3], xln[3][3], t0[3][3], t0n[3][3];

  // Extract deformational displacement from C0 to C0n configurations 
  extractDefDisp(node1,node2,node3, ns1,ns2,ns3, xl0,xln, t0,t0n, elDisp.data());

  double *staten = refState->getElemState(getGlNum()) + subNum*numStates() + gpmat->GetNumStates();
  double *statenp = curState.getElemState(getGlNum()) + subNum*numStates() + gpmat->GetNumStates();

  getVonMisesImpl(stress, weight, c0, elDisp, strIndex, surface, (double*)0, ylayer, zlayer, 0, staten, statenp);
}

void
FelippaShell::extractDeformations(GeomState &geomState, CoordSet &cs, double *vld,
                                  int &nlflag)
{
  if(type == 4) nlflag = 2;
  else Shell3Corotator::extractDeformations(geomState, cs, vld, nlflag);
}

void
FelippaShell::getVonMisesImpl(Vector &stress, Vector &weight, CoordSet &cs,
                              Vector &elDisp, int strInd, int surface,
                              double *ndTemps, double ylayer, double zlayer,
                              int flag, double *staten, double *statenp)
{
  weight = 1.0;
  if(strInd == -1) return;

  int strainFlg, offset;
  switch(strInd) {
    case 0 : case 1 : case 2 : case 3 : case 4 : case 5 : case 6 : {
      strainFlg = 0;
      offset = 0;
    } break;
    case 7 : case 8 : case 9 : case 10: case 11: case 12: case 13 : {
      strainFlg = 1;
      offset = 7;
    } break;
    case 18 : {
      strainFlg = 2;
      offset = 18;
    } break;
    case 19 : case 20 : case 21 : case 22 : case 23 : case 24 : {
      strainFlg = 3;
      offset = 19;
    } break;
    case 25 : case 26 : case 27 : case 28 : case 29 : case 30 : {
      strainFlg = 4;
      offset = 25;
    } break;
   case 17 : {
      strainFlg = 5;
      offset = 17;
    } break;
    default : {
      stress.zero();
      return;
    } 
  }

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);

  double x[3], y[3], z[3];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  int maxstr = 7;

  double elStress[3][7];

  double* disp = elDisp.data();

  Impl::andesvms(glNum+1, maxstr, prop->nu, x, y, z, disp,
                 (double*)elStress, type, nmat, strainFlg,
                 surface, sflg, ndTemps, flag, staten, statenp);

  stress[0] = elStress[0][strInd-offset];
  stress[1] = elStress[1][strInd-offset];
  stress[2] = elStress[2][strInd-offset];
}

void
FelippaShell::getAllStress(FullM &stress, Vector &weight, CoordSet &cs,
                           Vector &elDisp, int strInd, int surface,
                           double *ndTemps)
{
  getAllStressImpl(stress, weight, cs, elDisp, strInd, surface, ndTemps, 1);
}

void
FelippaShell::getNLAllStress(FullM &stress, Vector &weight, GeomState &curState,
                             GeomState *refState, CoordSet &c0, int strInd, int surface,
                             int measure)
{
  Vector elDisp(numDofs(), 0.);

  // Get Nodes original coordinates (C0 configuration)
  Node &node1 = c0.getNode( n1 );
  Node &node2 = c0.getNode( n2 );
  Node &node3 = c0.getNode( n3 );

  // Get Nodes current coordinates (C0n configuration)
  NodeState &ns1 = curState[ n1 ];
  NodeState &ns2 = curState[ n2 ];
  NodeState &ns3 = curState[ n3 ];

  double xl0[3][3], xln[3][3], t0[3][3], t0n[3][3];

  // Extract deformational displacement from C0 to C0n configurations 
  extractDefDisp(node1,node2,node3, ns1,ns2,ns3, xl0,xln, t0,t0n, elDisp.data());

  double *staten = refState->getElemState(getGlNum()) + subNum*numStates() + gpmat->GetNumStates();
  double *statenp = curState.getElemState(getGlNum()) + subNum*numStates() + gpmat->GetNumStates();

  getAllStressImpl(stress, weight, c0, elDisp, strInd, surface, (double*)0, 0, staten, statenp);
}

void
FelippaShell::getAllStressImpl(FullM &stress, Vector &weight, CoordSet &cs,
                               Vector &elDisp, int strInd, int surface,
                               double *ndTemps, int flag, double *staten, double *statenp)
{
  weight = 1.0;

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);

  double x[3], y[3], z[3];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  int maxstr = 7;

  double elStress[3][7];

  double* disp = elDisp.data();

  Impl::andesvms(glNum+1, maxstr, prop->nu, x, y, z, disp,
                 (double*)elStress, type, nmat, strInd,
                 surface, sflg, ndTemps, flag, staten, statenp);

  // Store all Stress or all Strain as defined by strInd
  int i,j;
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
FelippaShell::getGravityForce(CoordSet& cs, double *gravityAcceleration, 
                              Vector& gravityForce, int gravflg, GeomState *geomState)
{
  if (prop == NULL) {
    gravityForce.zero();
    return;
  }

  double x[3] = { cs[nn[0]]->x, cs[nn[1]]->x, cs[nn[2]]->x };
  double y[3] = { cs[nn[0]]->y, cs[nn[1]]->y, cs[nn[2]]->y };
  double z[3] = { cs[nn[0]]->z, cs[nn[1]]->z, cs[nn[2]]->z };

  double rhoh = gpmat->GetAreaDensity();

  Impl::andesgf(glNum+1, x, y, z, gravityForce.data(), gravityAcceleration, gravflg, rhoh);
}

double
FelippaShell::getMass(const CoordSet &cs) const
{ 
  if(prop == NULL) return 0.0;

  double x[3] = { cs[nn[0]]->x, cs[nn[1]]->x, cs[nn[2]]->x };
  double y[3] = { cs[nn[0]]->y, cs[nn[1]]->y, cs[nn[2]]->y };
  double z[3] = { cs[nn[0]]->z, cs[nn[1]]->z, cs[nn[2]]->z };
  double *gravityAcceleration = NULL, *grvfor = NULL;
  double totmas = 0.0;
  double rhoh = gpmat->GetAreaDensity();

  Impl::andesms(glNum+1, x, y, z, gravityAcceleration, grvfor, totmas, rhoh);

  return totmas;
}

FullSquareMatrix
FelippaShell::massMatrix(const CoordSet &cs, double *mel, int cmflg) const
{
  if(prop == NULL) {
    FullSquareMatrix ret(18,mel);
    ret.zero();
    return ret;
  }

  double x[3] = { cs[nn[0]]->x, cs[nn[1]]->x, cs[nn[2]]->x };
  double y[3] = { cs[nn[0]]->y, cs[nn[1]]->y, cs[nn[2]]->y };
  double z[3] = { cs[nn[0]]->z, cs[nn[1]]->z, cs[nn[2]]->z };
  double rhoh = gpmat->GetAreaDensity();
  int mflg = (solInfo.useMassAugmentation) ? 1 : 0;

  Impl::andesmm(glNum+1, x, y, z, mel, rhoh, mflg);

  FullSquareMatrix ret(18, mel);

  return ret;
}

void
FelippaShell::setProp(StructProp *p, bool myProp)
{
  Element::setProp(p, myProp);
  type = 0;
  if(p) {
    nmat = gpmat = new ShellMaterialType0<double>(p->E, p->eh, p->nu, p->rho, p->Ta, p->W);
  }
  else nmat = gpmat = 0; // phantom
}

void
FelippaShell::setCompositeData(int _type, int nlays, double *lData,
                               double *coefs, double *frame)
{
  type = _type;
  cFrame = frame;

  if(gpmat) delete gpmat; // delete default material instantiated in setProp
  switch(type) {

    case 1 :
      nmat = gpmat = new ShellMaterialType1<double>(coefs, frame, prop->rho, prop->eh, prop->Ta);
      break;

    case 2 :
      nmat = gpmat = new ShellMaterialTypes2And3<double>(nlays, lData, false, frame, prop->Ta, prop->rho);
      break;

    case 3 :
      nmat = gpmat = new ShellMaterialTypes2And3<double>(nlays, lData, true, frame, prop->Ta, prop->rho);
      break;

    case 5 :
      nmat = gpmat = new ShellMaterialType5<double>(coefs, frame, prop->rho, prop->eh, prop->Ta);
      break;

    default :
      nmat = gpmat = 0;
      throw std::runtime_error("\n"
          "*** FATAL ERROR in FelippaShell::setCompositeData ***\n"
          "*** Wrong type of constitutive law                ***\n"
          "*** Types allowed are:                            ***\n"
          "*** 1 = given constitutive law                    ***\n"
          "*** 2 = given layers properties                   ***\n"
          "***     (no coupling bending/membrane)            ***\n"
          "*** 3 = given layers properties                   ***\n"
          "***     (with coupling bending/membrane)          ***\n"
          "*** 5 = given orthotropic constitutive law        ***\n"
          "*** STOP ALL TREATMENTS RIGHT HERE                ***\n");

  }

}

#define PI 3.14159265358979

double *
FelippaShell::setCompositeData2(int _type, int nlays, double *lData,
                                double *coefs, CoordSet &cs, double theta)
{
  // PJSA: variant where cframe is not pre-defined but calculated from nodal coordinates and angle theta
  // theta is the angle in degrees between node1-node2 and the material x axis
  type = _type;

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

  // use Rodrigues' Rotation Formula to rotation x and y about z by an angle theta
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

  if(gpmat) delete gpmat;
  switch(type) {

    case 1 :
      nmat = gpmat = new ShellMaterialType1<double>(coefs, cFrame, prop->rho, prop->eh, prop->Ta);
      break;

    case 2 :
      nmat = gpmat = new ShellMaterialTypes2And3<double>(nlays, lData, false, cFrame, prop->Ta, prop->rho);
      break;

    case 3 :
      nmat = gpmat = new ShellMaterialTypes2And3<double>(nlays, lData, true, cFrame, prop->Ta, prop->rho);
      break;

    case 5 :
      nmat = gpmat = new ShellMaterialType5<double>(coefs, cFrame, prop->rho, prop->eh, prop->Ta);
      break;

    default :
      nmat = gpmat = 0;
      throw std::runtime_error("\n"
          "*** FATAL ERROR in FelippaShell::setCompositeData2 ***\n"
          "*** Wrong Type of constitutive law                 ***\n"
          "*** Types allowed are:                             ***\n"
          "*** 1 = given constitutive law                     ***\n"
          "*** 2 = given layers properties                    ***\n"
          "***     (no coupling bending/membrane)             ***\n"
          "*** 3 = given layers properties                    ***\n"
          "***     (with coupling bending/membrane)           ***\n"
          "*** 5 = given orthotropic constitutive law         ***\n"
          "*** STOP ALL TREATMENTS RIGHT HERE                 ***\n");

  }
 
  return cFrame;
}

void
FelippaShell::getCFrame(CoordSet &cs, double cFrame[3][3]) const
{
  if(FelippaShell::cFrame) {
    auto &nd1 = cs.getNode(nn[0]);
    auto &nd2 = cs.getNode(nn[1]);
    auto &nd3 = cs.getNode(nn[2]);

    double x[3], y[3], z[3];

    x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
    x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
    x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

    Impl::andesfrm(glNum+1, x, y, z, FelippaShell::cFrame, &cFrame[0][0]);
  }
  else {
    cFrame[0][0] = cFrame[1][1] = cFrame[2][2] = 1.;
    cFrame[0][1] = cFrame[0][2] = cFrame[1][0] = cFrame[1][2] = cFrame[2][0] = cFrame[2][1] = 0.;
  }
}

void
FelippaShell::setMaterial(NLMaterial *_mat)
{
  if(!prop) return; // phantom element
  ExpMat *expmat = dynamic_cast<ExpMat *>(_mat);
  if(expmat && expmat->optctv == 5) { // deprecated
    double E = expmat->ematpro[0], nu = expmat->ematpro[1];
    double lambda = E*nu/((1+nu)*(1-2*nu)), mu = E/(2*(1+nu));
    double sigmaY = expmat->ematpro[3], K = expmat->ematpro[4], H = expmat->ematpro[5];
    double tol = expmat->ematpro[6];
    double epsF = (expmat->ematpro[7] <= 0) ? std::numeric_limits<double>::infinity() : expmat->ematpro[7];
    IsotropicLinearElasticJ2PlasticPlaneStressMaterial *mat = new IsotropicLinearElasticJ2PlasticPlaneStressMaterial(lambda, mu, sigmaY, K, H, tol, epsF);
    if(expmat->ysst) {
      for(int i=0; i<expmat->ysst->getNumPoints(); ++i)
        mat->SetExperimentalCurveData(expmat->ysst->getT(i), expmat->ysst->getV(i));
    }
    if(expmat->yssrt) {
      for(int i=0; i<expmat->yssrt->getNumPoints(); ++i)
        mat->SetExperimentalCurveData2(expmat->yssrt->getT(i), expmat->yssrt->getV(i));
    }
    type = 4;
    if(gpmat) delete gpmat;
    gpmat = new ShellMaterialType4<double,IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(prop->eh, prop->nu, prop->rho, mat, 5, 3);
    nmat = new ShellMaterialType4<double,IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(prop->eh, prop->nu, prop->rho, mat, 3, 3);
    delete mat;
  }
  else {
    MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial> *mat 
      = dynamic_cast<MaterialWrapper<IsotropicLinearElasticJ2PlasticPlaneStressMaterial> *>(_mat);
    if(mat) {
      type = 4;
      if(gpmat) delete gpmat;
      gpmat = new ShellMaterialType4<double,IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(prop->eh, prop->nu, prop->rho, mat->getMaterial(), 5, 3);
      nmat = new ShellMaterialType4<double,IsotropicLinearElasticJ2PlasticPlaneStressMaterial>(prop->eh, prop->nu, prop->rho, mat->getMaterial(), 3, 3);
    }
    else {
      if(_mat->getGenStrainEvaluator() && !_mat->getGenStrainEvaluator()->isNonLinear()) {
        type = 4;
        double *cCoefs = (gpmat && gpmat->GetCoefOfConstitutiveLaw()) ? gpmat->GetCoefOfConstitutiveLaw() : 0;
        if(cCoefs) {
          NLMaterial *mat = _mat->clone();
          if(mat) {
            double C[6][6], alpha[6];
            // transform local constitutive matrix to global frame
            rotateConstitutiveMatrix(cCoefs, cFrame, C);
            mat->setTangentMaterial(C);
            // transform local coefficients of thermal expansion to global frame
            rotateVector(cCoefs+36, cFrame, alpha);
            mat->setThermalExpansionCoef(alpha);
          }
          else {
            fprintf(stderr, " *** ERROR: Unsupported material law specified for element %d.\n", getGlNum()+1);
            exit(-1);
          }
          if(gpmat) delete gpmat;
          gpmat = new ShellMaterialType6<double,NLMaterial>(mat, 5, 3);
          nmat = new ShellMaterialType6<double,NLMaterial>(mat, 3, 3);
        }
        else {
          if(gpmat) delete gpmat;
          gpmat = new ShellMaterialType6<double,NLMaterial>(_mat, 5, 3);
          nmat = new ShellMaterialType6<double,NLMaterial>(_mat, 3, 3);
        }
      }
      else {
        fprintf(stderr, " *** ERROR: Unsupported material law specified for element %d.\n", getGlNum()+1);
        exit(-1);
      }
    }
  }
}

int
FelippaShell::numStates()
{
  if(prop == NULL) return 0;
  return gpmat->GetNumStates() + nmat->GetNumStates();
}

FullSquareMatrix
FelippaShell::stiffness(const CoordSet &cs, double *d, int flg) const
{
  if(prop == NULL) {
    FullSquareMatrix ret(18,d);
    ret.zero();
    return ret;
  }

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);

  double x[3], y[3], z[3];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  double disp[18]; for(int i=0; i<18; ++i) disp[i] = 0;
  double *fint = NULL;

  Impl::andesstf(glNum+1, d, fint, prop->nu, x, y, z, disp, type, gpmat, flg, tflg);

  FullSquareMatrix ret(18,d);

  return ret;
}

Corotator*
FelippaShell::getCorotator(CoordSet& cs, double* kel, int _fitAlg, int)
{
  // Set Node Numbers
  Shell3Corotator::n1 = nn[0];
  Shell3Corotator::n2 = nn[1];
  Shell3Corotator::n3 = nn[2];

  if(type != 4) {

    // for linear material precompute stiffness matrix
    FullSquareMatrix myStiff = stiffness(cs, kel, 0);

    // Copy Element's stiffness matrix to origK (original stiffness matrix)
    for(int i = 0; i < 18; ++i)
      for(int j = 0; j < 18; ++j)
        Shell3Corotator::origK[i][j] = myStiff[i][j];
  }

  // Set fit algorithm
  Shell3Corotator::fitAlg = _fitAlg;

  return this;
}

extern "C" {

  void _FORTRAN(trirotation)(double *, double *);

}

void
FelippaShell::getStiffAndForce(GeomState *refState, GeomState &geomState, CoordSet &cs,
                               FullSquareMatrix &elK, double *f, double dt, double t)
{
  if(prop == NULL) {
    elK.zero();
    for(int i=0; i<18; ++i) f[i] = 0;
    return;
  }

  // Get Nodes original coordinates (C0 configuration)
  Node &node1 = cs.getNode( n1 );
  Node &node2 = cs.getNode( n2 );
  Node &node3 = cs.getNode( n3 );

  // Get Nodes current coordinates (C0n configuration)
  NodeState &ns1 = geomState[ n1 ];
  NodeState &ns2 = geomState[ n2 ];
  NodeState &ns3 = geomState[ n3 ];

  double xl0[3][3], xln[3][3], t0[3][3], t0n[3][3], vld[18], locF[18];

  // C0    = initial configuration
  // C0n   = nth configuration
  // xl0   = C0 local coordinates
  // xln   = C0n local coordinates
  // t0    = transformation matrix to C0
  // t0n   = transformation matrix to C0n
  // vld   = local deformation vector
  // locF  = local unprojected internal force
  // origK = original stiffness matrix

  // Extract deformational displacement from C0 to C0n configurations

  extractDefDisp(node1,node2,node3, ns1,ns2,ns3, xl0,xln, t0,t0n, vld );

  int i, j;

  // Form unprojected internal forces and initialize stiffness matrix

  if(type != 4) {
    for(i=0; i<18; ++i) {
      locF[i] = 0.0;
      for(j=0; j<18; ++j) elK[i][j] = origK[i][j];
    }

    // compute locF (local Force) as origK*vld
    _FORTRAN(dgemv)('N',18,18,1.0,(double *)origK,18,vld,1,0.0,locF,1);
  }
  else { // nonlinear material...

    double x[3], y[3], z[3];

    x[0] = node1.x; y[0] = node1.y; z[0] = node1.z;
    x[1] = node2.x; y[1] = node2.y; z[1] = node2.z;
    x[2] = node3.x; y[2] = node3.y; z[2] = node3.z;

    double *staten = refState->getElemState(getGlNum()) + subNum*numStates();
    double *statenp = geomState.getElemState(getGlNum()) + subNum*numStates();

    // get the nodal temperatures
    Vector ndTemps(3);
    int nodeNum[3] = { n1, n2, n3 };
    geomState.get_temperature(3, nodeNum, ndTemps, prop->Ta);

    Impl::andesstf(glNum+1, elK.data(), locF, prop->nu, x, y, z, vld, type, gpmat, 0,
                   tflg, ndTemps.data(), dt, staten, statenp);

    if(numStates() > 0) {
      double *state = geomState.getElemState(getGlNum()) + subNum*numStates();
      gpmat->GetState( state+0 );
      nmat->GetState ( state+gpmat->GetNumStates() );
    }
  }

  // Compute gradients of the nodal deformational pseudorotations
  // Correct element stiffness and internal force

  double rotvar[3][3][3];

  int inode;
  for(inode=0; inode<3; ++inode)
    pseudorot_var( vld+inode*6+3, rotvar[inode] );

  leftmult_rotvar( 3, 1, rotvar, elK);
  rightmult_rotvar( 3, 0, rotvar, elK);

  double fe[18];

  for(inode=0; inode<3; ++inode)
    for(i=0; i<3; ++i) {
      fe[6*inode+i]   = locF[6*inode+i];
      fe[6*inode+i+3] = rotvar[inode][0][i]*locF[6*inode+3] +
                        rotvar[inode][1][i]*locF[6*inode+4] +
                        rotvar[inode][2][i]*locF[6*inode+5];
    }

  // make copy of rotvar for later
  double rotvarCopy[3][3][3];
  for (inode=0;inode<3;inode++)
    for (j=0;j<3;j++)
      for (i=0;i<3;i++)
        rotvarCopy[inode][j][i] = rotvar[inode][j][i];

  // Add second variation of pseudorotations contracted with the
  // nodal moment to the diagonal blocks of the stiffness

  for(inode = 0; inode<3; ++inode) {
    pseudorot_2var(vld+inode*6+3, locF+inode*6+3, rotvar[inode]);
    for(i=0; i<3; ++i)
      for(j=0; j<3; ++j)
        elK[i+inode*6+3][j+inode*6+3] += rotvar[inode][i][j];
  }

  // Compute nonlinear projector matrix relative to deformed element
  // and correct stiffness and force

  double pmat[18][18], gmat[3][18];

  gradDefDisp(xl0, xln, pmat, gmat);

  double scrstiff[18][18];

  // Form: [K*] = [P'][K]

  _FORTRAN(dgemm)('N','T',18,18,18,1.0,elK.data(),18,
                  (double*)pmat,18,0.0,(double*)scrstiff,18);

  // Form: [K] = [K*][P]

  _FORTRAN(dgemm)('N','N',18,18,18,1.0,(double*)pmat,18,
                  (double*)scrstiff,18,0.0,elK.data(),18);

  // Form: {f} = [P']{fe}

  _FORTRAN(dgemv)('N',18,18,1.0,(double *)pmat,18,fe,1,0.0,f,1);

  // Form geometric stiffness from internal force and material stiffness

  double stiffGeo1[18][18], stiffGeo2[18][18];

  // Build correct geometric stiffness matrix
  // For a zero deformation, stiffGeo1 and stiffGeo2 are zero matrics.

  formCorrectGeometricStiffness(rotvarCopy,xln,pmat,gmat,f,
                                stiffGeo1,stiffGeo2,fe,t0n);

  // Assemble element tangent Stiffness Matrix K = K + Kg1 + Kg2 

  for(i=0; i<18; ++i)
    for(j=0; j<18; ++j)
      elK[i][j] += stiffGeo1[i][j] + stiffGeo2[i][j];

  _FORTRAN(trirotation)( elK.data(), (double*)t0n );

  // Transform internal force vector from local to global coordinates

  tran_force(f, t0n, 3);

  // The skew symmetric load stiffness matrix due to axial external moments is
  // added separately (see Domain::getFollowerForce in Driver.d/NLStatic.C)
  elK.symmetrize();
}

void
FelippaShell::getInternalForce(GeomState *refState, GeomState &geomState, CoordSet &cs,
                               FullSquareMatrix &, double *f, double dt, double t)
{
  if(prop == NULL) {
    for(int i=0; i<18; ++i) f[i] = 0;
    return;
  }

  // Get Nodes original coordinates (C0 configuration)
  Node &node1 = cs.getNode( n1 );
  Node &node2 = cs.getNode( n2 );
  Node &node3 = cs.getNode( n3 );

  // Get Nodes current coordinates (C0n configuration)
  NodeState &ns1 = geomState[ n1 ];
  NodeState &ns2 = geomState[ n2 ];
  NodeState &ns3 = geomState[ n3 ];

  double xl0[3][3], xln[3][3], t0[3][3], t0n[3][3], vld[18], locF[18];

  // C0    = initial configuration
  // C0n   = nth configuration
  // xl0   = C0 local coordinates
  // xln   = C0n local coordinates
  // t0n   = transformation matrix between C0n and C0
  // vld   = local deformation vector
  // locF  = local unprojected internal force
  // origK = original stiffness matrix

  // f = T' P' H' K v

  // Extract deformational displacement from C0 to C0n configurations

  extractDefDisp(node1,node2,node3, ns1,ns2,ns3, xl0,xln, t0,t0n, vld );

  // Form unprojected internal forces and initialize stiffness matrix

  if(type != 4) {

    // compute locF (local Force) as origK*vld

    _FORTRAN(dgemv)('N',18,18,1.0,(double *)origK,18,vld,1,0.0,locF,1);

  }
  else { // nonlinear material...

    double x[3], y[3], z[3];

    x[0] = node1.x; y[0] = node1.y; z[0] = node1.z;
    x[1] = node2.x; y[1] = node2.y; z[1] = node2.z;
    x[2] = node3.x; y[2] = node3.y; z[2] = node3.z;

    double *staten = refState->getElemState(getGlNum()) + subNum*numStates();
    double *statenp = geomState.getElemState(getGlNum()) + subNum*numStates();

    // get the nodal temperatures
    Vector ndTemps(3);
    int nodeNum[3] = { n1, n2, n3 };
    geomState.get_temperature(3, nodeNum, ndTemps, prop->Ta);

    Impl::andesstf(glNum+1, (double*)NULL, locF, prop->nu, x, y, z, vld, type, gpmat, 0,
                   tflg, ndTemps.data(), dt, staten, statenp);

    if(numStates() > 0) {
      double *state = geomState.getElemState(getGlNum()) + subNum*numStates();
      gpmat->GetState( state+0 );
      nmat->GetState ( state+gpmat->GetNumStates() ); 
    }
  }

  // Compute gradients of the nodal deformational pseudorotations
  // Correct element stiffness and internal force

  double rotvar[3][3][3];

  int inode,i;
  for(inode=0; inode<3; ++inode)
    pseudorot_var( vld+inode*6+3, rotvar[inode] );

  double fe[18];

  for(inode=0; inode<3; ++inode)
    for(i=0; i<3; ++i) {
      fe[6*inode+i]   = locF[6*inode+i];
      fe[6*inode+i+3] = rotvar[inode][0][i]*locF[6*inode+3] +
                        rotvar[inode][1][i]*locF[6*inode+4] +
                        rotvar[inode][2][i]*locF[6*inode+5];
    }

  // Compute nonlinear projector matrix relative to deformed element
  // and correct stiffness and force

  double pmat[18][18], gmat[3][18];

  gradDefDisp(xl0, xln, pmat, gmat);

  // Form: {f} = [P']{fe}

  _FORTRAN(dgemv)('N',18,18,1.0,(double *)pmat,18,fe,1,0.0,f,1);

  // Transform internal force vector from local to global coordinates

  tran_force(f, t0n, 3);
}

void
FelippaShell::updateStates(GeomState *refState, GeomState &geomState, CoordSet &cs, double dt)
{
  if(numStates() > 0) {

    // Get Nodes original coordinates (C0 configuration)
    Node &node1 = cs.getNode( n1 );
    Node &node2 = cs.getNode( n2 );
    Node &node3 = cs.getNode( n3 );

    // Get Nodes current coordinates (C0n configuration)
    NodeState &ns1 = geomState[ n1 ];
    NodeState &ns2 = geomState[ n2 ];
    NodeState &ns3 = geomState[ n3 ];

    double xl0[3][3], xln[3][3], t0[3][3], t0n[3][3], vld[18];

    // Extract deformational displacement from C0 to C0n configurations
    extractDefDisp(node1,node2,node3, ns1,ns2,ns3, xl0,xln, t0,t0n, vld);

    double *staten;
    if(refState) {
      staten = refState->getElemState(getGlNum()) + subNum*numStates();
    }
    else {
      staten = new double[numStates()];
      for(int i = 0; i < numStates(); ++i) staten[i] = 0.0;
    }

    gpmat->SetState( staten+0 );
    nmat->SetState ( staten+gpmat->GetNumStates() );

    double x[3], y[3], z[3];

    x[0] = node1.x; y[0] = node1.y; z[0] = node1.z;
    x[1] = node2.x; y[1] = node2.y; z[1] = node2.z;
    x[2] = node3.x; y[2] = node3.y; z[2] = node3.z;

    double *statenp = geomState.getElemState(getGlNum()) + subNum*numStates();

    // get the nodal temperatures
    Vector ndTemps(3);
    int nodeNum[3] = { n1, n2, n3 };
    geomState.get_temperature(3, nodeNum, ndTemps, prop->Ta);

    Impl::andesups(glNum+1, staten, statenp, x, y, z, vld, gpmat, nmat, sflg,
                   tflg, ndTemps.data(), dt);

    if(!refState) delete [] staten;
  }
}

bool
FelippaShell::checkElementDeletion(GeomState &geomState)
{
  double *state = geomState.getElemState(getGlNum()) + subNum*numStates();
  return (prop && type == 4 && gpmat->CheckFailure(state));
}

void
FelippaShell::initStates(double *staten)
{
  if(numStates() > 0) {
    gpmat->SetState( staten+0 );
    nmat->SetState ( staten+gpmat->GetNumStates() );
  }
}

double
FelippaShell::getDissipatedEnergy(GeomState &curState, CoordSet &cs)
{
  if(type != 4) {

    return 0.0;
  }
  else { // nonlinear material...
    
    // Get Nodes original coordinates (C0 configuration)
    Node &node1 = cs.getNode( n1 );
    Node &node2 = cs.getNode( n2 );
    Node &node3 = cs.getNode( n3 );

    double x[3], y[3], z[3];
    x[0] = node1.x; y[0] = node1.y; z[0] = node1.z;
    x[1] = node2.x; y[1] = node2.y; z[1] = node2.z;
    x[2] = node3.x; y[2] = node3.y; z[2] = node3.z;

    // AN: following part updated. Now send state to andesden
    double *statenp = curState.getElemState(getGlNum()) + subNum*numStates();;

    double D;
    Impl::andesden(glNum+1, x, y, z, gpmat, statenp, D);

    return D;
  }
}

int
FelippaShell::numNodes() const
{
  return 3;
}

int*
FelippaShell::nodes(int *p) const
{
  if(p == 0) p = new int[3];
  p[0] = nn[0];
  p[1] = nn[1];
  p[2] = nn[2];
  return p;
}

int
FelippaShell::numDofs() const
{
  return 18;
}

int*
FelippaShell::dofs(DofSetArray &dsa, int *p) const
{
  if(p == 0) p = new int[18];

  dsa.number(nn[0],DofSet::XYZdisp | DofSet::XYZrot, p  );
  dsa.number(nn[1],DofSet::XYZdisp | DofSet::XYZrot, p+6);
  dsa.number(nn[2],DofSet::XYZdisp | DofSet::XYZrot, p+12);

  return p;
}

void
FelippaShell::markDofs(DofSetArray &dsa) const
{
  dsa.mark(nn, 3, DofSet::XYZdisp | DofSet::XYZrot);
}

void
FelippaShell::computeDisp(CoordSet&, State &state, const InterpPoint &ip,
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
FelippaShell::getFlLoad(CoordSet &, const InterpPoint &ip, double *flF, 
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

int
FelippaShell::getTopNumber() const
{
  return 108;
}

void
FelippaShell::computePressureForce(CoordSet& cs, Vector& elPressureForce,
                                   GeomState *geomState, int cflg, double time)
{
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

        double lmy = -pressureForce*length/8.0;
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

    // compute centroid
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

    double fmf = 8; // or 12?

    elPressureForce[0]  = 0;
    elPressureForce[1]  = 0;
    elPressureForce[2]  = pressureForce;
    elPressureForce[3]  = cflg*pressureForce*( xl0[2][1] - xl0[0][1]
                                              +xl0[1][1] - xl0[0][1])/fmf;
    elPressureForce[4]  = cflg*pressureForce*( xl0[0][0] - xl0[2][0]
                                              +xl0[0][0] - xl0[1][0])/fmf;
    elPressureForce[5]  = 0;

    elPressureForce[6]  = 0;
    elPressureForce[7]  = 0;
    elPressureForce[8]  = pressureForce;
    elPressureForce[9]  = cflg*pressureForce*( xl0[0][1] - xl0[1][1]
                                              +xl0[2][1] - xl0[1][1])/fmf;
    elPressureForce[10] = cflg*pressureForce*( xl0[1][0] - xl0[0][0]
                                              +xl0[1][0] - xl0[2][0])/fmf;
    elPressureForce[11] = 0;

    elPressureForce[12] = 0;
    elPressureForce[13] = 0;
    elPressureForce[14] = pressureForce;
    elPressureForce[15] = cflg*pressureForce*( xl0[1][1] - xl0[2][1]
                                              +xl0[0][1] - xl0[2][1])/fmf;
    elPressureForce[16] = cflg*pressureForce*( xl0[2][0] - xl0[1][0]
                                              +xl0[2][0] - xl0[0][0])/fmf;
    elPressureForce[17] = 0;

  }
}

void
FelippaShell::getThermalForce(CoordSet& cs, Vector& ndTemps, Vector &elThermalForce, int glflag, GeomState *geomState)
{
  // force is returned in global coordinates if glflag = 0
  // for nonlinear analyses with nonlinear material law, the thermal load for this element is now computed in getStiffAndForce
  if(prop == NULL || (geomState && type == 4)) {
    elThermalForce.zero();
    return;
  }

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);

  double x[3], y[3], z[3];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  double disp[18]; for(int i=0; i<18; ++i) disp[i] = 0;

  int flg = (glflag == 0) ? 1 : 0;
  Impl::andesstf(glNum+1, (double *)NULL, elThermalForce.data(), prop->nu,
                 x, y, z, disp, type, gpmat, flg, tflg, ndTemps.data());
  elThermalForce *= -1;
}

double
FelippaShell::getMassThicknessSensitivity(CoordSet &cs)
{ 
  if(prop == NULL) return 0.0;
  if(type == 1 || type == 2 || type == 3) {
     throw std::runtime_error("\n"
          "*** FATAL ERROR in FelippaShell::getMassThicknessSensitivity ***\n"
          "*** Selected composite constitutive law is not               ***\n"
          "*** suitable for thickness sensitivity analysis              ***\n"
          "*** STOP ALL TREATMENTS RIGHT HERE                           ***\n");
  }

  double x[3] = { cs[nn[0]]->x, cs[nn[1]]->x, cs[nn[2]]->x };
  double y[3] = { cs[nn[0]]->y, cs[nn[1]]->y, cs[nn[2]]->y };
  double z[3] = { cs[nn[0]]->z, cs[nn[1]]->z, cs[nn[2]]->z };
  double *gravityAcceleration = NULL, *grvfor = NULL;
  double totmas = 0.0;
  double rhoh = gpmat->GetAreaDensity();

  // note: For types other than 2 or 3, the contribution of the "non-structural" mass (which is presumably independent of the thickness)
  //       to the sensitivity cannot be accounted for correctly because it is combined with the "structural" mass in the input file.
  //       This can be fixed by adding to the input file an option to define the proportion of the total mass which is non-structural. 
  double nsm = (type == 2 || type == 3) ? prop->rho : 0;

  Impl::andesms(glNum+1, x, y, z, gravityAcceleration, grvfor, totmas, rhoh-nsm);

  return totmas/gpmat->GetShellThickness();
}

void
FelippaShell::getWeightNodalCoordinateSensitivity(Vector &dwdx, CoordSet& cs, double *gravityAcceleration)
{
  if(prop == NULL || gravityAcceleration == NULL) {
    dwdx.zeroAll();
    return;
  }

  double x[3] = { cs[nn[0]]->x, cs[nn[1]]->x, cs[nn[2]]->x };
  double y[3] = { cs[nn[0]]->y, cs[nn[1]]->y, cs[nn[2]]->y };
  double z[3] = { cs[nn[0]]->z, cs[nn[1]]->z, cs[nn[2]]->z };

  double rhoh = gpmat->GetAreaDensity();

  Impl::andesmsWRTcoord(glNum+1, x, y, z, dwdx.data(), rhoh);

  double gravAccNorm = sqrt(gravityAcceleration[0]*gravityAcceleration[0] +
                            gravityAcceleration[1]*gravityAcceleration[1] +
                            gravityAcceleration[2]*gravityAcceleration[2]);

  dwdx *= gravAccNorm;
}

void
FelippaShell::getGravityForceThicknessSensitivity(CoordSet& cs, double *gravityAcceleration, 
                                                  Vector& dGfdthick, int gravflg, GeomState*)
{
  if(prop == NULL) {
    dGfdthick.zero();
    return;
  }
  if(type == 1 || type == 2 || type == 3) {
     throw std::runtime_error("\n"
          "*** FATAL ERROR in FelippaShell::getGravityForceThicknessSensitivity ***\n"
          "*** Selected composite constitutive law is not                       ***\n"
          "*** suitable for thickness sensitivity analysis                      ***\n"
          "*** STOP ALL TREATMENTS RIGHT HERE                                   ***\n");
  }

  double x[3] = { cs[nn[0]]->x, cs[nn[1]]->x, cs[nn[2]]->x };
  double y[3] = { cs[nn[0]]->y, cs[nn[1]]->y, cs[nn[2]]->y };
  double z[3] = { cs[nn[0]]->z, cs[nn[1]]->z, cs[nn[2]]->z };

  double rhoh = gpmat->GetAreaDensity();

  // note: For types other than 2 or 3, the contribution of the "non-structural" mass (which is presumably independent of the thickness)
  //       to the sensitivity cannot be accounted for correctly because it is combined with the "structural" mass in the input file.
  //       This can be fixed by adding to the input file an option to define the proportion of the total mass which is non-structural. 
  double nsm = (type == 2 || type == 3) ? prop->rho : 0;

  Impl::andesgf(glNum+1, x, y, z, dGfdthick.data(),
                gravityAcceleration, gravflg, rhoh-nsm);

  dGfdthick /= gpmat->GetShellThickness();
}

void
FelippaShell::getGravityForceNodalCoordinateSensitivity(CoordSet& cs, double *gravityAcceleration,
                                                        GenFullM<double> &dGfdx, int gravflg, GeomState*)
{
  if(prop == NULL) {
    dGfdx.zero();
    return;
  }

  double x[3] = { cs[nn[0]]->x, cs[nn[1]]->x, cs[nn[2]]->x };
  double y[3] = { cs[nn[0]]->y, cs[nn[1]]->y, cs[nn[2]]->y };
  double z[3] = { cs[nn[0]]->z, cs[nn[1]]->z, cs[nn[2]]->z };

  double rhoh = gpmat->GetAreaDensity();

  Impl::andesgfWRTcoord(glNum+1, x, y, z, dGfdx.data(),
                        gravityAcceleration, gravflg, rhoh);
}

void 
FelippaShell::getStiffnessThicknessSensitivity(CoordSet &cs, FullSquareMatrix &dStiffdThick, int flg)
{
  if(prop == NULL) {
    dStiffdThick.zero();
    return;
  }
  if(type == 1 || type == 2 || type == 3) {
     throw std::runtime_error("\n"
          "*** FATAL ERROR in FelippaShell::getStiffnessThicknessSensitivity ***\n"
          "*** Selected composite constitutive law is not                    ***\n"
          "*** suitable for thickness sensitivity analysis                   ***\n"
          "*** STOP ALL TREATMENTS RIGHT HERE                                ***\n");
  }

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);

  double x[3], y[3], z[3];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  double disp[18]; for(int i=0; i<18; ++i) disp[i] = 0;
  Impl::andesstfWRTthic(glNum+1, dStiffdThick.data(), (double*)NULL, prop->nu,
                        x, y, z, disp, type, gpmat, flg, tflg);
}

void 
FelippaShell::getStiffnessNodalCoordinateSensitivity(FullSquareMatrix *&dStiffdx, CoordSet &cs)
{
  if(prop == NULL) {
    for(int i=0; i<9; ++i) dStiffdx[i].zero();
    return;
  }

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);

  double x[3], y[3], z[3];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  double *data[9];
  for(int i=0; i<9; ++i) data[i] = dStiffdx[i].data();

  int flg = 1;

  Impl::andesstfWRTcoord(glNum+1, data, prop->E, prop->nu,
                         prop->rho, prop->eh, prop->Ta, prop->W,
                         cFrame, x, y, z, type,
                         gpmat->GetCoefOfConstitutiveLaw(), flg, tflg);
}

void 
FelippaShell::getVonMisesThicknessSensitivity(Vector &dStdThick, Vector &weight, CoordSet &cs,
                                              Vector &elDisp, int, int surface, double *ndTemps,
                                              int avgnum, double, double)
{
  if(type == 1 || type == 2 || type == 3) {
     throw std::runtime_error("\n"
          "*** FATAL ERROR in FelippaShell::getVonMisesThicknessSensitivity ***\n"
          "*** Selected composite constitutive law is not                   ***\n"
          "*** suitable for thickness sensitivity analysis                  ***\n"
          "*** STOP ALL TREATMENTS RIGHT HERE                               ***\n");
  }

  weight = 1.0;

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);

  double x[3], y[3], z[3];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  Impl::andesvmsWRTthic(glNum+1, prop->nu, x, y, z, elDisp.data(),
                        dStdThick.getData(), type, nmat, surface,
                        sflg, ndTemps);
}

void 
FelippaShell::getVonMisesNodalCoordinateSensitivity(GenFullM<double> &dStdx, Vector &weight, CoordSet &cs,
                                                    Vector &elDisp, int strInd, int surface, double *ndTemps,
                                                    int avgnum, double, double)
{
  weight = 1.0;

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);

  double x[3], y[3], z[3];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  Impl::andesvmsWRTcoord(glNum+1, prop->E, prop->nu, prop->rho,
                         prop->eh, prop->Ta, prop->W, cFrame,
                         x, y, z, elDisp.data(), dStdx.getData(),
                         type, nmat->GetCoefOfConstitutiveLaw(),
                         surface, sflg, ndTemps);
}

void 
FelippaShell::getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, GenFullM<double> *dDispDisp,
                                                 CoordSet &cs, Vector &elDisp, int strInd, int surface, double *ndTemps,
                                                 int avgnum, double, double)
{
  weight = 1.0;

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);

  double x[3], y[3], z[3];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  double* disp = elDisp.data();

  Impl::andesvmsWRTdisp(glNum+1, prop->nu, x, y, z, elDisp.data(),
                        dStdDisp.getData(), type, nmat, surface,
                        sflg, ndTemps);  
  if(dDispDisp) dStdDisp ^= (*dDispDisp); 
 
}

void
FelippaShell::getInternalForceThicknessSensitivity(GeomState *refState, GeomState &geomState, CoordSet &cs, Vector &dFintdThick,
                                                   double dt, double t)
{
  if(prop == NULL) {
    for(int i=0; i<18; ++i) dFintdThick[i] = 0;
    return;
  }
  if(type == 1 || type == 2 || type == 3) {
     throw std::runtime_error("\n"
          "*** FATAL ERROR in FelippaShell::getInternalForceThicknessSensitivity ***\n"
          "*** Selected composite constitutive law is not                        ***\n"
          "*** suitable for thickness sensitivity analysis                       ***\n"
          "*** STOP ALL TREATMENTS RIGHT HERE                                    ***\n");
  }

  // Get Nodes original coordinates (C0 configuration)
  Node &node1 = cs.getNode( n1 );
  Node &node2 = cs.getNode( n2 );
  Node &node3 = cs.getNode( n3 );

  // Get Nodes current coordinates (C0n configuration)
  NodeState &ns1 = geomState[ n1 ];
  NodeState &ns2 = geomState[ n2 ];
  NodeState &ns3 = geomState[ n3 ];

  double xl0[3][3], xln[3][3], t0[3][3], t0n[3][3], vld[18], locF[18];

  // C0    = initial configuration
  // C0n   = nth configuration
  // xl0   = C0 local coordinates
  // xln   = C0n local coordinates
  // t0n   = transformation matrix between C0n and C0
  // vld   = local deformation vector
  // locF  = local unprojected internal force
  // origK = original stiffness matrix

  // f = T' P' H' K v

  // Extract deformational displacement from C0 to C0n configurations

  extractDefDisp(node1,node2,node3, ns1,ns2,ns3, xl0,xln, t0,t0n, vld);

  // Form unprojected internal forces and initialize stiffness matrix

  double x[3], y[3], z[3];

  x[0] = node1.x; y[0] = node1.y; z[0] = node1.z;
  x[1] = node2.x; y[1] = node2.y; z[1] = node2.z;
  x[2] = node3.x; y[2] = node3.y; z[2] = node3.z;

  Impl::andesstfWRTthic(glNum+1, (double*)NULL, locF, prop->nu, x, y, z, vld, type, gpmat, 0, tflg);

  // Compute gradients of the nodal deformational pseudorotations
  // Correct element stiffness and internal force

  double rotvar[3][3][3];

  int inode,i;
  for(inode=0; inode<3; ++inode)
    pseudorot_var( vld+inode*6+3, rotvar[inode] );

  double fe[18];

  for(inode=0; inode<3; ++inode)
    for(i=0; i<3; ++i) {
      fe[6*inode+i]   = locF[6*inode+i];
      fe[6*inode+i+3] = rotvar[inode][0][i]*locF[6*inode+3] +
                        rotvar[inode][1][i]*locF[6*inode+4] +
                        rotvar[inode][2][i]*locF[6*inode+5];
    }

  // Compute nonlinear projector matrix relative to deformed element
  // and correct stiffness and force

  double pmat[18][18], gmat[3][18];

  gradDefDisp(xl0, xln, pmat, gmat);

  // Form: {f} = [P']{fe}

  _FORTRAN(dgemv)('N',18,18,1.0,(double *)pmat,18,fe,1,0.0,dFintdThick.data(),1);

  // Transform internal force vector from local to global coordinates

  tran_force(dFintdThick.data(), t0n, 3);
}

#endif
