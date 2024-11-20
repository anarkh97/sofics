#ifdef USE_EIGEN3
#include <Element.d/Rigid.d/RigidFourNodeShell.h>
#include <Element.d/Rigid.d/RigidBeam.h>
#include <Corotational.d/GeomState.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/Vector.h>
#include <Utils.d/dbg_alloca.h>
#include <Utils.d/linkfc.h>

extern "C" {
  void _FORTRAN(elemaslbt)(int&, double*, double*, double*, double*);
  void _FORTRAN(elefbc3dbrkshl2)(int&, double*, double*, double*, double*);
}


RigidFourNodeShell::RigidFourNodeShell(int *_nn)
 : SuperElement(true)
{
  nnodes = 4;
  nn = new int[nnodes];
  for(int i = 0; i < nnodes; ++i) nn[i] = _nn[i];
  nSubElems = 3;
  subElems = new Element * [nSubElems];
  for(int i = 0; i < nSubElems; ++i) {
    int indices[2] = { i+1, 0 };
    subElems[i] = new RigidBeam(indices);
  }
  pbc = 0;
}

double
RigidFourNodeShell::getMass(const CoordSet& cs) const
{
  if (prop == NULL || prop->rho == 0 || prop->eh == 0) return 0.0;

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);
  auto &nd4 = cs.getNode(nn[3]);

  Vector r1(3), r2(3), r3(3), r4(3);

  r1[0] = nd1.x; r1[1] = nd1.y; r1[2] = 0.0;
  r2[0] = nd2.x; r2[1] = nd2.y; r2[2] = 0.0;
  r3[0] = nd3.x; r3[1] = nd3.y; r3[2] = 0.0;
  r4[0] = nd4.x; r4[1] = nd4.y; r4[2] = 0.0;

  Vector v1(3), v2(3), v3(3), v4(3), v5(3);

  v1 = r2 - r1;
  v2 = r3 - r1;
  v3 = r4 - r1;

  v4 = v1.cross(v2);
  v5 = v2.cross(v3);

  double area = 0.5*(v4.magnitude() + v5.magnitude());
  double mass = area*prop->rho*prop->eh;

  return mass;
}

void
RigidFourNodeShell::getGravityForce(CoordSet& cs, double *gravityAcceleration,
                                    Vector& gravityForce, int gravflg, GeomState *geomState)
{
  gravityForce.zero();
  if (prop == NULL || prop->rho == 0 || prop->eh == 0) return;

  double massPerNode = 0.25*getMass(cs);
  double fx = massPerNode*gravityAcceleration[0];
  double fy = massPerNode*gravityAcceleration[1];
  double fz = massPerNode*gravityAcceleration[2];

  for(int i = 0; i < 4; ++i) {
    gravityForce[6*i+0] = fx;
    gravityForce[6*i+1] = fy;
    gravityForce[6*i+2] = fz;
  }
}

FullSquareMatrix
RigidFourNodeShell::massMatrix(const CoordSet &cs, double *mel, int cmflg) const
{
  int nndof = 6, ndime = 3;
  FullSquareMatrix ret(numDofs(), mel);
  ret.zero();
  if (prop == NULL || prop->rho == 0 || prop->eh == 0) return ret;

  // Check for element which has no mass
  if(prop && prop->rho != 0 && prop->eh != 0) {
    double* ecord = new double[nnodes*ndime];
    double* edisp = new double[nnodes*nndof];
    for(int i = 0; i < nnodes; ++i) {
      ecord[i*ndime+0] = cs[nn[i]]->x;
      ecord[i*ndime+1] = cs[nn[i]]->y;
      ecord[i*ndime+2] = cs[nn[i]]->z;
      for(int j = 0; j < nndof; ++j) edisp[i*nndof+j] = 0; // constant mass matrix
    }

    double* emasl = (double*) dbg_alloca(sizeof(double)*nnodes*nndof);
    double ematpro[20];
    ematpro[2] = prop->rho; ematpro[19] = prop->eh;
    // get bt shell element lumped mass
    _FORTRAN(elemaslbt)(nndof, ematpro, ecord, edisp, emasl);
       // input : nndof,ematpro,ecord,edisp
       // output : emasl
    delete [] ecord;
    delete [] edisp;

    for(int i = 0; i < nnodes*nndof; ++i) ret[i][i] = emasl[i];
  }

  return ret;
}

void
RigidFourNodeShell::computePressureForce(CoordSet& cs, Vector& elPressureForce,
                                         GeomState *geomState, int cflg, double time)
{
  if(!pbc) { elPressureForce = 0.; return; }
  int opttrc = 0; // 0 : pressure
                  // 1 : traction
  int optele = 3, ndime = 3;
  double* ecord = (double*) dbg_alloca(sizeof(double)*nnodes*ndime);
  double* edisp = (double*) dbg_alloca(sizeof(double)*nnodes*ndime); // translations only
  int iloc;
  for(int i = 0; i < nnodes; ++i) {
    iloc = i*ndime;
    ecord[iloc+0] = cs[nn[i]]->x;
    ecord[iloc+1] = cs[nn[i]]->y;
    ecord[iloc+2] = cs[nn[i]]->z;
    edisp[iloc+0] = (geomState) ? (*geomState)[nn[i]].x - cs[nn[i]]->x : 0;
    edisp[iloc+1] = (geomState) ? (*geomState)[nn[i]].y - cs[nn[i]]->y : 0;
    edisp[iloc+2] = (geomState) ? (*geomState)[nn[i]].z - cs[nn[i]]->z : 0;
  }
  double pressure = pbc->val;
  // Check if Conwep is being used. If so, add the pressure from the blast loading function.
  if(pbc->conwep && pbc->conwepswitch) {
    pressure += BlastLoading::ComputeShellPressureLoad(ecord, time, *(pbc->conwep));
  }
  double trac[3] = { -pressure, 0, 0 };
  double *efbc = (double*) dbg_alloca(sizeof(double)*nnodes*ndime); // translations only

  _FORTRAN(elefbc3dbrkshl2)(opttrc, ecord, edisp, trac, efbc);

  for(int i = 0; i < nnodes; ++i)
    for(int j = 0; j < ndime; ++j)
      elPressureForce[6*i+j] = efbc[3*i+j];

  for(int i=nnodes*ndime; i<numDofs(); ++i) elPressureForce[i] = 0; // lagrange multiplier dofs, if any
}

#include <Element.d/State.h>
#include <Hetero.d/InterpPoint.h>

void
RigidFourNodeShell::computeDisp(CoordSet& cs, State& state, const InterpPoint& ip,
                                double *res, GeomState *gs)
{
  const double *gp = ip.xy;
  double xyz[4][6];
  state.getDV(nn[0], xyz[0], xyz[0]+3);
  state.getDV(nn[1], xyz[1], xyz[1]+3);
  state.getDV(nn[2], xyz[2], xyz[2]+3);
  state.getDV(nn[3], xyz[3], xyz[3]+3);

  int j;
  for(j=0; j<6; ++j)
    res[j] = (1-gp[0])*(1-gp[1])* xyz[0][j] +
             gp[0]*(1-gp[1])* xyz[1][j] +
             (1-gp[0])*gp[1]* xyz[3][j] +
             gp[0]*gp[1]*xyz[2][j];
}

void
RigidFourNodeShell::getFlLoad(CoordSet& cs, const InterpPoint& ip, double *flF,
                              double *resF, GeomState *gs)
{
  // PJSA 9/10/2010 reversed resF[12+i] and resF[18+i] to match 
  // FlExchanger::getQuadFlLoad in Xfem/Hetero.d/FlExchange.C
  const double *gp = ip.xy;
  for(int i = 0; i < 3; ++i) {
    resF[i]    = (1-gp[0])*(1-gp[1])* flF[i];
    resF[6+i]  = gp[0]*(1-gp[1])* flF[i];
    resF[12+i] = gp[0]*gp[1]* flF[i];
    resF[18+i] = (1-gp[0])*gp[1]* flF[i];
    resF[i+3]  = resF[i+9] = resF[i+15] = resF[i+21] = 0.0;
  }
}

#endif
