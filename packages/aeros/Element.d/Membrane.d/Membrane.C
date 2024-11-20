#include <cstdio>
#include <cmath>

#include <Element.d/Membrane.d/Membrane.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/pstress.h>
#include <Corotational.d/utilities.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/Vector.h>
#include <Math.d/matrix.h>
#include <Corotational.d/Shell3Corotator.h>

extern int verboseFlag;
extern "C" {
void _FORTRAN(trimem)(int&, double*, double*, double*, double&, double&,
                      double*, double*);
void _FORTRAN(sands19)(double*, double*, double*, double&, double&, double*, double*,
                       double*, const int&, const int&, const int&, const int&, const int&);
}

Membrane::Membrane(int* nodeNumbers)
{
        nn[0] = nodeNumbers[0];
        nn[1] = nodeNumbers[1];
        nn[2] = nodeNumbers[2];
}

Element *
Membrane::clone()
{
        return new Membrane(*this);
}

void
Membrane::renum(const int *table)
{
        nn[0] = table[nn[0]];
        nn[1] = table[nn[1]];
        nn[2] = table[nn[2]];
}

void
Membrane::renum(EleRenumMap& table)
{
        nn[0] = table[nn[0]];
        nn[1] = table[nn[1]];
        nn[2] = table[nn[2]];
}

void
Membrane::getVonMises(Vector& stress, Vector& weight, CoordSet &cs,
                      Vector& elDisp, int strInd, int surface,
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

        // Thickness
        h[0] = h[1] = h[2] = prop->eh;

        double elStress[3][7];

        int  msize = 1;
        int maxstr = 7;
        int maxgus = 3;
        int    elm = 1;

        // SET STRAIN FLAG IF USER WANTS STRAIN OUTPUT 
        int strainFlg;
        if(strInd > 6) strainFlg = 1; 

       _FORTRAN(sands19)(x,y,z,prop->E,prop->nu,h,elDisp.data(),(double*)elStress,
                         msize,maxstr,maxgus,elm,strainFlg);

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
Membrane::getAllStress(FullM& stress, Vector& weight, CoordSet &cs,
                       Vector& elDisp, int strInd, int surface,
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

        // Thickness
        h[0] = h[1] = h[2] = prop->eh;

        double elStress[3][7];

        int  msize = 1;
        int maxstr = 7;
        int maxgus = 3;
        int    elm = 1;

       _FORTRAN(sands19)(x,y,z,prop->E,prop->nu,h,elDisp.data(),(double*)elStress,
                         msize,maxstr,maxgus,elm,strInd);

// Store all Stress or all Strain as defined by strInd
        int i,j;
        for (i=0; i<3; ++i) {
          for (j=0; j<6; ++j) {
            stress[i][j] = elStress[i][j];
          }
        }

// Get Element Principals
        double svec[6], pvec[3];
        for (j=0; j<6; ++j) {
          svec[j] = stress[0][j];
        }
// Convert Engineering to Tensor Strains
        if(strInd != 0) {
          svec[3] /= 2;
          svec[4] /= 2;
          svec[5] /= 2;
        }
        pstress(svec,pvec);
        for (i=0; i<3; ++i) {
          for (j=0; j<3; ++j) {
            stress[i][j+6] = pvec[j];
          }
        }
}

double
Membrane::getMass(const CoordSet& cs) const
{
        if (prop == NULL) return 0.0;

        auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);

        double r1[3], r2[3], r3[3], v1[3], v2[3], v3[3];

        r1[0] = nd1.x; r1[1] = nd1.y; r1[2] = nd1.z;
        r2[0] = nd2.x; r2[1] = nd2.y; r2[2] = nd2.z;
        r3[0] = nd3.x; r3[1] = nd3.y; r3[2] = nd3.z;

        v1[0] = r3[0] - r1[0];
        v1[1] = r3[1] - r1[1];
        v1[2] = r3[2] - r1[2];

        v2[0] = r2[0] - r1[0];
        v2[1] = r2[1] - r1[1];
        v2[2] = r2[2] - r1[2];

        crossprod(v1, v2, v3);

        double area = 0.5*sqrt(v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2]);
        double mass = area*prop->rho*prop->eh;

        return mass;

}

void
Membrane::getGravityForce(CoordSet& cs, double *gravityAcceleration,
                          Vector& gravityForce, int gravflg, GeomState *geomState)
{
        double mass = getMass(cs);
        double massPerNode = mass/3.0;
        double fx, fy, fz;

        // Lumped
        if(gravflg != 2) {

          fx = massPerNode*gravityAcceleration[0];
          fy = massPerNode*gravityAcceleration[1];
          fz = massPerNode*gravityAcceleration[2];

        }
        // Consistent
        else {
          int i;
          auto &nd1 = cs.getNode(nn[0]);
          auto &nd2 = cs.getNode(nn[1]);
          auto &nd3 = cs.getNode(nn[2]);
          double x[3], y[3], z[3], localg[3];
          double T1[3],T2[3],T3[3];

          // Set the coordinates
          x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
          x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
          x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

          // Local X-axis
          T1[0] = x[1] - x[0];
          T1[1] = y[1] - y[0];
          T1[2] = z[1] - z[0];
          normalize( T1 );
          // 2nd Vector In Plane
          T2[0] = x[2] - x[0];
          T2[1] = y[2] - y[0];
          T2[2] = z[2] - z[0];
          normalize( T2 );
          // Local Z-axis as cross product of x-axis and in-plane vector
          crossprod( T1, T2, T3 );
          normalize( T3 );
          // Local Y-axis as cross product of x-axis and z-axis
          crossprod( T3, T1, T2 );
          normalize( T2 );

          for(i=0; i<3; ++i)
            localg[i] = 0.0;

          for(i=0; i<3; ++i) {
            localg[0] += T1[i]*gravityAcceleration[i];
            localg[1] += T2[i]*gravityAcceleration[i];
            localg[2] += T3[i]*gravityAcceleration[i];
          }
          double localf[3];
          localf[0] = massPerNode*localg[0];
          localf[1] = massPerNode*localg[1];
          localf[2] = massPerNode*localg[2];

          fx = (T1[0]*localf[0]) + (T2[0]*localf[1]) + (T3[0]*localf[2]);
          fy = (T1[1]*localf[0]) + (T2[1]*localf[1]) + (T3[1]*localf[2]);
          fz = (T1[2]*localf[0]) + (T2[2]*localf[1]) + (T3[2]*localf[2]);
        }

        gravityForce[0] =  fx;
        gravityForce[1] =  fy;
        gravityForce[2] =  fz;
        gravityForce[3] =  0.0;
        gravityForce[4] =  0.0;
        gravityForce[5] =  0.0;
        gravityForce[6] =  fx;
        gravityForce[7] =  fy;
        gravityForce[8] =  fz;
        gravityForce[9] =  0.0;
        gravityForce[10] =  0.0;
        gravityForce[11] =  0.0;
        gravityForce[12] =  fx;
        gravityForce[13] =  fy;
        gravityForce[14] =  fz;
        gravityForce[15] =  0.0;
        gravityForce[16] =  0.0;
        gravityForce[17] =  0.0;
}

FullSquareMatrix
Membrane::massMatrix(const CoordSet &cs,double *mel,int cmflg) const
{
        auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);

        double x[3], y[3], z[3];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

        double x21   = x[1] - x[0];
        double y21   = y[1] - y[0];
        double z21   = z[1] - z[0];
        double x32   = x[2] - x[1];
        double y32   = y[2] - y[1];
        double z32   = z[2] - z[1];
        double x13   = x[0] - x[2];
        double y13   = y[0] - y[2];
        double z13   = z[0] - z[2];

        double rl[3];
        rl[0] = sqrt(x21*x21 + y21*y21 + z21*z21);
        rl[1] = sqrt(x32*x32 + y32*y32 + z32*z32);
        rl[2] = sqrt(x13*x13 + y13*y13 + z13*z13);

        double rmas = getMass(cs)/3.0; 

        FullSquareMatrix ret(18,mel);

        ret.zero();

        std::cerr << " *** WARNING: Membrane::massMatrix is not implemented\n";
/* XXX the following is not correct, mass matrix should be 18x18 and in global frame
        ret[0][0] = rmas;
        ret[1][1] = rmas;
        ret[2][2] = (rl[0]*rl[0]+rl[2]*rl[2])/420.0*rmas;
        ret[3][3] = rmas;
        ret[4][4] = rmas;
        ret[5][5] = (rl[1]*rl[1]+rl[0]*rl[0])/420.0*rmas;
        ret[6][6] = rmas;
        ret[7][7] = rmas;
        ret[8][8] = (rl[2]*rl[2]+rl[1]*rl[1])/420.0*rmas;
*/
        return ret;
}

FullSquareMatrix
Membrane::stiffness(const CoordSet &cs, double *d, int flg) const
{
        auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);

        double x[3], y[3], z[3], h[3];

        // Set the coordinates
        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

        // Set the thickness
        h[0] = h[1] = h[2] = prop->eh;

        if(h[0] <= 0.0)
          fprintf(stderr,"ERROR: Zero shell thickness (Membrane.C) %d %d %d\n",
                nn[0], nn[1], nn[2]);

        _FORTRAN(trimem)(flg, x, y, z, prop->E, prop->nu, h, (double *)d);

        FullSquareMatrix ret(18,d);

        return ret;
}

int
Membrane::numNodes() const
{
         return 3;
}

int*
Membrane::nodes(int *p) const
{
         if(p == 0) p = new int[3];
         p[0] = nn[0];
         p[1] = nn[1];
         p[2] = nn[2];
        return p;
}

int
Membrane::numDofs() const
{
         return 18;
}

int*
Membrane::dofs(DofSetArray &dsa, int *p) const
{
        if(p == 0) p = new int[18];

        dsa.number(nn[0], DofSet::XYZdisp | DofSet::XYZrot, p  );
        dsa.number(nn[1], DofSet::XYZdisp | DofSet::XYZrot, p+6);
        dsa.number(nn[2], DofSet::XYZdisp | DofSet::XYZrot, p+12);

        return p;
}

void
Membrane::markDofs(DofSetArray &dsa) const
{
        dsa.mark(nn, 3, DofSet::XYZdisp | DofSet::XYZrot);
}

int
Membrane::getTopNumber() const
{
  return 108;
}

Corotator *
Membrane::getCorotator(CoordSet &cs, double *kel, int fitAlg, int)
{
 int flag = 0; // signals stiffness routine to keep local matrix 
 FullSquareMatrix myStiff = stiffness(cs, kel, flag);
 return new Shell3Corotator(nn[0], nn[1], nn[2], myStiff, fitAlg);
}

#ifdef USE_EIGEN3
#include <Element.d/FelippaShell.d/ShellElementTemplate.hpp>
#include <Element.d/FelippaShell.d/EffMembraneTriangle.hpp>
#include <Element.d/FelippaShell.d/NoBendingTriangle.hpp>

typedef ShellElementTemplate<double,EffMembraneTriangle,NoBendingTriangle> Impl;

double
Membrane::getMassThicknessSensitivity(CoordSet& cs)
{
  if(prop == NULL) return 0.0;
  else return getMass(cs)/prop->eh;
}

void
Membrane::getWeightNodalCoordinateSensitivity(Vector &dwdx, CoordSet& cs, double *gravityAcceleration, int)
{
  if(prop == NULL || gravityAcceleration == NULL) {
    dwdx.zeroAll();
    return;
  }

  double x[3] = { cs[nn[0]]->x, cs[nn[1]]->x, cs[nn[2]]->x };
  double y[3] = { cs[nn[0]]->y, cs[nn[1]]->y, cs[nn[2]]->y };
  double z[3] = { cs[nn[0]]->z, cs[nn[1]]->z, cs[nn[2]]->z };

  double rhoh = prop->rho*prop->eh;

  Impl::andesmsWRTcoord(glNum+1, x, y, z, dwdx.data(), rhoh);

  double gravAccNorm = sqrt(gravityAcceleration[0]*gravityAcceleration[0] +
                            gravityAcceleration[1]*gravityAcceleration[1] +
                            gravityAcceleration[2]*gravityAcceleration[2]);

  dwdx *= gravAccNorm;
}

void
Membrane::getGravityForceThicknessSensitivity(CoordSet& cs, double *gravityAcceleration, int,
                                              Vector& dGfdthick, int gravflg, GeomState *geomState)
{
  if(prop == NULL) {
    dGfdthick.zero();
    return;
  }

  getGravityForce(cs, gravityAcceleration, dGfdthick, gravflg, geomState);

  dGfdthick /= prop->eh;
}

void
Membrane::getGravityForceNodalCoordinateSensitivity(CoordSet& cs, double *gravityAcceleration, int,
                                                    GenFullM<double> &dGfdx, int gravflg, GeomState*)
{
  std::cerr << " *** WARNING: Membrane::getGravityForceNodalCoordinateSensitivity is not implemented\n";
  dGfdx.zero();
}

void
Membrane::getStiffnessThicknessSensitivity(CoordSet &cs, FullSquareMatrix &dStiffdThick, int flg, int)
{
  if(prop == NULL) {
    dStiffdThick.zero();
    return;
  }

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);

  double x[3], y[3], z[3];

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
  x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;

  ShellMaterial<double> *mat = new ShellMaterialType0<double>(prop->E, prop->eh, prop->nu,
                                                              prop->rho, prop->Ta, prop->W); 

  double disp[18]; for(int i=0; i<18; ++i) disp[i] = 0;
  Impl::andesstfWRTthic(glNum+1, dStiffdThick.data(), (double*)NULL, prop->nu,
                        x, y, z, disp, 0, mat, flg);
  delete mat;
}

void
Membrane::getStiffnessNodalCoordinateSensitivity(FullSquareMatrix *&dStiffdx, CoordSet &cs, int)
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
                         (double*)NULL, x, y, z, 0, (double*)NULL,
                         flg);
}

void
Membrane::getVonMisesThicknessSensitivity(Vector &dStdThick, Vector &weight, CoordSet &cs,
                                          Vector &elDisp, int, int surface, int, double *ndTemps,
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

  ShellMaterial<double> *mat = new ShellMaterialType0<double>(prop->E, prop->eh, prop->nu, 
                                                              prop->rho, prop->Ta, prop->W);
  int sflg = 0;

  Impl::andesvmsWRTthic(glNum+1, prop->nu, x, y, z, elDisp.data(),
                        dStdThick.getData(), 0, mat, surface, sflg,
                        ndTemps);
  delete mat;
}

void
Membrane::getVonMisesNodalCoordinateSensitivity(GenFullM<double> &dStdx, Vector &weight, CoordSet &cs,
                                                Vector &elDisp, int, int surface, int, double *ndTemps,
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

  int sflg = 0;

  Impl::andesvmsWRTcoord(glNum+1, prop->E, prop->nu, prop->rho,
                         prop->eh, prop->Ta, prop->W, (double*)NULL,
                         x, y, z, elDisp.data(), dStdx.getData(),
                         0, (double*)NULL, surface, sflg, ndTemps);
}

void
Membrane::getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, GenFullM<double> *dDispDisp,
                                             CoordSet &cs, Vector &elDisp, int, int surface, int, double *ndTemps,
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

  ShellMaterial<double> *mat = new ShellMaterialType0<double>(prop->E, prop->eh, prop->nu, 
                                                              prop->rho, prop->Ta, prop->W);
  int sflg = 0;

  Impl::andesvmsWRTdisp(glNum+1, prop->nu, x, y, z, elDisp.data(),
                        dStdDisp.getData(), 0, mat, surface, sflg,
                        ndTemps);
  if(dDispDisp) dStdDisp ^= (*dDispDisp);

  delete mat;
}
#endif
