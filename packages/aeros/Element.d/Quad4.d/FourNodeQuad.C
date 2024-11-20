#include 	<Element.d/Quad4.d/FourNodeQuad.h>
#include 	<Element.d/Quad4.d/FourNodeQuadStressWRTDisplacementSensitivity.h>
#include 	<Element.d/Function.d/SpaceDerivatives.h>
#include        <Math.d/Vector.h>
#include 	<Math.d/FullSquareMatrix.h>
#include 	<Math.d/matrix.h>
#include 	<Utils.d/dofset.h>
#include 	<Element.d/State.h>
#include 	<Utils.d/linkfc.h>
#include 	<Utils.d/pstress.h>
#include 	<Hetero.d/InterpPoint.h>

extern int verboseFlag;
extern "C"      {
void    _FORTRAN(getcmt)(double&, double&, double&, double* );
void    _FORTRAN(quad4m)(double*, double*, double*, double*, const int&,
                         double*, const int&, double&);
void    _FORTRAN(q4dmas)(const int&, double&, double*, const int&, 
                         double*,double*,
                         double*, double*,double*,int&,double&,int&);
void _FORTRAN(sands2)(char*,double*,double*,double*,double*,double*,double*,
            const int&, const int&, const int&, const int&, const int&,
            const int&, double&, double&, double*);
void _FORTRAN(quad2d)(double*, double*, double&, const int*, 
                      double[4][8], const int*);

void _FORTRAN(qgauss)(int &, int &, int &, int &,
                      double &,  double &, double &);

void _FORTRAN(q4shpe)(double &, double &, double *, double *,
                      double *, double *, double *, double &);
}

FourNodeQuad::FourNodeQuad(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
	nn[3] = nodenums[3];
}

void
FourNodeQuad::spy()
{
 fprintf(stderr, "Spy %p with %d %d %d %d\n",
             this, nn[0],nn[1],nn[2],nn[3]);
}

Element *
FourNodeQuad::clone()
{
 return new FourNodeQuad(*this);
}

void
FourNodeQuad::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
}

void
FourNodeQuad::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
}

void
FourNodeQuad::getVonMises(Vector& stress,Vector& weight,CoordSet &cs, 
                          Vector& elDisp, int strInd,int,double *ndTemps,
                          double ylayer, double zlayer, int avgnum)
{
         
        // NOTE: SIGMAZZ, SIGMAYZ, SIGMAXZ, STRAINZZ, STRAINYZ & STRAINXZ
        // ARE EQUAL TO ZERO FOR A 4 NODE QUAD ELEMENT

        if(strInd == 2 || strInd == 4  || strInd == 5  ||
           strInd == 11 || strInd == 12   ) {

           weight = 0.0; stress = 0.0;

           return;
        }

        // ELSE CALCULATE SIGMAXX, SIGMAYY, SIGMAXY, STRAINXX, STRAINYY, 
        // STRAINXY AND VONMISES STRESS

        weight = 1.0;
        if(strInd == -1) return;

        auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);
        auto &nd4 = cs.getNode(nn[3]);

        double x[4], y[4], c[9];

        x[0] = nd1.x; y[0] = nd1.y;
        x[1] = nd2.x; y[1] = nd2.y;
        x[2] = nd3.x; y[2] = nd3.y;
        x[3] = nd4.x; y[3] = nd4.y;

       _FORTRAN(getcmt)(prop->A, prop->E, prop->nu, c);

       int maxgus = 4;
       int maxstr = 7;
       int elm    = 1;
       int numel  = 1;

       int vmflg = 0;  // PJSA init to 0
       if(strInd == 6) vmflg = 1;

       int strainFlg = 0; // PJSA init to 0
       if(strInd == 13) strainFlg = 1;

       double elStrain[4][7], elStress[4][7];
       double alpha = prop->W;
       double E     = prop->E;
       double nu    = prop->nu;
       double Tref  = prop->Ta;

       double tc = E*alpha/(1-nu);


       //char ESCM[7] = "DIRECT"; // ... DIRECT STRESS VALUE CALCULATION
       char ESCM[7] = "EXTRAP";   // ... STRESS EXTRAPOLATION FROM GAUSS POINTS

      _FORTRAN(sands2)(ESCM,x,y,c,elDisp.data(),(double*)elStress,
                       (double*)elStrain, maxgus,maxstr,elm,numel,
                       vmflg,strainFlg,tc,Tref,ndTemps); 

// if strInd <= 6, you are retrieving a stress value:
// if strInd >  6, you are retrieving a strain value:

      if(strInd <= 6) {
        stress[0] = elStress[0][strInd];
        stress[1] = elStress[1][strInd]; 
        stress[2] = elStress[2][strInd]; 
        stress[3] = elStress[3][strInd]; 
      }
      else {
        stress[0] = elStrain[0][strInd-7];
        stress[1] = elStrain[1][strInd-7]; 
        stress[2] = elStrain[2][strInd-7]; 
        stress[3] = elStrain[3][strInd-7]; 
      } 

        // ... TO COMPUTE OUT OF PLANE STRAIN
        // ... WHICH IS NOT NECESSARILY ZERO
	if(strInd == 9) {
	  double nu  = prop->nu;
	  double ezz = -nu/(1.0-nu)*(elStrain[0][0] + elStrain[0][1]); 
	  stress     = ezz;
	}

}

void
FourNodeQuad::getAllStress(FullM& stress,Vector& weight,CoordSet &cs,
                          Vector& elDisp, int strInd,int,double *ndTemps)
{

        // CALCULATE SIGMAXX, SIGMAYY, SIGMAXY, STRAINXX, STRAINYY,
        // STRAINXY AND VONMISES STRESS

        weight = 1.0;

        auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);
        auto &nd4 = cs.getNode(nn[3]);

        double x[4], y[4], c[9];

        x[0] = nd1.x; y[0] = nd1.y;
        x[1] = nd2.x; y[1] = nd2.y;
        x[2] = nd3.x; y[2] = nd3.y;
        x[3] = nd4.x; y[3] = nd4.y;

	_FORTRAN(getcmt)(prop->A, prop->E, prop->nu, c);

       int maxgus = 4;
       int maxstr = 7;
       int elm    = 1;
       int numel  = 1;

       int vmflg,strainFlg;
       vmflg  = 0;
       strainFlg = 0;

       double elStrain[4][7], elStress[4][7];
       double alpha = prop->W;
       double E     = prop->E;
       double nu    = prop->nu;
       double Tref  = prop->Ta;

       double tc = E*alpha/(1-nu);

       //char ESCM[7] = "DIRECT"; // ... DIRECT STRESS VALUE CALCULATION
       char ESCM[7] = "EXTRAP";   // ... STRESS EXTRAPOLATION FROM GAUSS POINTS

      _FORTRAN(sands2)(ESCM,x,y,c,elDisp.data(),(double*)elStress,
                       (double*)elStrain, maxgus,maxstr,elm,numel,
                       vmflg,strainFlg,tc,Tref,ndTemps);

// Store all Stress or all Strain as defined by strInd
        int i,j;
        if(strInd == 0) {
          for (i=0; i<4; ++i) {
            for (j=0; j<6; ++j) {
              stress[i][j] = elStress[i][j];
            }
          }
        }
        else {
          double nu  = prop->nu;
          double ezz = -nu/(1.0-nu)*(elStrain[0][0] + elStrain[0][1]);
          for (i=0; i<4; ++i) {
            for (j=0; j<6; ++j) {
              stress[i][j] = elStrain[i][j];
            }
            stress[i][2] = ezz;
          }
        }

// Get Element Principals
        double svec[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
        double pvec[3] = {0.0,0.0,0.0};
        for (j=0; j<6; ++j) {
          for (i=0; i<4; ++i) {
            svec[j] += stress[i][j];
          }
          svec[j] /= 4;
        }
// Convert Engineering to Tensor Strains
        if(strInd != 0) {
          svec[3] /= 2;
          svec[4] /= 2;
          svec[5] /= 2;
        }
        pstress(svec,pvec);
        for (i=0; i<4; ++i) {
          for (j=0; j<3; ++j) {
            stress[i][j+6] = pvec[j];
          }
        }
}

double
FourNodeQuad::getMass(const CoordSet& cs) const
{
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

double
FourNodeQuad::getMassSensitivityWRTthickness(CoordSet& cs)
{
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
        double massWRTthic = area*prop->rho;

        return massWRTthic;
}

double
FourNodeQuad::weight(CoordSet& cs, double *gravityAcceleration)
{
  if (prop == NULL) {
    return 0.0;
  }

  double _mass = getMass(cs);
  double gravAccNorm = sqrt(gravityAcceleration[0]*gravityAcceleration[0] + 
                            gravityAcceleration[1]*gravityAcceleration[1] +
                            gravityAcceleration[2]*gravityAcceleration[2]);
  return _mass*gravAccNorm;
}

double
FourNodeQuad::weightDerivativeWRTthickness(CoordSet& cs, double *gravityAcceleration, int senMethod)
{
  if (prop == NULL) {
    return 0.0;
  }
 
  if (senMethod == 0) {
    double _weight = weight(cs, gravityAcceleration);
    double thick = prop->eh;
    return _weight/thick;
  } else {
    fprintf(stderr," ... Error: FourNodeQuad::weightDerivativeWRTthickness for automatic differentiation and finite difference is not implemented\n");
    exit(-1);
  }
}

void
FourNodeQuad::getGravityForce(CoordSet& cs,double *gravityAcceleration, 
                              Vector& gravityForce, int gravflg, GeomState *geomState)
{

  // Lumped
  if (gravflg != 2) {
    double massPerNode = 0.25*getMass(cs);

    double fx = massPerNode*gravityAcceleration[0];
    double fy = massPerNode*gravityAcceleration[1];

    gravityForce[0] = fx;
    gravityForce[1] = fy;
    gravityForce[2] = fx;
    gravityForce[3] = fy;
    gravityForce[4] = fx;
    gravityForce[5] = fy;
    gravityForce[6] = fx;
    gravityForce[7] = fy;
  }
  // Consistent
  else {

    int i;
    int numgauss = 2;
    auto &nd1 = cs.getNode(nn[0]);
    auto &nd2 = cs.getNode(nn[1]);
    auto &nd3 = cs.getNode(nn[2]);
    auto &nd4 = cs.getNode(nn[3]);

    double x[4], y[4];
    x[0] = nd1.x; y[0] = nd1.y;
    x[1] = nd2.x; y[1] = nd2.y;
    x[2] = nd3.x; y[2] = nd3.y;
    x[3] = nd4.x; y[3] = nd4.y;
    double rho = prop->rho;
    double h = prop->eh;

    double lforce[4];
    for (i = 0; i < 4; ++i)
      lforce[i] = 0.0;

    int fortran = 1;  // fortran routines start from index 1
    int pt1, pt2;
    for (pt1 = 0 + fortran; pt1 < numgauss + fortran; pt1++)  {
      for (pt2 = 0 + fortran; pt2 < numgauss + fortran; pt2++)  {
        // get gauss point
        double xi, eta, wt;
        _FORTRAN(qgauss)(numgauss, pt1, numgauss, pt2, xi, eta, wt);

        //compute shape functions
        double shapeFunc[4], shapeGradX[4], shapeGradY[4];
        double detJ;  //det of jacobian

        _FORTRAN(q4shpe)(xi, eta, x, y,
                         shapeFunc, shapeGradX, shapeGradY, detJ);

        for (i = 0; i < 4; ++i)
          lforce[i] += wt*shapeFunc[i]*detJ;
      }
    }
    for(i=0; i<4; ++i) {
      gravityForce[2*i+0] = lforce[i]*h*rho*gravityAcceleration[0];
      gravityForce[2*i+1] = lforce[i]*h*rho*gravityAcceleration[1];
    }
  }
}

void
FourNodeQuad::getGravityForceSensitivityWRTthickness(CoordSet& cs, double *gravityAcceleration, int senMethod,
                                                     Vector& gravityForceSensitivity, int gravflg, GeomState *geomState)
{

  // Lumped
  if (gravflg != 2) {
    double massPerNodePerThickness = 0.25*getMass(cs)/prop->eh;

    double fx = massPerNodePerThickness*gravityAcceleration[0];
    double fy = massPerNodePerThickness*gravityAcceleration[1];

    gravityForceSensitivity[0] = fx;
    gravityForceSensitivity[1] = fy;
    gravityForceSensitivity[2] = fx;
    gravityForceSensitivity[3] = fy;
    gravityForceSensitivity[4] = fx;
    gravityForceSensitivity[5] = fy;
    gravityForceSensitivity[6] = fx;
    gravityForceSensitivity[7] = fy;
  }
  // Consistent
  else {

    int i;
    int numgauss = 2;
    auto &nd1 = cs.getNode(nn[0]);
    auto &nd2 = cs.getNode(nn[1]);
    auto &nd3 = cs.getNode(nn[2]);
    auto &nd4 = cs.getNode(nn[3]);

    double x[4], y[4];
    x[0] = nd1.x; y[0] = nd1.y;
    x[1] = nd2.x; y[1] = nd2.y;
    x[2] = nd3.x; y[2] = nd3.y;
    x[3] = nd4.x; y[3] = nd4.y;
    double rho = prop->rho;
    double h = prop->eh;

    double lforce[4];
    for (i = 0; i < 4; ++i)
      lforce[i] = 0.0;

    int fortran = 1;  // fortran routines start from index 1
    int pt1, pt2;
    for (pt1 = 0 + fortran; pt1 < numgauss + fortran; pt1++)  {
      for (pt2 = 0 + fortran; pt2 < numgauss + fortran; pt2++)  {
        // get gauss point
        double xi, eta, wt;
        _FORTRAN(qgauss)(numgauss, pt1, numgauss, pt2, xi, eta, wt);

        //compute shape functions
        double shapeFunc[4], shapeGradX[4], shapeGradY[4];
        double detJ;  //det of jacobian

        _FORTRAN(q4shpe)(xi, eta, x, y,
                         shapeFunc, shapeGradX, shapeGradY, detJ);

        for (i = 0; i < 4; ++i)
          lforce[i] += wt*shapeFunc[i]*detJ;
      }
    }
    for(i=0; i<4; ++i) {
      gravityForceSensitivity[2*i+0] = lforce[i]*rho*gravityAcceleration[0];
      gravityForceSensitivity[2*i+1] = lforce[i]*rho*gravityAcceleration[1];
    }
  }
}

FullSquareMatrix
FourNodeQuad::massMatrix(const CoordSet &cs,double *mel,int cmflg) const
{
        auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);
        auto &nd4 = cs.getNode(nn[3]);

        double x[4], y[4], h[4];

        x[0] = nd1.x; y[0] = nd1.y;
        x[1] = nd2.x; y[1] = nd2.y;
        x[2] = nd3.x; y[2] = nd3.y;
        x[3] = nd4.x; y[3] = nd4.y;

        h[0] = h[1] = h[2] = h[3] = prop->eh;

        double *gravityAcceleration = 0;
        double *grvfor = 0;
        double totmas = 0.0;

        int grvflg = 0, masflg = 0;

        const int numgauss = 2;
        const int numdof   = 8;

       _FORTRAN(q4dmas)(numgauss,prop->rho,(double*)mel,numdof,h,x,y,
                        gravityAcceleration,grvfor,grvflg,totmas,masflg);

        FullSquareMatrix ret(8,mel);

        return ret;
}

FullSquareMatrix
FourNodeQuad::stiffness(const CoordSet &cs, double *d, int flg) const
{
	auto &nd1 = cs.getNode(nn[0]);
	auto &nd2 = cs.getNode(nn[1]);
	auto &nd3 = cs.getNode(nn[2]);
	auto &nd4 = cs.getNode(nn[3]);

	double c[9];

	double x[4], y[4], h[4];

	x[0] = nd1.x; y[0] = nd1.y; 
	x[1] = nd2.x; y[1] = nd2.y;
	x[2] = nd3.x; y[2] = nd3.y; 
	x[3] = nd4.x; y[3] = nd4.y;

	h[0] = h[1] = h[2] = h[3] = prop->eh;

	_FORTRAN(getcmt)(prop->A, prop->E, prop->nu, c);

	const int numgauss = 2;
        const int numdof   = 8;

	_FORTRAN(quad4m)(x, y, h, c, numgauss, static_cast<double *>(d), numdof, prop->A);

        FullSquareMatrix ret(8,d);

        return ret;
}

int
FourNodeQuad::numNodes() const
{
 	return 4;
}

int*
FourNodeQuad::nodes(int *p) const
{
 	if(p == 0) p = new int[4];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
 	p[3] = nn[3];
	return p;
}

int
FourNodeQuad::numDofs() const
{
 	return 8;
}

int*
FourNodeQuad::dofs(DofSetArray &dsa, int *p) const
{
 	if(p == 0) p = new int[8];

	dsa.number(nn[0],DofSet::Xdisp | DofSet::Ydisp, p);
	dsa.number(nn[1],DofSet::Xdisp | DofSet::Ydisp, p+2);
	dsa.number(nn[2],DofSet::Xdisp | DofSet::Ydisp, p+4);
	dsa.number(nn[3],DofSet::Xdisp | DofSet::Ydisp, p+6);

	return p;
}

void
FourNodeQuad::markDofs(DofSetArray &dsa) const
{
        dsa.mark(nn, 4, DofSet::Xdisp | DofSet::Ydisp);
}

int
FourNodeQuad::getTopNumber() const
{
  return 102;
}


void
FourNodeQuad::computeDisp(CoordSet&cs, State &state, const InterpPoint &ip,
                          double*res, GeomState *gs)
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
FourNodeQuad::getFlLoad(CoordSet &cs, const InterpPoint &ip, double *flF, 
                        double *resF, GeomState *gs)
{
 const double *gp = ip.xy;
//i=2: Fx and Fy for each node
 int i;
 for(i = 0; i < 2; ++i) {
   resF[i]    = (1-gp[0])*(1-gp[1])* flF[i];
   resF[2+i]  = gp[0]*(1-gp[1])* flF[i];
   resF[6+i]  = (1-gp[0])*gp[1]* flF[i];
   resF[4+i]  = gp[0]*gp[1]* flF[i];
  }
// fprintf(stderr,"Forces : %d %f %f %f %f\n", i, resF[i], resF[2+i], resF[6+i], resF[4+i]);

}

void
FourNodeQuad::getThermalForce(CoordSet &cs, Vector &ndTemps,
                            Vector &elementThermalForce, int glflag, GeomState *geomState)
{
// Computes the thermal-mechanical coupling force C*theta
// The matrix C[8x4] is computed in quad2d.f
// Called from Dynam.C, since it is a time-dependant RHS force.

        auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);
        auto &nd4 = cs.getNode(nn[3]);

        double x[4], y[4], elC[4][8];
// BE CAREFUL!!! Fortran=(8,4)--> C++ = (4X8)

        x[0] = nd1.x; y[0] = nd1.y;
        x[1] = nd2.x; y[1] = nd2.y;
        x[2] = nd3.x; y[2] = nd3.y;
        x[3] = nd4.x; y[3] = nd4.y;

        const int numgauss = 2;
        const int numdof   = 4;

        double h = prop ->eh;
        double alpha = prop->W;
        double E     = prop->E;
        double nu    = prop->nu;
        double Tref  = prop->Ta;

//        double c = E*alpha*h/(1-2.0*nu);
        double c = E*alpha*h/(1-nu);

//        ndTemps.print("-----------");

        _FORTRAN(quad2d)(x, y, c, &numgauss, elC, &numdof);

// For Quads, C is a 8x4 matrix

        int i, j;
        for (i=0; i<8; i++){
          elementThermalForce[i] = 0.0;
          for(j=0; j<4; j++)
           elementThermalForce[i] += elC[j][i]*(ndTemps[j] - Tref);
         }
//         elementThermalForce.print("FourNodeQuad");
}

void
FourNodeQuad::getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, GenFullM<double> *dDispDisp, CoordSet &cs, Vector &elDisp, int strInd, int surface,
                                                 int senMethod, double *ndTemps, int avgnum, double ylayer, double zlayer)
{
#ifdef USE_EIGEN3
  if(strInd != 6) {
    std::cerr << " ... Error: strInd must be 6 in FourNodeQuad::getVonMisesDisplacementSensitivity\n";
    exit(-1);
  }
  if(dStdDisp.numRow() != 4 || dStdDisp.numCol() !=8) {
    std::cerr << " ... Error: dimenstion of sensitivity matrix is wrong\n";
    exit(-1);
  }
  weight = 1;
  // scalar parameters
  Eigen::Array<double,17,1> dconst;
  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);
  auto &nd4 = cs.getNode(nn[3]);

  double x[4], y[4];

  x[0] = nd1.x; y[0] = nd1.y; 
  x[1] = nd2.x; y[1] = nd2.y; 
  x[2] = nd3.x; y[2] = nd3.y; 
  x[3] = nd4.x; y[3] = nd4.y; 

  dconst[0] = nd1.x; dconst[1] = nd2.x; dconst[2] = nd3.x; dconst[3] = nd4.x; // x coordinates
  dconst[4] = nd1.y; dconst[5] = nd2.y; dconst[6] = nd3.y; dconst[7] = nd4.y; // y coordinates
  dconst[8] = prop->E;
  dconst[9] = prop->A;
  dconst[10] = prop->nu;
  dconst[11] = prop->W;
  dconst[12] = prop->Ta;
  if(ndTemps) {
    dconst[13] = ndTemps[0];
    dconst[14] = ndTemps[1];
    dconst[15] = ndTemps[2];
    dconst[16] = ndTemps[3];
  } else {
    dconst[13] = 0;
    dconst[14] = 0;
    dconst[15] = 0;
    dconst[16] = 0;
  }

  // integer parameters
  Eigen::Array<int,1,1> iconst;
  iconst[0] = avgnum;
  // inputs
  Eigen::Matrix<double,8,1> q = Eigen::Map<Eigen::Matrix<double,8,1> >(elDisp.data()).segment(0,8); // displacements

  //Jacobian evaluation
  Eigen::Matrix<double,4,8> dStressdDisp;
  Eigen::Matrix<double,7,3> stress;

  if(senMethod == 1) { // via automatic differentiation
#ifndef AEROS_NO_AD 
    Simo::Jacobian<double,FourNodeQuadStressWRTDisplacementSensitivity> dSdu(dconst,iconst);
    dStressdDisp = dSdu(q, 0);
    dStdDisp.copy(dStressdDisp.data());
#else
  std::cerr << " ... Error: AEROS_NO_AD is defined in FourNodeQuad::getVonMisesDisplacementSensitivity\n";   exit(-1);
#endif
  }
 
  if(senMethod == 0) { // analytic
    dStressdDisp.setZero();
    char escm[7] = "extrap";
    int numel = 1;
    int elm = 1;
    int maxgus = 4;
    int maxstr = 7;
    bool vmflg = true;
    bool strainFlg = false;
    double c[9];
    getcmt(prop->A, prop->E, prop->nu, c);
    double tc = prop->E*prop->W/(1.0-prop->nu);
      
    Eigen::Matrix<double,4,1> ndtemps = Eigen::Map<Eigen::Matrix<double,17,1> >(dconst.data()).segment(13,4); // extract eframe
    vms2WRTdisp(escm, x, y, c, q.data(), 
                dStressdDisp.data(), 0, 
                maxgus, maxstr, elm, numel, vmflg, 
                strainFlg, tc, prop->Ta, ndtemps.data());
    dStdDisp.copy(dStressdDisp.data());
  }

  if(senMethod == 2) { // via finite difference
    FourNodeQuadStressWRTDisplacementSensitivity<double> foo(dconst,iconst);
    double h = 1.0e-6;
    for(int j=0; j<8; ++j) {
      Eigen::Matrix<double,8,1> q_plus(q);
      Eigen::Matrix<double,8,1> q_minus(q);
      q_plus[j] += h;  q_minus[j] -= h;
      Eigen::Matrix<double,4,1> S_plus = foo(q_plus,0);   
      Eigen::Matrix<double,4,1> S_minus = foo(q_minus,0);
      Eigen::Matrix<double,4,1> dS = (S_plus-S_minus)/(2*h);
      dStressdDisp(0,j) = dS[0];
      dStressdDisp(1,j) = dS[1];
      dStressdDisp(2,j) = dS[2];
      dStressdDisp(3,j) = dS[3];
    }
    dStdDisp.copy(dStressdDisp.data());
  }
#else
  std::cerr << " ... ERROR! FourNodeQuad::getVonMisesDisplacementSensitivity needs Eigen library.\n";
  exit(-1);
#endif
  if(dDispDisp) dStdDisp ^= (*dDispDisp);
}
