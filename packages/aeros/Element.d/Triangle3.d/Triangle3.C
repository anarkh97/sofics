#include <Element.d/Triangle3.d/Triangle3.h>
#include <Element.d/Triangle3.d/Triangle3StressWRTDisplacementSensitivity.h>
#include <Element.d/Function.d/SpaceDerivatives.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/matrix.h>
#include <Math.d/Vector.h>
#include <Utils.d/dofset.h>
#include <Utils.d/pstress.h>
#include <cmath>
#include <iostream>

extern int verboseFlag;

Triangle3::Triangle3(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
}

Element *
Triangle3::clone()
{
 return new Triangle3(*this);
}

void
Triangle3::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
}

void
Triangle3::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
}

void
Triangle3::getVonMises(Vector& stress,Vector& weight,CoordSet &cs, 
                       Vector& elDisp, int strInd,int,double *,
		       double ylayer, double zlayer, int avgnum)
{

        // NOTE: SIGMAZZ, SIGMAYZ, SIGMAXZ, STRAINZZ, STRAINYZ & STRAINXZ
        // ARE EQUAL TO ZERO FOR A 3 NODE TRIANGLE ELEMENT

        if(strInd == 2  || strInd == 4  || strInd == 5  ||
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

        double x[3], y[3];

        x[0] = nd1.x; y[0] = nd1.y;
        x[1] = nd2.x; y[1] = nd2.y;
        x[2] = nd3.x; y[2] = nd3.y;

        double area2 = ((x[1]*y[2]-x[2]*y[1])+
                        (x[2]*y[0]-x[0]*y[2])+
                        (x[0]*y[1]-x[1]*y[0]));

        double x21 = x[1] - x[0];
        double x32 = x[2] - x[1];
        double x13 = x[0] - x[2];

        double y12 = y[0] - y[1];
        double y23 = y[1] - y[2];
        double y31 = y[2] - y[0];

        double ux1 = elDisp[0];
        double uy1 = elDisp[1];
        double ux2 = elDisp[2];
        double uy2 = elDisp[3];
        double ux3 = elDisp[4];
        double uy3 = elDisp[5];

        double coef = 1.0/area2;
        double exx  = coef*(y23*ux1 + y31*ux2 + y12*ux3); 
        double eyy  = coef*(x32*uy1 + x13*uy2 + x21*uy3); 
        double exy  = coef*(x32*ux1 + y23*uy1 + x13*ux2 + 
                            y31*uy2 + x21*ux3 + y12*uy3); 

        double elStress[3][7];

        if(strInd > 6) {
          elStress[0][0] = exx;
          elStress[1][0] = exx;
          elStress[2][0] = exx;

          elStress[0][1] = eyy;
          elStress[1][1] = eyy;
          elStress[2][1] = eyy;

          elStress[0][2] = 0.0;
          elStress[1][2] = 0.0;
          elStress[2][2] = 0.0;

          elStress[0][3] = exy;
          elStress[1][3] = exy;
          elStress[2][3] = exy;

          elStress[0][4] = 0.0;
          elStress[1][4] = 0.0;
          elStress[2][4] = 0.0;

          elStress[0][5] = 0.0;
          elStress[1][5] = 0.0;
          elStress[2][5] = 0.0;

          stress[0] = elStress[0][strInd-7];
          stress[1] = elStress[1][strInd-7];
          stress[2] = elStress[2][strInd-7];

          return;
        }

        double E  = prop->E;
        double nu = prop->nu;

        double c[3][3];

        // Constitutive matrix
        c[0][0] =  E / (1.0-nu*nu);
        c[1][1] =  c[0][0];
        c[2][2] =  0.5*c[0][0]*(1.0-nu);
        c[0][1] =  c[0][0]*nu;
        c[1][0] =  c[0][1];
        c[0][2] =  0.0;
        c[1][2] =  0.0;
        c[2][0] =  0.0;
        c[2][1] =  0.0;

        double sxx = c[0][0]*exx + c[0][1]*eyy + c[0][2]*exy;
        double syy = c[1][0]*exx + c[1][1]*eyy + c[1][2]*exy;
        double sxy = c[2][0]*exx + c[2][1]*eyy + c[2][2]*exy;

        elStress[0][0] = sxx;
        elStress[1][0] = sxx;
        elStress[2][0] = sxx;

        elStress[0][1] = syy;
        elStress[1][1] = syy;
        elStress[2][1] = syy;

        elStress[0][2] = 0.0;
        elStress[1][2] = 0.0;
        elStress[2][2] = 0.0;

        elStress[0][3] = sxy;
        elStress[1][3] = sxy;
        elStress[2][3] = sxy;

        elStress[0][4] = 0.0;
        elStress[1][4] = 0.0;
        elStress[2][4] = 0.0;

        elStress[0][5] = 0.0;
        elStress[1][5] = 0.0;
        elStress[2][5] = 0.0;

        double von = 0.0;

        double stress_dxy = elStress[0][0] - elStress[0][1];
        double stress_dyz = elStress[0][1] - elStress[0][2];
        double stress_dxz = elStress[0][0] - elStress[0][2];

        // compute von mises stress!
        von = sqrt( 0.5*(stress_dxy*stress_dxy + stress_dyz*stress_dyz + stress_dxz*stress_dxz) +
                      3*(elStress[0][3]*elStress[0][3] + elStress[0][4]*elStress[0][4] + elStress[0][5]*elStress[0][5]) );

        elStress[0][6] = von;
        elStress[1][6] = von;
        elStress[2][6] = von;

        stress[0] = elStress[0][strInd];
        stress[1] = elStress[1][strInd];
        stress[2] = elStress[2][strInd];
	
}

void
Triangle3::getAllStress(FullM& stress,Vector& weight,CoordSet &cs, 
                       Vector& elDisp, int strInd,int,double *)
{

        // CALCULATE SIGMAXX, SIGMAYY, SIGMAXY, STRAINXX, STRAINYY, 
        // STRAINXY AND VONMISES STRESS

 	weight = 1.0;

	auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);

        double x[3], y[3];

        x[0] = nd1.x; y[0] = nd1.y;
        x[1] = nd2.x; y[1] = nd2.y;
        x[2] = nd3.x; y[2] = nd3.y;

        double area2 = ((x[1]*y[2]-x[2]*y[1])+
                        (x[2]*y[0]-x[0]*y[2])+
                        (x[0]*y[1]-x[1]*y[0]));

        double x21 = x[1] - x[0];
        double x32 = x[2] - x[1];
        double x13 = x[0] - x[2];

        double y12 = y[0] - y[1];
        double y23 = y[1] - y[2];
        double y31 = y[2] - y[0];

	double ux1 = elDisp[0];
	double uy1 = elDisp[1];
	double ux2 = elDisp[2];
	double uy2 = elDisp[3];
	double ux3 = elDisp[4];
	double uy3 = elDisp[5];

	double coef = 1.0/area2;
	double exx  = coef*(y23*ux1 + y31*ux2 + y12*ux3); 
	double eyy  = coef*(x32*uy1 + x13*uy2 + x21*uy3); 
	double exy  = coef*(x32*ux1 + y23*uy1 + x13*ux2 + 
                            y31*uy2 + x21*ux3 + y12*uy3); 

	double elStress[3][7];

	if(strInd != 0) {
	  elStress[0][0] = exx;
          elStress[1][0] = exx;
          elStress[2][0] = exx;

          elStress[0][1] = eyy;
          elStress[1][1] = eyy;
          elStress[2][1] = eyy;

          elStress[0][2] = 0.0;
          elStress[1][2] = 0.0;
          elStress[2][2] = 0.0;

          elStress[0][3] = exy;
          elStress[1][3] = exy;
          elStress[2][3] = exy;

          elStress[0][4] = 0.0;
          elStress[1][4] = 0.0;
          elStress[2][4] = 0.0;

          elStress[0][5] = 0.0;
          elStress[1][5] = 0.0;
          elStress[2][5] = 0.0;

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
          svec[3] /= 2;
          svec[4] /= 2;
          svec[5] /= 2;
          pstress(svec,pvec);
          for (i=0; i<3; ++i) {
            for (j=0; j<3; ++j) {
              stress[i][j+6] = pvec[j];
            }
          }
	  return;
        }

	double E  = prop->E;
	double nu = prop->nu;

        double c[3][3];

        // Constitutive matrix
        c[0][0] =  E / (1.0-nu*nu);
        c[1][1] =  c[0][0];
        c[2][2] =  0.5*c[0][0]*(1.0-nu);
        c[0][1] =  c[0][0]*nu;
        c[1][0] =  c[0][1];
        c[0][2] =  0.0;
        c[1][2] =  0.0;
        c[2][0] =  0.0;
        c[2][1] =  0.0;

	double sxx = c[0][0]*exx + c[0][1]*eyy + c[0][2]*exy;
	double syy = c[1][0]*exx + c[1][1]*eyy + c[1][2]*exy;
	double sxy = c[2][0]*exx + c[2][1]*eyy + c[2][2]*exy;

	elStress[0][0] = sxx;
	elStress[1][0] = sxx;
	elStress[2][0] = sxx;

	elStress[0][1] = syy;
	elStress[1][1] = syy;
	elStress[2][1] = syy;

	elStress[0][2] = 0.0;
	elStress[1][2] = 0.0;
	elStress[2][2] = 0.0;

	elStress[0][3] = sxy;
	elStress[1][3] = sxy;
	elStress[2][3] = sxy;

	elStress[0][4] = 0.0;
	elStress[1][4] = 0.0;
	elStress[2][4] = 0.0;

	elStress[0][5] = 0.0;
	elStress[1][5] = 0.0;
	elStress[2][5] = 0.0;

	double von = 0.0;

	// compute von mises stress!
	elStress[0][6] = von;
	elStress[1][6] = von;
	elStress[2][6] = von;

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
        pstress(svec,pvec);
        for (i=0; i<3; ++i) {
          for (j=0; j<3; ++j) {
            stress[i][j+6] = pvec[j];
          }
        }
}

double
Triangle3::getMass(const CoordSet& cs) const
{
        auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);

        double x[3], y[3];

        x[0] = nd1.x; y[0] = nd1.y;
        x[1] = nd2.x; y[1] = nd2.y;
        x[2] = nd3.x; y[2] = nd3.y;

	double area = 0.5*((x[1]*y[2]-x[2]*y[1])+
                           (x[2]*y[0]-x[0]*y[2])+
                           (x[0]*y[1]-x[1]*y[0]));

	double density = prop->rho;
	double t       = prop->eh;

	double mass = area*t*density;

        return mass;

}

double
Triangle3::getMassThicknessSensitivity(CoordSet& cs)
{
        auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);

        double x[3], y[3];

        x[0] = nd1.x; y[0] = nd1.y;
        x[1] = nd2.x; y[1] = nd2.y;
        x[2] = nd3.x; y[2] = nd3.y;

	double area = 0.5*((x[1]*y[2]-x[2]*y[1])+
                           (x[2]*y[0]-x[0]*y[2])+
                           (x[0]*y[1]-x[1]*y[0]));

	double density = prop->rho;
	double massWRTthic = area*density;
  return massWRTthic;
}

void
Triangle3::getGravityForce(CoordSet& cs,double *gravityAcceleration, 
                           Vector& gravityForce, int gravflg, GeomState *geomState)
{
        double massPerNode = getMass(cs)/3.0;

        double fx = massPerNode*gravityAcceleration[0];
        double fy = massPerNode*gravityAcceleration[1];

	gravityForce[0] = fx;
	gravityForce[1] = fy;
	gravityForce[2] = fx;
	gravityForce[3] = fy;
	gravityForce[4] = fx;
	gravityForce[5] = fy;
}

void
Triangle3::getGravityForceThicknessSensitivity(CoordSet& cs,double *gravityAcceleration,
                                               Vector& gravityForceSensitivity, int gravflg, GeomState *geomState)
{
  double massPerNode = getMass(cs)/3.0;
  
  double fx = massPerNode*gravityAcceleration[0];
  double fy = massPerNode*gravityAcceleration[1];

  gravityForceSensitivity[0] = fx;
  gravityForceSensitivity[1] = fy;
  gravityForceSensitivity[2] = fx;
  gravityForceSensitivity[3] = fy;
  gravityForceSensitivity[4] = fx;
  gravityForceSensitivity[5] = fy;
}

FullSquareMatrix
Triangle3::massMatrix(const CoordSet &cs,double *mel,int cmflg) const
{
	double mass = getMass(cs);
	double massPerNode = mass/3.0;

        FullSquareMatrix ret(6,mel);

	ret.zero();

	int i;
        for(i=0; i<6; ++i)
          ret[i][i] = massPerNode;

        return ret;
}

FullSquareMatrix
Triangle3::stiffness(const CoordSet &cs, double *d, int flg) const
{
	auto &nd1 = cs.getNode(nn[0]);
	auto &nd2 = cs.getNode(nn[1]);
	auto &nd3 = cs.getNode(nn[2]);

	double x[3], y[3];
	x[0] = nd1.x; y[0] = nd1.y; 
	x[1] = nd2.x; y[1] = nd2.y;
	x[2] = nd3.x; y[2] = nd3.y; 

	double area2 = ((x[1]*y[2]-x[2]*y[1])+
                        (x[2]*y[0]-x[0]*y[2])+
                        (x[0]*y[1]-x[1]*y[0]));

	double x21 = x[1] - x[0];
	double x32 = x[2] - x[1];
	double x13 = x[0] - x[2];

        double y12 = y[0] - y[1];
        double y23 = y[1] - y[2];
        double y31 = y[2] - y[0];

	double t  = prop->eh;
	double E  = prop->E;
	double nu = prop->nu;

	double c[3][3];

	// Constitutive matrix
        c[0][0] =  E / (1.0-nu*nu);
        c[1][1] =  c[0][0];
        c[2][2] =  0.5*c[0][0]*(1.0-nu);
        c[0][1] =  c[0][0]*nu;
        c[1][0] =  c[0][1];
        c[0][2] =  0.0;
        c[1][2] =  0.0;
        c[2][0] =  0.0;
        c[2][1] =  0.0;


        FullSquareMatrix K(6,d);

	K.zero();

	double ke = t/(4.0*area2);

        K[0][0] = ke*(c[2][2]*x32*x32 + c[0][0]*y23*y23);
        K[0][1] = ke*(c[0][1]+c[2][2])*x32*y23;
        K[0][2] = ke*(c[2][2]*x13*x32 + c[0][0]*y31*y23);
        K[0][3] = ke*(c[2][2]*x32*y31 + c[0][1]*x13*y23);
        K[0][4] = ke*(c[2][2]*x21*x32 + c[0][0]*y12*y23);
        K[0][5] = ke*(c[2][2]*x32*y12 + c[0][1]*x21*y23);

        K[1][1] = ke*(c[1][1]*x32*x32 + c[2][2]*y23*y23);
        K[1][2] = ke*(c[2][2]*x13*y23 + c[0][1]*x32*y31);
        K[1][3] = ke*(c[1][1]*x13*x32 + c[2][2]*y23*y31);
        K[1][4] = ke*(c[2][2]*x21*y23 + c[0][1]*x32*y12);
        K[1][5] = ke*(c[1][1]*x21*x32 + c[2][2]*y23*y12);

        K[2][2] = ke*(c[2][2]*x13*x13 + c[0][0]*y31*y31);
        K[2][3] = ke*(c[0][1]+c[2][2])*x13*y31;
        K[2][4] = ke*(c[2][2]*x13*x21 + c[0][0]*y12*y31);
        K[2][5] = ke*(c[2][2]*x13*y12 + c[0][1]*x21*y31);

        K[3][3] = ke*(c[1][1]*x13*x13 + c[2][2]*y31*y31);
        K[3][4] = ke*(c[2][2]*x21*y31 + c[0][1]*x13*y12);
        K[3][5] = ke*(c[1][1]*x13*x21 + c[2][2]*y12*y31);

        K[4][4] = ke*(c[2][2]*x21*x21 + c[0][0]*y12*y12);
        K[4][5] = ke*(c[0][1]+c[2][2])*x21*y12;

        K[5][5] = ke*(c[1][1]*x21*x21 + c[2][2]*y12*y12);

        K[1][0] = K[0][1];

        K[2][0] = K[0][2];
        K[2][1] = K[1][2];

        K[3][0] = K[0][3];
        K[3][1] = K[1][3];
        K[3][2] = K[2][3];

        K[4][0] = K[0][4];
        K[4][1] = K[1][4];
        K[4][2] = K[2][4];
        K[4][3] = K[3][4];

        K[5][0] = K[0][5];
        K[5][1] = K[1][5];
        K[5][2] = K[2][5];
        K[5][3] = K[3][5];
        K[5][4] = K[4][5];

        return K;
}

int
Triangle3::numNodes() const
{
 	return 3;
}

int*
Triangle3::nodes(int *p) const
{
 	if(p == 0) p = new int[3];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
	return p;
}

int
Triangle3::numDofs() const
{
 	return 6;
}

int*
Triangle3::dofs(DofSetArray &dsa, int *p) const
{
 	if(p == 0) p = new int[6];

	dsa.number(nn[0], DofSet::Xdisp | DofSet::Ydisp, p);
	dsa.number(nn[1], DofSet::Xdisp | DofSet::Ydisp, p+2);
	dsa.number(nn[2], DofSet::Xdisp | DofSet::Ydisp, p+4);

	return p;
}

void
Triangle3::markDofs(DofSetArray &dsa) const
{
        dsa.mark(nn, numNodes(), DofSet::Xdisp | DofSet::Ydisp);
}

int
Triangle3::getTopNumber() const
{
  return 104;//4;
}

void
Triangle3::getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, GenFullM<double> *dDispDisp,
                                              CoordSet &cs, Vector &elDisp, int strInd, int surface,
                                              double *ndTemps, int avgnum, double ylayer, double zlayer)
{ 
#ifdef USE_EIGEN3
  if(strInd != 6) {
    std::cerr << " ... Error: strInd must be 6 in FourNodeQuad::getVonMisesDisplacementSensitivity\n";
    exit(-1); 
  } 
  if(dStdDisp.numRow() != 3 || dStdDisp.numCol() !=6) {
    std::cerr << " ... Error: dimenstion of sensitivity matrix is wrong\n";
    exit(-1);
  } 
  weight = 1;
  // scalar parameters
  Eigen::Array<double,8,1> dconst;
  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);
  auto &nd3 = cs.getNode(nn[2]);
  
  double x[3], y[3];
  
  x[0] = nd1.x; y[0] = nd1.y; 
  x[1] = nd2.x; y[1] = nd2.y;
  x[2] = nd3.x; y[2] = nd3.y;

  dconst[0] = nd1.x; dconst[1] = nd2.x; dconst[2] = nd3.x; // x coordinates
  dconst[3] = nd1.y; dconst[4] = nd2.y; dconst[5] = nd3.y; // y coordinates
  dconst[6] = prop->E;
  dconst[7] = prop->nu;

  // integer parameters
  Eigen::Array<int,1,1> iconst;
  iconst[0] = avgnum;
  // inputs
  Eigen::Matrix<double,6,1> q = Eigen::Map<Eigen::Matrix<double,6,1> >(elDisp.data()).segment(0,6); // displacements

  //Jacobian evaluation
  Eigen::Matrix<double,3,6> dStressdDisp;
  dStressdDisp.setZero();
  vms4WRTdisp(x, y, q.data(), 
              dStressdDisp.data(), 
              prop->E, prop->nu); 
  dStdDisp.copy(dStressdDisp.data());
#else
  std::cerr << " ... Error! Triangle3::getVonMisesDisplacementSensitivity needs Eigen library\n";
  exit(-1);
#endif
  if(dDispDisp) dStdDisp ^= (*dDispDisp);
}
