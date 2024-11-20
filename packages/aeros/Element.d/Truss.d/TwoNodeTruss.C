#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include <Element.d/Truss.d/TwoNodeTruss.h>
#include <Element.d/Truss.d/TwoNodeTrussStiffnessWRTNodalCoordinateSensitivity.h>
#include <Element.d/Truss.d/TwoNodeTrussStressWRTNodalCoordinateSensitivity.h>
#include <Element.d/Function.d/SpaceDerivatives.h>
#include <Corotational.d/BarCorotator.h>
#include <Corotational.d/GeomState.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/matrix.h>
#include <Corotational.d/utilities.h>
#include <Utils.d/dofset.h>
#include <Utils.d/linkfc.h>

extern int verboseFlag;

extern "C" {
 void _FORTRAN(transform)(double*, double*, double*, double*, double*, double*, double*);
}

TwoNodeTruss::TwoNodeTruss(int* nodenums)
{
        nn[0] = nodenums[0];
        nn[1] = nodenums[1];
        preload = 0.0;
}

void 
TwoNodeTruss::setPreLoad(std::vector<double> &load)
{ 
        preload = load[0];
}

std::vector<double>
TwoNodeTruss::getPreLoad()
{
        std::vector<double> ret;
        ret.push_back(preload);
        return ret;
}

Element *
TwoNodeTruss::clone()
{
	return new TwoNodeTruss(*this);
}

void
TwoNodeTruss::renum(const int *table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
}

void
TwoNodeTruss::renum(EleRenumMap& table)
{
	nn[0] = table[nn[0]];
	nn[1] = table[nn[1]];
}

void
TwoNodeTruss::getIntrnForce(Vector& elForce, CoordSet& cs,
			    double *elDisp, int forceIndex, double *ndTemps)
{
        using std::sqrt; 

        // ... BARS ONLY CARRY AXIAL FORCES
        if(forceIndex > 0) {
          elForce[0] = 0.0;
          elForce[1] = 0.0;
          return;
        }

        auto &nd1 = cs.getNode( nn[0] );
        auto &nd2 = cs.getNode( nn[1] );

	double dx = nd2.x - nd1.x;
	double dy = nd2.y - nd1.y;
	double dz = nd2.z - nd1.z;

	double length = sqrt(dx*dx + dy*dy + dz*dz);

	// scale dx, dy, and dz by the length
	dx /= length;
	dy /= length;
	dz /= length;

	// Compute the change in length of the element
	double dq = dx*(elDisp[3]-elDisp[0]) 
                  + dy*(elDisp[4]-elDisp[1]) 
                  + dz*(elDisp[5]-elDisp[2]);

	// Compute axial strain
	double exx = dq/length;

	// Compute force
	double f = prop->A*prop->E*exx;

        // Add Preload
        f  += preload;

        // Compute thermal force
        double coefficient = prop->A*prop->E*prop->W;
        double Tref = prop->Ta;

        double fth1 = coefficient*(ndTemps[0]-Tref);
        double fth2 = coefficient*(ndTemps[1]-Tref);

        // If We take the mean temperature:
//      double fth = coefficient*((ndTemps[0]+ndTemps[1])/2-Tref);

	// return nodal forces
        elForce[0] = -f + fth1;
        elForce[1] =  f - fth2;
}

double
TwoNodeTruss::getMass(const CoordSet& cs) const
{
        using std::sqrt;

        auto &nd1 = cs.getNode( nn[0] );
        auto &nd2 = cs.getNode( nn[1] );

        double x[2], y[2], z[2];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

        double dx = x[1] - x[0];
        double dy = y[1] - y[0];
        double dz = z[1] - z[0];

        double length = sqrt( dx*dx + dy*dy + dz*dz );
        double   mass = (length*prop->A*prop->rho);

        return mass;
}

void
TwoNodeTruss::getMassNodalCoordinateSensitivity(CoordSet &cs, Vector &dMassdx)
{
        using std::sqrt;

        if(dMassdx.size() != 6) {
          std::cerr << " ... Error: dimension of sensitivity matrix is wrong\n";
          exit(-1);
        }

        auto &nd1 = cs.getNode( nn[0] );
        auto &nd2 = cs.getNode( nn[1] );

        double x[2], y[2], z[2];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

        double dx = x[1] - x[0];
        double dy = y[1] - y[0];
        double dz = z[1] - z[0];

        double length = sqrt( dx*dx + dy*dy + dz*dz );
        double prefix = prop->A*prop->rho/length;
        dMassdx[0] = -prefix*dx;
        dMassdx[1] = -prefix*dy;
        dMassdx[2] = -prefix*dz;
        dMassdx[3] = prefix*dx;
        dMassdx[4] = prefix*dy;
        dMassdx[5] = prefix*dz;       
}

void
TwoNodeTruss::getLengthNodalCoordinateSensitivity(CoordSet &cs, Vector &dLengthdx)
{
        using std::sqrt;

        if(dLengthdx.size() != 6) {
          std::cerr << " ... Error: dimension of sensitivity matrix is wrong\n";
          exit(-1);
        }

        auto &nd1 = cs.getNode( nn[0] );
        auto &nd2 = cs.getNode( nn[1] );

        double x[2], y[2], z[2];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

        double dx = x[1] - x[0];
        double dy = y[1] - y[0];
        double dz = z[1] - z[0];

        double length = sqrt( dx*dx + dy*dy + dz*dz );
        dLengthdx[0] = -dx/length;
        dLengthdx[1] = -dy/length;
        dLengthdx[2] = -dz/length;
        dLengthdx[3] = dx/length;
        dLengthdx[4] = dy/length;
        dLengthdx[5] = dz/length;        
}

void
TwoNodeTruss::getWeightNodalCoordinateSensitivity(Vector &dwdx, CoordSet& cs, double *gravityAcceleration)
{
  if(dwdx.size() != 6) {
     std::cerr << " ... Error: dimension of sensitivity matrix is wrong\n";
     exit(-1);
  }
  Vector dMassdx(6);
  getMassNodalCoordinateSensitivity(cs, dMassdx);
  double gravAccNorm = sqrt(gravityAcceleration[0]*gravityAcceleration[0] +
                            gravityAcceleration[1]*gravityAcceleration[1] +
                            gravityAcceleration[2]*gravityAcceleration[2]);
  dwdx = gravAccNorm*dMassdx;
}

void
TwoNodeTruss::getGravityForceNodalCoordinateSensitivity(CoordSet& cs, double *gravityAcceleration,
                                                        GenFullM<double> &dGfdx, int gravflg, GeomState *geomState)
{
       Vector dMassdx(6);
       getMassNodalCoordinateSensitivity(cs, dMassdx);
       dGfdx[0][0] = 0.5*gravityAcceleration[0]*dMassdx[0];
       dGfdx[0][1] = 0.5*gravityAcceleration[1]*dMassdx[0];
       dGfdx[0][2] = 0.5*gravityAcceleration[2]*dMassdx[0];
       dGfdx[0][3] = dGfdx[0][0]; 
       dGfdx[0][4] = dGfdx[0][1];
       dGfdx[0][5] = dGfdx[0][2];

       dGfdx[1][0] = 0.5*gravityAcceleration[0]*dMassdx[1];
       dGfdx[1][1] = 0.5*gravityAcceleration[1]*dMassdx[1];
       dGfdx[1][2] = 0.5*gravityAcceleration[2]*dMassdx[1];
       dGfdx[1][3] = dGfdx[1][0]; 
       dGfdx[1][4] = dGfdx[1][1];
       dGfdx[1][5] = dGfdx[1][2];
      
       dGfdx[2][0] = 0.5*gravityAcceleration[0]*dMassdx[2];
       dGfdx[2][1] = 0.5*gravityAcceleration[1]*dMassdx[2];
       dGfdx[2][2] = 0.5*gravityAcceleration[2]*dMassdx[2];
       dGfdx[2][3] = dGfdx[2][0]; 
       dGfdx[2][4] = dGfdx[2][1];
       dGfdx[2][5] = dGfdx[2][2];
      
       dGfdx[3][0] = 0.5*gravityAcceleration[0]*dMassdx[3];
       dGfdx[3][1] = 0.5*gravityAcceleration[1]*dMassdx[3];
       dGfdx[3][2] = 0.5*gravityAcceleration[2]*dMassdx[3];
       dGfdx[3][3] = dGfdx[3][0]; 
       dGfdx[3][4] = dGfdx[3][1];
       dGfdx[3][5] = dGfdx[3][2];

       dGfdx[4][0] = 0.5*gravityAcceleration[0]*dMassdx[4];
       dGfdx[4][1] = 0.5*gravityAcceleration[1]*dMassdx[4];
       dGfdx[4][2] = 0.5*gravityAcceleration[2]*dMassdx[4];
       dGfdx[4][3] = dGfdx[4][0]; 
       dGfdx[4][4] = dGfdx[4][1];
       dGfdx[4][5] = dGfdx[4][2];

       dGfdx[5][0] = 0.5*gravityAcceleration[0]*dMassdx[5];
       dGfdx[5][1] = 0.5*gravityAcceleration[1]*dMassdx[5];
       dGfdx[5][2] = 0.5*gravityAcceleration[2]*dMassdx[5];
       dGfdx[5][3] = dGfdx[5][0]; 
       dGfdx[5][4] = dGfdx[5][1];
       dGfdx[5][5] = dGfdx[5][2];
}

void
TwoNodeTruss::getGravityForce(CoordSet& cs, double *gravityAcceleration, Vector &gravityForce,
                              int gravflg, GeomState *geomState)
{
        // Consistent and lumped mass lead to same expression
        const double massPerNode = 0.5 * getMass(cs);

        const double fx = massPerNode * gravityAcceleration[0];
        const double fy = massPerNode * gravityAcceleration[1];
        const double fz = massPerNode * gravityAcceleration[2];

        gravityForce[0] = fx;
        gravityForce[1] = fy;
        gravityForce[2] = fz;
        gravityForce[3] = fx;
        gravityForce[4] = fy;
        gravityForce[5] = fz;
}

FullSquareMatrix
TwoNodeTruss::massMatrix(const CoordSet &cs, double *mel, int cmflg) const
{
        FullSquareMatrix elementMassMatrix(6, mel);
        elementMassMatrix.zero();

        const double mass = getMass(cs);

        if (cmflg) {
                // Consistent mass matrix
                const double outDiagMass = mass / 6.0;
                const double diagMass = 2.0 * outDiagMass;

                for (int i = 0; i < 6; ++i) {
                  elementMassMatrix[i][i] = diagMass;
                }
                for (int i = 0; i < 3; ++i) {
                  const int j = i + 3;
                  elementMassMatrix[i][j] = outDiagMass;
                  elementMassMatrix[j][i] = outDiagMass;
                }
        } else {
                // Lumped mass matrix
                const double massPerNode = 0.5 * mass;
                for (int i = 0; i < 6; ++i) {
                  elementMassMatrix[i][i] = massPerNode;
                }
        }

        return elementMassMatrix;
}

void 
TwoNodeTruss::getStiffnessNodalCoordinateSensitivity(FullSquareMatrix *&dStiffdx, CoordSet &cs)
{
#ifdef USE_EIGEN3
        auto &nd1 = cs.getNode( nn[0] );
        auto &nd2 = cs.getNode( nn[1] );

        Eigen::Array<double,3,1> dconst;
        Eigen::Array<int,0,1> iconst;
        Eigen::Matrix<double,6,1> q;

        dconst[0] = prop->E;
        dconst[1] = prop->A;
        dconst[2] = preload;

        q << nd1.x, nd1.y, nd1.z, nd2.x, nd2.y, nd2.z;
        
        Eigen::Array<Eigen::Matrix<double,6,6>,1,6> dStiffnessdx; 
  
        Simo::FirstPartialSpaceDerivatives<double, TwoNodeTrussStiffnessWRTNodalCoordinateSensitivity> dKdx(dconst,iconst); 
        dStiffnessdx = dKdx(q, 0);
        for(int i=0; i<6; ++i) dStiffdx[i].copy(dStiffnessdx[i].data());
#endif
}

FullSquareMatrix
TwoNodeTruss::stiffness(const CoordSet &cs, double *k, int flg) const
{
        using std::sqrt;

        auto &nd1 = cs.getNode( nn[0] );
        auto &nd2 = cs.getNode( nn[1] );

        double x[2], y[2], z[2];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;

        double dx = x[1] - x[0];
        double dy = y[1] - y[0];
        double dz = z[1] - z[0];

        double length = sqrt( dx*dx + dy*dy + dz*dz );

        if(length == 0.0) {
          fprintf(stderr,"ERROR: truss has zero length %d %d\n",nn[0],nn[1]);
          fprintf(stderr," ... exiting fem program ...\n");
          exit(-1);
        }

        double c1[3];
        
        c1[0] = dx/length;
        c1[1] = dy/length;
        c1[2] = dz/length;

        FullSquareMatrix ret(6,k);

        // Check for negative or zero area and zero modulus
        if(prop->A <= 0.0)
          fprintf(stderr,"ERROR: truss has zero area %d %d\n",nn[0],nn[1]);
        if(prop->E <= 0.0)
          fprintf(stderr,"ERROR: truss has zero modulus %d %d\n",nn[0],nn[1]);

        double elementK = prop->E*prop->A/length;

        int i,j;
        for(i=0; i < 3; ++i)
          for(j=0; j < 3; ++j) {
             ret[i][j]     = elementK*c1[i]*c1[j];
             ret[i+3][j+3] = elementK*c1[i]*c1[j];
             ret[i+3][j]   = -ret[i][j];
             ret[i][j+3]   = -ret[i][j];
          }

        if (preload != 0.0)  {
          for(i=0; i < 3; ++i) {
             ret[i][i]     += preload/length;
             ret[i+3][i+3] += preload/length;
             ret[i+3][i]   = -ret[i][i];
             ret[i][i+3]   = -ret[i][i];
          }
        }
        return ret;
}

int
TwoNodeTruss::numNodes() const
{
        return 2;
}

int *
TwoNodeTruss::nodes(int *p) const
{
        if(p == 0) p = new int[2];
        p[0] = nn[0];
        p[1] = nn[1];
        return p;
}

int
TwoNodeTruss::numDofs() const
{
        return 6;
}

int *
TwoNodeTruss::dofs(DofSetArray &dsa, int *p) const
{
        if(p == 0) p = new int[6];

        dsa.number(nn[0],DofSet::XYZdisp, p  );
        dsa.number(nn[1],DofSet::XYZdisp, p+3);

        return p;
}

void
TwoNodeTruss::markDofs(DofSetArray& dsa) const
{
        dsa.mark(nn, 2, DofSet::XYZdisp);
}

// For Corotational Non Linear problems

Corotator *
TwoNodeTruss::getCorotator(CoordSet &cs, double *kel, int, int)
{
  return new BarCorotator(nn[0], nn[1], prop, preload, cs);
}

int
TwoNodeTruss::getTopNumber() const
{
  return 101;
}

void
TwoNodeTruss::getThermalForce(CoordSet &cs, Vector &ndTemps, Vector &elementThermalForce,
                              int glflag, GeomState *)
{
  // Computes the thermal-mechanical coupling force C*theta
  // A = cross section, W = dilatation coeff
    using std::sqrt;
  double x[2], y[2], z[2], elC[6][2];
  double dx,dy,dz,length;
  int i, j;
 
  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);

  x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
  x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;	
     
  dx = x[1] - x[0];
  dy = y[1] - y[0];
  dz = z[1] - z[0];

  length = sqrt( dx*dx + dy*dy + dz*dz ); 

  if(glflag == 1) { // compute force in local frame
    dx = length;
    dy = 0;
    dz = 0;
  }
     
  double Tref  = prop->Ta;
  double coeff = prop->E*prop->W*prop->A/(2*length);

  //  Coupling matrix:

  elC[0][0] = -coeff*dx;
  elC[0][1] = -coeff*dx;
  elC[1][0] = -coeff*dy;
  elC[1][1] = -coeff*dy;
  elC[2][0] = -coeff*dz;
  elC[2][1] = -coeff*dz;

  elC[3][0] = coeff*dx;
  elC[3][1] = coeff*dx;
  elC[4][0] = coeff*dy;
  elC[4][1] = coeff*dy;
  elC[5][0] = coeff*dz;
  elC[5][1] = coeff*dz;

  for(i=0; i<6; ++i) {
    elementThermalForce[i] = 0.0;
    for(j=0; j<2; ++j)
      elementThermalForce[i] += elC[i][j]*(ndTemps[j]-Tref);
  }
}

void
TwoNodeTruss::getVonMises(Vector& stress, Vector& weight, CoordSet& cs,
                          Vector& elDisp, int strInd, int surface, 
                          double *ndTemps, double ylayer, double zlayer, int avgnum)
{
#ifndef SALINAS
   using std::abs;
   using std::sqrt;
 
   weight = 1.0;
   if(strInd == -1) return;

   auto &nd1 = cs.getNode( nn[0] );
   auto &nd2 = cs.getNode( nn[1] );

   double dx = nd2.x - nd1.x;
   double dy = nd2.y - nd1.y;
   double dz = nd2.z - nd1.z;

   double length = sqrt(dx*dx + dy*dy + dz*dz);

   // scale dx, dy, and dz by the length
   dx /= length;
   dy /= length;
   dz /= length;

   // Compute the change in length of the element
   double dq = dx*(elDisp[3]-elDisp[0])
             + dy*(elDisp[4]-elDisp[1])
             + dz*(elDisp[5]-elDisp[2]);

   // Compute axial strain
   double exx = dq/length;


   switch (avgnum) {

      case 0:
      {
        if (strInd == 0 || strInd == 6) {
           // Compute axial force
           double f = prop->A*prop->E*exx;
	
           // Add Preload
           f  += preload;

           // Compute thermal force
           double coefficient = prop->E*prop->A*prop->W;
           double Tref = prop->Ta;
 
           double fth1(0); 
           double fth2(0); 
           if(ndTemps) {
              fth1 = coefficient*(ndTemps[0]-Tref);
              fth2 = coefficient*(ndTemps[1]-Tref);
           }

           // compute stresses 
           double elForce[2]={0.0,0.0};
           elForce[0] = -f + fth1;
           elForce[1] =  f - fth2;
           if(strInd == 0) {    
              stress[0] = -elForce[0]/prop->A;
              stress[1] =  elForce[1]/prop->A;
           } else {
              stress[0] = std::abs(elForce[0]);
              stress[1] = std::abs(elForce[1]);
           }
        }
        else if (strInd == 7) {
           stress[0] =  exx;
           stress[1] =  exx;
        }
        else {
           stress[0] = 0.0;
           stress[1] = 0.0;
        }
        break;
      }

      case 1: case 3:
      { 
        double xl[3][3];
        double xg[3][3] = {{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
        buildBarFrame(cs, xg, xl);

        if (strInd < 6) {
           // Compute axial force
           double f = prop->A*prop->E*exx;

           // Add Preload
           f  += preload;

           // Compute thermal force
           double coefficient = prop->E*prop->A*prop->W;
           double Tref = prop->Ta;
           double fth1(0);
           double fth2(0);
           if(ndTemps) {
             fth1 = coefficient*(ndTemps[0]-Tref);
             fth2 = coefficient*(ndTemps[1]-Tref);
           }

           // return stresses
           double elForce[2]={0.0,0.0};
           elForce[0] = -f + fth1;
           elForce[1] =  f - fth2;
           double tmpStr1[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
           double tmpStr2[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
           tmpStr1[0] = -elForce[0]/prop->A;
           tmpStr2[0] =  elForce[1]/prop->A;
          _FORTRAN(transform)(xl[0], xl[1], xl[2], xg[0], xg[1], xg[2], tmpStr1);
          _FORTRAN(transform)(xl[0], xl[1], xl[2], xg[0], xg[1], xg[2], tmpStr2);
           stress[0] = tmpStr1[strInd];
           stress[1] = tmpStr2[strInd];
        } else if (strInd == 6) {
           // Compute von Mises stress resultant
           double f = prop->A*prop->E*exx;
           f += preload;

           // Compute thermal force
           double coefficient = prop->E*prop->A*prop->W;
           double Tref = prop->Ta;
           double fth1(0);
           double fth2(0);
           if(ndTemps) {
             fth1 = coefficient*(ndTemps[0]-Tref);
             fth2 = coefficient*(ndTemps[1]-Tref);
           }
          stress[0] = std::abs(-f + fth1);
          stress[1] = std::abs( f - fth1);      
        } else if (strInd > 6 && strInd < 13) {
           double tmpStr[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
           tmpStr[0] = exx;
           _FORTRAN(transform)(xl[0], xl[1], xl[2], xg[0], xg[1], xg[2], tmpStr);
           stress[0] = tmpStr[strInd-7];
           stress[1] = tmpStr[strInd-7];
        }
        else {
           stress[0] = 0.0;
           stress[1] = 0.0;
        }
        break;
      }

      case 2:
      {
        weight = 0.0;
        stress[0] = 0.0;
        stress[1] = 0.0;
        break;
      }

      default:
        std::cerr << "avgnum = " << avgnum << " is not a valid number\n";
    }
#endif
}

void      
TwoNodeTruss::getVonMisesNodalCoordinateSensitivity(GenFullM<double> &dStdx, Vector &weight, CoordSet &cs, Vector &elDisp, int strInd, int surface,
                                                    double *ndTemps, int avgnum, double ylayer, double zlayer)
{ 
  using std::sqrt;

  if(strInd != 6) {
    std::cerr << " ... Error: strInd must be 6 in TwoNodeTruss::getVonMisesNodalCoordinateSensitivity\n";
    exit(-1);
  }
  if(dStdx.numRow() != 6 || dStdx.numCol() != 2) {
    std::cerr << " ... Error: dimension of sensitivity matrix is wrong\n";
    exit(-1);
  }
  weight = 1.0;

  auto &nd1 = cs.getNode(nn[0]);
  auto &nd2 = cs.getNode(nn[1]);

  double dx = nd2.x - nd1.x;
  double dy = nd2.y - nd1.y;
  double dz = nd2.z - nd1.z;
 
  double lengthsquare = dx*dx + dy*dy + dz*dz;
  double length = sqrt(dx*dx + dy*dy + dz*dz);

  // Compute the change in length of the element
  double dq = dx*(elDisp[3]-elDisp[0])
            + dy*(elDisp[4]-elDisp[1])
            + dz*(elDisp[5]-elDisp[2]);

  switch (avgnum) {

    case 1: case 3:
      { 
        if (strInd == 6) {
          double AE = prop->A*prop->E;
          dStdx[0][0] = AE/lengthsquare*(elDisp[3]-elDisp[0]) - 2*AE/(lengthsquare*lengthsquare)*dx*dq;
          dStdx[1][0] = AE/lengthsquare*(elDisp[4]-elDisp[1]) - 2*AE/(lengthsquare*lengthsquare)*dy*dq;
          dStdx[2][0] = AE/lengthsquare*(elDisp[5]-elDisp[2]) - 2*AE/(lengthsquare*lengthsquare)*dz*dq;
          dStdx[3][0] = -AE/lengthsquare*(elDisp[3]-elDisp[0]) + 2*AE/(lengthsquare*lengthsquare)*dx*dq;
          dStdx[4][0] = -AE/lengthsquare*(elDisp[4]-elDisp[1]) + 2*AE/(lengthsquare*lengthsquare)*dy*dq;
          dStdx[5][0] = -AE/lengthsquare*(elDisp[5]-elDisp[2]) + 2*AE/(lengthsquare*lengthsquare)*dz*dq;
          dStdx[0][1] = dStdx[0][0];
          dStdx[1][1] = dStdx[1][0];
          dStdx[2][1] = dStdx[2][0];
          dStdx[3][1] = dStdx[3][0];
          dStdx[4][1] = dStdx[4][0];
          dStdx[5][1] = dStdx[5][0];
          
          // scale dx, dy, and dz by the length
          dx /= length;
          dy /= length;
          dz /= length;

          // Compute the change in length of the element
          dq = dx*(elDisp[3]-elDisp[0])
             + dy*(elDisp[4]-elDisp[1])
             + dz*(elDisp[5]-elDisp[2]);

          // Compute axial strain
          double exx = dq/length;
          // Compute von Mises stress resultant
          double f = prop->A*prop->E*exx;
          f += preload;

          // Compute thermal force
          double coefficient = prop->E*prop->A*prop->W;
          double Tref = prop->Ta;
          double fth1(0);
          double fth2(0);
          if(ndTemps) {
             fth1 = coefficient*(ndTemps[0]-Tref);
             fth2 = coefficient*(ndTemps[1]-Tref);
          }
   
          if(-f+fth1<0) dStdx *= -1; 

        } 
        break;
      }

    case 2:
      {
        weight = 0.0;
        dStdx.zero();
        break;
      }

    default:
      std::cerr << "avgnum = " << avgnum << " is not a valid number\n";
  }
}

void
TwoNodeTruss::getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, GenFullM<double> *dDispDisp,
                                                 CoordSet &cs, Vector &elDisp, int strInd, int surface,
                                                 double *ndTemps, int avgnum, double ylayer, double zlayer)
{
  using std::sqrt;
  using std::abs;

  if(strInd != 6) {
    std::cerr << " ... Error: strInd must be 6 in TwoNodeTruss::getVonMisesDisplacementSensitivity\n";
    exit(-1);
  }
  if(dStdDisp.numRow() != 6 || dStdDisp.numCol() != 2) {
    std::cerr << " ... Error: dimension of sensitivity matrix is wrong\n";
    exit(-1);
  }
  weight = 1.0;

  auto &nd1 = cs.getNode( nn[0] );
  auto &nd2 = cs.getNode( nn[1] );

  double dx = nd2.x - nd1.x;
  double dy = nd2.y - nd1.y;
  double dz = nd2.z - nd1.z;

  double length = sqrt(dx*dx + dy*dy + dz*dz);

  // scale dx, dy, and dz by the length
  dx /= length;
  dy /= length;
  dz /= length;

  // Compute the change in length of the element
  double dq = dx*(elDisp[3]-elDisp[0])
            + dy*(elDisp[4]-elDisp[1])
            + dz*(elDisp[5]-elDisp[2]);

  // Compute axial strain
  double exx = dq/length;

  switch (avgnum) {

    case 0:
    case 1:
    case 3:
      {
        // Compute axial force
        double f = prop->A*prop->E*exx;
	
        // Add Preload
        f  += preload;

        // Compute thermal force
        double coefficient = prop->E*prop->A*prop->W;
        double Tref = prop->Ta;
        double fth1, fth2;

        if(ndTemps) {
          fth1 = coefficient*(ndTemps[0]-Tref);
          fth2 = coefficient*(ndTemps[1]-Tref);
        } else {
          fth1 = 0.0;
          fth2 = 0.0;
        }

        // compute stresses
        double stress[2]; 
        double elForce[2]={0.0,0.0};
        elForce[0] = -f + fth1;
        elForce[1] =  f - fth2;
        stress[0] = std::abs(elForce[0]);
        stress[1] = std::abs(elForce[1]);
        double f1s1 = elForce[0]/stress[0];
        double f2s2 = elForce[1]/stress[1];
   
        // replace automatic differentiation routine with analytic one
        dStdDisp[0][0] =  f1s1*dx;   dStdDisp[1][0] =  f1s1*dy;   dStdDisp[2][0] =  f1s1*dz;   
        dStdDisp[3][0] = -f1s1*dx;   dStdDisp[4][0] = -f1s1*dy;   dStdDisp[5][0] = -f1s1*dz; 
        dStdDisp[0][1] = -f2s2*dx;   dStdDisp[1][1] = -f2s2*dy;   dStdDisp[2][1] = -f2s2*dz;   
        dStdDisp[3][1] =  f2s2*dx;   dStdDisp[4][1] =  f2s2*dy;   dStdDisp[5][1] =  f2s2*dz; 
         
        dStdDisp *= (prop->A*prop->E/length);
        break;
      }

    case 2:
      {
        weight = 0.0;
        dStdDisp.zero();
        break;
      }

    default:
      std::cerr << "avgnum = " << avgnum << " is not a valid number\n";
  }
  if(dDispDisp) dStdDisp ^= (*dDispDisp);
}

void
TwoNodeTruss::buildBarFrame(CoordSet& cs, double xg[3][3], double xl[3][3])
{

        auto &nd1 = cs.getNode( nn[0] );
        auto &nd2 = cs.getNode( nn[1] );

        xl[0][0] = nd2.x - nd1.x;
        xl[0][1] = nd2.y - nd1.y;
        xl[0][2] = nd2.z - nd1.z;

        // local x-axis
        normalize( xl[0] );

        // local y-axis
        double tmp1[3],tmp2[3];

        crossprod(xg[1], xl[0], tmp1); // cross product of global y-axis with local x-axis
        crossprod(xg[2], xl[0], tmp2); // cross product of global z-axis with local x-axis

        if(magnitude(tmp1) > magnitude(tmp2)) {
           xl[1][0] = tmp1[0];
           xl[1][1] = tmp1[1];
           xl[1][2] = tmp1[2];
        } else {
           xl[1][0] = tmp2[0];
           xl[1][1] = tmp2[1];
           xl[1][2] = tmp2[2];
        }
        normalize( xl[1]);

        // local z-axis
        crossprod(xl[0], xl[1], xl[2]);
        normalize( xl[2] );

}

