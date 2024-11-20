#include <Element.d/Shear.d/ShearPanel.h>
#include <Element.d/Shear.d/ShearPanelTemplate.cpp>
#include <Element.d/Shear.d/ShearPanelStressWRTDisplacementSensitivity.h>
#include <Element.d/Function.d/SpaceDerivatives.h>
#include <Utils.d/dbg_alloca.h>
#include <Math.d/Vector.h>
#include <Math.d/FullSquareMatrix.h>
#include <Math.d/matrix.h>
#include <Utils.d/dofset.h>
#include <Element.d/State.h>
#include <Utils.d/linkfc.h>
#include <Utils.d/pstress.h>
#include <Corotational.d/utilities.h>
#include <iostream>

extern int verboseFlag;
extern "C" {
void _FORTRAN(shearpanel)(double*, double*, double*, double&, double&,
                          double&, double&, double&, const int&, double*, 
                          const int&, double&);

void _FORTRAN(shearmass) (double*, double*, double*, double&, double&,
                          const int&, double*, const int&, double&,
                          double*, double*, int&, double&, int&);

void _FORTRAN(spstress) (double*, double*, double*, double*, double&, 
                         double&, double&, double&, double*, double*,
                         double&, double&);

void _FORTRAN(qgauss)(int &, int &, int &, int &,
                      double &,  double &, double &);

void _FORTRAN(q4shpe)(double &, double &, double *, double *,
                      double *, double *, double *, double &);
}

ShearPanel::ShearPanel(int* nodenums)
{
	nn[0] = nodenums[0];
	nn[1] = nodenums[1];
	nn[2] = nodenums[2];
	nn[3] = nodenums[3];
}

Element *
ShearPanel::clone()
{
 return new ShearPanel(*this);
}

void
ShearPanel::renum(const int *table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
}

void
ShearPanel::renum(EleRenumMap& table)
{
  nn[0] = table[nn[0]];
  nn[1] = table[nn[1]];
  nn[2] = table[nn[2]];
  nn[3] = table[nn[3]];
}

void
ShearPanel::getVonMises(Vector& stress,Vector& weight,CoordSet &cs, 
                        Vector& elDisp, int strInd,int,double *,
			double ylayer, double zlayer, int avgnum)
{
 	weight = 1.0;
  if(strInd == -1) return;

	auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);
        auto &nd4 = cs.getNode(nn[3]);

        double x[4], y[4], z[4];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
        x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;

       double elStrain[4][6], elStress[4][6];

       double E  = prop->E;
       double nu = prop->nu;
       double G = E/(2.0*(1.0+nu));

       double F1 = prop->Ixx;
       double F2 = prop->Iyy;

       double vmssig, vmseps;
      _FORTRAN(spstress)(x,y,z,elDisp.data(),G,E,F1,F2,
                        (double*)elStress, (double*)elStrain, vmssig, vmseps);

// if strInd <= 6, you are retrieving a stress value:
// if strInd >  6, you are retrieving a strain value:

      if(strInd <= 6) {
        if(strInd == 6) {
          stress[0] = vmssig; 
          stress[1] = vmssig;
          stress[2] = vmssig;
          stress[3] = vmssig;
        } else {
          stress[0] = elStress[0][strInd];
          stress[1] = elStress[1][strInd]; 
          stress[2] = elStress[2][strInd]; 
          stress[3] = elStress[3][strInd]; 
        }
      } else {
        if(strInd == 13) {
          stress[0] = vmseps;
          stress[1] = vmseps;
          stress[2] = vmseps;
          stress[3] = vmseps;
        } else {
          stress[0] = elStrain[0][strInd-7];
          stress[1] = elStrain[1][strInd-7]; 
          stress[2] = elStrain[2][strInd-7]; 
          stress[3] = elStrain[3][strInd-7]; 
        }
      } 
}

void
ShearPanel::getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, GenFullM<double> *dDispDisp,
                                               CoordSet &cs, Vector &elDisp, int strInd, int surface,
                                               double *ndTemps, int avgnum, double ylayer, double zlayer)
{
#ifdef USE_EIGEN3
   if(strInd != 6) {
     std::cerr << " ... Error: strInd must be 6 in ShearPanel::getVonMisesDisplacementSensitivity\n";
     exit(-1);
   }
   if(dStdDisp.numRow() != 4 || dStdDisp.numCol() !=12) {
     std::cerr << " ... Error: dimenstion of sensitivity matrix is wrong\n";
     exit(-1);
   }
   weight = 1.0;

   auto &nd1 = cs.getNode(nn[0]);
   auto &nd2 = cs.getNode(nn[1]);
   auto &nd3 = cs.getNode(nn[2]);
   auto &nd4 = cs.getNode(nn[3]);

   double vmssig;

  // scalar parameters
  Eigen::Array<double,14,1> dconst;

  dconst[0] = nd1.x; dconst[1] = nd2.x; dconst[2] = nd3.x;  dconst[3] = nd4.x;  // x coordinates
  dconst[4] = nd1.y; dconst[5] = nd2.y; dconst[6] = nd3.y;  dconst[7] = nd4.y;  // y coordinates
  dconst[8] = nd1.z; dconst[9] = nd2.z; dconst[10] = nd3.z; dconst[11] = nd4.z; // z coordinates
  dconst[12] = prop->E;     // E
  dconst[13] = prop->E/(2.0*(1.0+prop->nu));   // G

  Eigen::Array<double,4,1> globalx = dconst.segment<4>(0).cast<double>();
  Eigen::Array<double,4,1> globaly = dconst.segment<4>(4).cast<double>();
  Eigen::Array<double,4,1> globalz = dconst.segment<4>(8).cast<double>();
   
  // integer parameters
  Eigen::Array<int,0,1> iconst;
  // inputs
  Eigen::Matrix<double,12,1> q = Eigen::Map<Eigen::Matrix<double,12,1> >(elDisp.data()).segment(0,12); //displacements

  // Jacobian evaluation
  Eigen::Matrix<double,4,12> dStressdDisp;
  dStressdDisp.setZero();
  vmssWRTdisp(globalx.data(), globaly.data(), globalz.data(), elDisp.data(),
              dconst[13], prop->E, dStressdDisp.data(), vmssig);   
  dStdDisp.copy(dStressdDisp.data());
#else
  std::cerr << " ... Error! ShearPanel::getVonMisesDisplacementSensitivity needs Eigen library.\n";
  exit(-1);
#endif
  if(dDispDisp) dStdDisp ^= (*dDispDisp);
}

void
ShearPanel::getAllStress(FullM& stress,Vector& weight,CoordSet &cs, 
                        Vector& elDisp, int strInd,int,double *)
{
 	weight = 1.0;

	auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);
        auto &nd4 = cs.getNode(nn[3]);

        double x[4], y[4], z[4];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
        x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;

       double elStrain[4][6], elStress[4][6];

       double E  = prop->E;
       double nu = prop->nu;
       double G = E/(2.0*(1.0+nu));

       double F1 = prop->Ixx;
       double F2 = prop->Iyy;

       double vmssig, vmseps;
      _FORTRAN(spstress)(x,y,z,elDisp.data(),G,E,F1,F2,
                        (double*)elStress, (double*)elStrain, vmssig,
                         vmseps);

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
          for (i=0; i<4; ++i) {
            for (j=0; j<6; ++j) {
              stress[i][j] = elStrain[i][j];
            }
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
ShearPanel::getMass(const CoordSet& cs) const
{
        auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);
        auto &nd4 = cs.getNode(nn[3]);

        double x[4], y[4], z[4];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
        x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;

        double *gravityAcceleration = 0;
        double *grvfor = 0;
        double totmas = 0.0;

        int grvflg = 0, masflg = 1;

        const int numgauss = 2;
        double volume;
        double *mel= (double *) dbg_alloca(sizeof(double)*12*12);

       _FORTRAN(shearmass)(x,y,z,prop->eh, prop->rho, numgauss, (double*)mel,
                           numDofs(), volume, gravityAcceleration,
                           grvfor, grvflg, totmas, masflg);

        return totmas;
}

double
ShearPanel::getMassThicknessSensitivity(CoordSet& cs)
{
        auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);
        auto &nd4 = cs.getNode(nn[3]);

        double x[4], y[4], z[4];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
        x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;

        double *gravityAcceleration = 0;
        double *grvfor = 0;
        double totmas = 0.0;

        int grvflg = 0, masflg = 1;

        const int numgauss = 2;
        double volume;
        double *mel= (double *) dbg_alloca(sizeof(double)*12*12);

       _FORTRAN(shearmass)(x,y,z,prop->eh, prop->rho, numgauss, (double*)mel,
                           numDofs(), volume, gravityAcceleration,
                           grvfor, grvflg, totmas, masflg);

        return totmas/prop->eh;
}

void
ShearPanel::getGravityForce(CoordSet& cs,double *gravityAcceleration, 
                              Vector& gravityForce, int gravflg, GeomState *geomState)
{
      // Lumped
      if(gravflg != 2) {
        double massPerNode = 0.25*getMass(cs);

        double fx = massPerNode*gravityAcceleration[0];
        double fy = massPerNode*gravityAcceleration[1];
        double fz = massPerNode*gravityAcceleration[2];

        gravityForce[0] = fx;
        gravityForce[1] = fy;
        gravityForce[2] = fz;
        gravityForce[3] = fx;
        gravityForce[4] = fy;
        gravityForce[5] = fz;
        gravityForce[6] = fx;
        gravityForce[7] = fy;
        gravityForce[8] = fz;
        gravityForce[9]  = fx;
        gravityForce[10] = fy;
        gravityForce[11] = fz;
      }
      
      // Consistent
      else {

        int i;
        int numgauss = 2;
        auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);
        auto &nd4 = cs.getNode(nn[3]);

        double x[4], y[4], z[4];
        double xl[4], yl[4];
        double localg[2];
        double T1[3],T2[3],T3[3],V[3];
        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
        x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;
        double rho = prop->rho;
        double h = prop->eh;

        //  Shift Origin to Node #1
        for (i = 0; i < 4; ++i) {
          x[i] -= x[0];
          y[i] -= y[0];
          z[i] -= z[0];
        }
        // Local X-axis (Node 1->2)
        T1[0] = x[1] - x[0];
        T1[1] = y[1] - y[0];
        T1[2] = z[1] - z[0];
        normalize( T1 );
        // Vector 2 from Node 4->2
        T2[0] = x[1] - x[3];
        T2[1] = y[1] - y[3];
        T2[2] = z[1] - z[3];
        normalize( T2 );
        // Vector 3 from Node 1->3
        T3[0] = x[2] - x[0];
        T3[1] = y[2] - y[0];
        T3[2] = z[2] - z[0];
        normalize( T3 );
        // Perpindicular as cross between v2 and v3
        crossprod( T2, T3, V );
        normalize( V );
        // Local Y-axis as cross between X and V
        crossprod( V, T1, T2 );
        normalize( T2);
        // Local Z-axis as cross between X and Y
        crossprod( T1, T2, T3 );
        normalize( T3);

        // Compute Local "In-plane" Coordinates (required by shape routines)
        for (i = 0; i < 4; ++i) {
          xl[i] = 0.0;
          yl[i] = 0.0;
        }
        for (i = 0; i < 4; ++i) {
          xl[i] = (T1[0]*x[i]) + (T1[1]*y[i]) + (T1[2]*z[i]);
          yl[i] = (T2[0]*x[i]) + (T2[1]*y[i]) + (T2[2]*z[i]);
        }

        // Compute Gravity in Local Frame
        for (i = 0; i < 2; ++i)
          localg[i] = 0.0;
        for (i = 0; i < 3; ++i) {
          localg[0] += T1[i]*gravityAcceleration[i];
          localg[1] += T2[i]*gravityAcceleration[i];
        }

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

            _FORTRAN(q4shpe)(xi, eta, xl, yl,
                             shapeFunc, shapeGradX, shapeGradY, detJ);

            for (i = 0; i < 4; ++i)
              lforce[i] += wt*shapeFunc[i]*detJ;
          }
        }
        double lfx, lfy;
        for(i=0; i<4; ++i) {
          lfx = lforce[i]*h*rho*localg[0];
          lfy = lforce[i]*h*rho*localg[1];
          gravityForce[3*i+0] = (T1[0]*lfx) + (T2[0]*lfy);
          gravityForce[3*i+1] = (T1[1]*lfx) + (T2[1]*lfy);
          gravityForce[3*i+2] = (T1[2]*lfx) + (T2[2]*lfy);
        }
      }

}

void
ShearPanel::getGravityForceThicknessSensitivity(CoordSet& cs,double *gravityAcceleration,
                                                Vector& gravityForceSensitivity, int gravflg, GeomState *geomState)
{
      // Lumped
      if(gravflg != 2) {
        double massPerNodePerThick = 0.25*getMass(cs)/prop->eh;

        double fx = massPerNodePerThick*gravityAcceleration[0];
        double fy = massPerNodePerThick*gravityAcceleration[1];
        double fz = massPerNodePerThick*gravityAcceleration[2];

        gravityForceSensitivity[0] = fx;
        gravityForceSensitivity[1] = fy;
        gravityForceSensitivity[2] = fz;
        gravityForceSensitivity[3] = fx;
        gravityForceSensitivity[4] = fy;
        gravityForceSensitivity[5] = fz;
        gravityForceSensitivity[6] = fx;
        gravityForceSensitivity[7] = fy;
        gravityForceSensitivity[8] = fz;
        gravityForceSensitivity[9]  = fx;
        gravityForceSensitivity[10] = fy;
        gravityForceSensitivity[11] = fz;
      }
      
      // Consistent
      else {

        int i;
        int numgauss = 2;
        auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);
        auto &nd4 = cs.getNode(nn[3]);

        double x[4], y[4], z[4];
        double xl[4], yl[4];
        double localg[2];
        double T1[3],T2[3],T3[3],V[3];
        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
        x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;
        double rho = prop->rho;
        double h = prop->eh;

        //  Shift Origin to Node #1
        for (i = 0; i < 4; ++i) {
          x[i] -= x[0];
          y[i] -= y[0];
          z[i] -= z[0];
        }
        // Local X-axis (Node 1->2)
        T1[0] = x[1] - x[0];
        T1[1] = y[1] - y[0];
        T1[2] = z[1] - z[0];
        normalize( T1 );
        // Vector 2 from Node 4->2
        T2[0] = x[1] - x[3];
        T2[1] = y[1] - y[3];
        T2[2] = z[1] - z[3];
        normalize( T2 );
        // Vector 3 from Node 1->3
        T3[0] = x[2] - x[0];
        T3[1] = y[2] - y[0];
        T3[2] = z[2] - z[0];
        normalize( T3 );
        // Perpindicular as cross between v2 and v3
        crossprod( T2, T3, V );
        normalize( V );
        // Local Y-axis as cross between X and V
        crossprod( V, T1, T2 );
        normalize( T2);
        // Local Z-axis as cross between X and Y
        crossprod( T1, T2, T3 );
        normalize( T3);

        // Compute Local "In-plane" Coordinates (required by shape routines)
        for (i = 0; i < 4; ++i) {
          xl[i] = 0.0;
          yl[i] = 0.0;
        }
        for (i = 0; i < 4; ++i) {
          xl[i] = (T1[0]*x[i]) + (T1[1]*y[i]) + (T1[2]*z[i]);
          yl[i] = (T2[0]*x[i]) + (T2[1]*y[i]) + (T2[2]*z[i]);
        }

        // Compute Gravity in Local Frame
        for (i = 0; i < 2; ++i)
          localg[i] = 0.0;
        for (i = 0; i < 3; ++i) {
          localg[0] += T1[i]*gravityAcceleration[i];
          localg[1] += T2[i]*gravityAcceleration[i];
        }

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

            _FORTRAN(q4shpe)(xi, eta, xl, yl,
                             shapeFunc, shapeGradX, shapeGradY, detJ);

            for (i = 0; i < 4; ++i)
              lforce[i] += wt*shapeFunc[i]*detJ;
          }
        }
        double lfx, lfy;
        for(i=0; i<4; ++i) {
          lfx = lforce[i]*rho*localg[0];
          lfy = lforce[i]*rho*localg[1];
          gravityForceSensitivity[3*i+0] = (T1[0]*lfx) + (T2[0]*lfy);
          gravityForceSensitivity[3*i+1] = (T1[1]*lfx) + (T2[1]*lfy);
          gravityForceSensitivity[3*i+2] = (T1[2]*lfx) + (T2[2]*lfy);
        }
      }

}

FullSquareMatrix
ShearPanel::massMatrix(const CoordSet &cs,double *mel,int cmflg) const
{
        auto &nd1 = cs.getNode(nn[0]);
        auto &nd2 = cs.getNode(nn[1]);
        auto &nd3 = cs.getNode(nn[2]);
        auto &nd4 = cs.getNode(nn[3]);

        double x[4], y[4], z[4];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
        x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;

        double *gravityAcceleration = 0;
        double *grvfor = 0;
        double totmas = 0.0;

        int grvflg = 0, masflg = 0;

        const int numgauss = 2;
        double volume;
       _FORTRAN(shearmass)(x,y,z,prop->eh, prop->rho, numgauss, (double*)mel,
                           numDofs(), volume, gravityAcceleration,
                           grvfor, grvflg, totmas, masflg);

        FullSquareMatrix ret(numDofs(),mel);

        return ret;
}

FullSquareMatrix
ShearPanel::stiffness(const CoordSet &cs, double *d, int flg) const
{
	auto &nd1 = cs.getNode(nn[0]);
	auto &nd2 = cs.getNode(nn[1]);
	auto &nd3 = cs.getNode(nn[2]);
	auto &nd4 = cs.getNode(nn[3]);

        double x[4], y[4], z[4];

        x[0] = nd1.x; y[0] = nd1.y; z[0] = nd1.z;
        x[1] = nd2.x; y[1] = nd2.y; z[1] = nd2.z;
        x[2] = nd3.x; y[2] = nd3.y; z[2] = nd3.z;
        x[3] = nd4.x; y[3] = nd4.y; z[3] = nd4.z;

        double E  = prop->E;
        double nu = prop->nu;
        double G = E/(2.0*(1.0+nu));

        double F1 = prop->Ixx;
        double F2 = prop->Iyy;

        const int numgauss = 2;
        double volume;
        _FORTRAN(shearpanel)(x,y,z,prop->eh,G,E,F1,F2,numgauss,(double*)d,
                             numDofs(), volume);

        FullSquareMatrix ret(numDofs(),d);

        return ret;
}

int
ShearPanel::numNodes() const
{
 	return 4;
}

int*
ShearPanel::nodes(int *p) const
{
 	if(p == 0) p = new int[4];
 	p[0] = nn[0];
 	p[1] = nn[1];
 	p[2] = nn[2];
 	p[3] = nn[3];
	return p;
}

int
ShearPanel::numDofs() const
{
 	return 12;
}

int*
ShearPanel::dofs(DofSetArray &dsa, int *p) const
{
 	if(p == 0) p = new int[numDofs()];

	dsa.number(nn[0], DofSet::XYZdisp, p);
	dsa.number(nn[1], DofSet::XYZdisp, p+3);
	dsa.number(nn[2], DofSet::XYZdisp, p+6);
	dsa.number(nn[3], DofSet::XYZdisp, p+9);

	return p;
}

void
ShearPanel::markDofs(DofSetArray &dsa) const
{
        dsa.mark(nn, 4, DofSet::XYZdisp);
}

int
ShearPanel::getTopNumber() const
{
  return 118; // 2
}

