#ifdef USE_EIGEN3
#ifndef _TRIANGLE3ELEMENTTEMPLATE_CPP_
#define _TRIANGLE3ELEMENTTEMPLATE_CPP_

#include <cmath>
#include <iostream>
#include <Element.d/Triangle3.d/Triangle3ElementTemplate.hpp>
#include <Eigen/Dense>
#include <Eigen/Geometry>

/*
C=AUTHOR  Y. CHOI, Jan. 2014 
C     Given the node displacements, SANDS4 computes corner
C     point stresses on a four-node quadrilateral in plane stress
C
C
C     The calling sequence is
C
C       CALL   SANDS4 (X, Y, V, STRESS, STRAIN, E, nu)
C
C     where the input arguments are
C
C       X         (3 x 1) array of x coordinates of quadrilateral nodes
C       Y         (3 x 1) array of y coordinates of quadrilateral nodes
C       V         (6 x 1) array of element node displacements arranged
C                    ux1,uy1, ux2,uy2 ...  uy3
C       E         Young's Modulus
C       nu        Poisson's ratio
C
C     The outputs are:
C
C       STRESS    (NUMEL x 3 x 3) array of corner node stresses arranged
C                  sigxx1,sigyy1,tauxy1, sigxx2,sigyy2,tauxy2, ... tauxy3
C       STRAIN    (NUMEL x 3 x 3) array of corner node strains arranged just like
C                  the stresses
C
C=END USAGE
*/

template<typename doublereal>
void
Triangle3ElementTemplate<doublereal>
::sands4(doublereal *_x, doublereal *_y, doublereal *_v, 
         doublereal *_stress, doublereal *_strain, 
         doublereal E, doublereal nu)
{
      Eigen::Map<Eigen::Matrix<doublereal,3,1> > x(_x), y(_y);
      Eigen::Map<Eigen::Matrix<doublereal,6,1> > v(_v);
      Eigen::Map<Eigen::Matrix<doublereal,Eigen::Dynamic,Eigen::Dynamic> > stress(_stress,3,7);
      Eigen::Map<Eigen::Matrix<doublereal,Eigen::Dynamic,Eigen::Dynamic> > strain(_strain,3,7);
      Eigen::Matrix<doublereal, 3, 3> c;

      doublereal area2 = ((x[1]*y[2]-x[2]*y[1])+
                          (x[2]*y[0]-x[0]*y[2])+
                          (x[0]*y[1]-x[1]*y[0]));

      doublereal x21 = x[1] - x[0];
      doublereal x32 = x[2] - x[1];
      doublereal x13 = x[0] - x[2];

      doublereal y12 = y[0] - y[1];
      doublereal y23 = y[1] - y[2];
      doublereal y31 = y[2] - y[0];

      doublereal ux1 = v[0];
      doublereal uy1 = v[1];
      doublereal ux2 = v[2];
      doublereal uy2 = v[3];
      doublereal ux3 = v[4];
      doublereal uy3 = v[5];

      doublereal coef = 1.0/area2;
      doublereal exx  = coef*(y23*ux1 + y31*ux2 + y12*ux3);
      doublereal eyy  = coef*(x32*uy1 + x13*uy2 + x21*uy3);
      doublereal exy  = coef*(x32*ux1 + y23*uy1 + x13*ux2 +
                            y31*uy2 + x21*ux3 + y12*uy3);

      strain(0,0) = exx;
      strain(1,0) = exx;
      strain(2,0) = exx;

      strain(0,1) = eyy;
      strain(1,1) = eyy;
      strain(2,1) = eyy;

      strain(0,2) = 0.0;
      strain(1,2) = 0.0;
      strain(2,2) = 0.0;

      strain(0,3) = exy;
      strain(1,3) = exy;
      strain(2,3) = exy;

      strain(0,4) = 0.0;
      strain(1,4) = 0.0;
      strain(2,4) = 0.0;

      strain(0,5) = 0.0;
      strain(1,5) = 0.0;
      strain(2,5) = 0.0;

//.... Constitutive matrix
      c(0,0) =  E / (1.0-nu*nu);
      c(1,1) =  c(0,0);
      c(2,2) =  0.5*c(0,0)*(1.0-nu);
      c(0,1) =  c(0,0)*nu;
      c(1,0) =  c(0,1);
      c(0,2) =  0.0;
      c(1,2) =  0.0;
      c(2,0) =  0.0;
      c(2,1) =  0.0;

      doublereal sxx = c(0,0)*exx + c(0,1)*eyy + c(0,2)*exy;
      doublereal syy = c(1,0)*exx + c(1,1)*eyy + c(1,2)*exy;
      doublereal sxy = c(2,0)*exx + c(2,1)*eyy + c(2,2)*exy;

      stress(0,0) = sxx;
      stress(1,0) = sxx;
      stress(2,0) = sxx;

      stress(0,1) = syy;
      stress(1,1) = syy;
      stress(2,1) = syy;

      stress(0,2) = 0.0;
      stress(1,2) = 0.0;
      stress(2,2) = 0.0;

      stress(0,3) = sxy;
      stress(1,3) = sxy;
      stress(2,3) = sxy;

      stress(0,4) = 0.0;
      stress(1,4) = 0.0;
      stress(2,4) = 0.0;

      stress(0,5) = 0.0;
      stress(1,5) = 0.0;
      stress(2,5) = 0.0;

//.... COMPUTE THE VON MISES STRESS IN THE PLATE

      doublereal dsxy = stress(0,0) - stress(0,1);
      doublereal dsyz = stress(0,1) - stress(0,2);
      doublereal dsxz = stress(0,0) - stress(0,2);

      doublereal von = sqrt( 0.5*(dsxy*dsxy + dsyz*dsyz + dsxz*dsxz) +
                               3*(stress(0,3)*stress(0,3) + stress(0,4)*stress(0,4) + stress(0,5)*stress(0,5)) );

      stress(0,6) = von;
      stress(1,6) = von;
      stress(2,6) = von;

//.... COMPUTE THE VON MISES STRAIN IN THE PLATE

      dsxy = stress(0,0) - stress(0,1);
      dsyz = stress(0,1) - stress(0,2);
      dsxz = stress(0,0) - stress(0,2);

      von = sqrt( 0.5*(dsxy*dsxy + dsyz*dsyz + dsxz*dsxz) +
                    3*(strain(0,3)*strain(0,3) + strain(0,4)*strain(0,4) + strain(0,5)*strain(0,5)) );

      strain(0,6) = von;
      strain(1,6) = von;
      strain(2,6) = von;
}

template<typename doublereal>
void
Triangle3ElementTemplate<doublereal>
::vms4WRTdisp(doublereal *_x, doublereal *_y, doublereal *_v, 
              doublereal *_vmsWRTdisp,  
              doublereal E, doublereal nu)
{
      Eigen::Map<Eigen::Matrix<doublereal,3,1> > x(_x), y(_y);
      Eigen::Map<Eigen::Matrix<doublereal,6,1> > v(_v);
      Eigen::Map<Eigen::Matrix<doublereal,Eigen::Dynamic,Eigen::Dynamic> > vmsWRTdisp(_vmsWRTdisp,3,6);
      Eigen::Matrix<doublereal,3,7> stress;
      Eigen::Matrix<doublereal,3,3> c;

      Eigen::Matrix<doublereal,3,6> depsilondv;
      Eigen::Matrix<doublereal,3,3> dsigmadepsilon;
      Eigen::Matrix<doublereal,3,3> dvmsdstress;

      depsilondv.setZero();
      dsigmadepsilon.setZero();
      dvmsdstress.setZero();

      doublereal area2 = ((x[1]*y[2]-x[2]*y[1])+
                          (x[2]*y[0]-x[0]*y[2])+
                          (x[0]*y[1]-x[1]*y[0]));

      doublereal x21 = x[1] - x[0];
      doublereal x32 = x[2] - x[1];
      doublereal x13 = x[0] - x[2];

      doublereal y12 = y[0] - y[1];
      doublereal y23 = y[1] - y[2];
      doublereal y31 = y[2] - y[0];

      doublereal ux1 = v[0];
      doublereal uy1 = v[1];
      doublereal ux2 = v[2];
      doublereal uy2 = v[3];
      doublereal ux3 = v[4];
      doublereal uy3 = v[5];

      doublereal coef = 1.0/area2;
      doublereal exx  = coef*(y23*ux1 + y31*ux2 + y12*ux3);
      doublereal eyy  = coef*(x32*uy1 + x13*uy2 + x21*uy3);
      doublereal exy  = coef*(x32*ux1 + y23*uy1 + x13*ux2 +
                            y31*uy2 + x21*ux3 + y12*uy3);

      depsilondv(0,0) = y23;   depsilondv(0,2) = y31;   depsilondv(0,4) = y12;
      depsilondv(1,1) = x32;   depsilondv(1,3) = x13;   depsilondv(1,5) = x21;
      depsilondv(2,0) = x32;   depsilondv(2,1) = y23;   depsilondv(2,2) = x13;
      depsilondv(2,3) = y31;   depsilondv(2,4) = x21;   depsilondv(2,5) = y12;
     

//.... Constitutive matrix
      c(0,0) =  E / (1.0-nu*nu);
      c(1,1) =  c(0,0);
      c(2,2) =  0.5*c(0,0)*(1.0-nu);
      c(0,1) =  c(0,0)*nu;
      c(1,0) =  c(0,1);
      c(0,2) =  0.0;
      c(1,2) =  0.0;
      c(2,0) =  0.0;
      c(2,1) =  0.0;

      doublereal sxx = c(0,0)*exx + c(0,1)*eyy + c(0,2)*exy;
      doublereal syy = c(1,0)*exx + c(1,1)*eyy + c(1,2)*exy;
      doublereal sxy = c(2,0)*exx + c(2,1)*eyy + c(2,2)*exy;

      dsigmadepsilon = c;

      stress.setZero();
      stress(0,0) = sxx;
      stress(1,0) = sxx;
      stress(2,0) = sxx;

      stress(0,1) = syy;
      stress(1,1) = syy;
      stress(2,1) = syy;

      stress(0,3) = sxy;
      stress(1,3) = sxy;
      stress(2,3) = sxy;

//.... COMPUTE THE VON MISES STRESS IN THE PLATE

      doublereal dsxy = stress(0,0) - stress(0,1);
      doublereal dsyz = stress(0,1) - stress(0,2);
      doublereal dsxz = stress(0,0) - stress(0,2);

      doublereal von = sqrt( 0.5*(dsxy*dsxy + dsyz*dsyz + dsxz*dsxz) +
                               3*(stress(0,3)*stress(0,3) + stress(0,4)*stress(0,4) + stress(0,5)*stress(0,5)) );

      stress(0,6) = von;
      stress(1,6) = von;
      stress(2,6) = von;

      dvmsdstress(2,0) = dvmsdstress(1,0) = dvmsdstress(0,0) = (2.*stress(0,0)-stress(0,1))/(2.*stress(0,6));
      dvmsdstress(2,1) = dvmsdstress(1,1) = dvmsdstress(0,1) = (2.*stress(0,1)-stress(0,0))/(2.*stress(0,6));
      dvmsdstress(2,2) = dvmsdstress(1,2) = dvmsdstress(0,2) = (3.*stress(0,3))/stress(0,6);

      vmsWRTdisp = dvmsdstress * dsigmadepsilon * depsilondv;
}

#endif
#endif
