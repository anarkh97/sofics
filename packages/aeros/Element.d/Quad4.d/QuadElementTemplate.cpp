#ifdef USE_EIGEN3
#ifndef _QUADELEMENTTEMPLATE_CPP_
#define _QUADELEMENTTEMPLATE_CPP_

#include <cmath>
#include <iostream>
#include <Element.d/Quad4.d/QuadElementTemplate.hpp>
#include <Eigen/Dense>
#include <Eigen/Geometry>

extern int verboseFlag;
/*
C=AUTHOR C. A. Felippa, May 1967
C=REVISED P. R. STERN, MARCH 1990
C=REVISED K. H. PIERSON, MARCH 1997
C=TRANSFORMED to C++ Y. CHOI, Jan. 2014 
C     Given the node displacements, SANDS2 computes corner
C     point stresses on a four-node quadrilateral in plane stress
C
C
C     The calling sequence is
C
C       CALL   SANDS2 (ESCM, X, Y, C, V, STRESS, STRAIN, MAXGUS, MAXSTR, ELM, NUMEL, VMF, TC, TREF, NDTEMPS)
C
C     where the input arguments are
C
C       ESCM      Character string defining stress computation method:
C                 DIRECT  direct stress evaluation at corners
C                 EXTRAP  extrapolate from 2 x 2 Gauss points
C       X         (4 x 1) array of x coordinates of quadrilateral nodes
C       Y         (4 x 1) array of y coordinates of quadrilateral nodes
C       C         (3 x 3) stress-strain constitutive matrix
C       V         (8 x 1) array of element node displacements arranged
C                    ux1,uy1, ux2,uy2 ...  uy4
C       ELM       Element Number
C       VMFLAG    VonMises Stress Flag
C       STRAINFLG VonMises Strain Flag
C       TC        TC=E*alpha/(1-nu)
C       TREF      Reference Temperature
C       NDTEMPS   The temperatures of the quad
C
C
C     The outputs are:
C
C       STRESS    (NUMEL x 3 x 4) array of corner node stresses arranged
C                  sigxx1,sigyy1,tauxy1, sigxx2,sigyy2,tauxy2, ... tauxy4
C       STRAIN    (NUMEL x 3 x 4) array of corner node strains arranged just like
C                  the stresses
C
C=END USAGE
C=BLOCK FORTRAN
*/

template<typename doublereal>
void
QuadElementTemplate<doublereal>
::sands2(char *escm, doublereal *_x, doublereal *_y, doublereal *_c, doublereal *_v, 
         doublereal *_stress, doublereal *_strain, 
         int maxgus, int maxstr, int elm, int msize, bool vmflg, 
         bool strainFlg, doublereal tc, doublereal tref, doublereal *_ndtemps)
{
      Eigen::Map<Eigen::Matrix<doublereal,4,1> > x(_x), y(_y);
      Eigen::Map<Eigen::Matrix<doublereal,3,3> > c(_c);
      Eigen::Map<Eigen::Matrix<doublereal,8,1> > v(_v);
      Eigen::Map<Eigen::Matrix<doublereal,4,1> > ndtemps(_ndtemps);
      Eigen::Map<Eigen::Matrix<doublereal,Eigen::Dynamic,Eigen::Dynamic> > stress(_stress,maxstr,maxgus);
      Eigen::Map<Eigen::Matrix<doublereal,Eigen::Dynamic,Eigen::Dynamic> > strain(_strain,maxstr,maxgus);

      Eigen::Matrix<doublereal,4,1> q, qx, qy, xinod, etanod, tl;
      Eigen::Matrix<doublereal,3,1> sigauss;
      Eigen::Matrix<doublereal,4,4> cext;
      doublereal xi, eta, det, epsxx, epsyy, gamxy, tgp, eptxo, eptyo;
      tl.setZero();      

      xinod << -1.0, 1.0, 1.0,-1.0;
      etanod << -1.0,-1.0, 1.0, 1.0;
      cext << 1.866025404,        -0.5, 0.133974596,        -0.5,
                     -0.5, 1.866025404,        -0.5, 0.133974596,
              0.133974596,        -0.5, 1.866025404,        -0.5,
                     -0.5, 0.133974596,        -0.5, 1.866025404;

//.... COMPUTE THE THERMAL FIELD

      if(_ndtemps) for (int i=0; i<4; ++i) tl[i] = ndtemps[i] - tref;

//                   L O G I C

      if (escm[1] == 'D') { 
        for (int n = 0; n<4; ++n) {
          xi  = xinod[n];
          eta = etanod[n];
          q4shpe(xi, eta, x.data(), y.data(), q.data(), qx.data(), qy.data(), det);
          if (det <= 0.0) {
            std::cerr << " ... Error: Negative Jacobian determinant\n";
            if (det == 0.0) {
              std::cerr << " ... Error: Zero Jacobian determinant\n";
            }
            exit(-1); 
          }

//.... COMPUTE THE THERMAL STRESS

          tgp   =  q[0]*tl[0]+q[1]*tl[1]+q[2]*tl[2]+q[3]*tl[3];
          eptxo = tc*tgp;
          eptyo = tc*tgp;

//.... COMPUTE THE TOTAL STRAIN FROM THE COUPLED SOLVE FOR DISPLACEMENTS

          epsxx = qx[0]*v[0] + qx[1]*v[2] + qx[2]*v[4] + qx[3]*v[6];
          epsyy = qy[0]*v[1] + qy[1]*v[3] + qy[2]*v[5] + qy[3]*v[7];
          gamxy = qy[0]*v[0] + qy[1]*v[2] + qy[2]*v[4] + qy[3]*v[6]
                + qx[0]*v[1] + qx[1]*v[3] + qx[2]*v[5] + qx[3]*v[7];


//.... COMPUTE THE TOTAL STRESS

          stress(0,n) = c(0,0)*epsxx + c(0,1)*epsyy + c(0,2)*gamxy - eptxo;
          strain(0,n) = epsxx; 
          stress(1,n) = c(1,0)*epsxx + c(1,1)*epsyy + c(1,2)*gamxy - eptyo;
          strain(1,n) = epsyy;
          stress(3,n) = c(2,0)*epsxx + c(2,1)*epsyy + c(2,2)*gamxy;
          strain(3,n) = gamxy;
        }

      } else { //.... EXTRAPOLATE FROM THE GAUSS POINTS
        for(int n =0; n<4; ++n ) {
          stress(0,n) = 0.0;
          stress(1,n) = 0.0;
          stress(2,n) = 0.0;
          stress(3,n) = 0.0;
          strain(0,n) = 0.0;   
          strain(1,n) = 0.0;  
          strain(2,n) = 0.0;    
          strain(3,n) = 0.0;    
        }
        for(int i = 0; i<4; ++i) {
          xi =     xinod[i]*0.577350269;
          eta =    etanod[i]*0.577350269;
          q4shpe(xi, eta, x.data(), y.data(), q.data(), qx.data(), qy.data(), det);
          if (det <= 0.0) {
            std::cerr << " ... Error: Negative Jacobian determinant, " << det << "\n";
            if (det == 0.0) {
              std::cerr << " ... Error: Zero Jacobian determinant\n";
            }
            exit(-1); 
          }

//.... COMPUTE THE THERMAL STRESS

          tgp   = q[0]*tl[0]+q[1]*tl[1]+q[2]*tl[2]+q[3]*tl[3];
          eptxo = tc*tgp;
          eptyo = tc*tgp;


//.... COMPUTE THE TOTAL STRAIN FROM THE COUPLED SOLVE FOR DISPLACEMENTS

          epsxx = qx[0]*v[0] + qx[1]*v[2] + qx[2]*v[4] + qx[3]*v[6];
          epsyy = qy[0]*v[1] + qy[1]*v[3] + qy[2]*v[5] + qy[3]*v[7];
          gamxy = qy[0]*v[0] + qy[1]*v[2] + qy[2]*v[4] + qy[3]*v[6]
                + qx[0]*v[1] + qx[1]*v[3] + qx[2]*v[5] + qx[3]*v[7];

//.... COMPUTE THE TOTAL STRESS

          sigauss[0] = c(0,0)*epsxx + c(0,1)*epsyy + c(0,2)*gamxy - eptxo;
          sigauss[1] = c(1,0)*epsxx + c(1,1)*epsyy + c(1,2)*gamxy - eptyo;
          sigauss[2] = c(2,0)*epsxx + c(2,1)*epsyy + c(2,2)*gamxy;

          for(int n=0; n<4; ++n) {
            stress(0,n) =  stress(0,n) + cext(i,n)*sigauss[0];
            strain(0,n) =  epsxx;
            stress(1,n) =  stress(1,n) + cext(i,n)*sigauss[1];
            strain(1,n) =  epsyy; 
            stress(3,n) =  stress(3,n) + cext(i,n)*sigauss[2];
            strain(3,n) =  gamxy; 
          } 
        }
      }

//     set remain stress/strain components to zero

      for(int n = 0; n < 4; ++n) {
        for(int i = 4; i<maxstr; ++i) {
          stress(i,n) = 0.0;
          strain(i,n) = 0.0;    
        }
      }

//.... COMPUTE THE VON MISES STRESS IN THE PLATE

      if (vmflg) {
        vmelmv(stress.data(),maxgus,maxstr,msize,elm,4);
      }

//.... COMPUTE THE VON MISES STRAIN IN THE PLATE

      if (strainFlg) {
        strainvm(strain.data(),maxgus,maxstr,msize,4);
      }

}

template<typename doublereal>
void
QuadElementTemplate<doublereal>
::vms2WRTdisp(char *escm, doublereal *_x, doublereal *_y, doublereal *_c, doublereal *_v, 
              doublereal *_vmsWRTdisp, doublereal *_stressWRTdisp, 
              int maxgus, int maxstr, int elm, int msize, bool vmflg, 
              bool stressFlg, doublereal tc, doublereal tref, doublereal *_ndtemps)
{
      Eigen::Map<Eigen::Matrix<doublereal,4,1> > x(_x), y(_y);
      Eigen::Map<Eigen::Matrix<doublereal,3,3> > c(_c);
      Eigen::Map<Eigen::Matrix<doublereal,8,1> > v(_v);
      Eigen::Map<Eigen::Matrix<doublereal,4,1> > ndtemps(_ndtemps);
      Eigen::Map<Eigen::Matrix<doublereal,Eigen::Dynamic,Eigen::Dynamic> > vmsWRTdisp(_vmsWRTdisp,maxgus,8);

      Eigen::Matrix<doublereal,4,1> q, qx, qy, xinod, etanod, tl;
      Eigen::Matrix<doublereal,3,1> sigauss;
      Eigen::Matrix<doublereal,4,4> cext;
      Eigen::Matrix<doublereal,7,4> stress;
      Eigen::Matrix<doublereal,7,4> strain;
      Eigen::Matrix<doublereal,12,8> depsilondv;
      Eigen::Matrix<doublereal,12,12> dsigaussdepsilon;
      Eigen::Matrix<doublereal,12,12> dstressdsigauss;
      Eigen::Matrix<doublereal,4,12> dvmsdstress;
      Eigen::Map<Eigen::Matrix<doublereal,12,8> > stressWRTdisp(_stressWRTdisp);
   
      doublereal xi, eta, det, epsxx, epsyy, gamxy, tgp, eptxo, eptyo;

//.... INITIALIZE DERIVATIVE VARIABLES
      depsilondv.setZero();
      dsigaussdepsilon.setZero();
      dvmsdstress.setZero(); 
      dstressdsigauss.setZero();
      tl.setZero();

      xinod << -1.0, 1.0, 1.0,-1.0;
      etanod << -1.0,-1.0, 1.0, 1.0;
      cext << 1.866025404,        -0.5, 0.133974596,        -0.5,
                     -0.5, 1.866025404,        -0.5, 0.133974596,
              0.133974596,        -0.5, 1.866025404,        -0.5,
                     -0.5, 0.133974596,        -0.5, 1.866025404;

//.... COMPUTE THE THERMAL FIELD

      if(_ndtemps) for (int i=0; i<4; ++i) tl[i] = ndtemps[i] - tref;

//                   L O G I C

      if (escm[1] == 'D') { 
        for (int n = 0; n<4; ++n) {
          xi  = xinod[n];
          eta = etanod[n];
          q4shpe(xi, eta, x.data(), y.data(), q.data(), qx.data(), qy.data(), det);
          if (det <= 0.0) {
            std::cerr << " ... Error: Negative Jacobian determinant\n";
            if (det == 0.0) {
              std::cerr << " ... Error: Zero Jacobian determinant\n";
            }
            exit(-1); 
          }

//.... COMPUTE THE THERMAL STRESS

          tgp   =  q[0]*tl[0]+q[1]*tl[1]+q[2]*tl[2]+q[3]*tl[3];
          eptxo = tc*tgp;
          eptyo = tc*tgp;

//.... COMPUTE THE TOTAL STRAIN FROM THE COUPLED SOLVE FOR DISPLACEMENTS

          epsxx = qx[0]*v[0] + qx[1]*v[2] + qx[2]*v[4] + qx[3]*v[6];
          epsyy = qy[0]*v[1] + qy[1]*v[3] + qy[2]*v[5] + qy[3]*v[7];
          gamxy = qy[0]*v[0] + qy[1]*v[2] + qy[2]*v[4] + qy[3]*v[6]
                + qx[0]*v[1] + qx[1]*v[3] + qx[2]*v[5] + qx[3]*v[7];


//... COMPUTE THE SENSITIVITY OF TOTAL STRAIN WRT DISPLACEMENTS
 
          depsilondv(3*n,0) = qx[0];   depsilondv(3*n,2) = qx[1];   depsilondv(3*n,4) = qx[2];   depsilondv(3*n,6) = qx[3];
          depsilondv(3*n+1,1) = qy[0]; depsilondv(3*n+1,3) = qy[1]; depsilondv(3*n+1,5) = qy[2]; depsilondv(3*n+1,7) = qy[3];
          depsilondv(3*n+2,0) = qy[0]; depsilondv(3*n+2,1) = qx[0]; depsilondv(3*n+2,2) = qy[1]; depsilondv(3*n+2,3) = qx[1];
          depsilondv(3*n+2,4) = qy[2]; depsilondv(3*n+2,5) = qx[2]; depsilondv(3*n+2,6) = qy[3]; depsilondv(3*n+2,7) = qx[3];


//.... COMPUTE THE TOTAL STRESS

          stress(0,n) = c(0,0)*epsxx + c(0,1)*epsyy + c(0,2)*gamxy - eptxo;
          strain(0,n) = epsxx; 
          stress(1,n) = c(1,0)*epsxx + c(1,1)*epsyy + c(1,2)*gamxy - eptyo;
          strain(1,n) = epsyy;
          stress(3,n) = c(2,0)*epsxx + c(2,1)*epsyy + c(2,2)*gamxy;
          strain(3,n) = gamxy;

//.... COMPUTE THE SENSITIVITY OF TOTAL STRESS WRT DISPLACEMENTS

          dsigaussdepsilon(3*n,3*n) = c(0,0);    dsigaussdepsilon(3*n,3*n+1) = c(0,1);    dsigaussdepsilon(3*n,3*n+2) = c(0,2);
          dsigaussdepsilon(3*n+1,3*n) = c(1,0);  dsigaussdepsilon(3*n+1,3*n+1) = c(1,1);  dsigaussdepsilon(3*n+1,3*n+2) = c(1,2);
          dsigaussdepsilon(3*n+2,3*n) = c(2,0);  dsigaussdepsilon(3*n+2,3*n+1) = c(2,1);  dsigaussdepsilon(3*n+2,3*n+2) = c(2,2);

          dstressdsigauss.Identity();
        }

      } else { //.... EXTRAPOLATE FROM THE GAUSS POINTS
        for(int n =0; n<4; ++n ) {
          stress(0,n) = 0.0;
          stress(1,n) = 0.0;
          stress(2,n) = 0.0;
          stress(3,n) = 0.0;
          strain(0,n) = 0.0;   
          strain(1,n) = 0.0;  
          strain(2,n) = 0.0;    
          strain(3,n) = 0.0;    
        }
        for(int i = 0; i<4; ++i) {
          xi =     xinod[i]*0.577350269;
          eta =    etanod[i]*0.577350269;
          q4shpe(xi, eta, x.data(), y.data(), q.data(), qx.data(), qy.data(), det);
          if (det <= 0.0) {
            std::cerr << " ... Error: Negative Jacobian determinant\n";
            if (det == 0.0) {
              std::cerr << " ... Error: Zero Jacobian determinant\n";
            }
            exit(-1); 
          }

//.... COMPUTE THE THERMAL STRESS

          tgp   = q[0]*tl[0]+q[1]*tl[1]+q[2]*tl[2]+q[3]*tl[3];
          eptxo = tc*tgp;
          eptyo = tc*tgp;


//.... COMPUTE THE TOTAL STRAIN FROM THE COUPLED SOLVE FOR DISPLACEMENTS

          epsxx = qx[0]*v[0] + qx[1]*v[2] + qx[2]*v[4] + qx[3]*v[6];
          epsyy = qy[0]*v[1] + qy[1]*v[3] + qy[2]*v[5] + qy[3]*v[7];
          gamxy = qy[0]*v[0] + qy[1]*v[2] + qy[2]*v[4] + qy[3]*v[6]
                + qx[0]*v[1] + qx[1]*v[3] + qx[2]*v[5] + qx[3]*v[7];

//... COMPUTE THE SENSITIVITY OF TOTAL STRAIN WRT DISPLACEMENTS
 
          depsilondv(3*i,0) = qx[0];   depsilondv(3*i,2) = qx[1];   depsilondv(3*i,4) = qx[2];   depsilondv(3*i,6) = qx[3];
          depsilondv(3*i+1,1) = qy[0]; depsilondv(3*i+1,3) = qy[1]; depsilondv(3*i+1,5) = qy[2]; depsilondv(3*i+1,7) = qy[3];
          depsilondv(3*i+2,0) = qy[0]; depsilondv(3*i+2,1) = qx[0]; depsilondv(3*i+2,2) = qy[1]; depsilondv(3*i+2,3) = qx[1];
          depsilondv(3*i+2,4) = qy[2]; depsilondv(3*i+2,5) = qx[2]; depsilondv(3*i+2,6) = qy[3]; depsilondv(3*i+2,7) = qx[3];

//.... COMPUTE THE TOTAL STRESS AND ITS SENSITIVITY WRT DISPLACEMENTS

          sigauss[0] = c(0,0)*epsxx + c(0,1)*epsyy + c(0,2)*gamxy - eptxo;
          sigauss[1] = c(1,0)*epsxx + c(1,1)*epsyy + c(1,2)*gamxy - eptyo;
          sigauss[2] = c(2,0)*epsxx + c(2,1)*epsyy + c(2,2)*gamxy;

          dsigaussdepsilon(3*i,3*i) = c(0,0);    dsigaussdepsilon(3*i,3*i+1) = c(0,1);    dsigaussdepsilon(3*i,3*i+2) = c(0,2);
          dsigaussdepsilon(3*i+1,3*i) = c(1,0);  dsigaussdepsilon(3*i+1,3*i+1) = c(1,1);  dsigaussdepsilon(3*i+1,3*i+2) = c(1,2);
          dsigaussdepsilon(3*i+2,3*i) = c(2,0);  dsigaussdepsilon(3*i+2,3*i+1) = c(2,1);  dsigaussdepsilon(3*i+2,3*i+2) = c(2,2);

          for(int n=0; n<4; ++n) {
            stress(0,n) =  stress(0,n) + cext(i,n)*sigauss[0];
            strain(0,n) =  epsxx;
            stress(1,n) =  stress(1,n) + cext(i,n)*sigauss[1];
            strain(1,n) =  epsyy; 
            stress(3,n) =  stress(3,n) + cext(i,n)*sigauss[2];
            strain(3,n) =  gamxy;
            dstressdsigauss(3*n,3*i) = cext(i,n);    
            dstressdsigauss(3*n+1,3*i+1) = cext(i,n);
            dstressdsigauss(3*n+2,3*i+2) = cext(i,n);  
          }
        }
      }

//     set remain stress/strain components to zero

      for(int n = 0; n < 4; ++n) {
        for(int i = 4; i<maxstr; ++i) {
          stress(i,n) = 0.0;
          strain(i,n) = 0.0;    
        }
      }

//.... COMPUTE THE VON MISES STRESS IN THE PLATE AND ITS SENSITIVITY WRT DISPLACEMENTS

      if (vmflg) {
        vmelmv(stress.data(),maxgus,maxstr,msize,elm,4);

        for(int n=0; n<4; ++n) {
          dvmsdstress(n,3*n) = (2.*stress(0,n)-stress(1,n))/(2.*stress(6,n));    
          dvmsdstress(n,3*n+1) = (2.*stress(1,n)-stress(0,n))/(2.*stress(6,n));    
          dvmsdstress(n,3*n+2) = (3.*stress(3,n))/stress(6,n);   
        }
//        if(verboseFlag) std::cerr << "printing quadstress\n" << stress << std::endl; 
        vmsWRTdisp = dvmsdstress * dstressdsigauss * dsigaussdepsilon * depsilondv;
//        if(verboseFlag) std::cerr << "printing vmsWRTdisp in quadelement template\n" << vmsWRTdisp << std::endl;
      }

      if (stressFlg) {
        stressWRTdisp = (dstressdsigauss * dsigaussdepsilon) * depsilondv;
/*        if(verboseFlag) {
          std::cerr << "printing depsilondv in quadelement template\n" << depsilondv << std::endl;
          std::cerr << "printing dsigaussdepsilon in quadelement template\n" << dsigaussdepsilon << std::endl;
          std::cerr << "printing dstressdsigauss in quadelement template\n" << dstressdsigauss << std::endl;
          std::cerr << "printing dstressdsigauss*dsigaussdepsilon in quadelement template\n" << dstressdsigauss * dsigaussdepsilon << std::endl;
          std::cerr << "printing stressWRTdisp in quadelement template\n" << stressWRTdisp << std::endl;
        }  */
      }
}

template<typename doublereal>
void
QuadElementTemplate<doublereal>
::q4shpe(doublereal xi, doublereal eta, doublereal* x, doublereal* y, 
         doublereal* s, doublereal* sx, doublereal* sy, doublereal &det)
{
      doublereal d1, d2, d3, d4, d1h, d2h, d3h, d4h, cdet, xd1, yd1, xd2, yd2;
      Eigen::Matrix<doublereal, 4, 1> s1, s2;

//                   L O G I C
//.... COMPUTE THE SHAPE FUNCTIONS

      d1 =     0.5 * (1.0+xi);
      d2 =     0.5 * (1.0+eta);
      d3 =     1.0 - d1;
      d4 =     1.0 - d2;
      s[0] =   d3 * d4;
      s[1] =   d4 * d1;
      s[2] =   d1 * d2;
      s[3] =   d2 * d3;

//.... COMPUTE THE SHAPE FUNCTION DERIVATIVES

      d1h =    0.5 * d1;
      d2h =    0.5 * d2;
      d3h =    0.5 * d3;
      d4h =    0.5 * d4;
      s1[0] =  -d4h;
      s1[1] =   d4h;
      s1[2] =   d2h;
      s1[3] =  -d2h;
      s2[0] =  -d3h;
      s2[1] =  -d1h;
      s2[2] =   d1h;
      s2[3] =   d3h;
      xd1 =     (x[1]-x[0])*d4h - (x[3]-x[2])*d2h;
      yd1 =     (y[1]-y[0])*d4h - (y[3]-y[2])*d2h;
      xd2 =     (x[2]-x[1])*d1h - (x[0]-x[3])*d3h;
      yd2 =     (y[2]-y[1])*d1h - (y[0]-y[3])*d3h;

//.... COMPUTE THE DETERMINANT OF THE JACOBIAN

      det =      xd1*yd2 - yd1*xd2;
      if (det == 0.0) return;
      cdet =    1.0 /det;
      for(int i = 0; i<4; ++i) {
        sx[i] =   cdet * ( yd2*s1[i] - yd1*s2[i]);
        sy[i] =   cdet * (-xd2*s1[i] + xd1*s2[i]);
      }
}

template<typename doublereal>
void
QuadElementTemplate<doublereal>
::vmelmv(doublereal *_stress, int maxgus, int maxstr, int msize, int elm, int nno)
{
        Eigen::Map<Eigen::Matrix<doublereal, Eigen::Dynamic, Eigen::Dynamic> > stress(_stress,maxstr,maxgus);

//.... DECLARE ALL LOCAL VARIABLES

        doublereal sxx,syy,szz,sxy,sxz,syz;
        doublereal dsxx,dsyy,dszz,dsxy,dsxz,dsyz;
        doublereal j2,comp,vms;

        vms = 0.0;
//.... LOOP OVER ALL NODES
        for(int n=0; n<nno; ++n) {

          sxx = stress(0,n);
          syy = stress(1,n);
          szz = stress(2,n);
          sxy = stress(3,n);
          syz = stress(4,n);
          sxz = stress(5,n);

//.... COMPUTE THE FIRST DEVEATORIC STRESSES

          comp = (sxx + syy + szz)/3.0;
          dsxx = sxx - comp;
          dsyy = syy - comp;
          dszz = szz - comp;
          dsxy = sxy;
          dsyz = syz;
          dsxz = sxz;

//.... COMPUTE THE SECOND DEVEATORIC STRESS

          j2 = ((dsxx*dsxx)+(dsyy*dsyy)+(dszz*dszz))/2.0 +
                (dsxy*dsxy)+(dsyz*dsyz)+(dsxz*dsxz);

          stress(6,n) = sqrt(3.0 * j2);
        }
}

template<typename doublereal>
void
QuadElementTemplate<doublereal>
::strainvm(doublereal *_strain, int maxgus, int maxstr, int msize, int numnod)
{
        Eigen::Map<Eigen::Matrix<doublereal, Eigen::Dynamic, Eigen::Dynamic> > strain(_strain,maxstr,maxgus);
         
        doublereal exx,eyy,ezz,exy,exz,eyz;
        doublereal dexx,deyy,dezz,dexy,dexz,deyz;
        doublereal j2,comp,vms;

        vms = 0.0;

//.... LOOP OVER ALL NODES

        for(int n=0; n<numnod; ++n) {

          exx = strain(0,n);
          eyy = strain(1,n);
          ezz = strain(2,n);
// convert to strain tensor from engineering strain
          exy = strain(3,n)/2.0;
          eyz = strain(4,n)/2.0;
          exz = strain(5,n)/2.0;


// ... COMPUTE THE MEAN HYDROSTATIC STRAIN

          comp = (exx + eyy + ezz)/3.0;

//.... COMPUTE THE FIRST DEVIATORIC STRAINS

          dexx = exx - comp;
          deyy = eyy - comp;
          dezz = ezz - comp;
          dexy = exy;
          deyz = eyz;
          dexz = exz;

//.... COMPUTE THE SECOND DEVIATORIC STRAIN

          j2 = ((dexx*dexx)+(deyy*deyy)+(dezz*dezz))/2.0 +
                (dexy*dexy)+(deyz*deyz)+(dexz*dexz);

//.... COMPUTE THE VON MISES STRAIN without averaging

          strain(6,n) = sqrt(3.0 * j2);

        }
}

template<typename doublereal>
void
QuadElementTemplate<doublereal>
::getcmt(doublereal rip, doublereal e, doublereal nu, doublereal *_c)
{
  Eigen::Map<Eigen::Matrix<doublereal,3,3> > c(_c);
  doublereal omn, om2n, opn;

  if(rip == 0) {
     c(0,0) =  e / (1.0-nu*nu);
     c(1,1) =  c(0,0);
     c(2,2) =  0.5*c(0,0)*(1.0-nu);
     c(0,1) =  c(0,0)*nu;
     c(1,0) =  c(0,1);
     c(0,2) =  0.0;
     c(1,2) =  0.0;
     c(2,0) =  0.0;
     c(2,1) =  0.0;
  }

  if(rip == 1) {
     omn  = 1.0 - nu;
     om2n = 1.0 - 2.0*nu;
     opn  = 1.0 + nu;

     c(0,0) =  e*omn/opn/om2n;
     c(1,1) =  c(0,0);
     c(2,2) =  e/opn/2.0;
     c(0,1) =  e*nu /opn/om2n;
     c(1,0) =  c(0,1);
     c(0,2) =  0.0;
     c(1,2) =  0.0;
     c(2,0) =  0.0;
     c(2,1) =  0.0;
  }
}

#endif
#endif
