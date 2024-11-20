#ifdef USE_EIGEN3
#ifndef _TETRAELEMENTTEMPLATE_CPP_
#define _TETRAELEMENTTEMPLATE_CPP_

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <Element.d/Tetra.d/TetraElementTemplate.hpp>
#include <Eigen/Dense>
#include <Eigen/Geometry>

template<typename doublereal>
void
TetraElementTemplate<doublereal>
::sands23(int elm, doublereal *_x, doublereal *_y, doublereal *_z,
          doublereal e, doublereal xnu, doublereal *_u, doublereal *_stress, doublereal *_strain, 
          int maxgus, int maxstr, int msize, int outerr, bool vmflg, bool strainFlg)
{
      int ndim(3),npo(4),npi(4),nno(4);

      Eigen::Map<Eigen::Matrix<doublereal,12,1> > u(_u);
      Eigen::Map<Eigen::Matrix<doublereal,7,4> > stress(_stress);
      Eigen::Map<Eigen::Matrix<doublereal,7,4> > strain(_strain);
      Eigen::Map<Eigen::Matrix<doublereal,4,1> > X(_x), Y(_y), Z(_z); 

      Eigen::Matrix<int,12,1> IJT;
      Eigen::Matrix<doublereal,72,4> sigmae, straine;
      Eigen::Matrix<doublereal,4,1> xint,yint,zint; 
      Eigen::Matrix<doublereal,9,1> elas;
      Eigen::Matrix<doublereal,4,1> delta;
      Eigen::Matrix<doublereal,3,3> df;
      Eigen::Matrix<doublereal,12,4> a2;
      Eigen::Matrix<doublereal,4,12> DP;
      Eigen::Matrix<doublereal,4,4> vp1;
      Eigen::Matrix<doublereal,3,3> dfinv;

      DP << -1.,-1.,-1.,1.,0.,0.,0.,1.,0.,0.,0.,1.,
            -1.,-1.,-1.,1.,0.,0.,0.,1.,0.,0.,0.,1.,   
            -1.,-1.,-1.,1.,0.,0.,0.,1.,0.,0.,0.,1.,   
            -1.,-1.,-1.,1.,0.,0.,0.,1.,0.,0.,0.,1.;
      IJT << 0, 3, 6, 9,  1, 4, 7, 10,  2, 5, 8, 11;
      vp1.setIdentity();

      doublereal XLAM,XMU;
      int INDICE(1),I,J,J1,L,N;
      for(L=0; L<npi; ++L) {

         df.setZero();
         for(I = 0; I<ndim; ++I) {
           for(N = 0; N<npo; ++N) {

             df(I,0) += DP(L,ndim*N+I) * X(N,0);
             df(I,1) += DP(L,ndim*N+I) * Y(N,0);
             df(I,2) += DP(L,ndim*N+I) * Z(N,0);
           }
         }
         if (ndim == 3) {
            dfinv(0,0) = df(1,1)*df(2,2)-df(1,2)*df(2,1);
            dfinv(1,1) = df(0,0)*df(2,2)-df(2,0)*df(0,2);
            dfinv(2,2) = df(0,0)*df(1,1)-df(1,0)*df(0,1);
            dfinv(0,1) = df(0,2)*df(2,1)-df(0,1)*df(2,2);
            dfinv(1,0) = df(1,2)*df(2,0)-df(1,0)*df(2,2);
            dfinv(0,2) = df(0,1)*df(1,2)-df(0,2)*df(1,1);
            dfinv(2,0) = df(1,0)*df(2,1)-df(1,1)*df(2,0);
            dfinv(1,2) = df(1,0)*df(0,2)-df(1,2)*df(0,0);
            dfinv(2,1) = df(2,0)*df(0,1)-df(2,1)*df(0,0);

            delta[L] = df(0,0) * dfinv(0,0) + df(1,0) * dfinv(0,1) + df(2,0) * dfinv(0,2);
         } else {
            delta[L] = sqrt( (df(0,0)*df(1,1)-df(0,1)*df(1,0))*(df(0,0)*df(1,1)-df(0,1)*df(1,0))
                     + (df(0,1)*df(1,2)-df(0,2)*df(1,1))*(df(0,1)*df(1,2)-df(0,2)*df(1,1))
                     + (df(0,2)*df(1,0)-df(0,0)*df(1,2))*(df(0,2)*df(1,0)-df(0,0)*df(1,2)));
         }
         if (INDICE > 0 && INDICE <= 2) {

           for(I = 0; I<3; ++I) {
             for(J = 0; J<nno; ++J) {
              a2(3*J+I,L) =  dfinv(I,0) * DP(L,ndim*J)  +
                             dfinv(I,1) * DP(L,ndim*J+1)  +
                             dfinv(I,2) * DP(L,ndim*J+2);
             }
           }
         }

         if (INDICE >= 2) {
              xint[L] = 0.0;
              yint[L] = 0.0;
              zint[L] = 0.0;
           for(I=0; I<npo; ++I) {
              xint[L] = xint[L] + vp1(I,L) * X[I];
              yint[L] = yint[L] + vp1(I,L) * Y[I];
              zint[L] = zint[L] + vp1(I,L) * Z[I];
           }
        }
      }
//      fobase(3,3,nno,npo,npi,vp1.data(),DP.data(),DP.data(),X.data(),Y.data(),Z.data(),
//             a2.data(),xint.data(),yint.data(),zint.data(),delta.data(),dfinv.data(),df.data(),INDICE);

//           --------------
      XMU   = e/(2.*(1.+xnu));
      XLAM  = e*xnu/((1.+xnu)* (1.-2.*xnu));
      elas[0] = XLAM + 2. * XMU;
      elas[2] = XLAM + 2. * XMU;
      elas[5] = XLAM + 2. * XMU;
      elas[1] = XLAM;
      elas[3] = XLAM;
      elas[4] = XLAM;
      elas[6] = XMU;
      elas[7] = XMU;
      elas[8] = XMU;

      for(L=0; L<npi; ++L) {

//        B = delta * (df)-1 * DP(L)

        for(J=0; J<nno; ++J) {

//   --- BLOC1

            J1=IJT[J];
            straine(6*J1,L) = a2(3*J,L) / delta[L];
            straine(6*J1+1,L) = 0.0;
            straine(6*J1+2,L) = 0.0;
//           straine(6*J1+1,L) = a2(3*J,L) / delta[L];
//           straine(6*J1+2,L) = a2(3*J,L) / delta[L];
            straine(6*J1+3,L) = a2(3*J+1,L) / delta[L];
            straine(6*J1+4,L) = 0.0;
            straine(6*J1+5,L) = a2(3*J+2,L) / delta[L];

            sigmae(6*J1,L) = elas[0]*a2(3*J,L) / delta[L];
            sigmae(6*J1+1,L) = elas[1]*a2(3*J,L) / delta[L];
            sigmae(6*J1+2,L) = elas[3]*a2(3*J,L) / delta[L];
            sigmae(6*J1+3,L) = elas[6]*a2(3*J+1,L) / delta[L];
            sigmae(6*J1+4,L) = 0.0;
            sigmae(6*J1+5,L) = elas[8]*a2(3*J+2,L) / delta[L];

//   --- BLOC2

            J1=IJT[J+nno];
            straine(6*J1,L) = 0.0;
//           straine(6*J1,L) = a2(3*J+1,L) / delta[L];
            straine(6*J1+1,L) = a2(3*J+1,L) / delta[L];
            straine(6*J1+2,L) = 0.0;
//           straine(6*J1+2,L) = a2(3*J+1,L) / delta[L];
            straine(6*J1+3,L) = a2(3*J,L) / delta[L];
            straine(6*J1+4,L) = a2(3*J+2,L) / delta[L];
            straine(6*J1+5,L) = 0.0;

            sigmae(6*J1,L) = elas[1]*a2(3*J+1,L) / delta(L,0);
            sigmae(6*J1+1,L) = elas[2]*a2(3*J+1,L) / delta(L,0);
            sigmae(6*J1+2,L) = elas[3]*a2(3*J+1,L) / delta(L,0);
            sigmae(6*J1+3,L) = elas[6]*a2(3*J,L) / delta(L,0);
            sigmae(6*J1+4,L) = elas[7]*a2(3*J+2,L) / delta(L,0);
            sigmae(6*J1+5,L) = 0.0;

//   --- BLOC3

            J1=IJT[J+nno*2];
            straine(6*J1,L) = 0.0;
            straine(6*J1+1,L) = 0.0;
//           straine(6*J1,L) = a2(3*J+2,L) / delta(L,0);
//           straine(6*J1+1,L) = a2(3*J+2,L) / delta(L,0);
            straine(6*J1+2,L) = a2(3*J+2,L) / delta(L,0);
            straine(6*J1+3,L) = 0.0;
            straine(6*J1+4,L) = a2(3*J+1,L) / delta(L,0);
            straine(6*J1+5,L) = a2(3*J,L) / delta(L,0);

            sigmae(6*J1,L) = elas[3]*a2(3*J+2,L) / delta(L,0);
            sigmae(6*J1+1,L) = elas[4]*a2(3*J+2,L) / delta(L,0);
            sigmae(6*J1+2,L) = elas[5]*a2(3*J+2,L) / delta(L,0);
            sigmae(6*J1+3,L) = 0.0;
            sigmae(6*J1+4,L) = elas[7]*a2(3*J+1,L) / delta(L,0);
            sigmae(6*J1+5,L) = elas[8]*a2(3*J,L) / delta(L,0);

        } // J

//       multiplier la contraine elementaire par le deplacement

        for(I=0; I<6; ++I) {

           stress(I,L) = 0.0;
           strain(I,L) = 0.0;
           for(J=0; J<3*nno; ++J) {
              stress(I,L) += sigmae(6*J+I,L)*u(J,0);
              strain(I,L) += straine(6*J+I,L)*u(J,0); 
           }
        }

    }  // L

//  les coordonnees des points de calcul sont xint,yint,zint

//.... COMPUTE THE VON MISES STRESS

        if (vmflg) 
          vmelmv(stress.data(),maxgus,maxstr,msize,elm,nno);


//.... COMPUTE THE VON MISES STRAIN

        if(strainFlg) 
          strainvm(strain.data(),maxgus,maxstr,msize,nno);
}

template<typename doublereal>
void
TetraElementTemplate<doublereal>
::vmelmv(doublereal *_stress, int maxgus, int maxstr, int msize, int elm, int nno)
{
        Eigen::Map<Eigen::Matrix<doublereal, 7, 4> > stress(_stress);

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
TetraElementTemplate<doublereal>
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

#endif
#endif
