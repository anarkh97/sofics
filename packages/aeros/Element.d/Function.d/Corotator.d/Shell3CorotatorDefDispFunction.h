#ifndef _SHELL3COROTATORDEFDISPFUNCTION_H_
#define _SHELL3COROTATORDEFDISPFUNCTION_H_

#include <Element.d/Function.d/Function.h>
#include <Element.d/Function.d/utilities.hpp>
#include <Element.d/Function.d/AutoDiffScalarPlugin.h>

#include <cmath>
#include <iostream>

namespace Simo {

template<typename Scalar>
class Shell3CorotatorDefDispFunction : public VectorValuedFunction<18,18,Scalar,36,1,double>
{
    Eigen::Matrix<double,3,3> x0;
    Eigen::Array<Eigen::Matrix3d,3,1> R;
    int fitAlg;

  public:
    Shell3CorotatorDefDispFunction(const Eigen::Array<double,36,1>& sconst, const Eigen::Array<int,1,1>& iconst)
    {
      x0 = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(const_cast<double*>(sconst.data())+0);
      R[0] = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(const_cast<double*>(sconst.data())+9);
      R[1] = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(const_cast<double*>(sconst.data())+18);
      R[2] = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(const_cast<double*>(sconst.data())+27);
      using std::abs;
      fitAlg = iconst[0];
    }

    Eigen::Matrix<Scalar,18,1> operator() (const Eigen::Matrix<Scalar,18,1>& q, Scalar)
    {
      // q[0] = x translation of node 1
      // q[1] = y translation of node 1
      // q[2] = z translation of node 1
      // q[3] = 1st component of axis-angle rotation vector of node 1
      // q[4] = 2nd component of axis-angle rotation vector of node 1
      // q[5] = 3rd component of axis-angle rotation vector of node 1
      // q[6] = x translation of node 2
      // q[7] = y translation of node 2
      // q[8] = z translation of node 2
      // q[9] = 1st component of axis-angle rotation vector of node 2
      // q[10] = 2nd component of axis-angle rotation vector of node 2
      // q[11] = 3rd component of axis-angle rotation vector of node 2
      // q[12] = x translation of node 3
      // q[13] = y translation of node 3
      // q[14] = z translation of node 3
      // q[15] = 1st component of axis-angle rotation vector of node 3
      // q[16] = 2nd component of axis-angle rotation vector of node 3
      // q[17] = 3rd component of axis-angle rotation vector of node 3

      // return value: element strain energy

      Eigen::Matrix<Scalar,3,3> xn;
      // Set current configuration in xn
      xn(0,0) = x0(0,0) + q[0]; // x coordinate of node 1
      xn(0,1) = x0(0,1) + q[1]; // y coordinate of node 1
      xn(0,2) = x0(0,2) + q[2]; // z coordinate of node 1

      xn(1,0) = x0(1,0) + q[6]; // x coordinate of node 2
      xn(1,1) = x0(1,1) + q[7]; // y coordinate of node 2
      xn(1,2) = x0(1,2) + q[8]; // z coordinate of node 2

      xn(2,0) = x0(2,0) + q[12]; // x coordinate of node 3
      xn(2,1) = x0(2,1) + q[13]; // y coordinate of node 3
      xn(2,2) = x0(2,2) + q[14]; // z coordinate of node 3

      Eigen::Matrix<double,3,3> t0, xl0;
      Eigen::Matrix<Scalar,3,3> t0n, xln;
      // TODO: don't cast x0 to Scalar
      localCoord(x0, xn, t0, t0n, xl0, xln, fitAlg);
      
      // Form the displacement from C0n to Cn
      Eigen::Matrix<Scalar,18,1> vld;

      // Translation for node 1
      vld[0]  = xln(0,0) - xl0(0,0);
      vld[1]  = xln(0,1) - xl0(0,1);
      vld[2]  = xln(0,2) - xl0(0,2);

      // Translation for node 2
      vld[6]  = xln(1,0) - xl0(1,0);
      vld[7]  = xln(1,1) - xl0(1,1);
      vld[8]  = xln(1,2) - xl0(1,2);

      // Translation for node 3
      vld[12] = xln(2,0) - xl0(2,0);
      vld[13] = xln(2,1) - xl0(2,1);
      vld[14] = xln(2,2) - xl0(2,2);

      // Create rotation part of the deformation vector
      Eigen::Matrix<Scalar,3,3> r0,dr;
      Eigen::Matrix<Scalar,3,1> w,dw;

      for(int inode = 0; inode < 3; ++inode) {

         w = q.template segment<3>(inode*6 + 3);
         vec_to_mat<Scalar>(w, r0);
         dr = t0n*r0*(R[inode]*t0.transpose()).template cast<Scalar>();

         mat_to_vec<Scalar>(dr, dw);
         vld.template segment<3>(inode*6 + 3) = dw;
      }

      // transform element displacement vector from local to global coordinates
      for(int i = 0; i < 6; ++i) {
         vld.template segment<3>(i*3) = t0n.transpose()*vld.template segment<3>(i*3);
      }

      return vld;
    }

  private:
    void localCoord(const Eigen::Matrix<double,3,3> &x0, const Eigen::Matrix<Scalar,3,3> &xn,
                    Eigen::Matrix<double,3,3> &t0, Eigen::Matrix<Scalar,3,3> &t0n,
                    Eigen::Matrix<double,3,3> &xl0, Eigen::Matrix<Scalar,3,3> &xln,
                    const int fitAlg) const
    /*****************************************************************
     *
     *  Purpose:
     *     Form the transformation matrix to local coordinate system
     *     for a 3 node element.
     *
     *  Method:
     *     Local z-axis is computed normal to the element both for t0 and t0n.
     *         Local x-axis is long side 1-2 of the x0 coordinates for t0
     *     transformation matrix and xl0 are the coordinates of the element
     *     in this local coordinate system.
     *         Local x-axis is rotated relative to side 1-2 for the deformed
     *     element given by xn coordinates.
     *     xln are the coordinates of the element in this local coordinate system.
     *         The x-axis is rotated so that xl0 are a best fit of the undeformed
     *     element on top of the deformed element given by xln. fitalg
     *     gives which algorithm to use for this rotation.
     *
     *
     *  Input;
     *     fitalg : == 1 x axis along side 1-2, both for t0 and t0n
     *              == 2 x axis rotated from side 1-2 in so that
     *                   the sum of angles between xl0 and xln side edges
     *                   are zero.
     *              == 3 x axis rotated so that xl0 is rotated relative
     *                   to xln with the continuum mechanics definition
     *                   of rotations
     * (conf = configuration)
     *     x0(i,j) : global coordinate componenet j of node i in undeformed conf.
     *     xn(i,j) : global coordinate componenet j of node i in deformed conf.
     *
     *  Output:
     *     t0       : transformation matrix to local undeformed coordinates
     *     t0n      : transformation matrix to local deformed coordinates
     *     xl0(i,j): local coordinate componenet j of node i in undeformed conf.
     *                also best fit of undef. element on def element.
     *     xln(i,j): local coordinate componenet j of node i in deformed conf.
     *
     *  Coded by: Bjorn Haugen; adjusted for C++ by Teymour Manzouri
     *******************************************************************/
     {
       using std::asin;
       using std::sqrt;
       using std::atan;
       int    i,  nod;
       double l0;
       Scalar alpha, ln;
       Eigen::Matrix<double,3,1> xc0, s0;
       Eigen::Matrix<Scalar,3,1> xcn, sn;
       Eigen::Matrix<Scalar,3,3> rten;
    
       // Compute coordinates of centroid for undeformed and deformed configuration
       for( i=0; i<3; ++i ) {
          xc0[i] = ( x0(0,i) + x0(1,i) + x0(2,i) )/3.0;
          xcn[i] = ( xn(0,i) + xn(1,i) + xn(2,i) )/3.0;
       }
    
       int nodi = 0;
       int nodj = 1;
       int nodk = 2;
    
       // Compute t0 transformation matrix with x axis along side 1-2 
       for( i=0; i<3; i++ ) t0(0,i) = x0(nodj,i) - x0(nodi,i);
       t0.row(0).normalize();
    
       // local y axis
       for( i=0; i<3; i++ ) t0(1,i) = x0(nodk,i) - x0(nodi,i);
       t0.row(2) = t0.row(0).cross(t0.row(1));
       t0.row(2).normalize();
    
       // local z axis
       t0.row(1) = t0.row(2).cross(t0.row(0));
       t0.row(1).normalize();
    
       // Compute local coordinates of undeformed element
       for( nod=0; nod<3; nod++) {
          for( i=0; i<3; i++ ) {
             xl0(nod,i) = t0(i,0)*(x0(nod,0) - xc0[0])
                         +t0(i,1)*(x0(nod,1) - xc0[1])
                         +t0(i,2)*(x0(nod,2) - xc0[2]);
         }
       }
    
       // Compute t0n transformation matrix with x axis along side 1-2
       for( i=0; i<3; i++ ) t0n(0,i) = xn(nodj,i) - xn(nodi,i);
       t0n.row(0).normalize();
    
       for( i=0; i<3; i++ ) t0n(1,i) = xn(nodk,i) - xn(nodi,i);
       t0n.row(2) = t0n.row(0).cross(t0n.row(1));
       t0n.row(2).normalize();
    
       t0n.row(1) = t0n.row(2).cross(t0n.row(0));
       t0n.row(1).normalize();
    
       // Compute local coordinates of deformed element
       for( nod=0; nod<3; nod++) {
          for( i=0; i<3; i++ ) {
             xln(nod,i) = t0n(i,0)*(xn(nod,0) - xcn[0])
                         +t0n(i,1)*(xn(nod,1) - xcn[1])
                         +t0n(i,2)*(xn(nod,2) - xcn[2]);
         }
       }
    
       switch( fitAlg ) {
         case 1 : {       // *********** fitalg == 1 ************
           alpha = Scalar(0.0);
         } break;
    
         default :
         case 2 : {  // *********** fitalg == 2 ************
           // Compute angle between side 2-3 for deformed and undeformed
           Scalar innerVal;
           for( i=0; i<2; i++ ) {
             s0[i] = xl0(2,i) - xl0(1,i);
             sn[i] = xln(2,i) - xln(1,i);
           }
           l0 = sqrt ( s0[0]*s0[0] + s0[1]*s0[1] );
           ln = sqrt ( sn[0]*sn[0] + sn[1]*sn[1] );
           for( i=0; i<2; i++ ) {
             s0[i] /= l0;
             sn[i] /= ln;
           }
           innerVal = s0[0]*sn[0] + s0[1]*sn[1];
           if (innerVal < 0.0){
             std::cerr << "WARNING OBTUSE ANGLE 2-3 In localCoord\n";
             alpha = M_PI - asin(sn[0]*s0[1] - s0[0]*sn[1]);
           }
           else alpha = asin(sn[0]*s0[1] - s0[0]*sn[1]);
    
           // Compute angle between side 3-1 for deformed and undeformed 
           for( i=0; i<2; i++ ) {
             s0[i] = xl0(0,i) -xl0(2,i);
             sn[i] = xln(0,i) -xln(2,i);
           }
           l0 = sqrt ( s0[0]*s0[0] + s0[1]*s0[1] );
           ln = sqrt ( sn[0]*sn[0] + sn[1]*sn[1] );
           for( i=0; i<2; i++ ) {
             s0[i] /= l0;
             sn[i] /= ln;
           }
           innerVal = s0[0]*sn[0] + s0[1]*sn[1];
           if (innerVal < 0.0){
             std::cerr << "WARNING OBTUSE ANGLE 3-1 In localCoord\n";
             alpha += M_PI - asin( sn[0]*s0[1] - s0[0]*sn[1] );
           }
           else alpha += asin( sn[0]*s0[1] - s0[0]*sn[1] );
           alpha = -alpha/3.0;
         } break;
    
         case 3 : {  // *********** fitalg == 3 ************
    
           Scalar u3  = (xln(nodk,0) - xln(nodi,0)) - (xl0(nodk,0) - xl0(nodi,0));
    
           Scalar u3p = ((xln(nodj,0) - xln(nodi,0)) - (xl0(nodj,0) - xl0(nodi,0)))
                        *(xl0(nodk,0) - xl0(nodi,0))/(xln(nodj,0) - xln(nodi,0));
    
           Scalar h3  = (xln(nodk,1) - xln(nodi,1));
    
           alpha = (0.5*atan(-(u3 - u3p)/h3));
         } break;
    
         case 4 : { // minimization of local nodal displacements
    
           Scalar p = Scalar(0.0), q = Scalar(0.0);
           for(int i=0; i<3; ++i) { p += (xln(i,1)*xl0(i,0)-xln(i,0)*xl0(i,1));
                                    q += (xln(i,0)*xl0(i,0)+xln(i,1)*xl0(i,1)); }
           alpha = atan(p/q);
         } break;
    
         case 5 : { // zero local spin at the centroid of the element
    
           int index[5] = { 0,1,2,0,1 };
           Scalar p = Scalar(0.0), q = Scalar(0.0);
           for(int l=0; l<3; ++l) { int i = index[l], j = index[l+1], k = index[l+2]; 
                                    p += (xl0(i,1)-xl0(j,1))*xln(k,1) - (xl0(j,0)-xl0(i,0))*xln(k,0);
                                    q += (xl0(j,0)-xl0(i,0))*xln(k,1) + (xl0(i,1)-xl0(j,1))*xln(k,0); }
           alpha = atan(p/q);
         } break;
       }
    
       // Compute rotation vector and rotation tensor for t0n matrix
       for( i=0; i<3; i++ )
         sn[i] = alpha*t0n(2,i);
  
       vec_to_mat<Scalar>(sn, rten);
    
       // Rotate x-axis angle alpha about z-axis for t0n
       for( i=0; i<3; i++ ) sn[i] = rten(i,0)*t0n(0,0)
                                   +rten(i,1)*t0n(0,1)
                                   +rten(i,2)*t0n(0,2);
       for( i=0; i<3; i++ )
         t0n(0,i) = sn[i];
    
       t0n.row(0).normalize();
 
       // Compute y-axis as cross product of z-axis and x-axis
       t0n.row(1) = t0n.row(2).cross(t0n.row(0));
       t0n.row(1).normalize();
 
       // Recompute local coordinates for deformed configuration
       for( nod=0; nod<3; nod++) {
          for( i=0; i<3; i++ ) {
             xln(nod,i) = t0n(i,0)*(xn(nod,0) - xcn[0])
                         +t0n(i,1)*(xn(nod,1) - xcn[1])
                         +t0n(i,2)*(xn(nod,2) - xcn[2]);
         }
       }
    
    }

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

} // namespace Simo

#endif
