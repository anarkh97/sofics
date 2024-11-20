C-----------------------------------------------------------------------
C     HEXA20SHAPE computes the value of the shape functions for a
C     twenty-noded isoparametric hexahedron and its
C     x-y-z derivatives, at a sample  point given by its hexahedron
C     coordinates (xi,eta,mu)
C
C     Input arguments:
C
C       XI,ETA,MU  Hexahedral (iso-P) coordinates of given point
C       X         (20 x 1) array of x coordinates of hexahedron corners
C       Y         (20 x 1) array of y coordinates of hexahedron corners
C       Z         (20 x 1) array of z coordinates of hexahedron corners
C
C     Outputs arguments:
C
C       SF        (20 x 1) array of shape function values
C       SX        (20 x 1) array of shape function x-derivatives
C       SY        (20 x 1) array of shape function y-derivatives
C       SZ        (20 x 1) array of shape function z-derivatives 
C       DET        Value of Jacobian determinant
C-----------------------------------------------------------------------
C HEXA20 nodes/shape fcts numbering:
C                     16
C             5+-------+-------+8
C             /|              /|             
C            / |             / |             
C         13+  |          15+  |               
C          / 17+           /   +20             
C         /    | 14       /    |    t          
C       6+-------+-------+7    |    |          
C        |     |      12 |     |    |          
C        |    1+-------+-|-----+ 4  +---- s    
C        |    /          |    /    /           
C      18+   /         19+   /    /           
C        | 9+            |  +11  r       
C        | /             | /            
C        |/              |/                   
C       2+-------+-------+3                  
C               10                  
C-----------------------------------------------------------------------
      subroutine  h20shpe(xi, eta, mu, x, y, z, sf, sx, sy, sz, det)
C
C                   A R G U M E N T S
C
      double precision  xi, eta, mu, x(20), y(20), z(20)
      double precision  sf(20), sx(20), sy(20), sz(20), det
C
C                   L O C A L   V A R I A B L E S
C
      integer            i
      double precision   r, s, t, rp, rm, sp, sm, tp, tm, rrm, ssm, ttm
      double precision   xd1, xd2, xd3, yd1, yd2, yd3, zd1, zd2, zd3
      double precision   a11, a12, a13, a21, a22, a23, a31, a32, a33
      double precision   cdet
      double precision   sd1(20), sd2(20), sd3(20)
C
C-------------------------------------------- FORM BASIC FUNCTION VALUES
      r  =  xi
      s  =  eta
      t  =  mu
      rp =  1.d0 + r
      rm =  1.d0 - r
      sp =  1.d0 + s
      sm =  1.d0 - s
      tp =  1.d0 + t
      tm =  1.d0 - t
      rrm = 1.d0 - r*r
      ssm = 1.d0 - s*s
      ttm = 1.d0 - t*t
C     +-----------------------------------------------------------------
C     I Q U A D R A T I C  SHAPE FUNCTIONS AND THEIR NATURAL DERIVATIVES
C     I20-NODED BRICK ELEMENT 
C     +-----------------------------------------------------------------
      sf(1)  = 0.125d0*rm*sm*tm*(rm+sm+tm-5.d0)
      sf(2)  = 0.125d0*rp*sm*tm*(rp+sm+tm-5.d0)
      sf(3)  = 0.125d0*rp*sp*tm*(rp+sp+tm-5.d0)
      sf(4)  = 0.125d0*rm*sp*tm*(rm+sp+tm-5.d0)
      sf(5)  = 0.125d0*rm*sm*tp*(rm+sm+tp-5.d0)
      sf(6)  = 0.125d0*rp*sm*tp*(rp+sm+tp-5.d0)
      sf(7)  = 0.125d0*rp*sp*tp*(rp+sp+tp-5.d0)
      sf(8)  = 0.125d0*rm*sp*tp*(rm+sp+tp-5.d0)
      sf(9)  = 0.25d0*rrm*sm*tm
      sf(10) = 0.25d0*rp*ssm*tm
      sf(11) = 0.25d0*rrm*sp*tm
      sf(12) = 0.25d0*rm*ssm*tm
      sf(17) = 0.25d0*rm*sm*ttm
      sf(18) = 0.25d0*rp*sm*ttm
      sf(19) = 0.25d0*rp*sp*ttm
      sf(20) = 0.25d0*rm*sp*ttm
      sf(13) = 0.25d0*rrm*sm*tp
      sf(14) = 0.25d0*rp*ssm*tp
      sf(15) = 0.25d0*rrm*sp*tp
      sf(16) = 0.25d0*rm*ssm*tp
C--------------------------------------- DERIVATIVE EVALUATION
      sd1(1)  = -0.125d0*sm*tm*(2.d0*rm+sm+tm-5.d0)
      sd1(2)  =  0.125d0*sm*tm*(2.d0*rp+sm+tm-5.d0)
      sd1(3)  =  0.125d0*sp*tm*(2.d0*rp+sp+tm-5.d0)
      sd1(4)  = -0.125d0*sp*tm*(2.d0*rm+sp+tm-5.d0)
      sd1(5)  = -0.125d0*sm*tp*(2.d0*rm+sm+tp-5.d0)
      sd1(6)  =  0.125d0*sm*tp*(2.d0*rp+sm+tp-5.d0)
      sd1(7)  =  0.125d0*sp*tp*(2.d0*rp+sp+tp-5.d0)
      sd1(8)  = -0.125d0*sp*tp*(2.d0*rm+sp+tp-5.d0)
      sd1(9)  = -0.5d0*r*sm*tm
      sd1(10) =  0.25d0*ssm*tm
      sd1(11) = -0.5d0*r*sp*tm
      sd1(12) = -sd1(10)
      sd1(18) =  0.25d0*sm*ttm
      sd1(17) = -sd1(18)
      sd1(19) =  0.25d0*sp*ttm
      sd1(20) = -sd1(19)
      sd1(13) = -0.5d0*r*sm*tp
      sd1(14) =  0.25d0*ssm*tp
      sd1(15) = -0.5d0*r*sp*tp
      sd1(16) = -sd1(14)
      sd2(1)  = -0.125d0*tm*rm*(2.d0*sm+tm+rm-5.d0)
      sd2(2)  = -0.125d0*tm*rp*(2.d0*sm+tm+rp-5.d0)
      sd2(3)  =  0.125d0*tm*rp*(2.d0*sp+tm+rp-5.d0)
      sd2(4)  =  0.125d0*tm*rm*(2.d0*sp+tm+rm-5.d0)
      sd2(5)  = -0.125d0*tp*rm*(2.d0*sm+tp+rm-5.d0)
      sd2(6)  = -0.125d0*tp*rp*(2.d0*sm+tp+rp-5.d0)
      sd2(7)  =  0.125d0*tp*rp*(2.d0*sp+tp+rp-5.d0)
      sd2(8)  =  0.125d0*tp*rm*(2.d0*sp+tp+rm-5.d0)
      sd2(10) = -0.5d0*s*tm*rp
      sd2(11) =  0.25d0*rrm*tm
      sd2(9)  = -sd2(11)
      sd2(12) = -0.5d0*s*tm*rm
      sd2(18) = -0.25d0*ttm*rp
      sd2(19) = -sd2(18)
      sd2(20) =  0.25d0*ttm*rm
      sd2(17) = -sd2(20)
      sd2(15) =  0.25d0*rrm*tp
      sd2(13) = -sd2(15)
      sd2(14) = -0.5d0*s*tp*rp
      sd2(16) = -0.5d0*s*tp*rm
      sd3(1)  = -0.125d0*rm*sm*(2.d0*tm+rm+sm-5.d0)
      sd3(2)  = -0.125d0*rp*sm*(2.d0*tm+rp+sm-5.d0)
      sd3(3)  = -0.125d0*rp*sp*(2.d0*tm+rp+sp-5.d0)
      sd3(4)  = -0.125d0*rm*sp*(2.d0*tm+rm+sp-5.d0)
      sd3(5)  =  0.125d0*rm*sm*(2.d0*tp+rm+sm-5.d0)
      sd3(6)  =  0.125d0*rp*sm*(2.d0*tp+rp+sm-5.d0)
      sd3(7)  =  0.125d0*rp*sp*(2.d0*tp+rp+sp-5.d0)
      sd3(8)  =  0.125d0*rm*sp*(2.d0*tp+rm+sp-5.d0)
      sd3(9)  = -0.25d0*rrm*sm
      sd3(10) = -0.25d0*ssm*rp
      sd3(11) = -0.25d0*rrm*sp
      sd3(12) = -0.25d0*ssm*rm
      sd3(17) = -0.5d0*t*rm*sm
      sd3(18) = -0.5d0*t*rp*sm
      sd3(19) = -0.5d0*t*rp*sp
      sd3(20) = -0.5d0*t*rm*sp
      sd3(14) = -sd3(10)
      sd3(15) = -sd3(11)
      sd3(16) = -sd3(12)
      sd3(13) = -sd3(9)
      xd1   =    0.0
      xd2   =    0.0
      xd3   =    0.0
      yd1   =    0.0
      yd2   =    0.0
      yd3   =    0.0
      zd1   =    0.0
      zd2   =    0.0
      zd3   =    0.0
      do 1500   i = 1,20
         xd1  =  xd1 +  x(i) * sd1(i)
         xd2  =  xd2 +  x(i) * sd2(i)
         xd3  =  xd3 +  x(i) * sd3(i)
         yd1  =  yd1 +  y(i) * sd1(i)
         yd2  =  yd2 +  y(i) * sd2(i)
         yd3  =  yd3 +  y(i) * sd3(i)
         zd1  =  zd1 +  z(i) * sd1(i)
         zd2  =  zd2 +  z(i) * sd2(i)
         zd3  =  zd3 +  z(i) * sd3(i)
 1500  continue
       a11 = yd2*zd3 - yd3*zd2
       a12 = yd3*zd1 - yd1*zd3
       a13 = yd1*zd2 - yd2*zd1
       a21 = xd3*zd2 - xd2*zd3
       a22 = xd1*zd3 - zd1*xd3
       a23 = xd2*zd1 - xd1*zd2
       a31 = xd2*yd3 - xd3*yd2
       a32 = yd1*xd3 - xd1*yd3
       a33 = xd1*yd2 - yd1*xd2
       det = xd1*a11 + yd1*a21 + zd1*a31
       if ( det .eq. 0.0)  then
            return
       end if
       cdet  =  1.0/det
       do 2000  i = 1,20
          sx(i) =  cdet * ( a11 * sd1(i) + a12 * sd2(i) + a13 * sd3(i))
          sy(i) =  cdet * ( a21 * sd1(i) + a22 * sd2(i) + a23 * sd3(i))
          sz(i) =  cdet * ( a31 * sd1(i) + a32 * sd2(i) + a33 * sd3(i))
 2000    continue

       return
       end
