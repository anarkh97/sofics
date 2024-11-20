******************************************************************
C     H8SHPE computes the value of the shape functions for a
C     eight-noded isoparametric hexahedron and its
C     x-y-z derivatives, at a sample  point given by its hexahedron
C     coordinates (xi,eta,emu)
C
C     Input arguments:
C
C       XI, ETA,EMU   Hexahedral  coordinates of given point
C       X         (8 x 1) array of x coordinates of hexahedron corners
C       Y         (8 x 1) array of y coordinates of hexahedron corners
C       Z         (8 x 1) array of z coordinates of hexahedron corners
C
C     Outputs arguments:
C
C       S         (8 x 1) array of shape function values
C       SX        (8 x 1) array of shape function x-derivatives
C       SY        (8 x 1) array of shape function y-derivatives
C       SZ        (8 x 1) array of shape function z-derivatives 
C       DET        Value of Jacobian determinant
C
C
      subroutine    h8shpe(xi,eta,emu,x,y,z,s,sx,sy,sz,det)
C
C                   A R G U M E N T S
C
      real*8  xi, eta, emu, x(8), y(8), z(8), s(8), sx(8)  
      real*8  sy(8),sz(8),det
C
C                   L O C A L   V A R I A B L E S
C
      integer            i
      real*8   d1, d2, d3, d4, d5, d6, d1h, d2h, d3h, d4h
      real*8   d5h, d6h, d5d6,d2d6,d5d3,d2d3
      real*8   s1(8), s2(8), s3(8), xd1, yd1, zd1, xd2
      real*8   yd2, zd2
      real*8   xd3, yd3, zd3
      real*8   a11, a12, a13, a21, a22, a23, a31, a32, a33
      real*8   cdet
C
C                   L O G I C
C
      d1 =     0.5 * (1.0 + xi)
      d2 =     0.5 * (1.0 + eta)
      d3 =     0.5 * (1.0 + emu)
      d4 =     1.0 - d1
      d5 =     1.0 - d2
      d6 =     1.0 - d3
      d5d6 = d5 * d6
      d2d6 = d2 * d6
      d5d3 = d5 * d3
      d2d3 = d2 * d3
      s(1) =   d4 * d5d6 
      s(2) =   d1 * d5d6
      s(3) =   d1 * d2d6
      s(4) =   d4 * d2d6
      s(5) =   d4 * d5d3
      s(6) =   d1 * d5d3
      s(7) =   d1 * d2d3
      s(8) =   d4 * d2d3
      d1h  =   0.5 * d1
      d2h  =   0.5 * d2
      d3h  =   0.5 * d3
      d4h  =   0.5 * d4
      d5h  =   0.5 * d5
      d6h  =   0.5 * d6

      s1(2) =    d5h * d6
      s1(1) =  - s1(2)
      s1(3) =    d2h * d6
      s1(4) =  - s1(3) 
      s1(6) =    d5h * d3
      s1(5) =  - s1(6)
      s1(7) =    d2h * d3
      s1(8) =  - s1(7) 

      s2(3) =    d1h * d6
      s2(4) =    d4h * d6
      s2(1) =  - s2(4)
      s2(2) =  - s2(3)

      s2(7) =    d1h * d3
      s2(8) =    d4h * d3
      s2(5) =  - s2(8)
      s2(6) =  - s2(7)

      s3(5) =    d4h * d5
      s3(6) =    d1h * d5
      s3(7) =    d1h * d2
      s3(8) =    d4h * d2
      s3(1) =  -s3(5) 
      s3(2) =  -s3(6) 
      s3(3) =  -s3(7) 
      s3(4) =  -s3(8) 
      xd1   =    0.0
      xd2   =    0.0
      xd3   =    0.0
      yd1   =    0.0
      yd2   =    0.0
      yd3   =    0.0
      zd1   =    0.0
      zd2   =    0.0
      zd3   =    0.0
        do   1000   i =  1 , 8
          xd1  =  xd1 +  x(i) * s1(i)
          xd2  =  xd2 +  x(i) * s2(i)
          xd3  =  xd3 +  x(i) * s3(i)
          yd1  =  yd1 +  y(i) * s1(i)
          yd2  =  yd2 +  y(i) * s2(i)
          yd3  =  yd3 +  y(i) * s3(i)
          zd1  =  zd1 +  z(i) * s1(i)
          zd2  =  zd2 +  z(i) * s2(i)
          zd3  =  zd3 +  z(i) * s3(i)
 1000  continue
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
         sx(i) =  0.0
         sy(i) =  0.0
         sz(i) =  0.0
         return
       end if
       cdet  =  1.0 / abs(det)
C      cdet  =  1.0 / det
       do   2000  i =  1 , 8
         sx(i) =  cdet * ( a11 * s1(i) + a12 * s2(i) + a13 * s3(i))
         sy(i) =  cdet * ( a21 * s1(i) + a22 * s2(i) + a23 * s3(i))
         sz(i) =  cdet * ( a31 * s1(i) + a32 * s2(i) + a33 * s3(i))
 2000  continue
       return
       end
