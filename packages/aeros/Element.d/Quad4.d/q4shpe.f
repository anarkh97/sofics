C
C     QUAD4SHAPE computes the value of the shape functions for a
C     four-noded isoparametric quadrilateral and its
C     x-y derivatives, at a sample  point given by its quadrilateral
C     coordinates (xi,eta)
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       call    q4shpe (xi, eta, x, y, s, sx, sy, det)
C
C     Input arguments:
C
C       XI, ETA   Quadrilateral coordinates of given point
C       X         (4 x 1) array of x coordinates of quadrilateral corners
C       Y         (4 x 1) array of y coordinates of quadrilateral corners
C
C     Outputs arguments:
C
C       S         (4 x 1) array of shape function values
C       SX        (4 x 1) array of shape function x-derivatives
C       SY        (4 x 1) array of shape function y-derivatives
C       DET        Value of Jacobian determinant
C
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine    q4shpe ( xi, eta, x, y, s, sx, sy, det)
C
C                   A R G U M E N T S
C
      double precision  xi, eta, x(*), y(*), s(*), sx(*), sy(*), det
C
C                   L O C A L   V A R I A B L E S
C
      integer            i
      double precision   d1, d2, d3, d4, d1h, d2h, d3h, d4h, cdet
      double precision   s1(4), s2(4), xd1, yd1, xd2, yd2
C
C                   L O G I C
C
C
C.... COMPUTE THE SHAPE FUNCTIONS
C
      d1 =     0.5 * (1.0+xi)
      d2 =     0.5 * (1.0+eta)
      d3 =     1.0 - d1
      d4 =     1.0 - d2
      s(1) =   d3 * d4
      s(2) =   d4 * d1
      s(3) =   d1 * d2
      s(4) =   d2 * d3
C
C.... COMPUTE THE SHAPE FUNCTION DERIVATIVES
C
      d1h =    0.5 * d1
      d2h =    0.5 * d2
      d3h =    0.5 * d3
      d4h =    0.5 * d4
      s1(1) =  -d4h
      s1(2) =   d4h
      s1(3) =   d2h
      s1(4) =  -d2h
      s2(1) =  -d3h
      s2(2) =  -d1h
      s2(3) =   d1h
      s2(4) =   d3h
      xd1 =     (x(2)-x(1))*d4h - (x(4)-x(3))*d2h
      yd1 =     (y(2)-y(1))*d4h - (y(4)-y(3))*d2h
      xd2 =     (x(3)-x(2))*d1h - (x(1)-x(4))*d3h
      yd2 =     (y(3)-y(2))*d1h - (y(1)-y(4))*d3h
C
C.... COMPUTE THE DETERMINANT OF THE JACOBIAN
C
      det =      xd1*yd2 - yd1*xd2
      if (det .eq. 0.0)       then
        return
      end if
      cdet =    1.0 /det
      do 2000  i = 1,4
        sx(i) =   cdet * ( yd2*s1(i) - yd1*s2(i))
        sy(i) =   cdet * (-xd2*s1(i) + xd1*s2(i))
 2000   continue
      return
      end
C=END FORTRAN
