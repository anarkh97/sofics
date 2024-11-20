C=DECK QUAD8SHAPE
C=PURPOSE Compute shape functions and x-y derivatives of 8-node quad
C=AUTHOR C. A. Felippa, April 1969
C=VERSION May 1982 (Fortran 77)
C=EQUIPMENT Machine independent
C=KEYWORDS four node quadrilateral shape functions
C=BLOCK ABSTRACT
C
C     QUAD8SHAPE computes the value of the shape functions for a
C     eight-noded isoparametric quadrilateral and its
C     x-y derivatives, at a sample point given by its quadrilateral
C     coordinates (xi,eta)
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       CALL      QUAD8SHAPE (HF, XI, ETA, X, Y, S, SX, SY, DET)
C
C     Input arguments:
C
C       HF        Option flag, presently ignored
C       XI, ETA   Quadrilateral coordinates of given point
C       X         (8 x 1) array of x coordinates of quadrilateral corners
C       Y         (8 x 1) of y coordinates of quadrilateral corners
C
C     Outputs arguments:
C
C       S         (8 x 1) array of shape function values
C       SX        (8 x 1) array of shape function x-derivatives
C       SY        (8 x 1) array of shape function y-derivatives
C       DET        Value of Jacobian determinant
C
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine    QUAD8SHAPE
     $             (hf, xi, eta, x, y, s, sx, sy, det)
C
C                   A R G U M E N T S
C
      character*1       hf
      double precision  xi, eta, x(8), y(8), s(8), sx(8), sy(8), det
C
C                   L O C A L   V A R I A B L E S
C
      integer            i
      double precision   d1, d2, d3, d4, d1h, d2h, d3h, d4h, cdet
      double precision   d13, d24
      double precision   s1(8), s2(8), xd1, yd1, xd2, yd2
C
C                   L O G I C
C
      d1 =       0.5 * (1.0+xi)
      d2 =       0.5 * (1.0+eta)
      d3 =       1.0 - d1
      d4 =       1.0 - d2
      d1h =      0.5 * d1
      d2h =      0.5 * d2
      d3h =      0.5 * d3
      d4h =      0.5 * d4
      s(1) =     d3*d4*(-xi-eta-1.0)
      s(2) =     d4*d1*(+xi-eta-1.0)
      s(3) =     d1*d2*(+xi+eta-1.0)
      s(4) =     d2*d3*(-xi+eta-1.0)
      s1(1) =   -d4h*(-xi-eta-1.0) - d3*d4
      s1(2) =    d4h*(+xi-eta-1.0) + d4*d1
      s1(3) =    d2h*(+xi+eta-1.0) + d1*d2
      s1(4) =   -d2h*(-xi+eta-1.0) - d2*d3
      s2(1) =   -d3h*(-xi-eta-1.0) - d3*d4
      s2(2) =   -d1h*(+xi-eta-1.0) - d4*d1
      s2(3) =    d1h*(+xi+eta-1.0) + d1*d2
      s2(4) =    d3h*(-xi+eta-1.0) + d2*d3
      d13 =      d1 * d3
      d24 =      d2 * d4
      s(5) =     4.0*d1*d3*d4
      s(6) =     4.0*d2*d4*d1
      s(7) =     4.0*d3*d1*d2
      s(8) =     4.0*d4*d2*d3
      s1(5) =   -2.0*xi*d4
      s1(6) =    2.0*d24
      s1(7) =   -2.0*xi*d2
      s1(8) =   -2.0*d24
      s2(5) =   -2.0*d13
      s2(6) =   -2.0*eta*d1
      s2(7) =    2.0*d13
      s2(8) =   -2.0*eta*d3
      xd1 =     0.0
      yd1 =     0.0
      xd2 =     0.0
      yd2 =     0.0
      do 1500  i = 1,8
        xd1 =     xd1 + x(i)*s1(i)
        yd1 =     yd1 + y(i)*s1(i)
        xd2 =     xd2 + x(i)*s2(i)
        yd2 =     yd2 + y(i)*s2(i)
 1500   continue
      det =      xd1*yd2 - yd1*xd2
      if (det .eq. 0.0)       then
        return
      end if
      cdet =    1.0 /det
      do 2000  i = 1,8
        sx(i) =   cdet * ( yd2*s1(i) - yd1*s2(i))
        sy(i) =   cdet * (-xd2*s1(i) + xd1*s2(i))
 2000   continue
      return
      end
C=END FORTRAN
