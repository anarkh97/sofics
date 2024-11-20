C=DECK TRIG6SHAPE
C=PURPOSE Form shape functions and derivatives of 4-node quad
C=AUTHOR C. A. Felippa, January 1967
C=VERSION May 1982 (Fortran 77)
C=EQUIPMENT Machine independent
C=KEYWORDS six node triangle shape functions
C=BLOCK ABSTRACT
C
C     TRIG6SHAPE computes the value of the shape functions for a
C     six-noded isoparametric triangle and its
C     x-y derivatives, at a sample point given by its triangle
C     coordinates (zeta1, zeta2, zeta3)
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       CALL  TRIG6SHAPE (HF, ZETA1, ZETA2, ZETA3, X, Y, S, SX, SY, DET)
C
C     Input arguments:
C
C       HF         A dummy character argument.
C      ZETA1,ZETA2,ZETA3  Triangular coordinates of given point
C       X         (3 x 1) array of x coordinates of triangle corners
C       Y         (3 x 1) array of y coordinates of triangle corners
C
C     Outputs arguments:
C
C       S         (3 x 1) array of shape function values
C       SX        (3 x 1) array of shape function x-derivatives
C       SY        (3 x 1) array of shape function y-derivatives
C       DET        Value of Jacobian determinant
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine    TRIG6SHAPE
     $       (hf, zeta1, zeta2, zeta3, x, y, s, sx, sy, det)
C
C                   A R G U M E N T S
C
      character*(*)     hf
      double precision  zeta1, zeta2, zeta3, x(6), y(6)
      double precision  s(6), sx(6), sy(6), det
C
C                   L O C A L   V A R I A B L E S
C
      double precision   x21, x32, x13, y12, y23, y31, cdet
      double precision   jx1, jx2, jx3, jy1, jy2, jy3
C
C                   L O G I C
C
      s(1) =     zeta1*(2.0*zeta1-1.0)
      s(2) =     zeta2*(2.0*zeta2-1.0)
      s(3) =     zeta3*(2.0*zeta3-1.0)
      s(4) =     4.0*zeta1*zeta2
      s(5) =     4.0*zeta2*zeta3
      s(6) =     4.0*zeta3*zeta1
      jx1  =     4.0*(x(1)*zeta1+x(4)*zeta2+x(6)*zeta3) - x(1)
      jx2  =     4.0*(x(2)*zeta2+x(5)*zeta3+x(4)*zeta1) - x(2)
      jx3  =     4.0*(x(3)*zeta3+x(6)*zeta1+x(5)*zeta2) - x(3)
      jy1  =     4.0*(y(1)*zeta1+y(4)*zeta2+y(6)*zeta3) - y(1)
      jy2  =     4.0*(y(2)*zeta2+y(5)*zeta3+y(4)*zeta1) - y(2)
      jy3  =     4.0*(y(3)*zeta3+y(6)*zeta1+y(5)*zeta2) - y(3)
      x21 =      jx2 - jx1
      x32 =      jx3 - jx2
      x13 =      jx1 - jx3
      y12 =      jy1 - jy2
      y23 =      jy2 - jy3
      y31 =      jy3 - jy1
      det =      x21*y31 - y12*x13
      if (det .le. 0.0)       then
        return
      end if
      cdet =      1.0 /det
      sx(1) =     cdet * (4.0*zeta1-1.0) * y23
      sx(2) =     cdet * (4.0*zeta2-1.0) * y31
      sx(3) =     cdet * (4.0*zeta3-1.0) * y12
      sy(1) =     cdet * (4.0*zeta1-1.0) * x32
      sy(2) =     cdet * (4.0*zeta2-1.0) * x13
      sy(3) =     cdet * (4.0*zeta3-1.0) * x21
      sx(4) =     cdet * 4.0*(zeta2*y23+zeta1*y31)
      sx(5) =     cdet * 4.0*(zeta3*y31+zeta2*y12)
      sx(6) =     cdet * 4.0*(zeta1*y12+zeta3*y23)
      sy(4) =     cdet * 4.0*(zeta2*x32+zeta1*x13)
      sy(5) =     cdet * 4.0*(zeta3*x13+zeta2*x21)
      sy(6) =     cdet * 4.0*(zeta1*x21+zeta3*x32)
      return
      end
C=END FORTRAN
C=DECK TRIGGAUSSQ
C=PURPOSE Get abscissas and weight factors for Gauss triangle rule
C=AUTHOR C. A. Felippa, April 1967
C=VERSION May 1982 (Fortran 77)
C=EQUIPMENT Machine independent
C=KEYWORDS triangle Gauss integration rule absissae weigth
C=BLOCK ABSTRACT
C
C     TRIGGAUSSQ returns the triangular coordinates of
C     sample points and weights for a Gauss-type integration
C     rule over a triangle. Rules are identified by total
C     number of points.
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       call      TRIGGAUSSQ  (P, I, ZETA1, ZETA2, ZETA3, WEIGHT)
C
C     Input arguments:
C
C       ABS(P)    Total number of Gauss points in rule:
C                 1, 3 or 7.  If none of these, assume 1.
C                 Sign of P is used to discriminate between two
C                 rules with same number of points.  For ABS(P)<8
C                 the only conflict occurs at 3 points:
C                   P=-3  use midpoint (1/2,1/2,0) rule
C                   P=+3  use (2/3,1/6,1/6) rule.
C       I         Index of sample point in  rule.
C
C     Outputs arguments:
C
C       ZETA1,ZETA2,ZETA3   Triangular coordinates of sample point
C       WEIGHT    Weight factor
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine    TRIGGAUSSQ
     $             (p, i, zeta1, zeta2, zeta3, weight)
C
C                   A R G U M E N T S
C
      integer           p, i
      double precision  zeta1, zeta2, zeta3, weight
C
C                   L O C A L   V A R I A B L E S
C
      double precision  zeta(3), sqrt15
      integer           pp
C
C                   L O G I C
C
      pp =   abs(p)
      if (pp .eq. 1)             then
        zeta(1) =   1.0D0/3.
        zeta(2) =   zeta(1)
        zeta(3) =   zeta(2)
        weight =  1.0D0
      else if (pp .eq. 3)        then
        if (p .lt. 0)            then
          zeta(1) =  0.5D0
          zeta(2) =  zeta(1)
          zeta(3) =  zeta(2)
          zeta(i) =  0.0
          weight =   1.0D0/3.
        else
          zeta(1) =  1.D0/6.
          zeta(2) =  zeta(1)
          zeta(3) =  zeta(2)
          zeta(i) =  2.D0/3.
          weight =   1.0D0/3.
        end if
      else if (pp .eq. 7)        then
        sqrt15 =   sqrt(15.D0)
        if (i .eq. 1)            then
          zeta(1) =  1.D0/3.
          zeta(2) =  zeta(1)
          zeta(3) =  zeta(2)
          weight =   9.D0/40.
        else if (i .le. 4)       then
          zeta(1) =  (6.D0-sqrt15)/21.
          zeta(2) =  zeta(1)
          zeta(3) =  zeta(2)
          zeta(i-1) = (9.D0+2.D0*sqrt15)/21.
          weight =   (155.D0-sqrt15)/1200.
        else
          zeta(1) =  (6.D0+sqrt15)/21.
          zeta(2) =  zeta(1)
          zeta(3) =  zeta(2)
          zeta(i-4) = (9.D0-2.D0*sqrt15)/21.
          weight =   (31.D0/120.)-((155.D0-sqrt15)/1200.)
        end if
      else
        zeta(1) =   1.0D0/3.
        zeta(2) =   zeta(1)
        zeta(3) =   zeta(2)
        weight =  1.0D0
      end if
      zeta1 =  zeta(1)
      zeta2 =  zeta(2)
      zeta3 =  zeta(3)
      return
      end
C=END FORTRAN
