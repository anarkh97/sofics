C
C     TRI3SHAPE computes the value of the shape functions for a
C     three-noded isoparametric triangle and its
C     x-y derivatives, at a sample  point given by its parent
C     coordinates (xi,eta)
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       call    q3shpe (xi, eta, x, y, s, sx, sy, det)
C
C     Input arguments:
C
C       XI, ETA   Quadrilateral coordinates of given point
C       X         (3 x 1) array of x coordinates of quadrilateral corners
C       Y         (3 x 1) array of y coordinates of quadrilateral corners
C
C     Outputs arguments:
C
C       S         (3 x 1) array of shape function values
C       SX        (3 x 1) array of shape function x-derivatives
C       SY        (3 x 1) array of shape function y-derivatives
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
      double precision   cdet
      double precision   s1(3), s2(4), xd1, yd1, xd2, yd2
C
C                   L O G I C
C
C
C.... COMPUTE THE SHAPE FUNCTIONS
C
      s(1) =   xi
      s(2) =   eta
      s(3) =   1-xi-eta
C
C.... COMPUTE THE SHAPE FUNCTION DERIVATIVES
C
      xd1 =     (y(2)-y(3))
      yd1 =    -(y(1)-y(3))
      xd2 =    -(x(2)-x(3))
      yd2 =     (x(1)-x(3))
      s1(1) =   1 
      s1(2) =   0
      s1(3) =  -1
      s2(1) =   0
      s2(2) =   1
      s2(3) =  -1
C
C.... COMPUTE THE DETERMINANT OF THE JACOBIAN
C
      det =     (xd1*yd2 - yd1*xd2)/(2*area)
      if (det .eq. 0.0)       then
        return
      end if
      cdet =    1.0 /det
  do 2000  i = 1,4
        sx(i) =   cdet * s1(i)
        sy(i) =   cdet * s2(i)
 2000   continue
      return
      end
C=END FORTRAN
