C=DECK TRIG6STIF1
C=PURPOSE Form stiffness of 6-node plane-stress triangle
C=AUTHOR C. A. Felippa,  February 1967
C=AUTHOR Antonini Puppin-Macedo (modif. for 1dof per node), June 1999
C=VERSION May 1982 (Fortran 77)
C=EQUIPMENT Machine independent
C=KEYWORDS plane stress six node triangle
C=KEYWORDS finite element stiffness matrix
C=BLOCK ABSTRACT
C
C     TRIG6MSTIF forms the element stiffness matrix of a
C     six-node triangle in plane stress.
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       CALL   TRIG6MSTIF (X, Y, P, SM, M)
C
C     where the input arguments are
C
C       X         (3 x 1) array of x coordinates of triangle nodes
C       Y         (3 x 1) array of y coordinates of triangle nodes
C       P         Identifies quadrature rule by value and sign
C                 (see TRIGGAUSSQ)
C       M         First dimension of SM in calling program.
C
C     The outputs are:
C
C       SM        (6 x 6) computed element stiffness matrix.  The
C                 arrangement of rows and columns pertains to node
C                 displacements arranged in the order
C                  (v1, v2, ... v6)
C                 This particular ordering is obtained by setting array
C                 LS to 1,3,5,7, ... 12 (see DATA statement below)
C       STATUS    Status character variable.  Blank if no error detected.
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine    TRIG6STIF1
     $          (x, y, p, sm, m)
C
C                   A R G U M E N T S
C
      integer           p, m
      double precision  x(6), y(6)
      double precision  sm(m,m)
C
C                   L O C A L   V A R I A B L E S
C
      double precision  q(6), qx(6), qy(6)
      double precision  zeta1, zeta2, zeta3, det, w, weight
      double precision  c1x, c3x
      integer           i, ix, j, jx, k, ls(6)
      data              ls /1,2,3,4,5,6/
C     data              ls /1,3,5,7,9,11, 2,4,6,8,10,12/
C
C                   L O G I C
C
C      status =   ' '
      do 1200  j = 1,6
        do 1100  i = 1,6
          sm(i,j) = 0.0
 1100     continue
 1200   continue
C
      do 3000  k = 1,abs(p)
        call   TRIGGAUSSQ (p, k, zeta1, zeta2, zeta3, weight)
        call   TRIG6SHAPE (' ', zeta1,zeta2,zeta3, x,y, q,qx,qy, det)
        if (det .le. 0.0)        then
C          status = 'Negative Jacobian determinant'
          if (det .eq. 0.0)      then
C            status = 'Zero Jacobian determinant'
          end if
          return
        end if
        w =    weight *(0.5*det)* (q(1)+q(2)+q(3)
     $         +q(4)+q(5)+q(6))
C
        do 2000  j = 1,6
          jx =    ls(j)
          c1x =  qx(j) * w
          c3x =  qy(j) * w
          do 1500  i = j,6
            ix =     ls(i)
            sm(ix,jx) =  sm(ix,jx) + qx(i)*c1x + qy(i)*c3x
            sm(jx,ix) =  sm(ix,jx)
 1500       continue
 2000     continue
 3000   continue
C
      return
      end
C=END FORTRAN
