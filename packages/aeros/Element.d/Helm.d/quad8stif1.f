C=DECK QUAD8STIF1
C=PURPOSE Form stiffness of 8-node quad in plane stress
C=AUTHOR C. A. Felippa,  February 1967
C=AUTHOR Antonini Puppin-Macedo (adapted for 1 dof per node), June 1999
C=VERSION May 1982 (Fortran 77)
C=EQUIPMENT Machine independent
C=KEYWORDS plane stress eight node quadrilateral
C=KEYWORDS finite element stiffness matrix
C=BLOCK ABSTRACT
C
C     QUAD8MSTIF forms the element stiffness matrix of a
C     eight-node quadrilateral in plane stress.
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       CALL   QUAD8STIF1 (X, Y, P, SM, M)
C
C     where the input arguments are
C
C       X         (8 x 1) array of x coordinates of quadrilateral nodes
C       Y         (8 x 1) array of y coordinates of quadrilateral nodes
C       P         Gauss quadrature rule (no. of points)
C       M         First dimension of SM in calling program.
C
C     The outputs are:
C
C       SM        (16 x 16) computed element stiffness matrix.  The
C                 arrangement of rows and columns pertains to node
C                 displacements arranged in the order
C                  (vx1, vy1, vx2, ... vy8)
C                 This particular ordering is obtained by setting array
C                 LS to 1,3,5, ... 16  (see DATA statement below)
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine    QUAD8STIF1
     $          (x, y, p, sm, m)
C
C                   A R G U M E N T S
C
      integer           p, m
      double precision  x(8), y(8)
      double precision  sm(m,m)
C
C                   L O C A L   V A R I A B L E S
C
      double precision  q(8), qx(8), qy(8)
      double precision  xi, eta, det, w, weight
      double precision  c1x, c3x
      integer           i, ix, j, jx, k, l
      integer           ls(8)
C
C                   D A T A
C
      data    ls /1,2,3,4,5,6,7,8/
C      data    ls /1,3,5,7,9,11,13,15, 2,4,6,8,10,12,14,16/
C
C                   L O G I C
C
C      status =   ' '
      do 1200  j = 1,8
        do 1100  i = 1,8
          sm(i,j) = 0.0
 1100     continue
 1200   continue
C
      do 3000  k = 1,p
        do 2500  l = 1,p
          call     QUADGAUSSQ (p, k, p, l, xi, eta, weight)
          call     QUAD8SHAPE (' ', xi, eta, x, y, q, qx, qy, det)
          if (det .le. 0.0)        then
C            status = 'Negative Jacobian determinant'
            if (det .eq. 0.0)      then
C              status = 'Zero Jacobian determinant'
            end if
            return
          end if
          w =    weight * det *
     $          (q(1)+q(2)+q(3)+q(4)
     $          +q(5)+q(6)+q(7)+q(8))
C
          do 2000  j = 1,8
            jx =    ls(j)
            c1x =  qx(j) * w
            c3x =  qy(j) * w
            do 1500  i = j,8
              ix =     ls(i)
              sm(ix,jx) =  sm(ix,jx) + qx(i)*c1x + qy(i)*c3x
              sm(jx,ix) =  sm(ix,jx)
 1500         continue
 2000       continue
 2500     continue
 3000   continue
C
      return
      end
C=END FORTRAN
