C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     This routine forms the element consistent mass matrix of a
C     four-node quadrilateral for 1 dof per node 
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       CALL   QUAD1DOFMASS ( X, Y, P, MM, M )
C
C     where the input arguments are
C
C       X         (4 x 1) array of x coordinates of quadrilateral nodes
C       Y         (4 x 1) array of y coordinates of quadrilateral nodes
C       P         Gauss quadrature rule (no. of points)
C       M         First dimension of MM in calling program.
C
C     The outputs are:
C
C       MM        (4 x 4) computed element mass matrix.
C                 As there is only one dof per node, we set  
C                 LS to 1,2,3,4  (see DATA statement below)
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine    q4d1dofmas(x, y, p, mm, m )
C
C                   A R G U M E N T S
C
      integer           p, m
      double precision  x(*), y(*)
      double precision  mm(m,*)
C
C                   L O C A L   V A R I A B L E S
C
      double precision  q(4), qx(4), qy(4)
      double precision  xi, eta, det, w, weight, c1x
      integer           i, ix, j, jx, k, l
      integer           ls(4)
C
C                   D A T A
C
      data              ls /1,2,3,4/
C
C                   L O G I C
C
      do 1200  j = 1,4
        do 1100  i = 1,4
          mm(i,j) = 0.0
 1100     continue
 1200   continue

C
      do 3000  k = 1,p
        do 2500  l = 1,p
          call     QGAUSS (p, k, p, l, xi, eta, weight)
          call     Q4SHPE (xi, eta, x, y, q, qx, qy, det)
          if (det .le. 0.0)        then
*           write(6,*)  'Negative Jacobian determinant'
            if (det .eq. 0.0)      then
*             write(6,*) 'Zero Jacobian determinant'
            end if
            stop 
          end if
          w =    weight * det *
     $          (q(1)+q(2)+q(3)+q(4))
C
          do 2000  j = 1,4
            jx =    ls(j)
            c1x =    q(j) * w
            do 1500  i = j,4
              ix =     ls(i)
              mm(ix,jx) =  mm(ix,jx) + q(i)*c1x 
              mm(jx,ix) =  mm(ix,jx)
 1500         continue
 2000       continue
 2500     continue
 3000   continue
 
C
      return
      end
