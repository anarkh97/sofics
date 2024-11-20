C
C     BARSLOSHFS forms the element mass matrix of a
C     two-node line for sloshing
C
C **********************************************************************
C
C     The calling sequence is
C
C       CALL   barsloshfs( XL, YL, H, SM, P, M)
C
C     where the input arguments are
C
C       XL,YL    (2 x 1) array of local coordinates of quadrilateral nodes
C       H        thicknesses at quadrilateral nodes (average) 
C       P        Gauss quadrature rule (no. of points)
C       M        Number of degrees of freedom of the element
C
C     The outputs are:
C
C       SM        (2 x 2) computed element stiffness matrix.
C
C=END USAGE
C=BLOCK FORTRAN

      subroutine    barsloshfs(xl, yl, h, sm, p, m)
C
C                   A R G U M E N T S
C
      integer   p, m
      real*8    xl(m), yl(m), h
      real*8    sm(m,m)
C
C                   L O C A L   V A R I A B L E S
C
      real*8  q(2)
      real*8  xi, det, w, weight
      integer           i, j, k
C
C                   L O G I C
C
      do 100  j = 1,2
        do 200  i = 1,2
          sm(i,j) = 0.0
200     continue
100   continue
C
      do 300  k = 1,p
          call LGAUSS (p, k, xi, weight)
C
C.... COMPUTE THE SHAPE FUNCTIONS
C
          q(1) =     0.5 * (1.0+xi)
          q(2) =     1.0 - q(1)
C
C.... COMPUTE THE DETERMINANT OF THE JACOBIAN
C
          det =     ((xl(2)-xl(1))**2 + (yl(2)-yl(1))**2) **0.5/2.0
          if (det .le. 0.0) then
            write(6,*)  'Negative Jacobian determinant'
            if (det .eq. 0.0)      then
              write(6,*) 'Zero Jacobian determinant'
            end if
            stop 
          end if
          w = weight * det * (q(1)+q(2)) * h

C... the sum of the shape functions is one, here we keep
C    the sum just in case we want to incorporate different 
C    thicknesses for each node
C
C.... COMPUTE THE THERMAL CONDUCTIVITY IN THE ELEMENT
C
          do 500  j = 1,2
	    do 600 i = j,2
              sm(i,j) = sm(i,j) + (q(i)*q(j))*w
              sm(j,i) = sm(i,j)
600         continue
500       continue
300   continue
C
      return
      end
