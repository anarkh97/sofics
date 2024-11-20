C
C     CONVECQUAD forms the element convective matrix of a
C     four-node quadrilateral for thermal bricks
C
C **********************************************************************
C
C     The calling sequence is
C
C       CALL   convecquad( XL, YL, C, SM, P, M)
C
C     where the input arguments are
C
C       XL,YL    (4 x 1) array of local coordinates of quadrilateral nodes
C       C        heat convection coefficients @ lateral faces
C       P        Gauss quadrature rule (no. of points)
C       M        Number of degrees of freedom of the element
C
C     The outputs are:
C
C       SM        (4 x 4) computed element stiffness matrix.
C
C=END USAGE
C=BLOCK FORTRAN

      subroutine    convecquad(xl, yl, c, sm, p, m)
C
C                   A R G U M E N T S
C
      integer   p, m
      real*8    xl(m), yl(m), c
      real*8    sm(m,m)
C
C                   L O C A L   V A R I A B L E S
C
      real*8  q(4), qx(4), qy(4)
      real*8  xi, eta, det, weight
      integer i, j, k, l
C
C                   L O G I C
C
      do 100  j = 1,4
        do 200  i = 1,4
          sm(i,j) = 0.0
200     continue
100   continue
C
      do 300  k = 1,p
        do 400  l = 1,p
          call QGAUSS (p, k, p, l, xi, eta, weight)
          call Q4SHPE (xi, eta, xl, yl, q, qx, qy, det)
          if (det .le. 0.0) then
            write(6,*)  'Negative Jacobian determinant'
            if (det .eq. 0.0)      then
              write(6,*) 'Zero Jacobian determinant'
            end if
            stop 
          end if

C... the sum of the shape functions is one, here we keep
C    the sum just in case we want to incorporate different 
C    thicknesses for each node
C
C.... COMPUTE THE CONVECTIVE MATRIX
C
          do 700 j = 1,4
            do 800 i = 1,4
              sm(i,j) = sm(i,j) + (q(i) * q(j)) *  c * weight * det
800         continue
700       continue
C
400     continue
300   continue
C
      return
      end
