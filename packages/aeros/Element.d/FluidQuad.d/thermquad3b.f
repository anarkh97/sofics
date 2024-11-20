C
C     THERMQUAD3B forms the element stiffness matrix of a
C     four-node quadrilateral for heat conduction
C
C **********************************************************************
C
C     The calling sequence is
C
C       CALL   thermquad3b( XL, YL, KO, C, H, SM, P, M)
C
C     where the input arguments are
C
C       XL,YL    (4 x 1) array of local coordinates of quadrilateral nodes
C       H        thicknesses at quadrilateral nodes (average) 
C       KO        heat conduction coefficient
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

      subroutine    thermquad3b(xl, yl, ko, c, h, sm, p, m)
C
C                   A R G U M E N T S
C
      integer   p, m
      real*8    xl(m), yl(m), h, ko, c
      real*8    sm(m,m)
C
C                   L O C A L   V A R I A B L E S
C
      real*8  q(4), qx(4), qy(4)
      real*8  xi, eta, det, w, weight
      integer           i, j, k, l
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
          w = weight * det * 
     &          (q(1)+q(2)+q(3)+q(4)) * h

C... the sum of the shape functions is one, here we keep
C    the sum just in case we want to incorporate different 
C    thicknesses for each node
C
C.... COMPUTE THE THERMAL CONDUCTIVITY IN THE ELEMENT
C
          do 500  j = 1,4
	    do 600 i = j,4
              sm(i,j) = sm(i,j) + ko*(qx(i)*qx(j) + qy(i)*qy(j))*w
              sm(j,i) = sm(i,j)
600         continue
500       continue
C
C.... COMPUTE THE CONVECTIVE HEAT TRANFER AT THE LATERAL FACES
C
C          do 700 j = 1,4
C            do 800 i = 1,4
C              sm(i,j) = sm(i,j) + (q(i) * q(j)) *  c * weight * det
C 800         continue
C 700       continue
C
400     continue
300   continue
C
      return
      end
