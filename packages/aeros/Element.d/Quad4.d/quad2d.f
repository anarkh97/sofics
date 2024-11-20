C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C MODIFIED BY PAUL STERN MARCH 7 1990
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     QUAD2D forms the coupling  matrix of a
C     four-node quadrilateral. 
C
C
C     The calling sequence is
C
C       CALL   QUAD2D ( X, Y, C, P, SM, M )
C
C     where the input arguments are
C
C       X         (4 x 1) array of x coordinates of quadrilateral nodes
C       Y         (4 x 1) array of y coordinates of quadrilateral nodes
C       C         element coupling coefficient
C       P         Gauss quadrature rule (no. of points)
C       M         First dimension of SM in calling program.
C
C     The outputs are:
C
C       SM        (8 x 4) computed element coupling  matrix.
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine    quad2d(x, y, c, p, cm, m )
C
C                   A R G U M E N T S
C
      integer           p, m
      real*8  x(*), y(*)
      real*8  c
      real*8  cm(m,*)
C
C                   L O C A L   V A R I A B L E S
C
      real*8  q(4), qx(4), qy(4)
      real*8  xi, eta, det, weight, w
      integer           i, j, k, l
C
C                   L O G I C
C
      do 100  j = 1,4
        do 200  i = 1,8
          cm(i,j) = 0.0
200     continue
100   continue

C
C.... COMPUTE THE COUPLING ELEMENT MATRIX            
C
      do 210 k = 1, p
        do 110 l = 1, p
C
C.... COMPUTE THE SHAPE FUNCTIONS & DERIVATIVES
C.... NOTE THE THERMAL AND MECHANICS USE THE SAME SFS
C
          call qgauss(p,k,p,l,xi,eta,weight)
          call q4shpe(xi,eta,x,y,q,qx,qy,det)
C
C.... CHECK THE DETERMINANT TO SEE IF THE QUAD IS TO DISTORTED
C
	  if (det .le. 0.0)        then
            write(6,*)  'Negative Jacobian determinant'
            if (det .eq. 0.0)      then
              write(6,*) 'Zero Jacobian determinant'
            end if
            stop
          end if
          w =    weight * det
C
C.... ASSEMBLE THE ELEMENT COUPLING MATRIX
C
          cm(1,1) = cm(1,1) + c*qx(1)*q(1)*w
          cm(1,2) = cm(1,2) + c*qx(1)*q(2)*w
          cm(1,3) = cm(1,3) + c*qx(1)*q(3)*w
          cm(1,4) = cm(1,4) + c*qx(1)*q(4)*w
C
          cm(2,1) = cm(2,1) + c*qy(1)*q(1)*w
          cm(2,2) = cm(2,2) + c*qy(1)*q(2)*w
          cm(2,3) = cm(2,3) + c*qy(1)*q(3)*w
          cm(2,4) = cm(2,4) + c*qy(1)*q(4)*w
C
          cm(3,1) = cm(3,1) + c*qx(2)*q(1)*w
          cm(3,2) = cm(3,2) + c*qx(2)*q(2)*w
          cm(3,3) = cm(3,3) + c*qx(2)*q(3)*w
          cm(3,4) = cm(3,4) + c*qx(2)*q(4)*w
C
          cm(4,1) = cm(4,1) + c*qy(2)*q(1)*w
          cm(4,2) = cm(4,2) + c*qy(2)*q(2)*w
          cm(4,3) = cm(4,3) + c*qy(2)*q(3)*w
          cm(4,4) = cm(4,4) + c*qy(2)*q(4)*w
C
          cm(5,1) = cm(5,1) + c*qx(3)*q(1)*w
          cm(5,2) = cm(5,2) + c*qx(3)*q(2)*w
          cm(5,3) = cm(5,3) + c*qx(3)*q(3)*w
          cm(5,4) = cm(5,4) + c*qx(3)*q(4)*w
C
          cm(6,1) = cm(6,1) + c*qy(3)*q(1)*w
          cm(6,2) = cm(6,2) + c*qy(3)*q(2)*w
          cm(6,3) = cm(6,3) + c*qy(3)*q(3)*w
          cm(6,4) = cm(6,4) + c*qy(3)*q(4)*w
C
          cm(7,1) = cm(7,1) + c*qx(4)*q(1)*w
          cm(7,2) = cm(7,2) + c*qx(4)*q(2)*w
          cm(7,3) = cm(7,3) + c*qx(4)*q(3)*w
          cm(7,4) = cm(7,4) + c*qx(4)*q(4)*w
C
          cm(8,1) = cm(8,1) + c*qy(4)*q(1)*w
          cm(8,2) = cm(8,2) + c*qy(4)*q(2)*w
          cm(8,3) = cm(8,3) + c*qy(4)*q(3)*w
          cm(8,4) = cm(8,4) + c*qy(4)*q(4)*w
C
110     continue
210   continue
C
      return
      end
