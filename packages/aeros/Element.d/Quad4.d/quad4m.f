      subroutine    quad4m(x, y, h, c, p, sm, m, rip)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     QUAD4MSTIF forms the element stiffness matrix of a
C     four-node quadrilateral in plane stress or plane strain.
C
C       X         (4 x 1) array of x coordinates of quadrilateral nodes
C       Y         (4 x 1) array of y coordinates of quadrilateral nodes
C       H         (4 x 1) array of thicknesses at quadrilateral nodes
C       CS        (4 x *) constitutive material matrix (not
C                         integrated thorugh the thickness)
C       P         Gauss quadrature rule (no. of points)
C       M         First dimension of SM in calling program
C       RIP       Plane stress - strain:  flag = 0 - plane stress
C                                         flag = 1 - plane strain
C
C     The outputs are:
C
C       SM        (8 x 8) computed element stiffness matrix.  The
C                 arrangement of rows and columns pertains to node
C                 displacements arranged in the order
C                  (vx1, vy1, vx2, ... vy4)
C                 This particular ordering is obtained by setting array
C                 LS to 1,3,5,7,2,4,6,8  (see DATA statement below)
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C                   A R G U M E N T S
C
      integer           p, m
      double precision  rip
      double precision  x(4), y(4), h(4), c(3,3)
      double precision  sm(8,8)
C
C                   L O C A L   V A R I A B L E S
C
C     double precision  efac
      double precision  q(4), qx(4), qy(4)
      double precision  xi, eta, det, w, weight
      double precision  c1x, c1y, c2x, c2y, c3x, c3y
      integer           i, ix, iy, j, jx, jy, k, l
      integer           ls(8)
C
C                   D A T A
C
      data              ls /1,3,5,7,2,4,6,8/
C     data              ids /1,2,4/
c
c     initialize stiffness matrix
C
      do 1200  j = 1,8
        do 1100  i = 1,8
          sm(i,j) = 0.0
 1100     continue
 1200   continue
c
c     Gauss quadrature loop
C
      do 3000  k = 1,p
        do 2500  l = 1,p
c	
          call QGAUSS (p, k, p, l, xi, eta, weight)
          call Q4SHPE (xi, eta, x, y, q, qx, qy, det)
C
          if (det .le. 0.0)        then
            write(6,*)  'Negative Jacobian determinant in 4 node quad'
            if (det .eq. 0.0)      then
              write(6,*) 'Zero Jacobian determinant in 4 node quad'
            end if
            stop 
          end if
C
          w =    weight * det *
     $          (h(1)*q(1)+h(2)*q(2)+h(3)*q(3)+h(4)*q(4))
C
          do 2000  j = 1,4
            jx =    ls(j)
            jy =    ls(j+4)
            c1x =  (c(1,1)*qx(j)+c(1,3)*qy(j)) * w
            c1y =  (c(1,3)*qx(j)+c(1,2)*qy(j)) * w
            c2x =  (c(1,2)*qx(j)+c(2,3)*qy(j)) * w
            c2y =  (c(2,3)*qx(j)+c(2,2)*qy(j)) * w
            c3x =  (c(1,3)*qx(j)+c(3,3)*qy(j)) * w
            c3y =  (c(3,3)*qx(j)+c(2,3)*qy(j)) * w
            do 1500  i = j,4
              ix =     ls(i)
              iy =     ls(i+4)
              sm(ix,jx) =  sm(ix,jx) + qx(i)*c1x + qy(i)*c3x
              sm(jx,ix) =  sm(ix,jx)
              sm(iy,jy) =  sm(iy,jy) + qx(i)*c3y + qy(i)*c2y
              sm(jy,iy) =  sm(iy,jy)
              sm(ix,jy) =  sm(ix,jy) + qx(i)*c1y + qy(i)*c3y
              sm(iy,jx) =  sm(iy,jx) + qx(i)*c3x + qy(i)*c2x
              sm(jy,ix) =  sm(ix,jy)
              sm(jx,iy) =  sm(iy,jx)
 1500         continue
 2000       continue
 2500     continue
 3000   continue
 
      return
      end
