C=PURPOSE Form stiffness of 20-node iso-P hexahedron
C
C     The input arguments are
C
C       X         (20 x 1) array of x coordinates of hexahedron nodes
C       Y         (20 x 1) array of y coordinates of hexahedron nodes
C       Z         (20 x 1) array of z coordinates of hexahedron nodes
C       C         (6 x 6) constitutive material matrix 
C       P         Gauss quadrature rule (no. of points in each dir)  Usually p=3.
C
C     The outputs are:
C
C       SM        (60 x 60) computed element stiffness matrix.  The
C                 arrangement of rows and columns pertains to node
C                 displacements arranged in the order
C                  (ux1, uy1, uz1, ux2, ... uz20)
C                 This particular ordering is obtained by setting array
C                 LS as shown by the DATA statement below
C
C       STATUS    if error was detected.
C-----------------------------------------------------------------------------
      subroutine  brik20v(x, y, z, c, p, sm, status)
C
C                   A R G U M E N T S
C
      integer           p,status
      double precision  x(20), y(20), z(20), c(6,6)
      double precision  sm(60,60)
C
C                   L O C A L   V A R I A B L E S
C
      double precision  q(20), qx(20), qy(20),qz(20)
      double precision  xi, eta, mu, det, w, weight
      double precision  c1x, c1y, c1z, c2x, c2y, c2z, c3x, c3y, c3z
      double precision  c4x, c4y, c4z, c5x, c5y, c5z, c6x, c6y, c6z
      integer           i, ix, iy, iz, j, jx, jy, jz, k, l, m
      integer           ls(60)
C
C                   D A T A
C
      data              ls 
     $    /1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,
     $     2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,
     $     3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60 /
C
C                   L O G I C
C
      status =  0
      do 1200  j = 1,60
        do 1100  i = 1,60
          sm(i,j) = 0.0d0
 1100     continue
 1200   continue
C
      do 3000  k = 1,p
        do 2500  l = 1,p
          do 2400  m = 1,p
            call hxgaus20 (p, k, p, l, p, m, xi, eta, mu, weight)
            call h20shpe (xi, eta, mu, x, y, z, q, qx, qy, qz, det)
            if (det .le. 0.0d0) then
               write(6,*) 'Negative Jacobian determinant in brik20v.f'
               status = -1
               return
            end if
C
            w =    weight * det 
C
            do 2000  j = 1,20
              jx =    ls(j)
              jy =    ls(j+20)
              jz =    ls(j+40)
              c1x = (c(1,1)*qx(j) + c(1,4)*qy(j) + c(1,6)*qz(j))*w
              c1y = (c(1,4)*qx(j) + c(1,2)*qy(j) + c(1,5)*qz(j))*w
              c1z = (c(1,6)*qx(j) + c(1,5)*qy(j) + c(1,3)*qz(j))*w
              c2x = (c(2,1)*qx(j) + c(2,4)*qy(j) + c(2,6)*qz(j))*w
              c2y = (c(2,4)*qx(j) + c(2,2)*qy(j) + c(2,5)*qz(j))*w
              c2z = (c(2,6)*qx(j) + c(2,5)*qy(j) + c(2,3)*qz(j))*w
              c3x = (c(3,1)*qx(j) + c(3,4)*qy(j) + c(3,6)*qz(j))*w
              c3y = (c(3,4)*qx(j) + c(3,2)*qy(j) + c(3,5)*qz(j))*w
              c3z = (c(3,6)*qx(j) + c(3,5)*qy(j) + c(3,3)*qz(j))*w
              c4x = (c(4,1)*qx(j) + c(4,4)*qy(j) + c(4,6)*qz(j))*w
              c4y = (c(4,4)*qx(j) + c(4,2)*qy(j) + c(4,5)*qz(j))*w
              c4z = (c(4,6)*qx(j) + c(4,5)*qy(j) + c(4,3)*qz(j))*w
              c5x = (c(5,1)*qx(j) + c(5,4)*qy(j) + c(5,6)*qz(j))*w
              c5y = (c(5,4)*qx(j) + c(5,2)*qy(j) + c(5,5)*qz(j))*w
              c5z = (c(5,6)*qx(j) + c(5,5)*qy(j) + c(5,3)*qz(j))*w
              c6x = (c(6,1)*qx(j) + c(6,4)*qy(j) + c(6,6)*qz(j))*w
              c6y = (c(6,4)*qx(j) + c(6,2)*qy(j) + c(6,5)*qz(j))*w
              c6z = (c(6,6)*qx(j) + c(6,5)*qy(j) + c(6,3)*qz(j))*w
              do 1500  i = j,20
                ix =     ls(i)
                iy =     ls(i+20)
                iz =     ls(i+40)
                sm(ix,jx) =sm(ix,jx)+qx(i)*c1x + qy(i)*c4x + qz(i)*c6x
                sm(jx,ix) =sm(ix,jx)
                sm(iy,jy) =sm(iy,jy)+qx(i)*c4y + qy(i)*c2y + qz(i)*c5y
                sm(jy,iy) =sm(iy,jy)
                sm(iz,jz) =sm(iz,jz)+qx(i)*c6z + qy(i)*c5z + qz(i)*c3z
                sm(jz,iz) =sm(iz,jz)
                sm(ix,jy) =sm(ix,jy)+qx(i)*c1y + qy(i)*c4y + qz(i)*c6y
                sm(iy,jx) =sm(iy,jx)+qx(i)*c4x + qy(i)*c2x + qz(i)*c5x
                sm(jy,ix) =sm(ix,jy)
                sm(jx,iy) =sm(iy,jx)
                sm(ix,jz) =sm(ix,jz)+qx(i)*c1z + qy(i)*c4z + qz(i)*c6z
                sm(iz,jx) =sm(iz,jx)+qx(i)*c6x + qy(i)*c5x + qz(i)*c3x
                sm(jz,ix) =sm(ix,jz)
                sm(jx,iz) =sm(iz,jx)
                sm(iy,jz) =sm(iy,jz)+qx(i)*c4z + qy(i)*c2z + qz(i)*c5z
                sm(iz,jy) =sm(iz,jy)+qx(i)*c6y + qy(i)*c5y + qz(i)*c3y
                sm(jz,iy) =sm(iy,jz)
                sm(jy,iz) =sm(iz,jy)
 1500           continue
 2000         continue
 2400       continue
 2500     continue
 3000   continue
C
      return
      end
