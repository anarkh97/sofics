C=PURPOSE Form stiffness of 8-node hexahedron in 3-D  stress
C
C     HEXA8STIF forms the element stiffness matrix of a
C     eight-node hexahedron in 3-D stress.
C
C
C     where the input arguments are
C
C       OPT       Option letter argument, presently ignored.
C       X         (8 x 1) array of x coordinates of hexahedron nodes
C       Y         (8 x 1) array of y coordinates of hexahedron nodes
C       Z         (8 x 1) array of z coordinates of hexahedron nodes
C       C         (6 x 6) constitutive material matrix 
C       P         Gauss quadrature rule (no. of points)
C
C     The outputs are:
C
C       SM        (24 x 24) computed element stiffness matrix.  The
C                 arrangement of rows and columns pertains to node
C                 displacements arranged in the order
C                  (vx1, vy1, vz1, vx2, ... vz8)
C                 This particular ordering is obtained by setting array
C                 LS to 1,4,7,10,13,16,19,22,2,5,8,11,14,
C                 17,20,23,3,6,9,12,15,18,21,24
C                 ( see DATA statement below)
C
C       STATUS    Status integer variable.  Zero if no error
C                 detected.
C
      subroutine    brik8v(x, y, z, c, p, sm, status)
C
C                   A R G U M E N T S
C
      integer p,status
      real*8  x(8), y(8), z(8), c(6,6)
      real*8  sm(24,24)
C
C                   L O C A L   V A R I A B L E S
C
      real*8  q(8), qx(8), qy(8),qz(8)
      real*8  xi, eta, emu, det, w, weight
      real*8  c1x, c1y, c1z, c2x, c2y, c2z, c3x, c3y, c3z
      real*8  c4x, c4y, c4z, c5x, c5y, c5z, c6x, c6y, c6z
      integer           i, ix, iy, iz, j, jx, jy, jz, k, l
      integer           ls(24), jj
      double precision jSign
C
C                   D A T A
C
      data              ls /1,4,7,10,13,16,19,22,2,5,8,11,14,17,20,23,
     $                      3,6,9,12,15,18,21,24/
C
C                   L O G I C
C
C
C     Initialize no error
      status = 0
      jSign = 0
C
      do 1200  j = 1,24
        do 1100  i = 1,24
          sm(i,j) = 0.0
 1100     continue
 1200   continue
C
      do 3000  k = 1,p
        do 2500  l = 1,p
          do 2400  jj = 1,p
            call hxgaus(p,k,p,l,p,jj,xi,eta,emu,weight)
            call h8shpe(xi,eta,emu,x,y,z,q,qx,qy,qz,det)
C
            if(det .eq. 0.0) then
              write(6,*) 'Zero Jacobian determinant in brik8v.f'
              status = -1
              return
            endif

            if (det .lt. 0.0) then
              if(jSign .gt. 0.0) then
                write(6,*)
     &            'Sign changing Jacobian determinant in brik8v.f'
                status = -1
                return
              endif
C
              det = -det 
              jSign = -1.0
C
            else
              if (det .gt. 0.0) then
                if(jSign .lt. 0.0) then
                  write(6,*)
     &              'Sign changing Jacobian determinant in brik8v.f'
                  status = -1
                  return
                endif
C
                jSign = +1.0
C
              end if
	    end if
C
            w = weight * det 
C
            do 2000  j = 1,8
              jx =    ls(j)
              jy =    ls(j+8)
              jz =    ls(j+16)
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
              do 1500  i = j,8
                ix =     ls(i)
                iy =     ls(i+8)
                iz =     ls(i+16)
                sm(ix,jx) = sm(ix,jx)+qx(i)*c1x + qy(i)*c4x + qz(i)*c6x
                sm(jx,ix) = sm(ix,jx)
                sm(iy,jy) = sm(iy,jy)+qx(i)*c4y + qy(i)*c2y + qz(i)*c5y
                sm(jy,iy) = sm(iy,jy)
                sm(iz,jz) = sm(iz,jz)+qx(i)*c6z + qy(i)*c5z + qz(i)*c3z
                sm(jz,iz) = sm(iz,jz)
                sm(ix,jy) = sm(ix,jy)+qx(i)*c1y + qy(i)*c4y + qz(i)*c6y
                sm(iy,jx) = sm(iy,jx)+qx(i)*c4x + qy(i)*c2x + qz(i)*c5x
                sm(jy,ix) = sm(ix,jy)
                sm(jx,iy) = sm(iy,jx)
                sm(ix,jz) = sm(ix,jz)+qx(i)*c1z + qy(i)*c4z + qz(i)*c6z
                sm(iz,jx) = sm(iz,jx)+qx(i)*c6x + qy(i)*c5x + qz(i)*c3x
                sm(jz,ix) = sm(ix,jz)
                sm(jx,iz) = sm(iz,jx)
                sm(iy,jz) = sm(iy,jz)+qx(i)*c4z + qy(i)*c2z + qz(i)*c5z
                sm(iz,jy) = sm(iz,jy)+qx(i)*c6y + qy(i)*c5y + qz(i)*c3y
                sm(jz,iy) = sm(iy,jz)
                sm(jy,iz) = sm(iz,jy)
 1500           continue
 2000         continue
 2400       continue
 2500     continue
 3000   continue
 
      return
      end
