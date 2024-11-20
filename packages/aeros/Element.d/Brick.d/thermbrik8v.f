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
C       P         Gauss quadrature rule (no. of points)
C       M         First dimension of SM in calling program.
C
C     The outputs are:
C
C       SM        (24 x 24) computed element stiffness matrix.  The
C                 arrangement of rows and columns pertains to node
C                 displacements arranged in the order
C                  (vx1, vy1, vz1, vx2, ... vz8)
C                 This particular ordering is obtained by setting array
C                 LS to 1,2,3,4,5,6,7,8
C                 ( see DATA statement below)
C
C       STATUS    Status character variable.  Blank if no error
C                 detected.
C
      subroutine    thermbrik8v(x, y, z, p, sm, m, outerr)
C
C                   A R G U M E N T S
C
      integer   p, m,outerr
      real*8    x(8), y(8), z(8)
      real*8    sm(m,m)
C
C                   L O C A L   V A R I A B L E S
C
      real*8  q(8), qx(8), qy(8),qz(8)
      real*8  xi, eta, emu, det, w, weight
      real*8  cx, cy, cz
      integer i, j, k, l, jj
C
C
C                   L O G I C
C
C

      do 1200  j = 1,8
        do 1100  i = 1,8
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
            if (det .le. 0.0)        then
              write(outerr,*)'Negative Jacobian determinant'
              write(6,*)'Negative Jacobian determinant'
              if (det .eq. 0.0)      then
                write(outerr,*)'Negative Jacobian determinant'
                write(6,*)'Negative Jacobian determinant'
              end if
              return
            end if
C
            w =    weight * det 
C
            do 2000  j = 1,8
             cx = qx(j) * w
             cy = qy(j) * w                
             cz = qz(j) * w
              do 1500  i = j,8
               sm(i,j) = sm(i,j) + qx(i)*cx + qy(i)*cy + qz(i)*cz
               sm(j,i) = sm(i,j)
 1500           continue
 2000         continue
 2400       continue
 2500     continue
 3000   continue
 
C
      return
      end
