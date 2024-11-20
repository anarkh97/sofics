C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Gregrory W. Brown
C January, 2000
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     SHEARPANEL forms the element stiffness matrix of a
C     four-node shear panel element.
C
C     Element is built as a plane element, then rotated to 3D
C
C     WARNING: THIS COULD RESULT IN A MECHANISM IF ELEMENT
C             IS NOT ATTACHED TO NODES CONTAINING FULL STIFFNESS
C     The calling sequence is
C
C       call   shearpanel (xg,yg,zg,t,G,E,F1,F2,p,Kshear,m,vol)
C
C     where the input arguments are
C
C       xg        (4 x 1) array of global x coordinates
C       yg        (4 x 1) array of global y coordinates
C       zg        (4 x 1) array of global y coordinates
C       t         thickness of quadrilateral nodes
C       G         Shear Modulus of element
C       E         Youngs Modulus of element
C       F1        Extensional factor 1
C       F2        Extensional factor 2
C       p         Gauss quadrature rule (no. of points)
C       m         First dimension of SM in calling program.
C
C     The outputs are:
C
C       Kshear    (12x12) computed element stiffness matrix.  The
C                 arrangement of rows and columns pertains to node
C                 displacements arranged in the order
C                  (vx1, vy1, yz1, vx2, ... vz4)
C                 This particular ordering is obtained by setting array
C                 LSG to  1,4,7,10,2,5,8,11,3,6,9,12 
C                 (see DATA statement below)
C
C       vol       Volume of element
C
      subroutine    shearpanel(xg,yg,zg,t,G,E,F1,F2,p,Kshear,m,vol)
C
C                   A R G U M E N T S
C
      integer p, m
      real*8  xg(*), yg(*), zg(*), t, G, E, F1, F2
      real*8  Kshear(m,*), vol
C
C                   L O C A L   V A R I A B L E S
C
C Local Coordinates
      real*8  x(4), y(4), z(4)
C Working Vectors
      real*8  node(3,4), local(3,4) ,v1(3),v2(3),v3(3),v4(3), len
C Consitutive Matrix and Thickness
      real*8  c(3,3), h(4)
C In-plane Stiffness Matrix, Extensional Stiffness Matrix
      real*8  sm(8,8), ex(8,8)
C Effective Extentional Areas
      real*8  area(4)
C Rod Coordinates and Stiffness
      real*8  xr(2), yr(2), zr(2)
      real*8  rodstif(6,6)
      integer rnod(2,4)
C Output Control
      logical output
C Other Variables
      real*8  q(4), qx(4), qy(4)
      real*8  xi, eta, det, w, weight
      real*8  c1x, c1y, c2x, c2y, c3x, c3y
      integer i, ix, iy, j, jx, jy, k, l
      integer ls(8), irx, jrx, iry, jry
      integer lsg(12), igx, jgx, igy, jgy, igz, jgz
C
C                   D A T A
C
      data    ls  /1,3,5,7,2,4,6,8/
      data    lsg /1,4,7,10,2,5,8,11,3,6,9,12/
C
C Set Output Status
      output = .false.
C
C                   L O G I C
C
      do 1200  j = 1,12
        do 1100  i = 1,12
          Kshear(i,j) = 0.0
 1100   continue
 1200 continue
      do 1201  j = 1,8
        do 1101  i = 1,8
          sm(i,j) = 0.0
          ex(i,j) = 0.0
 1101   continue
 1201 continue
C
C... Output Properties of Element
      if (output) then
        write(6,*)
        write(6,*) 'SHEARPANEL.F'
        write(6,*)
        write(6,*) 'Node Coordinates (global)'
        write(6,'(1x,A4,1x,A15,1x,A15,1x,A15)') 'Node',
     $    ' X Coordinate  ',' Y Coordinate  ',' Z Coordinate  '
        write(6,'(1x,A4,1x,A15,1x,A15,1x,A15)') '----',
     $    '---------------','---------------','--------------'
        do 100 i=1,4
          write(6,'(3x,I1,2x,E15.8,1x,E15.8,1x,E15.8)')
     $      i,xg(i),yg(i),zg(i)
 100    continue
        write(6,*)
        write(6,'(1x,A17,E15.8)') 'Youngs Modulus = ',E
        write(6,'(1x,A17,E15.8)') 'Shear Modulus  = ',G
        write(6,'(1x,A17,E15.8)') 'Thickness      = ',t
      endif
C
C     COMPUTE LOCAL COORDINATE SYSTEM 
C  (works best if element is plane)
C
C      Y                             
C     /|\                            
C      |  4            3             
C      |  /------------|             
C      | /             |             
C      |/              |             
C      /---------------|-----> X     
C      1               2             
C
C... store nodal coordinates in "node"
C... shift origin to node 1
C
      do 1300  i = 1,4
        node(1,i) = xg(i)-xg(1)
        node(2,i) = yg(i)-yg(1)
        node(3,i) = zg(i)-zg(1)
 1300 continue
C
C... Vector from 1->2 (Local X axis = v1)
      call unitv(node(1,1),node(1,2),v1)
C
C... Vector from 4->2
      call unitv(node(1,4),node(1,2),v2)
C
C... Vector from 1->3
      call unitv(node(1,1),node(1,3),v3)
C
C... Define perpindicular as cross between v2 and v3
C
      call cross(v2,v3,v4)
      call normve(v4)
      call length(v4,len)
      if (len .eq. 0.0)        then
        write(6,*)  'ERROR: subroutine shearpanel.f'
        write(6,*)  'Element has zero area.'
        stop
      endif
C
C... Cross v4 and v1 (Local Y axis = v2)
C
      call cross(v4,v1,v2)
      call normve(v2)
C
C... Cross v1 and v2 (Local Z axis = v3)
C
      call cross(v1,v2,v3)
      call normve(v3)
C
C... Compute Local Coordinates: store in "local"
C
      do 1310 i = 1,4
        x(i) = 0.0
        y(i) = 0.0
        z(i) = 0.0
        do 1315 j = 1,3
          local(j,i) = 0.0
 1315   continue
 1310 continue
C
C... Output Local Frame Vectors
      if (output) then
        write(6,*)
        write(6,*) 'Local Frame Vectors'
        write(6,'(1x,A4,1x,A15,1x,A15,1x,A15)') 'Dir ',
     $    ' X Coordinate  ',' Y Coordinate  ',' Z Coordinate  '
        write(6,'(1x,A4,1x,A15,1x,A15,1x,A15)') '----',
     $    '---------------','---------------','--------------'
        write(6,'(3x,I1,2x,E15.8,1x,E15.8,1x,E15.8)')
     $    1,v1(1),v1(2),v1(3)
        write(6,'(3x,I1,2x,E15.8,1x,E15.8,1x,E15.8)')
     $    2,v2(1),v2(2),v2(3)
        write(6,'(3x,I1,2x,E15.8,1x,E15.8,1x,E15.8)')
     $    3,v3(1),v3(2),v3(3)
      endif
C
C... v1,v2,v3 form transformation matrix
C... local = T*node
C
      do 1320 i = 1,4
        do 1325 j = 1,3
          local(1,i) = local(1,i) + v1(j)*node(j,i)
          local(2,i) = local(2,i) + v2(j)*node(j,i)
          local(3,i) = local(3,i) + v3(j)*node(j,i)
 1325   continue
 1320   continue
C
C... Store Local Coordinates in x,y,z
C
      do 1330 i = 1,4
        x(i) = local(1,i)
        y(i) = local(2,i)
        z(i) = local(3,i)
 1330 continue
C
C... Output Local Coordinates
      if (output) then
        write(6,*)
        write(6,*) 'Node Coordinates (local)'
        write(6,'(1x,A4,1x,A15,1x,A15,1x,A15)') 'Node',
     $  ' X Coordinate  ',' Y Coordinate  ',' Z Coordinate  '
        write(6,'(1x,A4,1x,A15,1x,A15,1x,A15)') '----',
     $  '---------------','---------------','--------------'
        do 110 i=1,4
          write(6,'(3x,I1,2x,E15.8,1x,E15.8,1x,E15.8)')
     $    i,x(i),y(i),z(i)
 110    continue
      endif
C
C... Check amount of out-of-plane behaviore
C    Use x coordinate of node 2 as representative length
C
      do 1340 i = 1,4
        len = abs(local(3,i)/local(1,2))
        if (len.gt.0.1) then
          write(6,*)  'WARNING: subroutine shearpanel.f'
          write(6,*)  'Element has warp > 10%'
          write(6,*)  'Warp = ',len
        endif
 1340 continue
C
C
C... Create consitutive matrix relevant to Shear element
C    The only non-zero entry is for the shear c(3,3) = G
C
      do 1400 i = 1,3
        do 1410 j = 1,3
          c(i,j) = 0.0
 1410   continue
 1400 continue
      c(3,3) = G
C
C... Element uses uniform thickness
C
      do 1420 i = 1,4
        h(i) = t
 1420 continue
C
C... Compute stiffness matrix for shear effects
C    This is identical to the in-plane quad
C
      vol = 0.0
      do 3000  k = 1,p
        do 2500  l = 1,p
          call     QGAUSS (p, k, p, l, xi, eta, weight)
          call     Q4SHPE (xi, eta, x, y, q, qx, qy, det)
          if (det .le. 0.0)        then
            write(6,*)  'ERROR: subroutine shearpanel.f'
            write(6,*)  'Negative Jacobian determinant in 4 node quad'
            if (det .eq. 0.0)      then
              write(6,*)  'ERROR: subroutine shearpanel.f'
              write(6,*) 'Zero Jacobian determinant in 4 node quad'
            end if
            stop 
          end if
C
          vol = vol + det
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
 1500       continue
 2000     continue
 2500   continue
 3000 continue
C
C... Determine the volume of the element
      vol = vol*t
      if (output) then
        write(6,*)
        write(6,'(1x,A17,E15.8)') 'Element Volume = ',vol
      endif
C
C... Output Shear Stiffness
      if (output) then
        write(6,*)
        write(6,*) 'Shear Stiffness (local)'
        write(6,'(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11)') 'Row',
     $  ' Column 1  ',' Column 2  ',' Column 3  ',' Column 4  '
        write(6,'(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11)') '---',
     $  '-----------','-----------','-----------','-----------'
        do 201 i=1,4
          write(6,'(2x,I1,2x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    i,(sm(i,j),j=1,4)
 201    continue
        write(6,'(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11)') 'Row',
     $  ' Column 5  ',' Column 6  ',' Column 7  ',' Column 8  '
        write(6,'(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11)') '---',
     $  '-----------','-----------','-----------','-----------'
        do 202 i=1,4
          write(6,'(2x,I1,2x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    i,(sm(i,j),j=5,8)
 202    continue
        write(6,'(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11)') 'Row',
     $  ' Column 1  ',' Column 2  ',' Column 3  ',' Column 4  '
        write(6,'(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11)') '---',
     $  '-----------','-----------','-----------','-----------'
        do 203 i=5,8
          write(6,'(2x,I1,2x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    i,(sm(i,j),j=1,4)
 203    continue
        write(6,'(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11)') 'Row',
     $  ' Column 5  ',' Column 6  ',' Column 7  ',' Column 8  '
        write(6,'(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11)') '---',
     $  '-----------','-----------','-----------','-----------'
        do 204 i=5,8
          write(6,'(2x,I1,2x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    i,(sm(i,j),j=5,8)
 204    continue
      endif
C
C... Compute Extensional Stiffness
C    This is defined as rods along the edges
C
C    Areas are based on info for NASTRAN CSHEAR element
C
C... Sides 1-2 and 3-4, based on F1
C    Compute "height" of panel
      len = (local(2,4)-local(2,1))+(local(2,3)-local(2,2))
      len = abs(len/2.0)
C
      if (F1.ne.0) then
C         
C    Compute rod area
        if (F1.le.(1.01)) then
          area(1) = 0.5*F1*t*len
        else
          area(1) = 0.5*F1*t*t
        endif
C
      else
        area(1) = 0.0
      endif
      area(2) = area(1)
C
C... Sides 1-4 and 2-3, based on F2
C    Compute "width" of panel
      len = (local(1,3)-local(1,4))+(local(1,2)-local(1,1))
      len = abs(len/2.0)
C
      if (F2.ne.0) then
C         
C    Compute rod area
        if (F2.le.(1.01)) then
          area(3) = 0.5*F2*t*len
        else
          area(3) = 0.5*F2*t*t
        endif
C
      else
        area(3) = 0.0
      endif
      area(4) = area(3)
C
C Local Z contribution to rods is ignored 
      zr(1) = 0.0
      zr(2) = 0.0
C
C Node definition of each rod
C Side 1-2
      rnod(1,1) = 1
      rnod(2,1) = 2
C Side 3-4
      rnod(1,2) = 3
      rnod(2,2) = 4
C Side 1-4
      rnod(1,3) = 1
      rnod(2,3) = 4
C Side 2-3
      rnod(1,4) = 2
      rnod(2,4) = 3
C
C Begin Loop Over 4 Rods
      do 4000 k = 1,4
C
C Get nodes of rod
        do 4010 i = 1,2
          xr(i) = x(rnod(i,k))
          yr(i) = y(rnod(i,k))
 4010   continue
C
C Zero rod stiffness matrix
        do 4020 i = 1,6   
         do 4030 j = 1,6   
           rodstif(i,j)=0.0
 4030    continue
 4020   continue
C
C Get rod stiffness
        call b3dstf(E,rodstif,6,area(k),xr,yr,zr,w)
C
C Add rod stiffness to element stiffness
        do 4040 j = 1,2
          jx =    ls(rnod(j,k))
          jy =    ls(rnod(j,k)+4)
          jrx = (3*j)-2
          jry = (3*j)-1
          do 4050 i = 1,2
            ix =    ls(rnod(i,k))
            iy =    ls(rnod(i,k)+4)
            irx = (3*i)-2
            iry = (3*i)-1
            ex(ix,jx) = ex(ix,jx) + rodstif(irx,jrx)
            ex(ix,jy) = ex(ix,jy) + rodstif(irx,jry)
            ex(iy,jx) = ex(iy,jx) + rodstif(iry,jrx)
            ex(iy,jy) = ex(iy,jy) + rodstif(iry,jry)
 4050     continue
 4040   continue
C
C End Loop on Rods
 4000 continue
C
C Add to SM
      do 4200  j = 1,8
        do 4100  i = 1,8
          sm(i,j) = sm(i,j) + ex(i,j)
 4100   continue
 4200 continue
C
C... Output Rod Stiffness
      if (output.and.((F1.ne.0).or.(F2.ne.0))) then
        write(6,*)
        write(6,*) 'Extensional Stiffness (local)'
        write(6,'(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11)') 'Row',
     $  ' Column 1  ',' Column 2  ',' Column 3  ',' Column 4  '
        write(6,'(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11)') '---',
     $  '-----------','-----------','-----------','-----------'
        do 301 i=1,4
          write(6,'(2x,I1,2x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    i,(ex(i,j),j=1,4)
 301    continue
        write(6,'(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11)') 'Row',
     $  ' Column 5  ',' Column 6  ',' Column 7  ',' Column 8  '
        write(6,'(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11)') '---',
     $  '-----------','-----------','-----------','-----------'
        do 302 i=1,4
          write(6,'(2x,I1,2x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    i,(ex(i,j),j=5,8)
 302    continue
        write(6,'(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11)') 'Row',
     $  ' Column 1  ',' Column 2  ',' Column 3  ',' Column 4  '
        write(6,'(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11)') '---',
     $  '-----------','-----------','-----------','-----------'
        do 303 i=5,8
          write(6,'(2x,I1,2x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    i,(ex(i,j),j=1,4)
 303    continue
        write(6,'(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11)') 'Row',
     $  ' Column 5  ',' Column 6  ',' Column 7  ',' Column 8  '
        write(6,'(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11)') '---',
     $  '-----------','-----------','-----------','-----------'
        do 304 i=5,8
          write(6,'(2x,I1,2x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    i,(ex(i,j),j=5,8)
 304    continue
      endif
C
C... Output Assembled Stiffness
      if (output) then
        write(6,*)
        write(6,*) 'Assembled Stiffness (local)'
        write(6,'(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11)') 'Row',
     $  ' Column 1  ',' Column 2  ',' Column 3  ',' Column 4  '
        write(6,'(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11)') '---',
     $  '-----------','-----------','-----------','-----------'
        do 401 i=1,4
          write(6,'(2x,I1,2x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    i,(sm(i,j),j=1,4)
 401    continue
        write(6,'(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11)') 'Row',
     $  ' Column 5  ',' Column 6  ',' Column 7  ',' Column 8  '
        write(6,'(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11)') '---',
     $  '-----------','-----------','-----------','-----------'
        do 402 i=1,4
          write(6,'(2x,I1,2x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    i,(sm(i,j),j=5,8)
 402    continue
        write(6,'(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11)') 'Row',
     $  ' Column 1  ',' Column 2  ',' Column 3  ',' Column 4  '
        write(6,'(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11)') '---',
     $  '-----------','-----------','-----------','-----------'
        do 403 i=5,8
          write(6,'(2x,I1,2x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    i,(sm(i,j),j=1,4)
 403    continue
        write(6,'(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11)') 'Row',
     $  ' Column 5  ',' Column 6  ',' Column 7  ',' Column 8  '
        write(6,'(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11)') '---',
     $  '-----------','-----------','-----------','-----------'
        do 404 i=5,8
          write(6,'(2x,I1,2x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    i,(sm(i,j),j=5,8)
 404    continue
      endif
C
C Transform to Global Frame
C
C     Kshear = T^t sm T
C
      do 5000  j = 1,4
        jx  = ls(j)
        jy  = ls(j+4)
        jgx = lsg(j)
        jgy = lsg(j+4)
        jgz = lsg(j+8)
        do 5010  i = 1,4
          ix  = ls(i)
          iy  = ls(i+4)
          igx = lsg(i)
          igy = lsg(i+4)
          igz = lsg(i+8)
C
          Kshear(jgx,igx) = Kshear(jgx,igx) 
     $      + sm(jx,ix)*v1(1)*v1(1)
     $      + sm(jy,ix)*v2(1)*v1(1)
     $      + sm(jx,iy)*v1(1)*v2(1)
     $      + sm(jy,iy)*v2(1)*v2(1)
          Kshear(jgx,igy) = Kshear(jgx,igy) 
     $      + sm(jx,ix)*v1(1)*v1(2)
     $      + sm(jy,ix)*v2(1)*v1(2)
     $      + sm(jx,iy)*v1(1)*v2(2)
     $      + sm(jy,iy)*v2(1)*v2(2)
          Kshear(jgx,igz) = Kshear(jgx,igz) 
     $      + sm(jx,ix)*v1(1)*v1(3)
     $      + sm(jy,ix)*v2(1)*v1(3)
     $      + sm(jx,iy)*v1(1)*v2(3)
     $      + sm(jy,iy)*v2(1)*v2(3)
          Kshear(jgy,igx) = Kshear(jgy,igx) 
     $      + sm(jx,ix)*v1(2)*v1(1)
     $      + sm(jy,ix)*v2(2)*v1(1)
     $      + sm(jx,iy)*v1(2)*v2(1)
     $      + sm(jy,iy)*v2(2)*v2(1)
          Kshear(jgy,igy) = Kshear(jgy,igy) 
     $      + sm(jx,ix)*v1(2)*v1(2)
     $      + sm(jy,ix)*v2(2)*v1(2)
     $      + sm(jx,iy)*v1(2)*v2(2)
     $      + sm(jy,iy)*v2(2)*v2(2)
          Kshear(jgy,igz) = Kshear(jgy,igz) 
     $      + sm(jx,ix)*v1(2)*v1(3)
     $      + sm(jy,ix)*v2(2)*v1(3)
     $      + sm(jx,iy)*v1(2)*v2(3)
     $      + sm(jy,iy)*v2(2)*v2(3)
          Kshear(jgz,igx) = Kshear(jgz,igx) 
     $      + sm(jx,ix)*v1(3)*v1(1)
     $      + sm(jy,ix)*v2(3)*v1(1)
     $      + sm(jx,iy)*v1(3)*v2(1)
     $      + sm(jy,iy)*v2(3)*v2(1)
          Kshear(jgz,igy) = Kshear(jgz,igy) 
     $      + sm(jx,ix)*v1(3)*v1(2)
     $      + sm(jy,ix)*v2(3)*v1(2)
     $      + sm(jx,iy)*v1(3)*v2(2)
     $      + sm(jy,iy)*v2(3)*v2(2)
          Kshear(jgz,igz) = Kshear(jgz,igz) 
     $      + sm(jx,ix)*v1(3)*v1(3)
     $      + sm(jy,ix)*v2(3)*v1(3)
     $      + sm(jx,iy)*v1(3)*v2(3)
     $      + sm(jy,iy)*v2(3)*v2(3)
C
 5010   continue
 5000 continue
C
C... Output Final Stiffness
      if (output) then
        write(6,*)
        write(6,*) 'Assembled Stiffness (global)'
        write(6,
     $  '(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11)')
     $  'Row',' Column 1  ',' Column 2  ',' Column 3  '
     $       ,' Column 4  ',' Column 5  ',' Column 6  '
        write(6,
     $  '(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11)')
     $  '---','-----------','-----------','-----------'
     $       ,'-----------','-----------','-----------'
        do 501 i=1,6
          write(6,
     $'(1x,I2,2x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    i,(Kshear(i,j),j=1,6)
 501    continue
        write(6,
     $  '(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11)')
     $  'Row',' Column 7  ',' Column 8  ',' Column 9  '
     $       ,' Column 10 ',' Column 11 ',' Column 12 '
        write(6,
     $  '(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11)')
     $  '---','-----------','-----------','-----------'
     $       ,'-----------','-----------','-----------'
        do 502 i=1,6
          write(6,
     $'(1x,I2,2x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    i,(Kshear(i,j),j=7,12)
 502    continue
        write(6,
     $  '(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11)')
     $  'Row',' Column 1  ',' Column 2  ',' Column 3  '
     $       ,' Column 4  ',' Column 5  ',' Column 6  '
        write(6,
     $  '(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11)')
     $  '---','-----------','-----------','-----------'
     $       ,'-----------','-----------','-----------'
        do 503 i=7,12
          write(6,
     $'(1x,I2,2x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    i,(Kshear(i,j),j=1,6)
 503    continue
        write(6,
     $  '(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11)')
     $  'Row',' Column 7  ',' Column 8  ',' Column 9  '
     $       ,' Column 10 ',' Column 11 ',' Column 12 '
        write(6,
     $  '(1x,A3,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11,1x,A11)')
     $  '---','-----------','-----------','-----------'
     $       ,'-----------','-----------','-----------'
        do 504 i=7,12
          write(6,
     $'(1x,I2,2x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    i,(Kshear(i,j),j=7,12)
 504    continue
      endif
      return
      end
C
