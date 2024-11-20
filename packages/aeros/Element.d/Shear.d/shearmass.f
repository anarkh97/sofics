C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Gregrory W. Brown
C January, 2000
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     SHEARMASS forms the element mass matrix of a
C     four-node shear panel element.
C
C     The calling sequence is
C
C       call   shearmass (xg,yg,zg,t,rho,p,Kmass,m,vol,gamma,
C                         grvfor,grvflg,totmas,masflg)
C
C     where the input arguments are
C
C       xg        (4 x 1) array of global x coordinates
C       yg        (4 x 1) array of global y coordinates
C       zg        (4 x 1) array of global y coordinates
C       t         thickness of quadrilateral nodes
C       rho       density of element
C       p         Gauss quadrature rule (no. of points)
C       m         First dimension of SM in calling program.
C       gamma     Accelerations due to gravity vector
C       grvflg    Flag to indicate gravity forces needed
C       masflg    Flag for subdomain mass computation

C
C     The outputs are:
C
C       Kshear    (12x12) computed element stiffness matrix.  The
C                 arrangement of rows and columns pertains to node
C                 displacements arranged in the order
C                  (vx1, vy1, yz1, vx2, ... vz4)
C                 This particular ordering is obtained by setting array
C                 LST to  1,4,7,10,2,5,8,11,3,6,9,12 
C                 (see DATA statement below)
C
C       vol       Volume of element
C       grvfor    Body forces due to gravity
C       totmas    Accumulated total mass
C
C The local control flag global determines whether the
C Mass matrix is built in the local or global frame
C
C if global = .false.
C   The mass matrix is built in plane and rotated to the global
C   frame.  This avoids mass in the directions in which the
C   shear panel has no stiffness, but could result in different
C   mass entries in each of the global coordinate directions.
C   This also can create off-diagonal entries, which is
C   incompatible with a lumped mass assumption
C   Until the software supports consistant mass, this must
C   be left set to .true.  
C if global = .true.  <------- Must be this way.
C   The mass matrix is built in the global frame.
C   frame.  This avoids different mass entries in each of the global
C   coordinate directions, but could result in mass in directions
C   in which the shear panel has no stiffness.
C
      subroutine    shearmass (xg,yg,zg,t,rho,p,Kmass,m,vol,gamma,
     $                         grvfor,grvflg,totmas,masflg)
C
C                   A R G U M E N T S
C
      integer p, m
      real*8  xg(*), yg(*), zg(*), t, rho
      real*8  Kmass(m,*), vol
      real*8  totmas, grvfor(*), gamma(*)
      logical grvflg, masflg 
C
C                   L O C A L   V A R I A B L E S
C
C Frame in which to build mass
      logical global
C Local Coordinates
      real*8  x(4), y(4), z(4)
C Working Vectors
      real*8  node(3,4), local(3,4) ,v1(3),v2(3),v3(3),v4(3), len
C In-plane Mass Matrix
      real*8  mm(8)
C Output Control
      logical output
C Other Variables
      real*8  nmass, emass
      real*8  q(4), qx(4), qy(4)
      real*8  xi, eta, det, weight
      integer i, j, jx, jy, k, l, ls(8)
      integer lsg(12), jgx, jgy, jgz
C
C                   D A T A
C
      data    ls  /1,3,5,7,2,4,6,8/
      data    lsg /1,4,7,10,2,5,8,11,3,6,9,12/
C
C Set Output Status
      output = .false.
C
C Set Frame in Which to Build Mass
C LEAVE THIS SET TO TRUE
      global = .true.
C
C                   L O G I C
C
      do 1200  j = 1,12
        do 1100  i = 1,12
          Kmass(i,j) = 0.0
 1100   continue
 1200 continue
      do 1201  j = 1,8
        mm(j) = 0.0
 1201 continue
C
C... Output Properties of Element
      if (output) then
        write(6,*)
        write(6,*) 'SHEARMASS.F'
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
        write(6,'(1x,A17,E15.8)') 'Density        = ',rho
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
        write(6,*)  'ERROR: subroutine shearmass.f'
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
          write(6,*)  'WARNING: subroutine shearmass.f'
          write(6,*)  'Element has warp > 10%'
          write(6,*)  'Warp = ',len
        endif
 1340 continue
C
C... Compute Volume of Element
C
      vol = 0.0
      do 3000  k = 1,p
        do 2500  l = 1,p
          call     QGAUSS (p, k, p, l, xi, eta, weight)
          call     Q4SHPE (xi, eta, x, y, q, qx, qy, det)
          if (det .le. 0.0)        then
            write(6,*)  'ERROR: subroutine shearmass.f'
            write(6,*)  'Negative Jacobian determinant in 4 node quad'
            if (det .eq. 0.0)      then
              write(6,*)  'ERROR: subroutine shearmass.f'
              write(6,*) 'Zero Jacobian determinant in 4 node quad'
            end if
            stop
          end if
          vol = vol + det
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
C... Determine total mass of the element
      emass = vol*rho
      if (output) then
        write(6,*)
        write(6,'(1x,A17,E15.8)') 'Element Mass   = ',emass
      endif
C
C... Determine nodal mass
      nmass = emass/4.0
C
C... Build Mass Matrix in Local Frame
      if (.not.global) then
        do 3500 i = 1,8
          mm(i) = nmass
 3500   continue
      endif
C
C... Output Mass Matrix
      if (output.and.(.not.global)) then
        write(6,*)
        write(6,*) 'Mass Matrix (local)'
        write(6,*) '(Diagonal Lumped Matrix)'
        write(6,'(1x,A3,1x,A3,1x,A11)') 'Row','Col','   Value   '
        write(6,'(1x,A3,1x,A3,1x,A11)') '---','---','-----------'
        do 201 i = 1,8
          write(6,'(1x,I3,1x,I3,1x,E11.4)') i,i,mm(i)
 201    continue
      endif
C
C if (global)
C Build in Global Frame
C
      if (global) then
        do 4000 i = 1,12
          Kmass(i,i) = nmass
 4000   continue
C
C else
C Transform to Global Frame
C
C     Kmass = T^t mm T
C
      else
        do 5000  j = 1,4
          jx  = ls(j)
          jy  = ls(j+4)
          jgx = lsg(j)
          jgy = lsg(j+4)
          jgz = lsg(j+8)
C
          Kmass(jgx,jgx) = Kmass(jgx,jgx) 
     $      + mm(jx)*v1(1)*v1(1)
     $      + mm(jy)*v2(1)*v2(1)
          Kmass(jgx,jgy) = Kmass(jgx,jgy) 
     $      + mm(jx)*v1(1)*v1(2)
     $      + mm(jy)*v2(1)*v2(2)
          Kmass(jgx,jgz) = Kmass(jgx,jgz) 
     $      + mm(jx)*v1(1)*v1(3)
     $      + mm(jy)*v2(1)*v2(3)
          Kmass(jgy,jgx) = Kmass(jgy,jgx) 
     $      + mm(jx)*v1(2)*v1(1)
     $      + mm(jy)*v2(2)*v2(1)
          Kmass(jgy,jgy) = Kmass(jgy,jgy) 
     $      + mm(jx)*v1(2)*v1(2)
     $      + mm(jy)*v2(2)*v2(2)
          Kmass(jgy,jgz) = Kmass(jgy,jgz) 
     $      + mm(jx)*v1(2)*v1(3)
     $      + mm(jy)*v2(2)*v2(3)
          Kmass(jgz,jgx) = Kmass(jgz,jgx) 
     $      + mm(jx)*v1(3)*v1(1)
     $      + mm(jy)*v2(3)*v2(1)
          Kmass(jgz,jgy) = Kmass(jgz,jgy) 
     $      + mm(jx)*v1(3)*v1(2)
     $      + mm(jy)*v2(3)*v2(2)
          Kmass(jgz,jgz) = Kmass(jgz,jgz) 
     $      + mm(jx)*v1(3)*v1(3)
     $      + mm(jy)*v2(3)*v2(3)
 5000   continue
      endif
C
C... Output Final Mess
      if (output) then
        write(6,*)
        write(6,*) 'Assembled Mass (global)'
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
     $    i,(Kmass(i,j),j=1,6)
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
     $    i,(Kmass(i,j),j=7,12)
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
     $    i,(Kmass(i,j),j=1,6)
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
     $    i,(Kmass(i,j),j=7,12)
 504    continue
      endif
C
C.... Compute the body force due to gravity if needed
C
      if (grvflg) then
        grvfor(1) = emass*gamma(1)
        grvfor(2) = emass*gamma(2)
        grvfor(3) = emass*gamma(3)
      endif
C
C... Compute the Subdomain total mass
C
      if (masflg) then
        totmas = totmas + emass
      endif
C
      return
      end
C
