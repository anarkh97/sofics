C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C Gregrory W. Brown
C January, 2000
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     SPSTRESS calculates stresses and strains for the
C     four-node shear panel element.
C
C            call spstress (xg,yg,zg,v,G,E,F1,F2,
C    $                      stress,strain,vmssig,vmseps)
C
C     where the input arguments are
C
C       xg        (4 x 1) array of global x coordinates
C       yg        (4 x 1) array of global y coordinates
C       zg        (4 x 1) array of global y coordinates
C       v         (12x 1) array of node dispalcements
C                 arranged, vx1,vy1,vz1,...,vz4
C       G         Shear Modulus of element
C       E         Youngs Modulus of element
C       F1        Extensional factor 1
C       F2        Extensional factor 2
C
C     The outputs are:
C
C       stress    (6x4) corner node stresses arranged
C                  (sigxx, sigyy, sigzz, tauxy, tauyz, tauxz)
C       strain    (6x4) corner node strains 
C       vmssig    Vommises stress
C       vmseps    Vommises strain
C
      subroutine  spstress (xg,yg,zg,v,G,E,F1,F2,
     $                      stress,strain,vmssig,vmseps)
C
C                   A R G U M E N T S
C
      real*8  xg(*), yg(*), zg(*), v(*), G, E, F1, F2
      real*8  stress(6,4), strain(6,4), vmssig, vmseps
C
C                   L O C A L   V A R I A B L E S
C
C Local Coordinates
      real*8  x(4), y(4), z(4), vl(8), zdl(4)
C Working Vectors
      real*8  node(3,4), local(3,4) ,v1(3),v2(3),v3(3),v4(3), len
C Consitutive Matrix
      real*8  c(3,3)
C For sands2.f routine
      character escm
      real*8    quadstress(1,7,4),quadstrain(1,7,4)
      integer   four, three, one, seven
      logical   falseflag
C Effective Extentional Areas
      real*8  area(4)
C Rod Coordinates and Stiffness
      real*8  epsrod,  rodstrain(6,4)
      real*8  xr(2), yr(2), zr(2), vrx(2), vry(2)
      integer rnod(2,4) 
      real*8  l1, l2, q1, q2, dq
C For vms routine
      real*8  vmsstress(1,7,4), vmsstrain(1,7,4)
C Output Control
      logical output
C Other Variables
      integer i, j, k
      real*8  avg
C
C Set Output Status
      output = .false.
C
C                   L O G I C
C
C... Output Properties of Element
      if (output) then
        write(6,*)
        write(6,*) 'SPSTRESS.F'
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
        write(6,*) 'Node Displacements (global)'
        write(6,'(1x,A4,1x,A15,1x,A15,1x,A15)') 'Node',
     $    'X Displacement ','Y Displacement ','Z Displacement '
        write(6,'(1x,A4,1x,A15,1x,A15,1x,A15)') '----',
     $    '---------------','---------------','--------------'
        do 101 i=1,4
          j = (3*i)-3
          write(6,'(3x,I1,2x,E15.8,1x,E15.8,1x,E15.8)')
     $      i,v(j+1),v(j+2),v(j+3)
 101    continue
        write(6,*)
        write(6,'(1x,A17,E15.8)') 'Youngs Modulus = ',E
        write(6,'(1x,A17,E15.8)') 'Shear Modulus  = ',G
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
        write(6,*)  'ERROR: subroutine spstress.f'
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
        zdl(i) = 0.0
        do 1315 j = 1,3
          local(j,i) = 0.0
 1315   continue
 1310 continue
      do 1311 i = 1,8
        vl(i) = 0.0
 1311 continue
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
 1320 continue
C
C... Transform displamcents
C... vl = T*v
C
      do 1321 i = 1,4
        do 1326 j =1,3
          vl((2*i)-1) = vl((2*i)-1) + v1(j)*v((3*i)-3+j)
          vl((2*i))   = vl((2*i))   + v2(j)*v((3*i)-3+j)
          zdl(i)      = zdl(i)      + v3(j)*v((3*i)-3+j)
 1326   continue
 1321 continue
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
        write(6,*)
        write(6,*) 'Node Displacements (local)'
        write(6,'(1x,A4,1x,A15,1x,A15,1x,A15)') 'Node',
     $    'X Displacement ','Y Displacement ','Z Displacement '
        write(6,'(1x,A4,1x,A15,1x,A15,1x,A15)') '----',
     $    '---------------','---------------','--------------'
        do 111 i=1,4
          j = (2*i)-2
          write(6,'(3x,I1,2x,E15.8,1x,E15.8,1x,E15.8)')
     $      i,vl(j+1),vl(j+2),zdl(i)
 111    continue
        write(6,*) '(Z Components ignored in stress calculation)'
      endif
C
C... Check amount of out-of-plane behaviore
C    Use x coordinate of node 2 as representative length
C
      do 1340 i = 1,4
        len = abs(local(3,i)/local(1,2))
        if (len.gt.0.1) then
          write(6,*)  'WARNING: subroutine spstress.f'
          write(6,*)  'Element has warp > 10%'
          write(6,*)  'Warp = ',len
        endif
 1340 continue
C
C... Set Parameters for sands2.f
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
      escm = 'E' 
C
      seven = 7
      four  = 4
      three = 3
      one   = 1
      falseflag = .false.
      do 3000  i = 1,7
        do 2500  j = 1,4
          quadstress(1,i,j) = 0.
          quadstrain(1,i,j) = 0.
 2500   continue
 3000 continue
C
C... Get Stress for Shear Portion
      call sands2(escm,x,y,c,vl,quadstress,quadstrain,four,
     &                  seven,one,one,falseflag,falseflag,0,0,0)
C
C... Output Stresses for Shear Portion
      if (output) then
        write(6,*)
        write(6,*) 'Stresses from Shear (local)'
        write(6,'(1x,A6,1x,A11,1x,A11,1x,A11,1x,A11)') 'Stress',
     $  '  Node 1   ','  Node 2   ','  Node 3   ','  Node 4   '
        write(6,'(1x,A6,1x,A11,1x,A11,1x,A11,1x,A11)') '------',
     $  '-----------','-----------','-----------','-----------'
          write(6,'(1x,A6,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    'sigxx ',(quadstress(1,1,j),j=1,4)
          write(6,'(1x,A6,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    'sigyy ',(quadstress(1,2,j),j=1,4)
          write(6,'(1x,A6,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    'tauxy ',(quadstress(1,4,j),j=1,4)
        write(6,*)
        write(6,*) 'Strains from Shear (local)'
        write(6,'(1x,A6,1x,A11,1x,A11,1x,A11,1x,A11)') 'Strain',
     $  '  Node 1   ','  Node 2   ','  Node 3   ','  Node 4   '
        write(6,'(1x,A6,1x,A11,1x,A11,1x,A11,1x,A11)') '------',
     $  '-----------','-----------','-----------','-----------'
          write(6,'(1x,A6,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    'exx   ',(quadstrain(1,1,j),j=1,4)
          write(6,'(1x,A6,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    'eyy   ',(quadstrain(1,2,j),j=1,4)
          write(6,'(1x,A6,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    'exy   ',(quadstrain(1,4,j),j=1,4)
      endif
C
C... Compute Extensional Stresses
C    This is defined as rods along the edges
C
      avg = 1.0
      do 4520 i = 1,4
        area(i) = 0.
 4520  continue
      if (F1.gt.0.) then
        avg = avg+1.0
        area(1) = 1.
        area(2) = 1.
      endif
      if (F2.gt.0.) then
        avg = avg+1.0
        area(3) = 1.
        area(4) = 1.
      endif
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
      do 4500 i = 1,6
        do 4510 j = 1,4
          rodstrain(i,j) = 0.
 4510   continue
 4500 continue
C
C Begin Loop Over 4 Rods
      do 4000 k = 1,4
        if (area(k).ne.0.) then
C
C Get nodes of rod
C and put origin at 1st node
          do 4010 i = 1,2
            xr(i) = x(rnod(i,k)) - x(rnod(1,k))
            yr(i) = y(rnod(i,k)) - y(rnod(1,k))
 4010     continue
C
C Length of Rod
          len =  xr(2)**2 + yr(2)**2
          if (len.le.0.) then
            write(6,*)  'ERROR: subroutine spstress.f'
            write(6,*)  'Rod has zero or negative length.'
            stop
          endif
          len = abs(sqrt(len))
C
C Get deformations of rod
          do 4020 i = 1,2
            vrx(i) = vl(2*rnod(i,k)-1)
            vry(i) = vl(2*rnod(i,k))
 4020     continue
C
C Stress Computation copied from sands1.f by P.R. Stern
C
C.... Calculate the defromation difference
C
          l1 = xr(2)/len
          l2 = yr(2)/len
C
          q1 = l1*vrx(1) + l2*vry(1)
          q2 = l1*vrx(2) + l2*vry(2)
C
          dq = q2-q1
C
C.... Calculate the Cauchy Strains
          epsrod = dq/len
C
C.... Add to strain matrix for each node
C     Transformation Matrix is  T = |  cos  sin |
C                                   | -sin  cos |
C    strain = T^t * | epsrod  0 | *T
C                   |   0     0 |   
C cos is stored in l1
C sin is stored in l1
C
          rodstrain(1,rnod(1,k)) = rodstrain(1,rnod(1,k)) + epsrod*l1*l1
          rodstrain(2,rnod(1,k)) = rodstrain(2,rnod(1,k)) + epsrod*l2*l2
          rodstrain(4,rnod(1,k)) = rodstrain(4,rnod(1,k)) + epsrod*l1*l2
          rodstrain(1,rnod(2,k)) = rodstrain(1,rnod(2,k)) + epsrod*l1*l1
          rodstrain(2,rnod(2,k)) = rodstrain(2,rnod(2,k)) + epsrod*l2*l2
          rodstrain(4,rnod(2,k)) = rodstrain(4,rnod(2,k)) + epsrod*l1*l2
C
        endif
C
C End Loop on Rods
 4000 continue
C
C Compute Stresses From Strains
C Add Stresses/Strains to Total from Shear
C
      do 4530 i = 1,6
        do 4540 j = 1,4
          quadstress(1,i,j) =  quadstress(1,i,j) + E*rodstrain(i,j)
          quadstrain(1,i,j) =  quadstrain(1,i,j) + rodstrain(i,j)
 4540   continue
 4530 continue
C
C Average Stresses
      do 4550 i = 1,6
        do 4560 j = 1,4
          quadstress(1,i,j) =  quadstress(1,i,j)/avg
          quadstrain(1,i,j) =  quadstrain(1,i,j)/avg
 4560   continue
 4550 continue
C
C... Output Stresses
      if (output) then
        write(6,*)
        write(6,*) 'Net Stresses (local)'
        write(6,'(1x,A6,1x,A11,1x,A11,1x,A11,1x,A11)') 'Stress',
     $  '  Node 1   ','  Node 2   ','  Node 3   ','  Node 4   '
        write(6,'(1x,A6,1x,A11,1x,A11,1x,A11,1x,A11)') '------',
     $  '-----------','-----------','-----------','-----------'
          write(6,'(1x,A6,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    'sigxx ',(quadstress(1,1,j),j=1,4)
          write(6,'(1x,A6,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    'sigyy ',(quadstress(1,2,j),j=1,4)
          write(6,'(1x,A6,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    'tauxy ',(quadstress(1,4,j),j=1,4)
        write(6,*)
        write(6,*) 'Net Strains (local)'
        write(6,'(1x,A6,1x,A11,1x,A11,1x,A11,1x,A11)') 'Strain',
     $  '  Node 1   ','  Node 2   ','  Node 3   ','  Node 4   '
        write(6,'(1x,A6,1x,A11,1x,A11,1x,A11,1x,A11)') '------',
     $  '-----------','-----------','-----------','-----------'
          write(6,'(1x,A6,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    'exx   ',(quadstrain(1,1,j),j=1,4)
          write(6,'(1x,A6,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    'eyy   ',(quadstrain(1,2,j),j=1,4)
          write(6,'(1x,A6,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    'exy   ',(quadstrain(1,4,j),j=1,4)
      endif
C
C Transform to Global Frame
C
C     stress = T^t stress_local T
C
C each part of quadstress holds [sigxx,sigyy,tauxy]
C each part of quadstress holds [sigxx,sigyy,sigzz,tauxy,tauyz,tauxz]
C
      do 5000  j = 1,4
C
C Global XX Strain
        strain(1,j) =
     $      + quadstrain(1,1,j)*v1(1)*v1(1)
     $      + quadstrain(1,4,j)*v2(1)*v1(1)
     $      + quadstrain(1,4,j)*v1(1)*v2(1)
     $      + quadstrain(1,2,j)*v2(1)*v2(1)
C Global YY Strain
        strain(2,j) =
     $      + quadstrain(1,1,j)*v1(2)*v1(2)
     $      + quadstrain(1,4,j)*v2(2)*v1(2)
     $      + quadstrain(1,4,j)*v1(2)*v2(2)
     $      + quadstrain(1,2,j)*v2(2)*v2(2)
C Global ZZ Strain
        strain(3,j) =
     $      + quadstrain(1,1,j)*v1(3)*v1(3)
     $      + quadstrain(1,4,j)*v2(3)*v1(3)
     $      + quadstrain(1,4,j)*v1(3)*v2(3)
     $      + quadstrain(1,2,j)*v2(3)*v2(3)
C Global XY Strain
        strain(4,j) =
     $      + quadstrain(1,1,j)*v1(1)*v1(2)
     $      + quadstrain(1,4,j)*v2(1)*v1(2)
     $      + quadstrain(1,4,j)*v1(1)*v2(2)
     $      + quadstrain(1,2,j)*v2(1)*v2(2)
C Global YZ Strain
        strain(5,j) =
     $      + quadstrain(1,1,j)*v1(2)*v1(3)
     $      + quadstrain(1,4,j)*v2(2)*v1(3)
     $      + quadstrain(1,4,j)*v1(2)*v2(3)
     $      + quadstrain(1,2,j)*v2(2)*v2(3)
C Global XZ Strain
        strain(6,j) =
     $      + quadstrain(1,1,j)*v1(1)*v1(3)
     $      + quadstrain(1,4,j)*v2(1)*v1(3)
     $      + quadstrain(1,4,j)*v1(1)*v2(3)
     $      + quadstrain(1,2,j)*v2(1)*v2(3)
C
C Global XX Stress
        stress(1,j) =
     $      + quadstress(1,1,j)*v1(1)*v1(1)
     $      + quadstress(1,4,j)*v2(1)*v1(1)
     $      + quadstress(1,4,j)*v1(1)*v2(1)
     $      + quadstress(1,2,j)*v2(1)*v2(1)
C Global YY Stress
        stress(2,j) =
     $      + quadstress(1,1,j)*v1(2)*v1(2)
     $      + quadstress(1,4,j)*v2(2)*v1(2)
     $      + quadstress(1,4,j)*v1(2)*v2(2)
     $      + quadstress(1,2,j)*v2(2)*v2(2)
C Global ZZ Stress
        stress(3,j) =
     $      + quadstress(1,1,j)*v1(3)*v1(3)
     $      + quadstress(1,4,j)*v2(3)*v1(3)
     $      + quadstress(1,4,j)*v1(3)*v2(3)
     $      + quadstress(1,2,j)*v2(3)*v2(3)
C Global XY Stress
        stress(4,j) =
     $      + quadstress(1,1,j)*v1(1)*v1(2)
     $      + quadstress(1,4,j)*v2(1)*v1(2)
     $      + quadstress(1,4,j)*v1(1)*v2(2)
     $      + quadstress(1,2,j)*v2(1)*v2(2)
C Global YZ Stress
        stress(5,j) =
     $      + quadstress(1,1,j)*v1(2)*v1(3)
     $      + quadstress(1,4,j)*v2(2)*v1(3)
     $      + quadstress(1,4,j)*v1(2)*v2(3)
     $      + quadstress(1,2,j)*v2(2)*v2(3)
C Global XZ Stress
        stress(6,j) =
     $      + quadstress(1,1,j)*v1(1)*v1(3)
     $      + quadstress(1,4,j)*v2(1)*v1(3)
     $      + quadstress(1,4,j)*v1(1)*v2(3)
     $      + quadstress(1,2,j)*v2(1)*v2(3)
C
 5000 continue
C
C... Output Stresses
      if (output) then
        write(6,*)
        write(6,*) 'Net Stresses (global)'
        write(6,'(1x,A6,1x,A11,1x,A11,1x,A11,1x,A11)') 'Stress',
     $  '  Node 1   ','  Node 2   ','  Node 3   ','  Node 4   '
        write(6,'(1x,A6,1x,A11,1x,A11,1x,A11,1x,A11)') '------',
     $  '-----------','-----------','-----------','-----------'
          write(6,'(1x,A6,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    'sigxx ',(stress(1,j),j=1,4)
          write(6,'(1x,A6,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    'sigyy ',(stress(2,j),j=1,4)
          write(6,'(1x,A6,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    'sigzz ',(stress(3,j),j=1,4)
          write(6,'(1x,A6,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    'tauxy ',(stress(4,j),j=1,4)
          write(6,'(1x,A6,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    'tauyz ',(stress(5,j),j=1,4)
          write(6,'(1x,A6,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    'tauxz ',(stress(6,j),j=1,4)
        write(6,*)
        write(6,*) 'Net Strains (global)'
        write(6,'(1x,A6,1x,A11,1x,A11,1x,A11,1x,A11)') 'Strain',
     $  '  Node 1   ','  Node 2   ','  Node 3   ','  Node 4   '
        write(6,'(1x,A6,1x,A11,1x,A11,1x,A11,1x,A11)') '------',
     $  '-----------','-----------','-----------','-----------'
          write(6,'(1x,A6,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    'epsxx ',(strain(1,j),j=1,4)
          write(6,'(1x,A6,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    'epsyy ',(strain(2,j),j=1,4)
          write(6,'(1x,A6,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    'epszz ',(strain(3,j),j=1,4)
          write(6,'(1x,A6,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    'epsxy ',(strain(4,j),j=1,4)
          write(6,'(1x,A6,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    'epsyz ',(strain(5,j),j=1,4)
          write(6,'(1x,A6,1x,E11.4,1x,E11.4,1x,E11.4,1x,E11.4)')
     $    'epsxz ',(strain(6,j),j=1,4)
      endif
C
C... Compute the Von Mises Stress
C
      do 6000 j = 1,4
        vmsstress(1,7,j) = 0.
        vmsstrain(1,7,j) = 0.
        do 6010 i = 1,6
          vmsstress(1,i,j) = stress(i,j)
          vmsstrain(1,i,j) = strain(i,j)
 6010   continue
 6000 continue
      call   vmelmv(vmsstress,four,seven,one,one,four)
      call strainvm(vmsstrain,four,seven,one,four)
      vmssig = vmsstress(1,7,1)
      vmseps = vmsstrain(1,7,1)

C... Output Stresses
      if (output) then
        write(6,*)
        write(6,*) 'Von Mises Stress = ',vmssig
        write(6,*) 'Von Mises Strain = ',vmseps
      endif
C
      return
      end
C
