C--------------------------------------------------------------------
C    Given the node displacement,computes
C    corner stresses on a eight-node hexahedron.
C
C    The stresses at the mid-nodes are interpolated
C
C    The input arguments are
C
C      x      (20x1)  array of x co-ordinates of hexahedron nodes.
C      y      (20x1)  array of y co-ordinates of hexahedron nodes.
C      z      (20x1)  array of z co-ordinates of hexahedron nodes.
C      c      (6x6)  stress-strain constitutive matrix.
C      u      (60x1) array of element node displacement arranged in
C             ux1,uy1,uz1,ux2,uy2,uz2, ............., ux20,uy20,uz20
C
C   The outputs are
C
C
C      ESIG   (6x8) array of corner node stresses arranged
C             sigxx1, sigyy1, sigzz1, tauxy1, tauyz1, tauxz1,
C             sigxx2, sigyy2, sigzz2, tauxy2, tauyz2, tauxz2,
C              ..    ..   ..    ..   .. ...    ...   ...  ..
C              ..    ..   ..    ..   ..       tauyz8, tauxz8
C---------------------------------------------------------------------
          subroutine sands20(elm,x,y,z,c,u,stress,strain,
     &                       maxgus,maxstr,maxsze,vmflg,strainFlg)
C
C                A R G U M E N T S
C
C        implicit none
         double precision  x(20), y(20), z(20), c(6,6), u(60)
         real*8  stress(maxsze,maxstr,maxgus) 
         real*8  strain(maxsze,maxstr,maxgus)
         integer elm,maxgus,maxstr,maxsze
         logical vmflg,strainFlg
C
C                L O C A L  V A R I A B L E S
C
         double precision  q(20), qx(20), qy(20), qz(20)
         double precision  xicorn(8), etacorn(8), mucorn(8)
         double precision  xi, eta, mu, det, epsxx, epsyy
         double precision  epszz, gamxy, gamyz, gamxz
         double precision  half 
C          
         integer           ipolm(12),ipoll(12),ipolr(12)
         integer           i, n, nm, nl, nr
C
C                D A T A
C
         data    xicorn  /-1.0,  1.0,  1.0, -1.0, -1.0,  1.0, 1.0, -1.0/
         data    etacorn /-1.0, -1.0,  1.0,  1.0, -1.0, -1.0, 1.0,  1.0/
         data    mucorn  /-1.0, -1.0, -1.0, -1.0,  1.0,  1.0, 1.0,  1.0/
c
         data    ipolm   / 9,10,11,12,13,14,15,16,17,18,19,20/
         data    ipoll   / 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4/
         data    ipolr   / 2, 3, 4, 1, 6, 7, 8, 5, 5, 6, 7, 8/
c
         half = 0.5d0
C
C                L O G I C
C
C
         do  2000  n = 1,8
            xi  = xicorn(n)
            eta = etacorn(n)
            mu  = mucorn(n)
            call h20shpe (xi, eta, mu, x, y, z, q, qx, qy, qz, det)
            if ( det .le. 0.0) then
              write(6,*) 'Negative  Jacobian determinant'
            if ( det .eq. 0.0) then
              write(6,*) 'Zero Jacobian determinant'
            endif
          return
         endif
         epsxx = 0.0d0
         epsyy = 0.0d0
         epszz = 0.0d0
         gamxy = 0.0d0
         gamyz = 0.0d0
         gamxz = 0.0d0
         do 1500 i = 1,20
            epsxx =   epsxx + qx(i)*u(3*i-2)
            epsyy =   epsyy + qy(i)*u(3*i-1) 
            epszz =   epszz + qz(i)*u(3*i  )  
            gamxy =   gamxy + qy(i)*u(3*i-2) + qx(i)*u(3*i-1)
            gamyz =   gamyz + qz(i)*u(3*i-1) + qy(i)*u(3*i  )
            gamxz =   gamxz + qx(i)*u(3*i  ) + qz(i)*u(3*i-2)
 1500       continue

         stress(elm,1,n) = c(1,1)*epsxx + c(1,2)*epsyy + c(1,3)*epszz
     $                   + c(1,4)*gamxy + c(1,5)*gamyz + c(1,6)*gamxz
         stress(elm,2,n) = c(2,1)*epsxx + c(2,2)*epsyy + c(2,3)*epszz
     $                   + c(2,4)*gamxy + c(2,5)*gamyz + c(2,6)*gamxz
         stress(elm,3,n) = c(3,1)*epsxx + c(3,2)*epsyy + c(3,3)*epszz
     $                   + c(3,4)*gamxy + c(3,5)*gamyz + c(3,6)*gamxz
         stress(elm,4,n) = c(4,1)*epsxx + c(4,2)*epsyy + c(4,3)*epszz
     $                   + c(4,4)*gamxy + c(4,5)*gamyz + c(4,6)*gamxz
         stress(elm,5,n) = c(5,1)*epsxx + c(5,2)*epsyy + c(5,3)*epszz
     $                   + c(5,4)*gamxy + c(5,5)*gamyz + c(5,6)*gamxz
         stress(elm,6,n) = c(6,1)*epsxx + c(6,2)*epsyy + c(6,3)*epszz
     $                   + c(6,4)*gamxy + c(6,5)*gamyz + c(6,6)*gamxz

         strain(elm,1,n) = epsxx
         strain(elm,2,n) = epsyy 
         strain(elm,3,n) = epszz 
         strain(elm,4,n) = gamxy 
         strain(elm,5,n) = gamyz 
         strain(elm,6,n) = gamxz 
C
 2000    continue
C
C.... COMPUTE THE VON MISES STRESS
C
        if (vmflg) then
          call vmelmv(stress,maxgus,maxstr,maxsze,elm,8)
        endif
C
C.... COMPUTE EQUIVALENT STRAIN
C
        if(strainFlg) then
          call strainvm(strain,maxgus,maxstr,maxsze,8)
        endif
C
C.... Interpolate stresses
C
      do i=1,12
c
        nm=ipolm(i)
        nl=ipoll(i)
        nr=ipolr(i)
c
        stress(elm,1,nm) = half*(stress(elm,1,nl)+stress(elm,1,nr))
        stress(elm,2,nm) = half*(stress(elm,2,nl)+stress(elm,2,nr))
        stress(elm,3,nm) = half*(stress(elm,3,nl)+stress(elm,3,nr))
        stress(elm,4,nm) = half*(stress(elm,4,nl)+stress(elm,4,nr))
        stress(elm,5,nm) = half*(stress(elm,5,nl)+stress(elm,5,nr))
        stress(elm,6,nm) = half*(stress(elm,6,nl)+stress(elm,6,nr))
        stress(elm,7,nm) = half*(stress(elm,7,nl)+stress(elm,7,nr))
c
        strain(elm,1,nm) = half*(strain(elm,1,nl)+strain(elm,1,nr))
        strain(elm,2,nm) = half*(strain(elm,2,nl)+strain(elm,2,nr))
        strain(elm,3,nm) = half*(strain(elm,3,nl)+strain(elm,3,nr))
        strain(elm,4,nm) = half*(strain(elm,4,nl)+strain(elm,4,nr))
        strain(elm,5,nm) = half*(strain(elm,5,nl)+strain(elm,5,nr))
        strain(elm,6,nm) = half*(strain(elm,6,nl)+strain(elm,6,nr))
        strain(elm,7,nm) = half*(strain(elm,7,nl)+strain(elm,7,nr))
c
      enddo
c
      return
      end
C
