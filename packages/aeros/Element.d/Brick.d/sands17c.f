C=PURPOSE Compute corner complex stresses of 8-node hexa
C
C    Given the node displacement, SANDS17 computes
C    corner stresses on a eight-node hexahedron.
C
C
C   The  calling sequence is
C
C   where the input arguments are
C
C      x      (8x1)  array of x co-ordinates of hexahedron nodes.
C      y      (8x1)  array of y co-ordinates of hexahedron nodes.
C      z      (8x1)  array of z co-ordinates of hexahedron nodes.
C      c      (6x6)  stress-strain constitutive matrix.
C      v      (24x1) array of element node displacement arranged in
C             ux1,uy1,uz1,ux2,uy2,uz2, ............., uz8
C
C   The outputs are
C
C
C  STRESS     sigxx1, sigyy1, sigzz1, tauxy1, tauyz1, tauxz1,
C             sigxx2, sigyy2, sigzz2, tauxy2, tauyz2, tauxz2,
C              ..    ..   ..    ..   .. ...    ...   ...  ..
C              ..    ..   ..    ..   ..       tauyz8, tauxz8
C
         subroutine sands17c(elm,x,y,z,c,v,stress,strain,
     &                     maxgus,maxstr,maxsze,vmflg,strainFlg)
C
C                A R G U M E N T S
C
	 integer maxsze,maxstr,maxgus,elm
         real*8  x(8), y(8), z(8), c(6,6) 
         double complex v(24)
         double complex stress(maxsze,maxstr,maxgus) 
         double complex strain(maxsze,maxstr,maxgus) 
         logical vmflg,strainFlg
C
C            L O C A L  V A R I A B L E S
C
         real*8  q(8), qx(8), qy(8), qz(8)
         real*8  xinod(8), etanod(8), emunod(8)
         real*8  xi, eta, emu, det
         double complex  epsxx, epsyy, epszz, gamxy, gamyz, gamxz
         integer n
C
C               D   A   T   A
C
         data    xinod/-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0/
         data    etanod/-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0/
         data    emunod/-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0/
C
C               L   O   G   I    C
C
         do 2000  n = 1, 8
           xi  = xinod(n)
           eta = etanod(n)
           emu = emunod(n)
           call  h8shpe(xi,eta,emu,x,y,z,q,qx,qy,qz,det)
C
           if (det .le. 0.0) then
C             write(6,*)'Negative Jacobian determinant in sands17.f'
C             if (det .eq. 0.0) then
C               write(6,*)'Zero jacobian determinant'
C             endif
C             stop  
           endif
C
         epsxx = qx(1)*v(1)  + qx(2)*v(4)  + qx(3)*v(7)  + qx(4)*v(10)
     $         + qx(5)*v(13) + qx(6)*v(16) + qx(7)*v(19) + qx(8)*v(22)

         epsyy = qy(1)*v(2)  + qy(2)*v(5)  + qy(3)*v(8)  + qy(4)*v(11)
     $         + qy(5)*v(14) + qy(6)*v(17) + qy(7)*v(20) + qy(8)*v(23)

         epszz = qz(1)*v(3)  + qz(2)*v(6)  + qz(3)*v(9)  + qz(4)*v(12)
     $         + qz(5)*v(15) + qz(6)*v(18) + qz(7)*v(21) + qz(8)*v(24)
         gamxy = qy(1)*v(1)  + qy(2)*v(4)  + qy(3)*v(7)  + qy(4)*v(10)
     $         + qy(5)*v(13) + qy(6)*v(16) + qy(7)*v(19) + qy(8)*v(22)
     $         + qx(1)*v(2)  + qx(2)*v(5)  + qx(3)*v(8)  + qx(4)*v(11)
     $         + qx(5)*v(14) + qx(6)*v(17) + qx(7)*v(20) + qx(8)*v(23)
         gamyz = qy(1)*v(3)  + qy(2)*v(6)  + qy(3)*v(9)  + qy(4)*v(12)
     $         + qy(5)*v(15) + qy(6)*v(18) + qy(7)*v(21) + qy(8)*v(24)
     $         + qz(1)*v(2)  + qz(2)*v(5)  + qz(3)*v(8)  + qz(4)*v(11)
     $         + qz(5)*v(14) + qz(6)*v(17) + qz(7)*v(20) + qz(8)*v(23)
         gamxz = qz(1)*v(1)  + qz(2)*v(4)  + qz(3)*v(7)  + qz(4)*v(10)
     $         + qz(5)*v(13) + qz(6)*v(16) + qz(7)*v(19) + qz(8)*v(22)
     $         + qx(1)*v(3)  + qx(2)*v(6)  + qx(3)*v(9)  + qx(4)*v(12)
     $         + qx(5)*v(15) + qx(6)*v(18) + qx(7)*v(21) + qx(8)*v(24)
C
C.... ENGINEERING STRESS COMPUTATION
C
         stress(elm,1,n) = c(1,1)*epsxx + c(1,2)*epsyy + c(1,3)*epszz
     $             + c(1,4)*gamxy + c(1,5)*gamyz + c(1,6)*gamxz
         stress(elm,2,n) = c(2,1)*epsxx + c(2,2)*epsyy + c(2,3)*epszz
     $             + c(2,4)*gamxy + c(2,5)*gamyz + c(2,6)*gamxz
         stress(elm,3,n) = c(3,1)*epsxx + c(3,2)*epsyy + c(3,3)*epszz
     $             + c(3,4)*gamxy + c(3,5)*gamyz + c(3,6)*gamxz
         stress(elm,4,n) = c(4,1)*epsxx + c(4,2)*epsyy + c(4,3)*epszz
     $             + c(4,4)*gamxy + c(4,5)*gamyz + c(4,6)*gamxz
         stress(elm,5,n) = c(5,1)*epsxx + c(5,2)*epsyy + c(5,3)*epszz
     $             + c(5,4)*gamxy + c(5,5)*gamyz + c(5,6)*gamxz
         stress(elm,6,n) = c(6,1)*epsxx + c(6,2)*epsyy + c(6,3)*epszz
     $             + c(6,4)*gamxy + c(6,5)*gamyz + c(6,6)*gamxz
C
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
          call vmelmvc(stress,maxgus,maxstr,maxsze,elm,8)
        endif

C
C.... COMPUTE EQUIVALENT STRAIN
C
        if(strainFlg) then
          call strainvm(strain,maxgus,maxstr,maxsze,8)
        endif
C
      return
      end
