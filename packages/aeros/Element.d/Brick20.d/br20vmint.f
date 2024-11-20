C=PURPOSE Compute integral of (Sigma_vm / Sigma_bar) of 8-node hexa
C
C     input arguments are
C
C       X         (20 x 1) array of x coordinates of hexahedron nodes
C       Y         (20 x 1) array of y coordinates of hexahedron nodes
C       Z         (20 x 1) array of z coordinates of hexahedron nodes
C       C         (6 x 6) constitutive material matrix 
C       P         Gauss quadrature rule (no. of points)
C       v         (24x1) array of element node displacement arranged in
C                 ux1,uy1,uz1,ux2,uy2,uz2, ............., uz8
C       sigbar    given stress value
C
C     The outputs are:
C
C       VMINT     integral of (Von Mises stresse / stress_bar) 
C       VOL       volume of element
C       STATUS    Status integer variable.  Zero if no error
C                 detected.
C
      subroutine  br20vmint(x,y,z,v,c,p,vmint,vol,sigbar,
     &                      fac,iarea,status)
C
C                   A R G U M E N T S
C
C      implicit none
      integer p,status
      real*8  x(20), y(20), z(20), c(6,6), v(*)
      real*8  vmint, sigbar, fac, vol
C
C                   L O C A L   V A R I A B L E S
C
      real*8  q(20), qx(20), qy(20),qz(20)
      real*8  xi, eta, emu, det, w, weight, a, b, pow, comp
      real*8  epsxx, epsyy,epszz, gamxy, gamyz, gamxz
      real*8  stress(6),dsxx, dsyy, dszz, dsxy, dsxz, dsyz, j2
      integer k, l, jj, i
C
C
C                   L O G I C
C
      pow(a,b) = exp(b*log(a))
C
C     Initialize 
      status = 0
      vmint  = 0.0d0
      vol    = 0.0d0
C
      do 3000  k = 1,p
        do 2500  l = 1,p
          do 2400  jj = 1,p
            call hxgaus20(p,k,p,l,p,jj,xi,eta,emu,weight)
            call h20shpe(xi,eta,emu,x,y,z,q,qx,qy,qz,det)
C
            if (det .le. 0.0d0) then
              write(6,*) 'Negative Jacobian determinant in brintvm.f'
              status = -1
              return
            end if
C
            w = weight * det
C
            epsxx = 0.0d0
            epsyy = 0.0d0
            epszz = 0.0d0
            gamxy = 0.0d0
            gamyz = 0.0d0
            gamxz = 0.0d0

            do 1500 i = 1,20

              epsxx =epsxx +qx(i)*v(3*i-2)
              epsyy =epsyy +qy(i)*v(3*i-1) 
              epszz =epszz +qz(i)*v(3*i  )  
              gamxy =gamxy +qy(i)*v(3*i-2) +qx(i)*v(3*i-1)
              gamyz =gamyz +qz(i)*v(3*i-1) +qy(i)*v(3*i  )
              gamxz =gamxz +qx(i)*v(3*i  ) +qz(i)*v(3*i-2)
 1500       continue
C
C.... ENGINEERING STRESS COMPUTATION
C
            stress(1) = c(1,1)*epsxx + c(1,2)*epsyy + c(1,3)*epszz
     $                + c(1,4)*gamxy + c(1,5)*gamyz + c(1,6)*gamxz
            stress(2) = c(2,1)*epsxx + c(2,2)*epsyy + c(2,3)*epszz
     $                + c(2,4)*gamxy + c(2,5)*gamyz + c(2,6)*gamxz
            stress(3) = c(3,1)*epsxx + c(3,2)*epsyy + c(3,3)*epszz
     $                + c(3,4)*gamxy + c(3,5)*gamyz + c(3,6)*gamxz
            stress(4) = c(4,1)*epsxx + c(4,2)*epsyy + c(4,3)*epszz
     $                + c(4,4)*gamxy + c(4,5)*gamyz + c(4,6)*gamxz
            stress(5) = c(5,1)*epsxx + c(5,2)*epsyy + c(5,3)*epszz
     $                + c(5,4)*gamxy + c(5,5)*gamyz + c(5,6)*gamxz
            stress(6) = c(6,1)*epsxx + c(6,2)*epsyy + c(6,3)*epszz
     $                + c(6,4)*gamxy + c(6,5)*gamyz + c(6,6)*gamxz
C
C.... COMPUTE THE FIRST DEVEATORIC STRESSES
C
            comp = (stress(1) + stress(2) + stress(3))/3.0d0
            dsxx = stress(1) - comp
            dsyy = stress(2) - comp
            dszz = stress(3) - comp
            dsxy = stress(4)
            dsyz = stress(5)
            dsxz = stress(6)
C
C.... COMPUTE THE SECOND DEVEATORIC STRESS
C
             j2 = ((dsxx*dsxx)+(dsyy*dsyy)+(dszz*dszz))/2.0d0+
     &             (dsxy*dsxy)+(dsyz*dsyz)+(dsxz*dsxz)
             j2 = dsqrt(3.0d0*j2)

C
C.... COMPUTE THE VON MISES STRESS AND VOLUME
C
             if (iarea.eq.0) w = 1.0d0

             vmint = vmint + w*pow((j2/sigbar),fac)
             vol   = vol + w
C
 2400       continue
 2500     continue
 3000   continue

C
C.... COMPUTE MAXIMUM STRESS VALUE AT GAUSS POINT; RETURN MAXVAL^FAC
C
       
      if (iarea.eq.0) then 
        vmint  = vmint/vol
        vol    = 1.0d0
      endif

      return
      end
