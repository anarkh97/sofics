C     The calling sequence is
C
C  CALL SLSAS2 (ESCM, X, Y, K, T, TMPGRD, MAXGUS, MAXSTR, ELM, NUMEL)
C
C     where the input arguments are
C
C       ESCM      Character string defining flux computation method:
C                 DIRECT  direct flux evaluation at corners
C                 EXTRAP  extrapolate from 2 x 2 Gauss points
C       X         (4 x 1) array of x coordinates of quadrilateral nodes
C       Y         (4 x 1) array of y coordinates of quadrilateral nodes
C       K         Should be 1.0 (Use to be Thermal Conduction Coefficient)
C       T         (4 x 1) array of element node temperatures          
C       ELM       Element Number
C
C
C     The outputs are:
C
C TMPGRD   (NUMEL x 2 x 4) array of corner node temperature gradients:
C                          potgx1,potgy1,potgx2,potgy2....,potgy4 
C
C=END USAGE
C=BLOCK FORTRAN

C     subroutine  slsas2(escm, x, y, k, t, potgrd, flux, maxgus, maxstr, 
C    &                              elm, numel)
      subroutine  slsas2(escm, x, y, k, t, potgrd, maxgus, maxstr, 
     &                              elm, numel)
C
C                   A R G U M E N T S
C
      integer           elm,numel,maxgus,maxstr
      character*(*)     escm
      real*8  x(4), y(4), k, t(4)
C     real*8  potgrd(numel,maxstr,maxgus), flux(numel,maxstr,maxgus) 
      real*8  potgrd(numel,maxstr,maxgus)
C
C                   L O C A L   V A R I A B L E S
C
      real*8  q(4), qx(4), qy(4)
C     real*8  xinod(4), etanod(4), sigauss(3), cext(4,4)
      real*8  xinod(4), etanod(4)
      real*8  xi, eta, det, potgx, potgy
      integer           i, n
C
C                   D A T A
C
      data          xinod  /-1.0, 1.0, 1.0,-1.0/
      data          etanod /-1.0,-1.0, 1.0, 1.0/
C     data          cext / 1.866025404, -0.5, 0.133974596, -0.5,
C    $                     -0.5, 1.866025404, -0.5, 0.133974596,
C    $                     0.133974596, -0.5, 1.866025404, -0.5,
C    $                     -0.5, 0.133974596, -0.5, 1.866025404 /
C
C                   L O G I C
C
C
      if (escm(1:1) .eq. 'D')                 then
        do 2000  n = 1,4
          xi =     xinod (n)
          eta =    etanod(n)
          call     Q4SHPE (xi, eta, x, y, q, qx, qy, det)
          if (det .le. 0.0)        then
            write(6,*) 'Negative Jacobian determinant'
            if (det .eq. 0.0)      then
              write(6,*)  'Zero Jacobian determinant'
            end if
            stop 
          end if
          potgx = qx(1)*t(1) + qx(2)*t(2) + qx(3)*t(3) + qx(4)*t(4)
          potgy = qy(1)*t(1) + qy(2)*t(2) + qy(3)*t(3) + qy(4)*t(4)

C
C         flux(elm,1,n)   = -k*potgx
C         flux(elm,2,n)   = -k*potgy
          potgrd(elm,1,n) = potgx 
          potgrd(elm,2,n) = potgy 
C
 2000     continue
C
      else

        do 2200  n = 1,4
C         flux(elm,1,n) = 0.0
C         flux(elm,2,n) = 0.0
 2200     continue

        do 3000  i = 1,4
          xi =     xinod (i)*0.577350269
          eta =    etanod(i)*0.577350269
          call     Q4SHPE (xi, eta, x, y, q, qx, qy, det)
          if (det .le. 0.0)        then
            write(6,*) 'Negative Jacobian determinant'
            if (det .eq. 0.0)      then
              write(6,*) 'Zero Jacobian determinant'
            end if
            stop 
          end if

          potgx = qx(1)*t(1) + qx(2)*t(2) + qx(3)*t(3) + qx(4)*t(4)
          potgy = qy(1)*t(1) + qy(2)*t(2) + qy(3)*t(3) + qy(4)*t(4)
C         sigauss(1) = -k*potgx 
C         sigauss(2) = -k*potgy 

          do 2500  n = 1,4
C           flux(elm,1,n)   =  potgrd(elm,1,n) + cext(i,n)*sigauss(1)
C           flux(elm,2,n)   =  potgrd(elm,2,n) + cext(i,n)*sigauss(2)
            potgrd(elm,1,n) =  potgx 
            potgrd(elm,2,n) =  potgy 
 2500       continue
 3000     continue
      end if
C
C
      return
      end
