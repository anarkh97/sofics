C     The calling sequence is
C
C  CALL HTSAS2 (ESCM, X, Y, K, T, TMPGRD, FLUX, MAXGUS, MAXSTR, ELM, NUMEL)
C
C     where the input arguments are
C
C       ESCM      Character string defining flux computation method:
C                 DIRECT  direct flux evaluation at corners
C                 EXTRAP  extrapolate from 2 x 2 Gauss points
C       X         (4 x 1) array of x coordinates of quadrilateral nodes
C       Y         (4 x 1) array of y coordinates of quadrilateral nodes
C       K         Thermal Conduction Coefficient 
C       T         (4 x 1) array of element node temperatures          
C       ELM       Element Number
C
C
C     The outputs are:
C
C TMPGRD   (NUMEL x 2 x 4) array of corner node temperature gradients:
C                          tmpgx1,tmpgy1,tmpgx2,tmpgy2....,tmpgy4 
C FLUX   (NUMEL x 2 x 4) array of heat fluxes:
C                          fluxx1,fluxy1,fluxx2,fluxy2....,fluxy4 
C
C=END USAGE
C=BLOCK FORTRAN

      subroutine  htsas2(escm, x, y, k, t, tmpgrd, flux, maxgus, maxstr, 
     &                              elm, numel)
C
C                   A R G U M E N T S
C
      integer           elm,numel,maxgus,maxstr
      character*(*)     escm
      real*8  x(4), y(4), k, t(4)
      real*8  tmpgrd(numel,maxstr,maxgus), flux(numel,maxstr,maxgus) 
C
C                   L O C A L   V A R I A B L E S
C
      real*8  q(4), qx(4), qy(4)
      real*8  xinod(4), etanod(4), sigauss(3), cext(4,4)
      real*8  xi, eta, det, tmpgx, tmpgy
      integer           i, n
C
C                   D A T A
C
      data          xinod  /-1.0, 1.0, 1.0,-1.0/
      data          etanod /-1.0,-1.0, 1.0, 1.0/
      data          cext / 1.866025404, -0.5, 0.133974596, -0.5,
     $                     -0.5, 1.866025404, -0.5, 0.133974596,
     $                     0.133974596, -0.5, 1.866025404, -0.5,
     $                     -0.5, 0.133974596, -0.5, 1.866025404 /
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
          tmpgx = qx(1)*t(1) + qx(2)*t(2) + qx(3)*t(3) + qx(4)*t(4)
          tmpgy = qy(1)*t(1) + qy(2)*t(2) + qy(3)*t(3) + qy(4)*t(4)

C
          flux(elm,1,n)   = -k*tmpgx
          flux(elm,2,n)   = -k*tmpgy
          tmpgrd(elm,1,n) = tmpgx 
          tmpgrd(elm,2,n) = tmpgy 
C
 2000     continue
C
      else

        do 2200  n = 1,4
          flux(elm,1,n) = 0.0
          flux(elm,2,n) = 0.0
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

          tmpgx = qx(1)*t(1) + qx(2)*t(2) + qx(3)*t(3) + qx(4)*t(4)
          tmpgy = qy(1)*t(1) + qy(2)*t(2) + qy(3)*t(3) + qy(4)*t(4)
          sigauss(1) = -k*tmpgx 
          sigauss(2) = -k*tmpgy 

          do 2500  n = 1,4
            flux(elm,1,n)   =  tmpgrd(elm,1,n) + cext(i,n)*sigauss(1)
            flux(elm,2,n)   =  tmpgrd(elm,2,n) + cext(i,n)*sigauss(2)
            tmpgrd(elm,1,n) =  tmpgx 
            tmpgrd(elm,2,n) =  tmpgy 
 2500       continue
 3000     continue
      end if
C
C
      return
      end
