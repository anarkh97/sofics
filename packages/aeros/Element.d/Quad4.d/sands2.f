C=AUTHOR C. A. Felippa, May 1967
C=REVISED P. R. STERN, MARCH 1990
C=REVISED K. H. PIERSON, MARCH 1997 
C     Given the node displacements, SANDS2 computes corner
C     point stresses on a four-node quadrilateral in plane stress
C
C
C     The calling sequence is
C
C       CALL   SANDS2 (ESCM, X, Y, C, V, STRESS, STRAIN, MAXGUS, MAXSTR, ELM, NUMEL, VMF, TC, TREF, NDTEMPS)
C
C     where the input arguments are
C
C       ESCM      Character string defining stress computation method:
C                 DIRECT  direct stress evaluation at corners
C                 EXTRAP  extrapolate from 2 x 2 Gauss points
C       X         (4 x 1) array of x coordinates of quadrilateral nodes
C       Y         (4 x 1) array of y coordinates of quadrilateral nodes
C       C         (3 x 3) stress-strain constitutive matrix
C       V         (8 x 1) array of element node displacements arranged
C                    ux1,uy1, ux2,uy2 ...  uy4
C       ELM       Element Number
C       VMFLAG    VonMises Stress Flag
C       STRAINFLG VonMises Strain Flag
C       TC        TC=E*alpha/(1-nu)
C       TREF      Reference Temperature
C       NDTEMPS   The temperatures of the quad
C
C
C     The outputs are:
C
C       STRESS    (NUMEL x 3 x 4) array of corner node stresses arranged
C                  sigxx1,sigyy1,tauxy1, sigxx2,sigyy2,tauxy2, ... tauxy4
C       STRAIN    (NUMEL x 3 x 4) array of corner node strains arranged just like
C                  the stresses
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine  sands2(escm, x, y, c, v, stress, strain, maxgus,
     &                   maxstr, elm, msize, vmflg, strainFlg, 
     &                   tc, tref, ndtemps)
C
C                   A R G U M E N T S
C
      integer           elm,msize,maxgus,maxstr
      character*(*)     escm
      real*8  x(4), y(4), c(3,3), v(8)
      real*8  stress(msize,maxstr,maxgus),strain(msize,maxstr,maxgus) 
      real*8  ndtemps(4), tref, tc
      logical vmflg,strainFlg
C
C                   L O C A L   V A R I A B L E S
C
      real*8  q(4), qx(4), qy(4)
      real*8  xinod(4), etanod(4), sigauss(3), cext(4,4)
      real*8  xi, eta, det, epsxx, epsyy, gamxy
      real*8  tl(4), tgp, eptxo, eptyo
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

C
C.... COMPUTE THE THERMAL FIELD
C
       do i = 1,4
          tl(i) = ndtemps(i) - tref
        enddo
    

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

C
C.... COMPUTE THE THERMAL STRESS
C
          tgp   =  q(1)*tl(1)+q(2)*tl(2)+q(3)*tl(3)+q(4)*tl(4)
          eptxo = tc*tgp
          eptyo = tc*tgp

C
C.... COMPUTE THE TOTAL STRAIN FROM THE COUPLED SOLVE FOR DISPLACEMENTS
C
          epsxx = qx(1)*v(1) + qx(2)*v(3) + qx(3)*v(5) + qx(4)*v(7)
          epsyy = qy(1)*v(2) + qy(2)*v(4) + qy(3)*v(6) + qy(4)*v(8)
          gamxy = qy(1)*v(1) + qy(2)*v(3) + qy(3)*v(5) + qy(4)*v(7)
     $          + qx(1)*v(2) + qx(2)*v(4) + qx(3)*v(6) + qx(4)*v(8)

C
C.... COMPUTE THE TOTAL STRESS
C
          stress(elm,1,n) = c(1,1)*epsxx + c(1,2)*epsyy + c(1,3)*gamxy
     &                      - eptxo
          strain(elm,1,n) = epsxx 
          stress(elm,2,n) = c(2,1)*epsxx + c(2,2)*epsyy + c(2,3)*gamxy
     &                      - eptyo
          strain(elm,2,n) = epsyy 
          stress(elm,4,n) = c(3,1)*epsxx + c(3,2)*epsyy + c(3,3)*gamxy
          strain(elm,4,n) = gamxy 
 2000     continue
C
C.... EXTRAPOLATE FROM THE GAUSS POINTS
c
      else
        do 2200  n = 1,4
          stress(elm,1,n) = 0.0
          stress(elm,2,n) = 0.0
          stress(elm,3,n) = 0.0
          stress(elm,4,n) = 0.0
          strain(elm,1,n) = 0.0    
          strain(elm,2,n) = 0.0    
          strain(elm,3,n) = 0.0    
          strain(elm,4,n) = 0.0    
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

C
C.... COMPUTE THE THERMAL STRESS
C
          tgp   =  q(1)*tl(1)+q(2)*tl(2)+q(3)*tl(3)+q(4)*tl(4)
          eptxo = tc*tgp
          eptyo = tc*tgp

C
C.... COMPUTE THE TOTAL STRAIN FROM THE COUPLED SOLVE FOR DISPLACEMENTS
C
          epsxx = qx(1)*v(1) + qx(2)*v(3) + qx(3)*v(5) + qx(4)*v(7)
          epsyy = qy(1)*v(2) + qy(2)*v(4) + qy(3)*v(6) + qy(4)*v(8)
          gamxy = qy(1)*v(1) + qy(2)*v(3) + qy(3)*v(5) + qy(4)*v(7)
     $          + qx(1)*v(2) + qx(2)*v(4) + qx(3)*v(6) + qx(4)*v(8)

C
C.... COMPUTE THE TOTAL STRESS
C
          sigauss(1) = c(1,1)*epsxx + c(1,2)*epsyy + c(1,3)*gamxy
     &                 - eptxo
          sigauss(2) = c(2,1)*epsxx + c(2,2)*epsyy + c(2,3)*gamxy
     &                 - eptyo
          sigauss(3) = c(3,1)*epsxx + c(3,2)*epsyy + c(3,3)*gamxy

          do 2500  n = 1,4
            stress(elm,1,n) =  stress(elm,1,n) + cext(i,n)*sigauss(1)
            strain(elm,1,n) =  epsxx 
            stress(elm,2,n) =  stress(elm,2,n) + cext(i,n)*sigauss(2)
            strain(elm,2,n) =  epsyy 
            stress(elm,4,n) =  stress(elm,4,n) + cext(i,n)*sigauss(3)
            strain(elm,4,n) =  gamxy 
C            print*,stress(elm,5,n), stress(elm,6,n)
 2500       continue
 3000     continue
      end if
C
C     set remain stress/strain components to zero
C
      do  n = 1,4
        do i = 5,maxstr
          stress(elm,i,n) = 0.0
          strain(elm,i,n) = 0.0    
        enddo
      enddo
C
C.... COMPUTE THE VON MISES STRESS IN THE PLATE
C
      if (vmflg) then
        call vmelmv(stress,maxgus,maxstr,msize,elm,4)
      endif
C
C.... COMPUTE THE VON MISES STRAIN IN THE PLATE
C
      if (strainFlg) then
        call strainvm(strain,maxgus,maxstr,msize,4)
      endif

      return
      end
