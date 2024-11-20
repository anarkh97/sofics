C====================================================================C
      subroutine cfjacobi( ak , am , xx , eigv , nsmax , tol , n)
C====================================================================C
C                                                                    C
C     Performs =   Solves the mass + stiffness generalized eigen     C
C     ----------   system with Jacobi rotations.  Matrices are full  C
C                  and overwriten by the Jacobi reduction.           C
C                  Tolerances are initialized at the beginning with  C
C                  respect to orders of magnitude of mass and        C
C                  stiffness matrices entries.                       C
C                                                                    C
C                                                                    C
C     Inputs/Outputs =                                               C
C     ----------------                                               C
C     AK     <input/output>  full stiffness matrix into subspace     C
C     AM     <input/output>  full mass matrix into subspace          C
C     XX     <output>        eigenvectors                            C
C     EIGV   <output>        eigenvalues                             C
C     NSMAX  <input>         maximum number of Jacobi sweeps         C
C     TOL    <input>         tolerance for convergence               C
C     N      <input>         dimension of the system                 C
C     OUT    <input>         logical unit number for error-messages  C
C                                                                    C
C                                                                    C
C     Warning =   Works only for positive definite systems.          C
C     ---------                                                      C
C                                                                    C
C                                                                    C
C     Outputs =   none, except warning and error-messages.           C
C     ---------                                                      C
C                                                                    C
C====================================================================C
C=Author =Charbel FARHAT-Modified by Francois HEMEZ                  C
C=Date   =October,15,1991                                            C
C=Version=2.0                                                        C
C=Comment=                                                           C
C====================================================================C
C
C.....DECLARES GLOBAL VARIABLES
C
      integer        out , n , nsmax
      real*8         tol , ak(*) , am(*) , xx(n,n) , eigv(n)
C
C.....DECLARES LOCAL VARIABLES
C
      real*8         zero , two , eps , epsam , epsak , d1 , d2
      real*8         den , akk , ajj , ab , xxj , xxk , tr
      real*8         check , sqch , ca , cg , akj , amj , amk
      real*8         ze1eps , ze2eps , xmax , xmin
      integer        nsweep , nr , iik , kk , jk , ii , jj
      integer        jp1 , jm1 , kp1 , km1 , ij , ik , ji , ki
      integer        i , j , k , nt , iconv
*     character*20   routin
*     character*60   messag
C
C--------------------------------------------------------------------C
C     L O G I C                                                      C
C--------------------------------------------------------------------C
C
*     routin = 'JACOBI'
      zero   = 0.0D0
      two    = 2.0D0
      out    = 6
C
C.....INITIALIZES TOLERANCES 'ZE1EPS' AND 'ZE2EPS'
C
      if ( am(1).eq.zero ) go to 100
      xmax = max(abs(ak(1)),abs(am(1)))
      xmin = min(abs(ak(1)),abs(am(1)))
      do 1000 i=2,n
         ii   = i*(i+1)/2
         if ( am(ii).eq.zero ) go to 100
         xmax = max(xmax,max(abs(ak(ii)),abs(am(ii))))
         xmin = min(xmin,min(abs(ak(ii)),abs(am(ii))))
 1000 continue
      ze1eps = 1.0D-35*0.5D0*(xmax+xmin)
C
      xmin = abs(am(1))
      xmax = abs(am(1))
      do 1100 i=2,(n*(n+1)/2)
         xmin = min(xmin,abs(am(i)))
         xmax = max(xmax,abs(am(i)))
 1100 continue
      ze2eps = 1.0D-08*0.5D0*(xmax+xmin)
C
C.....INITIALIZES EIGEN VECTORS 'XX'
C
      do 2000 i=1,n
         ii = i*(i-1)/2+i
         if ( am(ii).eq.zero ) go to 100
         if (abs(ak(ii)).lt.ze1eps.or.abs(am(ii)).lt.ze1eps) go to 110
         eigv(i) = ak(ii)/am(ii)
         do 2010 j=1,n
            xx(i,j) = zero
 2010    continue
         xx(i,i) = 1.0
 2000 continue
C
C.....SETS THE COUNTER
C
      nsweep = 0
      nr     = n-1
C
C-------------------------C
C     JACOBI ITERATIONS   C
C-------------------------C
C
   10 continue
      nsweep = nsweep+1
C
C.....CHECKS IF ZEROING IS REQUIRED
C
      eps = 0.01D0**nsweep
      eps = eps*eps
C
      do 3000 j=1,nr
         iik = j+1
         do 3100 k=iik,n
            jj    = j*(j-1)/2+j
            kk    = k*(k-1)/2+k
            jk    = k*(k-1)/2+j
            epsak = (ak(jk)*ak(jk))/(ak(jj)*ak(kk))
            epsam = (am(jk)*am(jk))/(am(jj)*am(kk))
            if ( epsak.ge.eps.or.epsam.ge.eps ) then
C
C.....CALCULATES ROTATION ELEMENTS
C
               akk   = ak(kk)*am(jk)-am(kk)*ak(jk)
               ajj   = ak(jj)*am(jk)-am(jj)*ak(jk)
               ab    = (ak(jj)*am(kk)-ak(kk)*am(jj))/two
               check = ab*ab+akk*ajj
               if ( check.lt.zero ) then
                  write(6,*)'Check & ze2eps',check,ze2eps
                  if ( abs(check).lt. ze2eps ) then
                     check = abs(check)
                  else
                     go to 200
                  endif
               endif
               sqch = dsqrt(check)
               d1   = ab+sqch
               d2   = ab-sqch
               den  = d1
               if ( abs(d2).gt.abs(d1) ) den = d2
               if ( den.eq.zero ) then
                  ca = zero
                  cg = -ak(jk)/ak(kk)
               else
                  ca = akk/den
                  cg = -ajj/den
               endif
C
C.....PERFORMS THE GENERALIZED ROTATION
C
               jp1 = j+1
               jm1 = j-1
               kp1 = k+1
               km1 = k-1
               if ( jm1.ge.1 ) then
                  do 3101 i=1,jm1
                     ij     = j*jm1/2+i
                     ik     = k*km1/2+i
                     akj    = ak(ij)
                     akk    = ak(ik)
                     amj    = am(ij)
                     amk    = am(ik)
                     ak(ij) = akj+cg*akk
                     am(ij) = amj+cg*amk
                     ak(ik) = akk+ca*akj
                     am(ik) = amk+ca*amj
 3101             continue
               endif
               if ( kp1.le.n ) then
                  do 3102 i=kp1,n
                     ji     = i*(i-1)/2+j
                     ki     = i*(i-1)/2+k
                     akj    = ak(ji)
                     amj    = am(ji)
                     akk    = ak(ki)
                     amk    = am(ki)
                     ak(ji) = akj+cg*akk
                     am(ji) = amj+cg*amk
                     ak(ki) = akk+ca*akj
                     am(ki) = amk+ca*amj
 3102             continue
               endif
               if ( jp1.le.km1 ) then
                  do 3103 i=jp1,km1
                     ji     = i*(i-1)/2+j
                     ik     = k*(k-1)/2+i
                     akj    = ak(ji)
                     amj    = am(ji)
                     akk    = ak(ik)
                     amk    = am(ik)
                     ak(ji) = akj+cg*akk
                     am(ji) = amj+cg*amk
                     ak(ik) = akk+ca*akj
                     am(ik) = amk+ca*amj
 3103             continue
               endif
               akk    = ak(kk)
               amk    = am(kk)
               akj    = ak(jj)
               amj    = am(jj)
               ak(kk) = akk+two*ca*ak(jk)+ca*ca*akj
               am(kk) = amk+two*ca*am(jk)+ca*ca*amj
               ak(jj) = akj+two*cg*ak(jk)+cg*cg*akk
               am(jj) = amj+two*cg*am(jk)+cg*cg*amk
               ak(jk) = zero
               am(jk) = zero
C
C.....UPDATES EIGEN VECTOR FOR THIS ROTATION
C
               do 3104 i=1,n
                  xxj     = xx(i,j)
                  xxk     = xx(i,k)
                  xx(i,j) = xxj+cg*xxk
                  xx(i,k) = xxk+ca*xxj
 3104          continue
            endif
 3100    continue
 3000 continue
C
C------------------------------C
C     UPDATES EIGEN VALUES     C
C     CHECKS FOR CONVERGENCE   C
C------------------------------C
C
      nt    = n*(n+1)/2
      iconv = 0
      do 4000 i=1,n
         ii = i*(i-1)/2+i
         if ( am(ii).eq.zero ) go to 100
         if (abs(ak(ii)).lt.ze1eps.or.abs(am(ii)).lt.ze1eps) go to 300
         tr      = ak(ii)/am(ii)
         den     = (tr-eigv(i))/tr
         eigv(i) = tr
         if ( abs(den).gt.tol ) iconv = 1
 4000 continue
C
      if ( iconv.eq.1 ) go to 20
C
C.....CHECKS OFF-DIAGONAL ENTRIES
C
      eps = tol*tol
C
      do 5000 j=1,nr
         iik = j+1
         do 5100 k=iik,n
            jj    = j*(j-1)/2+j
            kk    = k*(k-1)/2+k
            jk    = k*(k-1)/2+j
            epsak = (ak(jk)*ak(jk))/(ak(jj)*ak(kk))
            epsam = (am(jk)*am(jk))/(am(jj)*am(kk))
            if ( epsak.ge.eps.or.epsam.ge.eps ) go to 20
 5100    continue
 5000 continue
C
C.....CONVERGENCE IS ASSUMED - GO TO POST-TREATMENT
C
      go to 30
C
C.....KEEPS ON ITERATING AS LONG AS THE MAXIMUM
C.....NUMBER OF SWEEPS HAS NOT BEEN REACHED YET
C
   20 continue
      if ( nsweep.le.nsmax ) then
         go to 10
      else
         write(out,9) nsmax
         write(*  ,9) nsmax
      endif
C
C-------------------------------------------------------C
C     SCALES THE EIGEN VECTORS AND ENDS THE TREATMENT   C
C-------------------------------------------------------C
C
   30 continue
      do 6000 i=1,n
         ii = i*(i-1)/2+i
         if ( am(ii).lt.zero ) then
            if ( abs(am(ii)).lt.ze2eps ) then
               akk = dsqrt(-am(ii))
            else
               go to 400
            endif
         else
            akk = dsqrt(am(ii))
         endif
         akk = 1.0D0/akk
         do 6100 j=1,n
            xx(j,i) = akk*xx(j,i)
 6100    continue
 6000 continue
C
C.....RETURN
C
      return
C
C.....FORMATS
C
    1 format('*** ERROR in JACOBI routine')
   11 format('*** Zero diagonal mass entry')
   12 format('*** Equation ',I5,' AM( ',I5,' ) = ',E16.9)
   21 format('*** Matrix was found not definite (initialization)')
   22 format('*** Precision for this test is ZE1EPS = ',E16.9)
   31 format('*** Equation ',I5,' Kii = ',E16.9,' Mii = ',E16.9)
   32 format('*** Ratio stiffness/mass = ',E16.9)
    4 format('*** Problem occurs when computing a rotation')
    5 format('*** Can not take the square root of CHECK = ',E16.9)
   60 format('*** Matrix was found not definite')
   70 format('*** Can not take the square root')
   71 format('*** Precision for this test is ZE2EPS = ',E16.9)
    8 format('*** Equation ',I5,' AM( ',I5,' ) = ',E16.9)
    9 format('*** WARNING in JACOBI : TOLJAC could not be reached',/,
     $       '***                     even after ',I5,' sweeps')
C
C.....ERROR-MESSAGE FOR A ZERO MASS DIAGONAL ENTRY
C
  100 continue
*      write(out, 1)
*      write(out,11)
*      write(out,12) i,ii,am(ii)
      write(*  , 1) 
      write(*  ,11)
      write(*  ,12) i,ii,am(ii)
      stop
C
C.....ERROR-MESSAGE IF MATRIX IS NOT DEFINITE (INITIALIZATION)
C
 110  continue
      write(out, 1)
      write(out,21)
      write(out,22) ze1eps
      write(out,31) i,ak(ii),am(ii)
      write(out,32) ak(ii)/am(ii)
      write(*  , 1)
      write(*  ,21)
      write(*  ,22) ze1eps
      write(*  ,31) i,ak(ii),am(ii)
      write(*  ,32) ak(ii)/am(ii)
      stop
C
C.....ERROR-MESSAGE IF 'CHECK' IS NEGATIVE
C
  200 continue
      write(out,1)
      write(out,4)
      write(out,5) check
      write(*  ,1)
      write(*  ,4)
      write(*  ,5) check
      stop
C
C.....ERROR-MESSAGE IF MATRIX IS NOT DEFINITE (UPDATING)
C
  300 continue
      write(out, 1)
      write(out,60)
      write(out,22) ze1eps
      write(out,31) i,ak(ii),am(ii)
      write(out,32) ak(ii)/am(ii)
      write(*  , 1)
      write(*  ,60)
      write(*  ,22) ze1eps
      write(*  ,31) i,ak(ii),am(ii)
      write(*  ,32) ak(ii)/am(ii)
      stop
C
C.....ERROR-MESSAGE IF SCALING CAN NOT BE PERFORMED
C
  400 continue
      write(out, 1)
      write(out,70)
      write(out,71) ze2eps
      write(out, 8) i,ii,am(ii)
      write(*  , 1)
      write(*  ,70)
      write(*  ,71) ze2eps
      write(*  , 8) i,ii,am(ii)
      stop
C
      end
C=end of routine 'JACOBI'
C========================
