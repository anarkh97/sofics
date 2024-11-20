C=====================================================================C
      subroutine mass7(elm,emass,A,rho,x,y,z,gamma,grvfor
     .,                                     grvflg,totmas,masflg)
C=====================================================================C
C                                                                     C
C     Assemble the Elemental Mass Matrix of a 3D Timoschenko          C
C     Beam Finite Element. Lumping is Assumed Here.                   C
C                                                                     C
C     The Output Mass Matrix is a Block 12 by 12 Stored in the        C
C     Upper Left Corner of [emass] (First 12 Rows and Columns).       C
C                                                                     C
C     Francois M. Hemez - July 11th 1994 - Version 1.0                C
C                                                                     C
C     Original Parallel C Code By J.C. Chiou.                         C
C     See Directory: /mars/limbo/chiou/Beam.d/Timo.d.                 C
C     Gravity Body Force Modification November 1994 P.R. Stern        C
C                                                                     C
C=====================================================================C
C
C     ------------
C     DECLARATIONS
C     ------------
C
C.....GLOBAL VARIABLES
C
      integer   elm
      real*8  A,rho,emass(12,12),x(*),y(*),z(*)
      real*8  gamma(*),grvfor(*),totmas
      logical grvflg,masflg
C
C.....LOCAL VARIABLES
C
      integer   i , j
      real*8    zero , massdis , massrot, c1
      real*8    dx , dy , dz , length , two , twentyfour
C
C     ----
C     DATA
C     ----
C
      data zero       /0.000000D+00/
      data two        /0.200000D+01/
      data twentyfour /0.240000D+02/
C
C     -----
C     LOGIC
C     -----
C
C.....CLEAR THE OUTPUT MASS MATRIX
C
      do 1001 j=1,12
         do 1002 i=1,12
            emass(i,j) = zero
 1002    continue
 1001 continue
C
C.....COMPUTE THE LENGTH OF THE BEAM ELEMENT
C
      dx     = x(2) - x(1)
      dy     = y(2) - y(1)
      dz     = z(2) - z(1)
      length = dsqrt( dx*dx + dy*dy + dz*dz )
C
      if ( length.eq.zero ) go to 100
C
C.....INITIALIZE THE MASS COEFFICIENT FOR TRANSLATIONAL DOFs
C
      massdis = rho*A*length/two
C
C.....INITIALIZE THE MASS COEFFICIENT FOR ROTATIONAL DOFs
C
      massrot = rho*A*length*length*length/twentyfour
C
C.....ASSEMBLE THE (LUMPED) OUTPUT ELEMENTAL MASS MATRIX
C
      emass( 1, 1) = massdis
      emass( 2, 2) = massdis
      emass( 3, 3) = massdis
C
      emass( 4, 4) = massrot
      emass( 5, 5) = massrot
      emass( 6, 6) = massrot
C
      emass( 7, 7) = massdis
      emass( 8, 8) = massdis
      emass( 9, 9) = massdis
C
      emass(10,10) = massrot
      emass(11,11) = massrot
      emass(12,12) = massrot
C
C
C..... COMPUTE THE BODY FORCE DUE TO GRAVITY IF NEEDED
C
      if (grvflg) then
        c1 = 2.0d0*massdis
        grvfor(1) = c1*gamma(1)
        grvfor(2) = c1*gamma(2)
        grvfor(3) = c1*gamma(3)
      endif
C
C.... ACCUMULATE THE SUBDOMAIN MASS
C
        if (masflg) then
          totmas = totmas + 2.0d0*massdis
        endif
C     ------
C     RETURN
C     ------
C
      return
C
C     ---------------
C     ERROR-TREATMENT
C     ---------------
C
C.....ERROR-MESSAGE IF THE BEAM ELEMENT HAS ZERO LENGTH
C
  100 continue
      write(*,*) "*** FATAL ERROR in Routine MASS7  ***"
      write(*,*) "*** The Timoschenko Beam Element  ***"
      write(*,*) "*** Has Zero Length!              ***"
      write(*,*) "*** ... All Treatments Terminated ***"
      stop
C
C     ---
C     END
C     ---
C
      end
C=end of routine "MASS7"
C======================C
