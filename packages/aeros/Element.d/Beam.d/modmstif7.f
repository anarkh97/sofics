C=====================================================================C
C        1         2         3         4         5         6         7C
C234567890123456789012345678901234567890123456789012345678901234567890C
C=====================================================================C
      subroutine modmstif7(       estif , A      , E  ,
     $                   eframe , Ix    , Iy     , Iz , alphay ,
     $                   alphaz , C1    , nu     , x  , y , z,
     $                   flag  )
C=====================================================================C
C                                                                     C
C     Assemble the Elemental Stiffness Matrix of a 3D Timoschenko     C
C     Beam Finite Element that Couples Axial Deformation, Bending     C
C     Deformation and Transverse Shear Deformation.                   C
C     Cross-sectional Areas Normal to the Neutral Axis in the         C
C     Undeformed Configuration May Experience Small Angles With       C
C     Respect To the Neutral Axis in the Deformed Configuration.      C
C     However, They Remain Planar.                                    C
C                                                                     C
C     The Output Stiffness Matrix is a Block 12 by 12 Stored in the   C
C     Upper Left Corner of [estif] (First 12 Rows and Columns).       C
C                                                                     C
C     Original Parallel C Code By J.C. Chiou.                         C
C     See Directory: /mars/limbo/chiou/Beam.d/Timo.d.                 C
C     Element Derived From Timoshenko's "Strength of Materials" Book. C
C                                                                     C
C=====================================================================C
C                                                                     C
C     Francois M. Hemez - July 11th 1994 - Version 1.0                C
C                                                                     C
C=====================================================================C
C                                                                     C
C     Input/Output:                                                   C
C     -------------                                                   C
C                                                                     C
C     ELM     <input>   finite element number                         C
C     ESTIF   <output>  elemental stiffness matrix                    C
C     A       <input>   cross-sectional area                          C
C     E       <input>   Young modulus                                 C
C     EFRAME  <input>   elemental frame coordinates                   C
C     IX      <input>   moment of inertia around local [x] axis       C
C     IY      <input>   moment of inertia around local [y] axis       C
C     IZ      <input>   moment of inertia around local [z] axis       C
C     ALPHAY  <input>   shear deflection constant associated to [IY]  C
C     ALPHAZ  <input>   shear deflection constant associated to [IZ]  C
C     C1      <input>   non-uniform torsion constant                  C
C     NU      <input>   Poisson's ratio                               C
C     X       <input>   global X-coordinates of the two nodal points  C
C     Y       <input>   global Y-coordinates of the two nodal points  C
C     Z       <input>   global Z-coordinates of the two nodal points  C
C     FLAG    <input>   integer specifying whether to transform from  C
C                       local to global frame                         C
C                                                                     C
C=====================================================================C
C
C     ------------
C     DECLARATIONS
C     ------------
C
C.....GLOBAL VARIABLES
C
      real*8    A , E , nu , C1 , alphay , alphaz , Ix , Iy , Iz
      real*8    eframe(9) , estif(12,12) , x(2) , y(2) , z(2)
C
C.....LOCAL VARIABLES
C
      integer   i , j , k , l, flag
      real*8    T(9) , zero , half , Le(12,12)
      real*8    JJ , dx , dy , dz , length , G , one , two
      real*8    c11 , c22 , c33 , c12 , c13 , c23 , eps
      real*8    b , gammay , gammaz , twelve , six , xKe
      real*8    four , bendy , bendz , bendCy , bendCz
      real*8    locke(12,12) , s11 , s22 , s33 , s26 , s35
      real*8    s44 , s55 , s66 , s17 , s28 , s39 , s59 , s68
      real*8    s77 , s88 , s99 , s212 , s311 , s410 , s511
      real*8    s612 , s812 , s911 , s1010 , s1111 , s1212
      logical   ortho , specialcasey , specialcasez
C
C.....LOCAL VARIABLES FOR FRAME COMPUTATION
C
***   integer   index
***   real*8    normT , Id(3,3) , xx , yy , zz , xnorm
C
C     ----
C     DATA
C     ----
C
      data zero    /0.000000D+00/
      data half    /0.500000D+00/
      data one     /1.000000D+00/
      data two     /2.000000D+00/
      data four    /4.000000D+00/
      data six     /6.000000D+00/
      data twelve  /1.200000D+01/
C
C.....ZERO-MACHINE FOR FRAME ORTHONORMALITY CHECK
C
      data eps     /1.000000D-10/
C
C     -----
C     LOGIC
C     -----
C
C.....CLEAR THE OUTPUT STIFFNESS MATRIX
C
      do 1001 j=1,12
         do 1002 i=1,12
            estif(i,j) = zero
 1002    continue
 1001 continue
C
C.....CLEAR THE LOCAL ARRAYS
C
      do 1003 i=1,9
         T(i) = zero
 1003 continue
C
      do 1004 j=1,12
         do 1005 i=1,12
            locke(i,j) = zero
 1005    continue
 1004 continue
C
      do 1006 j=1,12
         do 1007 i=1,12
            Le(i,j) = zero
 1007    continue
 1006 continue
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
C.....CHECK IF MATERIAL PROPERTIES ARE POSITIVE OR ZERO
C
      if (      E.le.zero ) go to 200
      if (      A.lt.zero ) go to 200
      if (     nu.lt.-1 ) go to 200
      if (     Ix.lt.zero ) go to 200
      if (     Iy.lt.zero ) go to 200
      if (     Iz.lt.zero ) go to 200
      if ( alphay.lt.zero ) go to 200
      if ( alphaz.lt.zero ) go to 200
      if (     C1.lt.zero ) go to 200
C
C.....EXTRACT THE ROTATION MATRIX FROM [EFRAME]
C.....       [ x_X y_X z_X ]
C..... [T] = [ x_Y y_Y z_Y ]
C.....       [ x_Z y_Z z_Z ]
C
      do 2001 i=1,9
         T(i) = eframe(i)
 2001 continue
C
C.....CALCULATE THE ROTATION MATRIX IF THE INPUT ONE IS NILL
C.....ARCHTUNG!!! ASSUMES THAT THE BEAM ELEMENT IS SYMMETRIC
C.....(THAT IS: [Iy] = [Iz]) WITH CIRCULAR CROSS-SECTION SO
C.....THAT THE LOCAL [y] AND [z] AXES MAY BE ORIENTED
C.....ARBITRARILY IN THE PLANE NORMAL TO THE BEAM'S AXIS [x]
C
***   normT = zero
C
***   do 2002 i=1,9
***      normT = normT + T(i)*T(i)
*2002 continue
C
***   if ( ( normT.eq.zero ).and.( Iy.eq.Iz ) ) then
C
C.....INITIALIZE THE IDENTITY MATRIX
C
***   do 2003 i=1,3
***      do 2004 j=1,3
***         Id(i,j) = zero
*2004    continue
***      Id(i,i) = one
*2003 continue
C
C.....COMPUTE DIRECTION [x] STORED IN FIRST 3 ENTRIES OF [T]
C
***   T(1) = dx/length
***   T(2) = dy/length
***   T(3) = dz/length
C
C.....FIND WHICH AXIS HAS A MAXIMUM ANGLE WITH DIRECTION OF THE BEAM
C
***   index = 1
C
***   if ( abs(T(2)).lt.abs(T(1)).and.abs(T(2)).ne.zero ) then
***      index = 2
***   endif
C
***   if ( abs(T(3)).lt.abs(T(1)) ) then
***      if ( abs(T(3)).lt.abs(T(2)).and.abs(T(3)).ne.zero ) then
***         index = 3
***      endif
***   endif
C
C.....COMPUTE THE CROSS-PRODUCT [axis^x]
C
***   xx = Id(2,index)*T(3) - Id(3,index)*T(2)
***   yy = Id(3,index)*T(1) - Id(1,index)*T(3)
***   zz = Id(1,index)*T(2) - Id(2,index)*T(1)
C
C.....COMPUTE ITS NORM
C
***   xnorm = dsqrt( (xx*xx) + (yy*yy) + (zz*zz) )
C
C.....TREATMENT FOR A ZERO-NORM
C
***   if ( xnorm.eq.zero ) then
C
C.....THE BEAM IS IN THE [+/-X]-DIRECTION AND [y]=[Y] IS GIVEN
C
***      if ( index.eq.1 ) then
***         xx    = zero
***         yy    = one
***         zz    = zero
***         xnorm = one
***      endif
C
C.....THE BEAM IS IN THE [+/-Y]-DIRECTION AND [y]=[Z] IS GIVEN
C
***      if ( index.eq.2 ) then
***         xx    = zero
***         yy    = zero
***         zz    = one
***         xnorm = one
***      endif
C
C.....THE BEAM IS IN THE [+/-Z]-DIRECTION AND [y]=[X] IS GIVEN
C
***      if ( index.eq.3 ) then
***         xx    = one
***         yy    = zero
***         zz    = zero
***         xnorm = one
***      endif
C
C.....END OF EXCEPTION TREATMENT (IF [xnorm] IS ZERO)
C
***   endif
C
C.....NORMALIZE [y] AND STORE IN ENTRIES 4 TO 6 OF [T]
C
***   T(4) = xx/xnorm
***   T(5) = yy/xnorm
***   T(6) = zz/xnorm
C
C.....COMPUTE [z=x^y] AND STORE IN ENTRIES 7 TO 9 OF [T]
C
***   T(7) = T(2)*T(6) - T(3)*T(5)
***   T(8) = T(3)*T(4) - T(1)*T(6)
***   T(9) = T(1)*T(5) - T(2)*T(4)
C
C.....CHECK ORTHOGONALITY OF ROTATION MATRIX [T]
C
      c11 = T(1)*T(1) + T(2)*T(2) + T(3)*T(3)
      c22 = T(4)*T(4) + T(5)*T(5) + T(6)*T(6)
      c33 = T(7)*T(7) + T(8)*T(8) + T(9)*T(9)
      c12 = T(1)*T(4) + T(2)*T(5) + T(3)*T(6)
      c13 = T(1)*T(7) + T(2)*T(8) + T(3)*T(9)
      c23 = T(4)*T(7) + T(5)*T(8) + T(6)*T(9)
C
C
      ortho = ( abs(one-c11).le.eps.and.abs(one-c22).le.eps.and.
     $          abs(one-c33).le.eps.and.abs(    c12).le.eps.and.
     $          abs(    c13).le.eps.and.abs(    c23).le.eps    )
C
C.....END OF TREATMENT FOR ROTATION MATRIX OF A SYMMETRIC BEAM
C
***   endif
C
C.....ERROR-MESSAGE IF THE FRAME IS NOT AVAILABLE FOR A
C.....NON-SYMMETRIC TIMOSHENKO BEAM ELEMENT
C
***   if ( ( normT.eq.zero ).and.( Iy.ne.Iz ) ) go to 400
C
C.....ASSEMBLE THE TRANFORMATION MATRIX
C.....
C.....        [ [T]  0   0   0  ]
C.....        [  0  [T]  0   0  ]
C..... [Le] = [  0   0  [T]  0  ]
C.....        [  0   0   0  [T] ]
C
      do 2005 k=1,4
         do 2006 j=1,3
            do 2007 i=1,3
               Le(3*(k-1)+i,3*(k-1)+j) = T(3*(j-1)+i)
 2007       continue
 2006    continue
 2005 continue
C
C.....INITIALIZE THE VARIABLE [JJ] FOR AN ARBITRARY CROSS SECTION
C
      JJ = Ix
C
C.....INITIALIZE THE VARIABLE [JJ] FOR A CIRCULAR CROSS SECTION
C
***   JJ = half*A*A/acos(-one)
C
C.....INITIALIZE THE TRANSVERSE SHEAR MODULUS FOR ISOTROPIC MATERIAL
C
      G = E/(two*(one+nu))
C
C.....INITIALIZE THE ENTRIES OF THE LOCAL STIFFNESS MATRIX
C
      s11   = zero
      s22   = zero
      s33   = zero
      s26   = zero
      s35   = zero
      s44   = zero
      s55   = zero
      s66   = zero
      s17   = zero
      s28   = zero
      s39   = zero
      s59   = zero
      s68   = zero
      s77   = zero
      s88   = zero
      s99   = zero
      s212  = zero
      s311  = zero
      s410  = zero
      s511  = zero
      s612  = zero
      s812  = zero
      s911  = zero
      s1010 = zero
      s1111 = zero
      s1212 = zero
C
C.....INITIALIZE LOGICALS FOR PARTICULAR CASES
C
      specialcasey = .false.
      specialcasez = .false.
C
C.....INITIALIZE IF [ALPHA_Y] IS ZERO
C
      if ( alphay.eq.zero ) then
         s33          =  twelve*E*Iy/(length*length*length)
         s35          =    -six*E*Iy/(length*length)
         s55          =    four*E*Iy/length
         s39          = -twelve*E*Iy/(length*length*length)
         s59          =     six*E*Iy/(length*length)
         s99          =  twelve*E*Iy/(length*length*length)
         s311         =    -six*E*Iy/(length*length)
         s511         =     two*E*Iy/length
         s911         =     six*E*Iy/(length*length)
         s1111        =    four*E*Iy/length
         specialcasey = .true.
      endif
C
C.....INITIALIZE IF [ALPHA_Z] IS ZERO
C
      if ( alphaz.eq.zero ) then
         s22          =  twelve*E*Iz/(length*length*length)
         s26          =     six*E*Iz/(length*length)
         s66          =    four*E*Iz/length
         s28          = -twelve*E*Iz/(length*length*length)
         s68          =    -six*E*Iz/(length*length)
         s88          =  twelve*E*Iz/(length*length*length)
         s212         =     six*E*Iz/(length*length)
         s612         =     two*E*Iz/length
         s812         =    -six*E*Iz/(length*length)
         s1212        =    four*E*Iz/length
         specialcasez = .true.
      endif
C
C.....INITIALIZE IF [ALPHA_Y] IS NONZERO BUT [A] IS ZERO
C
      if ( (alphay.ne.zero).and.(A.eq.zero) ) then
         s55          =  E*Iy/length
         s511         = -E*Iy/length
         s1111        =  E*Iy/length
         specialcasey = .true.
      endif
C
C.....INITIALIZE IF [ALPHA_Z] IS NONZERO BUT [A] IS ZERO
C
      if ( (alphaz.ne.zero).and.(A.eq.zero) ) then
         s66          =  E*Iz/length
         s612         = -E*Iz/length
         s1212        =  E*Iz/length
         specialcasez = .true.
      endif
C
C.....INITIALIZE IF [I_Y] IS ZERO
C
      if ( Iy.eq.zero ) then
         specialcasey = .true.
      endif
C
C.....INITIALIZE IF [I_Z] IS ZERO
C
      if ( Iz.eq.zero ) then
         specialcasez = .true.
      endif
C
C.....INITIALIZE THE STIFFNESS ENTRIES IN ALL OTHER CASES
C
      s11 =  E*A/length
      s17 = -E*A/length
      s77 =  E*A/length
C
      if ( JJ.eq.zero ) then
         s44   = zero
         s410  = zero
         s1010 = zero
      else
         if ( C1.eq.zero ) then
            s44   =  G*JJ/length
            s410  = -G*JJ/length
            s1010 =  G*JJ/length
         else
            b     =  half*length*dsqrt(G*JJ/C1)
            s44   =  (G*JJ)/(length*(one-(tanh(b)/b)))
            s410  = -(G*JJ)/(length*(one-(tanh(b)/b)))
            s1010 =  (G*JJ)/(length*(one-(tanh(b)/b)))
         endif
      endif
C
      if ( .not.specialcasey ) then
         gammay = (G*A*length*length)/(twelve*E*Iy)
         bendy  = (alphay+four*gammay)/(twelve*gammay*(alphay+gammay))
         bendCy = (-alphay+two*gammay)/(twelve*gammay*(alphay+gammay))
         s33    =  G*A/(length*(alphay+gammay))
         s35    = -G*A/(two*(alphay+gammay))
         s55    =  bendy*G*A*length
         s39    = -G*A/(length*(alphay+gammay))
         s59    =  G*A/(two*(alphay+gammay))
         s99    =  G*A/(length*(alphay+gammay))
         s311   = -G*A/(two*(alphay+gammay))
         s511   =  bendCy*G*A*length
         s911   =  G*A/(two*(alphay+gammay))
         s1111  =  bendy*G*A*length
      endif
C
      if ( .not.specialcasez ) then
         gammaz = (G*A*length*length)/(twelve*E*Iz)
         bendz  = (alphaz+four*gammaz)/(twelve*gammaz*(alphaz+gammaz))
         bendCz = (-alphaz+two*gammaz)/(twelve*gammaz*(alphaz+gammaz))
         s22    =  G*A/(length*(alphaz+gammaz))
         s26    =  G*A/(two*(alphaz+gammaz))
         s66    =  bendz*G*A*length
         s28    = -G*A/(length*(alphaz+gammaz))
         s68    = -G*A/(two*(alphaz+gammaz))
         s88    =  G*A/(length*(alphaz+gammaz))
         s212   =  G*A/(two*(alphaz+gammaz))
         s612   =  bendCz*G*A*length
         s812   = -G*A/(two*(alphaz+gammaz))
         s1212  =  bendz*G*A*length
      endif
C
C.....INITIALIZE THE LOCAL STIFFNESS MATRIX
C
      locke( 1, 1) = s11
      locke( 2, 2) = s22
      locke( 3, 3) = s33
      locke( 2, 6) = s26
      locke( 3, 5) = s35
      locke( 4, 4) = s44
      locke( 5, 5) = s55
      locke( 6, 6) = s66
      locke( 1, 7) = s17
      locke( 2, 8) = s28
      locke( 3, 9) = s39
      locke( 5, 9) = s59
      locke( 6, 8) = s68
      locke( 7, 7) = s77
      locke( 8, 8) = s88
      locke( 9, 9) = s99
      locke( 2,12) = s212
      locke( 3,11) = s311
      locke( 4,10) = s410
      locke( 5,11) = s511
      locke( 6,12) = s612
      locke( 8,12) = s812
      locke( 9,11) = s911
      locke(10,10) = s1010
      locke(11,11) = s1111
      locke(12,12) = s1212
C
      do 3001 j=1,11
         do 3002 i=(j+1),12
            locke(i,j) = locke(j,i)
 3002    continue
 3001 continue
C
C.....ASSEMBLE THE OUTPUT ELEMENTAL STIFFNESS MATRIX
C.....
C..... [ESTIF] = [Le] * [LOCKE] * [Le]^T
C.....
C
C if flag equals 1, return the transformed element stiffness matrix
C else, return the local element stiffness matrix
C
      if(flag .eq. 1) then
        do 4001 l=1,12
           do 4002 k=1,12
              xKe = locke(k,l)
              do 4003 j=1,12
                 do 4004 i=1,12
                    estif(i,j) = estif(i,j) + Le(i,k)*xKe*Le(j,l)
 4004            continue
 4003         continue
 4002      continue
 4001   continue
      else
        do i=1,12
          do j=1,12
             estif(i,j) = locke(i,j)
          enddo
        enddo
      endif
C
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
      write(*,*) "*** FATAL ERROR in Routine MSTIF7 ***"
      write(*,*) "*** The Timoschenko Beam Element  ***"
      write(*,*) "*** Has Zero Length!              ***"
      write(*,*) "*** ... All Treatments Terminated ***"
      stop
C
C.....ERROR-MESSAGE IF A MATERIAL/GEOMETRICAL PROPERTY IS WRONG
C
  200 continue
      write(*,*) "*** FATAL ERROR in Routine MSTIF7      ***"
      write(*,*) "*** A Material or Geometrical Property ***"
      write(*,*) "*** is Not Positive: Check Input Data  ***"
      write(*,*) "*** ... All Treatments Terminated      ***"
      stop
C
C.....ERROR-MESSAGE IF THE ROTATION MATRIX IS NOT ORTHOGONAL
C
  300 continue
      write(*,*) "*** FATAL ERROR in Routine MSTIF7      ***"
      write(*,*) "*** The Rotation Matrix Computed For a ***"
      write(*,*) "*** Timoshenko Beam is Not Orthogonal! ***"
      write(*,*) "*** ... All Treatments Terminated      ***"
      stop
C
C.....ERROR-MESSAGE IF THE ROTATION MATRIX IS NOT AVAILABLE
C.....FOR A BEAM WITH NON-SYMMETRIC CROSS-SECTIONAL AREA
C
  400 continue
      write(*,*) "*** FATAL ERROR in Routine MSTIF7    ***"
      write(*,*) "*** There is No Rotation Matrix      ***"
      write(*,*) "*** Available For a Timoshenko Beam  ***"
      write(*,*) "*** With Non-Symmetric Cross-section ***"
      write(*,*) "*** ... All Treatments Terminated    ***"
      stop
C
C     ---
C     END
C     ---
C
      end
C=end of routine "MSTIF7"
C=======================C
