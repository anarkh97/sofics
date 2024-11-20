C=====================================================================C
C        1         2         3         4         5         6         7C
C234567890123456789012345678901234567890123456789012345678901234567890C
C=====================================================================C
      subroutine complay( nlayer , idlayer , mtlayer , x      , y     ,
     $                    z      , type    , ilayer  , cstbb  , cstmm ,
     $                    cstbm  , cstmb   , eframe  , aframe         )
C=====================================================================C
C                                                                     C
C     Perform =   Assembles the Four 3 by 3 Constitutive Matrices for C
C     ---------   a given Layer of a Multi-layer Composite Laminate.  C
C                                                                     C
C                                                                     C
C     Input/Output =                                                  C
C     --------------                                                  C
C     NLAYER  <input>  number of layers of the composite element      C
C     IDLAYER <input>  identificators for each layer                  C
C     MTLAYER <input>  material properties of each layer              C
C     X       <input>  triangular coordinates in local x-direction    C
C     Y       <input>  triangular coordinates in local y-direction    C
C     Z       <input>  triangular coordinates in local z-direction    C
C     TYPE    <input>  type of constitutive law                       C
C     ILAYER  <input>  layer number                                   C
C     CSTBB   <output> 3 by 3  bending- bending constitutive matrix   C
C     CSTMM   <output> 3 by 3 membrane-membrane constitutive matrix   C
C     CSTBM   <output> 3 by 3  bending-membrane constitutive matrix   C
C     CSTMB   <output> 3 by 3 membrane- bending constitutive matrix   C
C     EFRAME  <input>  coordinates of the elemental frame system      C
C     AFRAME  <input>  frame system for orienting the fibers          C
C                                                                     C
C                                                                     C
C     Computations =                                                  C
C     --------------                                                  C
C     Four different matrices [D] are assembled for the given layer   C
C     number [ilayer] of the composite shell element:                 C
C                                                                     C
C                             [ [D_mm]  [D_mb] ]                      C
C     [Constitutive_Matrix] = [                ]                      C
C            6 by 6           [ [D_bm]  [D_bb] ]                      C
C                                                                     C
C     where "b" and "m" stand for bending and membrane, respectively. C
C                                                                     C
C     The constitutive matrix [D_bb] relates the element's curvatures C
C     to the bending moments and the constitutive matrix [D_mm]       C
C     relates the element's normal efforts to the strains. Similarly, C
C     the constitutive matrices [D_bm] and [D_mb] couple the bending  C
C     and membrane effects:                                           C
C                                                                     C
C     [ M_x  ]            [ k_x  ]                                    C
C     [ M_y  ] = [D_bb] * [ k_y  ]                                    C
C     [ M_xy ]            [ k_xy ]                                    C
C                                                                     C
C     [ N_x  ]            [ e_x  ]                                    C
C     [ N_y  ] = [D_mm] * [ e_y  ]                                    C
C     [ N_xy ]            [ e_xy ]                                    C
C                                                                     C
C     [ M_x  ]            [ e_x  ]                                    C
C     [ M_y  ] = [D_bm] * [ e_y  ]                                    C
C     [ M_xy ]            [ e_xy ]                                    C
C                                                                     C
C     [ N_x  ]            [ k_x  ]                                    C
C     [ N_y  ] = [D_mb] * [ k_y  ]                                    C
C     [ N_xy ]            [ k_xy ]                                    C
C                                                                     C
C     The symmetry constraints are:                                   C
C                                                                     C
C     [D_bb]^T = [D_bb]                                               C
C     [D_mm]^T = [D_mm]                                               C
C     [D_bm]^T = [D_mb]                                               C
C                                                                     C
C                                                                     C
C     Caution =   It is assumed that the layer has a constant         C
C     ---------   thickness so that no numerical interpolation is     C
C                 required.                                           C
C                                                                     C
C=====================================================================C
C=Author  = Francois M. Hemez                                         C
C=Date    = June 10th, 1995                                           C
C=Version = 2.0                                                       C
C=Comment =                                                           C
C=====================================================================C
C
C     ------------
C     DECLARATIONS
C     ------------
C
C.....Global Variables
C
      integer   type , nlayer , idlayer(5,nlayer) , ilayer
      real*8    cstbb(3,3) , cstmm(3,3) , cstbm(3,3) , cstmb(3,3)
      real*8    x(3) , y(3) , z(3) , mtlayer(8,nlayer)
      real*8    eframe(3,3) , aframe(3,3)
C
C.....Local Variables
C
      integer   i , j , layernumber
      real*8    zero , one , pi , twopi , intthick
      real*8    E1 , E2 , nu12 , G12 , mu1 , mu2
      real*8    zsup , zinf , thetaF , thetaD , theta
      real*8    S11 , S12 , S13 , S22 , S23 , S33 , detS
      real*8    Q(3,3) , Qbar(3,3) , T(3,3) , R(3)
      real*8    invT(3,3) , costheta , sintheta
      real*8    qt1 , qt2 , qt3 , refvec(3)
      real*8    norm1 , norm2 , normref , proj1 , proj2
      real*8    cosine1 , cosine2 , orifiber(3)
C
C     ----
C     DATA
C     ----
C
      data zero /0.000000D+00/
      data one  /1.000000D+00/
C
C.....INITIALIZE THE MAIN DIAGONAL OF REUTER'S MATRIX
C
      data R    /1.000000D+00, 1.000000D+00, 0.500000D+00/
C
C     -----
C     LOGIC
C     -----
C
      pi    = acos(-one)
      twopi = 2.000000D+00*pi
C
C.....CHECK IF TYPE OF CONTITUTIVE LAW IS CORRECT
C
      if ( (type.ne.2).and.(type.ne.3) ) go to 100
C
C.....CLEAR THE 3 BY 3 CONSTITUTIVE MATRICES
C
      do 1001 j=1,3
         do 1002 i=1,3
            cstbb(i,j) = zero
 1002    continue
         do 1003 i=1,3
            cstmm(i,j) = zero
 1003    continue
         do 1004 i=1,3
            cstbm(i,j) = zero
 1004    continue
         do 1005 i=1,3
            cstmb(i,j) = zero
 1005    continue
 1001 continue
C
C.....INITIALIZE THE LAYER MATERIAL PROPERTIES
C
      layernumber = idlayer(3,ilayer)
C
      E1          = mtlayer(1,ilayer)
      E2          = mtlayer(2,ilayer)
      nu12        = mtlayer(3,ilayer)
      G12         = mtlayer(4,ilayer)
      mu1         = mtlayer(5,ilayer)
      mu2         = mtlayer(6,ilayer)
      thetaF      = mtlayer(8,ilayer)
C
C.....CHECK FOR OBVIOUS ERRORS IN THE INPUT DATA
C
      if ( (layernumber.lt.1).or.(layernumber.gt.nlayer) ) go to 200
C
      if (   E1.le.zero ) go to 300
      if (   E2.le.zero ) go to 300
C     if ( nu12.le.zero ) go to 300
      if (  G12.le.zero ) go to 300
      if (  mu1.lt.zero ) go to 300
      if (  mu2.lt.zero ) go to 300
C
C.....TRANSFORM ANGLE IN THE RANGE BETWEEN 0-360      
C
      if ( (thetaF.lt.zero).or.(thetaF.gt.360.00D+00) ) then 
        irot = thetaF / 360.0D0
        thetaF = thetaF - real(irot) * 360.0D0
        if (thetaF.lt.0.0D0) thetaF = 360.0D0 - thetaF
      endif
C
C.....TRANSFORM FROM DEGREE TO RADIAN THE ANGLE BETWEEN THE
C.....REFERENCE ORIENTATION VECTOR AND THE ORIENTATION OF THE FIBERS
C
      thetaF = (pi*thetaF)/180.00D+00
C
C.....SET THE REFERENCE VECTOR FOR ORIENTING THE FIBERS OF THE LAYER
C.....(ALWAYS TAKE THE FIRST VECTOR OF THE FRAME)
C
      refvec(1) = aframe(1,1)
      refvec(2) = aframe(2,1)
      refvec(3) = aframe(3,1)
C
C.....INITIALIZE THE UPPER AND LOWER [z] COORDINATES FOR THE LAYER
C
      zinf = -0.50D+00*mtlayer(7,ilayer)
      zsup =  0.50D+00*mtlayer(7,ilayer)
C
C.....PROJECT THE REFERENCE VECTOR INTO THE PLANE OF
C.....THE ELEMENT TO GET THE FIBER ORIENTATION VECTOR
C
      orifiber(1) = zero
      orifiber(2) = zero
      orifiber(3) = zero
C
      norm1   = zero
      norm2   = zero
      normref = zero
      proj1   = zero
      proj2   = zero
C
      do 2001 i=1,3
         norm1   =   norm1 + eframe(i,1)*eframe(i,1)
         norm2   =   norm2 + eframe(i,2)*eframe(i,2)
         normref = normref +   refvec(i)*refvec(i)
         proj1   =   proj1 + eframe(i,1)*refvec(i)
         proj2   =   proj2 + eframe(i,2)*refvec(i)
 2001 continue
C
      norm1   = sqrt(norm1)
      norm2   = sqrt(norm2)
      normref = sqrt(normref)
C
      if ( normref.eq.zero ) then
         cosine1 = one
         cosine2 = zero
      else
         cosine1 = proj1/(norm1*normref)
         cosine2 = proj2/(norm2*normref)
      endif
C
      do 2002 i=1,3
         orifiber(i) = cosine1*eframe(i,1) + cosine2*eframe(i,2)
 2002 continue
C
C.....CALCULATE THE ANGLE FROM THE [x] FRAME OF THE LOCAL
C.....COORDINATE SYSTEM TO THE REFERENCE ORIENTATION VECTOR
C
      proj1  = zero
      proj2  = zero
      thetaD = zero
C
      do 2003 i=1,3
         proj1 = proj1 + eframe(i,1)*orifiber(i)
         proj2 = proj2 + eframe(i,2)*orifiber(i)
 2003 continue
C
      if ( proj1.eq.zero ) then
         if ( proj2.eq.zero ) go to 500
         if ( proj2.gt.zero ) thetaD = 0.50D+00*pi
         if ( proj2.lt.zero ) thetaD = 1.50D+00*pi
      else
         thetaD = atan(proj2/proj1)
      endif
C
      if ( thetaD.lt.zero ) then
         thetaD = thetaD + twopi
      endif
C
      if ( (thetaD.lt.zero).or.(thetaD.gt.twopi) ) go to 600
C
C.....CALCULATE THE ANGLE FROM THE [x] FRAME OF THE LOCAL
C.....COORDINATE SYSTEM TO THE DIRECTION OF THE FIBER
C
      theta = thetaD + thetaF
C
      if ( theta.gt.twopi ) then
         theta = theta - twopi
      endif
C
C.....CALCULATE THE COMPLIANCE MATRIX [S] WHICH RELATES THE STRESSES
C.....[s1], [s2] AND [s12] TO THE STRAINS [e1], [e2] AND [e12] IN
C.....THE COORDINATE SYSTEM {1;2} ASSOCIATED WITH THE FIBER ORIENTATION
C
      S11 = one/E1
      S12 = -nu12/E1
      S13 = mu1/G12
      S22 = one/E2
      S23 = mu2/G12
      S33 = 1/G12
C
C.....CALCULATE THE DETERMINANT OF THE COMPLIANCE MATRIX
C
      detS = S33*(S11*S22-S12*S12) - S11*S23*S23 - S22*S13*S13
      detS = detS + 2.00D+00*S12*S13*S23
C
      if ( detS.eq.zero ) go to 700
C
C.....CALCULATE THE INVERSE OF THE COMPLIANCE MATRIX WHICH RELATES
C.....THE STRAINS [e1], [e2] AND [e12] TO THE STRESSES [s1], [s2]
C.....AND [s12] IN THE COORDINATE SYSTEM {1;2} OF THE FIBER ORIENTATION
C
      do 2202 j=1,3
         do 2203 i=1,3
            Q(i,j) = zero
 2203    continue
 2202 continue
C
      Q(1,1) = S22*S33 - S23*S23
      Q(1,2) = S13*S23 - S12*S33
      Q(1,3) = S12*S23 - S13*S22
C
      Q(2,1) = Q(1,2)
      Q(2,2) = S11*S33 - S13*S13
      Q(2,3) = S12*S13 - S11*S23
C
      Q(3,1) = Q(1,3)
      Q(3,2) = Q(2,3)
      Q(3,3) = S11*S22 - S12*S12
C
      do 2204 j=1,3
         do 2205 i=1,3
            Q(i,j) = Q(i,j)/detS
 2205    continue
 2204 continue
C
C.....INITIALIZE THE ROTATION MATRIX FROM THE FIBER COORDINATE
C.....SYSTEM {1;2} TO THE ELEMENT TRIANGULAR SYSTEM {x;y}
C
      do 2206 j=1,3
         do 2207 i=1,3
            T(i,j) = zero
 2207    continue
 2206 continue
C
      costheta = cos(theta)
      sintheta = sin(theta)
C
      T(1,1)   =  costheta*costheta
      T(1,2)   =  sintheta*sintheta
      T(1,3)   =  2.00D+00*costheta*sintheta
C
      T(2,1)   =  sintheta*sintheta
      T(2,2)   =  costheta*costheta
      T(2,3)   = -2.00D+00*costheta*sintheta
C
      T(3,1)   = -costheta*sintheta
      T(3,2)   =  costheta*sintheta
      T(3,3)   =  costheta*costheta - sintheta*sintheta
C
C.....COMPUTE THE INVERSE OF [T]:
C.....[invT] = inverse(diag[R]) * [T]^t * diag[R]
C
      do 2208 j=1,3
         do 2209 i=1,3
            invT(i,j) = zero
 2209    continue
 2208 continue
C
      do 2210 j=1,3
         do 2211 i=1,3
            invT(i,j) = T(j,i)*(R(j)/R(i))
 2211    continue
 2210 continue
C
C.....CALCULATE THE INVERSE OF THE COMPLIANCE MATRIX WHICH RELATES
C.....THE STRAINS [ex], [ey] AND [exy] TO THE STRESSES [sx], [sy]
C.....AND [sxy] IN THE TRIANGULAR COORDINATE SYSTEM {x;y}:
C.....[Qbar] = [invT] * [Q] * [invT]^t
C
      do 2212 j=1,3
         do 2213 i=1,3
            Qbar(i,j) = zero
 2213    continue
 2212 continue
C
      do 2214 j=1,3
C
         qt1 = Q(1,1)*invT(j,1) + Q(1,2)*invT(j,2) + Q(1,3)*invT(j,3)
         qt2 = Q(2,1)*invT(j,1) + Q(2,2)*invT(j,2) + Q(2,3)*invT(j,3)
         qt3 = Q(3,1)*invT(j,1) + Q(3,2)*invT(j,2) + Q(3,3)*invT(j,3)
C
         do 2215 i=1,3
            Qbar(i,j) = qt1*invT(i,1) + qt2*invT(i,2) + qt3*invT(i,3)
 2215    continue
C
 2214 continue
C
C     ------------------------------------------------
C      COMPOSITE MATERIAL WITH KNOWN LAYER PROPERTIES
C     NO COUPLING BETWEEN BENDING AND MEMBRANE EFFECTS
C      (NUMERICAL INTEGRATION THROUGH THE THICKNESS)
C     ------------------------------------------------
C
      if ( type.eq.2 ) then
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE BENDING
C
      intthick = (zsup*zsup*zsup - zinf*zinf*zinf)/3.00D+00
C
      do 3001 j=1,3
         do 3002 i=1,3
            cstbb(i,j) = cstbb(i,j) + Qbar(i,j)*intthick
 3002    continue
 3001 continue
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE MEMBRANE
C
      intthick = zsup - zinf
C
      do 3003 j=1,3
         do 3004 i=1,3
            cstmm(i,j) = cstmm(i,j) + Qbar(i,j)*intthick
 3004    continue
 3003 continue
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING BENDING-MEMBRANE
C
      do 3005 j=1,3
         do 3006 i=1,3
            cstbm(i,j) = zero
 3006    continue
 3005 continue
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING MEMBRANE-BENDING
C
      do 3007 j=1,3
         do 3008 i=1,3
            cstmb(i,j) = zero
 3008    continue
 3007 continue
C
C.....END OF TREATMENT FOR TYPE-2 CONSTITUTIVE LAW
C
      endif
C
C     --------------------------------------------------
C       COMPOSITE MATERIAL WITH KNOWN LAYER PROPERTIES
C     WITH COUPLING BETWEEN BENDING AND MEMBRANE EFFECTS
C       (NUMERICAL INTEGRATION THROUGH THE THICKNESS)
C     --------------------------------------------------
C
      if ( type.eq.3 ) then
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE BENDING
C
      intthick = (zsup*zsup*zsup - zinf*zinf*zinf)/3.00D+00
C
      do 4001 j=1,3
         do 4002 i=1,3
            cstbb(i,j) = cstbb(i,j) + Qbar(i,j)*intthick
 4002    continue
 4001 continue
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE MEMBRANE
C
      intthick = zsup - zinf
C
      do 4003 j=1,3
         do 4004 i=1,3
            cstmm(i,j) = cstmm(i,j) + Qbar(i,j)*intthick
 4004    continue
 4003 continue
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING BENDING-MEMBRANE
C.....(THE TRANSPOSED OF MATRIX [Qbar] IS TAKEN TO GET BEND-MEMB COUPLING)
C
      intthick = 0.50D+00*(zsup*zsup - zinf*zinf)
C
      do 4005 j=1,3
         do 4006 i=1,3
            cstbm(i,j) = cstbm(i,j) + Qbar(j,i)*intthick
 4006    continue
 4005 continue
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING MEMBRANE-BENDING
C
      intthick = 0.50D+00*(zsup*zsup - zinf*zinf)
C
      do 4007 j=1,3
         do 4008 i=1,3
            cstmb(i,j) = cstmb(i,j) + Qbar(i,j)*intthick
 4008    continue
 4007 continue
C
C.....END OF TREATMENT FOR TYPE-3 CONSTITUTIVE LAW
C
      endif
C
C     ------
C     RETURN
C     ------
C
      return
C
C     ------
C     FORMAT
C     ------
C
   91 format("*** Type Used is: ",I10,12x," ***")
C
C     ---------------
C     ERROR TREATMENT
C     ---------------
C
C.....ERROR-MESSAGE IF THE CONSTITUTIVE LAW IS INCORRECT
C
  100 continue
      write(*,*) "*** FATAL ERROR in Routine COMPLAY       ***"
      write(*,*) "*** Wrong Type of Constitutive Law       ***"
      write(*,91) type
      write(*,*) "*** Types Allowed are:                   ***"
      write(*,*) "*** 2 = given layers properties          ***"
      write(*,*) "***     (no coupling bending/membrane)   ***"
      write(*,*) "*** 3 = given layers properties          ***"
      write(*,*) "***     (with coupling bending/membrane) ***"
      write(*,*) "*** STOP ALL TREATMENTS RIGHT HERE       ***"
      stop
C
C.....ERROR-MESSAGE IF A LAYER IDENTIFICATION NUMBER IS NOT CORRECT
C
  200 continue
      write(*,*) "*** FATAL ERROR in Routine COMPLAY  ***"
      write(*,*) "*** The Local Layer Number is Not   ***"
      write(*,*) "*** Correct or Out of Bounds        ***"
      write(*,*) "*** STOP ALL TREATMENTS RIGHT HERE  ***"
      stop
C
C.....ERROR-MESSAGE IF A LAYER PROPERTY IS NOT CORRECT
C
  300 continue
      write(*,*) "*** FATAL ERROR in Routine COMPLAY     ***"
      write(*,*) "*** One of the Layer Material Property ***"
      write(*,*) "*** is Negative or Zero!               ***"
      write(*,*) "*** STOP ALL TREATMENTS RIGHT HERE     ***"
      stop
C
C.....ERROR-MESSAGE IF THE FIBER ANGLE IS OUT-OF-BOUNDS
C
  400 continue
      write(*,*) "*** FATAL ERROR in Routine COMPLAY    ***"
      write(*,*) "*** The Angle From the Reference      ***"
      write(*,*) "*** Direction to the Direction of     ***"
      write(*,*) "*** Fibers is Out-of-Bounds: it Must  ***"
      write(*,*) "*** be Within the Range 0-360 Degrees ***"
      write(*,*) "*** STOP ALL TREATMENTS RIGHT HERE    ***"
      stop
C
C.....ERROR-MESSAGE IF THE REFERENCE ORIENTATION IS BUGGY
C
  500 continue
      write(*,*) "*** FATAL ERROR in Routine COMPLAY   ***"
      write(*,*) "*** The Reference Orientation Vector ***"
      write(*,*) "*** is Parallel to the Two Inplane   ***"
      write(*,*) "*** and Orthogonal Local Frames!     ***"
      write(*,*) "*** STOP ALL TREATMENTS RIGHT HERE   ***"
      stop
C
C.....ERROR-MESSAGE IF THE ORIENTATION ANGLE IS OUT-OF-BOUNDS
C
  600 continue
      write(*,*) "*** FATAL ERROR in Routine COMPLAY    ***"
      write(*,*) "*** The Angle From the Local [x]      ***"
      write(*,*) "*** Axis of the Triangular Coordinate ***"
      write(*,*) "*** System to the Reference Direction ***"
      write(*,*) "*** is Out-of-Bounds: it Must be      ***"
      write(*,*) "*** Within the Range 0-2pi Radians    ***"
      write(*,*) "*** STOP ALL TREATMENTS RIGHT HERE    ***"
      stop
C
C.....ERROR-MESSAGE IF THE COMPLIANCE MATRIX IS SINGULAR
C
  700 continue
      write(*,*) "*** FATAL ERROR in Routine COMPLAY    ***"
      write(*,*) "*** The Compliance Matrix is Singular ***"
      write(*,*) "*** ... Check Material Properties ... ***"
      write(*,*) "*** STOP ALL TREATMENTS RIGHT HERE    ***"
      stop
C
C     ---
C     END
C     ---
C
      end
C=end of routine "COMPLAY"
C========================C
