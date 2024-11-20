C=====================================================================C
C        1         2         3         4         5         6         7C
C234567890123456789012345678901234567890123456789012345678901234567890C
C=====================================================================C
      subroutine compvms( elm       , type     , numel     , maxstr   ,
     $                    maxgus    , E        , nu        , h        ,
     $                    globalX   , globalY  , globalZ   , globalU  ,
     $                    stress    , laysigm  , nttco     , nttly    ,
     $                    ncmpfr    , cmpco    , idlay     , mtlay    ,
     $                    cmpfr     , maxgly   , laysgid   , iatt     ,
     $                    ctyp      , catt     , cfrm      , nfely    ,
     $                    msize     , strainFlg, surface   , thrmStr1 ,
     $                    thrmStr2                                    )
C=====================================================================C
C                                                                     C
C     -----------------                                               C
C     V A R I A B L E S                                               C
C     -----------------                                               C
C                                                                     C
C     elm      <input>   Finite Element Number                        C
C     type     <input>   Finite Element Type                          C
C     numel    <input>   Number of Finite Elements in the Mesh        C
C     maxstr   <input>   Maximum Number of Stresses                   C
C     maxgus   <input>   Maximum Number of Gauss Points for Stresses  C
C     E        <input>   Young Modulus (for an Isotropic Element)     C
C     nu       <input>   Poisson's Ratio (for an Isotropic Element)   C
C     h        <input>   Constant Thickness of the Element            C
C     globalX  <input>   X- Nodal Coordinates                         C
C     globalY  <input>   Y- Nodal Coordinates                         C
C     globalZ  <input>   Z- Nodal Coordinates                         C
C     globalU  <input>   Global Displacements at the Nodal Joints     C
C     stress   <output>  Stresses (Von Mises Stress) of the Element   C
C     laysigm  <output>  Stresses (Von Mises) per Layers of the Shell C
C     nttco    <input>   Number of Attributes of Type-1               C
C     nttly    <input>   Total Number of Layers for the Elements      C
C     ncmpfr   <input>   Total Number of Composite Frames             C
C     cmpco    <input>   Constitutive Coefficients                    C
C     idlay    <input>   Layer Identificators                         C
C     mtlay    <input>   Material Properties of the Layers            C
C     cmpfr    <input>   Composite Frames                             C
C     maxgly   <input>   Maximum Number of Gauss Points (Layers)      C
C     laysgid  <input>   Identificators for the Stresses per Layer    C
C     iatt     <input>   Attribute Number of the Element              C
C     ctyp     <input>   Type of Constitutive Law (0, 1, 2, or 3)     C
C     catt     <input>   Addressing for Composite Material Properties C
C     cfrm     <input>   Composite Frame Number                       C
C     nfely    <input>   Total Number of Composite Layers in the Mesh C
C     msize    <input>   Maximum Length of Stress Output Array        C
C                                                                     C
C=====================================================================C
C=Author   = Francois M. Hemez                                        C
C=Date     = June 10th, 1995                                          C
C=Version  = 2.0                                                      C
C=Modified = K. H. Pierson                                            C
C=Date     = April 11, 1997                                           C
C=Reason   = Added stress calculations for sigmaxx, sigmayy, sigmaxy  C
C            and von mises stress at top, median and bottom surfaces  C
C            Also added strain calculations for epsilonxx, epsilonyy, C
C            epsilonzz, epsilonxy and an equivalent strain at top,    C
C            median and bottom surfaces.                              C
C=====================================================================C
C
C     -----------
C     DECLARATION
C     -----------
C
C.....GLOBAL VARIABLES
C
      integer    elm , type , numel , maxstr
      integer    maxgus , maxgly , nttco , nttly , ncmpfr
      integer    idlay(5,nttly) , iatt , ctyp , catt , cfrm,surface
      integer    laysgid(5,nfely) , nfely , msize, strainFlg
C
      real*8     ebar,epsxx,epsyy,epszz,epsxy,t2
      real*8     globalX(3) , globalY(3) , globalZ(3) , globalU(18)
      real*8     cmpco(36,nttco) , mtlay(12,nttly) , cmpfr(9,ncmpfr)
      real*8     stress(msize,maxstr,maxgus) , E , nu , h(3) 
      real*8     str(6),xp(3),yp(3),zp(3)
      real*8     xg(3),yg(3),zg(3)
      real*8     laysigm(nfely,maxstr,maxgly)
      real*8     thrmStr1, thrmStr2
C
C.....SET THE MAXIMUM NUMBER OF LAYERS OF THE ELEMENT
C
      integer    maxlayer
      parameter( maxlayer = 1000 )
C
C.....LOCAL VARIABLES
C
      integer    i , j , row , dimp , idcmp23(5,maxlayer)
      integer    rowb(9) , colb(9) , rowm(9) , colm(9) , nlayer
      integer    ilayer , layerpos
C
      real*8     area , twicearea , factor , alpha , clr , cqr
      real*8     llr(9,3) , lqr(9,3) , L(18,3) , P(18,3)
      real*8     reducedL(9,3) , reducedP(9,3)
      real*8     x(3) , y(3) , z(3) , rot(6,6) , localU(18)
      real*8     x21 , x32 , x13 , y21 , y32 , y13
      real*8     x12 , x23 , x31 , y12 , y23 , y31
      real*8     x0 , y0 , dist12 , dist23 , dist31
      real*8     c12 , c23 , c31 , s12 , s23 , s31
      real*8     cc12 , cc23 , cc31 , ss12 , ss23 , ss31
      real*8     cs12 , cs23 , cs31 , thick , x1 , x2
      real*8     zero , one , two , three , six , xkj , xsj
      real*8     cstbb(3,3) , cstmm(3,3) , eleM(3)
      real*8     cstbm(3,3) , cstmb(3,3) , eleN(3)
      real*8     mtcmp23(8,maxlayer) , cstcoef(36)
      real*8     elecrv(3) , elestr(3) , vonmises
      real*8     ups(3) , mds(3), lws(3) , upr , lwr , ups1 , ups2
      real*8     mdr,mds1,mds2,mdvms
      real*8     lws1 , lws2 , upvms , lwvms , appxh2
      real*8     eframe(3,3) , aframe(3,3)
C
      logical    detecterror , fastcal
C
C     ----
C     DATA
C     ----
C
      data zero  /0.000000D+00/
      data one   /1.000000D+00/
      data two   /2.000000D+00/
      data three /3.000000D+00/
      data six   /6.000000D+00/
      data xg    /1.0,0.0,0.0/
      data yg    /0.0,1.0,0.0/
      data zg    /0.0,0.0,1.0/
C
C.....INITIALIZE THE WEIGHTS FOR PURE BENDING CONTRIBUTION
C
      data clr   /0.000000D+00/
      data cqr   /1.000000D+00/
C
C.....INITIALIZE THE WEIGHT FOR PURE MEMBRANE CONTRIBUTION
C
      data alpha /1.500000D+00/
C
C.....INITIALIZE THE LOGICAL FOR ENFORCING FASTER COMPUTATIONS
C
      data fastcal /.true./
C
C     -----
C     LOGIC
C     -----
C
C.....CLEAR THE LOCAL MATRICES
C
C      do 1001 i=1,9
C         rowb(i) = 0
C         rowm(i) = 0
C         colb(i) = 0
C         colm(i) = 0
C 1001 continue
C
      do 1002 j=1,6
         do 1003 i=1,6
            rot(i,j) = zero
 1003    continue
 1002 continue
C
      do 1004 j=1,3
         do 1005 i=1,9
            llr(i,j) = zero
 1005    continue
         do 1006 i=1,9
            lqr(i,j) = zero
 1006    continue
         do 1007 i=1,9
            reducedL(i,j) = zero
 1007    continue
         do 1008 i=1,9
            reducedP(i,j) = zero
 1008    continue
         do 1009 i=1,18
            L(i,j) = zero
 1009    continue
         do 1010 i=1,18
            P(i,j) = zero
 1010    continue
 1004 continue
C
      do 1011 i=1,3
         x(i) = zero
         y(i) = zero
         z(i) = zero
 1011 continue
C
      do 1012 j=1,3
         do 1013 i=1,3
            cstbb(i,j) = zero
 1013    continue
         do 1014 i=1,3
            cstmm(i,j) = zero
 1014    continue
         do 1015 i=1,3
            cstbm(i,j) = zero
 1015    continue
         do 1016 i=1,3
            cstmb(i,j) = zero
 1016    continue
 1012 continue
C
      do 1017 i=1,36
         cstcoef(i) = zero
 1017 continue
C
      nlayer = 0
C
      do 1018 j=1,maxlayer
         do 1019 i=1,5
            idcmp23(i,j) = 0
 1019    continue
         do 1020 i=1,8
            mtcmp23(i,j) = zero
 1020    continue
 1018 continue
C
      do 1021 i=1,18
         localU(i) = zero
 1021 continue
C
      do 1022 i=1,3
         elecrv(i) = zero
 1022 continue
C
      do 1023 i=1,3
         elestr(i) = zero
 1023 continue
C
      do 1024 i=1,3
         eleM(i) = zero
 1024 continue
C
      do 1025 i=1,3
         eleN(i) = zero
 1025 continue
C
      do 1026 i=1,3
         ups(i) = zero
 1026 continue
C
      do 1027 i=1,3
         lws(i) = zero
 1027 continue
C
      do 1028 j=1,3
         do 1029 i=1,3
            aframe(i,j) = zero
 1029    continue
 1028 continue
C
      do 1030 j=1,3
         do 1031 i=1,3
            eframe(i,j) = zero
 1031    continue
 1030 continue
C
C.....CHECK THE TYPE OF CONSTITUTIVE LAW
C
      if (      (ctyp.ne.0)
     $     .and.(ctyp.ne.1)
     $     .and.(ctyp.ne.2)
     $     .and.(ctyp.ne.3) ) go to 100
C
C.....CHECK THE ADDRESSING IN STORAGE [CMPCO]
C
      if ( ctyp.eq.1 ) then
         if ( (catt.lt.1).or.(catt.gt.nttco) ) go to 200
      endif
C
C.....CHECK THE ADDRESSING IN ARRAYS [IDLAY] AND [MTLAY]
C
      if ( (ctyp.eq.2).or.(ctyp.eq.3) ) then
         if ( (catt.lt.1).or.(catt.gt.nttly) ) go to 300
      endif
C
C.....CHECK THE ADDRESSING IN ARRAY [CMPFR]
C
      if ( (ctyp.eq.1).or.(ctyp.eq.2).or.(ctyp.eq.3) ) then
         if ( (cfrm.lt.0).or.(cfrm.gt.ncmpfr) ) go to 400
      endif
C
C.....INITIALIZE THE CONSTITUTIVE COEFFICIENTS IN CASE OF A TYPE-1 LAW
C
      if ( ctyp.eq.1 ) then
         do 2001 i=1,36
            cstcoef(i) = cmpco(i,catt)
 2001    continue
      endif
C
C.....INITIALIZE THE LAYER PROPERTIES FOR TYPE-2 AND TYPE-3 LAWS
C
      if ( (ctyp.eq.2).or.(ctyp.eq.3) ) then
         nlayer = idlay(2,catt)
         if ( nlayer.gt.maxlayer ) go to 500
         do 2002 i=catt,(catt+nlayer-1)
            ilayer = idlay(3,i)
            if ( ilayer.gt.nlayer ) go to 600
            idcmp23(1,ilayer) = idlay(1,i)
            idcmp23(2,ilayer) = idlay(2,i)
            idcmp23(3,ilayer) = idlay(3,i)
            idcmp23(4,ilayer) = idlay(4,i)
            idcmp23(5,ilayer) = idlay(5,i)
            mtcmp23(1,ilayer) = mtlay(1,i)
            mtcmp23(2,ilayer) = mtlay(2,i)
            mtcmp23(3,ilayer) = mtlay(3,i)
            mtcmp23(4,ilayer) = mtlay(4,i)
            mtcmp23(5,ilayer) = mtlay(5,i)
            mtcmp23(6,ilayer) = mtlay(6,i)
            mtcmp23(7,ilayer) = mtlay(8,i)
            mtcmp23(8,ilayer) = mtlay(9,i)
 2002    continue
      endif
C
C.....CHECK CONSISTENCY OF BASIC PARAMETERS FOR TYPES 2 AND 3
C
      if ( (ctyp.eq.2).or.(ctyp.eq.3) ) then
         do 2003 i=1,nlayer
            if ( idcmp23(1,i).ne.idcmp23(1,1) ) go to 700
            if ( idcmp23(2,i).ne.nlayer       ) go to 700
 2003    continue
         do 2004 i=(nlayer+1),maxlayer
            if ( idcmp23(1,i).ne.0 ) go to 700
            if ( idcmp23(2,i).ne.0 ) go to 700
 2004    continue
      endif
C
C.....GET THE ARBITRARY FRAME FOR DEFINITION OF THE CONSTITUTIVE LAW
C
      if ( cfrm.eq.0 ) then
C
C.....INITIALIZE WITH THE IDENTITY IF THE FRAME NUMBER IS ZERO
C
         aframe(1,1) = one
         aframe(2,1) = zero
         aframe(3,1) = zero
         aframe(1,2) = zero
         aframe(2,2) = one
         aframe(3,2) = zero
         aframe(1,3) = zero
         aframe(2,3) = zero
         aframe(3,3) = one
C
      else
C
         aframe(1,1) = cmpfr(1,cfrm)
         aframe(2,1) = cmpfr(2,cfrm)
         aframe(3,1) = cmpfr(3,cfrm)
         aframe(1,2) = cmpfr(4,cfrm)
         aframe(2,2) = cmpfr(5,cfrm)
         aframe(3,2) = cmpfr(6,cfrm)
         aframe(1,3) = cmpfr(7,cfrm)
         aframe(2,3) = cmpfr(8,cfrm)
         aframe(3,3) = cmpfr(9,cfrm)
C
      endif
C
C.....INITIALIZE THE THICKNESS FOR A TYPE-0 CONSTITUTIVE LAW
C.....IT IS ASSUMED CONSTANT HERE
C
      if ( ctyp.eq.0 ) then
         thick = h(1)
      endif
C
C.....INITIALIZE THE THICKNESS FOR A TYPE-1 CONSTITUTIVE LAW
C.....TAKE THE DEFAULT THICKNESS ASSUMED CONSTANT
C.....IF ZERO, ESTIMATE THE CONSTANT THICKNESS USING THE
C.....FIRST COEFFICIENTS OF EXTENTIONAL AND BENDING STIFFNESSES
C
      if ( ctyp.eq.1 ) then
C
         thick = h(1)
C
         if ( thick.eq.zero ) then
            if ( cstcoef(1).eq.zero ) go to 800
            appxh2 = three*cstcoef(22)/cstcoef(1)
            if ( appxh2.le.zero ) go to 900
            thick = sqrt(appxh2)
         endif
C
      endif
C
C.....INITIALIZE THE THICKNESS FOR TYPE-2 AND TYPE-3 CONSTITUTIVE LAWS
C.....IT IS ASSUMED CONSTANT AND EQUAL TO THE SUM OF EACH LAYER'S THICKNESS
C
      if ( (ctyp.eq.2).or.(ctyp.eq.3) ) then
C
         thick = zero
C
         do 2005 ilayer=1,nlayer
            thick = thick + mtcmp23(7,ilayer)
 2005    continue
C
      endif
C
C     ----------------------------------
C     STEP 1
C     COMPUTE THE TRIANGULAR COORDINATES
C     ----------------------------------
C
C.....GET THE ELEMENT TRIANGULAR COORDINATES
C.....GET THE ROTATION MATRIX
C.....GET THE DEGREE OF FREEDOM POINTERS
C
      call compcrd2( elm  , ctyp , globalX , globalY , globalZ ,
     $               rot  , x    , y       , z       , rowb    ,
     $               colb , rowm , colm    , xp      , yp      ,
     $               zp                                        )
C
C.....GET THE ELEMENT LEVEL FRAME
C
      eframe(1,1) = rot(1,1)
      eframe(2,1) = rot(2,1)
      eframe(3,1) = rot(3,1)
      eframe(1,2) = rot(1,2)
      eframe(2,2) = rot(2,2)
      eframe(3,2) = rot(3,2)
      eframe(1,3) = rot(1,3)
      eframe(2,3) = rot(2,3)
      eframe(3,3) = rot(3,3)
C
C     ---------------------------------------------------
C     STEP 2
C     COMPUTE THE INTEGRATED CURVATURE-NODAL DISPLACEMENT
C     MATRIX FOR PURE BENDING (BASIC STIFFNESS MATRIX)
C     ---------------------------------------------------
C
C.....CHECK IF [CLR] AND [CQR] SATISFY THE CONSTRAINT [CLR]+[CQR]=1
C
      if ( (clr+cqr).ne.one ) go to 910
C
C.....GET THE DISTANCES BETWEEN NODAL POINT X- AND Y- COORDINATES
C
      x21 =  x(2) - x(1)
      x12 = -x21
      x32 =  x(3) - x(2)
      x23 = -x32
      x13 =  x(1) - x(3)
      x31 = -x13
      y21 =  y(2) - y(1)
      y12 = -y21
      y32 =  y(3) - y(2)
      y23 = -y32
      y13 =  y(1) - y(3)
      y31 = -y13
C
C.....CALCULATE TWICE THE AREA OF THE TRIANGLE
C
      twicearea = y21*x13 - x21*y13
C
C.....CALCULATE THE AREA OF THE TRIANGLE
C
      area = twicearea/two
C
C.....CHECK THE AREA (ERROR IF NOT POSITIVE)
C
      if ( twicearea.le.zero ) go to 920
C
C.....GET THE COORDINATES OF THE CENTROID OF THE TRIANGLE
C
      x0 = ( x(1) + x(2) + x(3) )/three
      y0 = ( y(1) + y(2) + y(3) )/three
C
C.....GET THE DISTANCES BETWEEN NODES 1-2, 2-3 AND 3-1
C
      dist12 = sqrt( x12*x12 + y12*y12 )
      dist23 = sqrt( x23*x23 + y23*y23 )
      dist31 = sqrt( x31*x31 + y31*y31 )
C
C.....ASSEMBLE THE LOCAL MATRIX [LLR] W/ SHAPE FUNCTION DERIVATIVES
C
      if ( clr.ne.zero ) then
C
         llr(3,1) =  y32/two
         llr(6,1) =  y13/two
         llr(9,1) =  y21/two
C
         llr(2,2) =  x32/two
         llr(5,2) =  x13/two
         llr(8,2) =  x21/two
C
         llr(2,3) = -y32/two
         llr(3,3) = -x32/two
         llr(5,3) = -y13/two
         llr(6,3) = -x13/two
         llr(8,3) = -y21/two
         llr(9,3) = -x21/two
C
      endif
C
C.....ASSEMBLE THE LOCAL MATRIX [LQR] W/ SHAPE FUNCTION DERIVATIVES
C
      if ( cqr.ne.zero ) then
C
         c12      = y21/dist12
         s12      = x12/dist12
         c23      = y32/dist23
         s23      = x23/dist23
         c31      = y13/dist31
         s31      = x31/dist31
C
         cc12     = c12*c12
         cc23     = c23*c23
         cc31     = c31*c31
         ss12     = s12*s12
         ss23     = s23*s23
         ss31     = s31*s31
         cs12     = c12*s12
         cs23     = c23*s23
         cs31     = c31*s31
C
         lqr(1,1) =  cs12 - cs31
         lqr(1,2) = -lqr(1,1)
         lqr(1,3) =  (cc31-ss31) - (cc12-ss12)
C
         lqr(2,1) =  ( cc12*x12 + cc31*x31 )/two
         lqr(2,2) =  ( ss12*x12 + ss31*x31 )/two
         lqr(2,3) =  ss12*y21 + ss31*y13
C
         lqr(3,1) = -( cc12*y21 + cc31*y13 )/two
         lqr(3,2) = -lqr(2,3)/two
         lqr(3,3) = -two*lqr(2,1)
C
         lqr(4,1) =  cs23 - cs12
         lqr(4,2) = -lqr(4,1)
         lqr(4,3) =  (cc12-ss12) - (cc23-ss23)
C
         lqr(5,1) =  ( cc12*x12 + cc23*x23 )/two
         lqr(5,2) =  ( ss12*x12 + ss23*x23 )/two
         lqr(5,3) =  ss12*y21 + ss23*y32
C
         lqr(6,1) = -( cc12*y21 + cc23*y32 )/two
         lqr(6,2) = -lqr(5,3)/two
         lqr(6,3) = -two*lqr(5,1)
C
         lqr(7,1) =  cs31 - cs23
         lqr(7,2) = -lqr(7,1)
         lqr(7,3) =  (cc23-ss23) - (cc31-ss31)
C
         lqr(8,1) =  ( cc23*x23 + cc31*x31 )/two
         lqr(8,2) =  ( ss23*x23 + ss31*x31 )/two
         lqr(8,3) =  ss23*y32 + ss31*y13
C
         lqr(9,1) = -( cc23*y32 + cc31*y13 )/two
         lqr(9,2) = -lqr(8,3)/two
         lqr(9,3) = -two*lqr(8,1)
C
      endif
C
C.....ASSEMBLE THE LOCAL MATRIX [L] AS [CLR]*[LLR] + [CQR]*[LQR]
C
      if ( clr.eq.zero ) then
         do 3001 j=1,3
            do 3002 i=1,9
               reducedL(i,j) = lqr(i,j)
 3002       continue
 3001    continue
      endif
C
      if ( cqr.eq.zero ) then
         do 3003 j=1,3
            do 3004 i=1,9
               reducedL(i,j) = llr(i,j)
 3004       continue
 3003    continue
      endif
C
      if ( (clr.ne.zero).and.(cqr.ne.zero) ) then
         do 3005 j=1,3
            do 3006 i=1,9
               reducedL(i,j) = clr*llr(i,j) + cqr*lqr(i,j)
 3006       continue
 3005    continue
      endif
C
C.....CHANGED THIS TO DIVIDE [L] BY THE AREA AS IS STANDARD FOR BENDING.
C     NOTE THAT WHEN COMPUTING THE STIFFNESS MATRIX L IS DEFINED TO 
C     BE [lqr] DIVIDED BY THE SQUARE ROOT OF THE AREA BECAUSE IN THIS 
C     CASE WE USE L TO COMPUTE K = 1/area*lqr*D*lqr^T = L*D*L^T
C
      factor = one/area
C     factor = sqrt(one/area)
C
      do 3007 j=1,3
         do 3008 i=1,9
            reducedL(i,j) = factor*reducedL(i,j)
 3008    continue
 3007 continue
C
C.....TRANSFORM MATRIX [L] IN FEM-LIKE DOF NUMBERING
C
      do 3009 i=1,9
         row      = rowb(i)
         L(row,1) = reducedL(i,1)
         L(row,2) = reducedL(i,2)
         L(row,3) = reducedL(i,3)
 3009 continue
C
C     -------------------------------------------------
C     STEP 3
C     COMPUTE THE INTEGRATED STRAIN-NODAL DISPLACEMENT
C     MATRIX FOR PURE MEMBRANE (BASIC STIFFNESS MATRIX)
C     -------------------------------------------------
C
C.....CHECK IF FACTOR [ALPHA] IS POSITIVE
C
      if ( alpha.le.zero ) go to 930
C
C.....ASSEMBLE THE MATRIX [P] W/ SHAPE FUNCTION DERIVATIVES
C
      reducedP(1,1) = y23
      reducedP(2,1) = zero
      reducedP(3,1) = y31
      reducedP(4,1) = zero
      reducedP(5,1) = y12
      reducedP(6,1) = zero
C
      reducedP(1,2) = zero
      reducedP(2,2) = x32
      reducedP(3,2) = zero
      reducedP(4,2) = x13
      reducedP(5,2) = zero
      reducedP(6,2) = x21
C
      reducedP(1,3) = x32
      reducedP(2,3) = y23
      reducedP(3,3) = x13
      reducedP(4,3) = y31
      reducedP(5,3) = x21
      reducedP(6,3) = y12
C
      dimp          = 6
C
      if ( alpha.ne.zero ) then
C
         reducedP(7,1) = y23*(y13-y21)*alpha/six
         reducedP(7,2) = x32*(x31-x12)*alpha/six
         reducedP(7,3) = (x31*y13-x12*y21)*alpha/three
C
         reducedP(8,1) = y31*(y21-y32)*alpha/six
         reducedP(8,2) = x13*(x12-x23)*alpha/six
         reducedP(8,3) = (x12*y21-x23*y32)*alpha/three
C
         reducedP(9,1) = y12*(y32-y13)*alpha/six
         reducedP(9,2) = x21*(x23-x31)*alpha/six
         reducedP(9,3) = (x23*y32-x31*y13)*alpha/three
C
         dimp          = 9
C
      endif
C
C.....CHANGED THIS TO DIVIDE [P] BY TWICE THE AREA AS IS
C     STANDARD FOR MEMBRANE.  THIS CHANGE ALLOWS THE
C     COMPOSITE SHELL, WHEN USED FOR ISOTROPIC MATERIALS,
C     TO GIVE CONSISTENT RESULTS WITH THE STANDARD SHELL ELEMENT.
C     Gregory W. Brown  12/04/00
C
      factor = (one/(twicearea))
C      factor = sqrt(one/(two*twicearea))
C
      do 4001 j=1,3
         do 4002 i=1,dimp
            reducedP(i,j) = factor*reducedP(i,j)
 4002    continue
 4001 continue
C
C.....TRANSFORM MATRIX [L] IN FEM-LIKE DOF NUMBERING
C
      do 4003 i=1,9
         row      = rowm(i)
         P(row,1) = reducedP(i,1)
         P(row,2) = reducedP(i,2)
         P(row,3) = reducedP(i,3)
 4003 continue
C
C     -------------------------------------------
C     STEP 4
C     ROTATE THE NODAL DISPLACEMENTS TO THE LOCAL
C     FRAME SYSTEM (LOCAL TO THE SHELL ELEMENT)
C     -------------------------------------------
C
      do 5001 i=1,6
         do 5002 j=1,6
            localU(i) = localU(i) + rot(j,i)*globalU(j)
 5002    continue
         do 5003 j=1,6
            localU(i+6) = localU(i+6) + rot(j,i)*globalU(j+6)
 5003    continue
         do 5004 j=1,6
            localU(i+12) = localU(i+12) + rot(j,i)*globalU(j+12)
 5004    continue
 5001 continue
C
C     --------------------------------------------------
C     STEP 5
C     COMPUTE THE ELEMENTAL STRAIN AND CURVATURE VECTORS
C     --------------------------------------------------
C
C.....ELEMENTAL CURVATURE COMPUTATION (1/radius)
C
      do 6001 i=1,3
         do 6002 j=1,18
            elecrv(i) = elecrv(i) + L(j,i)*localU(j)
 6002    continue
 6001 continue
C
C.....ELEMENTAL STRAIN COMPUTATION
C
      do 6003 i=1,3
         do 6004 j=1,18
            elestr(i) = elestr(i) + P(j,i)*localU(j)
 6004    continue
 6003 continue
C
C COMPUTE EQUIVALENT STRAIN (VON MISES)
C
        if(strainFlg .eq. 1) then
          t2 = 0.5*thick

          elecrv(1) = -t2*elecrv(1)
          elecrv(2) = -t2*elecrv(2)
          elecrv(3) = -t2*elecrv(3)

          if(surface .eq. 2) then
            epsxx = elestr(1)
            epsyy = elestr(2)
            epsxy = 0.5*(elestr(3))
          else if(surface .eq. 3) then
            epsxx = elestr(1) - elecrv(1)
            epsyy = elestr(2) - elecrv(2)
            epsxy = 0.5*(elestr(3) - elecrv(3))
          else
            epsxx = elestr(1) + elecrv(1)
            epsyy = elestr(2) + elecrv(2)
            epsxy = 0.5*(elestr(3) + elecrv(3))
          endif

          epszz = -nu/(1.0-nu)*(epsxx + epsyy)

          str(1) = epsxx
          str(2) = epsyy
          str(3) = epszz
          str(4) = epsxy
          str(5) = 0.0
          str(6) = 0.0

          call transform(xp,yp,zp,xg,yg,zg,str)

          str(4) = 2.0D0*str(4)
          str(5) = 2.0D0*str(5)
          str(6) = 2.0D0*str(6)

          do 101 i=1,6
            stress(elm,i,1) = str(i)
            stress(elm,i,2) = str(i)
            stress(elm,i,3) = str(i)
101	  continue

          call straineq(elecrv(1),elecrv(2),elecrv(3),
     &     elestr(1),elestr(2),elestr(3),nu,surface,ebar)

          stress(elm,7,1) = ebar
          stress(elm,7,2) = ebar
          stress(elm,7,3) = ebar

          return
        end if
C
C     Subtract off thermal strain portions
      elestr(1) = elestr(1) - thrmStr1
      elestr(2) = elestr(2) - thrmStr2
C
C
C     ----------------------------------------------
C     STEP 6
C     COMPUTE THE ELEMENTAL MOMENT AND FORCE
C     PER UNIT LENGTH VECTORS AND COMPUTE THE
C     CENTROIDAL VON MISES STRESS RESULTANT
C     INTEGRATED OVER THE THICKNESS OF THE COMPOSITE
C     ----------------------------------------------
C
C.....GET THE CONSTITUTIVE MATRIX FOR PURE BENDING
C
      call compcst( E       , thick   , nu     , cstcoef , nlayer ,
     $              idcmp23 , mtcmp23 , x      , y       , z      ,
     $              cstbb   , ctyp    , eframe , aframe  , "BB"   )
C
C.....GET THE CONSTITUTIVE MATRIX FOR PURE MEMBRANE
C
      call compcst( E       , thick   , nu     , cstcoef , nlayer ,
     $              idcmp23 , mtcmp23 , x      , y       , z      ,
     $              cstmm   , ctyp    , eframe , aframe  , "MM"   )
C
C.....GET THE CONSTITUTIVE MATRIX FOR BENDING-MEMBRANE COUPLING
C
      call compcst( E       , thick   , nu     , cstcoef , nlayer ,
     $              idcmp23 , mtcmp23 , x      , y       , z      ,
     $              cstbm   , ctyp    , eframe , aframe  , "BM"   )
C
C.....GET THE CONSTITUTIVE MATRIX FOR MEMBRANE-BENDING COUPLING
C
      if ( fastcal ) then
C
      cstmb(1,1) = cstbm(1,1)
      cstmb(2,1) = cstbm(1,2)
      cstmb(3,1) = cstbm(1,3)
      cstmb(1,2) = cstbm(2,1)
      cstmb(2,2) = cstbm(2,2)
      cstmb(3,2) = cstbm(2,3)
      cstmb(1,3) = cstbm(3,1)
      cstmb(2,3) = cstbm(3,2)
      cstmb(3,3) = cstbm(3,3)
C
      else
C
      call compcst( E       , thick   , nu     , cstcoef , nlayer ,
     $              idcmp23 , mtcmp23 , x      , y       , z      ,
     $              cstmb   , ctyp    , eframe , aframe  , "MB"   )
C
      endif
C
C.....CHECK THE CONSTITUTIVE MATRIX
C.....(COMMENTED OUT HERE - USED FOR DEBUGGING ONLY)
C
***   call compchk( elm , cstbb , cstmm , cstbm , cstmb )
C
C.....ESTIMATE THE ELEMENT'S MOMENTS PER UNIT LENGTH
C
      do 6101 j=1,3
         xkj = elecrv(j)
         xsj = elestr(j)
         do 6102 i=1,3
            eleM(i) = eleM(i) + cstbb(i,j)*xkj + cstbm(i,j)*xsj
 6102    continue
 6101 continue
C
C.....ESTIMATE THE ELEMENT'S FORCES PER UNIT LENGTH
C
      do 6103 j=1,3
         xkj = elecrv(j)
         xsj = elestr(j)
         do 6104 i=1,3
            eleN(i) = eleN(i) + cstmb(i,j)*xkj + cstmm(i,j)*xsj
 6104    continue
 6103 continue
C
C.....INITIALIZE VON MISES STRESS RESULTANT TO ZERO
C
      vonmises = zero
C
C.....ESTIMATE THE STRESSES ON THE UPPER SURFACE
C
      do 6105 i=1,3
         ups(i) = (eleN(i)/thick) - (six*eleM(i)/(thick*thick))
 6105 continue
C
C.....STORE SIGMAXX, SIGMAYY, SIGMAXY (UPPER)
C
      do 10 i=1,maxgus
        stress(elm,1,i) = ups(1)
        stress(elm,2,i) = ups(2)
        stress(elm,4,i) = ups(3)
 10   continue
C
C.....CALCULATE THE RADIUS OF MOHR CIRCLE ON THE UPPER SURFACE
C
      x1  = ((ups(1)-ups(2))/two)*((ups(1)-ups(2))/two)
      x2  = ups(3)*ups(3)
      upr = sqrt( x1 + x2 )
C
C.....CALCULATE THE PRINCIPAL STRESSES ON THE UPPER SURFACE
C
      ups1 = ((ups(1)+ups(2))/two) + upr
      ups2 = ((ups(1)+ups(2))/two) - upr
C
C.....CALCULATE VON MISES STRESS OF THE UPPER LAYER
C
      x1    = ups1*ups1 + ups2*ups2
      x2    = (ups1-ups2)*(ups1-ups2)
      upvms = sqrt( x1 + x2 )/sqrt(two)
C
C.....ESTIMATE THE STRESSES ON THE LOWER SURFACE
C
      do 6106 i=1,3
         lws(i) = (eleN(i)/thick) + (six*eleM(i)/(thick*thick))
 6106 continue
C
C.....CALCULATE THE RADIUS OF MOHR CIRCLE ON THE LOWER SURFACE
C
      x1  = ((lws(1)-lws(2))/two)*((lws(1)-lws(2))/two)
      x2  = lws(3)*lws(3)
      lwr = sqrt( x1 + x2 )
C
C.....CALCULATE THE PRINCIPAL STRESSES ON THE LOWER SURFACE
C
      lws1 = ((lws(1)+lws(2))/two) + lwr
      lws2 = ((lws(1)+lws(2))/two) - lwr
C
C.....CALCULATE VON MISES STRESS OF THE LOWER SURFACE
C
      x1    = lws1*lws1 + lws2*lws2
      x2    = (lws1-lws2)*(lws1-lws2)
      lwvms = sqrt( x1 + x2 )/sqrt(two)
C
C.....THE VON MISES STRESS IS THE MAXIMUM OF THE TWO
C
      vonmises = max(upvms,lwvms)
C
C.....IF UPPER SURFACE IS REQUESTED, RETURN upvms
C
      if(surface .eq. 1) vonmises = upvms
C
C.....IF lOWER SURFACE IS REQUESTED, RETURN lwvms
C
      if(surface .eq. 3) then
        vonmises = lwvms
        do 11 i=1,maxgus
          stress(elm,1,i) = lws(1)
          stress(elm,2,i) = lws(2)
          stress(elm,4,i) = lws(3)
 11     continue
      endif

C
C.....IF MEDIAN SURFACE IS REQUESTED, RETURN mds
C
      if(surface .eq. 2) then
        do i=1,3
          mds(i) = eleN(i)/thick
        end do

        do 12 i=1,maxgus
          stress(elm,1,i) = mds(1)
          stress(elm,2,i) = mds(2)
          stress(elm,4,i) = mds(3)
 12     continue

C
C.....CALCULATE THE RADIUS OF MOHR CIRCLE ON THE MEDIAN SURFACE
C
        x1  = ((mds(1)-mds(2))/two)*((mds(1)-mds(2))/two)
        x2  = mds(3)*mds(3)
        mdr = sqrt( x1 + x2 )
C
C.....CALCULATE THE PRINCIPAL STRESSES ON THE MEDIAN SURFACE
C
        mds1 = ((mds(1)+mds(2))/two) + mdr
        mds2 = ((mds(1)+mds(2))/two) - mdr
C
C.....CALCULATE VON MISES STRESS OF THE MEDIAN LAYER
C
        x1    = mds1*mds1 + mds2*mds2
        x2    = (mds1-mds2)*(mds1-mds2)
        mdvms = sqrt( x1 + x2 )/sqrt(two)

        vonmises = mdvms

      endif

C
C.....ROTATE LOCAL STRESSES TO GLOBAL
C
        str(1) = stress(elm,1,1)
        str(2) = stress(elm,2,1)
        str(3) = 0.0
        str(4) = stress(elm,4,1)
        str(5) = 0.0
        str(6) = 0.0

        call transform(xp,yp,zp,xg,yg,zg,str)

        do 102 i=1,6
          stress(elm,i,1) = str(i)
          stress(elm,i,2) = str(i)
          stress(elm,i,3) = str(i)
102     continue

        
C
C.....STORE VON MISES STRESS FOR EACH ONE OF THE GAUSS POINT
C
      do 6107 i=1,maxgus
         stress(elm,7,i) = vonmises
 6107 continue
C
C     -------------------------------------------------
C     STEP 7
C     COMPUTE THE ELEMENTAL MOMENT AND FORCE
C     PER UNIT LENGTH VECTORS PER LAYER
C     COMPUTE THE CENTROIDAL VON MISES STRESS RESULTANT
C     LAYER PER LAYER IN THE COMPOSITE SHELL ELEMENT
C     -------------------------------------------------
C
      if ( (ctyp.eq.2).or.(ctyp.eq.3) ) then
C
C.....LOOP ON THE LAYERS OF THE COMPOSITE SHELL ELEMENT
C
      do 7001 ilayer=1,nlayer
C
C.....EXTRACT THE LAYER POSITION
C
      layerpos = 0
C
      detecterror = .false.
C
      do 7002 i=1,nfely
         if ( laysgid(1,i).eq.elm ) then
            if ( laysgid(5,i).eq.ilayer ) then
               if ( detecterror ) go to 940
               layerpos    = i
               detecterror = .true.
            endif
         endif
 7002 continue
C
C.....INITIALIZE THE LAYER'S CONSTANT THICKNESS
C
      thick = mtcmp23(7,ilayer)
C
C.....GET THE CONSTITUTIVE LAW FOR THE LAYER
C
      do 7003 j=1,3
         do 7004 i=1,3
            cstbb(i,j) = zero
 7004    continue
         do 7005 i=1,3
            cstmm(i,j) = zero
 7005    continue
         do 7006 i=1,3
            cstbm(i,j) = zero
 7006    continue
         do 7007 i=1,3
            cstmb(i,j) = zero
 7007    continue
 7003 continue
C
      call complay( nlayer , idcmp23 , mtcmp23 , x      , y      ,
     $              z      , ctyp    , ilayer  , cstbb  , cstmm  ,
     $              cstbm  , cstmb   , eframe  , aframe          )
C
C.....ESTIMATE THE LAYER'S MOMENTS PER UNIT LENGTH
C
      do 7101 i=1,3
         eleM(i) = zero
 7101 continue
C
      do 7102 j=1,3
         xkj = elecrv(j)
         xsj = elestr(j)
         do 7103 i=1,3
            eleM(i) = eleM(i) + cstbb(i,j)*xkj + cstbm(i,j)*xsj
 7103    continue
 7102 continue
C
C.....ESTIMATE THE LAYER'S FORCES PER UNIT LENGTH
C
      do 7201 i=1,3
         eleN(i) = zero
 7201 continue
C
      do 7202 j=1,3
         xkj = elecrv(j)
         xsj = elestr(j)
         do 7203 i=1,3
            eleN(i) = eleN(i) + cstmb(i,j)*xkj + cstmm(i,j)*xsj
 7203    continue
 7202 continue
C
C.....INITIALIZE VON MISES STRESS RESULTANT TO ZERO
C
      vonmises = zero
C
C.....ESTIMATE THE STRESSES ON THE LAYER'S UPPER SURFACE
C
      do 7301 i=1,3
         ups(i) = zero
 7301 continue
C
      do 7302 i=1,3
         ups(i) = (eleN(i)/thick) - (six*eleM(i)/(thick*thick))
 7302 continue
C
C.....CALCULATE THE RADIUS OF MOHR CIRCLE ON THE LAYER'S UPPER SURFACE
C
      x1  = ((ups(1)-ups(2))/two)*((ups(1)-ups(2))/two)
      x2  = ups(3)*ups(3)
      upr = sqrt( x1 + x2 )
C
C.....CALCULATE THE PRINCIPAL STRESSES ON THE LAYER'S UPPER SURFACE
C
      ups1 = ((ups(1)+ups(2))/two) + upr
      ups2 = ((ups(1)+ups(2))/two) - upr
C
C.....CALCULATE VON MISES STRESS OF THE LAYER'S UPPER LAYER
C
      x1    = ups1*ups1 + ups2*ups2
      x2    = (ups1-ups2)*(ups1-ups2)
      upvms = sqrt( x1 + x2 )/sqrt(two)
C
C.....ESTIMATE THE STRESSES ON THE LAYER'S LOWER SURFACE
C
      do 7401 i=1,3
         lws(i) = zero
 7401 continue
C
      do 7402 i=1,3
         lws(i) = (eleN(i)/thick) + (six*eleM(i)/(thick*thick))
 7402 continue
C
C.....CALCULATE THE RADIUS OF MOHR CIRCLE ON THE LAYER'S LOWER SURFACE
C
      x1  = ((lws(1)-lws(2))/two)*((lws(1)-lws(2))/two)
      x2  = lws(3)*lws(3)
      lwr = sqrt( x1 + x2 )
C
C.....CALCULATE THE PRINCIPAL STRESSES ON THE LAYER'S LOWER SURFACE
C
      lws1 = ((lws(1)+lws(2))/two) + lwr
      lws2 = ((lws(1)+lws(2))/two) - lwr
C
C.....CALCULATE VON MISES STRESS OF THE LAYER'S LOWER LAYER
C
      x1    = lws1*lws1 + lws2*lws2
      x2    = (lws1-lws2)*(lws1-lws2)
      lwvms = sqrt( x1 + x2 )/sqrt(two)
C
C.....THE LAYER'S VON MISES STRESS IS THE MAXIMUM OF THE TWO
C
      vonmises = max(upvms,lwvms)
C
C.....STORE THE LAYER'S VON MISES STRESS IN 7TH POSITION
C.....(POSITIONS 1-6 ARE FOR Sx Sy Sz Sxy Sxz Syz) FOR EACH
C.....ONE OF THE [maxgly] GAUSS POINTS
C
      do 7501 i=1,maxgly
         laysigm(layerpos,7,i) = vonmises
 7501 continue
C
C.....END OF LOOP ON THE LAYERS OF THE COMPOSITE SHELL ELEMENT
C
 7001 continue
C
      endif
C
C     ------
C     RETURN
C     ------
C
      return
C
C     ---------------
C     ERROR TREATMENT
C     ---------------
C
C.....ERROR-MESSAGE IF UNKNOWN TYPE OF CONSTITUTIVE LAW
C
  100 continue
      write(*,*) "*** FATAL ERROR in Routine COMPVMS ***"
      write(*,*) "*** The Type of Constitutive Law   ***"
      write(*,*) "*** Must Either be 0, 1, 2, or 3.  ***"
      write(*,*) "*** Execution Terminated Here      ***"
      stop
C
C.....ERROR-MESSAGE IF ADDRESSING IN [CMPCO] IS NOT CORRECT
C
  200 continue
      write(*,*) "*** FATAL ERROR in Routine COMPVMS ***"
      write(*,*) "*** The Address in Array [CMPCO]   ***"
      write(*,*) "*** is Out-of-Bounds.              ***"
      write(*,*) "*** Execution Terminated Here      ***"
      stop
C
C.....ERROR-MESSAGE IF ADDRESSING IN [IDLAY]/[MTLAY] IS NOT CORRECT
C
  300 continue
      write(*,*) "*** FATAL ERROR in Routine COMPVMS ***"
      write(*,*) "*** The Address in Arrays [IDLAY]  ***"
      write(*,*) "*** and [MTLAY] is Out-of-Bounds.  ***"
      write(*,*) "*** Execution Terminated Here      ***"
      stop
C
C.....ERROR-MESSAGE IF ADDRESSING IN [CMPFR] IS NOT CORRECT
C
  400 continue
      write(*,*) "*** FATAL ERROR in Routine COMPVMS ***"
      write(*,*) "*** The Address in Array [CMPFR]   ***"
      write(*,*) "*** is Out-of-Bounds.              ***"
      write(*,*) "*** Execution Terminated Here      ***"
      stop
C
C.....ERROR-MESSAGE IF THE MAXIMUM NUMBER OF LAYERS HAS BEEN EXCEEDED
C
  500 continue
      write(*,*) "*** FATAL ERROR in Routine COMPVMS        ***"
      write(*,*) "*** The Maximum Number of Layers Has Been ***"
      write(*,*) "*** Exceeded: Boost Parameter [MAXLAYER]. ***"
      write(*,*) "*** Execution Terminated Here             ***"
      stop
C
C.....ERROR-MESSAGE IF THE TOTAL NUMBER OF LAYERS IS NOT CONSISTENT
C
  600 continue
      write(*,*) "*** FATAL ERROR in Routine COMPVMS      ***"
      write(*,*) "*** A Layer Number Exceeds the Total    ***"
      write(*,*) "*** Number of Layers Stored in [IDLAY]. ***"
      write(*,*) "*** Execution Terminated Here           ***"
      stop
C
C.....ERROR-MESSAGE IF BASIC PARAMETERS FOR TYPES 2 & 3 ARE INCONSISTENT
C
  700 continue
      write(*,*) "*** FATAL ERROR in Routine COMPVMS        ***"
      write(*,*) "*** The Element Numbers and Total Number  ***"
      write(*,*) "*** of Composite Layers Do Not Check Out. ***"
      write(*,*) "*** Execution Terminated Here             ***"
      stop
C
C.....ERROR-MESSAGE IF THE FIRST COEFFICIENT OF EXTENTION IS ZERO
C
  800 continue
      write(*,*) "*** FATAL ERROR in routine COMPVMS       ***"
      write(*,*) "*** The First Coefficient of Extensional ***"
      write(*,*) "*** Stiffness Cbb(1,1) is Equal to Zero! ***"
      write(*,*) "*** Can Not Estimate the Thickness...    ***"
      write(*,*) "*** EXECUTION TERMINATED RIGHT HERE      ***"
      stop
C
C.....ERROR-MESSAGE IF THE RATIO BETWEEN FIRST COEFFICIENTS OF
C.....BENDING STIFFNESS OVER EXTENTIONAL STIFFNESS IS NEGATIVE
C
  900 continue
      write(*,*) "*** FATAL ERROR in routine COMPVMS            ***"
      write(*,*) "*** The Ratio Between the First Coefficient   ***"
      write(*,*) "*** of Bending Stiffness and the First        ***"
      write(*,*) "*** Coefficient of Extensional Stiffness is   ***"
      write(*,*) "*** Negative or Zero: Can Not Take the Square ***"
      write(*,*) "*** Root and Estimate the Shell Thickness.    ***"
      write(*,*) "*** EXECUTION TERMINATED RIGHT HERE           ***"
      stop
C
C.....ERROR-MESSAGE IF [CLR]+[CQR] IS DIFFERENT FROM ONE
C
  910 continue
      write(*,*) "*** FATAL ERROR in routine COMPVMS      ***"
      write(*,*) "*** The Factors [clr] and [cqr] Violate ***"
      write(*,*) "*** the Constraint [clr]+[cqr]=1:       ***"
      write(*,*) "*** Check the Calling Sequence.         ***"
      write(*,*) "*** EXECUTION TERMINATED RIGHT HERE     ***"
      stop
C
C.....ERROR-MESSAGE IF THE TRIANGLE'S AREA IS NEGATIVE OR ZERO
C
  920 continue
      write(*,*) "*** FATAL ERROR in routine COMPVMS         ***"
      write(*,*) "*** The Triangle Area is Found Negative or ***"
      write(*,*) "*** Zero: Check the Nodal Point Numbering  ***"
      write(*,*) "*** ... Counterclockwise?                  ***"
      write(*,*) "*** EXECUTION TERMINATED RIGHT HERE        ***"
      stop
C
C.....ERROR-MESSAGE IF [ALPHA] IS NEGATIVE
C
  930 continue
      write(*,*) "*** FATAL ERROR in routine COMPVMS   ***"
      write(*,*) "*** The Factor [alpha] is Negative   ***"
      write(*,*) "*** Check the Data Sequence: Factor  ***"
      write(*,*) "*** [alpha] Must be Positive or Zero ***"
      write(*,*) "*** EXECUTION TERMINATED RIGHT HERE  ***"
      stop
C
C.....ERROR-MESSAGE IF THE LAYER NUMBER IS NOT UNIQUE
C
  940 continue
      write(*,*) "*** FATAL ERROR in routine COMPVMS       ***"
      write(*,*) "*** The Same Layer Number Is Encountered ***"
      write(*,*) "*** More than Once for a Given Element   ***"
      write(*,*) "*** Can Not Happen: Check Input File.    ***"
      write(*,*) "*** EXECUTION TERMINATED RIGHT HERE      ***"
      stop
C
C     ---
C     END
C     ---
C
      end
C=end of routine "COMPVMS"
C========================C
