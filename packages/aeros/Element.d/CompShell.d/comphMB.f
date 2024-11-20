C=====================================================================C
C        1         2         3         4         5         6         7C
C234567890123456789012345678901234567890123456789012345678901234567890C
C=====================================================================C
      subroutine comphMB( elm , type , x    , y   , dmb     , fb   ,
     $                    fm  , rowm , colb , rot , Lh1     , Lh2  ,
     $                    Lh3 , Ph1  , Ph2  , Ph3 , fastcal , khMB )
C=====================================================================C
C
C     -----------
C     DECLARATION
C     -----------
C
C.....Global Variables
C
      integer   elm , type , rowm(9) , colb(9)
      real*8    x(3) , y(3) , dmb(3,3) , khMB(18,18) , rot(6,6)
      real*8    Lh1(9,3) , Lh2(9,3) , Lh3(9,3)
      real*8    fb , fm , Ph1(9,3) , Ph2(9,3) , Ph3(9,3)
      logical   fastcal
C
C.....Local Variables
C
      integer   i , j , col , row
      real*8    zero , x0 , y0 , x1 , x2 , x3 , y1 , y2 , y3
      real*8    x21 , x32 , x13 , y21 , y32 , y13
      real*8    area , twicearea , x2Ap3 , dist21 , dist32 , dist13
      real*8    al1 , al2 , al3 , bl1 , bl2 , bl3 , factor
      real*8    s1 , s2 , s3 , q1 , q2 , q3 , q4 , q5 , q6
      real*8    q(6,9) , sq(3,3) , rm1(3,6) , rm2(3,6) , rm3(3,6)
      real*8    rm1tsqt(6,3) , rm2tsqt(6,3) , rm3tsqt(6,3)
      real*8    gauss11 , gauss21 , gauss31 , invweight1
      real*8    gauss12 , gauss22 , gauss32 , invweight2
      real*8    gauss13 , gauss23 , gauss33 , invweight3
      real*8    mult11 , mult12 , mult13
      real*8    mult21 , mult22 , mult23
      real*8    mult31 , mult32 , mult33
      real*8    x12 , x23 , x31 , y12 , y23 , y31
      real*8    x10 , x20 , x30 , y10 , y20 , y30
      real*8    aa12 , aa23 , aa31 , ss12 , ss23 , ss31
      real*8    ca , caa12 , caa23 , caa31 , sum123 , sum456
      real*8    cax10 , cax20 , cax30 , cay10 , cay20 , cay30
      real*8    h1 , h2 , h3 , h4 , h5 , h6 , factor1 , factor2
      real*8    hmv(6,9) , hqv(6,9) , B(3,6) , factor3
      real*8    Z1(6) , Z2(6) , Z3(6)
      real*8    z1bt(6,3) , z2bt(6,3) , z3bt(6,3)
C
C     ----
C     DATA
C     ----
C
      data zero       /0.000000D+00/
C
C.....DEFINE THE FIRST SET OF GAUSS INTEGRATION POINTS (W/ MID-POINT RULE)
C
      data gauss11    /0.000000D+00/
      data gauss21    /0.500000D+00/
      data gauss31    /0.500000D+00/
      data invweight1 /3.000000D+00/
C
C.....DEFINE THE SECOND SET OF GAUSS INTEGRATION POINTS (W/ MID-POINT RULE)
C
      data gauss12    /0.500000D+00/
      data gauss22    /0.000000D+00/
      data gauss32    /0.500000D+00/
      data invweight2 /3.000000D+00/
C
C.....DEFINE THE THIRD SET OF GAUSS INTEGRATION POINTS (W/ MID-POINT RULE)
C
      data gauss13    /0.500000D+00/
      data gauss23    /0.500000D+00/
      data gauss33    /0.000000D+00/
      data invweight3 /3.000000D+00/
C
C     -----
C     LOGIC
C     -----
C
C     --------------------------------
C
C     BASIC CHECKS AND INITIALIZATIONS
C
C     --------------------------------
C
C.....CLEAR THE OUTPUT STIFFNESS MATRIX
C
      do 1001 j=1,18
         do 1002 i=1,18
            khMB(i,j) = zero
 1002    continue
 1001 continue
C
C.....CLEAR THE LOCAL MATRICES FOR BENDING
C
      do 1003 j=1,9
         do 1004 i=1,6
            q(i,j) = zero
 1004    continue
 1003 continue
C
      do 1005 j=1,3
         do 1006 i=1,3
            sq(i,j) = zero
 1006    continue
 1005 continue
C
      do 1007 j=1,6
         do 1008 i=1,3
            rm1(i,j) = zero
 1008    continue
         do 1009 i=1,3
            rm2(i,j) = zero
 1009    continue
         do 1010 i=1,3
            rm3(i,j) = zero
 1010    continue
 1007 continue
C
      do 1011 j=1,3
         do 1012 i=1,6
            rm1tsqt(i,j) = zero
 1012    continue
         do 1013 i=1,6
            rm2tsqt(i,j) = zero
 1013    continue
         do 1014 i=1,6
            rm3tsqt(i,j) = zero
 1014    continue
 1011 continue
C
C.....CLEAR THE LOCAL MATRICES FOR MEMBRANE
C
      do 1015 j=1,9
         do 1016 i=1,6
            hmv(i,j) = zero
 1016    continue
         do 1017 i=1,6
            hqv(i,j) = zero
 1017    continue
 1015 continue
C
      do 1018 j=1,6
         do 1019 i=1,3
            B(i,j) = zero
 1019    continue
 1018 continue
C
      do 1020 i=1,6
         Z1(i) = zero
 1020 continue
C
      do 1021 i=1,6
         Z2(i) = zero
 1021 continue
C
      do 1022 i=1,6
         Z3(i) = zero
 1022 continue
C
      do 1023 j=1,3
         do 1024 i=1,6
            z1bt(i,j) = zero
 1024    continue
         do 1025 i=1,6
            z2bt(i,j) = zero
 1025    continue
         do 1026 i=1,6
            z3bt(i,j) = zero
 1026    continue
 1023 continue
C
C.....INITIALIZE THE INTEGRATED STRAIN-TO-DISPLACEMENT AND
C.....CURVATURE-TO-DISPLACEMENT MATRICES TI ZERO OR SKIP
C.....THEIR ASSEMBLY IF THEY ARE AVAILABLE ALREADY
C
      if ( fastcal ) then
C
         go to 900
C
      else
C
         do 1027 j=1,3
            do 1028 i=1,9
               Lh1(i,j) = zero
 1028       continue
            do 1029 i=1,9
               Lh2(i,j) = zero
 1029       continue
            do 1030 i=1,9
               Lh3(i,j) = zero
 1030       continue
 1027    continue
C
         do 1031 j=1,3
            do 1032 i=1,9
               Ph1(i,j) = zero
 1032       continue
            do 1033 i=1,9
               Ph2(i,j) = zero
 1033       continue
            do 1034 i=1,9
               Ph3(i,j) = zero
 1034       continue
 1031    continue
C
      endif
C
C.....RETURN IF THE STIFFNESS FACTOR FOR BENDING IS ZERO
C
      if ( fb.eq.zero ) return
C
C.....RETURN IF THE STIFFNESS FACTOR FOR MEMBRANE IS ZERO
C
      if ( fm.eq.zero ) return
C
C.....CHECK IF THE STIFFNESS FACTOR [FB] FOR BENDING IS POSITIVE
C
      if ( fb.lt.zero ) go to 100
C
C.....CHECK IF THE STIFFNESS FACTOR [FM] FOR MEMBRANE IS POSITIVE
C
      if ( fm.lt.zero ) go to 200
C
C     ------------------------------------------------------
C
C     FORM THE LOCAL MOMENT-CURVATURE MATRICES FOR BENDING
C
C             AND FOR THE THREE INTEGRATION POINTS
C
C     ------------------------------------------------------
C
C.....GET THE COORDINATES OF THE CENTROID OF THE TRIANGLE
C
      x0 = ( x(1) + x(2) + x(3) )/3.00D+00
      y0 = ( y(1) + y(2) + y(3) )/3.00D+00
C
C.....GET THE DISTANCES BETWEEN NODES 1, 2 AND 3 AND THE CENTROID
C
      x1 = x(1) - x0
      x2 = x(2) - x0
      x3 = x(3) - x0
C
      y1 = y(1) - y0
      y2 = y(2) - y0
      y3 = y(3) - y0
C
C.....GET THE DISTANCES BETWEEN NODAL POINT X- AND Y- COORDINATES
C
      x21 = x2 - x1
      x32 = x3 - x2
      x13 = x1 - x3
C
      y21 = y2 - y1
      y32 = y3 - y2
      y13 = y1 - y3
C
      x12 = -x21
      x23 = -x32
      x31 = -x13
C
      y12 = -y21
      y23 = -y32
      y31 = -y13
C
C.....CALCULATE TWICE THE AREA OF THE TRIANGLE
C
      twicearea = y21*x13 - x21*y13
C
C.....CALCULATE THE AREA OF THE TRIANGLE
C
      area = 0.50D+00*twicearea
C
C.....CHECK THE AREA (ERROR IF NOT POSITIVE)
C
      if ( twicearea.le.zero ) go to 300
C
C.....GET THE DISTANCES BETWEEN NODES 1-2, 2-3 AND 3-1
C
      dist21 = sqrt( x21*x21 + y21*y21 )
      dist32 = sqrt( x32*x32 + y32*y32 )
      dist13 = sqrt( x13*x13 + y13*y13 )
C
C.....GET THE SIDE PROJECTIONS
C
      bl1 = ( (x2-x3)*(x1-x3) + (y2-y3)*(y1-y3) )/(dist13*dist13)
      bl2 = ( (x3-x1)*(x2-x1) + (y3-y1)*(y2-y1) )/(dist21*dist21)
      bl3 = ( (x1-x2)*(x3-x2) + (y1-y2)*(y3-y2) )/(dist32*dist32)
C
      al1 = 1.00D+00 - bl1
      al2 = 1.00D+00 - bl2
      al3 = 1.00D+00 - bl3
C
C.....FORM THE MATRIX [Q] W/ SHAPE FUNCTION DERIVATIVES
C
      q(1,1) =  6.00D+00
      q(1,2) = -2.00D+00*y13
      q(1,3) =  2.00D+00*x13
      q(1,7) = -6.00D+00
      q(1,8) = -4.00D+00*y13
      q(1,9) =  4.00D+00*x13
C
      q(2,1) = -6.00D+00
      q(2,2) =  4.00D+00*y13
      q(2,3) = -4.00D+00*x13
      q(2,7) =  6.00D+00
      q(2,8) =  2.00D+00*y13
      q(2,9) = -2.00D+00*x13
C
      q(3,1) = -6.00D+00
      q(3,2) = -4.00D+00*y21
      q(3,3) =  4.00D+00*x21
      q(3,4) =  6.00D+00
      q(3,5) = -2.00D+00*y21
      q(3,6) =  2.00D+00*x21
C
      q(4,1) =  6.00D+00
      q(4,2) =  2.00D+00*y21
      q(4,3) = -2.00D+00*x21
      q(4,4) = -6.00D+00
      q(4,5) =  4.00D+00*y21
      q(4,6) = -4.00D+00*x21
C
      q(5,4) = -6.00D+00
      q(5,5) = -4.00D+00*y32
      q(5,6) =  4.00D+00*x32
      q(5,7) =  6.00D+00
      q(5,8) = -2.00D+00*y32
      q(5,9) =  2.00D+00*x32
C
      q(6,4) =  6.00D+00
      q(6,5) =  2.00D+00*y32
      q(6,6) = -2.00D+00*x32
      q(6,7) = -6.00D+00
      q(6,8) =  4.00D+00*y32
      q(6,9) = -4.00D+00*x32
C
C.....GET THE MATRIX [SQ] THAT REPRESENTS THE INVERSE OF THE MATRIX
C.....RELATING INSIDE CURVATURES WITH BOUNDARY CURVATURES
C
      x2Ap3   = (twicearea*twicearea*twicearea)
C
      sq(1,1) = ( -x21*y21*y32*y32 + x32*y21*y21*y32 )/x2Ap3
      sq(1,2) = (  x13*y13*y32*y32 - x32*y13*y13*y32 )/x2Ap3
      sq(1,3) = (  x21*y21*y13*y13 - x13*y21*y21*y13 )/x2Ap3
C
      sq(2,1) = (  x21*x32*x32*y21 - x21*x21*x32*y32 )/x2Ap3
      sq(2,2) = ( -x13*x32*x32*y13 + x13*x13*x32*y32 )/x2Ap3
      sq(2,3) = ( -x21*x13*x13*y21 + x21*x21*x13*y13 )/x2Ap3
C
      sq(3,1) = (  x21*x21*y32*y32 - x32*x32*y21*y21 )/x2Ap3
      sq(3,2) = ( -x13*x13*y32*y32 + x32*x32*y13*y13 )/x2Ap3
      sq(3,3) = ( -x21*x21*y13*y13 + x13*x13*y21*y21 )/x2Ap3
C
C.....MULTIPLY MATRIX [SQ] BY THE SQUARE ROOT OF THE STIFFNESS FACTOR
C.....AND THE SQUARE ROOT OF THE AREA DIVIDED BY THE WEIGHT OF THE
C.....NUMERICAL INTEGRATION (EQUAL TO THREE W/ MID-POINT RULE HERE)
C
      factor = sqrt(fb*area/3.00D+00)
C
      do 2001 j=1,3
         do 2002 i=1,3
            sq(i,j) = factor*sq(i,j)
 2002    continue
 2001 continue
C
C.....ESTIMATE THE MATRIX [RM] AT THE FIRST GAUSS INTEGRATION POINT
C
      rm1(1,1) = gauss31 + al1*gauss21 - ( 1.00D+00 + al1 )/3.00D+00
      rm1(1,2) = gauss11 + bl1*gauss21 - ( 1.00D+00 + bl1 )/3.00D+00
C
      rm1(2,3) = gauss11 + al2*gauss31 - ( 1.00D+00 + al2 )/3.00D+00
      rm1(2,4) = gauss21 + bl2*gauss31 - ( 1.00D+00 + bl2 )/3.00D+00
C
      rm1(3,5) = gauss21 + al3*gauss11 - ( 1.00D+00 + al3 )/3.00D+00
      rm1(3,6) = gauss31 + bl3*gauss11 - ( 1.00D+00 + bl3 )/3.00D+00
C
C.....ESTIMATE THE MATRIX [RM] AT THE SECOND GAUSS INTEGRATION POINT
C
      rm2(1,1) = gauss32 + al1*gauss22 - ( 1.00D+00 + al1 )/3.00D+00
      rm2(1,2) = gauss12 + bl1*gauss22 - ( 1.00D+00 + bl1 )/3.00D+00
C
      rm2(2,3) = gauss12 + al2*gauss32 - ( 1.00D+00 + al2 )/3.00D+00
      rm2(2,4) = gauss22 + bl2*gauss32 - ( 1.00D+00 + bl2 )/3.00D+00
C
      rm2(3,5) = gauss22 + al3*gauss12 - ( 1.00D+00 + al3 )/3.00D+00
      rm2(3,6) = gauss32 + bl3*gauss12 - ( 1.00D+00 + bl3 )/3.00D+00
C
C.....ESTIMATE THE MATRIX [RM] AT THE THIRD GAUSS INTEGRATION POINT
C
      rm3(1,1) = gauss33 + al1*gauss23 - ( 1.00D+00 + al1 )/3.00D+00
      rm3(1,2) = gauss13 + bl1*gauss23 - ( 1.00D+00 + bl1 )/3.00D+00
C
      rm3(2,3) = gauss13 + al2*gauss33 - ( 1.00D+00 + al2 )/3.00D+00
      rm3(2,4) = gauss23 + bl2*gauss33 - ( 1.00D+00 + bl2 )/3.00D+00
C
      rm3(3,5) = gauss23 + al3*gauss13 - ( 1.00D+00 + al3 )/3.00D+00
      rm3(3,6) = gauss33 + bl3*gauss13 - ( 1.00D+00 + bl3 )/3.00D+00
C
C.....COMPUTE THE MATRIX-MATRIX PRODUCT [rm1]^T*[sq]^T INTO MATRIX [rm1tsqt]
C.....COMPUTE THE MATRIX-MATRIX PRODUCT [rm2]^T*[sq]^T INTO MATRIX [rm2tsqt]
C.....COMPUTE THE MATRIX-MATRIX PRODUCT [rm3]^T*[sq]^T INTO MATRIX [rm3tsqt]
C
      do 3001 j=1,3
C
         s1 = sq(j,1)
         s2 = sq(j,2)
         s3 = sq(j,3)
C
         do 3002 i=1,6
            rm1tsqt(i,j) = s1*rm1(1,i)+ s2*rm1(2,i) + s3*rm1(3,i)
 3002    continue
C
         do 3003 i=1,6
            rm2tsqt(i,j) = s1*rm2(1,i)+ s2*rm2(2,i) + s3*rm2(3,i)
 3003    continue
C
         do 3004 i=1,6
            rm3tsqt(i,j) = s1*rm3(1,i)+ s2*rm3(2,i) + s3*rm3(3,i)
 3004    continue
C
 3001 continue
C
C.....ASSEMBLE MATRICES [Lh1], [Lh2] AND [Lh3] FOR THE THREE GAUSS POINTS:
C.....[Lh1] = [q]^T * [rm1]^T * [sq]^T = [q]^T * [rm1tsqt]
C.....[Lh2] = [q]^T * [rm2]^T * [sq]^T = [q]^T * [rm2tsqt]
C.....[Lh3] = [q]^T * [rm3]^T * [sq]^T = [q]^T * [rm3tsqt]
C
      do 4001 i=1,9
C
         q1 = q(1,i)
         q2 = q(2,i)
         q3 = q(3,i)
         q4 = q(4,i)
         q5 = q(5,i)
         q6 = q(6,i)
C
         do 4002 j=1,3
            Lh1(i,j) = q1*rm1tsqt(1,j) + q2*rm1tsqt(2,j)
     $               + q3*rm1tsqt(3,j) + q4*rm1tsqt(4,j)
     $               + q5*rm1tsqt(5,j) + q6*rm1tsqt(6,j)
 4002    continue
C
         do 4003 j=1,3
            Lh2(i,j) = q1*rm2tsqt(1,j) + q2*rm2tsqt(2,j)
     $               + q3*rm2tsqt(3,j) + q4*rm2tsqt(4,j)
     $               + q5*rm2tsqt(5,j) + q6*rm2tsqt(6,j)
 4003    continue
C
         do 4004 j=1,3
            Lh3(i,j) = q1*rm3tsqt(1,j) + q2*rm3tsqt(2,j)
     $               + q3*rm3tsqt(3,j) + q4*rm3tsqt(4,j)
     $               + q5*rm3tsqt(5,j) + q6*rm3tsqt(6,j)
 4004    continue
C
 4001 continue
C
C     --------------------------------------------------------
C
C     FORM THE LOCAL NORMAL FORCE-STRAIN MATRICES FOR MEMBRANE
C
C              AND FOR THE THREE INTEGRATION POINTS
C
C     --------------------------------------------------------
C
C.....GET THE DISTANCES BETWEEN NODES 1, 2 AND 3 AND THE CENTROID
C
      x10 = x(1) - x0
      x20 = x(2) - x0
      x30 = x(3) - x0
C
      y10 = y(1) - y0
      y20 = y(2) - y0
      y30 = y(3) - y0
C
C.....CALCULATE BASIC COEFFICIENTS FOR THE ASSEMBLY OF MATRIX [HMV]
C
      aa12  = 2.25D+00*( x30*x30 + y30*y30 )
      aa23  = 2.25D+00*( x10*x10 + y10*y10 )
      aa31  = 2.25D+00*( x20*x20 + y20*y20 )
C
      caa12 = 15.00D+00/(128.00D+00*aa12)
      caa23 = 15.00D+00/(128.00D+00*aa23)
      caa31 = 15.00D+00/(128.00D+00*aa31)
C
      ss12  = x12*x12 + y12*y12
      ss23  = x23*x23 + y23*y23
      ss31  = x31*x31 + y31*y31
C
      ca    = 3.00D+00/(16.00D+00*area)
C
      cax10 = ca*x10
      cax20 = ca*x20
      cax30 = ca*x30
C
      cay10 = ca*y10
      cay20 = ca*y20
      cay30 = ca*y30
C
C.....CONSTRUCT LOCAL MATRIX [HMV] W/ SHAPE FUNCTION DERIVATIVES
C
      hmv(1,1) = cay30*x32
      hmv(1,2) = cay30*y32
      hmv(1,3) = cay30*x13
      hmv(1,4) = cay30*y13
      hmv(1,5) = cay30*x21
      hmv(1,6) = cay30*y21
      hmv(1,7) = ( ss23 - ss31 + (2.40D+00*aa12) )*y30*caa12
     $         + 4.00D+00*area*x30*caa12
      hmv(1,8) = (9.00D+00/16.00D+00)*y30 - hmv(1,7)
      hmv(1,9) = (3.00D+00/16.00D+00)*y30
C
      hmv(2,1) = cay10*x32
      hmv(2,2) = cay10*y32
      hmv(2,3) = cay10*x13
      hmv(2,4) = cay10*y13
      hmv(2,5) = cay10*x21
      hmv(2,6) = cay10*y21
      hmv(2,7) = (3.00D+00/16.00D+00)*y10
      hmv(2,8) = ( ss31 - ss12 + (2.40D+00*aa23) )*y10*caa23
     $         + 4.00D+00*area*x10*caa23
      hmv(2,9) = (9.00D+00/16.00D+00)*y10 - hmv(2,8)
C
      hmv(3,1) = cay20*x32
      hmv(3,2) = cay20*y32
      hmv(3,3) = cay20*x13
      hmv(3,4) = cay20*y13
      hmv(3,5) = cay20*x21
      hmv(3,6) = cay20*y21
      hmv(3,7) = ( ss23 - ss12 + (2.40D+00*aa31) )*y20*caa31
     $         - 4.00D+00*area*x20*caa31
      hmv(3,8) = (3.00D+00/16.00D+00)*y20
      hmv(3,9) = (9.00D+00/16.00D+00)*y20 - hmv(3,7)
C
      hmv(4,1) = -cax30*x32
      hmv(4,2) = -cax30*y32
      hmv(4,3) = -cax30*x13
      hmv(4,4) = -cax30*y13
      hmv(4,5) = -cax30*x21
      hmv(4,6) = -cax30*y21
      hmv(4,7) = ( ss31 - ss23 - (2.40D+00*aa12) )*x30*caa12
     $         + 4.00D+00*area*y30*caa12
      hmv(4,8) = -(9.00D+00/16.00D+00)*x30 - hmv(4,7)
      hmv(4,9) = -(3.00D+00/16.00D+00)*x30
C
      hmv(5,1) = -cax10*x32
      hmv(5,2) = -cax10*y32
      hmv(5,3) = -cax10*x13
      hmv(5,4) = -cax10*y13
      hmv(5,5) = -cax10*x21
      hmv(5,6) = -cax10*y21
      hmv(5,7) = -(3.00D+00/16.00D+00)*x10
      hmv(5,8) = ( ss12 - ss31 - (2.40D+00*aa23) )*x10*caa23
     $         + 4.00D+00*area*y10*caa23
      hmv(5,9) = -(9.00D+00/16.00D+00)*x10 - hmv(5,8)
C
      hmv(6,1) = -cax20*x32
      hmv(6,2) = -cax20*y32
      hmv(6,3) = -cax20*x13
      hmv(6,4) = -cax20*y13
      hmv(6,5) = -cax20*x21
      hmv(6,6) = -cax20*y21
      hmv(6,7) = ( ss12 - ss23 - (2.40D+00*aa31) )*x20*caa31
     $         - 4.00D+00*area*y20*caa31
      hmv(6,8) = -(3.00D+00/16.00D+00)*x20
      hmv(6,9) = -(9.00D+00/16.00D+00)*x20 - hmv(6,7)
C
C.....CONSTRUCT LOCAL MATRIX [HQV] FROM [HMV]
C
      do 5001 j=1,9
C
         sum123   = (2.00D+00/9.00D+00)*(hmv(1,j)+hmv(2,j)+hmv(3,j))
C
         hqv(1,j) = sum123 - (4.00D+00/3.00D+00)*hmv(1,j)
         hqv(2,j) = sum123 - (4.00D+00/3.00D+00)*hmv(2,j)
         hqv(3,j) = sum123 - (4.00D+00/3.00D+00)*hmv(3,j)
C
         sum456   = (2.00D+00/9.00D+00)*(hmv(4,j)+hmv(5,j)+hmv(6,j))
C
         hqv(4,j) = sum456 - (4.00D+00/3.00D+00)*hmv(4,j)
         hqv(5,j) = sum456 - (4.00D+00/3.00D+00)*hmv(5,j)
         hqv(6,j) = sum456 - (4.00D+00/3.00D+00)*hmv(6,j)
C
 5001 continue
C
C.....FORM THE LOCAL MATRIX [B]
C
      B(1,1) =  y30
      B(2,1) =  zero
      B(3,1) = -x30
C
      B(1,2) =  y10
      B(2,2) =  zero
      B(3,2) = -x10
C
      B(1,3) =  y20
      B(2,3) =  zero
      B(3,3) = -x20
C
      B(1,4) =  zero
      B(2,4) = -x30
      B(3,4) =  y30
C
      B(1,5) =  zero
      B(2,5) = -x10
      B(3,5) =  y10
C
      B(1,6) =  zero
      B(2,6) = -x20
      B(3,6) =  y20
C
C.....CALCULATE THE DIAGONAL TRIANGULAR COORDINATE MATRIX [Z1]
C.....AT THE FIRST GAUSS INTEGRATION POINT (W/ MID-POINT RULE)
C
      Z1(1) = gauss21 - gauss11
      Z1(2) = gauss31 - gauss21
      Z1(3) = gauss11 - gauss31
      Z1(4) = gauss21 - gauss11
      Z1(5) = gauss31 - gauss21
      Z1(6) = gauss11 - gauss31
C
C.....CALCULATE THE DIAGONAL TRIANGULAR COORDINATE MATRIX [Z2]
C.....AT THE SECOND GAUSS INTEGRATION POINT (W/ MID-POINT RULE)
C
      Z2(1) = gauss22 - gauss12
      Z2(2) = gauss32 - gauss22
      Z2(3) = gauss12 - gauss32
      Z2(4) = gauss22 - gauss12
      Z2(5) = gauss32 - gauss22
      Z2(6) = gauss12 - gauss32
C
C.....CALCULATE THE DIAGONAL TRIANGULAR COORDINATE MATRIX [Z3]
C.....AT THE THIRD GAUSS INTEGRATION POINT (W/ MID-POINT RULE)
C
      Z3(1) = gauss23 - gauss13
      Z3(2) = gauss33 - gauss23
      Z3(3) = gauss13 - gauss33
      Z3(4) = gauss23 - gauss13
      Z3(5) = gauss33 - gauss23
      Z3(6) = gauss13 - gauss33
C
C.....FORM THE PRODUCTS [Z1]*[B]^T, [Z2]*[B]^T AND [Z3]*[B]^T
C.....IN MATRICES [z1bt], [z2bt] AND [z3bt], RESPECTIVELY
C
      do 5002 j=1,3
C
         do 5003 i=1,6
            z1bt(i,j) = Z1(i)*B(j,i)
 5003    continue
C
         do 5004 i=1,6
            z2bt(i,j) = Z2(i)*B(j,i)
 5004    continue
C
         do 5006 i=1,6
            z3bt(i,j) = Z3(i)*B(j,i)
 5006    continue
C
 5002 continue
C
C.....MULTIPLY LOCAL MATRICES [z1bt], [z2bt] AND [z3bt] BY THE SQUARE
C.....ROOT OF ( 9.0 * [fm] * INTEGRATION_WEIGHT / AREA )
C
      factor1 = sqrt( 9.00D+00*fm/(area*invweight1) )
      factor2 = sqrt( 9.00D+00*fm/(area*invweight2) )
      factor3 = sqrt( 9.00D+00*fm/(area*invweight3) )
C
      do 5007 j=1,3
C
         do 5008 i=1,6
            z1bt(i,j) = factor1*z1bt(i,j)
 5008    continue
C
         do 5009 i=1,6
            z2bt(i,j) = factor2*z2bt(i,j)
 5009    continue
C
         do 5010 i=1,6
            z3bt(i,j) = factor3*z3bt(i,j)
 5010    continue
C
 5007 continue
C
C.....FORM THE LOCAL PRODUCTS:
C.....[Ph1] = [HQV]^T * [Z1] * [B]^T = [HQV]^T * [z1bt]
C.....[Ph2] = [HQV]^T * [Z2] * [B]^T = [HQV]^T * [z2bt]
C.....[Ph3] = [HQV]^T * [Z3] * [B]^T = [HQV]^T * [z3bt]
C
      do 5011 i=1,9
C
         h1 = hqv(1,i)
         h2 = hqv(2,i)
         h3 = hqv(3,i)
         h4 = hqv(4,i)
         h5 = hqv(5,i)
         h6 = hqv(6,i)
C
         do 5012 j=1,3
            Ph1(i,j) = h1*z1bt(1,j) + h2*z1bt(2,j) + h3*z1bt(3,j)
     $               + h4*z1bt(4,j) + h5*z1bt(5,j) + h6*z1bt(6,j)
 5012    continue
C
         do 5013 j=1,3
            Ph2(i,j) = h1*z2bt(1,j) + h2*z2bt(2,j) + h3*z2bt(3,j)
     $               + h4*z2bt(4,j) + h5*z2bt(5,j) + h6*z2bt(6,j)
 5013    continue
C
         do 5014 j=1,3
            Ph3(i,j) = h1*z3bt(1,j) + h2*z3bt(2,j) + h3*z3bt(3,j)
     $               + h4*z3bt(4,j) + h5*z3bt(5,j) + h6*z3bt(6,j)
 5014    continue
C
 5011 continue
C
C     --------------------------------------------------------------------
C
C     STIFFNESS MATRIX ASSEMBLY FOR HIGHER ORDER MEMBRANE-BENDING COUPLING
C
C     --------------------------------------------------------------------
C
  900 continue
C
C.....ASSEMBLE THE OUTPUT STIFFNESS SUCH THAT:
C.....[khMB] = [Ph1]*[dmb]*[Lh1]^T + [Ph2]*[dmb]*[Lh2]^T + [Ph3]*[dmb]*[Lh3]^T
C
      do 6001 i=1,9
C
      col    = colb(i)
C
      mult11 = dmb(1,1)*Lh1(i,1) + dmb(1,2)*Lh1(i,2) + dmb(1,3)*Lh1(i,3)
      mult12 = dmb(2,1)*Lh1(i,1) + dmb(2,2)*Lh1(i,2) + dmb(2,3)*Lh1(i,3)
      mult13 = dmb(3,1)*Lh1(i,1) + dmb(3,2)*Lh1(i,2) + dmb(3,3)*Lh1(i,3)
C
      mult21 = dmb(1,1)*Lh2(i,1) + dmb(1,2)*Lh2(i,2) + dmb(1,3)*Lh2(i,3)
      mult22 = dmb(2,1)*Lh2(i,1) + dmb(2,2)*Lh2(i,2) + dmb(2,3)*Lh2(i,3)
      mult23 = dmb(3,1)*Lh2(i,1) + dmb(3,2)*Lh2(i,2) + dmb(3,3)*Lh2(i,3)
C
      mult31 = dmb(1,1)*Lh3(i,1) + dmb(1,2)*Lh3(i,2) + dmb(1,3)*Lh3(i,3)
      mult32 = dmb(2,1)*Lh3(i,1) + dmb(2,2)*Lh3(i,2) + dmb(2,3)*Lh3(i,3)
      mult33 = dmb(3,1)*Lh3(i,1) + dmb(3,2)*Lh3(i,2) + dmb(3,3)*Lh3(i,3)
C
      do 6002 j=1,9
         row           = rowm(j)
         khMB(row,col) = khMB(row,col)   + mult11*Ph1(j,1)
     $                 + mult12*Ph1(j,2) + mult13*Ph1(j,3)
     $                 + mult21*Ph2(j,1) + mult22*Ph2(j,2)
     $                 + mult23*Ph2(j,3) + mult31*Ph3(j,1)
     $                 + mult32*Ph3(j,2) + mult33*Ph3(j,3)
 6002 continue
C
 6001 continue
C
C.....OUTPUT THE MATRIX PRIOR TO ROTATION (FOR DEBUGGING)
C
*     open(unit=90,file="khMB.m")
*     write(90,*) "khMB=["
*     do 991 i=1,18
* 991 write(90,9) (khMB(i,j),j=1,18)
*     write(90,*) "     ];"
*     close(90)
*   9 format(18(1x,E16.9))
C
C.....ROTATE THE OUTPUT STIFFNESS MATRIX
C
C      call compmrot( khMB , rot , rot , rot )
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
C.....ERROR-MESSAGE IF THE STIFFNESS FACTOR [FB] (BENDING) IS NEGATIVE
C
  100 continue
      write(*,*) "*** FATAL ERROR in routine COMPHMB          ***"
      write(*,*) "*** The Stiffness Factor [fb] for Bending   ***"
      write(*,*) "*** is Negative: Check the Calling Sequence ***"
      write(*,*) "*** Factor [fb] Must be Positive or Zero    ***"
      write(*,*) "*** EXECUTION TERNINATED RIGHT HERE         ***"
      stop
C
C.....ERROR-MESSAGE IF THE STIFFNESS FACTOR [FM] (MEMBRANE) IS NEGATIVE
C
  200 continue
      write(*,*) "*** FATAL ERROR in routine COMPHMB          ***"
      write(*,*) "*** The Stiffness Factor [fm] for Membrane  ***"
      write(*,*) "*** is Negative: Check the Calling Sequence ***"
      write(*,*) "*** Factor [fm] Must be Positive or Zero    ***"
      write(*,*) "*** EXECUTION TERNINATED RIGHT HERE         ***"
      stop
C
C.....ERROR-MESSAGE IF THE TRIANGLE'S AREA IS NEGATIVE OR ZERO
C
  300 continue
      write(*,*) "*** FATAL ERROR in routine COMPHMB         ***"
      write(*,*) "*** The Triangle Area is Found Negative or ***"
      write(*,*) "*** Zero: Check the Nodal Point Numbering  ***"
      write(*,*) "*** ... Counterclockwise?                  ***"
      write(*,*) "*** EXECUTION TERNINATED RIGHT HERE        ***"
      stop
C
C     ---
C     END
C     ---
C
      end
C=end of routine "COMPHMB"
C========================C
