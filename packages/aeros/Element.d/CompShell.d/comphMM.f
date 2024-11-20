C=====================================================================C
C        1         2         3         4         5         6         7C
C234567890123456789012345678901234567890123456789012345678901234567890C
C=====================================================================C
      subroutine comphMM( elm , type , x    , y   , dm  ,
     $                    f   , rowm , colm , rot , Ph1 ,
     $                    Ph2 , Ph3  , khMM             )
C=====================================================================C
C
C     -----------
C     DECLARATION
C     -----------
C
C.....Global Variables
C
      integer   elm , type , rowm(9) , colm(9)
      real*8    x(3) , y(3) , dm(3,3) , khMM(18,18) , f
      real*8    rot(6,6) , Ph1(9,3) , Ph2(9,3) , Ph3(9,3)
C
C.....Local Variables
C
      integer   i , j , col , row
      real*8    zero , twicearea , area
      real*8    gauss11 , gauss21 , gauss31 , invweight1
      real*8    gauss12 , gauss22 , gauss32 , invweight2
      real*8    gauss13 , gauss23 , gauss33 , invweight3
      real*8    x12 , x21 , x23 , x32 , x31 , x13
      real*8    y12 , y21 , y23 , y32 , y31 , y13
      real*8    x0 , y0 , x10 , x20 , x30 , y10 , y20 , y30
      real*8    aa12 , aa23 , aa31 , ss12 , ss23 , ss31
      real*8    ca , caa12 , caa23 , caa31 , sum123 , sum456
      real*8    cax10 , cax20 , cax30 , cay10 , cay20 , cay30
      real*8    h1 , h2 , h3 , h4 , h5 , h6 , factor1 , factor2
      real*8    dmPh1t1 , dmPh1t2 , dmPh1t3
      real*8    dmPh2t1 , dmPh2t2 , dmPh2t3
      real*8    dmPh3t1 , dmPh3t2 , dmPh3t3
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
C.....CLEAR THE OUTPUT STIFFNESS MATRIX
C
      do 1001 j=1,18
         do 1002 i=1,18
            khMM(i,j) = zero
 1002    continue
 1001 continue
C
C.....CLEAR THE LOCAL MATRICES
C
      do 1003 j=1,9
         do 1004 i=1,6
            hmv(i,j) = zero
 1004    continue
         do 1005 i=1,6
            hqv(i,j) = zero
 1005    continue
 1003 continue
C
      do 1006 j=1,6
         do 1007 i=1,3
            B(i,j) = zero
 1007    continue
 1006 continue
C
      do 1008 i=1,6
         Z1(i) = zero
 1008 continue
C
      do 1009 i=1,6
         Z2(i) = zero
 1009 continue
C
      do 1010 i=1,6
         Z3(i) = zero
 1010 continue
C
      do 1011 j=1,3
         do 1012 i=1,6
            z1bt(i,j) = zero
 1012    continue
         do 1013 i=1,6
            z2bt(i,j) = zero
 1013    continue
         do 1014 i=1,6
            z3bt(i,j) = zero
 1014    continue
 1011 continue
C
      do 1015 j=1,3
         do 1016 i=1,9
            Ph1(i,j) = zero
 1016    continue
         do 1017 i=1,9
            Ph2(i,j) = zero
 1017    continue
         do 1018 i=1,9
            Ph3(i,j) = zero
 1018    continue
 1015 continue
C
C.....RETURN IF THE STIFFNESS FACTOR IS ZERO
C
      if ( f.eq.zero ) return
C
C.....CHECK IF THE STIFFNESS FACTOR [F] IS POSITIVE
C
      if ( f.lt.zero ) go to 100
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
      area = 0.50D+00*twicearea
C
C.....CHECK THE AREA (ERROR IF NOT POSITIVE)
C
      if ( twicearea.le.zero ) go to 200
C
C.....GET THE COORDINATES OF THE CENTROID OF THE TRIANGLE
C
      x0 = ( x(1) + x(2) + x(3) )/3.00D+00
      y0 = ( y(1) + y(2) + y(3) )/3.00D+00
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
      do 2001 j=1,9
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
 2001 continue
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
      do 3001 j=1,3
C
         do 3002 i=1,6
            z1bt(i,j) = Z1(i)*B(j,i)
 3002    continue
C
         do 3003 i=1,6
            z2bt(i,j) = Z2(i)*B(j,i)
 3003    continue
C
         do 3004 i=1,6
            z3bt(i,j) = Z3(i)*B(j,i)
 3004    continue
C
 3001 continue
C
C.....MULTIPLY LOCAL MATRICES [z1bt], [z2bt] AND [z3bt] BY THE SQUARE
C.....ROOT OF ( 9.0 * [f] * INTEGRATION_WEIGHT / AREA )
C
      factor1 = sqrt( 9.00D+00*f/(area*invweight1) )
      factor2 = sqrt( 9.00D+00*f/(area*invweight2) )
      factor3 = sqrt( 9.00D+00*f/(area*invweight3) )
C
      do 4001 j=1,3
C
         do 4002 i=1,6
            z1bt(i,j) = factor1*z1bt(i,j)
 4002    continue
C
         do 4003 i=1,6
            z2bt(i,j) = factor2*z2bt(i,j)
 4003    continue
C
         do 4004 i=1,6
            z3bt(i,j) = factor3*z3bt(i,j)
 4004    continue
C
 4001 continue
C
C.....FORM THE LOCAL PRODUCTS:
C.....[Ph1] = [HQV]^T * [Z1] * [B]^T = [HQV]^T * [z1bt]
C.....[Ph2] = [HQV]^T * [Z2] * [B]^T = [HQV]^T * [z2bt]
C.....[Ph3] = [HQV]^T * [Z3] * [B]^T = [HQV]^T * [z3bt]
C
      do 5001 i=1,9
C
         h1 = hqv(1,i)
         h2 = hqv(2,i)
         h3 = hqv(3,i)
         h4 = hqv(4,i)
         h5 = hqv(5,i)
         h6 = hqv(6,i)
C
         do 5002 j=1,3
            Ph1(i,j) = h1*z1bt(1,j) + h2*z1bt(2,j) + h3*z1bt(3,j)
     $               + h4*z1bt(4,j) + h5*z1bt(5,j) + h6*z1bt(6,j)
 5002    continue
C
         do 5003 j=1,3
            Ph2(i,j) = h1*z2bt(1,j) + h2*z2bt(2,j) + h3*z2bt(3,j)
     $               + h4*z2bt(4,j) + h5*z2bt(5,j) + h6*z2bt(6,j)
 5003    continue
C
         do 5004 j=1,3
            Ph3(i,j) = h1*z3bt(1,j) + h2*z3bt(2,j) + h3*z3bt(3,j)
     $               + h4*z3bt(4,j) + h5*z3bt(5,j) + h6*z3bt(6,j)
 5004    continue
C
 5001 continue
C
C.....ASSEMBLE THE OUTPUT STIFFNESS SUCH THAT:
C.....[khMM] = [Ph1]*[dm]*[Ph1]^T + [Ph2]*[dm]*[Ph2]^T + [Ph3]*[dm]*[Ph3]^T
C
      do 6001 i=1,9
C
      col     = colm(i)
C
      dmPh1t1 = dm(1,1)*Ph1(i,1) + dm(1,2)*Ph1(i,2) + dm(1,3)*Ph1(i,3)
      dmPh1t2 = dm(2,1)*Ph1(i,1) + dm(2,2)*Ph1(i,2) + dm(2,3)*Ph1(i,3)
      dmPh1t3 = dm(3,1)*Ph1(i,1) + dm(3,2)*Ph1(i,2) + dm(3,3)*Ph1(i,3)
C
      dmPh2t1 = dm(1,1)*Ph2(i,1) + dm(1,2)*Ph2(i,2) + dm(1,3)*Ph2(i,3)
      dmPh2t2 = dm(2,1)*Ph2(i,1) + dm(2,2)*Ph2(i,2) + dm(2,3)*Ph2(i,3)
      dmPh2t3 = dm(3,1)*Ph2(i,1) + dm(3,2)*Ph2(i,2) + dm(3,3)*Ph2(i,3)
C
      dmPh3t1 = dm(1,1)*Ph3(i,1) + dm(1,2)*Ph3(i,2) + dm(1,3)*Ph3(i,3)
      dmPh3t2 = dm(2,1)*Ph3(i,1) + dm(2,2)*Ph3(i,2) + dm(2,3)*Ph3(i,3)
      dmPh3t3 = dm(3,1)*Ph3(i,1) + dm(3,2)*Ph3(i,2) + dm(3,3)*Ph3(i,3)
C
      do 6002 j=1,i
         row           = rowm(j)
         khMM(row,col) = khMM(row,col)    + dmPh1t1*Ph1(j,1)
     $                 + dmPh1t2*Ph1(j,2) + dmPh1t3*Ph1(j,3)
     $                 + dmPh2t1*Ph2(j,1) + dmPh2t2*Ph2(j,2)
     $                 + dmPh2t3*Ph2(j,3) + dmPh3t1*Ph3(j,1)
     $                 + dmPh3t2*Ph3(j,2) + dmPh3t3*Ph3(j,3)
         khMM(col,row) = khMM(row,col)
 6002 continue
C
 6001 continue
C
C.....OUTPUT THE MATRIX PRIOR TO ROTATION (FOR DEBUGGING ONLY)
C
*     open(unit=90,file="khMM.m")
*     write(90,*) "khMM=["
*     do 991 i=1,18
* 991 write(90,9) (khMM(i,j),j=1,18)
*     write(90,*) "     ];"
*     close(90)
*   9 format(18(1x,E16.9))
C
C.....ROTATE THE OUTPUT STIFFNESS MATRIX
C
C      call compmrot( khMM , rot , rot , rot )
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
C.....ERROR-MESSAGE IF THE STIFFNESS FACTOR [F] IS NEGATIVE
C
  100 continue
      write(*,*) "*** FATAL ERROR in routine COMPHMM       ***"
      write(*,*) "*** The Stiffness Factor [f] is Negative ***"
      write(*,*) "*** Check the Calling Sequence:          ***"
      write(*,*) "*** Factor [f] Must be Positive or Zero  ***"
      write(*,*) "*** EXECUTION TERNINATED RIGHT HERE      ***"
      stop
C
C.....ERROR-MESSAGE IF THE TRIANGLE'S AREA IS NEGATIVE OR ZERO
C
  200 continue
      write(*,*) "*** FATAL ERROR in routine COMPHMM         ***"
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
C=end of routine "COMPHMM"
C========================C
