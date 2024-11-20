C=====================================================================C
C        1         2         3         4         5         6         7C
C234567890123456789012345678901234567890123456789012345678901234567890C
C=====================================================================C
      subroutine comphBB( elm , type , x    , y   , db  ,
     $                    f   , rowb , colb , rot , Lh1 ,
     $                    Lh2 , Lh3  , khBB             )
C=====================================================================C
C
C     -----------
C     DECLARATION
C     -----------
C
C.....Global Variables
C
      integer   elm , type , rowb(9) , colb(9)
      real*8    x(3) , y(3) , db(3,3) , khBB(18,18) , f
      real*8    Lh1(9,3) , Lh2(9,3) , Lh3(9,3) , rot(6,6)
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
      real*8    gauss11 , gauss21 , gauss31
      real*8    gauss12 , gauss22 , gauss32
      real*8    gauss13 , gauss23 , gauss33
      real*8    dbLh1t1 , dbLh1t2 , dbLh1t3
      real*8    dbLh2t1 , dbLh2t2 , dbLh2t3
      real*8    dbLh3t1 , dbLh3t2 , dbLh3t3
C
C     ----
C     DATA
C     ----
C
      data zero    /0.000000D+00/
C
C.....DEFINE THE FIRST SET OF GAUSS INTEGRATION POINTS (W/ MID-POINT RULE)
C
      data gauss11 /0.000000D+00/
      data gauss21 /0.500000D+00/
      data gauss31 /0.500000D+00/
C
C.....DEFINE THE SECOND SET OF GAUSS INTEGRATION POINTS (W/ MID-POINT RULE)
C
      data gauss12 /0.500000D+00/
      data gauss22 /0.000000D+00/
      data gauss32 /0.500000D+00/
C
C.....DEFINE THE THIRD SET OF GAUSS INTEGRATION POINTS (W/ MID-POINT RULE)
C
      data gauss13 /0.500000D+00/
      data gauss23 /0.500000D+00/
      data gauss33 /0.000000D+00/
C
C     -----
C     LOGIC
C     -----
C
C.....CLEAR THE OUTPUT STIFFNESS MATRIX
C
      do 1001 j=1,18
         do 1002 i=1,18
            khBB(i,j) = zero
 1002    continue
 1001 continue
C
C.....CLEAR THE LOCAL MATRICES
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
      do 1015 j=1,3
         do 1016 i=1,9
            Lh1(i,j) = zero
 1016    continue
         do 1017 i=1,9
            Lh2(i,j) = zero
 1017    continue
         do 1018 i=1,9
            Lh3(i,j) = zero
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
      factor = sqrt(f*area/3.00D+00)
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
C.....ASSEMBLE THE OUTPUT STIFFNESS SUCH THAT:
C.....[khBB] = [Lh1]*[db]*[Lh1]^T + [Lh2]*[db]*[Lh2]^T + [Lh3]*[db]*[Lh3]^T
C
      do 5001 i=1,9
C
      col     = colb(i)
C
      dbLh1t1 = db(1,1)*Lh1(i,1) + db(1,2)*Lh1(i,2) + db(1,3)*Lh1(i,3)
      dbLh1t2 = db(2,1)*Lh1(i,1) + db(2,2)*Lh1(i,2) + db(2,3)*Lh1(i,3)
      dbLh1t3 = db(3,1)*Lh1(i,1) + db(3,2)*Lh1(i,2) + db(3,3)*Lh1(i,3)
C
      dbLh2t1 = db(1,1)*Lh2(i,1) + db(1,2)*Lh2(i,2) + db(1,3)*Lh2(i,3)
      dbLh2t2 = db(2,1)*Lh2(i,1) + db(2,2)*Lh2(i,2) + db(2,3)*Lh2(i,3)
      dbLh2t3 = db(3,1)*Lh2(i,1) + db(3,2)*Lh2(i,2) + db(3,3)*Lh2(i,3)
C
      dbLh3t1 = db(1,1)*Lh3(i,1) + db(1,2)*Lh3(i,2) + db(1,3)*Lh3(i,3)
      dbLh3t2 = db(2,1)*Lh3(i,1) + db(2,2)*Lh3(i,2) + db(2,3)*Lh3(i,3)
      dbLh3t3 = db(3,1)*Lh3(i,1) + db(3,2)*Lh3(i,2) + db(3,3)*Lh3(i,3)
C
      do 5002 j=1,i
         row           = rowb(j)
         khBB(row,col) = khBB(row,col)    + dbLh1t1*Lh1(j,1)
     $                 + dbLh1t2*Lh1(j,2) + dbLh1t3*Lh1(j,3)
     $                 + dbLh2t1*Lh2(j,1) + dbLh2t2*Lh2(j,2)
     $                 + dbLh2t3*Lh2(j,3) + dbLh3t1*Lh3(j,1)
     $                 + dbLh3t2*Lh3(j,2) + dbLh3t3*Lh3(j,3)
         khBB(col,row) = khBB(row,col)
 5002 continue
C
 5001 continue
C
C.....OUTPUT THE MATRIX PRIOR TO ROTATION (FOR DEBUGGING ONLY)
C
*     open(unit=90,file="khBB.m")
*     write(90,*) "khBB=["
*     do 991 i=1,18
* 991 write(90,9) (khBB(i,j),j=1,18)
*     write(90,*) "     ];"
*     close(90)
*   9 format(18(1x,E16.9))
C
C.....ROTATE THE OUTPUT STIFFNESS MATRIX
C
C      call compmrot( khBB , rot , rot , rot )
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
      write(*,*) "*** FATAL ERROR in routine COMPHBB       ***"
      write(*,*) "*** The Stiffness Factor [f] is Negative ***"
      write(*,*) "*** Check the Calling Sequence:          ***"
      write(*,*) "*** Factor [f] Must be Positive or Zero  ***"
      write(*,*) "*** EXECUTION TERNINATED RIGHT HERE      ***"
      stop
C
C.....ERROR-MESSAGE IF THE TRIANGLE'S AREA IS NEGATIVE OR ZERO
C
  200 continue
      write(*,*) "*** FATAL ERROR in routine COMPHBB         ***"
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
C=end of routine "COMPHBB"
C========================C
