C=====================================================================C
C        1         2         3         4         5         6         7C
C234567890123456789012345678901234567890123456789012345678901234567890C
C=====================================================================C
      subroutine compbBB( elm  , type , x    , y    , db  , f ,
     $                    clr  , cqr  , rowb , colb , rot , L ,
     $                    kbBB                                )
C=====================================================================C
C
C     -----------
C     DECLARATION
C     -----------
C
C.....Global Variables
C
      integer   elm , type , rowb(9) , colb(9)
      real*8    x(3) , y(3) , db(3,3) , kbBB(18,18) , rot(6,6)
      real*8    f , clr , cqr , L(9,3)
C
C.....Local Variables
C
      integer   i , j , row , col
      real*8    zero , twicearea , llr(9,3) , lqr(9,3)
      real*8    x21 , x32 , x13 , y21 , y32 , y13
      real*8    x12 , x23 , x31 , y12 , y23 , y31
      real*8    x0 , y0 , dist12 , dist23 , dist31
      real*8    c12 , c23 , c31 , s12 , s23 , s31
      real*8    cc12 , cc23 , cc31 , ss12 , ss23 , ss31
      real*8    one , cs12 , cs23 , cs31 , area
      real*8    dbLt1 , dbLt2 , dbLt3 , factor
C
C     ----
C     DATA
C     ----
C
      data zero /0.000000D+00/
      data one  /1.000000D+00/
C
C     -----
C     LOGIC
C     -----
C
C.....CLEAR THE OUTPUT STIFFNESS MATRIX
C
      do 1001 j=1,18
         do 1002 i=1,18
            kbBB(i,j) = zero
 1002    continue
 1001 continue
C
C.....CLEAR THE LOCAL MATRICES
C
      do 1003 j=1,3
         do 1004 i=1,9
            llr(i,j) = zero
 1004    continue
         do 1005 i=1,9
            lqr(i,j) = zero
 1005    continue
         do 1006 i=1,9
            L(i,j) = zero
 1006    continue
 1003 continue
C
C.....RETURN IF THE STIFFNESS FACTOR IS ZERO
C
      if ( f.eq.zero ) return
C
C.....CHECK IF [CLR] AND [CQR] SATISFY THE CONSTRAINT [CLR]+[CQR]=1
C
      if ( (clr+cqr).ne.one ) go to 100
C
C.....CHECK IF THE STIFFNESS FACTOR [F] IS POSITIVE
C
      if ( f.lt.zero ) go to 200
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
      if ( twicearea.le.zero ) go to 300
C
C.....GET THE COORDINATES OF THE CENTROID OF THE TRIANGLE
C
      x0 = ( x(1) + x(2) + x(3) )/3.00D+00
      y0 = ( y(1) + y(2) + y(3) )/3.00D+00
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
         llr(3,1) =  0.50D+00*y32
         llr(6,1) =  0.50D+00*y13
         llr(9,1) =  0.50D+00*y21
C
         llr(2,2) =  0.50D+00*x32
         llr(5,2) =  0.50D+00*x13
         llr(8,2) =  0.50D+00*x21
C
         llr(2,3) = -0.50D+00*y32
         llr(3,3) = -0.50D+00*x32
         llr(5,3) = -0.50D+00*y13
         llr(6,3) = -0.50D+00*x13
         llr(8,3) = -0.50D+00*y21
         llr(9,3) = -0.50D+00*x21
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
         lqr(2,1) =  0.50D+00*( cc12*x12 + cc31*x31 )
         lqr(2,2) =  0.50D+00*( ss12*x12 + ss31*x31 )
         lqr(2,3) =  ss12*y21 + ss31*y13
C
         lqr(3,1) = -0.50D+00*( cc12*y21 + cc31*y13 )
         lqr(3,2) = -0.50D+00*lqr(2,3)
         lqr(3,3) = -2.00D+00*lqr(2,1)
C
         lqr(4,1) =  cs23 - cs12
         lqr(4,2) = -lqr(4,1)
         lqr(4,3) =  (cc12-ss12) - (cc23-ss23)
C
         lqr(5,1) =  0.50D+00*( cc12*x12 + cc23*x23 )
         lqr(5,2) =  0.50D+00*( ss12*x12 + ss23*x23 )
         lqr(5,3) =  ss12*y21 + ss23*y32
C
         lqr(6,1) = -0.50D+00*( cc12*y21 + cc23*y32 )
         lqr(6,2) = -0.50D+00*lqr(5,3)
         lqr(6,3) = -2.00D+00*lqr(5,1)
C
         lqr(7,1) =  cs31 - cs23
         lqr(7,2) = -lqr(7,1)
         lqr(7,3) =  (cc23-ss23) - (cc31-ss31)
C
         lqr(8,1) =  0.50D+00*( cc23*x23 + cc31*x31 )
         lqr(8,2) =  0.50D+00*( ss23*x23 + ss31*x31 )
         lqr(8,3) =  ss23*y32 + ss31*y13
C
         lqr(9,1) = -0.50D+00*( cc23*y32 + cc31*y13 )
         lqr(9,2) = -0.50D+00*lqr(8,3)
         lqr(9,3) = -2.00D+00*lqr(8,1)
C
      endif
C
C.....ASSEMBLE THE LOCAL MATRIX [L] AS [CLR]*[LLR] + [CQR]*[LQR]
C
      if ( clr.eq.zero ) then
         do 2001 j=1,3
            do 2002 i=1,9
               L(i,j) = lqr(i,j)
 2002       continue
 2001    continue
      endif
C
      if ( cqr.eq.zero ) then
         do 2003 j=1,3
            do 2004 i=1,9
               L(i,j) = llr(i,j)
 2004       continue
 2003    continue
      endif
C
      if ( (clr.ne.zero).and.(cqr.ne.zero) ) then
         do 2005 j=1,3
            do 2006 i=1,9
               L(i,j) = clr*llr(i,j) + cqr*lqr(i,j)
 2006       continue
 2005    continue
      endif
C
C.....MULTIPLY MATRIX [L] BY THE SQUARE ROOT OF THE STIFFNESS FACTOR
C.....AND DIVIDE BY THE SQUARE ROOT OF THE AREA
C
      factor = sqrt(f/area)
C
      do 3001 j=1,3
         do 3002 i=1,9
            L(i,j) = factor*L(i,j)
 3002    continue
 3001 continue
C
C.....ASSEMBLE THE OUTPUT STIFFNESS SUCH THAT [kbBB] = [L]*[db]*[L]^T
C
      do 4001 i=1,9
C
         col   = colb(i)
         dbLt1 = db(1,1)*L(i,1) + db(1,2)*L(i,2) + db(1,3)*L(i,3)
         dbLt2 = db(2,1)*L(i,1) + db(2,2)*L(i,2) + db(2,3)*L(i,3)
         dbLt3 = db(3,1)*L(i,1) + db(3,2)*L(i,2) + db(3,3)*L(i,3)
C
         do 4002 j=1,i
            row           = rowb(j)
            kbBB(row,col) = kbBB(row,col)
     $                    + dbLt1*L(j,1)
     $                    + dbLt2*L(j,2)
     $                    + dbLt3*L(j,3)
            kbBB(col,row) = kbBB(row,col)
 4002    continue
C
 4001 continue
C
C.....OUTPUT THE MATRIX PRIOR TO ROTATION (FOR DEBUGGING ONLY)
C
*     open(unit=90,file="kbBB.m")
*     write(90,*) "kbBB=["
*     do 991 i=1,18
* 991 write(90,9) (kbBB(i,j),j=1,18)
*     write(90,*) "     ];"
*     close(90)
*   9 format(18(1x,E16.9))
C
C.....ROTATE THE OUTPUT STIFFNESS MATRIX
C
C      call compmrot( kbBB , rot , rot , rot )
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
C.....ERROR-MESSAGE IF [CLR]+[CQR] IS DIFFERENT FROM ONE
C
  100 continue
      write(*,*) "*** FATAL ERROR in routine COMPBBB      ***"
      write(*,*) "*** The Factors [clr] and [cqr] Violate ***"
      write(*,*) "*** the Constraint [clr]+[cqr]=1:       ***"
      write(*,*) "*** Check the Calling Sequence          ***"
      write(*,*) "*** EXECUTION TERNINATED RIGHT HERE     ***"
      stop
C
C.....ERROR-MESSAGE IF THE STIFFNESS FACTOR [F] IS NEGATIVE
C
  200 continue
      write(*,*) "*** FATAL ERROR in routine COMPBBB       ***"
      write(*,*) "*** The Stiffness Factor [f] is Negative ***"
      write(*,*) "*** Check the Calling Sequence:          ***"
      write(*,*) "*** Factor [f] Must be Positive or Zero  ***"
      write(*,*) "*** EXECUTION TERNINATED RIGHT HERE      ***"
      stop
C
C.....ERROR-MESSAGE IF THE TRIANGLE'S AREA IS NEGATIVE OR ZERO
C
  300 continue
      write(*,*) "*** FATAL ERROR in routine COMPBBB         ***"
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
C=end of routine "COMPBBB"
C========================C
