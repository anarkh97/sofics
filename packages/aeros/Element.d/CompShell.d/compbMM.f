C=====================================================================C
C        1         2         3         4         5         6         7C
C234567890123456789012345678901234567890123456789012345678901234567890C
C=====================================================================C
      subroutine compbMM( elm , type , x    , y   , dm , alpha ,
     $                    f   , rowm , colm , rot , P  , kbMM  )
C=====================================================================C
C
C     -----------
C     DECLARATION
C     -----------
C
C.....Global Variables
C
      integer   elm , type , rowm(9) , colm(9)
      real*8    x(3) , y(3) , dm(3,3) , alpha , f
      real*8    kbMM(18,18) , rot(6,6) , P(9,3)
C
C.....Local Variables
C
      integer   i , j , row , col , dimp
      real*8    zero , three , six
      real*8    twicearea , factor
      real*8    x21, x32, x13, y21, y32, y13
      real*8    x12, x23, x31, y12, y23, y31
      real*8    dmPt1 , dmPt2 , dmPt3
C
C     ----
C     DATA
C     ----
C
      data zero  /0.000000D+00/
      data three /3.000000D+00/
      data six   /6.000000D+00/
C
C     -----
C     LOGIC
C     -----
C
C.....CLEAR THE OUTPUT STIFFNESS MATRIX
C
      do 1001 j=1,18
         do 1002 i=1,18
            kbMM(i,j) = zero
 1002    continue
 1001 continue
C
C.....CLEAR THE LOCAL MATRICES
C
      do 1003 j=1,3
         do 1004 i=1,9
            P(i,j) = zero
 1004    continue
 1003 continue
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
C.....CHECK THE AREA (ERROR IF NOT POSITIVE)
C
      if ( twicearea.le.zero ) go to 200
C
C.....ASSEMBLE THE MATRIX [P] W/ SHAPE FUNCTION DERIVATIVES
C
      P(1,1) = y23
      P(2,1) = zero
      P(3,1) = y31
      P(4,1) = zero
      P(5,1) = y12
      P(6,1) = zero
C
      P(1,2) = zero
      P(2,2) = x32
      P(3,2) = zero
      P(4,2) = x13
      P(5,2) = zero
      P(6,2) = x21
C
      P(1,3) = x32
      P(2,3) = y23
      P(3,3) = x13
      P(4,3) = y31
      P(5,3) = x21
      P(6,3) = y12
C
      dimp   = 6
C
      if ( alpha.ne.zero ) then
C
         P(7,1) = y23*(y13-y21)*alpha/six
         P(7,2) = x32*(x31-x12)*alpha/six
         P(7,3) = (x31*y13-x12*y21)*alpha/three
C
         P(8,1) = y31*(y21-y32)*alpha/six
         P(8,2) = x13*(x12-x23)*alpha/six
         P(8,3) = (x12*y21-x23*y32)*alpha/three
C
         P(9,1) = y12*(y32-y13)*alpha/six
         P(9,2) = x21*(x23-x31)*alpha/six
         P(9,3) = (x23*y32-x31*y13)*alpha/three
C
         dimp   = 9
C
      endif
C
C.....MULTIPLY MATRIX [P] BY THE SQUARE ROOT OF THE STIFFNESS FACTOR
C.....AND DIVIDE BY THE SQUARE ROOT OF FOUR TIMES THE AREA
C
      factor = sqrt(f/(2.00D+00*twicearea))
C
      do 2001 j=1,3
         do 2002 i=1,dimp
            P(i,j) = factor*P(i,j)
 2002    continue
 2001 continue
C
C.....ASSEMBLE THE OUTPUT STIFFNESS SUCH THAT [kbMM] = [P]*[dm]*[P]^T
C
      do 3001 i=1,dimp
C
         col   = colm(i)
         dmPt1 = dm(1,1)*P(i,1) + dm(1,2)*P(i,2) + dm(1,3)*P(i,3)
         dmPt2 = dm(2,1)*P(i,1) + dm(2,2)*P(i,2) + dm(2,3)*P(i,3)
         dmPt3 = dm(3,1)*P(i,1) + dm(3,2)*P(i,2) + dm(3,3)*P(i,3)
C
         do 3002 j=1,i
            row           = rowm(j)
            kbMM(row,col) = kbMM(row,col)
     $                    + dmPt1*P(j,1)
     $                    + dmPt2*P(j,2)
     $                    + dmPt3*P(j,3)
            kbMM(col,row) = kbMM(row,col)
 3002    continue
C
 3001 continue
C
C.....OUTPUT THE MATRIX PRIOR TO ROTATION (FOR DEBUGING ONLY)
C
*     open(unit=90,file="kbMM.m")
*     write(90,*) "kbMM=["
*     do 991 i=1,18
* 991 write(90,9) (kbMM(i,j),j=1,18)
*     write(90,*) "     ];"
*     close(90)
*   9 format(18(1x,E16.9))
C
C.....ROTATE THE OUTPUT STIFFNESS MATRIX
C
C      call compmrot( kbMM , rot , rot , rot )
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
      write(*,*) "*** FATAL ERROR in routine COMPBMM       ***"
      write(*,*) "*** The Stiffness Factor [f] is Negative ***"
      write(*,*) "*** Check the Calling Sequence:          ***"
      write(*,*) "*** Factor [f] Must be Positive or Zero  ***"
      write(*,*) "*** EXECUTION TERNINATED RIGHT HERE      ***"
      stop
C
C.....ERROR-MESSAGE IF THE TRIANGLE'S AREA IS NEGATIVE OR ZERO
C
  200 continue
      write(*,*) "*** FATAL ERROR in routine COMPBMM         ***"
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
C=end of routine "COMPBMM"
C========================C
