C=====================================================================C
C        1         2         3         4         5         6         7C
C234567890123456789012345678901234567890123456789012345678901234567890C
C=====================================================================C
      subroutine compcrd( elm  , type , X   , Y    , Z    , rot  ,
     $                    xlp  , ylp  , zlp , rowb , colb , rowm ,
     $                    colm                                   )
C=====================================================================C
C                                                                     C
C     Perform =    This subroutine computes basic quantities needed   C
C     ---------    for the assembly of the basic and higher order     C
C                  composite stiffness matrices.                      C
C                                                                     C
C                                                                     C
C     Inputs/Outputs =                                                C
C     ----------------                                                C
C     ELM   <input>   finite element number                           C
C     TYPE  <input>   type of constitutive law                        C
C     X     <input>   nodal coordinates in the X-direction            C
C     Y     <input>   nodal coordinates in the Y-direction            C
C     Z     <input>   nodal coordinates in the Z-direction            C
C     ROT   <output>  rotation matrix obtained from nodal points      C
C     XLP   <output>  triangular coordinates in the X-direction       C
C     YLP   <output>  triangular coordinates in the Y-direction       C
C     ZLP   <output>  triangular coordinates in the Z-direction       C
C     ROWB  <output>  row pointer for bending DOFs                    C
C     COLB  <output>  column pointer for bending DOFs                 C
C     ROWM  <output>  row pointer for membrane DOFs                   C
C     COLM  <output>  column pointer for membrane DOFs                C
C                                                                     C
C                                                                     C
C     Computations =                                                  C
C     --------------                                                  C
C                                                                     C
C                                                                     C
C     Caution =                                                       C
C     ---------                                                       C
C                                                                     C
C                                                                     C
C     Outputs = no output.                                            C
C     ---------                                                       C
C                                                                     C
C=====================================================================C
C=Author  = Francois M. Hemez                                         C
C=Date    = June 9th, 1994                                            C
C=Version = 1.0                                                       C
C=Comment =                                                           C
C=====================================================================C
C
C     -----------
C     DECLARATION
C     -----------
C
C.....Global Variables
C
      integer   elm , type , rowb(9) , colb(9) , rowm(9) , colm(9)
      real*8    X(3) , Y(3) , Z(3)
      real*8    rot(6,6) , xlp(3) , ylp(3) , zlp(3)
C
C.....Local Variables
C
      integer   i
C     integer   j
      real*8    zero , side21length , side32length , projection
      real*8    x21 , y21 , z21 , x13 , y13 , z13 , x32 , y32 , z32
      real*8    signedarea , area , lengthY , lengthZ , one
      real*8    xp(3) , yp(3) , zp(3) , xcg , ycg , zcg
      real*8    xlcg , ylcg , zlcg , v1n(3) , v2n(3) , v3n(3)
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
C.....CLEAR THE OUTPUT COORDINATES
C
      do 1001 i=1,3
         xlp(i) = zero
 1001 continue
C
      do 1002 i=1,3
         ylp(i) = zero
 1002 continue
C
      do 1003 i=1,3
         zlp(i) = zero
 1003 continue
C
C.....CLEAR THE DEGREE OF FREEDOM POINTERS
C
      do 1006 i=1,9
         rowb(i) = 0
 1006 continue
C
      do 1007 i=1,9
         colb(i) = 0
 1007 continue
C
      do 1008 i=1,9
         rowm(i) = 0
 1008 continue
C
      do 1009 i=1,9
         colm(i) = 0
 1009 continue
C
C.....COMPUTE THE NODAL COORDINATE DIFFERENCES
C
      x21 = X(2) - X(1)
      y21 = Y(2) - Y(1)
      z21 = Z(2) - Z(1)
C
      x13 = X(1) - X(3)
      y13 = Y(1) - Y(3)
      z13 = Z(1) - Z(3)
C
      x32 = X(3) - X(2)
      y32 = Y(3) - Y(2)
      z32 = Z(3) - Z(2)
C
C.....COMPUTE THE LENGTH OF SIDE 2-1
C
      side21length = sqrt(x21*x21 + y21*y21 + z21*z21)
C
C.....CHECK IF LENGTH 2-1 IS DIFFERENT FROM ZERO
C
      if ( side21length.eq.zero ) go to 100
C
C.....COMPUTE THE LENGTH OF SIDE 3-2
C
      side32length = sqrt(x32*x32 + y32*y32 + z32*z32)
C
C.....CHECK IF LENGTH 3-2 IS DIFFERENT FROM ZERO
C
      if ( side32length.eq.zero ) go to 200
C
C.....COMPUTE THE DISTANCE OF THE OPPOSING NODE 3 TO SIDE 2-1
C
      projection = abs(x21*x32 + y21*y32 + z21*z32)/side21length
C
C.....GET THE AREA OF THE TRIANGLE
C
      signedarea = side32length*side32length - projection*projection
C
      if ( signedarea.le.zero ) go to 300
C
      area = 0.50D+00*side21length*sqrt(signedarea)
C
C.....COMPUTE THE DIRECTION COSINES OF THE LOCAL SYSTEM
C.....DIRECTION [X] IS DIRECTED PARALLEL TO THE SIDE 2-1
C.....DIRECTION [Z] IS THE EXTERNAL NORMAL (COUNTERCLOCKWISE)
C.....DIRECTION [Y] IS COMPUTED AS [Z] x [X] (TENSORIAL PRODUCT)
C
      xp(1)   = x21/side21length
      xp(2)   = y21/side21length
      xp(3)   = z21/side21length
C
      zp(1)   = y21*z32 - z21*y32
      zp(2)   = z21*x32 - x21*z32
      zp(3)   = x21*y32 - y21*x32
C
      lengthZ = sqrt(zp(1)*zp(1) + zp(2)*zp(2) + zp(3)*zp(3))
C
      if ( lengthZ.eq.zero ) go to 400
C
      zp(1)   = zp(1)/lengthZ
      zp(2)   = zp(2)/lengthZ
      zp(3)   = zp(3)/lengthZ
C
      yp(1)   = zp(2)*xp(3) - zp(3)*xp(2)
      yp(2)   = zp(3)*xp(1) - zp(1)*xp(3)
      yp(3)   = zp(1)*xp(2) - zp(2)*xp(1)
C
      lengthY = sqrt(yp(1)*yp(1) + yp(2)*yp(2) + yp(3)*yp(3))
C
      if ( lengthY.eq.zero ) go to 400
C
      yp(1)   = yp(1)/lengthY
      yp(2)   = yp(2)/lengthY
      yp(3)   = yp(3)/lengthY
C
C.....COMPUTE THE COORDINATES FOR THE CENTER OF GRAVITY
C
      xcg = (X(1) + X(2) + X(3))/3.00D+00
      ycg = (Y(1) + Y(2) + Y(3))/3.00D+00
      zcg = (Z(1) + Z(2) + Z(3))/3.00D+00
C
C.....COMPUTE THE LOCAL COORDINATES
C
      do 2001 i=1,3
         xlcg   = X(i) - xcg
         ylcg   = Y(i) - ycg
         zlcg   = Z(i) - zcg
         xlp(i) = xp(1)*xlcg + xp(2)*ylcg + xp(3)*zlcg
         ylp(i) = yp(1)*xlcg + yp(2)*ylcg + yp(3)*zlcg
         zlp(i) = zp(1)*xlcg + zp(2)*ylcg + zp(3)*zlcg
 2001 continue
C
C.....COMPUTE THE NODAL ROTATION MATRIX
C
      v1n(1)=  one
      v1n(2)= zero
      v1n(3)= zero
C
      v2n(1)= zero
      v2n(2)=  one
      v2n(3)= zero
C
      v3n(1)= zero
      v3n(2)= zero
      v3n(3)=  one
C
      call compfrot( xp , yp , zp , v1n , v2n , v3n , rot )
C
C.....DEFINE ROW POINTER FOR BENDING CONTRIBUTIONS
C
      rowb(1) =  3
      rowb(2) =  4
      rowb(3) =  5
      rowb(4) =  9
      rowb(5) = 10
      rowb(6) = 11
      rowb(7) = 15
      rowb(8) = 16
      rowb(9) = 17
C
C.....DEFINE COLUMN POINTER FOR BENDING CONTRIBUTIONS
C
      colb(1) =  3
      colb(2) =  4
      colb(3) =  5
      colb(4) =  9
      colb(5) = 10
      colb(6) = 11
      colb(7) = 15
      colb(8) = 16
      colb(9) = 17
C
C.....DEFINE ROW POINTER FOR MEMBRANE CONTRIBUTION
C
      rowm(1) =  1
      rowm(2) =  2
      rowm(3) =  7
      rowm(4) =  8
      rowm(5) = 13
      rowm(6) = 14
      rowm(7) =  6
      rowm(8) = 12
      rowm(9) = 18
C
C.....DEFINE COLUMN POINTER FOR MEMBRANE CONTRIBUTION
C
      colm(1) =  1
      colm(2) =  2
      colm(3) =  7
      colm(4) =  8
      colm(5) = 13
      colm(6) = 14
      colm(7) =  6
      colm(8) = 12
      colm(9) = 18
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
C.....ERROR-MESSAGE IF SIDE 2-1 HAS ZERO-LENGTH
C
  100 continue
      write(*,*) "*** FATAL ERROR in routine COMPCRD          ***"
      write(*,*) "*** Side Between Nodes 1 and 2 Has 0-Length ***"
      write(*,*) "*** Check Coordinates and FE Topology       ***"
      write(*,*) "*** EXECUTION TERMINATED HERE               ***"
      stop
C
C.....ERROR-MESSAGE IF SIDE 3-2 HAS ZERO-LENGTH
C
  200 continue
      write(*,*) "*** FATAL ERROR in routine COMPCRD          ***"
      write(*,*) "*** Side Between Nodes 2 and 3 Has 0-Length ***"
      write(*,*) "*** Check Coordinates and FE Topology       ***"
      write(*,*) "*** EXECUTION TERMINATED HERE               ***"
      stop
C
C.....ERROR-MESSAGE IF THE AREA IS NEGATIVE OR ZERO
C
  300 continue
      write(*,*) "*** FATAL ERROR in routine COMPCRD    ***"
      write(*,*) "*** The Area is Negative or Zero      ***"
      write(*,*) "*** Check Coordinates and FE Topology ***"
      write(*,*) "*** EXECUTION TERMINATED HERE         ***"
      stop
C
C.....ERROR-MESSAGE IF A LOCAL FRAME VECTOR HAS ZERO-LENGTH
C
  400 continue
      write(*,*) "*** FATAL ERROR in routine COMPCRD     ***"
      write(*,*) "*** A Local Frame Vector has 0-Length  ***"
      write(*,*) "*** Check FE Topology and Computations ***"
      write(*,*) "*** EXECUTION TERMINATED HERE          ***"
      stop
C
C     ---
C     END
C     ---
C
      end
C=end of routine "COMPCRD"
C========================C
