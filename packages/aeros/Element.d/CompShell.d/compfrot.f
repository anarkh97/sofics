C=====================================================================C
C        1         2         3         4         5         6         7C
C234567890123456789012345678901234567890123456789012345678901234567890C
C=====================================================================C
      subroutine compfrot( v1o , v2o , v3o , v1n , v2n , v3n , rot )
C=====================================================================C
C                                                                     C
C     This Subroutine Forms the Rotation Matrix From Local Elemental  C
C     Frames to Global (Computational) Axes such that the Rotation    C
C     Matrix is:                                                      C
C                                                                     C
C           [ [rot]   0     0   ]                                     C
C     [R] = [   0   [rot]   0   ]   and:   [Ke] = [R] * [ke] * [R]^T  C
C           [   0     0   [rot] ]                                     C
C                                                                     C
C     where [Ke] and [ke] Represent the (Same) Elemental Stiffness    C
C     Matrix Expressed in the Global and Local Frame Systems, Resp.   C
C     (See Routine "compmrot.f".)                                     C
C                                                                     C
C     The Rotation Provides the Transformation Matrix [rot] for:      C
C                                                                     C
C     [v] = [rot] * [x]                                               C
C                                                                     C
C     where:                                                          C
C                                                                     C
C     [x] are the Coordinates in the Old System                       C
C     [v] are the Coordinates in the New System                       C
C                                                                     C
C     The Input [v1o], [v2o] and [v3o] Represent the Old Basis and    C
C     the Input [v1n], [v2n] and [v3n] Represent the New Basis Both   C
C     in Cartesian Coordinates.                                       C
C                                                                     C
C=====================================================================C
C
C     -----------
C     DECLARATION
C     -----------
C
C.....Global Variables
C
      real*8    v1o(3) , v2o(3) , v3o(3)
      real*8    v1n(3) , v2n(3) , v3n(3)
      real*8    rot(6,6)
C
C.....Local Variables
C
      integer   i , j
      real*8    zero
C
C     ----
C     DATA
C     ----
C
      data zero /0.000000D+00/
C
C     -----
C     LOGIC
C     -----
C
C.....CLEAR THE OUTPUT ROTATION MATRIX
C
      do 1001 j=1,6
         do 1002 i=1,6
            rot(i,j) = zero
 1002    continue
 1001 continue
C
C.....COMPUTE THE ROTATION MATRIX FOR TRANSLATIONAL DOFS
C
      do 2001 i=1,3
         rot(1,1) = rot(1,1) + v1n(i)*v1o(i)
         rot(1,2) = rot(1,2) + v1n(i)*v2o(i)
         rot(1,3) = rot(1,3) + v1n(i)*v3o(i)
         rot(2,1) = rot(2,1) + v2n(i)*v1o(i)
         rot(2,2) = rot(2,2) + v2n(i)*v2o(i)
         rot(2,3) = rot(2,3) + v2n(i)*v3o(i)
         rot(3,1) = rot(3,1) + v3n(i)*v1o(i)
         rot(3,2) = rot(3,2) + v3n(i)*v2o(i)
         rot(3,3) = rot(3,3) + v3n(i)*v3o(i)
 2001 continue
C
C.....FILL OUT THE ROTATION MATRIX FOR ROTATIONAL DOFS
C
      do 3001 j=1,3
         do 3002 i=1,3
            rot(i+3,j+3) = rot(i,j)
 3002    continue
 3001 continue
C
C     ------
C     RETURN
C     ------
C
      return
C
C     ---
C     END
C     ---
C
      end
C=end of routine "COMPFROT"
C=========================C
