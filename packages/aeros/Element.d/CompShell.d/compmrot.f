C=====================================================================C
C        1         2         3         4         5         6         7C
C234567890123456789012345678901234567890123456789012345678901234567890C
C=====================================================================C
      subroutine compmrot( ke , r1 , r2 , r3 )
C=====================================================================C
C                                                                     C
C     This Subroutine Rotates the Local Stiffness Matrix to the       C
C     Global (Computational) Axes Defined in Each Corner Node:        C
C                                                                     C
C     [ke] Becomes [R] * [ke] * [R]^T where:                          C
C                                                                     C
C           [ [r1]  0    0   ]                                        C
C     [R] = [  0   [r2]  0   ]                                        C
C           [  0    0   [r3] ]                                        C
C                                                                     C
C=====================================================================C
C
C     -----------
C     DECLARATION
C     -----------
C
C.....Global Variables
C
      real*8    ke(18,18) , r1(6,6) , r2(6,6) , r3(6,6)
C
C.....Local Variables
C
      integer   i , j
C
      real*8    rke(18,18) , zero
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
      do 1001 i=1,18
      do 1001 j=1,18
         rke(i,j) = zero
 1001 continue
C
      do i=1,6
        do j=1,18
          do k=1,6
            rke(i,j)    = rke(i,j)    + r1(i,k)*ke(k,j)
            rke(6+i, j) = rke(6+i,j)  + r2(i,k)*ke(6+k,j)
            rke(12+i,j) = rke(12+i,j) + r3(i,k)*ke(12+k,j)
          end do
        end do
      end do
C
      do i=1,18
        do j=1,6
          ke(i,j)    = 0.0
	  ke(i,j+6)  = 0.0
          ke(i,j+12) = 0.0
          do k=1,6
            ke(i,j)    = ke(i,j)    + rke(i,k)   * r1(j,k)
            ke(i,j+6)  = ke(i,j+6)  + rke(i,k+6) * r2(j,k)
            ke(i,j+12) = ke(i,j+12) + rke(i,k+12)* r3(j,k)
          end do
        end do
      end do
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
C=end of routine "COMPMROT"
C=========================C
