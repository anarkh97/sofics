C
C Authors: Michel Lesoinne and Kendall Pierson
C
C Modified routine: 3-31-98
C
C This routine provides the transformation
C matrix  R for:
C
C         vp= R x
C
C  x  are the coordinates in the old system
C  vp are the coordinates in the new system
C
C  Input:
C       v1o,v2o,v3o old basis described in cartesian
C                   components
C       v1n,v2n,v3n new basis
C
      subroutine rotation(v1o,v2o,v3o,v1n,v2n,v3n,r)
C
C Declarations
C
      real*8    v1o(3),v2o(3),v3o(3)
      real*8    v1n(3),v2n(3),v3n(3),r(3,3)
C
C Local Declarations
C
      integer i,j
C
C Zero rotation matrix
C
      do i=1,3
        do j=1,3
          r(i,j)=0.0d+00
        enddo
      enddo
C
C Compute rotation matrix
C
      do i=1,3
        r(1,1) = r(1,1) + v1n(i)*v1o(i)
        r(1,2) = r(1,2) + v1n(i)*v2o(i)
        r(1,3) = r(1,3) + v1n(i)*v3o(i)
        r(2,1) = r(2,1) + v2n(i)*v1o(i)
        r(2,2) = r(2,2) + v2n(i)*v2o(i)
        r(2,3) = r(2,3) + v2n(i)*v3o(i)
        r(3,1) = r(3,1) + v3n(i)*v1o(i)
        r(3,2) = r(3,2) + v3n(i)*v2o(i)
        r(3,3) = r(3,3) + v3n(i)*v3o(i)
      enddo 

      return
      end
