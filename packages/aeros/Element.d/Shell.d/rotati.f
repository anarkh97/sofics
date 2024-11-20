C this routine provides the transformation
C matrix  R for:
C
C         vp= R x
C
C  x are the coordinates in the old system
C  vp are the coordinates in the new system
C
C  Input:
C       v1o,v2o,v3o old basis described in cartesian
C                   components
C       v1n,v2n,v3n new basis
C
      subroutine rotati(v1o,v2o,v3o,v1n,v2n,v3n,r)
C
      real*8    v1o(3),v2o(3),v3o(3)
      real*8    v1n(3),v2n(3),v3n(3),r(6,6)
      integer i,j
C
      do 5  i=1,6
        do 5 j=1,6
   5      r(i,j)=0.0d+00

      do 10 i=1,3
        r(1,1) = r(1,1)+v1n(i)*v1o(i)
        r(1,2) = r(1,2)+v1n(i)*v2o(i)
        r(1,3) = r(1,3)+v1n(i)*v3o(i)
        r(2,1) = r(2,1)+v2n(i)*v1o(i)
        r(2,2) = r(2,2)+v2n(i)*v2o(i)
        r(2,3) = r(2,3)+v2n(i)*v3o(i)
        r(3,1) = r(3,1)+v3n(i)*v1o(i)
        r(3,2) = r(3,2)+v3n(i)*v2o(i)
        r(3,3) = r(3,3)+v3n(i)*v3o(i)
   10 continue
      do 20 i=1,3
        do 20 j=1,3
          r(i+3,j+3)= r(i,j)
   20 continue
      return
      end
