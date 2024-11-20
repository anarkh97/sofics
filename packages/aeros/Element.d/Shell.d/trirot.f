C************************************************
C this subroutine rotates the local stiffness
C matrix to global (or computational) axis defined
C in each corner node
C************************************************

      subroutine trirot(sm,r1,r2,r3)
      real*8 r1(6,6),r2(6,6),r3(6,6)
      real*8 sm(18,18), prod1(18)
      integer i, j, k
C premultiplication by rotation matrix
      do 10 j=1,18
        do 20 k=1,6
          prod1(k)   =0.0d+00
          prod1(k+6) =0.0d+00
          prod1(k+12)=0.0d+00
          do 20 i=1,6
            prod1(k)   =prod1(k)   +r1(k,i)*sm(i,j)
            prod1(k+6) =prod1(k+6) +r2(k,i)*sm(i+6,j)
            prod1(k+12)=prod1(k+12)+r3(k,i)*sm(i+12,j)
   20   continue
        do 30 k=1,18
          sm(k,j)=prod1(k)
   30   continue
   10 continue
C postmultiplication by rotation matrix transposed
      do 40 j=1,18
        do 50 k=1,6
          prod1(k)   =0.0d+00
          prod1(k+6) =0.0d+00
          prod1(k+12)=0.0d+00
          do 50 i=1,6
            prod1(k)   =prod1(k)   +sm(j,i)*r1(k,i)
            prod1(k+6) =prod1(k+6) +sm(j,i+6)*r2(k,i)
            prod1(k+12)=prod1(k+12)+sm(j,i+12)*r3(k,i)
   50   continue
        do 60 k=1,18
          sm(j,k)=prod1(k)
   60   continue
   40 continue
C final
      return
      end
