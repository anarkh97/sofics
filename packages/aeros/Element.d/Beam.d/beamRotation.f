C
C Authors: Michel Lesoinne and Kendall Pierson
C
C************************************************
C This subroutine rotates the local stiffness
C matrix to global (or computational) axis defined
C in each corner node for a 2 node Beam element
C************************************************
C
      subroutine beamRotation(sm,r1)
C
C Declarations
C
      real*8 r1(3,3), sm(12,12)
C
C Local Declarations
C
      real*8 prod1(12)
      integer i, j, k
C
C premultiplication by rotation matrix
C
      do 10 j=1,12
        do 20 k=1,3
          prod1(k)    = 0.0d+00
          prod1(k+3)  = 0.0d+00
          prod1(k+6)  = 0.0d+00
          prod1(k+9)  = 0.0d+00
          do 20 i=1,3
            prod1(k)    = prod1(k)    + r1(k,i)*sm(i   ,j)
            prod1(k+3)  = prod1(k+3)  + r1(k,i)*sm(i+3 ,j)
            prod1(k+6)  = prod1(k+6)  + r1(k,i)*sm(i+6 ,j)
            prod1(k+9)  = prod1(k+9)  + r1(k,i)*sm(i+9 ,j)
   20   continue

        do 30 k=1,12
          sm(k,j) = prod1(k)
   30   continue
   10 continue
C
C postmultiplication by rotation matrix transposed
C
      do 40 j=1,12
        do 50 k=1,3
          prod1(k)   =0.0d+00
          prod1(k+3) =0.0d+00
          prod1(k+6) =0.0d+00
          prod1(k+9) =0.0d+00
          do 50 i=1,3
            prod1(k)   =prod1(k)   +sm(j,i)*r1(k,i)
            prod1(k+3) =prod1(k+3) +sm(j,i+3)*r1(k,i)
            prod1(k+6) =prod1(k+6) +sm(j,i+6)*r1(k,i)
            prod1(k+9) =prod1(k+9) +sm(j,i+9)*r1(k,i)
   50   continue
C
        do 60 k=1,12
          sm(j,k) = prod1(k)
   60   continue
   40 continue
C
      return
      end
