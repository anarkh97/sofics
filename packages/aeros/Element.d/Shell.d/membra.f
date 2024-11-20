C=BLOCK FORTRAN
      subroutine    membra
     $             (x, y, alpha, le, q, status)
C
C                   T Y P E   &   D I M E N S I O N
C
C=BLOCK VAX
C     implicit      none
C=END VAX
      character*(*)  status
C=BLOCK REAL*8
      real*8  x(3), y(3), p(9,3), q(18,3)
      real*8  area2, c , alpha
      real*8  x21, x32, x13, y21, y32, y13
      real*8  x12, x23, x31, y12, y23, y31
C=ELSE
C=END REAL*8
      integer       i, j, kk
      integer       le(9), n
C
C                   L O G I C
C
      status =   ' '
      x21 =      x(2) - x(1)
      x12 =     -x21
      x32 =      x(3) - x(2)
      x23 =     -x32
      x13 =      x(1) - x(3)
      x31 =     -x13
      y21 =      y(2) - y(1)
      y12 =     -y21
      y32 =      y(3) - y(2)
      y23 =     -y32
      y13 =      y(1) - y(3)
      y31 =     -y13
      area2 =    y21*x13 - x21*y13
      if (area2 .le. 0.0)      then
        status = 'SM3MB: Zero area'
        if (area2 .eq. 0.0)   status = 'SM3MB: Zero area'
        return
      end if
      p(1,1) =   y23
      p(2,1) =   0.0
      p(3,1) =   y31
      p(4,1) =   0.0
      p(5,1) =   y12
      p(6,1) =   0.0
      p(1,2) =   0.0
      p(2,2) =   x32
      p(3,2) =   0.0
      p(4,2) =   x13
      p(5,2) =   0.0
      p(6,2) =   x21
      p(1,3) =   x32
      p(2,3) =   y23
      p(3,3) =   x13
      p(4,3) =   y31
      p(5,3) =   x21
      p(6,3) =   y12
      if (alpha .ne. 0.0)         then
        p(7,1) =  y23*(y13-y21)*alpha/6.
        p(7,2) =  x32*(x31-x12)*alpha/6.
        p(7,3) =  (x31*y13-x12*y21)*alpha/3.
        p(8,1) =  y31*(y21-y32)*alpha/6.
        p(8,2) =  x13*(x12-x23)*alpha/6.
        p(8,3) =  (x12*y21-x23*y32)*alpha/3.
        p(9,1) =  y12*(y32-y13)*alpha/6.
        p(9,2) =  x21*(x23-x31)*alpha/6.
        p(9,3) =  (x23*y32-x31*y13)*alpha/3.
        n =        9
      end if
      c =   1.0   /area2
      do 300 i=1,9
      kk=le(i)
        do 300 j=1,3
          q(kk,j)=p(i,j)*c
300   continue
      return
      end
C=END FORTRAN
