c=block fortran
      subroutine    momen(x,y,lb,l,status)
c
c                   a r g u m e n t s
c
      integer  lb(*)
      real*8   x(*),y(*),l(18,*)
      character*(*)      status
c
c                   t y p e   &   d i m e n s i o n
c
      real*8  lqr(9,3)
      real*8  x0, y0, cab, a1, a2, a3, b1, b2, b3
      real*8  x21, x32, x13, y21, y32, y13
      real*8  x12, x23, x31, y12, y23, y31
      real*8  xl12, xl23, xl31, c12, c23, c31, s12, s23, s31
      real*8  cc12, cc23, cc31, ss12, ss23, ss31
      real*8  cs12, cs23, cs31
      real*8  area2
      integer            i, j,  kk
c
c                   l o g i c
c
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
        status = 'basico: negative area'
        if (area2 .eq. 0.0)   status = 'basico: zero area'
        return
      end if
      x0 =      (x(1)+x(2)+x(3))/3.
      y0 =      (y(1)+y(2)+y(3))/3.
      cab =      3.0/area2
      a1 =      -cab*(y(3)-y0)
      a2 =      -cab*(y(1)-y0)
      a3 =      -cab*(y(2)-y0)
      b1 =       cab*(x(3)-x0)
      b2 =       cab*(x(1)-x0)
      b3 =       cab*(x(2)-x0)
      xl12 =    dsqrt(x12**2+y12**2)
      xl23 =    dsqrt(x23**2+y23**2)
      xl31 =    dsqrt(x31**2+y31**2)
      do 1200  j = 1,3
        do 1200  i = 1,9
          lqr(i,j) =  0.0
 1200   continue
c
c
        c12 =       y21/xl12
        s12 =       x12/xl12
        c23 =       y32/xl23
        s23 =       x23/xl23
        c31 =       y13/xl31
        s31 =       x31/xl31
        cc12 =      c12*c12
        cc23 =      c23*c23
        cc31 =      c31*c31
        ss12 =      s12*s12
        ss23 =      s23*s23
        ss31 =      s31*s31
        cs12 =      c12*s12
        cs23 =      c23*s23
        cs31 =      c31*s31
        lqr(1,1) =   cs12 - cs31
        lqr(1,2) =  -lqr(1,1)
        lqr(1,3) =  (cc31-ss31) - (cc12-ss12)
        lqr(2,1) =  (cc12*x12 + cc31*x31)*.5
        lqr(2,2) =  (ss12*x12 + ss31*x31)*.5
        lqr(2,3) =   ss12*y21 + ss31*y13
        lqr(3,1) = -(cc12*y21 + cc31*y13)*.5
        lqr(3,2) = -0.5*lqr(2,3)
        lqr(3,3) =  -2.*lqr(2,1)
        lqr(4,1) =  cs23 - cs12 
        lqr(4,2) =  -lqr(4,1)
        lqr(4,3) =  (cc12-ss12) - (cc23-ss23)
        lqr(5,1) =  (cc12*x12 + cc23*x23)*.5
        lqr(5,2) =  (ss12*x12 + ss23*x23)*.5
        lqr(5,3) =   ss12*y21 + ss23*y32
        lqr(6,1) = -(cc12*y21 + cc23*y32)*.5
        lqr(6,2) = -0.5*lqr(5,3)
        lqr(6,3) =  -2.*lqr(5,1)
        lqr(7,1) =  cs31 - cs23
        lqr(7,2) =  -lqr(7,1)
        lqr(7,3) =  (cc23-ss23) - (cc31-ss31)
        lqr(8,1) =  (cc23*x23 + cc31*x31)*.5
        lqr(8,2) =  (ss23*x23 + ss31*x31)*.5
        lqr(8,3) =   ss23*y32 + ss31*y13
        lqr(9,1) = -(cc23*y32 + cc31*y13)*.5
        lqr(9,2) = -0.5*lqr(8,3)
        lqr(9,3) =  -2.*lqr(8,1)
c
      do 1600  j = 1,9
        kk= lb(j)
        l(kk,1) = lqr(j,1)*dsqrt(2.00D+00/area2)
        l(kk,2) = lqr(j,2)*dsqrt(2.00D+00/area2)
        l(kk,3) = lqr(j,3)*dsqrt(2.00D+00/area2)
 1600   continue
      return
      end
c=end fortran
