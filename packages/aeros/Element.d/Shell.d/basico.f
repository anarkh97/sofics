c=deck basico basico fortran
c=purpose form basico bending stiffness of c1 triangle
c=author c. a. felippa, may 1984
c=version september 1989
c=equipment machine independent
c=keywords thin plate bending
c=keywords finite element triangle basico stiffness matrix
c=block abstract
c
c     basico forms the basico material element stiffness matrix
c     of a c1 triangle constructed with the generalized
c     ff/ans formulation
c
c=end abstract
c=block usage
c
c     the calling sequence is
c
c       call      basico (x,y,db,f,clr,cqr,ls,sm,m,status)
c
c     where the input arguments are
c
c       x         (3 x 1) array of x coordinates of triangle nodes
c       y         (3 x 1) array of y coordinates of triangle nodes
c       db        (5 x 5) stress resultant-strain matrix
c       f         factor by which stiffness entries will be multiplied.
c       clr,cqr   use clr*llr+cqr*lqr for l
c                 thus clr+cqr  must add up to 1.
c                 llr=linear rotation, lqr=quadratic rotation (kirchhoff)
c       ls        (9 x 1) array of stiffness location pointers
c                 (see output sm).
c       sm        incoming material stiffness array.
c       m         first dimension of sm in calling program.
c
c     the outputs are:
c
c       sm        output stiffness array with basico stiffness
c                 coefficients added in.  the (i,j)-th entry of the
c                 (9 by 9) element bending stiffness is added to
c                 sm(k,l), where k=ls(i) and l=ls(j).
c       status    status character variable.  blank if no error
c                 detected.
c
c=end usage
c=block fortran
      subroutine    basico( x, y, db, f, clr,cqr,ls,sm, m, status)
c
c                   a r g u m e n t s
c
      real*8             x(3), y(3), db(3,3), f, sm(m,m)
      real*8             clr, cqr
      integer            m, ls(9)
      character*(*)      status
c
c                   t y p e   &   d i m e n s i o n
c
      real*8             llr(9,3), lqr(9,3), l(9,3)
      real*8             db11, db12, db13, db22, db23, db33
      real*8             x0, y0, cab, a1, a2, a3, b1, b2, b3
      real*8             x21, x32, x13, y21, y32, y13
      real*8             x12, x23, x31, y12, y23, y31
      real*8             xl12, xl23, xl31, c12, c23, c31, s12, s23, s31
      real*8             cc12, cc23, cc31, ss12, ss23, ss31
      real*8             cs12, cs23, cs31, s1, s2, s3
      real*8             area2, c
      integer            i, j, ii, jj
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
C
      xl12 =    dsqrt(x12*x12+y12*y12)
      xl23 =    dsqrt(x23*x23+y23*y23)
      xl31 =    dsqrt(x31*x31+y31*y31)
C
      do 1200  j = 1,3
        do 1200  i = 1,9
          llr(i,j) =  0.0
          lqr(i,j) =  0.0
 1200   continue
c
      if (clr .ne. 0.0)         then
        llr(3,1) =   y32*.5
        llr(6,1) =   y13*.5
        llr(9,1) =   y21*.5
        llr(2,2) =   x32*.5
        llr(5,2) =   x13*.5
        llr(8,2) =   x21*.5
        llr(2,3) =  -y32*.5
        llr(3,3) =  -x32*.5
        llr(5,3) =  -y13*.5
        llr(6,3) =  -x13*.5
        llr(8,3) =  -y21*.5
        llr(9,3) =  -x21*.5
      end if
c
      if (cqr .ne. 0.0)         then
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
      end if

c     write(*,*) ((lqr(i,j),j=1,3),i=1,9)
c
      do 1600  j = 1,9
        l(j,1) =   clr*llr(j,1) + cqr*lqr(j,1)
        l(j,2) =   clr*llr(j,2) + cqr*lqr(j,2)
        l(j,3) =   clr*llr(j,3) + cqr*lqr(j,3)
 1600   continue
c
      c =        2.0*f/area2
      db11 =     c*db(1,1)
      db22 =     c*db(2,2)
      db33 =     c*db(3,3)
      db12 =     c*db(1,2)
      db13 =     c*db(1,3)
      db23 =     c*db(2,3)
      do 4000  j = 1,9
        jj =     ls(j)
        s1 =     db11*l(j,1) + db12*l(j,2) + db13*l(j,3)
        s2 =     db12*l(j,1) + db22*l(j,2) + db23*l(j,3)
        s3 =     db13*l(j,1) + db23*l(j,2) + db33*l(j,3)
        do 3500  i = 1,j
          ii =      ls(i)
          sm(jj,ii) =  sm(jj,ii) + (s1*l(i,1)+s2*l(i,2)+s3*l(i,3))
          sm(ii,jj) =  sm(jj,ii)
 3500     continue
 4000   continue
      return
      end
c=end fortran
