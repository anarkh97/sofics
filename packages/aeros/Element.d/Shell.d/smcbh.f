c=deck smcbh smcbh fortran
c=purpose form higher-order bending stiffness obtained
c          from curvatures over the sides
c=author c. militello, may 1989
c=version may 1989
c=equipment machine independent
c=keywords thin plate bending
c=keywords finite element triangle higher order stiffness matrix
c=block abstract
c
c     smcbh forms the higher order material stiffness matrix of a
c     9-dof thin-plate-bending triangle obtained by using linear
c     curvatures over the sides
c
c=end abstract
c=block usage
c
c     the calling sequence is
c
c       call      smcbh (x, y, db, f, ,ls, sm, m, status)
c
c     where the input arguments are
c
c       x         (3 x 1) array of x coordinates of triangle nodes
c       y         (3 x 1) array of y coordinates of triangle nodes
c       db        (3 x 3) moment-curvature matrix.
c       f         factor by which stiffness entries will be multiplied.
c       ls        (9 x 1) array of stiffness location pointers
c                 (see output sm).
c       sm        incoming material stiffness array.
c       m         first dimension of sm in calling program.
c
c     the outputs are:
c
c       sm        output stiffness array with higher order stiffness
c                 coefficients added in.  the (i,j)-th entry of the
c                 (9 by 9) element bending stiffness is added to
c                 sm(k,l), where k=ls(i) and l=ls(j).
c       status    status character variable.  blank if no error
c                 detected.
c
c=end usage
c=block fortran
      subroutine    smcbh(x, y, db, f, ls, sm, m, status)
c
c                   t y p e   &   d i m e n s i o n
c
c=block vax
c     implicit      none
c=end vax
      integer m
      character*(*)  status
      double precision x(3),y(3),db(3,3), pg(3,3), rsd(6,6)
      double precision sm(m,m), sq(3,3), sds(3,3), q(6,9)
      double precision rm(3,6), l1, l2, l3
      double precision f, x0, y0, x1, x2, x3, y1, y2, y3
      double precision x21, x32, x13, y21, y32, y13, area, area2
      double precision l21, l32, l13, bl2, al2, bl3, al3, bl1, al1
      double precision cc, d11, d22, d33, d12, d13, d23
      double precision s1, s2, s3, s4, s5, s6
      integer       ls(9)
      integer       i, j, k, l
      data pg /0.d0,0.5d0,0.5d0,0.5d0,0.d0,0.5d0,0.5d0,0.5d0,0.d0/
c
c                   l o g i c
c
      status =   ' '
c     cleaning
c
      do 100 i = 1,3
        do 100 j = 1,6
        rm(i,j) = 0.0d0
100   continue
      do 200 i = 1,6
        do 200 j = 1,9
        q(i,j) = 0.0d0
200   continue
      do 300 i=1,6
        do 300 j=1,6
300       rsd(i,j)=0.0
      do 400 i=1,3
        do 400 j=1,3
400       sds(i,j)=0.0
c     coordinates
      x0=(x(1)+x(2)+x(3))/3.
      y0=(y(1)+y(2)+y(3))/3.
      x1= x(1) -x0
      x2= x(2) -x0
      x3= x(3) -x0
      y1= y(1) -y0
      y2= y(2) -y0
      y3= y(3) -y0
      x21 =      x2 - x1
      x32 =      x3 - x2
      x13 =      x1 - x3
      y21 =      y2 - y1
      y32 =      y3 - y2
      y13 =      y1 - y3
      area2 =    y21*x13 - x21*y13
      if (area2 .le. 0.0)      then
        status = 'nega_area'
        if (area2 .eq. 0.0)   status = 'zero_area'
        return
      end if
c          side lenghts
        l21 =  dsqrt( x21**2+y21**2 )
        l32 =  dsqrt( x32**2+y32**2 )
        l13 =  dsqrt( x13**2+y13**2 )
c          side proyections
        bl2=((x3-x1)*(x2-x1)+(y3-y1)*(y2-y1))/l21**2
        al2=1.0d0-bl2
        bl3=((x1-x2)*(x3-x2)+(y1-y2)*(y3-y2))/l32**2
        al3=1.0d0-bl3
        bl1=((x2-x3)*(x1-x3)+(y2-y3)*(y1-y3))/l13**2
        al1=1.0d0-bl1
c          inverse of the matrix relating inside curvatures
c          xx,yy,xy with boundary curvatures
c
        cc=area2**3
        sq(1,1)=( -x21*y21*y32**2 + x32*y21**2*y32)/cc
        sq(1,2)=(  x13*y13*y32**2 - x32*y13**2*y32)/cc
        sq(1,3)=(  x21*y21*y13**2 - x13*y21**2*y13)/cc
        sq(2,1)=(  x21*x32**2*y21 - x21**2*x32*y32)/cc
        sq(2,2)=( -x13*x32**2*y13 + x13**2*x32*y32)/cc
        sq(2,3)=( -x21*x13**2*y21 + x21**2*x13*y13)/cc
        sq(3,1)=(  x21**2*y32**2  - x32**2*y21**2) /cc
        sq(3,2)=( -x13**2*y32**2  + x32**2*y13**2) /cc
        sq(3,3)=( -x21**2*y13**2  + x13**2*y21**2) /cc
c     print '(''area'',f8.3)',area
c     print '(''inver'',3f8.3)',((sq(i,j),j=1,3),i=1,3)
      d11 =      db(1,1)
      d22 =      db(2,2)
      d33 =      db(3,3)
      d12 =      db(1,2)
      d13 =      db(1,3)
      d23 =      db(2,3)
      area =     0.5*area2
      do 2000 j=1,3
        s1=   d11*sq(1,j) + d12*sq(2,j) + d13*sq(3,j)
        s2=   d12*sq(1,j) + d22*sq(2,j) + d23*sq(3,j)
        s3=   d13*sq(1,j) + d23*sq(2,j) + d33*sq(3,j)
        do 1500 i=1,j
          sds(i,j)= sds(i,j)+ (s1*sq(1,i)+s2*sq(2,i)+s3*sq(3,i))
          sds(j,i)= sds(i,j)
1500    continue
2000  continue
      do 2100 j=1,3
        do 2100 i=1,3
2100      sds(i,j)=sds(i,j)*area/3.0*f
c     print '(''sds'',3f8.3)',((sds(i,j),j=1,3),i=1,3)
      
c    
c         matrix q
      q(1,1)=  6.0d0
      q(1,2)= -2.0*y13
      q(1,3)=  2.0*x13
      q(1,7)= -6.0d0
      q(1,8)= -4.0*y13
      q(1,9)=  4.0*x13
      q(2,1)= -6.0d0
      q(2,2)=  4.0*y13
      q(2,3)= -4.0*x13
      q(2,7)=  6.0d0
      q(2,8)=  2.0*y13
      q(2,9)= -2.0*x13
      q(3,1)= -6.0d0
      q(3,2)= -4.0*y21
      q(3,3)=  4.0*x21
      q(3,4)=  6.0d0
      q(3,5)= -2.0*y21
      q(3,6)=  2.0*x21
      q(4,1)=  6.0d0
      q(4,2)=  2.0*y21
      q(4,3)= -2.0*x21
      q(4,4)= -6.0d0
      q(4,5)=  4.0*y21
      q(4,6)= -4.0*x21
      q(5,4)= -6.0d0
      q(5,5)= -4.0*y32
      q(5,6)=  4.0*x32
      q(5,7)=  6.0d0
      q(5,8)= -2.0*y32
      q(5,9)=  2.0*x32
      q(6,4)=  6.0d0
      q(6,5)=  2.0*y32
      q(6,6)= -2.0*x32
      q(6,7)= -6.0d0
      q(6,8)=  4.0*y32
      q(6,9)= -4.0*x32
c     print '(''q'',9f5.1)',((q(i,j),j=1,9),i=1,6)
c
c    numerical integration
c
      do 2500 k=1,3
        l1=pg(k,1)
        l2=pg(k,2)
        l3=pg(k,3)
c
c    compute rm in the integration point
c
        rm(1,1)= l3 +al1 *l2 -(1.+al1)/3.
        rm(1,2)= l1 +bl1 *l2 -(1.+bl1)/3.
        rm(2,3)= l1 +al2 *l3 -(1.+al2)/3.
        rm(2,4)= l2 +bl2 *l3 -(1.+bl2)/3.
        rm(3,5)= l2 +al3 *l1 -(1.+al3)/3.
        rm(3,6)= l3 +bl3 *l1 -(1.+bl3)/3.
        do 2400 j=1,6
          s1=sds(1,1)*rm(1,j)+sds(1,2)*rm(2,j)+sds(1,3)*rm(3,j)
          s2=sds(2,1)*rm(1,j)+sds(2,2)*rm(2,j)+sds(2,3)*rm(3,j)
          s3=sds(3,1)*rm(1,j)+sds(3,2)*rm(2,j)+sds(3,3)*rm(3,j)
          do 2300 i=1,j
            rsd(i,j)=rsd(i,j)+(s1*rm(1,i)+s2*rm(2,i)+
     $                         s3*rm(3,i))
            rsd(j,i)=rsd(i,j)
 2300     continue
 2400   continue
 2500 continue
c     print '(''int'',6f7.3)',((rsd(i,j),j=1,6),i=1,6)
c
c    computing the kh matrix
c
      do 2600 j=1,9
        k =ls(j)
        s1=rsd(1,1)*q(1,j)+rsd(1,2)*q(2,j)+rsd(1,3)*q(3,j)
     $    +rsd(1,4)*q(4,j)+rsd(1,5)*q(5,j)+rsd(1,6)*q(6,j)
        s2=rsd(2,1)*q(1,j)+rsd(2,2)*q(2,j)+rsd(2,3)*q(3,j)
     $    +rsd(2,4)*q(4,j)+rsd(2,5)*q(5,j)+rsd(2,6)*q(6,j)
        s3=rsd(3,1)*q(1,j)+rsd(3,2)*q(2,j)+rsd(3,3)*q(3,j)
     $    +rsd(3,4)*q(4,j)+rsd(3,5)*q(5,j)+rsd(3,6)*q(6,j)
        s4=rsd(4,1)*q(1,j)+rsd(4,2)*q(2,j)+rsd(4,3)*q(3,j)
     $    +rsd(4,4)*q(4,j)+rsd(4,5)*q(5,j)+rsd(4,6)*q(6,j)
        s5=rsd(5,1)*q(1,j)+rsd(5,2)*q(2,j)+rsd(5,3)*q(3,j)
     $    +rsd(5,4)*q(4,j)+rsd(5,5)*q(5,j)+rsd(5,6)*q(6,j)
        s6=rsd(6,1)*q(1,j)+rsd(6,2)*q(2,j)+rsd(6,3)*q(3,j)
     $    +rsd(6,4)*q(4,j)+rsd(6,5)*q(5,j)+rsd(6,6)*q(6,j)
        do 2800 i=1,j
          l      =ls(i)
          sm(l,k)=sm(l,k)+(s1*q(1,i)+s2*q(2,i)+s3*q(3,i)
     $                    +s4*q(4,i)+s5*q(5,i)+s6*q(6,i))
          sm(k,l)=sm(l,k)
 2800   continue
 2600   continue  
c      write(46,5000) ((sm(i,j),i=1,9),j=1,9)
c5000  format(9f10.3)
      return
      end
c=end fortran
