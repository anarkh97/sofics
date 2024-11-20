        subroutine trithmfr(x,y,z,t,a,e,nu,h,alpha,f,globflag)
c
c----------------------------------------------------------------------*
c
c This routine caclulates the membrane themally induced mechanical force  
c for the 18 dof three node triangle shell element. 
c Coded by Joe Pajot on 1/16/03
c
c input variables:
c	x  = x coordinates of nodes
c       y  = y coordinates of nodes
c	z  = z coordinates of nodes
c       t  = mean temperature difference in element 
c	a  = coefficent of thermal expansion
c	e  = Young's modulus
c	nu = Poisson's ratio
c       h  = shell thickness
c       globflag = flag to return to global coordinates
c
c output variables:
c	f  = thermally induced mechanical force
c
c----------------------------------------------------------------------
C
        double precision t,a,e,nu,h,alpha
        double precision x(*),y(*),z(*),f(*)
        double precision xp(3),yp(3),zp(3),xlp(3),ylp(3),zlp(3)
        double precision str(3),temp(3),ftemp(9)
        double precision dm(3,3),p(9,3),rot(3,3)
        double precision v1n(3),v2n(3),v3n(3)
        double precision x21,y21,z21,x32,y32,z32,x13,y13,z13
        double precision x12,y12,x23,y23,x31,y31
        double precision cb,rlr,area2,coef1,coef2
c       double precision rlb,bpr,area
        double precision xlcg,ylcg,zlcg,ylr,zlr,xcg,ycg,zcg
        integer i,j,globflag
        data v1n/1.0d0,0.0d0,0.0d0/
        data v2n/0.0d0,1.0d0,0.0d0/
        data v3n/0.0d0,0.0d0,1.0d0/
c        integer i,j,globflag
C
C dimension variables
C
        x21 = x(2) - x(1)
        y21 = y(2) - y(1)
        z21 = z(2) - z(1)
        x32 = x(3) - x(2)
        y32 = y(3) - y(2)
        z32 = z(3) - z(2)
        x13 = x(1) - x(3)
        y13 = y(1) - y(3)
        z13 = z(1) - z(3)
C triangle in space : we compute the length of one side and the distance of the
C opposing node to that side to compute the area
        rlr = dsqrt( x21*x21 + y21*y21 + z21*z21 )
c       rlb = dsqrt( x32*x32 + y32*y32 + z32*z32 )
c       bpr = dsqrt((x21 * x32 + y21 * y32 + z21 *z32 )**2)/rlr
c       area= rlr*(dsqrt(rlb**2-bpr**2))/2.0d+00
C direction cosines of the local system . X' is directed parallel 
C to the 2-1 side
C Z' is the external normal (counterclockwise). Y' computed as Z' x X'
        xp(1) = x21/rlr
        xp(2) = y21/rlr
        xp(3) = z21/rlr
C Z'
        zp(1) = y21 * z32 - z21 * y32
        zp(2) = z21 * x32 - x21 * z32
        zp(3) = x21 * y32 - y21 * x32
        zlr   = dsqrt( zp(1)*zp(1) + zp(2)*zp(2)+ zp(3)*zp(3) )
        zp(1) = zp(1)/zlr
        zp(2) = zp(2)/zlr
        zp(3) = zp(3)/zlr
C Y'
        yp(1) = zp(2) * xp(3) - zp(3) * xp(2)
        yp(2) = zp(3) * xp(1) - zp(1) * xp(3)
        yp(3) = zp(1) * xp(2) - zp(2) * xp(1)
        ylr   = dsqrt( yp(1)*yp(1) + yp(2)*yp(2) + yp(3)*yp(3) )
        yp(1) = yp(1)/ylr
        yp(2) = yp(2)/ylr
        yp(3) = yp(3)/ylr
C compute center of gravity
        xcg = (x(1) + x(2) + x(3))/3.0d0
        ycg = (y(1) + y(2) + y(3))/3.0d0
        zcg = (z(1) + z(2) + z(3))/3.0d0
C compute local coordinates 
        do 10, i=1,3
          xlcg   = x(i) - xcg
          ylcg   = y(i) - ycg
          zlcg   = z(i) - zcg
          xlp(i) = xp(1) * xlcg + xp(2) * ylcg + xp(3) * zlcg
          ylp(i) = yp(1) * xlcg + yp(2) * ylcg + yp(3) * zlcg
          zlp(i) = zp(1) * xlcg + zp(2) * ylcg + zp(3) * zlcg
10      continue
c        end do
c
c  dimension variables in local coordinates
c		
      x21 =      xlp(2) - xlp(1)
      x12 =     -x21
      x32 =      xlp(3) - xlp(2)
      x23 =     -x32
      x13 =      xlp(1) - xlp(3)
      x31 =     -x13
      y21 =      ylp(2) - ylp(1)
      y12 =     -y21
      y32 =      ylp(3) - ylp(2)
      y23 =     -y32
      y13 =      ylp(1) - ylp(3)
      y31 =     -y13
      area2 =    y21*x13 - x21*y13
c
c membrane elastic matrix
c
      cb=e*(h/2.0d0)/(1.0d0-(nu*nu))
      dm(1,1) =    cb
      dm(1,2) = nu*cb
      dm(1,3) = 0.0d0
      dm(2,1) = dm(1,2)
      dm(2,2) = cb
      dm(2,3) = 0.0d0
      dm(3,1) = 0.0d0
      dm(3,2) = 0.0d0
      dm(3,3) = ((1.0d0-nu)/2.0d0)*cb
c
c  create strain vector
c
      str(1) = a*t 
      str(2) = a*t 
      str(3) = 0.0d0
c
c  create strain-displacement matrix p (assign half and dV terms to cb)
c
      p(1,1) =   y23
      p(2,1) =   0.0d0
      p(3,1) =   y31
      p(4,1) =   0.0d0
      p(5,1) =   y12
      p(6,1) =   0.0d0
      p(1,2) =   0.0d0
      p(2,2) =   x32
      p(3,2) =   0.0d0
      p(4,2) =   x13
      p(5,2) =   0.0d0
      p(6,2) =   x21
      p(1,3) =   x32
      p(2,3) =   y23
      p(3,3) =   x13
      p(4,3) =   y31
      p(5,3) =   x21
      p(6,3) =   y12
      coef1  = alpha/6.0d0
      coef2  = alpha/3.0d0
      p(7,1) =  y23*(y13-y21)*coef1
      p(7,2) =  x32*(x31-x12)*coef1
      p(7,3) =  (x31*y13-x12*y21)*coef2
      p(8,1) =  y31*(y21-y32)*coef1
      p(8,2) =  x13*(x12-x23)*coef1
      p(8,3) =  (x12*y21-x23*y32)*coef2
      p(9,1) =  y12*(y32-y13)*coef1
      p(9,2) =  x21*(x23-x31)*coef1
      p(9,3) =  (x23*y32-x31*y13)*coef2
c
c create vector = p'*D*str
c            
      do 20, i=1,3
        temp(i) = dm(i,1)*str(1) + dm(i,2)*str(2) + dm(i,3)*str(3)
20    continue
      do 30, j=1,9
        ftemp(j) = p(j,1)*temp(1) + p(j,2)*temp(2) + p(j,3)*temp(3)
30    continue
c
      if (globflag .eq. 0) then
c transform to global coordinates	
c	first create the transformation matrix (rot),
c	next carry out this multiplication node by node which has 
c       forces [flx fly mz] 
        call rotation(xp,yp,zp,v1n,v2n,v3n,rot)
        do 40, i =1,3
           f(i*6-5) =  rot(1,1)*ftemp(i*2-1) + rot(1,2)*ftemp(i*2)
           f(i*6-4) =  rot(2,1)*ftemp(i*2-1) + rot(2,2)*ftemp(i*2)
           f(i*6-3) =  rot(3,1)*ftemp(i*2-1) + rot(3,2)*ftemp(i*2)
           f(i*6-2) =  rot(1,3)*ftemp(6+i)
           f(i*6-1) =  rot(2,3)*ftemp(6+i)
           f(i*6  ) =  rot(3,3)*ftemp(6+i)
40      continue
      else 
         do 50, i=1,3
           f(i*6-5) = ftemp(i*2-1)
           f(i*6-4) = ftemp(i*2  )
           f(i*6-3) = 0
           f(i*6-2) = 0
           f(i*6-1) = 0
           f(i*6  ) = ftemp(6+i)
50       continue
      end if 
      return
      end 

