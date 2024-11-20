	subroutine compthmfr(x,y,z,t,a1,a2,e1,e2,nu,h,phi,ta,numlay,
     &                       cmpfr,f,cfrm,globalflag)
c
c----------------------------------------------------------------------*
c
c This routine caclulates the membrane themally induced mechanical forces  
c for the 18 dof three node triangle composite shell element. 
c Coded by Joe Pajot on 2/11/03
c
c input variables:
c	x           = x coordinates of nodes
c       y           = y coordinates of nodes
c	z           = z coordinates of nodes
c       t           = mean temperature in element 
c	a1          = coefficent of thermal expansion 1(in array of length numlay)
c	a2          = coefficent of thermal expansion 2(in array of length numlay)
c	e1,2        = Young's modulus (in array of length numlay)
c       g12         = Shear modulus
c	nu          = Poisson's ratio (in array of length numlay)
c       h           = shell thickness (in array of length numlay)
c       phi         = fiber orientaion w.r.t. cmpfr (in array of length numlay)
c       ta          = reference temperature in each element
c       numlay      = number of layers in the composite
c       globalflag  = flag to return to local or global coordinates
c       cmpfr       = composite frame for orthotropic properties
c       cfrm        = frame number for definition of the shell     
c
c output variables:
c	f  = thermally induced mechanical force
c
c----------------------------------------------------------------------
C
	integer numlay,irot
	real*8 t
	real*8 a1(numlay),a2(numlay),e1(numlay),nu(numlay)
	real*8 h(numlay),phi(numlay),ta(numlay)
	real*8 x(3),y(3),z(3),f(18),cmpfr(9),e2(numlay)
	real*8 zvec(numlay+1),Qbar(3,3),Tr(3,3),invT(3,3)
	real*8 alpha,aframe1(3),aframe2(3),aframe3(3)
	real*8 qt1,qt2,qt3,thetaF,thetaD,theta
	real*8 xp(3),yp(3),zp(3),xlp(3),ylp(3),zlp(3)
	real*8 str(3),temp(3),ftemp(9),flocal(18)
	real*8 dm(3,3),dml(3,3),p(9,3),rot(3,3),R(3)
	real*8 v1n(3),v2n(3),v3n(3),refvec(3),orifiber(3)
	real*8 norm1 , norm2 , normref , proj1 , proj2
	real*8 x21,y21,z21,x32,y32,z32,x13,y13,z13
	real*8 x12,y12,z12,x23,y23,z23,x31,y31,z31
	real*8 cb,rlr,area2,coef1,coef2
	real*8 xlcg,ylcg,zlcg,ylr,zlr,xcg,ycg,zcg
	real*8 pi,twopi,cosine1,cosine2,costheta,sintheta
	integer i,j,k,m,globalflag,cfrm	
	data alpha/1.5D0/
	data v1n/1.0D0,0.0D0,0.0D0/
        data v2n/0.0D0,1.0D0,0.0D0/
        data v3n/0.0D0,0.0D0,1.0D0/
c	
	pi    = acos(-1.0D0)
        twopi = 2.000000D+00*pi
	thetaF = 0.0D+00
	thetaD = 0.0D+00
	theta  = 0.0D+00
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
        xcg = (x(1) + x(2) + x(3))/3.0d+00
        ycg = (y(1) + y(2) + y(3))/3.0d+00
        zcg = (z(1) + z(2) + z(3))/3.0d+00
C compute local coordinates 
        do  i=1,3
          xlcg   = x(i) - xcg
          ylcg   = y(i) - ycg
          zlcg   = z(i) - zcg
          xlp(i) = xp(1) * xlcg + xp(2) * ylcg + xp(3) * zlcg
          ylp(i) = yp(1) * xlcg + yp(2) * ylcg + yp(3) * zlcg
          zlp(i) = zp(1) * xlcg + zp(2) * ylcg + zp(3) * zlcg
  	end do
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
c  create strain-displacement matrix p (assign half and dV terms to cb)
c
      p(1,1) =   y23
      p(2,1) =   0.0D0
      p(3,1) =   y31
      p(4,1) =   0.0D0
      p(5,1) =   y12
      p(6,1) =   0.0D0
      p(1,2) =   0.0D0
      p(2,2) =   x32
      p(3,2) =   0.0D0
      p(4,2) =   x13
      p(5,2) =   0.0D0
      p(6,2) =   x21
      p(1,3) =   x32
      p(2,3) =   y23
      p(3,3) =   x13
      p(4,3) =   y31
      p(5,3) =   x21
      p(6,3) =   y12
      coef1  = alpha/6.0
      coef2  = alpha/3.0
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
c     create the z-vector
      zvec(1) = 0.0D0
      do i= 1,numlay
         zvec(1) = zvec(1) + h(i)/2.0D+00
      end do
      do i = 1,numlay
         zvec(1+i) = zvec(i) - h(i)
      end do		      
c
c     zero the local force vector 
c
      do i=1,18
         flocal(i) = 0.0D0
      end do 
c
c create the requisite frames
c
      if ( cfrm .eq. 0 ) then
C
C.....INITIALIZE WITH THE IDENTITY IF THE FRAME NUMBER IS ZERO
C
         aframe1(1) = 1.0D0
         aframe1(2) = 0.0D0
         aframe1(3) = 0.0D0
         aframe2(1) = 0.0D0
         aframe2(2) = 1.0D0
         aframe2(3) = 0.0D0
         aframe3(1) = 0.0D0
         aframe3(2) = 0.0D0
         aframe3(3) = 1.0D0
C
      else
C
         aframe1(1) = cmpfr(1)
         aframe1(2) = cmpfr(2)
         aframe1(3) = cmpfr(3)
         aframe2(1) = cmpfr(4)
         aframe2(2) = cmpfr(5)
         aframe2(3) = cmpfr(6)
         aframe3(1) = cmpfr(7)
         aframe3(2) = cmpfr(8)
         aframe3(3) = cmpfr(9)
C
      endif
c
      do i=1,3
        do j =1,3
         rot(i,j) = 0.0D+00
	end do
      end do
c     
      call rotation(xp,yp,zp,v1n,v2n,v3n,rot)
c
C ....SET THE REFERENCE VECTOR FOR ORIENTING THE FIBERS OF THE LAYER
C ....(ALWAYS TAKE THE FIRST VECTOR OF THE FRAME)
C
      refvec(1) = aframe1(1)
      refvec(2) = aframe1(2)
      refvec(3) = aframe1(3)
C
C ....PROJECT THE REFERENCE VECTOR INTO THE PLANE OF
C ....THE ELEMENT TO GET THE FIBER ORIENTATION VECTOR
C
      orifiber(1) = 0.0D+00
      orifiber(2) = 0.0D+00
      orifiber(3) = 0.0D+00
C
      norm1   = 0.0D+00
      norm2   = 0.0D+00
      normref = 0.0D+00
      proj1   = 0.0D+00
      proj2   = 0.0D+00
C
      do j=1,3
       norm1   =   norm1 + rot(j,1)*rot(j,1)
       norm2   =   norm2 + rot(j,2)*rot(j,2)
       normref = normref + refvec(j)*refvec(j)
       proj1   =   proj1 + rot(j,1)*refvec(j)
       proj2   =   proj2 + rot(j,2)*refvec(j)
      end do
C
      norm1   = sqrt(norm1)
      norm2   = sqrt(norm2)
      normref = sqrt(normref)
c             
c
      if ( normref.eq.0.0D+00 ) then
  	 cosine1 = 1.0D+00
  	 cosine2 = 0.0D+00
      else
  	 cosine1 = proj1/(norm1*normref)
  	 cosine2 = proj2/(norm2*normref)
      endif
C
      do j=1,3
  	 orifiber(j) = cosine1*rot(j,1) + cosine2*rot(j,2)
      end do
C
C ....CALCULATE THE ANGLE FROM THE [x] FRAME OF THE LOCAL
C ....COORDINATE SYSTEM TO THE REFERENCE ORIENTATION VECTOR
C
      proj1  = 0.0D+00
      proj2  = 0.0D+00
      thetaD = 0.0D+00
C
      do  j=1,3
  	 proj1 = proj1 + rot(j,1)*orifiber(j)
  	 proj2 = proj2 + rot(j,2)*orifiber(j)
      end do
C
      if ( proj1.eq.0.0D+00 ) then
  	 if ( proj2.eq.0.0D+00 ) go to 700
  	 if ( proj2.gt.0.0D+00 ) thetaD = 0.50D+00*pi
  	 if ( proj2.lt.0.0d+00 ) thetaD = 1.50D+00*pi
      else
  	 if ( proj2.eq.0.0D+00 ) then
  	    if ( proj1.eq.0.0D+00 ) go to 700
  	    if ( proj1.gt.0.0D+00 ) thetaD = 0.0D+00
  	    if ( proj1.lt.0.0D+00 ) thetaD = pi
  	 else
  	    thetaD = atan(proj2/proj1)
  	 endif
      endif
C
      if ( thetaD.lt.0.0D+00 ) then
  	 thetaD = thetaD + twopi
      endif
C
      if ( (thetaD.lt.0.0D+00).or.(thetaD.gt.twopi) ) go to 800
c
C
        do j=1,3
           do k=1,3
              invT(k,j) = 0.0D+00
	      dm(k,j) = 0.0D+00
	      Tr(k,j) = 0.0D+00
           end do
        end do
C
c     loop over all shell layers and compute forces and moments     
c
      do i=1,numlay
c	
c       TRANSFORM ANGLE IN THE RANGE BETWEEN 0-360      
C       
        thetaF = phi(i)
        if ( (thetaF.lt.0.0).or.(thetaF.gt.360.00D+00) ) then 
           irot = thetaF / 360.0D0
           thetaF = thetaF - real(irot) * 360.0D0
           if (thetaF.lt.0.0D0) thetaF = 360.0D0 + thetaF
        endif
C
C   .....TRANSFORM FROM DEGREE TO RADIAN THE ANGLE BETWEEN THE
C   .....REFERENCE ORIENTATION VECTOR AND THE ORIENTATION OF THE FIBERS
C
        thetaF = (pi*thetaF)/180.00D+00
C
C  .....CALCULATE THE ANGLE FROM THE [x] FRAME OF THE LOCAL
C.  ....COORDINATE SYSTEM TO THE DIRECTION OF THE FIBER
C
        theta = thetaD - thetaF
c
C  .....INITIALIZE THE ROTATION MATRIX FROM THE FIBER COORDINATE
C  .....SYSTEM {1;2} TO THE ELEMENT TRIANGULAR SYSTEM {x;y}
C
        costheta = cos(theta)
        sintheta = sin(theta)
C
        Tr(1,1)   =  costheta
        Tr(1,2)   =  1.0D+00*sintheta
        Tr(1,3)   =  0.0D+00
C
        Tr(2,1)   =  -1.0D+00*sintheta
        Tr(2,2)   =  costheta
        Tr(2,3)   = 0.0D+00
C
        Tr(3,1)   = 0.0D+00
        Tr(3,2)   = 0.0D+00
        Tr(3,3)   = 1.0D+00
C
        do j=1,3
           do k=1,3
              invT(k,j) = Tr(j,k)
           end do
        end do
C
c       Calculate stress tensor in fiber coordinates
c
        cb = (e2(i)*(h(i)/2))/(e1(i)-
     *        e2(i)*nu(i)*nu(i))
c     
        dm(1,1) = e1(i)*e1(i)*a1(i)*cb/e2(i)
	dm(1,2) = e1(i)*a2(i)*cb*nu(i)
	dm(1,3) = 0.0D+00
	dm(2,1) = e1(i)*cb*nu(i)*a1(i)
	dm(2,2) = e1(i)*cb*a2(i)
	dm(2,3) = 0.0D+00
	dm(3,1) = 0.0D+00
	dm(3,2) = 0.0D+00
	dm(3,3) = 0.0D+00
c
c 	create vector temp = (D*str)_global
c       
	do j=1,3
	  temp(j) = dm(j,1)*(t-ta(i)) + 
     *	            dm(j,2)*(t-ta(i)) + dm(j,3)*(t-ta(i))
	end do
c
c       write out stress tensor
        dm(1,1) = temp(1)
	dm(1,2) = temp(3)
	dm(1,3) = 0.0D+00
	dm(2,1) = temp(3)
	dm(2,2) = temp(2)
	dm(2,3) = 0.0D+00
	dm(3,1) = 0.0D+00
	dm(3,2) = 0.0D+00
	dm(3,3) = 0.0D+00
c
c       transform stress tensor to element local coordinates
c
        do j=1,3
C
          qt1 = dm(1,1)*invT(j,1) + dm(1,2)*invT(j,2) + 
     *	        dm(1,3)*invT(j,3)
          qt2 = dm(2,1)*invT(j,1) + dm(2,2)*invT(j,2) + 
     *          dm(2,3)*invT(j,3)
          qt3 = dm(3,1)*invT(j,1) + dm(3,2)*invT(j,2) + 
     *	        dm(3,3)*invT(j,3)
C
          do k=1,3
            Qbar(k,j) = qt1*invT(k,1) + qt2*invT(k,2) + qt3*invT(k,3)
          end do
C
        end do
c	
c 	create vector temp = (D*str)_global
c       
	temp(1) = Qbar(1,1) 
	temp(2) = Qbar(2,2)
	temp(3) = Qbar(2,1)
c
c	
	do j=1,9
	  ftemp(j) = p(j,1)*temp(1) + p(j,2)*temp(2) + 
     *               p(j,3)*temp(3)
	end do
c	
c	 add moment and force contributions from this layer to flocal				 
c	
         do j = 1,3
	 	 flocal(j*6-5) = flocal(j*6-5) + ftemp(2*j-1)
	 	 flocal(j*6-4) = flocal(j*6-4) + ftemp(2*j  )
	 	 flocal(j*6-2) = flocal(j*6-2) - 
     *	 	                 ftemp(2*j)*(zvec(i) + zvec(i+1))/2.0
	 	 flocal(j*6-1) = flocal(j*6-1) +		 
     *	 	                 ftemp(2*j-1)*(zvec(i) + zvec(i+1))/2.0
	 	 flocal(j*6  ) = flocal(j*6  ) + ftemp(6+j)
	
	 end do
c	 
      end do
c
      if (globalflag .eq. 0) then
c       transorm force local to force global
c	first create the transformation matrix which rot,
c	next carry out this multipication node by node which has 
c       forces [flx fly 0 mlx mly mlz] 
	do i =1,3
	   f(i*6-5) =  rot(1,1)*flocal(i*6-5) + rot(1,2)*flocal(i*6-4)
	   f(i*6-4) =  rot(2,1)*flocal(i*6-5) + rot(2,2)*flocal(i*6-4)
	   f(i*6-3) =  rot(3,1)*flocal(i*6-5) + rot(3,2)*flocal(i*6-4)
	   f(i*6-2) =  rot(1,1)*flocal(i*6-2) + rot(1,2)*flocal(i*6-1) + 
     *	               rot(1,3)*flocal(i*6)   
	   f(i*6-1) =  rot(2,1)*flocal(i*6-2) + rot(2,2)*flocal(i*6-1) +
     *                 rot(2,3)*flocal(i*6)   
	   f(i*6  ) =  rot(3,1)*flocal(i*6-2) + rot(3,2)*flocal(i*6-1) +
     *                 rot(3,3)*flocal(i*6)
	end do
      else 
         do i=1,3
	   f(i*6-5) = flocal(i*6-5)
	   f(i*6-4) = flocal(i*6-4)
	   f(i*6-3) = flocal(i*6-3)
	   f(i*6-2) = flocal(i*6-2)
	   f(i*6-1) = flocal(i*6-1)
	   f(i*6  ) = flocal(i*6  )
	 end do
      end if 
c
      return
c      
  700 continue
      write(*,*) "*** FATAL ERROR in Routine COMPTHMFR ***"
      write(*,*) "*** The Reference Orientation Vector ***"
      write(*,*) "*** is Parallel to the Two Inplane   ***"
      write(*,*) "*** and Orthogonal Local Frames!     ***"
      write(*,*) "*** STOP ALL TREATMENTS RIGHT HERE   ***"
      stop
c
C
  800 continue
      write(*,*) "*** FATAL ERROR in Routine COMPTHMFR  ***"
      write(*,*) "*** The Angle From the Local [x]      ***"
      write(*,*) "*** Axis of the Triangular Coordinate ***"
      write(*,*) "*** System to the Reference Direction ***"
      write(*,*) "*** is Out-of-Bounds: it Must be      ***"
      write(*,*) "*** Within the Range 0-2pi Radians    ***"
      write(*,*) "*** STOP ALL TREATMENTS RIGHT HERE    ***"
      stop
      
      end 
