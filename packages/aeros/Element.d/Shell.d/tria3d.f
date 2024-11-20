        subroutine tria3d(flag, xl,yl,zl,e,nu,h,rk)
c
c-------------------------------------------------------------------*
c This subroutine evaluates the stiffness matrix and mass vector
c for the spatial 18 d.o.f tree node triangle obtained as
c a combination of the aqr bending triangle plus the membrane
c with drilling d.o.f. developed by Felippa et al.
c
c        rkb will be used for the basic stiffness 
c	 rkm will be used for the higher order stiffness
c
c input variables:
c 	xl = x coordinates
c 	yl = y coordinates
c 	zl = z coordinates
c 	e  = elastic modulus
c 	nu = poisson's ratio
c 	h  = thickness
c
c output variables:
c	rk = stiffness matrix
c 
*-------------------------------------------------------------------*
*
*  SUBROUTINES CALLED: BASICO
*                      SM3MB
*                      SMCBH
*                      SM3MHEFF
*                      ROTATION
*                      TRIROT
*-------------------------------------------------------------------*
C
        real*8 e,nu,alpha
        real*8 xl(3),yl(3),zl(3),h(3),rk(18,18)
        real*8 db(3,3),dm(3,3)
        real*8 xp(3),yp(3),zp(3),xlp(3),ylp(3),zlp(3)
        real*8 r1(3,3),rkm(18,18),rkb(18,18)
        real*8 v1n(3),v2n(3),v3n(3)
        real*8 cb,x21,y21,z21,x32,y32,z32,x13,y13,z13
        real*8 rlr,rlb,bpr,area,ylr,zlr
        real*8 ycg,xcg,zcg,xlcg,ylcg,zlcg,f,clr,cqr,esp
        integer lb(9),le(9)
        integer i,j, flag
        character *10 status
        data lb/3,4,5,9,10,11,15,16,17/
        data le/1,2,7,8,13,14,6,12,18/
        data v1n/1.0,0.0,0.0/
        data v2n/0.0,1.0,0.0/
        data v3n/0.0,0.0,1.0/

        PARAMETER(f = 1.0d+00, clr = 0.0d+00, cqr = 1.0d+00)
        PARAMETER(alpha = 1.5d+00)
c
c h   = shell thickness
c e   = elastic modulus
c nu  = poisson's ratio
c
c flag = 0 do NOT perform transform from local to global
c flag = 1 perform transformation from local to global
c
c set thickness
c
        esp = h(1)
c
c bending elastic matrix
c
        cb      = e*esp**3/12.0/(1.0-(nu*nu))
        db(1,1) =    cb
        db(1,2) = nu*cb
        db(1,3) = 0.0
        db(2,1) = db(1,2)
        db(2,2) = cb
        db(2,3) = 0.0
        db(3,1) = 0.0
        db(3,2) = 0.0
        db(3,3) = ((1.0-nu)/2.0)*cb
c
c membrane elastic matrix
c
        cb=e*esp/(1.0-(nu*nu))
        dm(1,1) =    cb
        dm(1,2) = nu*cb
        dm(1,3) = 0.0
        dm(2,1) = dm(1,2)
        dm(2,2) = cb
        dm(2,3) = 0.0
        dm(3,1) = 0.0
        dm(3,2) = 0.0
        dm(3,3) = ((1.0-nu)/2.0)*cb
C
C dimension variables
C
        x21 = xl(2) - xl(1)
        y21 = yl(2) - yl(1)
        z21 = zl(2) - zl(1)
        x32 = xl(3) - xl(2)
        y32 = yl(3) - yl(2)
        z32 = zl(3) - zl(2)
        x13 = xl(1) - xl(3)
        y13 = yl(1) - yl(3)
        z13 = zl(1) - zl(3)
C triangle in space : we compute the length of one side and the distance of the
C opposing node to that side to compute the area
        rlr = dsqrt( x21*x21 + y21*y21 + z21*z21 )
        rlb = dsqrt( x32*x32 + y32*y32 + z32*z32 )
        bpr = dsqrt((x21 * x32 + y21 * y32 + z21 *z32 )**2)/rlr
        area= rlr*(dsqrt(rlb**2-bpr**2))/2.0d+00
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
        xcg = (xl(1) + xl(2) + xl(3))/3.0d+00
        ycg = (yl(1) + yl(2) + yl(3))/3.0d+00
        zcg = (zl(1) + zl(2) + zl(3))/3.0d+00
C compute local coordinates 
        do 43 i=1,3
          xlcg   = xl(i) - xcg
          ylcg   = yl(i) - ycg
          zlcg   = zl(i) - zcg
          xlp(i) = xp(1) * xlcg + xp(2) * ylcg + xp(3) * zlcg
          ylp(i) = yp(1) * xlcg + yp(2) * ylcg + yp(3) * zlcg
          zlp(i) = zp(1) * xlcg + zp(2) * ylcg + zp(3) * zlcg
  43    continue
c
c zero stiffness matrices
c
        do 15 i=1,18
          do 15 j=1,18
            rk(i,j)  = 0.0d+00
            rkm(i,j) = 0.0d+00
15          rkb(i,j) = 0.0d+00
c
c form local basic bending stiffness
c
      call basico(xlp,ylp,db,1.0d+00,clr,cqr,lb,rkb,18,status)
c
c form local basic membrane stiffness
c
      call sm3mb(xlp,ylp,dm,alpha,1.0d+00,le,rkm,18,status)
c
c form local higher order bending stiffness
c
      call smcbh(xlp,ylp,db,1.0d+00,lb,rkb,18,status)
c
c form local higher order membrane stiffness
c
      call sm3mhe(xlp,ylp,dm,0.32d+00,le,rkm,18,status)
C
C add bending stiffness and membrane stiffness
C
      do 16 i=1,18
        do 16 j=1,18
          rk(i,j) = rkb(i,j)+rkm(i,j)
  16  continue
C
C rotate stiffness matrix from local coordinate system to global
C coordinate system in the case of linear FEM. In the case of
C nonlinear FEM with corotational method, do not perform this 
C transformation as the corotational routines expect a stiffness
C matrix in local coordinates.
C
      if(flag .eq. 1) then
C
C     compute local to global rotation matrix
C
        call rotation(xp,yp,zp,v1n,v2n,v3n,r1)

        call trirotation(rk,r1)

      endif
C
      return
      end

