        subroutine trimem(flag, xl,yl,zl,e,nu,h,rk)
C
C THIS SUBROUTINE EVALUATES THE STIFFNESS MATRIX 
C FOR THE SPATIAL 9 D.O.F THREE NODE TRIANGLE MEMBRANE   
C
C.... DECLARE ALL GLOBAL VARAIBALES
C
C.... REAL CONSTANTS
C
        real*8 e, nu
C
C.... REAL ARRAYS
C
        real*8 xl(3),yl(3),zl(3),h(3),sm(9,9),rk(18,18)
C
C.... DECLARE ALL LOCAL VARIABLES
C
        real*8 dm(3,3),r1(3,3)
        real*8 xp(3),yp(3),zp(3),xlp(3),ylp(3),zlp(3)
        real*8 v1n(3),v2n(3),v3n(3)
        real*8 cb,x21,y21,z21,x32,y32,z32,x13,z13
        real*8 rx,ry,rz,bx,by,bz,rlr,rlb,bpr,area,ylr,zlr
        real*8 ycg,xcg,zcg,xlcg,ylcg,zlcg,esp
        integer le(9)
        integer ptr(9)
C
        integer i,j,flag
        character*10 status
        data le/1,2,7,8,13,14,6,12,18/
        data ptr/1,2,6,7,8,12,13,14,18/
        data v1n/1.0,0.0,0.0/
        data v2n/0.0,1.0,0.0/
        data v3n/0.0,0.0,1.0/
C
C.... thickness
C
        esp = h(1)
C
C.... elastic matrix
C
        cb=(e*esp)/(1.0-(nu*nu))
        dm(1,1) = cb
        dm(1,2) = nu*cb
        dm(1,3) = 0.0
        dm(2,1) = dm(1,2)
        dm(2,2) = cb
        dm(2,3) = 0.0
        dm(3,1) = 0.0
        dm(3,2) = 0.0
        dm(3,3) = ((1.0-nu)/2.0)*cb
C
C.... dimension variables
C
        x21 = xl(2) - xl(1)
        y21 = yl(2) - yl(1)
        z21 = zl(2) - zl(1)
        x32 = xl(3) - xl(2)
        y32 = yl(3) - yl(2)
        z32 = zl(3) - zl(2)
        x13 = xl(1) - xl(3)
        z13 = zl(1) - zl(3)
C
C.... triangle in space : we compute the length of one 
C.... side and the distance of the
C.... opposing node to that side to compute the area
C
        rx   = x21
        ry   = y21
        rz   = z21
        bx   = x32
        by   = y32
        bz   = z32
C
        rlr  = dsqrt( rx*rx + ry*ry + rz*rz )
        rlb  = dsqrt( bx*bx + by*by + bz*bz )
        bpr  = dsqrt((rx * bx + ry * by + rz *bz )**2)/rlr
        area = 0.5*rlr*(dsqrt(rlb*rlb-bpr*bpr))
C
C.... direction cosines of the local system . 
C.... X' is dircted parallel to the 2-1 side
C.... Z' is the external normal (counterclockwise). 
C.... Y' computed as Z'x X'
C
        xp(1) = x21/rlr
        xp(2) = y21/rlr
        xp(3) = z21/rlr
C
C.... Z'
C
        zp(1) = y21 * z32 - z21 * y32
        zp(2) = z21 * x32 - x21 * z32
        zp(3) = x21 * y32 - y21 * x32
        zlr   = dsqrt( zp(1)*zp(1) + zp(2)*zp(2) + zp(3)*zp(3) )
        zp(1) = zp(1)/zlr
        zp(2) = zp(2)/zlr
        zp(3) = zp(3)/zlr
C
C.... Y'
C
        yp(1) = zp(2) * xp(3) - zp(3) * xp(2)
        yp(2) = zp(3) * xp(1) - zp(1) * xp(3)
        yp(3) = zp(1) * xp(2) - zp(2) * xp(1)
        ylr   = dsqrt( yp(1)*yp(1) + yp(2)*yp(2) + yp(3)*yp(3) )
        yp(1) = yp(1)/ylr
        yp(2) = yp(2)/ylr
        yp(3) = yp(3)/ylr
C
C.... center of gravity
C
        xcg = (xl(1) + xl(2) + xl(3))/3.0d+00
        ycg = (yl(1) + yl(2) + yl(3))/3.0d+00
        zcg = (zl(1) + zl(2) + zl(3))/3.0d+00
C
C.... computing local coordinates 
C
        do 43 i=1,3
          xlcg   = xl(i)-xcg
          ylcg   = yl(i)-ycg
          zlcg   = zl(i)-zcg
          xlp(i) = xp(1) * xlcg + xp(2) * ylcg + xp(3) * zlcg
          ylp(i) = yp(1) * xlcg + yp(2) * ylcg + yp(3) * zlcg
          zlp(i) = zp(1) * xlcg + zp(2) * ylcg + zp(3) * zlcg
43      continue
C
C.... cleaning stiffness matrix
C
        do 15 i=1,18
           do 15 j=1,18
             rk(i,j) = 0.0d+00
15      continue
C
C.... forming local basic stiffness for membrane
C
        call sm3mb(xlp,ylp,dm,1.5d+00,1.0d+00,le,rk,18,status)
C
C.... forming local higher order stiffness for membrane
C
        call sm3mhe(xlp,ylp,dm,0.32d+00,le,rk,18,status)
C
C rotate stiffness matrix from local coordinate system to global
C coordinate system in the case of linear FEM. In the case of
C nonlinear FEM with corotational method, do not perform this 
C transformation as the corotational routines expect a stiffness
C matrix in local coordinates.
C
        if(flag .eq. 1) then
C
C.... computing nodal rotation matrices
C
        call rotation(xp,yp,zp,v1n,v2n,v3n,r1)
C
C.... rotate membrane stiffness
C
        call trirotation(rk,r1)

        endif
C
C        do 60 i=1,9
C          sm(i,1) = rk(ptr(i), 1)
C          sm(i,2) = rk(ptr(i), 2)
C          sm(i,3) = rk(ptr(i), 6)
C          sm(i,4) = rk(ptr(i), 7)
C          sm(i,5) = rk(ptr(i), 8)
C          sm(i,6) = rk(ptr(i),12)
C          sm(i,7) = rk(ptr(i),13)
C          sm(i,8) = rk(ptr(i),14)
C          sm(i,9) = rk(ptr(i),18)
C60      continue

        return
        end
