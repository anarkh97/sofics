        subroutine trianarea(x,y,z,area)
c
c-------------------------------------------------------------------*
c input variables:
c 	x = x coordinates
c 	y = y coordinates
c 	z = z coordinates
c
c output variable:
c	area = triangle area
c 
*-------------------------------------------------------------------*
*
        real*8 x(3),y(3),z(3)
        real*8 x21,y21,z21,x32,y32,z32,x13,y13,z13
        real*8 rlr,rlb,bpr,area

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
        rlb = dsqrt( x32*x32 + y32*y32 + z32*z32 )
        bpr = dsqrt((x21 * x32 + y21 * y32 + z21 *z32 )**2)/rlr
        area= rlr*(dsqrt(rlb**2-bpr**2))/2.0d+00
C
      return
      end

