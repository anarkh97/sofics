CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C This subroutine computes the local coordinates and the area for
C the 3 node shell element
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C
        subroutine localxy(x, y, z, xl, yl, localxyflag, area)
C --------------------------------------------------------------

        integer i,localxyflag
        real*8 x(3), y(3), z(3)
        real*8 xp(3), yp(3), zp(3), xl(3), yl(3), zl(3)
        real*8 x21,y21,z21,x32,y32,z32,x13,y13,z13
        real*8 ycg, xcg, zcg, xlcg, ylcg, zlcg
        real*8 rlr, area, ylr, zlr, v1, v2, v3

C dimension variables

          x21 = x(2)-x(1)
          y21 = y(2)-y(1)
          z21 = z(2)-z(1)
          x32 = x(3)-x(2)
          y32 = y(3)-y(2)
          z32 = z(3)-z(2)
          x13 = x(1)-x(3)
          y13 = y(1)-y(3)
          z13 = z(1)-z(3)

c vector product

         v1 = y21*z13 - z21*y13
         v2 = z21*x13 - x21*z13
         v3 = x21*y13 - y21*x13
          
c area of triangle = 0.5*norm(v)

         area = dsqrt(v1*v1 + v2*v2 + v3*v3)*0.5d+00

c   Compute local coordinates if localxyflag = 1

        IF (localxyflag .eq. 1) THEN

C direction cosines of the local system . X' is directed parallel
C to the 2-1 side
C Z' is the external normal (counterclockwise). Y' computed as Z' x X'

        rlr = dsqrt(x21**2+y21**2+z21**2)

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
        do 43 i=1,3
          xlcg   = x(i) - xcg
          ylcg   = y(i) - ycg
          zlcg   = z(i) - zcg
          xl(i) = xp(1) * xlcg + xp(2) * ylcg + xp(3) * zlcg
          yl(i) = yp(1) * xlcg + yp(2) * ylcg + yp(3) * zlcg
          zl(i) = zp(1) * xlcg + zp(2) * ylcg + zp(3) * zlcg
  43    continue
     
      ENDIF
C
      return
      end
