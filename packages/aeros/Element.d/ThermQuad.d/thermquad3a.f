        subroutine thermquad3a(x, y, z, xl, yl)

******************************************************************
* THIS  SUBROUTINE WILL COMPUTE THE CONDUCTION AND COVECTION COM-*
* PONENTS OF THE HEAT ELEMENT CONDUCTION MATRICIES.              *
* FOR THE THREE DIMENSIONAL HEAT ELEMENT. A 3-D QUADRILATERAL    *
* / HEAT SHELL WITH 4 NODES.                                     *
******************************************************************
*
*
*           DECLARE GLOBAL VARIABLES
*
*        X,Y,Z = ELEMENT COORDINATES
* 
*            RETURNED VARIABLES :
*
*        XL, YL = ELEMENT LOCAL COORDINATES
*
******************************************************************
C
C.... DECLARE ALL GLOBAL VARIABLES
C
C.... REAL ARRAYS
C
        real*8 x(*), y(*), z(*)
C
C.... LOCAL VARIABLES
C
        real*8 sr
        real*8 xl(4),yl(4),op(3)
        real*8 mx(4),my(4),mz(4),o(4),p(3),q(3),a(3,3)
C
C.... COMPUTE THE MIDPONTS OF THE QUADRILATERAL HEAT SHELL SIDES
C
        do 10 i=1,3
          mx(i) = (x(i+1) + x(i)) / 2.0
          my(i) = (y(i+1) + y(i)) / 2.0
          mz(i) = (z(i+1) + z(i)) / 2.0
10      continue
C
        mx(4) = (x(4) + x(1)) / 2.0
        my(4) = (y(4) + y(1)) / 2.0
        mz(4) = (z(4) + z(1)) / 2.0
C
C.... COMPUTE THE MEDIANS P & Q
C
        p(1) = mx(2) - mx(4)
        p(2) = my(2) - my(4)
        p(3) = mz(2) - mz(4)
C
        q(1) = mx(3) - mx(1)
        q(2) = my(3) - my(1)
        q(3) = mz(3) - mz(1)
C
C.... COMPUTE THE CROSS-PRODUCT OF P AND Q TO GET W(normal vector)
C
        call veccrs(p,q,a(1,3))
C
C.... COMPUTE THE MIDDLE OF THE ELEMENT
C
        o(1) = ( x(1) + x(2) + x(3) + x(4) ) /4.0
        o(2) = ( y(1) + y(2) + y(3) + y(4) ) /4.0
        o(3) = ( z(1) + z(2) + z(3) + z(4) ) /4.0
C
C.... COMPUTE THE U COORDINATE(parallel to horizonal)
C
        a(1,1) = mx(2) - o(1)
        a(2,1) = my(2) - o(2)
        a(3,1) = mz(2) - o(3)
C
C.... NORMALIZE U COORDINATE
C
        sr = dsqrt( a(1,1)**2 + a(2,1)**2 + a(3,1)**2 )
C
        do 20 i=1,3
          a(i,1) = a(i,1) / sr
20      continue
C
C.... COMPUTE THE V COORDINATE
C
        call veccrs(a(1,3),a(1,1),a(1,2))
C
C.... COMPUTE THE LOCAL COORDINATES
C
        do 30 i=1,4
C
C.... FIND VECTOR FROM CENTER TO EACH NODE
C
          op(1) = x(i) - o(1)
          op(2) = y(i) - o(2)
          op(3) = z(i) - o(3)
C
C.... COMPUTE THE DOT PRODUCT
C
          xl(i) = ( op(1)*a(1,1) + op(2)*a(2,1) + op(3)*a(3,1) )
          yl(i) = ( op(1)*a(1,2) + op(2)*a(2,2) + op(3)*a(3,2) )
C
30      continue
C

        return
        end
