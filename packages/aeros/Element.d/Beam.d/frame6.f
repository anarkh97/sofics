C=DECK  FRAME6
C=BLOCK FORTRAN
        subroutine        frame6(frame,u,v,w)
        real*8  frame(9), u(3), v(3), w(3)
        u(1) = frame(1)
        u(2) = frame(2)
        u(3) = frame(3)
        v(1) = frame(4)
        v(2) = frame(5)
        v(3) = frame(6)
        w(1) = frame(7)
        w(2) = frame(8)
        w(3) = frame(9)
        return
        end
C=END FORTRAN
