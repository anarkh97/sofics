        subroutine outprd(u1,u2,u3,mat,const)
C
        real*8 u1,u2,u3,mat(3,3),const
C
        mat(1,1) = const*u1*u1
        mat(1,2) = const*u1*u2
        mat(1,3) = const*u1*u3
        mat(2,1) = mat(1,2) 
        mat(2,2) = const*u2*u2
        mat(2,3) = const*u2*u3 
        mat(3,1) = mat(1,3)
        mat(3,2) = mat(2,3)
        mat(3,3) = const*u3*u3
C
        return
        end
