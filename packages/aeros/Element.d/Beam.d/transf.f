C=DECK TRANSF
C=BLOCK FORTRAN
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C*                                                           *
C*  Purpose :  Construct Element Stiffness Matrix for Beam   *
C*             element.  It stores only the upper triangular *
C*                                                           *
C*  Authors :  E. Pramono & Charbel Farhat                   *
C*                                                           *
C*            CSSC/ BOULDER 1990                             *
C*                                                           *
C*  Modified By : J. C. Chiou and P. Stern, 1/15/91          *
C*                                                           *
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
        subroutine	transf(u,v,w,tran,trant)
C
        real*8 u(3), v(3), w(3), tran(12,12), trant(12,12)
C
        integer i, j, k, ic
        real*8 t33(3,3)
C
        do 100 j = 1, 12
          do 100 i = 1, 12
            tran(i,j)  = 0.0
            trant(i,j) = 0.0
 100    continue
 
        t33(1,1) = u(1)
        t33(1,2) = u(2)
        t33(1,3) = u(3)
        t33(2,1) = v(1)
        t33(2,2) = v(2)
        t33(2,3) = v(3)
        t33(3,1) = w(1)
        t33(3,2) = w(2)
        t33(3,3) = w(3)
C
        do 200 k = 1, 4
          ic = 3*(k-1)
          do 150 i = 1, 3
            do 150 j = 1, 3
            tran(ic+i,ic+j) = t33(i,j)
 150    continue
 200    continue 
C
        do 300 j = 1, 12
          do 300 i = 1, 12
            trant(i,j) = tran(j,i)
 300    continue
C
        return
        end
C=END FORTRAN
