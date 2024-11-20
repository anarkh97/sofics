      subroutine     brkcmt(e, nu, c)
      real*8  e, nu, c(6,6)
      integer           i, j
      real*8 lam, mu
      do 20  i=1,6
        do 10  j=1,6
          c(i,j) = 0.0
 10     continue
 20   continue
      lam = nu*e/((1.0+nu)*(1.0-2.0*nu))
      mu  = e/(2.0*(1.0+nu))
      c(1,1) =  lam+2*mu 
      c(2,2) =  c(1,1)
      c(3,3) =  c(1,1)
      c(1,2) =  lam 
      c(2,1) =  c(1,2)
      c(1,3) =  c(1,2)
      c(2,3) =  c(1,2)
      c(3,1) =  c(1,2)
      c(3,2) =  c(1,2)
      c(4,4) =  mu 
      c(5,5) =  c(4,4)
      c(6,6) =  c(4,4)
      return
      end
C=END FORTRAN
