C=====================================================================C
C        1         2         3         4         5         6         7C
C234567890123456789012345678901234567890123456789012345678901234567890C
C=====================================================================C
      subroutine compjac( a , n , np , d , v , b , z , nrot )
C=====================================================================C
C                                                                     C
C     WARNING   The Original Matrix [a] is Overwritten!               C
C     -------                                                         C
C                                                                     C
C=====================================================================C
C
C     -----------
C     DECLARATION
C     -----------
C
C.....Global Variables
C
      integer    n , np , nrot
      real*8     a(np,np) , d(np) , v(np,np) , b(np) , z(np)
C
C.....Define the Maximum Number of Jacobi Iterations
C
      integer    maxjac
      parameter( maxjac = 1000 )
C
C.....Local Variables
C
      integer    i , j , k , ip , iq
      real*8     zero , one , sm , tresh , g , h , t
      real*8     p , c , s , tau , theta
C
C     ----
C     DATA
C     ----
C
      data zero  /0.000000D+00/
      data one   /1.000000D+00/
C
C     -----
C     LOGIC
C     -----
C
C.....OUTPUT THE MATRIX FOR CHECK
C
***   open( unit=66 , file="ke.mat" )
***   do 9001 i=1,n
***      write(66,1) (a(i,j),j=1,n)
*9001 continue
*** 1 format(50(1x,E13.6))
***   close(66)
C
C.....CLEAR ALL OUTPUT AND LOCAL STORAGES
C
      do 1001 iq=1,np
         do 1002 ip=1,np
            v(ip,iq) = zero
 1002    continue
 1001 continue
C
      do 1003 ip=1,np
         d(ip) = zero
 1003 continue
C
      do 1004 ip=1,np
         b(ip) = zero
 1004 continue
C
      do 1005 ip=1,np
         z(ip) = zero
 1005 continue
C
      nrot = 0
C
C.....INITIALIZE THE MODE SHAPE MATRIX TO UNITY
C
      do 2001 ip=1,n
         v(ip,ip) = one
 2001 continue
C
C.....INITIALIZE THE STORAGES [D] AND [B] W/ THE MAIN DIAGONAL
C.....OF THE INPUT MATRIX [A]
C
      do 2002 ip=1,n
         d(ip) = a(ip,ip)
         b(ip) = a(ip,ip)
 2002 continue
C
C.....JACOBI ITERATIONS
C
      do 3001 i=1,maxjac
C
      sm = zero
C
      do 3101 ip=1,(n-1)
         do 3102 iq=(ip+1),n
            sm = sm + abs(a(ip,iq))
 3102    continue
 3101 continue
C
      if ( sm.eq.zero ) go to 100
C
      if ( i.lt.4 ) then
         tresh = 0.200D+00*sm/(n*n)
      else
         tresh = zero
      endif
C
      do 3103 ip=1,n-1
      do 3104 iq=ip+1,n
C
         g = 100.00D+00*abs(a(ip,iq))
C
         if (      (i.gt.4)
     $        .and.((abs(d(ip))+g).eq.abs(d(ip)))
     $        .and.((abs(d(iq))+g).eq.abs(d(iq))) ) then
C
            a(ip,iq) = zero
C
         else
C
            if ( abs(a(ip,iq)).gt.tresh ) then
C
               h = d(iq) - d(ip)
C
               if ( (abs(h)+g).eq.abs(h) ) then
                  t = a(ip,iq)/h
               else
                  theta = 0.500D+00*h/a(ip,iq)
                  t     = one/(abs(theta)+sqrt(one+(theta*theta)))
                  if ( theta.lt.zero ) t = -t
               endif
C
               c        = one/sqrt(one+(t*t))
               s        = t*c
               tau      = s/(one+c)
               h        = t*a(ip,iq)
               z(ip)    = z(ip) - h
               z(iq)    = z(iq) + h
               d(ip)    = d(ip) - h
               d(iq)    = d(iq) + h
               a(ip,iq) = zero
C
               do 3105 j=1,(ip-1)
                  g       = a(j,ip)
                  h       = a(j,iq)
                  a(j,ip) = g - s*(h+(g*tau))
                  a(j,iq) = h + s*(g-(h*tau))
 3105          continue
C
               do 3106 j=(ip+1),(iq-1)
                  g       = a(ip,j)
                  h       = a(j,iq)
                  a(ip,j) = g - s*(h+(g*tau))
                  a(j,iq) = h + s*(g-(h*tau))
 3106          continue
C
               do 3107 j=(iq+1),n
                  g       = a(ip,j)
                  h       = a(iq,j)
                  a(ip,j) = g - s*(h+(g*tau))
                  a(iq,j) = h + s*(g-(h*tau))
 3107          continue
C
               do 3108 j=1,n
                  g       = v(j,ip)
                  h       = v(j,iq)
                  v(j,ip) = g - s*(h+(g*tau))
                  v(j,iq) = h + s*(g-(h*tau))
 3108          continue
C
               nrot = nrot + 1
C
            endif
C
         endif
C
 3104 continue
 3103 continue
C
      do 3109 ip=1,n
          b(ip) = b(ip) + z(ip)
          d(ip) = b(ip)
          z(ip) = zero
 3109 continue
C
C.....END OF JACOBI ITERATION
C
 3001 continue
C
C.....MAXIMUM NUMBER OF JACOBI ITERATIONS EXCEEDED => ERROR
C
      go to 200
C
C.....SORT THE EIGENVALUES AND EIGENVECTORS
C
  100 continue
C
      do 4001 i=1,(n-1)
C
         k = i
         p = d(i)
C
         do 4002 j=(i+1),n
            if ( d(j).le.p ) then
               k = j
               p = d(j)
            endif
 4002    continue
C
         if ( k.ne.i ) then
            d(k) = d(i)
            d(i) = p
            do 4003 j=1,n
               p      = v(j,i)
               v(j,i) = v(j,k)
               v(j,k) = p
 4003       continue
         endif
C
 4001 continue
C
C     ------
C     RETURN
C     ------
C
      return
C
C     ---------------
C     ERROR TREATMENT
C     ---------------
C
  200 continue
      write(*,*) "*** FATAL ERROR in Routine COMPJAC       ***"
      write(*,*) "*** Reached the Maximum Number of Jacobi ***"
      write(*,*) "*** Iterations: Should Never Happen ...  ***"
      write(*,*) "*** Boost Local Parameter [maxjac] And   ***"
      write(*,*) "*** Check the Input Matrix (Symmetric?)  ***"
      stop
C
C     ---
C     END
C     ---
C
      end
C=end of routine "COMPJAC"
C========================C
