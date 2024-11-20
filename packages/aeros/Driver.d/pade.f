C
C********************************************************************************
      SUBROUTINE PADE(c,ic,a,ia,b,ib,l,m,x,w,w1,w2,ik,iw)
C********************************************************************************
C
C     The real function PADE returns the value of the [L/M] Pade
C     approximants evaluated at x using f(x) and the (L+M)-th derivative
C     of f at x. The Pade numerator coefficients are stored in the 
C     first l+1 locations of the array a.
C
C     The denominator coefficients are stored in the first m+1
C     locations of the array b.
C
C     Double precision arrays w(iw,iw), w1(iw), w2(i2) and integer
C     array ik(iw) are required for work space.
C     The strict inequalities iw > 0 and iw.ge.m, ib > m, ia > l,
C     ic > l+m are necessary to fulfill the storage requirements
C     at entry.
C
C     c is a real array of dimension ic. On entry, it contains
C     the given series coefficients and is unchanged on exit.
C
C     a is a real array of dimension ia. On exit, a(1), a(2), ..., a(l+1)
C     contain the numerator parameters a0, a1, ..., a(l), respectively
C
C     b is a real array of dimension ib. On exit, b(1), b(2), ..., b(m+1)
C     contain the numerator parameters b0, b1, ..., bm, respectively.
C
C     l is an integer specifying the degree of the numerator.
C
C     m is an integer specifying the degree of the denominator.
C
C     w is a two-dimensional double precision real array of dimensions (iw,iw)
C     which must be provided for workspace.
C
C     w1, w2 are two double precision real arrays of dimension iw which must
C     be provided for worksapce.
C
C     ik is an integer array of dimension iw which must be provided for workspace.
C
C     ic, ia, ib, iw are four integers specifying various array dimensions, to be
C     set on entry and unchanged on exit
C
C     x is the real variable of the given power series which must be set with some
C     value on entry. x is unchanged on exit.
C
C     If the integers ic, ia, ib, and iw are inconsistently chosen (it is required
C     that l>= 0, m>=0, iw>=m, ib>m, ia>l, and ic>l+m), or the approximant is 
C     degenerate, an error message is printed, and control is returned to the
C     calling program. If the approximant is evaluated at a pole, an error message
C     is printed. The arrays a, b are set correctly, and the artificial value
C     pade( ) = 0 is returned.
C
C********************************************************************************
C
      real*8  c(ic),b(ib),a(ia),w(iw,iw),w1(iw),w2(iw)
      integer ik(iw)
      real*8  one,ought,det,t,rpade
      data one,ought/0.1d 01,0.0d 00/
C
C.....Next statement applies to computers with 12 figure precision
C
      eps = 0.1e-24
      lp1 = l+1
      b(1) = 1.0
      if(l) 40,1,1
1     if(m) 40,6,2
2     if(iw-m) 40,3,3
3     if(ib-m) 40,40,4
4     if(ia-l) 40,40,5
5     if(ic-l-m) 40,40,8
C
6     do 7 i = 1,lp1
7     a(i) = c(i)
C
      d = 1.0
C
      go to 35
C
8     do 13 i = 1,m
C
         do 12 j = 1,m
         if(l+i-j) 11,10,10
10       i1 = l + i - j + 1
C
      w(i,j) = c(i1)
C
      go to 12
C
11    w(i,j) = ought
12    continue
C
      i1 = l + i + 1
      w1(i) = -c(i1)
13    continue
C
C.....Solve denominator equations using Gauss-Jordan elimination
C.....with full pivoting and double precision
C
      det = one
C
      do 14 j = 1,m
14    ik(j) = 0
C
      do 30 i = 1,m
      t = ought
C
         do 19 j = 1,m
         if(ik(j) - 1) 15,19,15
C
15          do 18 k = 1,m
            if(ik(k) - 1) 16,18,30
16          if(dabs(t) - dabs(w(j,k))) 17,17,18
17          irow = j
            icol = k
            t = w(j,k)
18          continue
C
19       continue
C
      ik(icol) = ik(icol) + 1
      if(irow-icol) 20,22,20
20    det = -det
C
         do 21 n = 1,m
         t = w(irow,n)
         w(irow,n) = w(icol,n)
21       w(icol,n) = t
C
      t = w1(irow)
      w1(irow) = w1(icol)
      w1(icol) = t
22    w2(i) = w(icol,icol)
      im1 = i-1
      if(i.eq.1) go to 25
      t = dexp(dlog(dabs(det))/dble(float(im1)))
      if(dabs(w2(i)) - t*dble(float(i))*eps) 23,25,25
23    det = ought
24    format(1x,5hThe [,i4,1H/,i4,
     *45h] Pade approximant apparently is reducible or/
     *49h does not exist. If this is so, try using a lower/
     *19h order approximant.)
      write(6,24) l,m
      return
C
25    det = det*w2(i)
      w(icol,icol) = one
C
         do 26 n = 1,m
26       w(icol,n) = w(icol,n)/w2(i)
C
      w1(icol) = w1(icol)/w2(i)
C
         do 29 li = 1,m
         if(li-icol) 27,29,27
27       t = w(li,icol)
         w(li,icol) = ought
C
            do 28 n = 1,m
28          w(li,n) = w(li,n) - w(icol,n)*t
C
         w1(li) = w1(li) - w1(icol)*t
29       continue
C
30    continue
C
C.....Evaluate denominator
C
      do 31 i = 1,m
31    b(i+1) = w1(i)
C
      d = b(m+1)
C
      do 32 i = 1,m
      i1 = m + 1 -i
32    d = x*d + b(i1)
C
C.....Evaluate numerator
C
      do 34 i = 1,lp1
      k = m + 1
      if(k.gt.i) k = i
      y = 0.0
C
         do 33 j = 1,k
         i1 = i - j + 1
33       y = y + b(j)*c(i1)
C
      a(i) = y
34    continue
C
      y = abs(d)/m
      if(y.lt.eps) go to 42
35    if(l) 40,36,37
36    rpade = c(1)/d
      return
C
37    y = a(l+1)
C
      do 38 i = 1,l
      i1 = l - i + 1
38    y = y*x + a(i1)
C
      rpade = y/d
      return
C
39    format(4h ic=,i4,4h ia=,i4,4h ib=,i4,3h l=,i4,3h m=,i4,4h iw=,
     *i4/ 43h and there is an error among these integers)
40    write(6,39) ic,ia,ib,l,m,iw
      rpade = 0.0
      return
C
41    format(1x,e16.8,40h apparently is a pole of the approximant)
42    write(6,41) x
      rpade = 0.0
      return
      end
