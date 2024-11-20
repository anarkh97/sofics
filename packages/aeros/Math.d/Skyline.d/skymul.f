C=DECK SKYMUL
C=PURPOSE Postmultiply skyline-stored matrix by vector
C=AUTHOR C. A. Felippa, April 1969
C=VERSION March 1982 (Fortran 77)
C=EQUIPMENT Machine independent
C=BLOCK ABSTRACT
C
C     This subroutine postmultiplies skyline-stored matrix [A] by
C     vector [X]. The result [AX] is produced in double precision.
C     Further operations are controlled by input flag IOP :
C       IOP =  0    exit after obtaining [AX]
C       IOP =  1    form residual vector [R] = [B] - [AX] in R
C       IOP = -1    move [X] to R, and form [B] - [AX] in X.
C
C     Arrays AX (d.p.) and X may be equivalent in the calling program
C
C=END ABSTRACT
C=BLOCK FORTRAN
      subroutine SKYMUL (a, n, ld, x, ax, iop, b, r)
C
C                       A R G U M E N T S
C
      double precision        ax(*)
      double precision        a(*),     b(*),     x(*),     r(*)
      integer       iop, n, ld(*)
C
C                       T Y P E   &   D I M E N S I O N
C
      integer           i, j, k, ii, m
      double precision  aij, axi
C
C                        L O G I C
C
      do 2000  i = 1,n
        ii =     iabs(ld(i+1))
        ax(i) =  a(ii)*x(i)
        m =      ii - iabs(ld(i)) - 1
        if (m. eq. 0)                  go to 2000
        do 1500  k = 1,m
          j =      i - k
          aij =    a(ii-k)
          ax(i) =  ax(i) + aij*x(j)
          ax(j) =  ax(j) + aij*x(i)
 1500     continue
 2000   continue
C
      if (iop)                         2500,5000,3500
 2500 do 2800  i = 1,n
        axi = ax(i)
        r(i) = x(i)
        x(i) = b(i) - axi
        if (ld(i+1).le.0)              x(i) = 0.
 2800   continue
      go to 5000
 3500 do 3800  i = 1,n
        r(i) = b(i) - ax(i)
        if (ld(i+1).le.0)              r(i) = 0.
 3800   continue
C
 5000 return
      end
C=END FORTRAN

