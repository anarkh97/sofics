C=DECK SM3MB
C=PURPOSE Form basic membrane stiffness of 9-dof triangle
C=AUTHOR C. A. Felippa, June 1984
C=VERSION June 1984
C=EQUIPMENT Machine independent
C=KEYWORDS finite element membrane plane stress
C=KEYWORDS basic material stiffness matrix
C=BLOCK ABSTRACT
C
C     SM3MB forms the material element stiffness matrix associated with
C     the basic displacement modes (rigid modes + constant strain
C     modes) of a 9-dof plane-stress triangle based on the
C     free formulation of Bergan and Nygard.
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       CALL      SM3MB (X, Y, DM, ALPHA, F, LS, SM, M, STATUS)
C
C     where the input arguments are
C
C       X         (3 x 1) array of x coordinates of triangle nodes.
C       Y         (3 x 1) array of y coordinates of triangle nodes.
C       DM        (3 x 3) matrix relating in-plane forces to strains.
C       ALPHA     Rotational lumping factor; if zero form CST.
C       F         Factor by which stiffness entries will be multiplied.
C       LS        (9 x 1) array of stiffness location pointers
C                 (see Output SM).
C       SM        Incoming material stiffness array.
C       M         First dimension of SM in calling program.
C
C     The outputs are:
C
C       SM        Output stiffness array with basic stiffness
C                 coefficients added in.  The (i,j)-th entry of the
C                 basic element stiffness is added to SM(K,L),
C                 where K=LS(I) and L=LS(J).
C       STATUS    Status character variable.  Blank if no error
C                 detected.
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine    SM3MB
     $             (x, y, dm, alpha, f, ls, sm, m, status)
C
C                   T Y P E   &   D I M E N S I O N
C
C=BLOCK VAX
C     implicit      none
C=END VAX
      character*(*)  status
      integer       m, ls(9)
C=BLOCK DOUBLE
      double precision   x(3), y(3), dm(3,3), alpha, f, p(9,3), sm(m,m)
      double precision   area2, c
      double precision   d11, d12, d13, d22, d23, d33
      double precision   x21, x32, x13, y21, y32, y13
      double precision   x12, x23, x31, y12, y23, y31
      double precision   s1, s2, s3
C=ELSE
C=END DOUBLE
      integer       i, j, k, l, n
C
C                   L O G I C
C
      status =   ' '
      x21 =      x(2) - x(1)
      x12 =     -x21
      x32 =      x(3) - x(2)
      x23 =     -x32
      x13 =      x(1) - x(3)
      x31 =     -x13
      y21 =      y(2) - y(1)
      y12 =     -y21
      y32 =      y(3) - y(2)
      y23 =     -y32
      y13 =      y(1) - y(3)
      y31 =     -y13
      area2 =    y21*x13 - x21*y13
      if (area2 .le. 0.0)      then
        status = 'SM3MB: Zero area'
        if (area2 .eq. 0.0)   status = 'SM3MB: Zero area'
        return
      end if
      p(1,1) =   y23
      p(2,1) =   0.0
      p(3,1) =   y31
      p(4,1) =   0.0
      p(5,1) =   y12
      p(6,1) =   0.0
      p(1,2) =   0.0
      p(2,2) =   x32
      p(3,2) =   0.0
      p(4,2) =   x13
      p(5,2) =   0.0
      p(6,2) =   x21
      p(1,3) =   x32
      p(2,3) =   y23
      p(3,3) =   x13
      p(4,3) =   y31
      p(5,3) =   x21
      p(6,3) =   y12
      n =        6
      if (alpha .ne. 0.0)         then
	coef1  = alpha/6.0
	coef2  = alpha/3.0
        p(7,1) =  y23*(y13-y21)*coef1
        p(7,2) =  x32*(x31-x12)*coef1
        p(7,3) =  (x31*y13-x12*y21)*coef2
        p(8,1) =  y31*(y21-y32)*coef1
        p(8,2) =  x13*(x12-x23)*coef1
        p(8,3) =  (x12*y21-x23*y32)*coef2
        p(9,1) =  y12*(y32-y13)*coef1
        p(9,2) =  x21*(x23-x31)*coef1
        p(9,3) =  (x23*y32-x31*y13)*coef2
        n = 9
      end if
      c =       0.5*f/area2
      d11 =     c * dm(1,1)
      d22 =     c * dm(2,2)
      d33 =     c * dm(3,3)
      d12 =     c * dm(1,2)
      d13 =     c * dm(1,3)
      d23 =     c * dm(2,3)
      do 3000  j = 1,n
        l =      ls(j)
        s1 =     d11*p(j,1) + d12*p(j,2) + d13*p(j,3)
        s2 =     d12*p(j,1) + d22*p(j,2) + d23*p(j,3)
        s3 =     d13*p(j,1) + d23*p(j,2) + d33*p(j,3)
        do 2500  i = 1,j
          k =      ls(i)
          sm(k,l) =  sm(k,l) + (s1*p(i,1) + s2*p(i,2) + s3*p(i,3))
          sm(l,k) =  sm(k,l)
 2500     continue
 3000   continue
      return
      end
C=END FORTRAN
