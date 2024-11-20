C=DECK SM3MHE
C=PURPOSE Form high-order material stiffness of 9-dof EFF triangle
C=AUTHOR C. A. Felippa
C=VERSION June 1991
C=EQUIPMENT Machine independent
C=KEYWORDS finite element
C=KEYWORDS material stiffness matrix
C=KEYWORDS triangle membrane high-order extended free formulation
C=BLOCK ABSTRACT
C
C     SM3ME forms the higher order stiffness matrix of a 9-dof
C     membrane triangle based on the extended free formulation.
C     This implementation has alphah=5/4 hardwired, and is
C     optimized for maximum formation speed.
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       CALL      SM3MHE (X, Y, DM, F, LS, SM, M, STATUS)
C
C     The inputs are:
C
C       X         (3 x 1) array of x coordinates of triangle nodes
C       Y         (3 x 1) array of y coordinates of triangle nodes
C       DM        (3 x 3) matrix constitutive matrix already
C                 integrated through the thickness
C       F         Factor by which all stiffness entries will be multiplied.
C                 It is beta or 0.5*beta
C       SM        Incoming material stiffness array.
C       LS        (9 x 1) array of stiffness location pointers
C                 (see examples in SM3MB).
C                 three rotational DOF will appear at the end.
C       M         First dimension of SM in calling program.
C
C     The outputs are:
C
C       SM        Output stiffness array with higher order stiffness
C                 coefficients added in.
C                 The (i,j)-th entry of the basic element stiffness is added
C                 to SM(K,L), where K=LS(I) and L=LS(J).
C                 (Drilling freedoms are internally 7,8,9)
C
C       STATUS    Status character variable.  Blank if no error
C                 detected.
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine    SM3MHE(x, y, dm, f, ls, sm, m, status)
C
C                   A R G U M E N T S
C
      integer           ls(9), m
      double precision  x(3), y(3), dm(3,3), f, sm(m,m)
      character*(*)     status
C
C                   T Y P E   &   D I M E N S I O N
C
      double precision  x0,y0, x10,x20,x30, y10,y20,y30
      double precision  x12, x21, x23, x32, x31, x13
      double precision  y12, y21, y23, y32, y31, y13
      double precision  aa12,aa23,aa31,ss12,ss23,ss31,ss1,ss2,ss3
      double precision  caa12,caa23,caa31, sum
      double precision  ca,cax10,cax20,cax30,cay10,cay20,cay30
      double precision  area, area2, kfac
      double precision  kqh(6,6),hmt(6,3),hqt(6,3),kth(3,3)
      double precision  s(3),w(6),xyij(6)
      double precision  e11,e22,e33,e12,e13,e23
      integer           i,j,k,l
C
C                   L O G I C
C
      status =   ' '
      if (f .eq. 0.0)         return
      x12 =      x(1) - x(2)
      x21 =     -x12
      x23 =      x(2) - x(3)
      x32 =     -x23
      x31 =      x(3) - x(1)
      x13 =     -x31
      y12 =      y(1) - y(2)
      y21 =     -y12
      y23 =      y(2) - y(3)
      y32 =     -y23
      y31 =      y(3) - y(1)
      y13 =     -y31
      area2 =    x21*y31-x31*y21
      if (area2 .le. 0.0)      then
        status = 'SM3MBE: Negative area'
        if (area2 .eq. 0.0)   status = 'SM3MBE: Zero area'
        return
      end if
      area  =    0.5D0*area2
      x0 =       (x(1)+x(2)+x(3))/3.0
      y0 =       (y(1)+y(2)+y(3))/3.0
      x10 =      x(1) - x0
      x20 =      x(2) - x0
      x30 =      x(3) - x0
      y10 =      y(1) - y0
      y20 =      y(2) - y0
      y30 =      y(3) - y0
      aa12 =     2.25D0*(x30**2+y30**2)
      aa23 =     2.25D0*(x10**2+y10**2)
      aa31 =     2.25D0*(x20**2+y20**2)
      caa12 =    15.D0/(32.*aa12)
      caa23 =    15.D0/(32.*aa23)
      caa31 =    15.D0/(32.*aa31)
      ss12 =     x12**2+y12**2
      ss23 =     x23**2+y23**2
      ss31 =     x31**2+y31**2
      ss1 =      0.25D0*(ss12-ss31)
      ss2 =      0.25D0*(ss23-ss12)
      ss3 =      0.25D0*(ss31-ss23)
      cay10 =    0.1875D0*y10
      cay20 =    0.1875D0*y20
      cay30 =    0.1875D0*y30
      cax10 =    0.1875D0*x10
      cax20 =    0.1875D0*x20
      cax30 =    0.1875D0*x30
      hmt(1,1) = caa12*((-ss3+0.6D0*aa12)*y30+area*x30)
      hmt(1,2) =  3.*cay30 - hmt(1,1)
      hmt(1,3) = cay30
      hmt(2,1) = cay10
      hmt(2,2) = caa23*((-ss1+0.6D0*aa23)*y10+area*x10)
      hmt(2,3) =  3.*cay10 - hmt(2,2)
      hmt(3,1) = caa31*((ss2+0.6D0*aa31)*y20-area*x20)
      hmt(3,2) = cay20
      hmt(3,3) =  3.*cay20 - hmt(3,1)
      hmt(4,1) = caa12*((ss3-0.6D0*aa12)*x30+area*y30)
      hmt(4,2) = -3.*cax30 - hmt(4,1)
      hmt(4,3) = -cax30
      hmt(5,1) = -cax10
      hmt(5,2) = caa23*((ss1-0.6D0*aa23)*x10+area*y10)
      hmt(5,3) = -3.*cax10 - hmt(5,2)
      hmt(6,1) = caa31*((-ss2-0.6D0*aa31)*x20-area*y20)
      hmt(6,2) = -cax20
      hmt(6,3) = -3.*cax20 - hmt(6,1)
      do 2000  j = 1,3
        sum =    (2.D0/9.)*(hmt(1,j)+hmt(2,j)+hmt(3,j))
        hqt(1,j) =  sum - (4.D0/3.)*hmt(1,j)
        hqt(2,j) =  sum - (4.D0/3.)*hmt(2,j)
        hqt(3,j) =  sum - (4.D0/3.)*hmt(3,j)
        sum =    (2.D0/9.)*(hmt(4,j)+hmt(5,j)+hmt(6,j))
        hqt(4,j) =  sum - (4.D0/3.)*hmt(4,j)
        hqt(5,j) =  sum - (4.D0/3.)*hmt(5,j)
        hqt(6,j) =  sum - (4.D0/3.)*hmt(6,j)
 2000   continue
      kfac =     1.5D0*f/area2
      e11 =      kfac * dm(1,1)
      e22 =      kfac * dm(2,2)
      e33 =      kfac * dm(3,3)
      e12 =      kfac * dm(1,2)
      e13 =      kfac * dm(1,3)
      e23 =      kfac * dm(2,3)
      kqh(1,1) = 2*(e11*y30**2-2*e13*x30*y30+e33*x30**2)
      kqh(1,2) = ((e13*x10-e11*y10)*y30+(e13*y10-e33*x10)*x30)
      kqh(1,3) = ((e13*x20-e11*y20)*y30+(e13*y20-e33*x20)*x30)
      kqh(1,4) = 2*(e13*y30**2-(e33+e12)*x30*y30+e23*x30**2)
      kqh(1,5) = ((e12*x10-e13*y10)*y30+(e33*y10-e23*x10)*x30)
      kqh(1,6) = ((e12*x20-e13*y20)*y30+(e33*y20-e23*x20)*x30)
      kqh(2,1) = kqh(1,2)
      kqh(2,2) = 2*(e11*y10**2-2*e13*x10*y10+e33*x10**2)
      kqh(2,3) = ((e13*x10-e11*y10)*y20+(e13*y10-e33*x10)*x20)
      kqh(2,4) = ((e33*x10-e13*y10)*y30+(e12*y10-e23*x10)*x30)
      kqh(2,5) = 2*(e13*y10**2-(e33+e12)*x10*y10+e23*x10**2)
      kqh(2,6) = ((e33*x10-e13*y10)*y20+(e12*y10-e23*x10)*x20)
      kqh(3,1) = kqh(1,3)
      kqh(3,2) = kqh(2,3)
      kqh(3,3) = 2*(e11*y20**2-2*e13*x20*y20+e33*x20**2)
      kqh(3,4) = ((e33*x20-e13*y20)*y30+(e12*y20-e23*x20)*x30)
      kqh(3,5) = ((e12*x10-e13*y10)*y20+(e33*y10-e23*x10)*x20)
      kqh(3,6) = 2*(e13*y20**2-(e33+e12)*x20*y20+e23*x20**2)
      kqh(4,1) = kqh(1,4)
      kqh(4,2) = kqh(2,4)
      kqh(4,3) = kqh(3,4)
      kqh(4,4) = 2*(e33*y30**2-2*e23*x30*y30+e22*x30**2)
      kqh(4,5) = ((e23*x10-e33*y10)*y30+(e23*y10-e22*x10)*x30)
      kqh(4,6) = ((e23*x20-e33*y20)*y30+(e23*y20-e22*x20)*x30)
      kqh(5,1) = kqh(1,5)
      kqh(5,2) = kqh(2,5)
      kqh(5,3) = kqh(3,5)
      kqh(5,4) = kqh(4,5)
      kqh(5,5) = 2*(e33*y10**2-2*e23*x10*y10+e22*x10**2)
      kqh(5,6) = ((e23*x10-e33*y10)*y20+(e23*y10-e22*x10)*x20)
      kqh(6,1) = kqh(1,6)
      kqh(6,2) = kqh(2,6)
      kqh(6,3) = kqh(3,6)
      kqh(6,4) = kqh(4,6)
      kqh(6,5) = kqh(5,6)
      kqh(6,6) = 2*(e33*y20**2-2*e23*x20*y20+e22*x20**2)
      kth(1,1) =  0.0
      kth(1,2) =  0.0
      kth(2,2) =  0.0
      kth(1,3) =  0.0
      kth(2,3) =  0.0
      kth(3,3) =  0.0
      do 3500  j = 1,3
        do 3200  i = 1,6
          w(i) =    kqh(i,1)*hqt(1,j) + kqh(i,2)*hqt(2,j)
     $            + kqh(i,3)*hqt(3,j) + kqh(i,4)*hqt(4,j)
     $            + kqh(i,5)*hqt(5,j) + kqh(i,6)*hqt(6,j)
 3200     continue
        do 3300  i = 1,j
          kth(i,j) = kth(i,j) + hqt(1,i)*w(1) + hqt(2,i)*w(2)
     $                        + hqt(3,i)*w(3) + hqt(4,i)*w(4)
     $                        + hqt(5,i)*w(5) + hqt(6,i)*w(6)
          kth(j,i) = kth(i,j)
 3300     continue
 3500   continue
      s(1) =   kth(1,1) + kth(1,2) + kth(1,3)
      s(2) =   kth(2,1) + kth(2,2) + kth(2,3)
      s(3) =   kth(3,1) + kth(3,2) + kth(3,3)
      ca =       0.25D0/area
      xyij(1) =  ca*x32
      xyij(2) =  ca*y32
      xyij(3) =  ca*x13
      xyij(4) =  ca*y13
      xyij(5) =  ca*x21
      xyij(6) =  ca*y21
      do 4000  j = 1,9
        l =    ls(j)
        do 3600  i = 1,3
          if (j .le. 6) then
             w(i) =  s(i)*xyij(j)
          else
             w(i) =  kth(i,j-6)
          end if
 3600     continue
        sum =    w(1) + w(2) + w(3)
        do 3700  i = 1,j
          k =      ls(i)
          if (i .le. 6) then
             sm(k,l) =  sm(k,l) + sum*xyij(i)
          else
             sm(k,l) =  sm(k,l) + w(i-6)
          end if
          sm(l,k) =  sm(k,l)
 3700     continue
 4000   continue
      return
      end
C=END FORTRAN

