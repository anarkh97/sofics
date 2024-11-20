C=PURPOSE Tester to exercise 20-node hexahedron element stiffness routines
C=AUTHOR C. A. Felippa, June 1966
C=VERSION May 1980 (Fortran 77)
C=USAGE Self-contained data. Just execute 'as is'
C=BLOCK DISCUSSION
C
C     TESTHEXA20 is a driver that exercises the stiffness routines
C     of the iso-P 20-node hexahedron finite element for solid analysis.
C
C     The stiffness matrix formation subroutine tested is HEXA20STIF.
C     Subordinate to these are the shape function evaluation routine
C     HEXA20SHAPE and the Gauss-quadrature-information routine QUADGAUSSQ.
C
C     The element tested by TESTHEXA20 is a 20-node hexahedron.  Local
C     node numbering is shown in the Figure below.  Corners are numbered
C     1-8 and midpoints 9-20:
C
C                                   19
C                       8 o+++++++++o+++++++++++o 7
C                       +                     + +
C                  20 o                  18 o   +
C                   +       17            +     o15
C               5 o++++++++++o++++++++++o 6     +
C                 +                     +       +
C                 +                     +       o 3
C              13 o                     o14   +
C                 +                     +   o10
C                 +         9           + +
C               1 o+++++++++o+++++++++++o 2
C
C      ^z
C      |
C      |   /y
C      | /
C      o------> x
C
C     Element geometry may be arbitrary, but usually the tests are carried
C     out on cubes or parallelepipeds to help "eyeball" verification and
C     analytical verification.   If the element geometry is cubic-like,
C     [K] should be identical for Gauss quadrature rules 3x3x3 and 4x4x4
C
C     Besides node coordinates (arrays xe,ye,ze) other test data includes
C
C      e,nu         elastic modulus and Poisson's ratio for isotropic
C                   material (to test for more general materials would
C                   require more input data)
C      pmin,pmax    Gauss integration rule range to investigate
C
C     The test data is preset via DATA statements.
C
C     The test performed include:
C
C     1)  forming and printing the element stiffness matrices.
C     2)  checking symmetry.
C     3)  computing and printing forces associated with the
C         three rigid body motions (should be zero within roundoff)
C     4)  spectral test: all eigenvalues positive except for zero
C         eigenvalues (within roundoff) for rigid-body and possibly
C         zero-energy modes.  If the integration order is 2 or higher
C         there should be only 6 zero eigenvalues associated with the
C         rigid body modes if the stress-strain matrix has full rank.
C
C     A formation time evaluation may be enabled although of course
C     this does not test code correctness; only efficiency.
C
C=END DISCUSSION
C=BLOCK TEST
C     For a rectangular paralellepided element dimensioned 5 x 1 x 3 along
C     the x,y,z dimensions, elastic modulus E=90 and Poisson's ratio 0.20,
C     the stiffness matrix for p=2,3,4 should be (to 5-7 digits):
C
C=END TEST
C=BLOCK FORTRAN
      program       TESTHEXA20
C
C                   T Y P E   A N D   D I M E N S I O N
C
      double precision  e, nu, c(6,6), xe(20), ye(20), ze(20)
      double precision  xc(8), yc(8), zc(8)
C      double precision  esm(60,60), ve(60,6), pe(60,6)
      double precision  esm(60,60)
      double precision  eval(60), evec(60,60)
      integer           p, pmin, pmax
      integer           i, j, m, mij(2,12)
      character         status*60
C
C                   D A T A
C
      data          xc  /1.0,6.0,6.0,1.0,1.0,6.0,6.0,1.0/
      data          yc  /1.0,1.0,2.0,2.0,1.0,1.0,2.0,2.0/
      data          zc  /1.0,1.0,1.0,1.0,4.0,4.0,4.0,4.0/
*     data          xc/ -1.,1.,1.,-1., -1.,1.,1.,-1./
*     data          yc/ -2.,-2, 2.,2., -2.,-2.,2.,2./
*     data          zc/ -3.,-3.,-3.,-3., 3.,3.,3.,3./
      data          mij/ 1,2, 2,3, 3,4, 4,1, 
     $                   1,5, 2,6, 3,7, 4,8,
     $                   5,6, 6,7, 7,8, 8,5 /
      data          e   /9.d0/,  nu/0.2d0/
*     data          pmin,pmax /2,2/
      data          pmin,pmax /1,4/
C
C                   L O G I C
C
      do 1500  i = 1,8
        xe(i) =  xc(i)
        ye(i) =  yc(i)
        ze(i) =  zc(i)
 1500   continue
      do 1800  m = 1,12
        i  = mij(1,m)
        j  = mij(2,m)
        xe(m+8) = (xc(i)+xc(j))/2.
        ye(m+8) = (yc(i)+yc(j))/2.
        ze(m+8) = (zc(i)+zc(j))/2.
 1800   continue
      print'(/A)',' =================================================='
      print'( A)',' TESTHEXA20: Test 20-node hexahedron elem stiffness'
      print'(A/)',' =================================================='
      print'('' Corner x :'',8F8.4)', xc
      print'('' Corner y :'',8F8.4)', yc
      print'('' Corner z :'',8F8.4)', zc
      call       GETCMT3D (e, nu, c)
      print'('' Stress-strain matrix [C]:'')'
      call       MATRIXPRINT8 (c, 6, 6)
C
C                  Cycle over integration rules
C
      do 2000  p = pmin,pmax
        print '(/A)',    ' ===================================='
        print '( A,I4)', ' Testing for integration rule p =',p
        print '( A)',    ' ===================================='
        call   HEXA20STIF (' ', xe, ye, ze, c, p, esm, 60, status)
        if (status(1:1) .ne. ' ')   then
          print *, status
        else
          print '(/A)', ' Element stiffness matrix [K]:'
*         call   MATRIXPRINT8 (esm, 60, 60)
          call   TSYMM  (esm, 60)
*         call   TRIGID3D (xe, ye, ze, 20, esm, ve, pe, 60)
          call   TEIGEN (esm, 60, eval, evec)
C
C                  Next subroutine is machine dependent (in Fortran 77)
C                  because of system timer. Normally kept disabled
C
*         call   TIMESM (xe, ye, ze, 20, c, p, esm, 60)
        end if
 2000   continue
C
      stop
      end
C=END FORTRAN
C=DECK GETCMT3D
C=PURPOSE Form 3D constitutive matrix for an isotropic material
C=BLOCK FORTRAN
      subroutine     GETCMT3D (e, nu, c)
      double precision  e, nu, c(6,6)
      integer           i, j
      do 2000  i = 1,6
        do 1000  j = 1,6
          c(i,j) = 0.
 1000     continue
 2000   continue
      c(1,1) =  e * (1.-nu)/((1.-2.*nu)*(1.+nu))
      c(2,2) =  c(1,1)
      c(3,3) =  c(1,1)
      c(1,2) =  c(1,1)*nu/(1.-nu)
      c(2,1) =  c(1,2)
      c(1,3) =  c(1,2)
      c(2,3) =  c(1,2)
      c(3,1) =  c(1,2)
      c(3,2) =  c(1,2)
      c(4,4) =  e/(2.*(1.+nu))
      c(5,5) =  c(4,4)
      c(6,6) =  c(4,4)
      return
      end
C=END FORTRAN
C=DECK MATRIXPRINT6
C=PURPOSE Print (m x n) DP matrix as 2D array, 6 items/line
C=BLOCK FORTRAN
      subroutine   MATRIXPRINT6 (a, m, n)
      integer           m, n, i, j, jref
      double precision  a(m,n)
C
      do 2000  jref = 0,n-1,6
        print '(2X,6I12)',(j,j=jref+1,min(jref+6,n))
        do 1500  i = 1,m
          print '(I5,6F12.4)',i,(a(i,j),j=jref+1,min(jref+6,n))
 1500     continue
 2000   continue
C
      return
      end
C=END FORTRAN
C=DECK MATRIXPRINT8
C=PURPOSE Print (m x n) DP matrix as 2D array, 8 items/line
C=BLOCK FORTRAN
      subroutine   MATRIXPRINT8 (a, m, n)
      integer           m, n, i, j, jref
      double precision  a(m,n)
C
      do 2000  jref = 0,n-1,8
        print '(2X,8I9)',(j,j=jref+1,min(jref+8,n))
        do 1500  i = 1,m
          print '(I5,8F9.4)',i,(a(i,j),j=jref+1,min(jref+8,n))
 1500     continue
 2000   continue
C
      return
      end
C=END FORTRAN
C=DECK MATRIXPRODUCT
C=PURPOSE  Multiply two matrices
C=BLOCK FORTRAN
      subroutine   MATRIXPRODUCT (a, b, c, l, m, n)
      integer      l, m, n
      double precision a(l,m), b(m,n), c(l,n)
      integer      i, j, k
C
      do 2000  j = 1,n
        do 1800  i = 1,l
          c(i,j) = 0.0
          do 1600  k = 1,m
            c(i,j) = c(i,j) + a(i,k)*b(k,j)
 1600       continue
 1800     continue
 2000   continue
C
      return
      end
C=END FORTRAN
C=DECK RIGID3D
C=PURPOSE Construct rigid-body displacement modes of n-node 3D element
C=BLOCK FORTRAN
      subroutine  RIGID3D  (which, n, x, y, z, v)
      character*2       which
      integer           i, n
      double precision  x(*), y(*), z(*), v(*)
C
      do 1200  i = 1,3*n
        v(i) =  0.0
 1200   continue
C
      if (which .eq. 'TX')        then
        do 2000  i = 1,n
          v(3*i-2) = 1.0
 2000     continue
      else if (which .eq. 'TY')   then
        do 2200  i = 1,n
          v(3*i-1) = 1.0
 2200     continue
      else if (which .eq. 'TZ')   then
        do 2400  i = 1,n
          v(3*i) =   1.0
 2400     continue
      else if (which .eq. 'RX')   then
        do 2600  i = 1,n
          v(3*i-1) =   z(i)
          v(3*i) =    -y(i)
 2600     continue
      else if (which .eq. 'RY')   then
        do 2800  i = 1,n
          v(3*i) =     x(i)
          v(3*i-2) =  -z(i)
 2800     continue
      else if (which .eq. 'RZ')   then
        do 3000  i = 1,n
          v(3*i-2) =   y(i)
          v(3*i-1) =  -x(i)
 3000     continue
      else
        print '('' *RIGID* Illegal mode identifier '',A)', which
      end if
      return
      end
C=END FORTRAN
C=DECK TEIGEN
C=PURPOSE Test the eigenspectrum of the element stiffness matrix
C=BLOCK FORTRAN
      subroutine TEIGEN  (esm, m, eval, evec)
      double precision    esm(m,m), eval(*), evec(m,m), aux(100)
      integer      m
      call   CFJACOBI (esm, m, m, eval, .false., evec, aux)
      print '(/'' Computed stiffness matrix eigenvalues:'')'
      call   MATRIXPRINT8 (eval, 1, m)
*     print '(/'' Computed stiffness matrix eigenvectors:'')'
*     call   MATRIXPRINT8 (evec, m, m)
      return
      end
C=END FORTRAN
C=DECK TIMESM
C=PURPOSE Test element stiffness formation time
C=BLOCK FORTRAN
*     subroutine TIMESM  (xe, ye, ze, n, c, p, esm, m)
*     double precision  xe(n), ye(n), h(n), c(3,3), esm(m,m)
*     character    status*60
*     integer      n, m, p
*     integer      ntimes, ne
*     real         tform, tbeg, tend
C
*     ntimes =  100
*     call   SECNDS (tbeg)
*     do 3000  ne = 1,ntimes
*       call  HEXA20STIF (' ', xe, ye, ze, c, p, esm, m, status)
*3000   continue
*     call   SECNDS (tend)
*     tform =   1000.*(tend-tbeg)/ntimes
*     print '(/'' Element formation time (Sun/Unix):'',F10.2,
*    $         '' milliseconds'')', tform
*     return
*     end
C=END FORTRAN
C=DECK TRIGID
C=PURPOSE Test whether rigid body modes are correctly represented
C=BLOCK FORTRAN
      subroutine   TRIGID3D  (x, y, z, n, esm, vr, pr, m)
      integer           n, m
      double precision  x(*), y(*), z(*), esm(m,*), vr(m,*), pr(m,*)
C
      call    RIGID3D ('TX', n, x, y, z, vr(1,1))
      call    RIGID3D ('TY', n, x, y, z, vr(1,2))
      call    RIGID3D ('TZ', n, x, y, z, vr(1,3))
      call    RIGID3D ('RX', n, x, y, z, vr(1,4))
      call    RIGID3D ('RY', n, x, y, z, vr(1,5))
      call    RIGID3D ('RZ', n, x, y, z, vr(1,6))
      call    MATRIXPRODUCT (esm, vr, pr, m, m, 6)
      print '(/A)',
     $ ' Rigid body modes [vr] (x/y-translations, z-rotation):'
      call    MATRIXPRINT8 (vr, m, 6)
      print '(/A)', ' Forces [K][vr] (should vanish):'
      call    MATRIXPRINT8 (pr, m, 6)
      return
      end
C=END FORTRAN
C=DECK TSYMM
C=PURPOSE Test the symmetry of the element stiffness matrix
C=BLOCK FORTRAN
      subroutine   TSYMM  (esm, m)
      double precision  esm(m,m)
      integer           i, j, m
      double precision  snorm, sum, tol
      data              tol /1.0D-12/
      sum =    0.0
      snorm =  0.0
      do 1500  i = 1,m
        do 1500  j = 1,m
          sum =    sum + abs(esm(i,j)-esm(j,i))
          snorm =  snorm + esm(i,j)**2
 1500   continue
      if (sum .gt. tol*sqrt(snorm))      then
        print *, 'Element stiffness matrix is not symmetric!'
        print *, 'Mean unsymmetry detected:',sum/(m*m)
      else
        print '(/ '' Symmetry verified'')'
      end if
      return
      end
C=END FORTRAN
C=DECK SECNDS
C=PURPOSE Simulation of VAX/VMS's timer on Unix
C=BLOCK FORTRAN
      subroutine   SECNDS (tim)
      real         t(3), tim
*     call etime   (t)
      tim =        t(1)
      return
      end
C=END FORTRAN
C=DECK ERROR
C 
C  =====================  Element routines follow =======================
C

C=END FORTRAN
C=DECK HEXA20SIGM
C=PURPOSE Compute corner node stresses of 20-node iso-P hexahedrom
C=AUTHOR C. A. Felippa, August 1966
C=VERSION  May 1980 (Fortran 77)
C=EQUIPMENT Machine independent
C=KEYWORDS twenty node hexahedron
C=KEYWORDS finite element strains stresses
C=BLOCK ABSTRACT
C
C    Given the node displacement, HEXA20SIGM computes
C    corner stresses on a eight-node hexahedron.
C
C=END ABSTRACT
C=BLOCK usage
C
C   The  calling sequence is
C
C   CALL  HEXA20SIGM (ESCM, X, Y, Z, U, ESIG, STATUS)
C
C   where the input arguments are
C
C      ESCM   Character string defining stress computation method.
C             DIRECT  direct stress evalution at corners.
C             EXTRAP  extrapolate from 2x2x2 Gauss points.
C      x      (20x1)  array of x co-ordinates of hexahedron nodes.
C      y      (20x1)  array of y co-ordinates of hexahedron nodes.
C      z      (20x1)  array of z co-ordinates of hexahedron nodes.
C      c      (6x6)  stress-strain constitutive matrix.
C      u      (60x1) array of element node displacement arranged in
C             ux1,uy1,uz1,ux2,uy2,uz2, ............., ux20,uy20,uz20
C
C   The outputs are
C
C
C      ESIG   (6x8) array of corner node stresses arranged
C             sigxx1, sigyy1, sigzz1, tauxy1, tauyz1, tauxz1,
C             sigxx2, sigyy2, sigzz2, tauxy2, tauyz2, tauxz2,
C              ..    ..   ..    ..   .. ...    ...   ...  ..
C              ..    ..   ..    ..   ..       tauyz8, tauxz8
C
C      STATUS  Status character variable. Blank if no error detected.
C=END USAGE
C=BLOCK FORTRAN
          subroutine    HEXA20SIGM
     $                  (escm, x, y, z, c, u, esig, status )
C
C                A R G U M E N T S
C
         character*(*)  escm, status
         double precision  x(20), y(20), z(20), c(6,6), u(60)
         double precision  esig(6,8)
C
C                L O C A L  V A R I A B L E S
C
         double precision  q(20), qx(20), qy(20), qz(20)
         double precision  xicorn(8), etacorn(8), mucorn(8)
         double precision  xi, eta, mu, det, epsxx, epsyy
         double precision  epszz, gamxy, gamyz, gamxz
         integer           i, n
C
C                D A T A
C
         data    xicorn  /-1.0,  1.0,  1.0, -1.0, -1.0,  1.0, 1.0, -1.0/
         data    etacorn /-1.0, -1.0,  1.0,  1.0, -1.0, -1.0, 1.0,  1.0/
         data    mucorn  /-1.0, -1.0, -1.0, -1.0,  1.0,  1.0, 1.0,  1.0/
C
C                L O G I C
C
         status =  ' '
C
         do  2000  n = 1,8
            xi  = xicorn(n)
            eta = etacorn(n)
            mu = mucorn(n)
            call   HEXA20SHAPE (' ', xi, eta, mu, x, y, z, 
     $                          q, qx, qy, qz, det)
         if  ( det .le. 0.0) then
            status = 'Negative  Jacobian determinant'
           if  ( det .eq. 0.0) then
              status = 'Zero Jacobian determinant'
           endif
          return
         endif
         epsxx = 0.0
         epsyy = 0.0
         epszz = 0.0
         gamxy = 0.0
         gamyz = 0.0
         gamxz = 0.0
         do 1500 i = 1,20
            epsxx =   epsxx + qx(i)*u(3*i-2)
            epsyy =   epsyy + qy(i)*u(3*i-1) 
            epszz =   epszz + qz(i)*u(3*i  )  
            gamxy =   gamxy + qy(i)*u(3*i-2) + qx(i)*u(3*i-1)
            gamyz =   gamyz + qz(i)*u(3*i-1) + qy(i)*u(3*i  )
            gamxz =   gamxz + qx(i)*u(3*i  ) + qz(i)*u(3*i-2)
 1500       continue
         esig(1,n) = c(1,1)*epsxx + c(1,2)*epsyy + c(1,3)*epszz
     $             + c(1,4)*gamxy + c(1,5)*gamyz + c(1,6)*gamxz
         esig(2,n) = c(2,1)*epsxx + c(2,2)*epsyy + c(2,3)*epszz
     $             + c(2,4)*gamxy + c(2,5)*gamyz + c(2,6)*gamxz
         esig(3,n) = c(3,1)*epsxx + c(3,2)*epsyy + c(3,3)*epszz
     $             + c(3,4)*gamxy + c(3,5)*gamyz + c(3,6)*gamxz
         esig(4,n) = c(4,1)*epsxx + c(4,2)*epsyy + c(4,3)*epszz
     $             + c(4,4)*gamxy + c(4,5)*gamyz + c(4,6)*gamxz
         esig(5,n) = c(5,1)*epsxx + c(5,2)*epsyy + c(5,3)*epszz
     $             + c(5,4)*gamxy + c(5,5)*gamyz + c(5,6)*gamxz
         esig(6,n) = c(6,1)*epsxx + c(6,2)*epsyy + c(6,3)*epszz
     $             + c(6,4)*gamxy + c(6,5)*gamyz + c(6,6)*gamxz
 2000    continue
      return
      end
C
C=END FORTRAN
C=DECK HEXA20STIF
C=PURPOSE Form stiffness of 20-node iso-P hexahedron
C=AUTHOR C. A. Felippa, April 1966 (Fortran IV)
C=VERSION May 1980 (Fortran 77)
C=EQUIPMENT Machine independent
C=KEYWORDS twenty node hexahedron
C=KEYWORDS finite element stiffness matrix
C=BLOCK ABSTRACT
C
C     HEXA20STIF forms the element stiffness matrix of a
C     twenty-node hexahedron in 3-D stress.
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       CALL   HEXA20STIF (OPT, X, Y, Z, C, P, SM, NS, STATUS)
C
C     where the input arguments are
C
C       OPT       Option letter argument, presently ignored.
C       X         (20 x 1) array of x coordinates of hexahedron nodes
C       Y         (20 x 1) array of y coordinates of hexahedron nodes
C       Z         (20 x 1) array of z coordinates of hexahedron nodes
C       C         (6 x 6) constitutive material matrix 
C       P         Gauss quadrature rule (no. of points in each dir)  Usually p=3.
C       NS        First dimension of SM in calling program.
C
C     The outputs are:
C
C       SM        (60 x 60) computed element stiffness matrix.  The
C                 arrangement of rows and columns pertains to node
C                 displacements arranged in the order
C                  (ux1, uy1, uz1, ux2, ... uz20)
C                 This particular ordering is obtained by setting array
C                 LS as shown by the DATA statement below
C
C       STATUS    Status character variable.  Blank if no error
C                 detected.
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine    HEXA20STIF
     $          (opt, x, y, z, c, p, sm, ns, status)
C
C                   A R G U M E N T S
C
      character*(*)     opt, status
      integer           p, ns
      double precision  x(20), y(20), z(20), c(6,6)
      double precision  sm(ns,ns)
C
C                   L O C A L   V A R I A B L E S
C
      double precision  q(20), qx(20), qy(20),qz(20)
      double precision  xi, eta, mu, det, w, weight
      double precision  c1x, c1y, c1z, c2x, c2y, c2z, c3x, c3y, c3z
      double precision  c4x, c4y, c4z, c5x, c5y, c5z, c6x, c6y, c6z
      integer           i, ix, iy, iz, j, jx, jy, jz, k, l, m
      integer           ls(60)
C
C                   D A T A
C
      data              ls 
     $    /1,4,7,10,13,16,19,22,25,28,31,34,37,40,43,46,49,52,55,58,
     $     2,5,8,11,14,17,20,23,26,29,32,35,38,41,44,47,50,53,56,59,
     $     3,6,9,12,15,18,21,24,27,30,33,36,39,42,45,48,51,54,57,60 /
C
C                   L O G I C
C
      status =   ' '
      do 1200  j = 1,60
        do 1100  i = 1,60
          sm(i,j) = 0.0
 1100     continue
 1200   continue
C
      do 3000  k = 1,p
        do 2500  l = 1,p
          do 2400  m = 1,p
            call     HEXAGAUSSQ (p, k, p, l, p, m, xi, eta, mu,
     $                           weight)
            call     HEXA20SHAPE (' ', xi, eta, mu, x, y, z, q, qx,
     $                           qy, qz, det)
            if (det .le. 0.0)        then
              status = 'Negative Jacobian determinant'
              if (det .eq. 0.0)      then
                status = 'Zero Jacobian determinant'
              end if
              return
            end if
            w =    weight * det 
C
            do 2000  j = 1,20
              jx =    ls(j)
              jy =    ls(j+20)
              jz =    ls(j+40)
              c1x = (c(1,1)*qx(j) + c(1,4)*qy(j) + c(1,6)*qz(j))*w
              c1y = (c(1,4)*qx(j) + c(1,2)*qy(j) + c(1,5)*qz(j))*w
              c1z = (c(1,6)*qx(j) + c(1,5)*qy(j) + c(1,3)*qz(j))*w
              c2x = (c(2,1)*qx(j) + c(2,4)*qy(j) + c(2,6)*qz(j))*w
              c2y = (c(2,4)*qx(j) + c(2,2)*qy(j) + c(2,5)*qz(j))*w
              c2z = (c(2,6)*qx(j) + c(2,5)*qy(j) + c(2,3)*qz(j))*w
              c3x = (c(3,1)*qx(j) + c(3,4)*qy(j) + c(3,6)*qz(j))*w
              c3y = (c(3,4)*qx(j) + c(3,2)*qy(j) + c(3,5)*qz(j))*w
              c3z = (c(3,6)*qx(j) + c(3,5)*qy(j) + c(3,3)*qz(j))*w
              c4x = (c(4,1)*qx(j) + c(4,4)*qy(j) + c(4,6)*qz(j))*w
              c4y = (c(4,4)*qx(j) + c(4,2)*qy(j) + c(4,5)*qz(j))*w
              c4z = (c(4,6)*qx(j) + c(4,5)*qy(j) + c(4,3)*qz(j))*w
              c5x = (c(5,1)*qx(j) + c(5,4)*qy(j) + c(5,6)*qz(j))*w
              c5y = (c(5,4)*qx(j) + c(5,2)*qy(j) + c(5,5)*qz(j))*w
              c5z = (c(5,6)*qx(j) + c(5,5)*qy(j) + c(5,3)*qz(j))*w
              c6x = (c(6,1)*qx(j) + c(6,4)*qy(j) + c(6,6)*qz(j))*w
              c6y = (c(6,4)*qx(j) + c(6,2)*qy(j) + c(6,5)*qz(j))*w
              c6z = (c(6,6)*qx(j) + c(6,5)*qy(j) + c(6,3)*qz(j))*w
              do 1500  i = j,20
                ix =     ls(i)
                iy =     ls(i+20)
                iz =     ls(i+40)
                sm(ix,jx) = sm(ix,jx)+qx(i)*c1x + qy(i)*c4x + qz(i)*c6x
                sm(jx,ix) = sm(ix,jx)
                sm(iy,jy) = sm(iy,jy)+qx(i)*c4y + qy(i)*c2y + qz(i)*c5y
                sm(jy,iy) = sm(iy,jy)
                sm(iz,jz) = sm(iz,jz)+qx(i)*c6z + qy(i)*c5z + qz(i)*c3z
                sm(jz,iz) = sm(iz,jz)
                sm(ix,jy) = sm(ix,jy)+qx(i)*c1y + qy(i)*c4y + qz(i)*c6y
                sm(iy,jx) = sm(iy,jx)+qx(i)*c4x + qy(i)*c2x + qz(i)*c5x
                sm(jy,ix) = sm(ix,jy)
                sm(jx,iy) = sm(iy,jx)
                sm(ix,jz) = sm(ix,jz)+qx(i)*c1z + qy(i)*c4z + qz(i)*c6z
                sm(iz,jx) = sm(iz,jx)+qx(i)*c6x + qy(i)*c5x + qz(i)*c3x
                sm(jz,ix) = sm(ix,jz)
                sm(jx,iz) = sm(iz,jx)
                sm(iy,jz) = sm(iy,jz)+qx(i)*c4z + qy(i)*c2z + qz(i)*c5z
                sm(iz,jy) = sm(iz,jy)+qx(i)*c6y + qy(i)*c5y + qz(i)*c3y
                sm(jz,iy) = sm(iy,jz)
                sm(jy,iz) = sm(iz,jy)
 1500           continue
 2000         continue
 2400       continue
 2500     continue
 3000   continue
C
      return
      end
C=END FORTRAN
C=DECK HEXA20SHAPE
C=PURPOSE Compute shape functions and xyz derivatives of 20-node hex
C=AUTHOR C. A. Felippa, April 1966 (Fortran IV)
C=VERSION May 1999 - shape functions provided by Kurt Maute
C=EQUIPMENT Machine independent
C=KEYWORDS twenty node hexahedron shape functions
C=BLOCK ABSTRACT
C
C     HEXA20SHAPE computes the value of the shape functions for a
C     twenty-noded isoparametric hexahedron and its
C     x-y-z derivatives, at a sample  point given by its hexahedron
C     coordinates (xi,eta,mu)
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       CALL HEXA20SHAPE (HF, XI,ETA,MU, X,Y,Z, SF,SX,SY,SZ, DET)
C
C     Input arguments:
C
C       HF        A dummy character argument
C       XI,ETA,MU  Hexahedral (iso-P) coordinates of given point
C       X         (20 x 1) array of x coordinates of hexahedron corners
C       Y         (20 x 1) array of y coordinates of hexahedron corners
C       Z         (20 x 1) array of z coordinates of hexahedron corners
C
C     Outputs arguments:
C
C       SF        (20 x 1) array of shape function values
C       SX        (20 x 1) array of shape function x-derivatives
C       SY        (20 x 1) array of shape function y-derivatives
C       SZ        (20 x 1) array of shape function z-derivatives 
C       DET        Value of Jacobian determinant
C
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine    HEXA20SHAPE
     $             (hf, xi, eta, mu, x, y, z, sf, sx, sy, sz, det)
C
C                   A R G U M E N T S
C
      character*(*)     hf
      double precision  xi, eta, mu, x(20), y(20), z(20)
      double precision  sf(20), sx(20), sy(20), sz(20), det
C
C                   L O C A L   V A R I A B L E S
C
      integer            i
      double precision   r, s, t, rp, rm, sp, sm, tp, tm, rrm, ssm, ttm
      double precision   xd1, xd2, xd3, yd1, yd2, yd3, zd1, zd2, zd3
      double precision   a11, a12, a13, a21, a22, a23, a31, a32, a33
      double precision   cdet
      double precision   sd1(20), sd2(20), sd3(20)
C
C                   L O G I C
C
C-------------------------------------------- FORM BASIC FUNCTION VALUES
      r  =  xi
      s  =  eta
      t  =  mu
      rp =  1.d0 + r
      rm =  1.d0 - r
      sp =  1.d0 + s
      sm =  1.d0 - s
      tp =  1.d0 + t
      tm =  1.d0 - t
      rrm = 1.d0 - r*r
      ssm = 1.d0 - s*s
      ttm = 1.d0 - t*t
C     +-----------------------------------------------------------------
C     I Q U A D R A T I C  SHAPE FUNCTIONS AND THEIR NATURAL DERIVATIVES
C     I20-NODED BRICK ELEMENT 
C     +-----------------------------------------------------------------
      sf(1)  = 0.125d0*rm*sm*tm*(rm+sm+tm-5.d0)
      sf(2)  = 0.125d0*rp*sm*tm*(rp+sm+tm-5.d0)
      sf(3)  = 0.125d0*rp*sp*tm*(rp+sp+tm-5.d0)
      sf(4)  = 0.125d0*rm*sp*tm*(rm+sp+tm-5.d0)
      sf(5)  = 0.125d0*rm*sm*tp*(rm+sm+tp-5.d0)
      sf(6)  = 0.125d0*rp*sm*tp*(rp+sm+tp-5.d0)
      sf(7)  = 0.125d0*rp*sp*tp*(rp+sp+tp-5.d0)
      sf(8)  = 0.125d0*rm*sp*tp*(rm+sp+tp-5.d0)
      sf(9)  = 0.25d0*rrm*sm*tm
      sf(10) = 0.25d0*rp*ssm*tm
      sf(11) = 0.25d0*rrm*sp*tm
      sf(12) = 0.25d0*rm*ssm*tm
      sf(13) = 0.25d0*rm*sm*ttm
      sf(14) = 0.25d0*rp*sm*ttm
      sf(15) = 0.25d0*rp*sp*ttm
      sf(16) = 0.25d0*rm*sp*ttm
      sf(17) = 0.25d0*rrm*sm*tp
      sf(18) = 0.25d0*rp*ssm*tp
      sf(19) = 0.25d0*rrm*sp*tp
      sf(20) = 0.25d0*rm*ssm*tp
C--------------------------------------- DERIVATIVE EVALUATION
      sd1(1)  = -0.125d0*sm*tm*(2.d0*rm+sm+tm-5.d0)
      sd1(2)  =  0.125d0*sm*tm*(2.d0*rp+sm+tm-5.d0)
      sd1(3)  =  0.125d0*sp*tm*(2.d0*rp+sp+tm-5.d0)
      sd1(4)  = -0.125d0*sp*tm*(2.d0*rm+sp+tm-5.d0)
      sd1(5)  = -0.125d0*sm*tp*(2.d0*rm+sm+tp-5.d0)
      sd1(6)  =  0.125d0*sm*tp*(2.d0*rp+sm+tp-5.d0)
      sd1(7)  =  0.125d0*sp*tp*(2.d0*rp+sp+tp-5.d0)
      sd1(8)  = -0.125d0*sp*tp*(2.d0*rm+sp+tp-5.d0)
      sd1(9)  = -0.5d0*r*sm*tm
      sd1(10) =  0.25d0*ssm*tm
      sd1(11) = -0.5d0*r*sp*tm
      sd1(12) = -sd1(10)
      sd1(14) =  0.25d0*sm*ttm
      sd1(13) = -sd1(14)
      sd1(15) =  0.25d0*sp*ttm
      sd1(16) = -sd1(15)
      sd1(17) = -0.5d0*r*sm*tp
      sd1(18) =  0.25d0*ssm*tp
      sd1(19) = -0.5d0*r*sp*tp
      sd1(20) = -sd1(18)
      sd2(1)  = -0.125d0*tm*rm*(2.d0*sm+tm+rm-5.d0)
      sd2(2)  = -0.125d0*tm*rp*(2.d0*sm+tm+rp-5.d0)
      sd2(3)  =  0.125d0*tm*rp*(2.d0*sp+tm+rp-5.d0)
      sd2(4)  =  0.125d0*tm*rm*(2.d0*sp+tm+rm-5.d0)
      sd2(5)  = -0.125d0*tp*rm*(2.d0*sm+tp+rm-5.d0)
      sd2(6)  = -0.125d0*tp*rp*(2.d0*sm+tp+rp-5.d0)
      sd2(7)  =  0.125d0*tp*rp*(2.d0*sp+tp+rp-5.d0)
      sd2(8)  =  0.125d0*tp*rm*(2.d0*sp+tp+rm-5.d0)
      sd2(10) = -0.5d0*s*tm*rp
      sd2(11) =  0.25d0*rrm*tm
      sd2(9)  = -sd2(11)
      sd2(12) = -0.5d0*s*tm*rm
      sd2(14) = -0.25d0*ttm*rp
      sd2(15) = -sd2(14)
      sd2(16) =  0.25d0*ttm*rm
      sd2(13) = -sd2(16)
      sd2(19) =  0.25d0*rrm*tp
      sd2(17) = -sd2(19)
      sd2(18) = -0.5d0*s*tp*rp
      sd2(20) = -0.5d0*s*tp*rm
      sd3(1)  = -0.125d0*rm*sm*(2.d0*tm+rm+sm-5.d0)
      sd3(2)  = -0.125d0*rp*sm*(2.d0*tm+rp+sm-5.d0)
      sd3(3)  = -0.125d0*rp*sp*(2.d0*tm+rp+sp-5.d0)
      sd3(4)  = -0.125d0*rm*sp*(2.d0*tm+rm+sp-5.d0)
      sd3(5)  =  0.125d0*rm*sm*(2.d0*tp+rm+sm-5.d0)
      sd3(6)  =  0.125d0*rp*sm*(2.d0*tp+rp+sm-5.d0)
      sd3(7)  =  0.125d0*rp*sp*(2.d0*tp+rp+sp-5.d0)
      sd3(8)  =  0.125d0*rm*sp*(2.d0*tp+rm+sp-5.d0)
      sd3(9)  = -0.25d0*rrm*sm
      sd3(10) = -0.25d0*ssm*rp
      sd3(11) = -0.25d0*rrm*sp
      sd3(12) = -0.25d0*ssm*rm
      sd3(13) = -0.5d0*t*rm*sm
      sd3(14) = -0.5d0*t*rp*sm
      sd3(15) = -0.5d0*t*rp*sp
      sd3(16) = -0.5d0*t*rm*sp
      sd3(18) = -sd3(10)
      sd3(19) = -sd3(11)
      sd3(20) = -sd3(12)
      sd3(17) = -sd3(9)
      xd1   =    0.0
      xd2   =    0.0
      xd3   =    0.0
      yd1   =    0.0
      yd2   =    0.0
      yd3   =    0.0
      zd1   =    0.0
      zd2   =    0.0
      zd3   =    0.0
      do 1500   i = 1,20
         xd1  =  xd1 +  x(i) * sd1(i)
         xd2  =  xd2 +  x(i) * sd2(i)
         xd3  =  xd3 +  x(i) * sd3(i)
         yd1  =  yd1 +  y(i) * sd1(i)
         yd2  =  yd2 +  y(i) * sd2(i)
         yd3  =  yd3 +  y(i) * sd3(i)
         zd1  =  zd1 +  z(i) * sd1(i)
         zd2  =  zd2 +  z(i) * sd2(i)
         zd3  =  zd3 +  z(i) * sd3(i)
 1500  continue
       a11 = yd2*zd3 - yd3*zd2
       a12 = yd3*zd1 - yd1*zd3
       a13 = yd1*zd2 - yd2*zd1
       a21 = xd3*zd2 - xd2*zd3
       a22 = xd1*zd3 - zd1*xd3
       a23 = xd2*zd1 - xd1*zd2
       a31 = xd2*yd3 - xd3*yd2
       a32 = yd1*xd3 - xd1*yd3
       a33 = xd1*yd2 - yd1*xd2
       det = xd1*a11 + yd1*a21 + zd1*a31
       if ( det .eq. 0.0)  then
            return
       end if
       cdet  =  1.0/det
       do 2000  i = 1,20
          sx(i) =  cdet * ( a11 * sd1(i) + a12 * sd2(i) + a13 * sd3(i))
          sy(i) =  cdet * ( a21 * sd1(i) + a22 * sd2(i) + a23 * sd3(i))
          sz(i) =  cdet * ( a31 * sd1(i) + a32 * sd2(i) + a33 * sd3(i))
 2000    continue
       return
       end
C=END FORTRAN
C=DECK HEXAGAUSSQ
C=PURPOSE Get s-point abscissas and weight for hexa product Gauss rule
C=AUTHOR C. A. Felippa, August 1966
C=VERSION May 1982 (Fortran 77)
C=EQUIPMENT Machine independent
C=KEYWORDS hexahedron Gauss integration rule abscissae weight
C=BLOCK ABSTRACT
C
C     HEXAGAUSSQ returns the hexahedron coordinates of
C     sample points and weights for a Gauss-product integration
C     rule over an isoparametric hexahedron.
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       CALL      HEXAGAUSSQ  (P1, I1, P2, I2, P3, I3,  XI, ETA, MU,
C                              WEIGHT)
C
C     Input arguments:
C
C       P1        Number of Gauss points in the XI direction
C       I1        Index of sample point in the XI direction
C       P2        Number of Gauss points in the ETA direction
C       I2        Index of sample point in the ETA direction
C       P3        Number of Gauss points in the MU direction
C       I3        Index of sample points in the MU direction
C        
C
C     Outputs arguments:
C
C       XI, ETA, MU   Hexahedral coordinates of sample point
C       WEIGHT    Weight factor
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine    HEXAGAUSSQ
     $             (p1, i1, p2, i2, p3, i3, xi, eta, mu, weight)
C
C                   A R G U M E N T S
C
      integer           p1, i1, p2, i2, p3, i3
      double precision  xi, eta, mu, weight
C
C                   L O C A L   V A R I A B L E S
C
      double precision  w1, w2, w3
C
C                   L O G I C
C
      call      LINEGAUSSQ (p1, i1, xi,  w1)
      call      LINEGAUSSQ (p2, i2, eta, w2)
      call      LINEGAUSSQ (p3, i3, mu,  w3)
      weight =  w1 * w2 * w3
      return
      end
C=END FORTRAN
