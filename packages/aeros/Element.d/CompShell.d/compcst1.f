C=====================================================================C
C        1         2         3         4         5         6         7C
C234567890123456789012345678901234567890123456789012345678901234567890C
C=====================================================================C
      subroutine compcst1( e       , thick   , nu     , coef   , nlayer,
     $                    idlayer , mtlayer , x      , y      , z     ,
     $                    d       , type    , eframe , aframe , effect)
C=====================================================================C
C                                                                     C
C     Perform =   Assembles the 3 by 3 Constitutive Matrix According  C
C     ---------   to the Type of Constitutive Law Requested.          C
C                                                                     C
C                                                                     C
C     Input/Output =                                                  C
C     --------------                                                  C
C     E       <input>  Young modulus                                  C
C     THICK   <input>  thickness (assumed constant over the element)  C
C     NU      <input>  Poisson's ratio                                C
C     COEF    <input>  coefficients of the constitutive law           C
C     NLAYER  <input>  number of layers of the composite element      C
C     IDLAYER <input>  identificators for each layer                  C
C     MTLAYER <input>  material properties of each layer              C
C     X       <input>  triangular coordinates in local x-direction    C
C     Y       <input>  triangular coordinates in local y-direction    C
C     Z       <input>  triangular coordinates in local z-direction    C
C     D       <output> 3 by 3 constitutive matrix                     C
C     TYPE    <input>  type of constitutive law                       C
C     EFRAME  <input>  element level 3x3 frame                        C
C     AFRAME  <input>  arbitrary 3x3 frame of the constitutive law    C
C     EFFECT  <input>  type of matrix [D] requested                   C
C                                                                     C
C=====================================================================C
C=Author  = Francois M. Hemez                                         C
C=Date    = June 9th, 1995                                            C
C=Version = 2.0                                                       C
C=Comment =                                                           C
C=====================================================================C
C
C     ------------
C     DECLARATIONS
C     ------------
C
C.....Global Variables
C
      integer     type , nlayer , idlayer(5,nlayer)
      real*8      e , nu , thick , coef(36) , d(3,3)
      real*8      x(3) , y(3) , z(3) , mtlayer(12,nlayer)
      real*8      eframe(3,3) , aframe(3,3)
      character   effect*2
C
C.....Local Variables
C
      integer     i , j
      real*8      zero , one , pi , twopi
      real*8      thetaD , theta
      real*8      Qbar(3,3) , T(3,3) , R(3)
      real*8      invT(3,3) , costheta , sintheta
      real*8      qt1 , qt2 , qt3
      real*8      refvec(3)
      real*8      norm1 , norm2 , normref , proj1 , proj2
      real*8      cosine1 , cosine2 , orifiber(3)
      logical     pureben , puremem , cbenmem , cmemben
c PJSA 5-11-2007 tensor rotation matrix and inverse
c      real*8      rtrmat(3,3), tmatinv(3,3)
c      real*8      c,s
C
C     ----
C     DATA
C     ----
C
      data zero   /0.000000D+00/
      data one    /1.000000D+00/
C
C.....INITIALIZE THE MAIN DIAGONAL OF REUTER'S MATRIX
C
      data R      /1.000000D+00, 1.000000D+00, 2.000000D+00/
C
C     -----
C     LOGIC
C     -----
C
      pi    = acos(-one)
      twopi = 2.000000D+00*pi
C
C.....CHECK IF THE ELEMENT IS A COMPOSITE SHELL
C
      if ( type.eq.-1 ) go to 100
C
C.....CHECK IF TYPE OF CONTITUTIVE LAW IS CORRECT
C
      if (type.ne.1) go to 200
C
C.....CHECK THE PHYSICAL EFFECT
C
      if (      (effect.ne."BB")
     $     .and.(effect.ne."MM")
     $     .and.(effect.ne."BM")
     $     .and.(effect.ne."MB") ) go to 300
C
C.....SET THE TYPE OF PHYSICAL EFFECT
C
      pureben = ( effect.eq."BB" )
      puremem = ( effect.eq."MM" )
      cbenmem = ( effect.eq."BM" )
      cmemben = ( effect.eq."MB" )
C
C.....CLEAR THE 3 BY 3 CONSTITUTIVE MATRIX
C
      do 1001 j=1,3
         do 1002 i=1,3
            d(i,j) = zero
 1002    continue
 1001 continue
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE BENDING
C
      if ( pureben ) then
C
         d(1,1) = coef(22)
         d(1,2) = coef(23)
         d(1,3) = coef(24)
         d(2,1) = coef(28)
         d(2,2) = coef(29)
         d(2,3) = coef(30)
         d(3,1) = coef(34)
         d(3,2) = coef(35)
         d(3,3) = coef(36)
C
      endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE MEMBRANE
C
      if ( puremem ) then
C
         d(1,1) = coef( 1)
         d(1,2) = coef( 2)
         d(1,3) = coef( 3)
         d(2,1) = coef( 7)
         d(2,2) = coef( 8)
         d(2,3) = coef( 9)
         d(3,1) = coef(13)
         d(3,2) = coef(14)
         d(3,3) = coef(15)
C
      endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING BENDING-MEMBRANE
C
      if ( cbenmem ) then
C
         d(1,1) = coef(19)
         d(1,2) = coef(20)
         d(1,3) = coef(21)
         d(2,1) = coef(25)
         d(2,2) = coef(26)
         d(2,3) = coef(27)
         d(3,1) = coef(31)
         d(3,2) = coef(32)
         d(3,3) = coef(33)
C
      endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING MEMBRANE-BENDING
C
      if ( cmemben ) then
C
         d(1,1) = coef( 4)
         d(1,2) = coef( 5)
         d(1,3) = coef( 6)
         d(2,1) = coef(10)
         d(2,2) = coef(11)
         d(2,3) = coef(12)
         d(3,1) = coef(16)
         d(3,2) = coef(17)
         d(3,3) = coef(18)
C
      endif
C
C.....SET THE REFERENCE VECTOR FOR ORIENTING THE FIBERS OF THE LAYER
C.....(ALWAYS TAKE THE FIRST VECTOR OF THE FRAME)
C
      refvec(1) = aframe(1,1)
      refvec(2) = aframe(2,1)
      refvec(3) = aframe(3,1)
C
C.....PROJECT THE REFERENCE VECTOR INTO THE PLANE OF
C.....THE ELEMENT TO GET THE FIBER ORIENTATION VECTOR
C
      orifiber(1) = zero
      orifiber(2) = zero
      orifiber(3) = zero
C
      norm1   = zero
      norm2   = zero
      normref = zero
      proj1   = zero
      proj2   = zero
C
      do 2004 i=1,3
         norm1   =   norm1 + eframe(i,1)*eframe(i,1)
         norm2   =   norm2 + eframe(i,2)*eframe(i,2)
         normref = normref +   refvec(i)*refvec(i)
         proj1   =   proj1 + eframe(i,1)*refvec(i)
         proj2   =   proj2 + eframe(i,2)*refvec(i)
 2004 continue
C
      norm1   = sqrt(norm1)
      norm2   = sqrt(norm2)
      normref = sqrt(normref)
C
      if ( normref.eq.zero ) then
         cosine1 = one
         cosine2 = zero
      else
         cosine1 = proj1/(norm1*normref)
         cosine2 = proj2/(norm2*normref)
      endif
C
      do 2005 i=1,3
         orifiber(i) = cosine1*eframe(i,1) + cosine2*eframe(i,2)
 2005 continue
C
C.....CALCULATE THE ANGLE FROM THE [x] FRAME OF THE LOCAL
C.....COORDINATE SYSTEM TO THE REFERENCE ORIENTATION VECTOR
C
      proj1  = zero
      proj2  = zero
      thetaD = zero
C
      do 2006 i=1,3
         proj1 = proj1 + eframe(i,1)*orifiber(i)
         proj2 = proj2 + eframe(i,2)*orifiber(i)
 2006 continue
C
      if ( proj1.eq.zero ) then
         if ( proj2.eq.zero ) go to 700
         if ( proj2.gt.zero ) thetaD = 0.50D+00*pi
         if ( proj2.lt.zero ) thetaD = 1.50D+00*pi
      else
         if ( proj2.eq.zero ) then
            if ( proj1.eq.zero ) go to 700
            if ( proj1.gt.zero ) thetaD = zero
            if ( proj1.lt.zero ) thetaD = pi
         else
            thetaD = atan(proj2/proj1)
         endif
      endif
C
      if ( thetaD.lt.zero ) then
         thetaD = thetaD + twopi
      endif
C
      if ( (thetaD.lt.zero).or.(thetaD.gt.twopi) ) go to 800
C
C.....CALCULATE THE ANGLE FROM THE [x] FRAME OF THE LOCAL
C.....COORDINATE SYSTEM TO THE DIRECTION OF THE FIBER
C
      theta = thetaD 
C
      if ( theta.gt.twopi ) then
         theta = theta - twopi
      endif
C
C.....INITIALIZE THE ROTATION MATRIX FROM THE FIBER COORDINATE
C.....SYSTEM {1;2} TO THE ELEMENT TRIANGULAR SYSTEM {x;y}
C
      do 2206 j=1,3
         do 2207 i=1,3
            T(i,j) = zero
 2207    continue
 2206 continue
C
      costheta = cos(theta)
      sintheta = sin(theta)
C
      T(1,1)   =  costheta*costheta
      T(1,2)   =  sintheta*sintheta
      T(1,3)   =  2.00D+00*costheta*sintheta
C
      T(2,1)   =  sintheta*sintheta
      T(2,2)   =  costheta*costheta
      T(2,3)   = -2.00D+00*costheta*sintheta
C
      T(3,1)   = -costheta*sintheta
      T(3,2)   =  costheta*sintheta
      T(3,3)   =  costheta*costheta - sintheta*sintheta
C
C.....COMPUTE THE INVERSE OF [T]:
C.....[invT] = inverse(diag[R]) * [T]^t * diag[R]
C
      do 2208 j=1,3
         do 2209 i=1,3
            invT(i,j) = zero
 2209    continue
 2208 continue
C
      do 2210 j=1,3
         do 2211 i=1,3
            invT(i,j) = T(j,i)*(R(j)/R(i))
 2211    continue
 2210 continue
C
C.....CALCULATE THE INVERSE OF THE COMPLIANCE MATRIX WHICH RELATES
C.....THE STRAINS [ex], [ey] AND [exy] TO THE STRESSES [sx], [sy]
C.....AND [sxy] IN THE TRIANGULAR COORDINATE SYSTEM {x;y}:
C.....[Qbar] = [invT] * [Q] * [invT]^t
C
      do 2212 j=1,3
         do 2213 i=1,3
            Qbar(i,j) = zero
 2213    continue
 2212 continue
C
      do 2214 j=1,3
C
         qt1 = d(1,1)*invT(j,1) + d(1,2)*invT(j,2) + d(1,3)*invT(j,3)
         qt2 = d(2,1)*invT(j,1) + d(2,2)*invT(j,2) + d(2,3)*invT(j,3)
         qt3 = d(3,1)*invT(j,1) + d(3,2)*invT(j,2) + d(3,3)*invT(j,3)
C
         do 2215 i=1,3
            Qbar(i,j) = qt1*invT(i,1) + qt2*invT(i,2) + qt3*invT(i,3)
 2215    continue
C
 2214 continue
C
      do 3001 j=1,3
         do 3002 i=1,3
            d(i,j) = Qbar(i,j)
 3002    continue
 3001 continue
C
C     ------
C     RETURN
C     ------
C
      return
C
C     ------
C     FORMAT
C     ------
C
   91 format("*** Type Used is: ",I10,12x," ***")
   92 format("*** Type Used is: ",5x,A2,9x," ***")
C
C     ---------------
C     ERROR TREATMENT
C     ---------------
C
C.....ERROR-MESSAGE IF THE ELEMENT IS NOT A COMPOSITE SHELL
C
  100 continue
      write(*,*) "*** FATAL ERROR in Routine COMPCST1       ***"
      write(*,*) "*** The Finite Element is not a Composite ***"
      write(*,*) "*** Shell Element: Inconsistency Detected ***"
      write(*,*) "*** STOP ALL TREATMENTS RIGHT HERE        ***"
      stop
C
C.....ERROR-MESSAGE IF THE CONSTITUTIVE LAW IS INCORRECT
C
  200 continue
      write(*,*) "*** FATAL ERROR in Routine COMPCST1      ***"
      write(*,*) "*** Wrong Type of Constitutive Law       ***"
      write(*,91) type
      write(*,*) "*** Types Allowed are:                   ***"
      write(*,*) "*** 1 = given constitutive law           ***"
      write(*,*) "*** STOP ALL TREATMENTS RIGHT HERE       ***"
      stop
C
C.....ERROR-MESSAGE IF THE PHYSICAL EFFECT IS INCORRECT
C
  300 continue
      write(*,*) "*** FATAL ERROR in Routine COMPCST1 ***"
      write(*,*) "*** Wrong Type of Physical Effect   ***"
      write(*,92) effect
      write(*,*) "*** Types Allowed are:              ***"
      write(*,*) "*** BB = pure bending               ***"
      write(*,*) "*** MM = pure membrane              ***"
      write(*,*) "*** BM = coupling bending-membrane  ***"
      write(*,*) "*** MB = coupling membrane-bending  ***"
      write(*,*) "*** STOP ALL TREATMENTS RIGHT HERE  ***"
      stop
C
C.....ERROR-MESSAGE IF THE REFERENCE ORIENTATION IS BUGGY
C
  700 continue
      write(*,*) "*** FATAL ERROR in Routine COMPCST1  ***"
      write(*,*) "*** The Reference Orientation Vector ***"
      write(*,*) "*** is Parallel to the Two Inplane   ***"
      write(*,*) "*** and Orthogonal Local Frames!     ***"
      write(*,*) "*** STOP ALL TREATMENTS RIGHT HERE   ***"
      stop
C
C.....ERROR-MESSAGE IF THE ORIENTATION ANGLE IS OUT-OF-BOUNDS
C
  800 continue
      write(*,*) "*** FATAL ERROR in Routine COMPCST1   ***"
      write(*,*) "*** The Angle From the Local [x]      ***"
      write(*,*) "*** Axis of the Triangular Coordinate ***"
      write(*,*) "*** System to the Reference Direction ***"
      write(*,*) "*** is Out-of-Bounds: it Must be      ***"
      write(*,*) "*** Within the Range 0-2pi Radians    ***"
      write(*,*) "*** STOP ALL TREATMENTS RIGHT HERE    ***"
      stop
C
C     ---
C     END
C     ---
C
      end
C=end of routine "COMPCST1
C========================C
