C   
C  Algorithm SPNNLS: SPARSE NONNEGATIVE LEAST SQUARES
C   
C  The original version of this code was developed by
C  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
C  1973 JUN 15, and published in the book
C  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
C  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C
C  Adapted by Julien Cortial at Stanford University
C
C     GIVEN AN M BY N MATRIX, A, AND AN M-VECTOR, B,  COMPUTE A
C     SPARSE N-VECTOR, X, THAT VERIFIES
C   
C         ||A * X - B||_2 <= RELTOL * ||B||_2  SUBJECT TO X .GE. 0   
C     ------------------------------------------------------------------
C                     Subroutine Arguments
C
C     A(),MDA,M,N     MDA IS THE FIRST DIMENSIONING PARAMETER FOR THE   
C                     ARRAY, A().   ON ENTRY A() CONTAINS THE M BY N    
C                     MATRIX, A.           ON EXIT A() CONTAINS 
C                     THE PRODUCT MATRIX, Q*A , WHERE Q IS AN   
C                     M BY M ORTHOGONAL MATRIX GENERATED IMPLICITLY BY  
C                     THIS SUBROUTINE.  
C     B()     ON ENTRY B() CONTAINS THE M-VECTOR, B.   ON EXIT B() CON- 
C             TAINS Q*B.
C     X()     ON ENTRY X() NEED NOT BE INITIALIZED.  ON EXIT X() WILL   
C             CONTAIN THE SOLUTION VECTOR. 
C     RELTOL  RELATIVE TOLERANCE
C             (STOPPING CRITERION: ||B - A*X||_2 < RELTOL * ||B||_2).
C     RNORM   ON EXIT RNORM CONTAINS THE EUCLIDEAN NORM OF THE  
C             RESIDUAL VECTOR.  
C     W()     AN N-ARRAY OF WORKING SPACE.  ON EXIT W() WILL CONTAIN    
C             THE DUAL SOLUTION VECTOR.   W WILL SATISFY W(I) = 0.  
C             FOR ALL I IN SET P  AND W(I) .LE. 0. FOR ALL I IN SET Z   
C     ZZ()     AN M-ARRAY OF WORKING SPACE.     
C     ZZ2()    AN N-ARRAY OF WORKING SPACE.     
C     INDEX()     AN INTEGER WORKING ARRAY OF LENGTH AT LEAST N.
C                 ON EXIT THE CONTENTS OF THIS ARRAY DEFINE THE SETS    
C                 P AND Z AS FOLLOWS..  
C   
C                 INDEX(1)   THRU INDEX(NSETP) = SET P.     
C                 INDEX(IZ1) THRU INDEX(IZ2)   = SET Z.     
C                 IZ1 = NSETP + 1 = NPP1
C                 IZ2 = N   
C     MODE    THIS IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING 
C             MEANINGS. 
C             1     THE SOLUTION HAS BEEN COMPUTED SUCCESSFULLY.
C             2     THE DIMENSIONS OF THE PROBLEM ARE BAD.  
C                   EITHER M .LE. 0 OR N .LE. 0.
C             3     ITERATION COUNT EXCEEDED.
C                   MORE THAN MAXITE*N ITERATIONS. 
C     MAXSZE  ALTERNATIVE STOPPING CRITERION BASED ON ACTIVE SET SIZE
C             NSETP <= MIN(M,MAXSZE*N)
C   
C     ------------------------------------------------------------------
      SUBROUTINE SPNNLS (A,MDA,M,N,B,X,RELTOL,RNORM,W,ZZ,ZZ2,INDEX,MODE,
     +                   PRTFLG,SCAFLG,MAXSZE,MAXITE,DTIME)
C     ------------------------------------------------------------------
      USE ISO_C_BINDING
      IMPLICIT NONE
C     include 'mpif.h'
      integer(kind=C_LONG) I, II, IP, ITER, ITMAX, IZ, IZ1, IZ2, IZMAX,
     +                     J, JJ, JZ, L, DDATE
      integer(kind=C_LONG) M, MDA, MODE, N, NPP1, NSETP, RTNKEY, SPMAX
      integer(kind=C_LONG) INDEX(*)  
      double precision A(MDA,*), B(*), W(*), X(*), ZZ(*), ZZ2(*) 
      double precision ALPHA, ASAVE, CC, DIFF, DUMMY, FACTOR, RNORM
      double precision ABSTOL,RELTOL,MAXSZE,MAXITE
      REAL(8)::T1,T2,DTIME
      double precision ONE, SM, SS, T, TEMP, TWO, UNORM, UP, WMAX
      double precision ZERO, ZTEST
      integer(kind=C_LONG) PRTFLG, SCAFLG
      parameter(FACTOR = 0.01d0)
      parameter(TWO = 2.0d0, ONE = 1.0d0, ZERO = 0.0d0)
C     ------------------------------------------------------------------
      MODE=1
      IF (M .le. 0 .or. N .le. 0) then
         MODE=2
         RETURN
      endif
      ITER=0
      ITMAX=INT(MAXITE*N) 
      SPMAX=MIN(M,INT(MAXSZE*N,C_LONG))
      DDATE=0
      DTIME=0.0
      T1=0.0
      T2=0.0
C   
C                    INITIALIZE THE ARRAYS INDEX() AND X(). 
C   
      DO 20 I=1,N   
         X(I)=ZERO     
         INDEX(I)=I    
   20 CONTINUE 
C   
      IZ2=N 
      IZ1=1 
      NSETP=0   
      NPP1=1
C                    INIT ZZ2 = tr(A^T * A)^{-1}
      DO 25 I = 1,N
         IF(SCAFLG.NE.0) THEN
            ZZ2(I) = ZERO
            DO 26 L = 1,M
               ZZ2(I) = ZZ2(I) + A(L,I)**2
   26       CONTINUE
            ZZ2(I) = ONE / sqrt(ZZ2(I))
         ELSE
            ZZ2(I) = ONE
         ENDIF
   25 CONTINUE
C
C                    INIT ABSTOL = RELTOL^2 * ||B||^2
      SM=ZERO   
      DO 28 I=1,M   
   28     SM=SM+B(I)**2
      ABSTOL = RELTOL * sqrt(SM)
C                             ******  MAIN LOOP BEGINS HERE  ******     
   30 CONTINUE  
C
C                     COMPUTE THE NORM^2 OF THE RESIDUAL VECTOR.    
      SM=ZERO   
      DO 35 I=NPP1,M   
  35     SM=SM+B(I)**2
      RNORM = sqrt(SM)
C         
      IF(PRTFLG.NE.0) THEN
        write (*,500) 'Iteration = ', ITER,
     +                'Downdate = ', DDATE,
     +                'Active set size = ', NSETP,
     +                'Residual norm = ', RNORM,
     +                'Target = ', ABSTOL 
      ENDIF
  500 FORMAT (A, I9, 4X, A, I9, 4X, A, I9, 4X, A, ES13.6, 4X, A, ES13.6)
C
C                  QUIT IF ALL COEFFICIENTS ARE ALREADY IN THE SOLUTION.
C                        OR IF M COLS OF A HAVE BEEN TRIANGULARIZED.    
C
      IF (IZ1.GT.IZ2.OR.NSETP.GE.SPMAX) GO TO 350
C 
C         STOPPING CRITERION
      IF (RNORM .LT. ABSTOL) GO TO 350
C         
C         COMPUTE COMPONENTS OF THE DUAL (NEGATIVE GRADIENT) VECTOR W().
C   
      DO 50 IZ=IZ1,IZ2  
         J=INDEX(IZ)   
         SM=ZERO   
         DO 40 L=NPP1,M
   40        SM=SM+A(L,J)*B(L)     
         W(J)=SM*ZZ2(J)   
   50 continue
C                                   FIND LARGEST POSITIVE W(J). 
   60 continue
      WMAX=ZERO 
      DO 70 IZ=IZ1,IZ2  
         J=INDEX(IZ)   
         IF (W(J) .gt. WMAX) then
            WMAX=W(J)     
            IZMAX=IZ  
         endif
   70 CONTINUE  
C   
C             IF WMAX .LE. 0. GO TO TERMINATION.
C             THIS INDICATES SATISFACTION OF THE KUHN-TUCKER CONDITIONS.
C   
      IF (WMAX .le. ZERO) go to 350
      IZ=IZMAX  
      J=INDEX(IZ)   
C   
C     THE SIGN OF W(J) IS OK FOR J TO BE MOVED TO SET P.    
C     BEGIN THE TRANSFORMATION AND CHECK NEW DIAGONAL ELEMENT TO AVOID  
C     NEAR LINEAR DEPENDENCE.   
C   
      ASAVE=A(NPP1,J)   
      CALL H12 (1,NPP1,NPP1+1,M,A(1,J),1,UP,DUMMY,1,1,0)    
      UNORM=ZERO
      IF (NSETP .ne. 0) then
          DO 90 L=1,NSETP   
   90       UNORM=UNORM+A(L,J)**2     
      endif
      UNORM=sqrt(UNORM) 
      IF (DIFF(UNORM+ABS(A(NPP1,J))*FACTOR,UNORM) .gt. ZERO) then
C   
C        COL J IS SUFFICIENTLY INDEPENDENT.  COPY B INTO ZZ, UPDATE ZZ
C        AND SOLVE FOR ZTEST ( = PROPOSED NEW VALUE FOR X(J) ).    
C   
         DO 120 L=1,M  
  120        ZZ(L)=B(L)    
         CALL H12 (2,NPP1,NPP1+1,M,A(1,J),1,UP,ZZ,1,1,1)   
         ZTEST=ZZ(NPP1)/A(NPP1,J)  
C   
C                                     SEE IF ZTEST IS POSITIVE  
C   
         IF (ZTEST .gt. ZERO) go to 140
      endif
C   
C     REJECT J AS A CANDIDATE TO BE MOVED FROM SET Z TO SET P.  
C     RESTORE A(NPP1,J), SET W(J)=0., AND LOOP BACK TO TEST DUAL
C     COEFFS AGAIN.     
C   
      A(NPP1,J)=ASAVE   
      W(J)=ZERO 
      GO TO 60  
C   
C     THE INDEX  J=INDEX(IZ)  HAS BEEN SELECTED TO BE MOVED FROM
C     SET Z TO SET P.    UPDATE B,  UPDATE INDICES,  APPLY HOUSEHOLDER  
C     TRANSFORMATIONS TO COLS IN NEW SET Z,  ZERO SUBDIAGONAL ELTS IN   
C     COL J,  SET W(J)=0.   
C   
  140 continue
      DO 150 L=1,M  
  150    B(L)=ZZ(L)    
C   
      INDEX(IZ)=INDEX(IZ1)  
      INDEX(IZ1)=J  
      IZ1=IZ1+1 
      NSETP=NPP1
      NPP1=NPP1+1   
C   
      IF (IZ1 .le. IZ2) then
         DO 160 JZ=IZ1,IZ2 
            JJ=INDEX(JZ)  
            CALL H12 (2,NSETP,NPP1,M,A(1,J),1,UP,A(1,JJ),1,MDA,1)
  160    continue
      endif
C   
      IF (NSETP .ne. M) then
         DO 180 L=NPP1,M   
  180       A(L,J)=ZERO   
      endif
C   
      W(J)=ZERO 
C                                SOLVE THE TRIANGULAR SYSTEM.   
C                                STORE THE SOLUTION TEMPORARILY IN ZZ().
      RTNKEY = 1
      GO TO 400 
  200 CONTINUE  
C   
C                       ******  SECONDARY LOOP BEGINS HERE ******   
C   
C                          ITERATION COUNTER.   
C 
  210 continue  
      ITER=ITER+1   
      IF (ITER .gt. ITMAX) then
         MODE=3
         GO TO 350 
      endif
C   
C                    SEE IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE.    
C                                  IF NOT COMPUTE ALPHA.    
C   
      ALPHA=TWO 
      DO 240 IP=1,NSETP 
         L=INDEX(IP)   
         IF (ZZ(IP) .le. ZERO) then
            T=-X(L)/(ZZ(IP)-X(L))     
            IF (ALPHA .gt. T) then
               ALPHA=T   
               JJ=IP 
            endif
         endif
  240 CONTINUE  
C   
C          IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE THEN ALPHA WILL   
C          STILL = 2.    IF SO EXIT FROM SECONDARY LOOP TO MAIN LOOP.   
C   
      IF (ALPHA.EQ.TWO) GO TO 330   
C     T1 = MPI_Wtime()
      DDATE = DDATE+1
C   
C          OTHERWISE USE ALPHA WHICH WILL BE BETWEEN 0. AND 1. TO   
C          INTERPOLATE BETWEEN THE OLD X AND THE NEW ZZ.    
C   
      DO 250 IP=1,NSETP 
         L=INDEX(IP)   
         X(L)=X(L)+ALPHA*(ZZ(IP)-X(L)) 
  250 continue
C   
C        MODIFY A AND B AND THE INDEX ARRAYS TO MOVE COEFFICIENT I  
C        FROM SET P TO SET Z.   
C   
      I=INDEX(JJ)   
C     write(*,*) 'removing index ',(I-1),'ZZ(JJ)=',ZZ(JJ)
  260 continue
      X(I)=ZERO 
C   
      IF (JJ .ne. NSETP) then
         JJ=JJ+1   
         DO 280 J=JJ,NSETP 
            II=INDEX(J)   
            INDEX(J-1)=II 
            CALL G1 (A(J-1,II),A(J,II),CC,SS,A(J-1,II))   
            A(J,II)=ZERO  
            DO 270 L=1,N  
               IF (L.NE.II) then
C
C                 Apply procedure G2 (CC,SS,A(J-1,L),A(J,L))  
C
                  TEMP = A(J-1,L)
                  A(J-1,L) = CC*TEMP + SS*A(J,L)
                  A(J,L)   =-SS*TEMP + CC*A(J,L)
               endif
  270       CONTINUE  
C
C                 Apply procedure G2 (CC,SS,B(J-1),B(J))   
C
            TEMP = B(J-1)
            B(J-1) = CC*TEMP + SS*B(J)    
            B(J)   =-SS*TEMP + CC*B(J)    
  280    continue
      endif
C
      NPP1=NSETP
      NSETP=NSETP-1     
      IZ1=IZ1-1 
      INDEX(IZ1)=I  
C   
C        SEE IF THE REMAINING COEFFS IN SET P ARE FEASIBLE.  THEY SHOULD
C        BE BECAUSE OF THE WAY ALPHA WAS DETERMINED.
C        IF ANY ARE INFEASIBLE IT IS DUE TO ROUND-OFF ERROR.  ANY   
C        THAT ARE NONPOSITIVE WILL BE SET TO ZERO   
C        AND MOVED FROM SET P TO SET Z. 
C   
      DO 300 JJ=1,NSETP 
         I=INDEX(JJ)   
         IF (X(I) .le. ZERO) go to 260
  300 CONTINUE  
C   
C         COPY B( ) INTO ZZ( ).  THEN SOLVE AGAIN AND LOOP BACK.
C   
      DO 310 I=1,M  
  310     ZZ(I)=B(I)    
      RTNKEY = 2
      GO TO 400 
  320 CONTINUE  
C     T2 = MPI_Wtime()
C     write(*,*) 'cpu time 1 = ',T1,' cpu time 2 = ',T2
C     DTIME = DTIME + T2 - T1
      GO TO 210 
C                      ******  END OF SECONDARY LOOP  ******
C   
  330 continue
      DO 340 IP=1,NSETP 
          I=INDEX(IP)   
  340     X(I)=ZZ(IP)   
C        ALL NEW COEFFS ARE POSITIVE.  LOOP BACK TO BEGINNING.  
      GO TO 30  
C   
C                        ******  END OF MAIN LOOP  ******   
C   
C                        COME TO HERE FOR TERMINATION.  
C                     COMPUTE THE NORM OF THE FINAL RESIDUAL VECTOR.    
C 
  350 continue  
      SM=ZERO   
      IF (NPP1 .le. M) then
         DO 360 I=NPP1,M   
  360       SM=SM+B(I)**2 
      else
         DO 380 J=1,N  
  380       W(J)=ZERO     
      endif
      RNORM=sqrt(SM)    
      RETURN
C   
C     THE FOLLOWING BLOCK OF CODE IS USED AS AN INTERNAL SUBROUTINE     
C     TO SOLVE THE TRIANGULAR SYSTEM, PUTTING THE SOLUTION IN ZZ().     
C   
  400 continue
      DO 430 L=1,NSETP  
         IP=NSETP+1-L  
         IF (L .ne. 1) then
            DO 410 II=1,IP
               ZZ(II)=ZZ(II)-A(II,JJ)*ZZ(IP+1)   
  410       continue
         endif
         JJ=INDEX(IP)  
         ZZ(IP)=ZZ(IP)/A(IP,JJ)    
  430 continue
      go to (200, 320), RTNKEY
      END   
