C=MODFIED: 1-19-98
C=MODULE  PFACT
C=PURPOSE Parallel/Vector solution of a symmetric positive
C=PURPOSE semi-definite system stored in skyline form
C=PURPOSE and using variable band information -
C=PURPOSE LDLt algorithm -
C=PURPOSE FACTOR is fully parallelized
C=PURPOSE and loop unrolls to level 4 -
C=PURPOSE FORWARD and BACKWARD are serial and rolled -
C=PURPOSE If not empty, Null Space (ZERO ENERGY MODES) is computed 
C=AUTHOR C. FARHAT and M. LESOINNE
C=AUTHOR March 1998
C=BLOCK Force
      subroutine pfact(COLVAL,LD,LACOL,B,W,PIVOT,
     &                  TOL,NEQ,FLAG,NOPS,NZEM,ZEM,PROCID,
     &                  NPROC, BAR, ll)
C----------------------------------------------------------------------------
C     SVBU4R =   skyline solver  using variable band information
C                with unrolling depth = 4 and returning eventual
C                ZERO ENERGY MODES.
C
C     Version 3:
C     ==========
C     FACTOR is fully parallelized and loop unrolls to level 4
C     FORWARD and BACKWARD are serial and rolled
C     Null Space (ZERO ENERGY MODES) is also computed.
C
C     Force routine for solution of K u = r, where K is symmetric
C     stored in a profile form. Zero pivots are flagged as PIVOT(I) = 0.
C     Null space (ZERO ENERGY MODES) is computed via a modified
C     Farhat's Lecture Notes method.
C
C     COLVAL   Stores sequentially the values in the column profile of K
C     LD       Locates the diagonal elements of K in COLVAL
C     B        Right hand side
C     W        Working vectors arranged in a matrix
C     PIVOT    stores the singularity information
C     TOL      user computed/defined tolerance for singularity
C     NEQ      Number of equations in the system
C     FLAG     = 1 Factor only (with Preprocessing of
C                               of eventual Null Space)
C              = 2 Retrieve Null Space (ZERO ENERGY MODES) only
C              = 3 Forward  solve only
C              = 4 Backward solve only
C              = 5 Forward and Backward solve only
C
C     NZEM     Number of ZERO ENERGY MODES 
C              --> with FLAG = 1, NZEM is computed and returned
C
C     ZEM      Stores the ZERO ENERGY MODES 
C----------------------------------------------------------------------------
C
C-----DECLARATIONS

      INTEGER LD(*),LACOL(*),PIVOT(*),NEQ,FLAG
      INTEGER NZEM,NOPS
      REAL*8  COLVAL(*),W(4,*),B(*),TOL
      REAL*8  ZEM(NEQ,*)
      INTEGER BAR,ll
      INTEGER NPROC
C
      INTEGER PROCID
C
C-----Local DECLARATIONS
C
      INTEGER K,J,ISTART,ISTOP,I,IW,KR
      REAL*8  M11,M21,M22,M31,M32,M33,F1,F2,F3,F4
      REAL*8  G1,G2,G3,G4
      INTEGER P1,P2,P3,P4, OFFSET
      REAL*8  INVPIV(4)
      INTEGER JREM
C      INTEGER DDD
      EXTERNAL ussetlock, usunsetlock
      INTEGER ussetlock, usunsetlock
C
      NOPS = 0
C
C-----INITIALIZE      
C
      IF(PROCID.EQ.0) THEN
        NZEM = 0
C
        DO IW  = 1, NEQ
          PIVOT(IW) = 1
          ZEM(IW,1) = ABS(COLVAL(LD(IW)))
        END DO
      END IF 
C
C-----FACTOR COLVAL
C
      DO 100 K = 1, NEQ - 3, 4
C
C     TRIANGULAR ZONE
C
C
C     Step k
C
      IF(PROCID.EQ.0) THEN
        P1 = LD(K)
        P2 = LD(K + 1)
        P3 = LD(K + 2)
        P4 = LD(K + 3)
C
C     Handle Singularities
C
        IF(ABS(COLVAL(P1)).LT.TOL*ZEM(K,1)) THEN
           PIVOT(K)      = 0
c          write(6,*) 'k = ',k,'  ',colval(P1)/ZEM(K,1)
           NZEM          = NZEM + 1
           COLVAL(P1)    = 0.0
           INVPIV(1)     = 0.0
        ELSE
           INVPIV(1)     = 1.0/COLVAL(P1)
        END IF
C
C
        M11 = COLVAL(P2 - 1)*INVPIV(1)
C
        COLVAL(P2)     = COLVAL(P2)     - M11*COLVAL(P2 - 1)
        COLVAL(P3 - 1) = COLVAL(P3 - 1) - M11*COLVAL(P3 - 2)
        COLVAL(P4 - 2) = COLVAL(P4 - 2) - M11*COLVAL(P4 - 3)
C
C     Step k + 1
C
C
C     Handle Singularities
C
        IF(ABS(COLVAL(P2)).LT.TOL*ZEM(K+1,1)) THEN
           PIVOT(K + 1) = 0
c          write(6,*)'k = ',k+1,'  ',colval(p2)/ZEM(K+1,1)
           NZEM         = NZEM + 1
           COLVAL(P2)   = 0.0
           INVPIV(2)    = 0.0
        ELSE
           INVPIV(2)    = 1.0/COLVAL(P2)
        END IF
C
        M21 = COLVAL(P3 - 2)*INVPIV(1)
        M22 = COLVAL(P3 - 1)*INVPIV(2)
C
        COLVAL(P3) = COLVAL(P3)- M21*COLVAL(P3 - 2) - M22*COLVAL(P3 - 1)
        COLVAL(P4 - 1) = COLVAL(P4 - 1) - M21*COLVAL(P4 - 3)
     .                                  - M22*COLVAL(P4 - 2)
C
        NOPS = NOPS + 17
C
C
C     Step k + 2
C
C
C     Handle Singularities
C
        IF(ABS(COLVAL(P3)).LT.TOL*ZEM(K+2,1)) THEN
           PIVOT(K + 2) = 0
c          write(6,*)'k = ',k+2,'  ',colval(p3)/ZEM(K+2,1)
           NZEM         = NZEM + 1
           COLVAL(P3)   = 0.0
           INVPIV(3)    = 0.0     
        ELSE
           INVPIV(3)    = 1.0/COLVAL(P3)
        END IF
C
        M31 = COLVAL(P4 - 3)*INVPIV(1)
        M32 = COLVAL(P4 - 2)*INVPIV(2)
        M33 = COLVAL(P4 - 1)*INVPIV(3)
C
        COLVAL(P4) = COLVAL(P4) - M31*COLVAL(P4 - 3) 
     .                          - M32*COLVAL(P4 - 2)
     .                          - M33*COLVAL(P4 - 1)
C
C     Handle Singularities
C
        IF(ABS(COLVAL(P4)).LT.TOL*ZEM(K+3,1)) THEN
c          write(6,*)'k = ',k+3,'  ',colval(p4)/ZEM(K+3,1)
           PIVOT(K + 3)  = 0
           NZEM          = NZEM + 1
           COLVAL(P4)    = 0.0
           INVPIV(4)     = 0.0
        ELSE
           INVPIV(4)     = 1.0/COLVAL(P4)
        END IF
C
        NOPS = NOPS + 9
C
C
C     ITS HORIZONTAL SHADOW
C
CVD$  NODEPCHK
           DO J = K + 4, LACOL(K)
             IF((LD(J) - LD(J - 1) - J + K).GT.0) THEN
                COLVAL(LD(J) - J + K + 1) = COLVAL(LD(J) - J + K + 1)
     .                                    - M11*COLVAL(LD(J) - J + K)
                COLVAL(LD(J) - J + K + 2) = COLVAL(LD(J) - J + K + 2)
     .                                  - M21*COLVAL(LD(J) - J + K)
     .                                  - M22*COLVAL(LD(J) - J + K + 1)
                COLVAL(LD(J) - J + K + 3) = COLVAL(LD(J) - J + K + 3)
     .                                  - M31*COLVAL(LD(J) - J + K)
     .                                  - M32*COLVAL(LD(J) - J + K + 1)
     .                                  - M33*COLVAL(LD(J) - J + K + 2)
                W(1,J - K - 3) = COLVAL(LD(J) - J + K)    *INVPIV(1)
                W(2,J - K - 3) = COLVAL(LD(J) - J + K + 1)*INVPIV(2)
                W(3,J - K - 3) = COLVAL(LD(J) - J + K + 2)*INVPIV(3)
                W(4,J - K - 3) = COLVAL(LD(J) - J + K + 3)*INVPIV(4)
C
                NOPS = NOPS + 19
C
             ELSE
                W(1,J - K - 3) = 0.0D0
                W(2,J - K - 3) = 0.0D0
                W(3,J - K - 3) = 0.0D0
                W(4,J - K - 3) = 0.0D0
             END IF
           END DO
         ENDIF
C
C
C     REMAINING REGION
C  Create a barrier
C        WRITE(6,*) 'Barrier',NPROC,PROCID
C        CALL FLUSH(6)
         CALL BARRIER(BAR, NPROC)
C
C         DO 400 J = K + 4, LACOL(K)
C Compute my first column
         JREM = MOD(K+3, 2*NPROC)
         J = K+4-JREM+2*PROCID
C Here we adjust back to do columns bigger than K+4. Notice that
C all the numbers will be done including K+4 thanks to the fact that K is
C always even
         IF(J.LT.K+4) J = J+2*NPROC
C        DDD = ussetlock(ll)
C        WRITE(0,*) J,K+4,procid
C        DDD = usunsetlock(ll)
         DO WHILE (J.LE.LACOL(K))
           IF((LD(J) - LD(J - 1) - J + K).GT.0) THEN
              ISTART   = LD(J) - J + K + 4
              ISTOP    = LD(J)
              F1 = COLVAL(ISTART - 4)
              F2 = COLVAL(ISTART - 3)
              F3 = COLVAL(ISTART - 2)
              F4 = COLVAL(ISTART - 1)
              IF(J.LT.LACOL(K) .AND.
     &           (LD(J+1) - LD(J) - J -1 + K).GT.0) THEN
                OFFSET = LD(J+1) - LD(J) - 1
                G1 = COLVAL(OFFSET + ISTART - 4)
                G2 = COLVAL(OFFSET + ISTART - 3)
                G3 = COLVAL(OFFSET + ISTART - 2)
                G4 = COLVAL(OFFSET + ISTART - 1)
                IW=1
C  Break the dependency check between COLVAL(I) and COLVAL(I+OFFSET)
CDIR$ IVDEP
                DO I  = ISTART, ISTOP
                   COLVAL(I) = COLVAL(I) - F1*W(1,IW) - F2*W(2,IW)
     &                                 - F3*W(3,IW) - F4*W(4,IW)
                   COLVAL(I+OFFSET) = COLVAL(I+OFFSET)
     &                                 - G1*W(1,IW) - G2*W(2,IW)
     &                                 - G3*W(3,IW) - G4*W(4,IW)
                   IW        = IW + 1
                ENDDO
                COLVAL(OFFSET+ISTOP+1) = COLVAL(OFFSET+ISTOP+1)
     &                                 - G1*W(1,IW) - G2*W(2,IW)
     &                                 - G3*W(3,IW) - G4*W(4,IW)
                J=J+2*NPROC
              ELSE
                IW = 1
CVD$ NODEPCHK
                   DO I  = ISTART, ISTOP
                   COLVAL(I) = COLVAL(I) - F1*W(1,IW) - F2*W(2,IW)
     .                                 - F3*W(3,IW) - F4*W(4,IW)
                   IW        = IW + 1
                   END DO
C
                NOPS = NOPS + 8*(ISTOP - ISTART + 1)
                J=J+2*NPROC
             ENDIF
           ELSE
             J = J+1
             IF(J.LE.LACOL(K).AND.
     .           (LD(J) - LD(J - 1) - J + K).GT.0) THEN
               ISTART   = LD(J) - J + K + 4
               ISTOP    = LD(J)
               F1 = COLVAL(ISTART - 4)
               F2 = COLVAL(ISTART - 3)
               F3 = COLVAL(ISTART - 2)
               F4 = COLVAL(ISTART - 1)
               IW = 1
CVD$ NODEPCHK
               DO I  = ISTART, ISTOP
                   COLVAL(I) = COLVAL(I) - F1*W(1,IW) - F2*W(2,IW)
     .                                 - F3*W(3,IW) - F4*W(4,IW)
                   IW        = IW + 1
               END DO
C
               NOPS = NOPS + 8*(ISTOP - ISTART + 1)
             ENDIF
             J=J-1+2*NPROC
           END IF
C
      END DO
      CALL BARRIER(BAR, NPROC)
C
100   CONTINUE
C
C     REMAINDER STEPS
C
C
      IF(PROCID.NE.0) RETURN
      KR = MOD(NEQ,4)
      IF(KR.EQ.3) THEN
C
C        Handle Singularities
C
         IF(ABS(COLVAL(LD(NEQ - 2))).LT.TOL*ZEM(NEQ-2,1)) THEN
            PIVOT(NEQ - 2)      = 0
c        write(6,*)'k = ',neq-2,'  ',colval(ld(neq-2))/ZEM(NEQ-2,1)
            NZEM                = NZEM + 1
            COLVAL(LD(NEQ - 2)) = 0.0
            INVPIV(1)           = 0.0
         ELSE
            INVPIV(1)           = 1.0/COLVAL(LD(NEQ - 2))
         END IF
C
         M11 = COLVAL(LD(NEQ - 1) - 1)*INVPIV(1)
         COLVAL(LD(NEQ - 1)) = COLVAL(LD(NEQ - 1))
     .                       - M11*COLVAL(LD(NEQ - 1) - 1)
         COLVAL(LD(NEQ) - 1) = COLVAL(LD(NEQ) - 1)
     .                       - M11*COLVAL(LD(NEQ) - 2)
C
         NOPS = NOPS + 5
C
C        Handle Singularities
C
         IF(ABS(COLVAL(LD(NEQ - 1))).LT.TOL*ZEM(NEQ-1,1)) THEN
            PIVOT(NEQ - 1)      = 0
c           write(6,*)'k = ',neq-1,'  ',colval(ld(neq-1))/ZEM(NEQ-1,1)
            NZEM                = NZEM + 1
            COLVAL(LD(NEQ - 1)) = 0.0
            INVPIV(2)           = 0.0
         ELSE
            INVPIV(2)           = 1.0/COLVAL(LD(NEQ - 1))
         END IF
C
         M21 = COLVAL(LD(NEQ) - 2)*INVPIV(1)
         M22 = COLVAL(LD(NEQ) - 1)*INVPIV(2)
         COLVAL(LD(NEQ)) = COLVAL(LD(NEQ))
     .                   - M21*COLVAL(LD(NEQ) - 2)
     .                   - M22*COLVAL(LD(NEQ) - 1)
C
         NOPS = NOPS + 6
C
C        Handle Singularities
C
         IF(ABS(COLVAL(LD(NEQ))).LT.TOL*ZEM(NEQ,1)) THEN
            PIVOT(NEQ)      = 0
c           write(6,*)'k = ',neq,'  ',colval(ld(neq))/ZEM(NEQ,1)
            NZEM            = NZEM + 1
            COLVAL(LD(NEQ)) = 0.0
         END IF
C
      ELSE IF (KR.EQ.2) THEN
C
C        Handle Singularities
C
         IF(ABS(COLVAL(LD(NEQ - 1))).LT.TOL*ZEM(NEQ-1,1)) THEN
            PIVOT(NEQ - 1)      = 0
c           write(6,*)'k = ',neq-1,'  ',colval(ld(neq-1))/ZEM(NEQ-1,1)
            NZEM                = NZEM + 1
            COLVAL(LD(NEQ - 1)) = 0.0
            INVPIV(1)           = 0.0
         ELSE
            INVPIV(1)           = 1.0/COLVAL(LD(NEQ - 1))
         END IF
C
         M11 = COLVAL(LD(NEQ) - 1)*INVPIV(1)
         COLVAL(LD(NEQ)) = COLVAL(LD(NEQ)) - M11*COLVAL(LD(NEQ) - 1)
C
         NOPS = NOPS + 3
C
C        Handle Singularities
C
         IF(ABS(COLVAL(LD(NEQ))).LT.TOL*ZEM(NEQ,1)) THEN
c           write(6,*)'k = ',neq,'  ',colval(ld(neq))/ZEM(NEQ,1)
            PIVOT(NEQ)      = 0
            NZEM            = NZEM + 1
            COLVAL(LD(NEQ)) = 0.0
         END IF
C
      ELSE IF (KR.EQ.1) THEN
C
C        Handle Singularities
C
         IF(ABS(COLVAL(LD(NEQ))).LT.TOL*ZEM(NEQ,1)) THEN
c           write(6,*)'k = ',neq,'  ',colval(ld(neq))/ZEM(NEQ,1)
            PIVOT(NEQ)      = 0
            NZEM            = NZEM + 1
            COLVAL(LD(NEQ)) = 0.0
         END IF
C
      END IF
C
      RETURN
      END

