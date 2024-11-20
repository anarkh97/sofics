C=MODULE  SVBU4GMB
C=PURPOSE Parallel/Vector solution of a symmetric positive
C=PURPOSE semi-definite system stored in skyline form
C=PURPOSE and using variable band information -
C=PURPOSE LDLt algorithm -
C=PURPOSE FACTOR is fully parallelized
C=PURPOSE and loop unrolls to level 4 -
C=PURPOSE FORWARD and BACKWARD are serial and rolled -
C=PURPOSE If not empty, Null Space (ZERO ENERGY MODES) is computed 
C=PURPOSE by the HYBRID GEOMETRIC-ALGEBRAIC METHOD
C=AUTHOR C. FARHAT 
C=AUTHOR June 1996
C=BLOCK Force
      subroutine SVBU4GMB(COLVAL,LD,LACOL,B,W,W3,W4,PIVOT,TOL,NEQ
     .,                         FLAG,NOPS,NZEM,ZEM,KSTOP,S1,S2,NGRBM,
     .                          SEQID) 
C----------------------------------------------------------------------------
C*************************************************************
C
C WARNING: This is a modified version of SVBU4G
C          Modification is related to the case FLAG = 11
C
C*************************************************************
C
C     SVBU4R =   skyline solver  using variable band information
C                with unrolling depth = 4 and returning eventual
C                ZERO ENERGY MODES.
C     
C
C     MEMORY SCHEME and KSTOP BUSINESS
C     --------------------------------
C              kstop = 4*((neq-10)/4)
C              if(kstop.lt.4) kstop = 4
C              call alocat(reaend,rlassz,s1)
C              call alocat(s1,(neq - kstop)*(neq - kstop),s2)
C              call alocat(s2,(neq - kstop)*(neq - kstop),reaend)
C              rlassz = (neq - kstop)*(neq - kstop)
C
C
C ****************
C  W A R N I N G
C ****************
C     This solver must be used only when the zero energy modes are to
C     be captured by the hybrid geometric-algebraic method. In particular
C     it requires that the equation numbering has been altered to ensure
C     that the ``small'' pivots associated with the NRBM rigid body modes
C     found previously by the geometric method are located AFTER the
C     equation corresponding to ``kstop''.
C     
C
C     Version 3:
C     ==========
C     FACTOR is fully parallelized and loop unrolls to level 4
C     FORWARD and BACKWARD are serial and rolled
C     Null Space (ZERO ENERGY MODES) is also computed.
C
C     Force routine for solution of K u = r, where K is symmetric
C     stored in a profile form. Zero pivots are flagged as PIVOT(I) = 0.
C
C     COLVAL   Stores sequentially the values in the column profile
C              of K
C     LD       Locates the diagonal elements of K in COLVAL
C     B        Right hand side
C     W        4 Working vectors stored in a matrix
C     S1---S2  2 Other Working Vectors
C     PIVOT    stores the singularity information
C     TOL      user computed/defined tolerance for singularity
C     NEQ      Number of equations in the system
C     KSTOP    Step at which the factorization process in interrupted to
C              shift to a full-pivoting version in another subroutine
C     FLAG     = 1  Factor only (if matrix is singular, factorization
C                                is incomplete)
C     FLAG     = 11 Complete factorization for a singular matrix
C                   and retrive the part of the null space 
C                   associated with mechanisms
C              = 2  Retrieve Complete Null Space (ZERO ENERGY MODES) 
C              = 3  Forward  solve only
C              = 4  Backward solve only
C              = 5  Forward and Backward solve only
C
C     NZEM     Number of ZERO ENERGY MODES
C              --> with FLAG = 1, NZEM returns the additional zero energy 
C                  modes that are not geometric rigid body modes
C              --> with FLAG = 2, NZEM must be input as the total number
C                  of zero energy modes
C
C     NGRBM    Number of GEOMETRIC RIGID BODY MODES
C
C     ZEM      Stores the ZERO ENERGY MODES 
C----------------------------------------------------------------------------
C
C-----DECLARATIONS

      INTEGER LD(*),LACOL(*),PIVOT(*),SEQID(*),NEQ,FLAG
      INTEGER NZEM,KSTOP,NGRBM
      REAL*8  COLVAL(*),W(4,*),W3(*),W4(*),B(*),TOL,NOPS
      REAL*8  ZEM(NEQ,*),S1(*),S2(*)
      REAL*8  CHDDOT
C
C-----Local DECLARATIONS
C
      INTEGER K,J,ISTART,ISTOP,I,IW,LEN,KR,I1,I2
      INTEGER IR,IZEM,ICOL,MNZEM
      REAL*8  M11,M21,M22,M31,M32,M33,F1,F2,F3,F4
      REAL*8  G1,G2,G3,G4
      INTEGER P1,P2,P3,P4,OFFSET
      REAL*8  INVPIV(4)
C
C-----OVERDIMENSION SEQID (Cheap and Effective)
C
C     INTEGER SEQID(100)
C
      NOPS = 0
C
      IF(FLAG.EQ.11) GOTO 11
      IF(FLAG.EQ.2)  GO TO 2
      IF(FLAG.EQ.3)  GO TO 3
      IF(FLAG.EQ.4)  GO TO 4
      IF(FLAG.EQ.5)  GO TO 3
C
C-----INITIALIZE      
C
      NZEM  = 0
C
      DO 10 IW  = 1, NEQ
      PIVOT(IW) = 1
10    CONTINUE
C
C-----FACTOR COLVAL
C
      DO 100 K = 1, KSTOP - 3, 4
C
C     TRIANGULAR ZONE
C
C
C     Step k
C
      P1 = LD(K)
      P2 = LD(K + 1)
      P3 = LD(K + 2)
      P4 = LD(K + 3)
C
C     Handle Singularities
C
      IF(ABS(COLVAL(P1)).LT.TOL) THEN
         PIVOT(K)      = 0
C        write(6,*)' mechanism      : k = ',k,' pivot = ',colval(p1)
         NZEM          = NZEM + 1
         SEQID(NZEM)   = K
         COLVAL(P1)    = 0.0
         INVPIV(1)     = 0.0
      ELSE
         INVPIV(1)     = 1.0/COLVAL(P1)
      END IF
C
      M11            = COLVAL(P2 - 1)*INVPIV(1)
      COLVAL(P2)     = COLVAL(P2)     - M11*COLVAL(P2 - 1)
      COLVAL(P3 - 1) = COLVAL(P3 - 1) - M11*COLVAL(P3 - 2)
      COLVAL(P4 - 2) = COLVAL(P4 - 2) - M11*COLVAL(P4 - 3)
C
C
C     Step k + 1
C
C
C     Handle Singularities
C
      IF(ABS(COLVAL(P2)).LT.TOL) THEN
         PIVOT(K + 1) = 0
C        write(6,*)' mechanism      : k = ',k,' pivot = ',colval(p2)
         NZEM         = NZEM + 1
         SEQID(NZEM)   = K + 1
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
     .                                - M22*COLVAL(P4 - 2)
C
      NOPS = NOPS + 17
C
C
C     Step k + 2
C
C
C     Handle Singularities
C
      IF(ABS(COLVAL(P3)).LT.TOL) THEN
         PIVOT(K + 2) = 0
C        write(6,*)' mechanism      : k = ',k+2,' pivot = ',colval(p3)
         NZEM         = NZEM + 1
         SEQID(NZEM)  = K + 2
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
      COLVAL(P4) = COLVAL(P4) - M31*COLVAL(P4 - 3) - M32*COLVAL(P4 - 2)
     .                        - M33*COLVAL(P4 - 1)
C
C     Handle Singularities
C
      IF(ABS(COLVAL(P4)).LT.TOL) THEN
C        write(6,*)' mechanism      : k = ',k+3,' pivot = ',colval(p4)
         PIVOT(K + 3)  = 0
         NZEM          = NZEM + 1
         SEQID(NZEM)  = K + 3
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
         DO 200 J = K + 4, LACOL(K)
         IF((LD(J) - LD(J - 1) - J + K).GT.0) THEN
            COLVAL(LD(J) - J + K + 1) = COLVAL(LD(J) - J + K + 1)
     .                                  - M11*COLVAL(LD(J) - J + K)
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
200      CONTINUE
C     REMAINING REGION
C
C         DO 400 J = K + 4, LACOL(K)
         J = K+4
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
     &                                   - F3*W(3,IW) - F4*W(4,IW)
                   COLVAL(I+OFFSET) = COLVAL(I+OFFSET)
     &                                   - G1*W(1,IW) - G2*W(2,IW)
     &                                   - G3*W(3,IW) - G4*W(4,IW)
                   IW        = IW + 1
                ENDDO
                COLVAL(OFFSET+ISTOP+1) = COLVAL(OFFSET+ISTOP+1)
     &                                   - G1*W(1,IW) - G2*W(2,IW)
     &                                   - G3*W(3,IW) - G4*W(4,IW)
                J=J+2
              ELSE
                IW = 1
CVD$ NODEPCHK
                   DO 410 I  = ISTART, ISTOP
                   COLVAL(I) = COLVAL(I) - F1*W(1,IW) - F2*W(2,IW)
     .                                   - F3*W(3,IW) - F4*W(4,IW)
                   IW        = IW + 1
410               CONTINUE
C
                NOPS = NOPS + 8*(ISTOP - ISTART + 1)
                J=J+1
             ENDIF
           ELSE
             J = J+1
           END IF
C
        ENDDO
CKHP
100   CONTINUE
C
C     REMAINDER STEPS
C
C
      KR = MOD(KSTOP,4)
      IF(KR.NE.0) THEN
         WRITE(6,*)'Internal Error in Subroutine SVBU4G'
         WRITE(6,*)'kstop must be a multiple of 4'
         WRITE(6,*)'STOP in SVBU4G'
         STOP
      ENDIF
C
      IF(FLAG.EQ.1) RETURN
C
C-----FORWARD SUBSTITUTE
C
3     CONTINUE
C
C
      IF(PIVOT(1).NE.0) THEN
         B(1) = B(1)/COLVAL(LD(1))
      ELSE
         B(1) = 0.0
      ENDIF
C
      NOPS = NOPS + 1
C
      DO 500 J = 2, NEQ
      IF(PIVOT(J).NE.0) THEN
         LEN  = LD(J) - LD(J - 1) - 1
         B(J) = (B(J) - CHDDOT(LEN,COLVAL(LD(J - 1) + 1),B(J - LEN)))
     .          / COLVAL(LD(J))
C
      NOPS = NOPS + 2*LEN
      ELSE
         B(J) = 0.0
      ENDIF
C
500   CONTINUE
C
      DO 600 K = 1, NEQ
      B(K)     = B(K)*COLVAL(LD(K))
600   CONTINUE
C
      NOPS = NOPS + NEQ
C
C
      IF(FLAG.EQ.3) RETURN
C
C-----BACKWARD SUBSTITUTE
C
4     CONTINUE
C
C
      DO 700 K = NEQ, 2, -1
C
C     Compute Solution Xk
C
C
C     Handle Singularities
C
      IF(PIVOT(K).EQ.0) THEN
         B(K) = 0.0
      ELSE
         B(K) = B(K)/COLVAL(LD(K))
C
         NOPS = NOPS + 1
C
      END IF
C
      IF(B(K).EQ.0.0) GO TO 700
C
C     Sweep
C
      I1 = LD(K - 1) + 1
      I2 = LD(K) - 1
      IF(I1.GT.I2) GO TO 700
      IW = K - (LD(K) - LD(K - 1)) + 1
CVD$ NODEPCHK
         DO 710 I = I1, I2
         B(IW)    = B(IW) - COLVAL(I)*B(K)
         IW       = IW + 1
710      CONTINUE
C
      NOPS = NOPS + 2*(LD(K) - LD(K - 1) - 1)
C
700   CONTINUE
C
C     Handle Singularities
C
      IF(PIVOT(1).EQ.0) THEN
         B(1) = 0.0
      ELSE
         B(1) = B(1)/COLVAL(LD(1))
C
         NOPS = NOPS + 1
C
      END IF
C
C
      IF(FLAG.EQ.4) RETURN
      IF(FLAG.EQ.5) RETURN
C
11    CONTINUE
2     CONTINUE
C
C-----RECOVER EVENTUAL ZERO ENERGY MODES 
C
C     Switch to a factorization algorithm with full pivoting
C     following the hybrid geometric-algebraic method.
C     Remember that this method requires that the equation numbering
C     be such that the rigid body mode pivots are in the
C     equations AFTER kstop. Note also that FFULPV preprocesses
C     the recovery of the zero energy modes
C
C      Note: I am passing a real array W4 whereas FFULPV needs an integer array.
C            Since I dont need to access the information stored in W4 by FFULPV
C            this is ok
C
      ICOL = NZEM - NGRBM + 1
      CALL  FFULPVM(NEQ,LD,COLVAL,(NEQ-KSTOP),NGRBM,PIVOT
     &                     ,SEQID(ICOL),ZEM(1,ICOL),S1,S2
     &                                        ,W3,W4,FLAG)
C
      IF(NZEM.EQ.0)  RETURN
C
C-----SCATTER AN IDENTITY SUBMATRIX 
C-----IN THE COLUMNS OF ZEM THAT HAVE NOT
C-----BEEN PREPROCESSED BY FFULPV
C     (Note: SVBU4G and  FFULPV do not use the
C      same technique for recovering singular vectors)
C
      IF(NZEM.GT.NGRBM) THEN
         IZEM = 1
         DO 760 I = 1, NEQ
            DO 770 IR = 1, NZEM - NGRBM
            ZEM(I,IR) = 0.0
770         CONTINUE
         IF(PIVOT(I).EQ.0) THEN
            ZEM(I,IZEM) = 1.0
            IZEM        = IZEM + 1
         ENDIF
760      CONTINUE
      ENDIF
C
C-----COMPUTE THE ZERO ENERGY MODES
C-----VIA BACKWARD SUBSTITUTION
C     (Note: FLAG 11 means that we will extract
C      only the zem due to mechanisms as the others
C      are provided by the geometric procedure)
C
      IF(FLAG.EQ.11) THEN
         MNZEM = NZEM - NGRBM
      ELSE
         MNZEM = NZEM
      ENDIF
C
      DO 800 IR = 1, MNZEM

         DO 810 K = NEQ, 2, -1
C
C        Compute Solution 
C
         IF(PIVOT(K).NE.0) THEN
            ZEM(K,IR) = ZEM(K,IR)/COLVAL(LD(K))
         ENDIF
C
         NOPS = NOPS + 1
C
         IF(ZEM(K,IR).EQ.0.0) GO TO 810
C
C        Sweep
C
         I1 = LD(K - 1) + 1
         I2 = LD(K) - 1
         IF(I1.GT.I2) GO TO 810
         IW = K - (LD(K) - LD(K - 1)) + 1
CVD$ NODEPCHK
           DO 820   I = I1, I2
           ZEM(IW,IR) = ZEM(IW,IR) - COLVAL(I)*ZEM(K,IR)
           IW       = IW + 1
820        CONTINUE
C
         NOPS = NOPS + 2*(LD(K) - LD(K - 1) - 1)
C
810      CONTINUE
C
         IF(PIVOT(1).NE.0) THEN
            ZEM(1,IR) = ZEM(1,IR)/COLVAL(LD(1))
         ENDIF
C
         NOPS = NOPS + 1
C
800   CONTINUE
C
C
      RETURN
      END
C
C
C
C=end force
