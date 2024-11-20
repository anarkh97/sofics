C=MODULE  FORBACKC4
C=PURPOSE Forward/backward substitution of a symmetric positive
C=PURPOSE semi-definite system decomposed and stored in skyline form
C=PURPOSE and using variable band information -
C=PURPOSE FORWARD and BACKWARD are serial and unrolled to level 4  -
C=AUTHOR  M. LESOINNE A. PUPPIN-MACEDO
C=AUTHOR  March 1998 
C=BLOCK 
      subroutine FORBACKC4(COLVAL,LD,W1,W2,W3,W4, PIVOT, NEQ)
C----------------------------------------------------------------------------
C     FORBACKC4 = forward/backward substitution of a symmetric skyline 
C                matrix  using variable band information
C                with unrolling depth = 4 
C
C     Version 1.1 (complex)
C     =====================
C
C     COLVAL   (complex) Stores sequentially the values in the column profile
C              of K
C     LD       (integer) Locates the diagonal elements of K in COLVAL
C     W1---W4  (complex) 4 RHS vectors
C     PIVOT    (integer) stores the singularity information
C     NEQ      (integer) Number of equations in the system
C
C----------------------------------------------------------------------------
C
C-----DECLARATIONS
C
      INTEGER LD(*),PIVOT(*),NEQ,NOPS
      COMPLEX*16  COLVAL(*),W1(*),W2(*),W3(*),W4(*)
      COMPLEX*16  XTEMP
C
C-----Local DECLARATIONS
C
      INTEGER K,J,I,IW,LEN,I1,I2
C
      NOPS = 0
C
C-----FORWARD SUBSTITUTE
C
      IF(PIVOT(1).NE.0) THEN
         XTEMP = COLVAL(LD(1))
         W1(1) = W1(1)/XTEMP
         W2(1) = W2(1)/XTEMP
         W3(1) = W3(1)/XTEMP
         W4(1) = W4(1)/XTEMP
      ELSE
         W1(1) = 0.0
         W2(1) = 0.0
         W3(1) = 0.0
         W4(1) = 0.0
      ENDIF
C
      NOPS = NOPS + 1
C
      DO 500 J = 2, NEQ
      IF(PIVOT(J).NE.0) THEN
         LEN  = LD(J) - LD(J - 1) - 1
C  Break the dependency check between W1, W2, W3 and W4
CDIR$ IVDEP
         DO K = 1, LEN
            XTEMP = COLVAL(LD(J - 1) + K)
            W1(J) = W1(J)- XTEMP*W1(J - LEN+K-1)
            W2(J) = W2(J)- XTEMP*W2(J - LEN+K-1)
            W3(J) = W3(J)- XTEMP*W3(J - LEN+K-1)
            W4(J) = W4(J)- XTEMP*W4(J - LEN+K-1)
C           W1(J) =  W1(J) - DDOTC(LEN, XTEMP, W1(J - LEN+K-1))
C           W2(J) =  W2(J) - DDOTC(LEN, XTEMP, W2(J - LEN+K-1))
C           W3(J) =  W3(J) - DDOTC(LEN, XTEMP, W3(J - LEN+K-1))
C           W4(J) =  W4(J) - DDOTC(LEN, XTEMP, W4(J - LEN+K-1))
         END DO
         XTEMP = COLVAL(LD(J))
         W1(J) = W1(J)/XTEMP
         W2(J) = W2(J)/XTEMP
         W3(J) = W3(J)/XTEMP
         W4(J) = W4(J)/XTEMP

C
      NOPS = NOPS + 8*LEN
      ELSE
         W1(J) = 0.0
         W2(J) = 0.0
         W3(J) = 0.0
         W4(J) = 0.0
      ENDIF
C
500   CONTINUE
C
      DO 600 K = 1, NEQ
      XTEMP     = COLVAL(LD(K))
      W1(K)     = W1(K)*XTEMP
      W2(K)     = W2(K)*XTEMP
      W3(K)     = W3(K)*XTEMP
      W4(K)     = W4(K)*XTEMP
600   CONTINUE
C
      NOPS = NOPS + NEQ
C
C
C
C-----BACKWARD SUBSTITUTE
C
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
         W1(K) = 0.0
         W2(K) = 0.0
         W3(K) = 0.0
         W4(K) = 0.0
      ELSE
         XTEMP = COLVAL(LD(K))
         W1(K) = W1(K)/XTEMP
         W2(K) = W2(K)/XTEMP
         W3(K) = W3(K)/XTEMP
         W4(K) = W4(K)/XTEMP
C
         NOPS = NOPS + 4
C
      END IF
C
C     Sweep
C
      I1 = LD(K - 1) + 1
      I2 = LD(K) - 1
      IF(I1.GT.I2) GO TO 700
      IW = K - (LD(K) - LD(K - 1)) + 1
CVD$ NODEPCHK
C  Break the dependency check between COLVAL(I) and COLVAL(I+OFFSET)
CDIR$ IVDEP
         DO 710 I = I1, I2
         XTEMP     = COLVAL(I)
         W1(IW)    = W1(IW) - XTEMP*W1(K)
         W2(IW)    = W2(IW) - XTEMP*W2(K)
         W3(IW)    = W3(IW) - XTEMP*W3(K)
         W4(IW)    = W4(IW) - XTEMP*W4(K)
         IW       = IW + 1
710      CONTINUE
C
      NOPS = NOPS + 8*(LD(K) - LD(K - 1) - 1)
C
700   CONTINUE
C
C     Handle Singularities
C
      IF(PIVOT(1).EQ.0) THEN
         W1(1) = 0.0
         W2(1) = 0.0
         W3(1) = 0.0
         W4(1) = 0.0
      ELSE
         XTEMP = COLVAL(LD(1))
         W1(1) = W1(1)/XTEMP
         W2(1) = W2(1)/XTEMP
         W3(1) = W3(1)/XTEMP
         W4(1) = W4(1)/XTEMP
C
         NOPS = NOPS + 1
C
      END IF
C
C
      RETURN
      END















