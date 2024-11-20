C=MODULE  FORBACKR2
C=PURPOSE Forward/backward substitution of a symmetric positive
C=PURPOSE semi-definite system decomposed and stored in skyline form
C=PURPOSE and using variable band information -
C=PURPOSE FORWARD and BACKWARD are serial 
C=AUTHOR  M. LESOINNE, A. PUPPIN-MACEDO, K. H. PIERSON
C=AUTHOR  May 1998
C=BLOCK
      subroutine FORBACKR2(COLVAL,LD,W1,W2,PIVOT,NEQ)
C----------------------------------------------------------------------------
C     FORBACKR2 = forward/backward substitution of a symmetric skyline
C                 matrix using variable band information
C                 with unrolling depth = 1
C
C     Version 1.1 (real)
C     =====================
C
C     COLVAL   (real) Stores sequentially the values in the column profile of K
C     LD       (integer) Locates the diagonal elements of K in COLVAL
C     W1---W2  (real) 2 RHS vectors
C     PIVOT    (integer) stores the singularity information
C     NEQ      (integer) Number of equations in the system
C
C----------------------------------------------------------------------------
C
C-----DECLARATIONS
C
      INTEGER LD(*),PIVOT(*),NEQ
      DOUBLE PRECISION  COLVAL(*),W1(*),W2(*),XTEMP
C
C-----Local DECLARATIONS
C
      INTEGER K,J,I,IW,LEN,I1,I2
C
C-----FORWARD SUBSTITUTE
C
      IF(PIVOT(1).EQ.0) THEN
         W1(1) = 0.0
         W2(1) = 0.0
      ELSE
         XTEMP = COLVAL(LD(1))
         W1(1) = W1(1)/XTEMP
         W2(1) = W2(1)/XTEMP
      ENDIF
C
      DO 500 J = 2, NEQ
      IF(PIVOT(J).EQ.0) THEN
         W1(J) = 0.0
         W2(J) = 0.0
      ELSE
           I1   = LD(J - 1)
           I2   = LD(J)
           LEN  = I2 - I1 - 1
C  Break the dependency check between W1 and W2
CDIR$ IVDEP
         DO K = 1, LEN
            XTEMP = COLVAL(I1 + K)
            W1(J) = W1(J) - XTEMP*W1(J - LEN + K - 1)
            W2(J) = W2(J) - XTEMP*W2(J - LEN + K - 1)
         END DO
         XTEMP = COLVAL(I2)
         W1(J) = W1(J)/XTEMP
         W2(J) = W2(J)/XTEMP
      ENDIF
C
500   CONTINUE
C
C-----DIAGONAL SCALING
C
      DO 600 K = 1, NEQ
        XTEMP   = COLVAL(LD(K))
        W1(K)   = W1(K)*XTEMP
        W2(K)   = W2(K)*XTEMP
600   CONTINUE
C
C-----BACKWARD SUBSTITUTE
C
      DO 700 K = NEQ, 2, -1
C
C     Handle Singularities
C
      IF(PIVOT(K).EQ.0) THEN
         W1(K) = 0.0
         W2(K) = 0.0
         GO TO 700
      ELSE
         XTEMP = COLVAL(LD(K))
         W1(K) = W1(K)/XTEMP
         W2(K) = W2(K)/XTEMP
      END IF
C
C     Sweep
C
        I1 = LD(K - 1) + 1

        I2 = LD(K) - 1

        IF(I1.GT.I2) GO TO 700

        IW = K - LD(K) + ( LD(K - 1) + 1 )

C  Break the dependency check between W1 and W2
CDIR$ IVDEP
         DO 710 I = I1, I2
         XTEMP     = COLVAL(I)
         W1(IW)    = W1(IW) - XTEMP*W1(K)
         W2(IW)    = W2(IW) - XTEMP*W2(K)
         IW        = IW + 1
710      CONTINUE
700   CONTINUE
C
C     Handle Singularities
C
      IF(PIVOT(1).EQ.0) THEN
         W1(1) = 0.0
         W2(1) = 0.0
      ELSE
         XTEMP = COLVAL(LD(1))
         W1(1) = W1(1)/XTEMP
         W2(1) = W2(1)/XTEMP
      END IF
      RETURN
      END
C
C=END FORTRAN
C
