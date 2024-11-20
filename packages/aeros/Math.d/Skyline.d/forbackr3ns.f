C=MODULE  FORBACKR3NS
C=PURPOSE Forward/backward substitution of a symmetric positive
C=PURPOSE semi-definite system decomposed and stored in skyline form
C=PURPOSE and using variable band information -
C=PURPOSE FORWARD and BACKWARD are serial 
C=AUTHOR  M. LESOINNE, A. PUPPIN-MACEDO, K. H. PIERSON
C=AUTHOR  June 1998
C=BLOCK
      subroutine FORBACKR3NS(COLVAL,LD,W1,W2,W3,NEQ)
C----------------------------------------------------------------------------
C     FORBACKR3NS = forward/backward substitution of a symmetric skyline
C                 matrix using variable band information
C                 with unrolling depth = 1
C
C     Version 1.1 (real)
C     =====================
C
C     COLVAL   (real) Stores sequentially the values in the column profile of K
C     LD       (integer) Locates the diagonal elements of K in COLVAL
C     W1---W3  (real) 3 RHS vectors
C     NEQ      (integer) Number of equations in the system
C
C----------------------------------------------------------------------------
C
C-----DECLARATIONS
C
      INTEGER LD(*),NEQ
      REAL*8  COLVAL(*),W1(*),W2(*),W3(*),XTEMP
C
C-----Local DECLARATIONS
C
      INTEGER K,J,I,IW,LEN,I1,I2
C
C-----FORWARD SUBSTITUTE
C
         XTEMP = COLVAL(LD(1))
         W1(1) = W1(1)/XTEMP
         W2(1) = W2(1)/XTEMP
         W3(1) = W3(1)/XTEMP
C
      DO 500 J = 2, NEQ
         LEN  = LD(J) - LD(J - 1) - 1
C  Break the dependency check between W1,W2,W3
CDIR$ IVDEP
         DO K = 1, LEN
            XTEMP = COLVAL(LD(J - 1) + K)
            W1(J) = W1(J) - XTEMP*W1(J - LEN+K-1)
            W2(J) = W2(J) - XTEMP*W2(J - LEN+K-1)
            W3(J) = W3(J) - XTEMP*W3(J - LEN+K-1)
         END DO
         XTEMP = COLVAL(LD(J))
         W1(J) = W1(J)/XTEMP
         W2(J) = W2(J)/XTEMP
         W3(J) = W3(J)/XTEMP
C
500   CONTINUE
C
C-----DIAGONAL SCALING
C
      DO 600 K = 1, NEQ
        XTEMP = COLVAL(LD(K))
        W1(K) = W1(K)*XTEMP
        W2(K) = W2(K)*XTEMP
        W3(K) = W3(K)*XTEMP
600   CONTINUE
C
C-----BACKWARD SUBSTITUTE
C
      DO 700 K = NEQ, 2, -1
        XTEMP = COLVAL(LD(K))
        W1(K) = W1(K)/XTEMP
        W2(K) = W2(K)/XTEMP
        W3(K) = W3(K)/XTEMP

        I1 = LD(K - 1) + 1

        I2 = LD(K) - 1

        IF(I1.GT.I2) GO TO 700

        IW = K - LD(K) + ( LD(K - 1) + 1 )

C  Break the dependency check between W1,W2,W3
CDIR$ IVDEP
        DO 710 I = I1, I2
            XTEMP     = COLVAL(I)
            W1(IW)    = W1(IW) - XTEMP*W1(K)
            W2(IW)    = W2(IW) - XTEMP*W2(K)
            W3(IW)    = W3(IW) - XTEMP*W3(K)
            IW        = IW + 1
710     CONTINUE
700   CONTINUE
C
      XTEMP = COLVAL(LD(1))
      W1(1) = W1(1)/XTEMP
      W2(1) = W2(1)/XTEMP
      W3(1) = W3(1)/XTEMP
C
      RETURN
      END
C
C=END FORTRAN
C
