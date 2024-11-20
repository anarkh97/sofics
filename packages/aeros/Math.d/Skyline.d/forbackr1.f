C=MODULE  FORBACKR1
C=PURPOSE Forward/backward substitution of a real symmetric positive
C=PURPOSE semi-definite system decomposed and stored in skyline form
C=PURPOSE and using variable band information
C=PURPOSE FORWARD and BACKWARD are serial
C=AUTHOR  K. H. PIERSON
C=AUTHOR  June 1998
C=BLOCK
      subroutine FORBACKR1CF(COLVAL, LD, B, PIVOT, NEQ)
C
C----------------------------------------------------------------------------
C
C     FORBACKR1 = forward/backward substitution of a real symmetric skyline
C                 matrix using variable band information
C                 with unrolling depth = 1
C
C     Version 1.1 (real)
C     =====================
C
C     COLVAL   (real)    Stores sequentially values in column profile of K
C     LD       (integer) Locates diagonal elements of K in COLVAL
C     B        (real)    RHS vector
C     PIVOT    (integer) stores singularity information
C     NEQ      (integer) Number of equations in system
C
C----------------------------------------------------------------------------
C
C-----DECLARATIONS
C
      INTEGER LD(*), PIVOT(*), NEQ
      REAL*8  COLVAL(*), B(*), XTEMP
C
C-----LOCAL DECLARATIONS
C
      INTEGER K,J,I,IW,LEN,I1,I2
      REAL*8  CHDDOT
C
C-----BEGIN
C
      IF(PIVOT(1) .EQ. 0) THEN
         B(1) = 0.0
      ELSE
         B(1) =  B(1) / COLVAL(LD(1))
      ENDIF
C
C-----FORWARD SUBSTITUTE
C
      DO 500 J = 2, NEQ
        IF(PIVOT(J) .EQ. 0) THEN
           B(J) = 0.0
        ELSE
           I1   = ( LD(J - 1) + 1 )
           I2   = LD(J)
           LEN  = ( I2 - I1 )
           B(J) = (B(J)-CHDDOT(LEN,COLVAL(I1),B(J - LEN)))/COLVAL(I2)
        ENDIF
500   CONTINUE
C
C-----DIAGONAL SCALING
C
      DO 600 K = 1, NEQ
        B(K) = B(K)*COLVAL(LD(K))
600   CONTINUE
C
C-----BACKWARD SUBSTITUTE
C
      DO 700 K = NEQ, 2, -1

         IF(B(K) .EQ. 0.0) GO TO 700

         IF(PIVOT(K) .EQ. 0) THEN
           B(K) = 0.0
           GO TO 700
         ELSE
           B(K) = B(K) / COLVAL(LD(K))
         END IF

         I1 = LD(K - 1) + 1

         I2 = LD(K) - 1

         IF(I1 .GT. I2) GO TO 700

         IW = K - LD(K) + ( LD(K - 1) + 1 )

         XTEMP = B(K)
         DO 710 I = I1, I2
           B(IW) = B(IW) - COLVAL(I)*XTEMP
           IW    = IW + 1
710      CONTINUE

700   CONTINUE

      IF(PIVOT(1) .EQ. 0) THEN
         B(1) = 0.0
      ELSE
         B(1) = B(1)/COLVAL(LD(1))
      END IF
C
      RETURN
      END
C
C
      subroutine FORR1(COLVAL, LD, B, PIVOT, NEQ)
C
C----------------------------------------------------------------------------
C
C     FORRR1 = forward substitution of a real symmetric skyline
C              matrix using variable band information
C              with unrolling depth = 1
C
C     Version 1.1 (real)
C     =====================
C
C     COLVAL   (real)    Stores sequentially values in column profile of K
C     LD       (integer) Locates diagonal elements of K in COLVAL
C     B        (real)    RHS vector
C     PIVOT    (integer) stores singularity information
C     NEQ      (integer) Number of equations in system
C
C----------------------------------------------------------------------------
C
C-----DECLARATIONS
C
      INTEGER LD(*), PIVOT(*), NEQ
      REAL*8  COLVAL(*), B(*)
C
C-----LOCAL DECLARATIONS
C
      INTEGER K,J,LEN,I1,I2
      REAL*8  CHDDOT
C
C-----BEGIN
C
      IF(PIVOT(1) .EQ. 0) THEN
         B(1) = 0.0
      ELSE
         B(1) =  B(1) / COLVAL(LD(1))
      ENDIF
C
C-----FORWARD SUBSTITUTE
C
      DO 500 J = 2, NEQ
        IF(PIVOT(J) .EQ. 0) THEN
           B(J) = 0.0
        ELSE
           I1   = ( LD(J - 1) + 1 )
           I2   = LD(J)
           LEN  = ( I2 - I1 )
           B(J) = (B(J)-CHDDOT(LEN,COLVAL(I1),B(J - LEN)))/COLVAL(I2)
        ENDIF
500   CONTINUE
C
C-----DIAGONAL SCALING
C
      DO 600 K = 1, NEQ
        B(K) = B(K)*SQRT(COLVAL(LD(K)))
c        B(K) = B(K)*COLVAL(LD(K))
600   CONTINUE
C
      RETURN
      END
C
C
      subroutine BACKR1(COLVAL, LD, B, PIVOT, NEQ)
C
C----------------------------------------------------------------------------
C
C     BACKR1 = backward substitution of a real symmetric skyline
C              matrix using variable band information
C              with unrolling depth = 1
C
C     Version 1.1 (real)
C     =====================
C
C     COLVAL   (real)    Stores sequentially values in column profile of K
C     LD       (integer) Locates diagonal elements of K in COLVAL
C     B        (real)    RHS vector
C     PIVOT    (integer) stores singularity information
C     NEQ      (integer) Number of equations in system
C
C----------------------------------------------------------------------------
C
C-----DECLARATIONS
C
      INTEGER LD(*), PIVOT(*), NEQ
      REAL*8  COLVAL(*), B(*), XTEMP
C
C-----LOCAL DECLARATIONS
C
      INTEGER K,I,IW,I1,I2
C
C-----DIAGONAL SCALING
C
      DO 600 K = 1, NEQ
        B(K) = B(K)*SQRT(COLVAL(LD(K)))
c        B(K) = B(K)*COLVAL(LD(K))
600   CONTINUE
C
C-----BACKWARD SUBSTITUTE
C
      DO 700 K = NEQ, 2, -1

         IF(B(K) .EQ. 0.0) GO TO 700

         IF(PIVOT(K) .EQ. 0) THEN
           B(K) = 0.0
           GO TO 700
         ELSE
           B(K) = B(K) / COLVAL(LD(K))
         END IF

         I1 = LD(K - 1) + 1

         I2 = LD(K) - 1

         IF(I1 .GT. I2) GO TO 700

         IW = K - LD(K) + ( LD(K - 1) + 1 )

         XTEMP = B(K)
         DO 710 I = I1, I2
           B(IW) = B(IW) - COLVAL(I)*XTEMP
           IW    = IW + 1
710      CONTINUE

700   CONTINUE

      IF(PIVOT(1) .EQ. 0) THEN
         B(1) = 0.0
      ELSE
         B(1) = B(1)/COLVAL(LD(1))
      END IF
C
      RETURN
      END

