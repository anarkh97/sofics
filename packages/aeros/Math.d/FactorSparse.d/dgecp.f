      SUBROUTINE  DGECP (
     &              N, A, LDA, TOL, ASDEF, NDEF, IPROW, IPCOL, INFO
     &          )
*
*  -- LAPACK-like routine --
*     Esmond G. Ng, Oak Ridge National Labroatory
*     March 18, 1998
*
*     .. Scalar Arguments ..
      INTEGER           INFO, LDA, N, NDEF, ASDEF
      DOUBLE PRECISION  TOL
*     ..
*     .. Array Arguments ..
      INTEGER           IPCOL(*), IPROW(*)
      DOUBLE PRECISION  A(LDA,*)
*     ..
*
*  Purpose
*  =======
*
*  DGECP computes the triangular factorization of a real general matrix
*  A using Gaussian elimination with complete pivoting.
*
*  The factorization has the form
*     A = P* L * U * Q' ,
*  where P and Q are permutation matrices, L is a unit lower triangular
*  matrix, and U is upper triangular,
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the matrix A.
*
*          On exit, if INFO = 0, the factors L and U from the triangular
*          factorization A = P*L*U*Q.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  TOL     (input) DOUBLE PRECISION
*          Tolerance for checking if a pivot is zero.
*
*  ASDEF   (input) the assumed deficiency of the matrix A
*
*  NDEF    (output) INTEGER
*          The additional deficiency of the matrix A.  NDEF <= N-ASDEF.
*
*  IPROW   (output) INTEGER array, dimension N.
*          The row permutation.
*
*  IPCOL   (output) INTEGER array, dimension N.
*          The column permutation.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*       .. Parameters ..
        DOUBLE PRECISION    ONE, ZERO
        PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*       ..
*       .. Local Scalars ..
        INTEGER             I, IROW, J, JCOL, K, LEN
        DOUBLE PRECISION    PIVOT, T
*       ..
*       .. External Functions ..
        INTEGER             IDAMAX
        EXTERNAL            IDAMAX
*       ..
*       .. External Subroutines ..
        EXTERNAL            DGER, DSCAL, DSWAP, XERBLA
*       ..
*       .. Intrinsic Functions ..
        INTRINSIC           ABS, MAX
*       ..
*       .. Executable Statements ..
*
*       Test the input parameters.
*
        INFO = 0
        IF  ( N .LT. 0 )  THEN
            INFO = -1
        ELSE IF  ( LDA .LT. MAX(1,N) )  THEN
            INFO = -3
        END IF
        IF  ( INFO .NE. 0 )  THEN
            CALL  XERBLA ( 'DGECP', -INFO )
            RETURN
        END IF
*
*       Quick return if possible
*
        IF  ( N .EQ. 0 )  RETURN
*
*       Gaussian elimination with complete pivoting ...
*
        DO  K = 1, N-ASDEF
*
*           Find pivot ...
*
            IROW = K
            JCOL = K
            PIVOT = ABS(A(K,K))
*
            LEN = N - K + 1
            DO  J = K, N
                I = (K - 1) + IDAMAX( LEN, A(K,J), 1 )
                T = ABS(A(I,J))
                IF  ( T .GT. PIVOT )  THEN
                    PIVOT = T
                    IROW = I
                    JCOL = J
                END IF
            END DO
*
*           Check for singularity ...
*
            IF  ( PIVOT .LE. TOL )  THEN
                NDEF = LEN
                DO  J = K, N
                    DO  I = K, N
                        A(I,J) = ZERO
                    END DO
                    A(J,J) = ONE
                    IPCOL(J) = J
                    IPROW(J) = J
                END DO
c                WRITE(6,*) "Exiting at the wrong place"
                RETURN
            END IF
*
            IPROW(K) = IROW
            IPCOL(K) = JCOL
*
*           Interchange rows ...
*
            IF  ( K .NE. IROW )  THEN
                CALL  DSWAP ( N, A(K,1), LDA, A(IROW,1), LDA )
            END IF
*
*           Interchange columns ...
*
            IF  ( K .NE. JCOL )  THEN
                CALL  DSWAP ( N, A(1,K), 1, A(1,JCOL), 1 )
            END IF
*
            IF  ( K .LT. N )  THEN
*
*               Compute multipliers ...
*
                CALL  DSCAL ( N-K, ONE/A(K,K), A(K+1,K), 1 )
*
*               Apply rank-1 update ...
*
                CALL  DGER ( N-K, N-K, -ONE, A(K+1,K), 1, A(K,K+1), LDA,
     &                       A(K+1,K+1), LDA )
            END IF
*
        END DO
*       Make identity with the assumed singular part of the matrix
        NDEF = ASDEF
        DO  J = K, N
            DO  I = K, N
                A(I,J) = ZERO
            END DO
            A(J,J) = ONE
            IPCOL(J) = J
            IPROW(J) = J
        END DO

        RETURN
*
*       End of DGECP
*
      END
