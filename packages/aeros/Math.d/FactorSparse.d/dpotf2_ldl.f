      SUBROUTINE  DPOTF2_LDL (
     &              UPLO, N, A, LDA, NDEF, IDEF, TOL, WORK, INFO
     &          )
*
*  -- LAPACK-like routine --
*     Esmond G. Ng, Oak Ridge National Laboratory
*     March 17, 1998
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N, NDEF
      DOUBLE PRECISION   TOL
*     ..
*     .. Array Arguments ..
      INTEGER            IDEF(*)
      DOUBLE PRECISION   A(LDA,*), WORK(*)
*     ..
*
*  Purpose
*  =======
*
*  DPOTF2_LDL computes the Cholesky factorization of a real symmetric
*  positive definite or semi-definite matrix A.
*
*  The factorization has the form
*     A = U' * D * U ,  if UPLO = 'U', or
*     A = L  * D * L',  if UPLO = 'L',
*  where U is a unit upper triangular matrix, L is unit lower
*  triangular, and D is diagonal.
*
*  This is the unblocked version of the algorithm, calling Level 2 BLAS.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          symmetric matrix A is stored.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          n by n upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n by n lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, if INFO = 0, the factor U+D or L+D from the Cholesky
*          factorization A = U'*D*U  or A = L*D*L'.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  NDEF    (output) INTEGER
*          The rank deficiency of the matrix A.  NDEF <= N.
*
*  IDEF    (output) INTEGER array, dimension NDEF
*          Indices of columns for which zero pivots are encountered.
*
*  TOL     (input) DOUBLE PRECISION
*          Tolerance for checking if a pivot is zero.
*
*  WORK    (input) DOUBLE PRECISION array, dimension N.
*          Temporary work array.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, the leading minor of order k is not
*               positive definite or semi-definite, and the
*               factorization could not be completed.
*
*  =====================================================================
*
*       .. Parameters ..
        DOUBLE PRECISION   ONE, ZERO
        PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*       ..
*       .. Local Scalars ..
        LOGICAL            UPPER
        INTEGER            I, J
        DOUBLE PRECISION   AJJ
*       ..
*       .. External Functions ..
        LOGICAL            LSAME
        DOUBLE PRECISION   DDOT
        EXTERNAL           LSAME, DDOT
*       ..
*       .. External Subroutines ..
        EXTERNAL           DGEMV, DSCAL, XERBLA
*       ..
*       .. Intrinsic Functions ..
        INTRINSIC          ABS, MAX
*       ..
*       .. Executable Statements ..
*
*       Test the input parameters.
*
        INFO = 0
        UPPER = LSAME(UPLO,'U')
        IF  ( .NOT. UPPER  .AND.  .NOT. LSAME(UPLO,'L') )  THEN
            INFO = -1
        ELSE IF  ( N .LT. 0 )  THEN
            INFO = -2
        ELSE IF  ( LDA .LT. MAX(1,N) )  THEN
            INFO = -4
        END IF
        IF  ( INFO .NE. 0 )  THEN
            CALL  XERBLA ( 'DPOTF2_LDL', -INFO )
            RETURN
        END IF
*
*       Quick return if possible
*
        IF  ( N .EQ. 0 )  RETURN
*
        IF  ( UPPER )  THEN
*
*           Compute the Cholesky factorization A = U'*D*U.
*
            NDEF = 0
            DO  J = 1, N
*
*               Compute U(J,J) and test for non-positive-definiteness.
*
                CALL  DCOPY ( J-1, A(1,J), 1, WORK, 1 )
                DO  I = 1, J-1
                    WORK(I) = WORK(I)*A(I,I)
                END DO
                AJJ = A(J,J) - DDOT( J-1, A(1,J), 1, WORK, 1 )
                IF   (ABS( AJJ) .GT. TOL )  THEN
                    A(J,J) = AJJ
*
*                   Compute elements J+1:N of row J.
*
                    IF  ( J .LT. N )  THEN
                        CALL  DGEMV ( 'Transpose', J-1, N-J,
     &                                -ONE, A(1,J+1), LDA, WORK, 1,
     &                                ONE, A(J,J+1), LDA )
                        CALL  DSCAL ( N-J, ONE/AJJ, A(J,J+1), LDA )
                    END IF
                ELSE IF  ( ABS(AJJ) .LE. TOL )  THEN
*
*                   Handle zero pivots ...
*
                    NDEF = NDEF + 1
                    IDEF(NDEF) = J
                    A(J,J) = ONE
                    DO  I = J+1, N
                        A(J,I) = ZERO
                    END DO
                ELSE
                    A(J,J) = AJJ
                    INFO = J
                    RETURN
                END IF
            END DO
*
        ELSE
*
*           Compute the Cholesky factorization A = L*D*L'.
*
            NDEF = 0
            DO  J = 1, N
*
*               Compute L(J,J) and test for non-positive-definiteness.
*
                CALL  DCOPY ( J-1, A(J,1), LDA, WORK, 1 )
                DO  I = 1, J-1
                    WORK(I) = WORK(I) * A(I,I)
                END DO
                AJJ = A(J,J) - DDOT( J-1, A(J,1), LDA, WORK, 1 )
                IF  ( ABS(AJJ) .GT. TOL )  THEN
                    A(J,J) = AJJ
*
*                   Compute elements J+1:N of column J.
*
                    IF  ( J .LT. N )  THEN
                        CALL  DGEMV ( 'No transpose', N-J, J-1,
     &                                -ONE, A(J+1,1), LDA, WORK, 1,
     &                                ONE, A(J+1,J), 1 )
                        CALL  DSCAL ( N-J, ONE/AJJ, A(J+1,J), 1 )
                    END IF
                ELSE IF  ( ABS(AJJ) .LE. TOL )  THEN
*
*                   Handle zero pivots ...
*
                    NDEF = NDEF + 1
                    IDEF(NDEF) = J
                    A(J,J) = ONE
                    DO  I = J+1, N
                        A(I,J) = ZERO
                    END DO
                ELSE
                    A(J,J) = AJJ
                    INFO = J
                    RETURN
                END IF
            END DO
        END IF
        RETURN
*
*       End of DPOTF2_LDL
*
      END
