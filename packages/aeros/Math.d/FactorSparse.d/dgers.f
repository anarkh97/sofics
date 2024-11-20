      SUBROUTINE  DGERS (
     &              N, NRHS, A, LDA, NDEF, IPROW, IPCOL, B, LDB, INFO
     &          )
*
*  -- LAPACK-like routine --
*     Esmond G. Ng, Oak Ridge National Labroatory
*     March 18, 1998
*
*     .. Scalar Arguments ..
      INTEGER           INFO, LDA, LDB, N, NDEF, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER           IPCOL(*), IPROW(*)
      DOUBLE PRECISION  A(LDA,*), B(LDB,*)
*     ..
*
*  Purpose
*  =======
*
*  DGERS solves a real general linear system using the triangular
*  factorization computed by DGECP.
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
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The triangular factors L and U from the factorization of
*          A = P*L*U*Q', as computed by DGECP.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  NDEF    (input) INTEGER
*          The rank deficiency of the matrix A.  NDEF <= N.
*
*  IPROW   (input) INTEGER array, dimension N.
*          The row permutation P.
*
*  IPCOL   (input) INTEGER array, dimension N.
*          The column permutation Q.
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*       .. Parameters ..
        DOUBLE PRECISION    ONE, ZERO
        PARAMETER           ( ONE = 1.0D+0 )
        PARAMETER           ( ZERO = 0.0D+0 )
*       ..
*       .. Local Scalars ..
        INTEGER             J, K, LEN
*       ..
*       .. External Subroutines ..
        EXTERNAL            DLASWP, DTRSM, XERBLA
*       ..
*       .. Intrinsic Functions ..
        INTRINSIC           MAX
*       ..
*       .. Executable Statements ..
*
*       Test the input parameters.
*
        INFO = 0
        IF  ( N .LT. 0 )  THEN
            INFO = -1
        ELSE IF  ( NRHS .LT. 0 )  THEN
            INFO = -2
        ELSE IF  ( LDA .LT. MAX(1,N) )  THEN
            INFO = -4
        ELSE IF  ( NDEF .LT. 0 )  THEN
            INFO = -5
        ELSE IF  ( LDB .LT. 0 )  THEN
            INFO = -9
        END IF
        IF  ( INFO .NE. 0 )  THEN
            CALL  XERBLA ( 'DGERS', -INFO )
            RETURN
        END IF
*
*       Quick return if possible
*
        IF  ( N .EQ. 0 )  RETURN
*
*       Apply row permutation to the right hand side matrix ...
*
        CALL  DLASWP ( NRHS, B, LDB, 1, N, IPROW, 1 )
*
*       Forward substitution ...
*
        CALL  DTRSM ( 'Left', 'Lower', 'No Transpose', 'Unit',
     &                N-NDEF, NRHS, ONE, A, LDA, B, LDB )
*
*       Handle rank-deficiency system ...
*
        IF  ( NDEF .GT. 0 )  THEN
            LEN = N - NDEF + 1
            DO  J = 1, NRHS
                DO  K = LEN, N
                    B(K,J) = ZERO
                END DO
            END DO
        END IF
*
*       Backward substitution ...
*
        CALL  DTRSM ( 'Left', 'Upper', 'No Transpose', 'Non-Unit',
     &                N-NDEF, NRHS, ONE, A, LDA, B, LDB )
*
*       Apply column permutation to the solution matrix ...
*
        CALL  DLASWP ( NRHS, B, LDB, 1, N, IPCOL, -1 )
*
        RETURN
*
*       End of DGERS
*
      END
