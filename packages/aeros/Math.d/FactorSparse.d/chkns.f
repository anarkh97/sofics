************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    Esmond G. Ng and Barry W. Peyton
*   Created:    09-23-1998 (Esmond G. Ng)
*
*   Computational Mathematics and Statistics Section
*   Oak Ridge National Laboratory
*
************************************************************************
************************************************************************
*****   CHKNS ... check validity of null space                     *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine performs a trivial check on the validity of
*     the computed null space.  The approach is to form the product
*     of the coefficient matrix and its null space, and look for the
*     component with the largest magnitude.
*
*   ----------
*   Arguments:
*   ----------
*
*     N         (input) integer
*               Number of equations.
*
*     COLPTR    (input) integer array, dimension N+1
*               Pointers to the column structure of the input matrix.
*
*     ROWIDX    (input) integer array, dimension COLPTR(N+1)-1
*               Row indices of nonzero entries in the input matrix,
*               stored by columns.
*
*     VALUES    (input) double precision array, dimension COLPTR(N+1)-1
*               Numerical values of nonzero entries in the input
*               matrix, stored by columns.
*
*     NDEF      (input) integer
*               Rank deficiency of the matrix.
*
*     NS        (input) double precision array, dimension (N,NDEF)
*               The null space of the matrix.
*
*     LDNS      (input) integer
*               Leading dimension of NS.
*
*     ERROR     (output) double precision
*               The absolute value of the component in the product
*               of the input matrix and its null space that has the
*               largest magnitude.
*
*     TEMP      (temporary) double precision array, dimension N
*               It holds the product of the input matrix and a column
*               in the null space.
*
************************************************************************
*
      SUBROUTINE  CHKNS   ( N     , COLPTR, ROWIDX, VALUES, NDEF  ,
     &                      NS    , LDNS  , ERROR , TEMP            )
*
************************************************************************
*
*       -----------------
*       Scalar Arguments:
*       -----------------
        INTEGER             LDNS  , N     , NDEF
        DOUBLE PRECISION    ERROR
*
*       ----------------
*       Array Arguments:
*       ----------------
        INTEGER             COLPTR(*), ROWIDX(*)
        DOUBLE PRECISION    VALUES(*)
        DOUBLE PRECISION    NS(LDNS,*)
        DOUBLE PRECISION    TEMP(*)
*
************************************************************************
*
*       -------------
*       Parameter ...
*       -------------
        DOUBLE PRECISION    ZERO
        PARAMETER       (   ZERO = 0.0D0 )
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             COL   , COLBEG, COLEND, I     , II    ,
     &                      J
        DOUBLE PRECISION    S
*
*       -----------------------
*       Intrinsic Functions ...
*       -----------------------
        INTRINSIC           ABS   , MAX
*
************************************************************************
*
        ERROR = ZERO
        DO  J = 1, NDEF
*           ----------------------------------
*           For each column of the null space,
*           form the product.
*           ----------------------------------
            DO  I = 1, N
                TEMP(I) = ZERO
            END DO
            COLBEG = COLPTR(1)
            DO  COL = 1, N
                COLEND = COLPTR(COL+1)
                S      = NS(COL,J)
                DO  II = COLBEG, COLEND-1
                    I       = ROWIDX(II)
                    TEMP(I) = TEMP(I) + S*VALUES(II)
                END DO
                COLBEG = COLEND
            END DO
*           ----------------------------
*           Determine the max magnitude.
*           ----------------------------
            DO  I = 1, N
                ERROR = MAX(ERROR,ABS(TEMP(I)))
            END DO
        END DO
        RETURN
*
*       -------------
*       End of CHKNS.
*       -------------
      END
