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
*****   GETRES ... compute residual in the linear system           *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine computes the residual in a linear system.  It is
*     assumed that the input matrix is stored by columns.
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
*     RHS       (input) double precision array, dimension N
*               The right hand side vector.
*
*     SOLN      (input) double precision array, dimension N
*               The solution vector.
*
*     RESID     (output) double precision
*               The L-infinity norm of the residual.
*
*     TEMP      (temporary) double precision array, dimension N
*               The residual vector.
*
************************************************************************
*
      SUBROUTINE  GETRES  ( N     , COLPTR, ROWIDX, VALUES, RHS   ,
     &                      SOLN  , RESID , TEMP                    )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             N
        DOUBLE PRECISION    RESID
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             COLPTR(*), ROWIDX(*)
        DOUBLE PRECISION    VALUES(*)
        DOUBLE PRECISION    RHS(*), SOLN(*)
        DOUBLE PRECISION    TEMP(*)
*
************************************************************************
*
*       --------------
*       Parameters ...
*       --------------
        DOUBLE PRECISION    ZERO
        PARAMETER       (   ZERO = 0.0D0 )
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             COLBEG, COLEND, I     , II    , J
        DOUBLE PRECISION    RMAX  , S
*
*       -----------------------
*       Intrinsic Functions ...
*       -----------------------
        INTRINSIC           ABS   , MAX
*
************************************************************************
*
*       ---------------
*       Initialization.
*       ---------------
        DO  J = 1, N
            TEMP(J) = ZERO
        END DO
*
        COLBEG = COLPTR(1)
        DO  J = 1, N
*           ------------------------------------------------
*           Form the matrix-vector product column by column.
*           ------------------------------------------------
            COLEND = COLPTR(J+1)
            S      = SOLN(J)
            DO  II = COLBEG, COLEND-1
                I       = ROWIDX(II)
                TEMP(I) = TEMP(I) + S*VALUES(II)
            END DO
            COLBEG = COLEND
        END DO
*       ------------------------------------------------------
*       Determine the relative residual (using infinity norm).
*       ------------------------------------------------------
        RMAX  = ZERO
        RESID = ZERO
        DO  J = 1, N
            TEMP(J) = RHS(J) - TEMP(J)
            RMAX    = MAX(RMAX,ABS(RHS(J)))
            RESID   = MAX(RESID,ABS(TEMP(J)))
        END DO
        RESID = RESID/RMAX
        RETURN
*
*       --------------
*       End of GETRES.
*       --------------
      END
