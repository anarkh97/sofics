************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    Esmond G. Ng and Barry W. Peyton
*   Created:    12-27-1994 (Barry W. Peyton)
*   Modified:   09-23-1998 (Esmond G. Ng)
*
*   Computational Mathematics and Statistics Section
*   Oak Ridge National Laboratory
*
************************************************************************
************************************************************************
*****   GETRHS ... generate right hand side vector                 *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine computes the right hand side vector from a known
*     solution.  The input matrix is assumed to be stored by columns.
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
*     SOLN      (input) double precision array, dimension N
*               The solution vector.
*
*     RHS       (output) double precision array, dimension N
*               The right hand side vector.
*
************************************************************************
*
      SUBROUTINE  GETRHS  ( N     , COLPTR, ROWIDX, VALUES, SOLN  ,
     &                      RHS                                     )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             N
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             COLPTR(*), ROWIDX(*)
        DOUBLE PRECISION    VALUES(*)
        DOUBLE PRECISION    RHS(*)   , SOLN(*)
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
        DOUBLE PRECISION    S
*
************************************************************************
*
*       ---------------
*       Initialization.
*       ---------------
        DO  J = 1, N
            RHS(J) = ZERO
        END DO
*
        COLBEG = COLPTR(1)
        DO  J = 1, N
*           -------------------
*           For each column ...
*           -------------------
            COLEND = COLPTR(J+1)
            S      = SOLN(J)
            DO  II = COLBEG, COLEND-1
                I      = ROWIDX(II)
                RHS(I) = RHS(I) + S*VALUES(II)
            END DO
            COLBEG = COLEND
        END DO
        RETURN
*
*       --------------
*       End of GETRHS.
*       --------------
      END
