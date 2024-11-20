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
*****   GETNRM ... compute the norm of the input matrix            *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutines computes the L-infinity norm of the input
*     (symmetric) matrix.  It is assumed that the matrix is stored
*     by columns.
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
*     TEMP      (temporary) double precision array, dimension N
*               TEMP(I) is the L-1 norm of row I of the input matrix.
*
*     ANORM     (output) double precision
*               L-infinity norm of the input matrix.
*           
************************************************************************
*
      SUBROUTINE  GETNRM  ( N     , COLPTR, ROWIDX, VALUES, TEMP  ,
     &                      ANORM                                   )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             N
        DOUBLE PRECISION    ANORM
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             COLPTR(*) , ROWIDX(*)
        DOUBLE PRECISION    VALUES(*)
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
        INTEGER             COLBEG, COLEND, J
        DOUBLE PRECISION    T
*
*       ----------------------
*       External Functions ...
*       ----------------------
        DOUBLE PRECISION    DASUM
        EXTERNAL            DASUM
*
*       -----------------------
*       Intrinsic Functions ...
*       -----------------------
        INTRINSIC           MAX
*
************************************************************************
*
        ANORM  = ZERO
        COLBEG = COLPTR(1)
        DO  J = 1, N
*           -------------------
*           For each column ...
*           -------------------
            COLEND = COLPTR(J+1)
            T      = DASUM( COLEND-COLBEG, VALUES(COLBEG), 1 )
            ANORM  = MAX(ANORM,T)
            COLBEG = COLEND
        END DO
        RETURN
*
*       --------------
*       End of GETNRM.
*       --------------
      END
