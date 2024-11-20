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
*****   INPNV ... input numerical values                           *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine inputs numerical values of a symmetric matrix
*     into sparse data structures that have been set up for Cholesky
*     factorization.  It is assumed that the input matrix is stored
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
*     PERM      (input) integer array, dimension N
*               The permutation vector.
*
*     INVP      (input) integer array, dimension N
*               The inverse of the permutation vector.
*   
*     NSUPER    (input) integer
*               Number of supernodes.
*
*     XSUPER    (input) integer array, dimension NSUPER+1
*               The supernodal partition.
*
*     XLINDX    (input) integer array, dimension NSUPER+1
*               Pointers to column structure of the Cholesky factor.
*
*     LINDX     (input) integer array, dimension XSUPER(NSUPER+1)-1
*               Row indices of nonzero entries in the Cholesky factor,
*               stored by columns, in a compressed representation.
*
*     XLNZ      (input) integer array, dimension N+1
*               Pointers to nonzero entries in the Cholesky factor.
*   
*     LNZ       (output) double precision array, dimension XLNZ(N+1)-1
*               Numerical values of nonzero entries in the Cholesky
*               factor, stored by columns.
*
*     OFFSET    (temporary) integer array, dimension N
*               Keep track of relative positions of nonzero entries.
*               
************************************************************************
*
      SUBROUTINE  INPNV   ( N     , COLPTR, ROWIDX, VALUES, PERM  ,
     &                      INVP  , NSUPER, XSUPER, XLINDX, LINDX ,
     &                      XLNZ  , LNZ   , OFFSET                  )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             N     , NSUPER
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             COLPTR(*), ROWIDX(*)
        DOUBLE PRECISION    VALUES(*)
        INTEGER             PERM(*), INVP(*)
        INTEGER             XSUPER(*)
        INTEGER             XLINDX(*), LINDX(*)
        INTEGER             XLNZ(*)
        DOUBLE PRECISION    LNZ(*)
        INTEGER             OFFSET(*)
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
        INTEGER             FSTCOL, II    , IROW  , JCOL  , JLEN  ,
     &                      JSUPER, LASTL , LSTCOL, LXBEG , LXEND ,
     &                      OLDJ
*
************************************************************************
*
*       ------------------------------
*       Initialize the data structure.
*       ------------------------------
        DO  II = 1, XLNZ(N+1)-1
            LNZ(II) = ZERO
        END DO
*
*       -----------------------------------------------
*       For each supernode JSUPER, do the following ...
*       -----------------------------------------------
*
        LXBEG  = XLINDX(1)
        FSTCOL = XSUPER(1)
        DO  JSUPER = 1, NSUPER
*           -----------------------------------------------
*           First get offset to facilitate numerical input.
*           -----------------------------------------------
            LXEND = XLINDX(JSUPER+1)
            JLEN  = LXEND - LXBEG
            DO  II = LXBEG, LXEND-1
                IROW         = LINDX(II)
                JLEN         = JLEN - 1
                OFFSET(IROW) = JLEN
            END DO
*
            LSTCOL = XSUPER(JSUPER+1)
            DO  JCOL = FSTCOL, LSTCOL-1
*               -----------------------------------
*               Next input the individual nonzeros.
*               -----------------------------------
                OLDJ  = PERM(JCOL)
                LASTL = XLNZ(JCOL+1) - 1
                DO  II = COLPTR(OLDJ), COLPTR(OLDJ+1)-1
                    IROW = INVP(ROWIDX(II))
                    IF  ( IROW .GE. FSTCOL )  THEN
                        LNZ(LASTL - OFFSET(IROW)) = VALUES(II)
                    END IF
                END DO
            END DO
            LXBEG  = LXEND
            FSTCOL = LSTCOL
        END DO
        RETURN
*
*       -------------
*       End of INPNV.
*       -------------
      END
