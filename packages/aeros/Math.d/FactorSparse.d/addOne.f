*
*   Version:    0.5
*   Modified:   01-14-1999 (Michel Lesoinne and Kendall H. Pierson)
*
*   Computational Mathematics and Statistics Section
*   Oak Ridge National Laboratory
*
************************************************************************
************************************************************************
*****   ADDONE ... input numerical values                           *****
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
*     ROWIDX    (input) integer, row index of input value
*
*     COLIDX    (input) integer, column index of input value
*
*     VALUE     (input) double precision value,
*               Numerical value of a nonzero entry in the input
*               matrix, stored by columns.
*
*     INVP      (input) integer array, dimension N
*               The inverse of the permutation vector.
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
*     INVSUPER  (temporary) integer array dimension N
*               Keep track of dof to super node map
*               
************************************************************************
*
      SUBROUTINE  ADDONE  ( ROWIDX, COLIDX, VALUE,
     &                      INVP  , XSUPER, XLINDX, LINDX ,
     &                      XLNZ  , LNZ   , INVSUPER       )
*
************************************************************************
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             ROWIDX, COLIDX
        DOUBLE PRECISION    VALUE
        INTEGER             INVP(*)
        INTEGER             XSUPER(*)
        INTEGER             XLINDX(*), LINDX(*)
        INTEGER             XLNZ(*)
        DOUBLE PRECISION    LNZ(*)
        INTEGER             INVSUPER(*)
*
************************************************************************
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             FSTCOL, JJ, IROW, LASTL, 
     &                      LXBEG , LXEND, CSUPER, OFFSET, P1
*
************************************************************************
*
        P1     = INVP(ROWIDX)
        CSUPER = INVSUPER(P1)
        FSTCOL = XSUPER(CSUPER)
        IROW   = INVP(COLIDX)
        IF ( IROW .GE. FSTCOL )  THEN
*
*           -----------------------------------------------
*           First get offset to facilitate numerical input.
*           -----------------------------------------------
*
          LXBEG  = XLINDX(CSUPER)
          LXEND  = XLINDX(CSUPER+1)
          OFFSET = LXEND - LXBEG
*
          DO JJ = LXBEG, LXEND-1
             OFFSET        = OFFSET - 1
             IF(LINDX(JJ) .EQ. IROW) GO TO 100
          END DO
100       LASTL  = XLNZ(P1+1) - 1 - OFFSET
*
          LNZ(LASTL) = LNZ(LASTL) + VALUE
*
        END IF
*
*       -------------
*       End of ADDONE.
*       -------------
      END
