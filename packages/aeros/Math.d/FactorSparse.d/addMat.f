************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    Esmond G. Ng and Barry W. Peyton
*   Created:    12-27-1994 (Barry W. Peyton)
*   Modified:   09-23-1998 (Esmond G. Ng)
*   Modified:   01-14-1999 (Michel Lesoinne and Kendall H. Pierson)
*
*   Computational Mathematics and Statistics Section
*   Oak Ridge National Laboratory
*
************************************************************************
************************************************************************
*****   ADDMAT ... input numerical values                           *****
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
*     INVSUPER  (temporary) integer array dimension N
*               Keep track of dof to super node map
*               
************************************************************************
*
      SUBROUTINE  ADDMAT  ( N     , ROWIDX, VALUES, PERM  ,
     &                      INVP  , NSUPER, XSUPER, XLINDX, LINDX ,
     &                      XLNZ  , LNZ   , OFFSET, INVSUPER       )
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
        INTEGER             ROWIDX(*)
        DOUBLE PRECISION    VALUES(*)
        INTEGER             PERM(*), INVP(*)
        INTEGER             XSUPER(*)
        INTEGER             XLINDX(*), LINDX(*)
        INTEGER             XLNZ(*)
        DOUBLE PRECISION    LNZ(*)
        INTEGER             OFFSET(*)
        INTEGER             INVSUPER(*)
*
************************************************************************
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             FSTCOL, II, JJ, IROW  , JLEN  ,
     &                      LASTL , LSTCOL, LXBEG , LXEND ,
     &                      CSUPER, LASTSUPER, P1
*
************************************************************************
*
        LASTSUPER = -1
        DO II=1, N
          P1 = INVP(ROWIDX(II))
          CSUPER = INVSUPER(P1)

          IF(CSUPER.NE.LASTSUPER) THEN
*
*           -----------------------------------------------
*           First get offset to facilitate numerical input.
*           -----------------------------------------------
*
            LXBEG = XLINDX(CSUPER)
            LXEND = XLINDX(CSUPER+1)
            JLEN  = LXEND - LXBEG
            DO  JJ = LXBEG, LXEND-1
                IROW         = LINDX(JJ)
                JLEN         = JLEN - 1
                OFFSET(IROW) = JLEN
            END DO
            FSTCOL = XSUPER(CSUPER)
            LSTCOL = XSUPER(CSUPER+1)
            LASTL  = XLNZ(P1+1) - 1
            DO JJ = 1, N
                  IROW = INVP(ROWIDX(JJ))
                  IF  ( IROW .GE. FSTCOL )  THEN
                    LNZ(LASTL - OFFSET(IROW)) = 
     &               LNZ(LASTL - OFFSET(IROW)) + VALUES(JJ)
                  END IF
            END DO
          END IF
 
        END DO
*
*       -------------
*       End of ADDMAT.
*       -------------
      END
