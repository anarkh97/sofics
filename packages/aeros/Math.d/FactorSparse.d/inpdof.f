************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    01-14-1999 (Michel Lesoinne and Kendall H. Pierson)
*
*   Modified from the inpnv.f routine.
*
************************************************************************
************************************************************************
*****   INPDOF ... input numerical values                           *****
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
*     INVSUPER  (temporary) integer array dimension N
*               Keep track of dof to super node map
*               
************************************************************************
*
      SUBROUTINE  INPDOF  ( N     , COLPTR, ROWIDX, VALUES, PERM  ,
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
        INTEGER             COLPTR(*), ROWIDX(*)
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
*       --------------
*       Parameters ...
*       --------------
*        DOUBLE PRECISION    ZERO
*        PARAMETER       (   ZERO = 0.0D0 )
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             FSTCOL, II, JJ, IROW  , JLEN  ,
     &                      LASTL , LSTCOL, LXBEG , LXEND ,
     &                      CSUPER, LASTSUPER, P1
*       INTEGER JSUPER, OLDJ
*
************************************************************************
*
*       ------------------------------
*       Initialize the data structure.
*       ------------------------------
*        DO  II = 1, XLNZ(N+1)-1
*            LNZ(II) = ZERO
*        END DO
*
*
*       ----------------------
*       Fill INVSUPER array
*       ----------------------
*        DO  JSUPER = 1, NSUPER
*          DO II=XSUPER(JSUPER), XSUPER(JSUPER+1) - 1
*              INVSUPER(II) = JSUPER
*          END DO
*        END DO
*
***********************************************************
*
        LASTSUPER = -1
        DO II=1, N
          P1 = INVP(II)
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
            DO JJ = COLPTR(II), COLPTR(II+1)-1
                  IROW = INVP(ROWIDX(JJ))
                  IF  ( IROW .GE. FSTCOL )  THEN
                      LNZ(LASTL - OFFSET(IROW)) = VALUES(JJ)
                  END IF
            END DO
          END IF
 
        END DO
*
*       -----------------------------------------------
*       For each supernode JSUPER, do the following ...
*       -----------------------------------------------
*
**        LXBEG  = XLINDX(1)
*
*       The span from FSTCOL to LSTCOL represent the 
*       degrees of freedom that belong to JSUPER node
*
**        FSTCOL = XSUPER(1)
*
**        DO  JSUPER = 1, NSUPER
*           -----------------------------------------------
*           First get offset to facilitate numerical input.
*           -----------------------------------------------
**            LXEND = XLINDX(JSUPER+1)
**            JLEN  = LXEND - LXBEG
**            DO  II = LXBEG, LXEND-1
**                IROW         = LINDX(II)
**                JLEN         = JLEN - 1
**                OFFSET(IROW) = JLEN
**            END DO
*
**            LSTCOL = XSUPER(JSUPER+1)
**            DO  JCOL = FSTCOL, LSTCOL-1
*               -----------------------------------
*               Next input the individual nonzeros.
*               -----------------------------------
**                OLDJ  = PERM(JCOL)
**                LASTL = XLNZ(JCOL+1) - 1
**                DO  II = COLPTR(OLDJ), COLPTR(OLDJ+1)-1
**                    IROW = INVP(ROWIDX(II))
**                    IF  ( IROW .GE. FSTCOL )  THEN
**                        LNZ(LASTL - OFFSET(IROW)) = VALUES(II)
**                    END IF
**                END DO
**            END DO
**            LXBEG  = LXEND
**            FSTCOL = LSTCOL
**        END DO
        RETURN
*
*       -------------
*       End of INPNV.
*       -------------
      END
