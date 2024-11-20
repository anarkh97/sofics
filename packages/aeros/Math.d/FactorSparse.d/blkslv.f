************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    Esmond G. Ng and Barry W. Peyton
*   Created:    05-06-1996 (Barry W. Peyton)
*   Modified:   09-23-1998 (Esmond G. Ng)
*
*   Computational Mathematics and Statistics Section
*   Oak Ridge National Laboratory
*
************************************************************************
************************************************************************
*****   BLKSLV ... block triangular solutions                      *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     Given the Cholesky factorization of a sparse symmetric matrix,
*     this subroutine performs the triangular solution.  It uses
*     output from BLKLDL.
*
*   ----------
*   Arguments:
*   ----------
*
*     NPANEL    (input) integer
*               Number of panels.  (N = XPANEL(NPANEL+1)-1.)
*
*     XPANEL    (input) integer array, dimension NPANEL+1
*               The panel partition.
*
*     XLINDX    (input) integer array, dimension NPANEL+1
*               Pointers to column structure of the Cholesky factor.
*
*     LINDX     (input) integer array, dimension XSUPER(NPANEL+1)-1
*               Row indices of nonzero entries in the Cholesky factor,
*               stored by columns, in a compressed representation.
*
*     XLNZ      (input) integer array, dimension N+1
*               Pointers to nonzero entries in the Cholesky factor.
*
*     LNZ       (input) double precision array, dimension XLNZ(N+1)-1
*               Numerical values of nonzero entries in the Cholesky
*               factor, stored by columns.
*
*     DEFBLK    (input) integer
*               The size of the last supernode, which contains the
*               known deficiency.
*
*     NDEF      (input) integer
*               Rank deficiency of the matrix.
*
*     LBDEF     (input) integer
*               If DEFBLK is nonzero, LBDEF is the deficiency of the
*               last block.
*
*     DEF       (input) integer array, dimension NDEF
*               It identifies the columns of the matrix that are
*               linearly dependent.
*
*     IPROW     (input) integer array, dimension DEFBLK
*               If DEFBLK is nonzero, IPROW contains the row pivoting
*               sequence for the last block.
*
*     IPCOL     (input) integer array, dimension DEFBLK
*               If DEFBLK is nonzero, IPCOL contains the column
*               pivoting sequence for the last block.
*
*     PERM      (input) integer array, dimension N
*               The permutation vector.
*
*     INVP      (input) integer array, dimension N
*               The inverse of the permutation vector.
*
*     RHS       (input) double precision array, dimension N
*               The right hand side vector.
*
*     SOLN      (output) double precision array, dimension N
*               The solution vector.
*
*     TEMP      (temporary) double precision array, dimension N
*               Temporary work array.
*
************************************************************************
*
      SUBROUTINE  BLKSLV  ( NPANEL, XPANEL, XLINDX, LINDX , XLNZ  ,
     &                      LNZ   , DEFBLK, NDEF  , LBDEF , DEF   ,
     &                      IPROW , IPCOL , PERM  , INVP  , RHS   ,
     &                      SOLN  , TEMP                            )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             DEFBLK, LBDEF , NDEF  , NPANEL
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             XPANEL(*)
        INTEGER             XLINDX(*), LINDX(*)
        INTEGER             XLNZ(*)
        DOUBLE PRECISION    LNZ(*)
        INTEGER             DEF(*)
        INTEGER             IPROW(*), IPCOL(*)
        INTEGER             PERM(*), INVP(*)
        DOUBLE PRECISION    RHS(*), SOLN(*)
        DOUBLE PRECISION    TEMP(*)
*
************************************************************************
*
*       --------------
*       Parameters ...
*       --------------
        DOUBLE PRECISION    ONE   , ZERO
        PARAMETER       (   ONE  = 1.0D0,
     &                      ZERO = 0.0D0 )
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             FJ    , INFO  , ISUB  , J     , JDPTR ,
     &                      JJ    , JLEN  , JLPNT , JPANEL, JXPNT ,
     &                      LJ    , LSTPNL, N     , NJ
*
*       ------------------------
*       External Subroutines ...
*       ------------------------
        EXTERNAL            DGEMV , DGERS , DTRSM
*
************************************************************************
*
        IF  ( NPANEL .LE. 0 )  RETURN
*
        N = XPANEL(NPANEL+1) - 1
        DO  J = 1, N
            TEMP(J) = RHS(PERM(J))
            SOLN(J) = ZERO
        END DO
*
*       ------------------------
*       Forward substitution ...
*       ------------------------
        LSTPNL= NPANEL
        IF  ( DEFBLK .NE. 0 )  LSTPNL = NPANEL - 1
*
        DO  JPANEL = 1, LSTPNL
*
            FJ    = XPANEL(JPANEL)
            NJ    = XPANEL(JPANEL+1) - FJ
            JLEN  = XLNZ(FJ+1) - XLNZ(FJ)
            JXPNT = XLINDX(JPANEL)
            JLPNT = XLNZ(FJ)
*
            CALL  DTRSM (
     &              'LEFT', 'LOWER', 'NO TRANSPOSE', 'UNIT',
     &              NJ, 1,
     &              ONE,
     &              LNZ(JLPNT), JLEN,
     &              TEMP(FJ), NJ
     &          )
*
            CALL  DGEMV (
     &              'NO TRANSPOSE',
     &              JLEN-NJ, NJ,
     &              -ONE,
     &              LNZ(JLPNT+NJ), JLEN,
     &              TEMP(FJ), 1,
     &              ZERO,
     &              SOLN, 1
     &          )
*
            JJ = JXPNT + NJ - 1
            DO  J = 1, JLEN-NJ
                JJ         = JJ + 1
                ISUB       = LINDX(JJ)
                TEMP(ISUB) = TEMP(ISUB) + SOLN(J)
                SOLN(J)    = ZERO
            END DO
*
        END DO
*
        IF  ( DEFBLK .NE. 0 )  THEN
            FJ    = XPANEL(NPANEL)
            LJ    = XPANEL(NPANEL+1) - 1
            NJ    = LJ - FJ + 1
            JLEN  = XLNZ(FJ+1) - XLNZ(FJ)
            JLPNT = XLNZ(FJ)
            CALL  DGERS ( NJ, 1, LNZ(JLPNT), JLEN, LBDEF,
     &                    IPROW, IPCOL, TEMP(FJ), NJ, INFO )
        END IF
*
        DO  JJ = 1, NDEF-LBDEF
            J       = DEF(JJ)
            TEMP(J) = ZERO
        END DO
        DO  JPANEL = 1, LSTPNL
            FJ    = XPANEL(JPANEL)
            LJ    = XPANEL(JPANEL+1) - 1
            JLEN  = XLNZ(FJ+1) - XLNZ(FJ)
            JDPTR = XLNZ(FJ)
            DO  J = FJ, LJ
                TEMP(J) = TEMP(J)/LNZ(JDPTR)
                JDPTR   = JDPTR + JLEN + 1
            END DO
        END DO
*
*       -------------------------
*       Backward substitution ...
*       -------------------------
        DO  JPANEL = LSTPNL, 1, -1
*
            FJ    = XPANEL(JPANEL)
            NJ    = XPANEL(JPANEL+1) - FJ
            JLEN  = XLNZ(FJ+1) - XLNZ(FJ)
            JXPNT = XLINDX(JPANEL)
            JLPNT = XLNZ(FJ)
*
            JJ = JXPNT + NJ - 1
********KHP
*            DO  J = 1, JLEN
            DO  J = 1, JLEN - NJ
                JJ      = JJ + 1
                ISUB    = LINDX(JJ)
                SOLN(J) = TEMP(ISUB)
            END DO
*
            IF  ( JLEN .GT. NJ )  THEN
                CALL  DGEMV (
     &                  'TRANSPOSE',
     &                  JLEN-NJ, NJ,
     &                  -ONE,
     &                  LNZ(JLPNT+NJ), JLEN,
     &                  SOLN, 1,
     &                  ONE,
     &                  TEMP(FJ), 1
     &              )
            END IF
*
            CALL  DTRSM (
     &              'LEFT', 'LOWER', 'TRANSPOSE', 'UNIT',
     &              NJ, 1,
     &              ONE,
     &              LNZ(JLPNT), JLEN,
     &              TEMP(FJ), NJ
     &          )
*
        END DO
*
        DO  J = 1, N
            SOLN(J) = TEMP(INVP(J))
        END DO
*
        RETURN
*
*       --------------
*       End of BLKSLV.
*       --------------
      END
