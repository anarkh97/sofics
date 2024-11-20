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
*****  ZBLKSLVP ... block triangular solutions (complex)           *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     Given the Cholesky factorization of a sparse symmetric matrix,
*     this subroutine performs the triangular solution for p vectors.
*     It uses the ouput from ZBLKLDL.
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
*     LNZ       (input) COMPLEX*16 array, dimension XLNZ(N+1)-1
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
*     NRHS      (input) integer
*               The number of right hand side vectors
*
*     R         (input) COMPLEX*16 array, dimension N
*               The right hand side vectors.
*
*     S         (output) COMPLEX*16 array, dimension N
*               The solution vectors.
*
*     T         (temporary) COMPLEX*16 array, dimension N
*               Temporary work arrays.
*
************************************************************************
*
      SUBROUTINE  ZBLKSLVP ( NPANEL, XPANEL, XLINDX, LINDX , XLNZ  ,
     &                      LNZ   , DEFBLK, NDEF  , LBDEF , DEF   ,
     &                      IPROW , IPCOL , PERM  , INVP  , NRHS,
     &                      R, LDR, S, T)
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
        COMPLEX*16          LNZ(*)
        INTEGER             DEF(*)
        INTEGER             IPROW(*), IPCOL(*)
        INTEGER             PERM(*), INVP(*)
        INTEGER             NRHS
        INTEGER             LDR
        COMPLEX*16          R(LDR,*)
        COMPLEX*16          S(LDR,*)
        COMPLEX*16          T(LDR,*)
*
************************************************************************
*
*       --------------
*       Parameters ...
*       --------------
        COMPLEX*16         ONE
        PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
        COMPLEX*16         ZERO
        PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             FJ    , INFO  , ISUB  , J     , JDPTR ,
     &                      JJ    , JLEN  , JLPNT , JPANEL, JXPNT ,
     &                      LJ    , LSTPNL, N     , NJ , P
*
*       ------------------------
*       External Subroutines ...
*       ------------------------
        EXTERNAL            ZGEMM , ZGERS , ZTRSM
*
************************************************************************
*
        IF  ( NPANEL .LE. 0 )  RETURN
*
        N = XPANEL(NPANEL+1) - 1
        DO  J = 1, N
          DO P = 1, NRHS
            T(J,P) = R(PERM(J),P)
            S(J,P) = ZERO
          END DO 
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
            CALL  ZTRSM ('LEFT', 'LOWER', 'NO TRANSPOSE', 'UNIT',
     &              NJ, NRHS, ONE, LNZ(JLPNT), JLEN, T(FJ,1), LDR)
*
            CALL  ZGEMM ('NO TRANSPOSE', 'NO TRANSPOSE', JLEN-NJ,
     &             NRHS, NJ, -ONE, LNZ(JLPNT+NJ), JLEN,
     &             T(FJ,1), LDR, ZERO, S, LDR)
*
            JJ = JXPNT + NJ - 1
            DO  J = 1, JLEN-NJ
                JJ         = JJ + 1
                ISUB       = LINDX(JJ)
                DO P = 1, NRHS
                  T(ISUB,P) = T(ISUB,P) + S(J,P)
                  S(J,P)    = ZERO
                END DO
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
            CALL  ZGERS ( NJ, NRHS, LNZ(JLPNT), JLEN, LBDEF,
     &                    IPROW, IPCOL, T(FJ,1), LDR, INFO )
        END IF
*
        DO  JJ = 1, NDEF-LBDEF
            J       = DEF(JJ)
            DO P = 1, NRHS
              T(J,P) = ZERO
            END DO
        END DO
        DO  JPANEL = 1, LSTPNL
            FJ    = XPANEL(JPANEL)
            LJ    = XPANEL(JPANEL+1) - 1
            JLEN  = XLNZ(FJ+1) - XLNZ(FJ)
            JDPTR = XLNZ(FJ)
            DO  J = FJ, LJ
              DO P = 1, NRHS
                T(J,P) = T(J,P)/LNZ(JDPTR)
              END DO
              JDPTR   = JDPTR + JLEN + 1
            END DO
        END DO
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
            DO  J = 1, JLEN
                JJ      = JJ + 1
                ISUB    = LINDX(JJ)
                DO P = 1, NRHS
                  S(J,P) = T(ISUB,P)
                END DO
            END DO
*
            IF  ( JLEN .GT. NJ )  THEN
*
            CALL  ZGEMM ('TRANSPOSE', 'NO TRANSPOSE', NJ,
     &             NRHS, JLEN-NJ, -ONE, LNZ(JLPNT+NJ), JLEN,
     &             S, LDR, ONE, T(FJ,1), LDR)
*
            END IF
*
            CALL  ZTRSM ('LEFT', 'LOWER', 'TRANSPOSE', 'UNIT',
     &              NJ, NRHS, ONE, LNZ(JLPNT), JLEN, T(FJ,1), LDR)
*
        END DO
*
        DO  J = 1, N
          DO P = 1, NRHS
            S(J,P) = T(INVP(J),P)
          END DO
        END DO
*
        RETURN
*
*       --------------
*       End of ZBLKSLV
*       --------------
      END
