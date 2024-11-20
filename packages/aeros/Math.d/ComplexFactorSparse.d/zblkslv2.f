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
*****  ZBLKSLV2 ... block triangular solutions (complex)           *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     Given the Cholesky factorization of a sparse symmetric matrix,
*     this subroutine performs the triangular solution for 2 vectors.
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
*     R1 - R2   (input) COMPLEX*16 array, dimension N
*               The right hand side vectors.
*
*     S1 - S2   (output) COMPLEX*16 array, dimension N
*               The solution vectors.
*
*     T1 - T2   (temporary) COMPLEX*16 array, dimension N
*               Temporary work arrays.
*
************************************************************************
*
      SUBROUTINE  ZBLKSLV2 ( NPANEL, XPANEL, XLINDX, LINDX , XLNZ  ,
     &                      LNZ   , DEFBLK, NDEF  , LBDEF , DEF   ,
     &                      IPROW , IPCOL , PERM  , INVP  , R1, R2,
     &                      S1, S2, T1, T2 )
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
        COMPLEX*16          R1(*), R2(*)
        COMPLEX*16          S1(*), S2(*)
        COMPLEX*16          T1(*), T2(*)
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
     &                      LJ    , LSTPNL, N     , NJ
*
*       ------------------------
*       External Subroutines ...
*       ------------------------
        EXTERNAL            ZGEMV , ZGERS , ZTRSM
*
************************************************************************
*
        IF  ( NPANEL .LE. 0 )  RETURN
*
        N = XPANEL(NPANEL+1) - 1
        DO  J = 1, N
            T1(J) = R1(PERM(J))
            S1(J) = ZERO
            T2(J) = R2(PERM(J))
            S2(J) = ZERO
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
     &              NJ, 1, ONE, LNZ(JLPNT), JLEN, T1(FJ), NJ)
            CALL  ZTRSM ('LEFT', 'LOWER', 'NO TRANSPOSE', 'UNIT',
     &              NJ, 1, ONE, LNZ(JLPNT), JLEN, T2(FJ), NJ)
*
            CALL  ZGEMV ('NO TRANSPOSE', JLEN-NJ, NJ, -ONE, 
     &             LNZ(JLPNT+NJ), JLEN, T1(FJ), 1, ZERO, S1, 1)
            CALL  ZGEMV ('NO TRANSPOSE', JLEN-NJ, NJ, -ONE, 
     &             LNZ(JLPNT+NJ), JLEN, T2(FJ), 1, ZERO, S2, 1)
*
            JJ = JXPNT + NJ - 1
            DO  J = 1, JLEN-NJ
                JJ         = JJ + 1
                ISUB       = LINDX(JJ)
                T1(ISUB) = T1(ISUB) + S1(J)
                S1(J)    = ZERO
                T2(ISUB) = T2(ISUB) + S2(J)
                S2(J)    = ZERO
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
            CALL  ZGERS ( NJ, 1, LNZ(JLPNT), JLEN, LBDEF,
     &                    IPROW, IPCOL, T1(FJ), NJ, INFO )
            CALL  ZGERS ( NJ, 1, LNZ(JLPNT), JLEN, LBDEF,
     &                    IPROW, IPCOL, T2(FJ), NJ, INFO )
        END IF
*
        DO  JJ = 1, NDEF-LBDEF
            J       = DEF(JJ)
            T1(J) = ZERO
            T2(J) = ZERO
        END DO
        DO  JPANEL = 1, LSTPNL
            FJ    = XPANEL(JPANEL)
            LJ    = XPANEL(JPANEL+1) - 1
            JLEN  = XLNZ(FJ+1) - XLNZ(FJ)
            JDPTR = XLNZ(FJ)
            DO  J = FJ, LJ
               T1(J) = T1(J)/LNZ(JDPTR)
               T2(J) = T2(J)/LNZ(JDPTR)
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
            DO  J = 1, JLEN
                JJ      = JJ + 1
                ISUB    = LINDX(JJ)
                S1(J) = T1(ISUB)
                S2(J) = T2(ISUB)
            END DO
*
            IF  ( JLEN .GT. NJ )  THEN
                CALL  ZGEMV ( 'TRANSPOSE',JLEN-NJ, NJ,-ONE,
     &                LNZ(JLPNT+NJ), JLEN, S1, 1, ONE,
     &                T1(FJ), 1)
                CALL  ZGEMV ( 'TRANSPOSE',JLEN-NJ, NJ,-ONE,
     &                LNZ(JLPNT+NJ), JLEN, S2, 1, ONE,
     &                T2(FJ), 1)
            END IF
*
            CALL  ZTRSM ('LEFT', 'LOWER', 'TRANSPOSE',
     &              'UNIT', NJ, 1, ONE, LNZ(JLPNT), JLEN,
     &              T1(FJ), NJ)
            CALL  ZTRSM ('LEFT', 'LOWER', 'TRANSPOSE',
     &              'UNIT', NJ, 1, ONE, LNZ(JLPNT), JLEN,
     &              T2(FJ), NJ)
*
        END DO
*
        DO  J = 1, N
            S1(J) = T1(INVP(J))
            S2(J) = T2(INVP(J))
        END DO
*
        RETURN
*
*       --------------
*       End of ZBLKSLV
*       --------------
      END
