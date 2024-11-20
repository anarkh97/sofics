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
*****  ZBLKNS ... block null space (complex matrices)              *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*     This subroutine computes the null space of a sparse symmetric
*     semi-definite matrix.  It is assumed that a LDL' factorization
*     of the matrix has been computed.
*
*   ----------
*   Arguments:
*   ----------
*
*     NPANEL    (input) integer
*               Number of panels.
*
*     XPANEL    (input) integer array, dimension NPANEL+1
*               The panel partition.
*
*     XLINDX    (input) integer array, dimension NPANEL+1
*               Pointers to column structure of the Cholesky factor.
*
*     LINDX     (input) integer array, dimension XPANEL(NPANEL+1)-1
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
*               The size of the last panel, which contains the known
*               deficiency.
*
*     NDEF      (input) integer
*               Rank deficiency of the matrix.
*
*     LBDEF     (output) integer
*               If DEFBLK is nonzero, LBDEF is the deficiency of the
*               last block.
*
*     DEF       (input) integer array, dimension NDEF
*               It identifies the columns of the matrix that are
*               linearly dependent.
*
*     IPCOL     (output) integer array, dimension DEFBLK
*               If DEFBLK is nonzero, IPCOL contains the column
*               pivoting sequence for the last block.
*
*     INVP      (input) integer array, dimension N
*               The inverse of the permutation vector.
*
*     NS        (input) COMPLEX*16 array, dimension (N,NDEF)
*               The null space of the matrix.
*
*     LDNS      (input) integer
*               Leading dimension of NS.
*
*     TEMP      (temporary) COMPLEX 16 array, dimension N
*               Temporary work space.
*
************************************************************************
*
      SUBROUTINE  ZBLKNS   ( NPANEL, XPANEL, XLINDX, LINDX , XLNZ  ,
     &                      LNZ   , DEFBLK, NDEF  , LBDEF , DEF   ,
     &                      IPCOL , INVP  , NS    , LDNS  , TEMP    )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             DEFBLK, LDNS  , NDEF  , NPANEL, LBDEF
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             XPANEL(*)
        INTEGER             XLINDX(*), LINDX(*)
        INTEGER             XLNZ(*)
        COMPLEX*16          LNZ(*)
        INTEGER             INVP(*)
        INTEGER             DEF(*)
        INTEGER             IPCOL(*)
        COMPLEX*16          NS(LDNS,*)
        COMPLEX*16          TEMP(*)
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
        INTEGER             FJ    , I     , ISUB  , J     , JCOL  ,
     &                      JJ    , JLEN  , JLPNT , JPANEL, JXPNT ,
     &                      LSTPNL, N     , NJ
*
*       ------------------------
*       External Subroutines ...
*       ------------------------
        EXTERNAL            ZGEMV , ZLASWP, ZTRSM
*
************************************************************************
*
        IF  ( NPANEL .LE. 0 )  RETURN
*
        IF  ( NDEF .EQ. 0 )  RETURN
*
*       ---------------
*       Initialization.
*       ---------------
        N = XPANEL(NPANEL+1) - 1
        DO  J = 1, NDEF
            DO  I = 1, N
                NS(I,J) = ZERO
            END DO
            I = DEF(J)
            NS(I,J) = ONE
        END DO
*
*       -------------------------
*       Backward substitution ...
*       -------------------------
        LSTPNL = NPANEL
        IF  ( DEFBLK .NE. 0 )  THEN
            LSTPNL = NPANEL - 1
            FJ     = XPANEL(NPANEL)
            NJ     = N - FJ + 1
            JLEN   = XLNZ(FJ+1) - XLNZ(FJ)
            JLPNT  = XLNZ(FJ)
            CALL  ZTRSM ( 'LEFT', 'UPPER', 'NO-TRANSPOSE', 'NON-UNIT',
     &                    NJ, NDEF, ONE, LNZ(JLPNT), JLEN,
     &                    NS(FJ,1), LDNS )
            CALL  ZLASWP ( NDEF, NS(FJ,1), LDNS, 1, NJ, IPCOL, -1 )
        END IF
*
        DO  JCOL = 1, NDEF
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
                    TEMP(J) = NS(ISUB,JCOL)
                END DO
*
                IF  ( JLEN .GT. NJ )  THEN
                    CALL  ZGEMV (
     &                      'TRANSPOSE',
     &                      JLEN-NJ, NJ,
     &                      -ONE,
     &                      LNZ(JLPNT+NJ), JLEN,
     &                      TEMP, 1,
     &                      ONE,
     &                      NS(FJ,JCOL), 1
     &                  )
                END IF
*
                CALL  ZTRSM (
     &                  'LEFT', 'LOWER', 'TRANSPOSE', 'UNIT',
     &                  NJ, 1,
     &                  ONE,
     &                  LNZ(JLPNT), JLEN,
     &                  NS(FJ,JCOL), NJ
     &          )
*
            END DO
        END DO
*
*       ---------------------
*       Apply the permutation
*       ---------------------
        DO  J = 1, NDEF
            DO  I = 1, N
                TEMP(I) = NS(I,J)
            END DO
            DO  I = 1, N
                NS(I,J) = TEMP(INVP(I))
            END DO
        END DO
*
        RETURN
*
*       -------------
*       End of ZBLKNS
*       -------------
      END
