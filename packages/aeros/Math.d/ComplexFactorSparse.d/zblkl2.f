************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    Esmond G. Ng and Barry W. Peyton
*   Created:    05-07-1998 (Esmond G. Ng)
*   Modified:   09-23-1998 (Esmond G. Ng)
*
*   Computational Mathematics and Statistics Section
*   Oak Ridge National Laboratory
*
************************************************************************
************************************************************************
** ZBLKL2 ... block general sparse Cholesky factorization (complex) ****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine computes the Cholesky factorization of a
*     sparse symmetric matrix.  It returns the LDL' factorization.
*     The matrix can be semi-definite.  The computation is
*     organized around kernels that perform panel-to-panel updates.
*     Level-3 BLAS subroutines are used whenever possible.  Each
*     panel is a dense lower trapezoidal submatrix within a panel
*     (portion of a supernode) of the Cholesky factor, but is stored
*     in a rectangular array.
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
*     PNODE     (input) integer array, dimension N
*               Membership of columns in the panel partition.
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
*     LNZ       (input/output) complex*16 array, dimension
*                   XLNZ(N+1)-1
*               Numerical values of nonzero entries in the Cholesky
*               factor, stored by columns.
*
*     DEFBLK    (input) integer
*               The size of the last supernode, which contains the
*               known deficiency.
*
*     NDEF      (output) integer
*               Rank deficiency of the matrix.
*
*     ASDEF     (input) the assumed deficiency of the matrix A
*
*     LBDEF     (output) integer
*               If DEFBLK is nonzero, LBDEF is the deficiency of the
*               last block.
*
*     DEF       (output) integer array, dimension NDEF
*               It identifies the columns of the matrix that are
*               linearly dependent.
*
*     TOL       (input) double precision
*               Tolerance for determining rank deficiency.
*
*     IPROW     (output) integer array, dimension DEFBLK
*               If DEFBLK is nonzero, IPROW contains the row pivoting
*               sequence for the last block.
*
*     IPCOL     (output) integer array, dimension DEFBLK
*               If DEFBLK is nonzero, IPCOL contains the column
*               pivoting sequence for the last block.
*
*     LINK      (temporary) integer array, dimension NPANEL
*               It links together the panels in a panel row.
*
*     LENGTH    (temporary) integer array, dimension NPANEL
*               It contains length of the active portion of each panel.
*
*     INDMAP    (temporary) integer array, dimension N
*               It is a vector into which the global indices are
*               scattered.
*
*     RELIND    (temporary) integer array, dimension N
*               It maps locations in the updating columns to the
*               corresponding locations in the updated columns.
*               (RELIND is gathered from INDMAP).
*
*     TMPSIZ    (input) integer
*               Size of work array TEMP.
*
*     TEMP      (temporary) complex*16 array, dimension TMPSIZ
*               Work space for ZBLKLDL for accumulating updates.  It
*               must accommodate all columns of a panel.
*
*     RWSIZE    (input) integer
*               Size of work array RWORK.
*
*     RWORK     (temporary) complex*16 array, dimension RWSIZE
*               Temporary work array; its size should be as big as
*               the largest panel in the matrix.
*
*     IFLAG     (output) integer
*               Error flag.
*
*               IFLAG =  0: successful factorization.
*               IFLAG = 31: insufficient work space in TEMP.
*               IFLAG = 33: insufficient work space in RWORK.
*
*   To Do:
*     Need to monitor work space requirement in RWORK?
*
************************************************************************
*
      SUBROUTINE  ZBLKL2   ( NPANEL, XPANEL, PNODE , XLINDX, LINDX ,
     &                      XLNZ  , LNZ   , DEFBLK, ASDEF,  NDEF  ,
     &                      LBDEF ,
     &                      DEF   , TOL   , IPROW , IPCOL , LINK  ,
     &                      LENGTH, INDMAP, RELIND, TMPSZE, TEMP  ,
     &                      RWSIZE, RWORK , IFLAG                   )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             DEFBLK, IFLAG , LBDEF , NDEF  , NPANEL,
     &                      RWSIZE, TMPSZE, ASDEF
        DOUBLE PRECISION    TOL
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             XPANEL(*)
        INTEGER             PNODE(*)
        INTEGER             XLINDX(*), LINDX(*)
        INTEGER             XLNZ(*)
        COMPLEX*16          LNZ(*)
        COMPLEX*16          TEMP(*)
        INTEGER             INDMAP(*), LENGTH(*), LINK(*), RELIND(*)
        INTEGER             IPROW(*), IPCOL(*)
        INTEGER             DEF(*)
        COMPLEX*16          RWORK(*)
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
        INTEGER             COL   , FJ    , FK    , I     , IDPTR ,
     &                      II    , ILPNT , ILPTR , INDDIF, INFO  ,
     &                      IWPTR , J     , JDEF  , JDPTR , JLEN  ,
     &                      JLPNT , JLPTR , JPANEL, JXPNT , K     ,
     &                      KFIRST, KLAST , KLEN  , KLPNT , KPANEL,
     &                      KPLEN , KXPNT , LJ    , LK    , N     ,
     &                      NJ    , NK    , NUPS  , NXKPAN, NXT   ,
     &                      NXTPAN, STORE
*
*       ------------------------
*       External Subroutines ...
*       ------------------------
        EXTERNAL            ZASSMB, ZCOPY, ZGECP, ZGEMM, ZPOTF2_LDL,
     &                      ZSCAL , ZTRSM , IGATHR, LDINDX, ZMMPYI
*
************************************************************************
*
        IFLAG = 0
        N     = XPANEL(NPANEL+1) - 1
*
*       -----------------------------------
*       Initialize empty row lists in LINK.
*       Zero out TEMP and RWORK.
*       -----------------------------------
        DO  JPANEL = 1, NPANEL
            LINK(JPANEL) = 0
        END DO
        DO  I = 1, TMPSZE
            TEMP(I) = ZERO
        END DO
        DO  I = 1, RWSIZE
            RWORK(I) = ZERO
        END DO
*
*       -------------------------
*       For each panel JPANEL ...
*       -------------------------
        DO  JPANEL = 1, NPANEL
*
*           ----------------------------------------------
*           FJ     :  first row/column of JPANEL.
*           LJ     :  last row/column of JPANEL.
*           NJ     :  number of rows/columns in JPANEL.
*           JLEN   :  length of JPANEL.
*           JXPNT  :  pointer to index of first nonzero
*                       in JPANEL.
*           JLPNT  :  pointer to location of first nonzero
*                       in column FJ.
*           ----------------------------------------------
            FJ    = XPANEL(JPANEL)
            LJ    = XPANEL(JPANEL+1) - 1
            NJ    = LJ - FJ + 1
            JLEN  = XLNZ(FJ+1) - XLNZ(FJ)
            JXPNT = XLINDX(JPANEL)
            JLPNT = XLNZ(FJ)
*
*           ----------------------------------------------------
*           Set up INDMAP to map the entries in update columns
*           to their corresponding positions in updated columns,
*           relative to the bottom of each updated column.
*           ----------------------------------------------------
            CALL  LDINDX ( JLEN, LINDX(JXPNT), INDMAP )
*
*           -----------------------------------------
*           For every panel KPANEL in row(JPANEL) ...
*           -----------------------------------------
            KPANEL = LINK(JPANEL)
  100       CONTINUE
            IF  ( KPANEL .GT. 0 )  THEN
                NXKPAN = LINK(KPANEL)
*
*               ----------------------------------------------
*               Get info about the cmod(JPANEL,KPANEL) update.
*
*               FK     :  first row/column of KPANEL.
*               NK     :  number of rows/columns in KPANEL.
*               KPLEN  :  length of KPANEL.
*               KLEN   :  length of active portion of KPANEL.
*               KXPNT  :  pointer to index of first nonzero
*                           in active portion of KPANEL.
*               KLPNT  :  pointer to location of first nonzero
*                           in active portion of column FK.
*               KUPNT  :  pointer to location of first nonzero
*                           in active portion of row FK.
*               ----------------------------------------------
                FK    = XPANEL(KPANEL)
                LK    = XPANEL(KPANEL+1) - 1
                NK    = LK - FK + 1
                KPLEN = XLNZ(FK+1) - XLNZ(FK)
                KLEN  = LENGTH(KPANEL)
                KXPNT = XLINDX(KPANEL+1) - KLEN
                KLPNT = XLNZ(FK+1) - KLEN
*
*               -----------------------------------------------
*               Perform cmod(JPANEL,KPANEL), with special cases
*               handled differently.
*
*               NUPS   :  number of rows/columns updated.
*               -----------------------------------------------
*
                IF  ( KLEN .EQ. JLEN )  THEN
*
*                   -------------------------------------------
*                   Dense cmod(JPANEL,KPANEL).
*                   JPANEL and KPANEL have identical structure.
*                   -------------------------------------------
                    IDPTR = XLNZ(FK)
                    ILPTR = KLPNT
                    IWPTR = 1
                    DO  K = 1, NK
                        CALL  ZCOPY ( NJ, LNZ(ILPTR), 1,
     &                                RWORK(IWPTR), 1 )
                        CALL  ZSCAL ( NJ, LNZ(IDPTR), RWORK(IWPTR), 1 )
                        IDPTR = IDPTR + KPLEN + 1
                        ILPTR = ILPTR + KPLEN
                        IWPTR = IWPTR + NJ
                    END DO
                    CALL  ZGEMM (
     &                      'NO TRANSPOSE', 'TRANSPOSE',
     &                      JLEN, NJ, NK,
     &                      -ONE,
     &                      LNZ(KLPNT), KPLEN,
     &                      RWORK, NJ,
     &                      ONE,
     &                      LNZ(JLPNT), JLEN
     &                  )
                    NUPS = NJ
                    IF  ( KLEN .GT. NJ )  THEN
                        NXT = LINDX(JXPNT+NJ)
                    END IF
*
                ELSE
*
*                   ---------------------------
*                   Sparse cmod(JPANEL,KPANEL).
*                   ---------------------------
*
*                   ------------------------------------
*                   Determine the number of rows/columns
*                   to be updated.
*                   ------------------------------------
                    DO  I = 0, KLEN-1
                        NXT = LINDX(KXPNT+I)
                        IF  ( NXT .GT. LJ )  GO TO 200
                    END DO
                    I = KLEN
  200               CONTINUE
                    NUPS = I
*
                    IF  ( NK .EQ. 1 )  THEN
*
*                       ----------------------------------------
*                       Updating target panel by a trivial panel
*                       (with one column).
*                       ----------------------------------------
                        CALL  ZCOPY ( NUPS, LNZ(KLPNT), 1,
     &                                RWORK, 1 )
                        CALL  ZSCAL ( NUPS, LNZ(XLNZ(FK)),
     &                                RWORK, 1 )
                        CALL  ZMMPYI (
     &                          KLEN, NUPS, LINDX(KXPNT), LINDX(KXPNT),
     &                          LNZ(KLPNT), RWORK,
     &                          XLNZ, LNZ, INDMAP
     &                      )
*
                    ELSE
*
*                       ------------------------------------------
*                       KFIRST :  first index of active portion of
*                                   panel KPANEL (first column to
*                                   be updated).
*                       KLAST  :  last index of active portion of
*                                   panel KPANEL.
*                       ------------------------------------------
*
                        KFIRST = LINDX(KXPNT)
                        KLAST  = LINDX(KXPNT+KLEN-1)
                        INDDIF = INDMAP(KFIRST) - INDMAP(KLAST)
*
                        IF  ( INDDIF .LT. KLEN )  THEN
*
*                           ---------------------------------------
*                           Dense cmod(JPANEL,KPANEL).
*
*                           ILPNT  :  pointer to first nonzero in
*                                       column KFIRST.
*                           ---------------------------------------
*
                            ILPNT = XLNZ(KFIRST) + ( KFIRST - FJ )
                            IDPTR = XLNZ(FK)
                            ILPTR = KLPNT
                            IWPTR = 1
                            DO  K = 1, NK
                                CALL  ZCOPY ( NUPS, LNZ(ILPTR), 1,
     &                                        RWORK(IWPTR), 1 )
                                CALL  ZSCAL ( NUPS, LNZ(IDPTR),
     &                                        RWORK(IWPTR), 1 )
                                IDPTR = IDPTR + KPLEN + 1
                                ILPTR = ILPTR + KPLEN
                                IWPTR = IWPTR + NUPS
                            END DO
                            CALL  ZGEMM (
     &                              'NO TRANSPOSE', 'TRANSPOSE',
     &                              KLEN, NUPS, NK,
     &                              -ONE,
     &                              LNZ(KLPNT), KPLEN,
     &                              RWORK, NUPS,
     &                              ONE,
     &                              LNZ(ILPNT), JLEN
     &                          )
*
                        ELSE
*
*                           -----------------------------------
*                           General sparse cmod(JPANEL,KPANEL).
*                           Compute cmod(JPANEL,KPANEL) update
*                           in work storage.
*                           -----------------------------------
                            STORE = KLEN * NUPS
                            IF  ( STORE .GT. TMPSZE )  THEN
                                IFLAG = 31
                                RETURN
                            END IF
*
*                           ---------------------------------
*                           Gather indices of KPANEL relative
*                           to JPANEL.
*                           ---------------------------------
                            CALL  IGATHR ( KLEN, LINDX(KXPNT),
     &                                     INDMAP, RELIND )
*
                            IDPTR = XLNZ(FK)
                            ILPTR = KLPNT
                            IWPTR = 1
                            DO  K = 1, NK
                                CALL  ZCOPY ( NUPS, LNZ(ILPTR), 1,
     &                                        RWORK(IWPTR), 1 )
                                CALL  ZSCAL ( NUPS, LNZ(IDPTR),
     &                                        RWORK(IWPTR), 1 )
                                IDPTR = IDPTR + KPLEN + 1
                                ILPTR = ILPTR + KPLEN
                                IWPTR = IWPTR + NUPS
                            END DO
                            CALL  ZGEMM (
     &                              'NO TRANSPOSE', 'TRANSPOSE',
     &                              KLEN, NUPS, NK,
     &                              -ONE,
     &                              LNZ(KLPNT), KPLEN,
     &                              RWORK, NUPS,
     &                              ONE,
     &                              TEMP, KLEN
     &                          )
*
*                           -----------------------------------------
*                           Incorporate the cmod(JPANEL,KPANEL) block
*                           update into the appropriate columns of L.
*                           -----------------------------------------
                            CALL  ZASSMB (
     &                              KLEN, NUPS, TEMP,
     &                              RELIND(1), RELIND(1),
     &                              XLNZ(FJ), LNZ, JLEN
     &                          )
*
                        END IF
*
                    END IF
*
                END IF
*
*               ----------------------------------------------
*               Link KPANEL into linked list of the next panel
*               it will update and decrement KPANEL'S active
*               length.
*               ----------------------------------------------
                IF  ( KLEN .GT. NUPS )  THEN
                    NXTPAN         = PNODE(NXT)
                    LINK(KPANEL)   = LINK(NXTPAN)
                    LINK(NXTPAN)   = KPANEL
                    LENGTH(KPANEL) = KLEN - NUPS
                ELSE
                    LENGTH(KPANEL) = 0
                END IF
*
*               -----------------------------
*               Next updating panel (KPANEL).
*               -----------------------------
                KPANEL = NXKPAN
                GO TO 100
*
            END IF
*
            IF  ( JPANEL .EQ. NPANEL  .AND.  DEFBLK .NE. 0 )  THEN
*
*               --------------------------------------------------
*               If this is the last block and if the special flag
*               is turned on, factor the last block using complete
*               pivoting.
*               --------------------------------------------------
                CALL  ZGECP ( NJ, LNZ(JLPNT), JLEN, TOL, ASDEF,
     &                        LBDEF,
     &                        IPROW, IPCOL, INFO )
                IF  ( INFO .NE. 0 )  THEN
                    IFLAG = 101
                    RETURN
                END IF
                IF  ( LBDEF .NE. 0 )  THEN
                    COL = N - LBDEF + 1
                    DO  I = NDEF+1, NDEF+LBDEF
                        DEF(I) = COL
                        COL    = COL + 1
                    END DO
                    NDEF = NDEF + LBDEF
                END IF
*
            ELSE
*
*               ---------------------------------------------
*               Otherwise, apply partial Cholesky to the diagonal
*               block.
*               ---------------------------------------------
                CALL  ZPOTF2_LDL ( 'LOWER', NJ, LNZ(JLPNT), JLEN,
     &                          JDEF, DEF(NDEF+1), TOL, RWORK, INFO )
                IF  ( INFO .NE. 0 )  THEN
                    write(6,*) '*** INFO = ',INFO
                    write(6,*) '*** RANK DEFICIENCY = ',NDEF
                    IFLAG = 102
                    RETURN
                END IF
*               ---------------
*               Update columns.
*               ---------------
                CALL  ZTRSM (
     &                  'RIGHT', 'LOWER', 'TRANSPOSE', 'UNIT',
     &                  JLEN-NJ, NJ,
     &                  ONE,
     &                  LNZ(JLPNT), JLEN,
     &                  LNZ(JLPNT+NJ), JLEN
     &              )
                JDPTR = XLNZ(FJ)
                JLPTR = JLPNT + NJ
                DO  J = FJ, LJ
                    CALL  ZSCAL ( JLEN-NJ, ONE/LNZ(JDPTR),
     &                          LNZ(JLPTR), 1 )
                    JDPTR = JDPTR + JLEN + 1
                    JLPTR = JLPTR + JLEN
                END DO
*
                IF  ( JDEF .NE. 0 )  THEN
                    DO  I = NDEF+1, NDEF+JDEF
                        COL    = DEF(I) + (FJ-1)
                        DEF(I) = COL
                        DO  II = XLNZ(COL)+NJ, XLNZ(COL+1)-1
                            LNZ(II) = ZERO
                        END DO
                    END DO
                    NDEF = NDEF + JDEF
                END IF
            END IF
*
*           ---------------------------------------------
*           Insert JPANEL into linked list of first panel
*           it will update.
*           ---------------------------------------------
            IF  ( JLEN .GT. NJ )  THEN
                NXT            = LINDX(JXPNT+NJ)
                NXTPAN         = PNODE(NXT)
                LINK(JPANEL)   = LINK(NXTPAN)
                LINK(NXTPAN)   = JPANEL
                LENGTH(JPANEL) = JLEN - NJ
            ELSE
                LENGTH(JPANEL) = 0
            END IF
*
        END DO
*
        RETURN
*
*       -------------
*       End of ZBLKL2
*       -------------
      END
