************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    Esmond G. Ng and Barry W. Peyton
*   Created:    03-06-1995 (Barry W. Peyton)
*   Modified:   09-23-1998 (Esmond G. Ng)
*
*   Computational Mathematics and Statistics Section
*   Oak Ridge National Laboratory
*
************************************************************************
************************************************************************
*****   BLKLDL ... block general sparse Cholesky factorization     *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine calls the block general sparse Cholesky
*     factorization subroutine, BLKL2, to compute a LDL'
*     factorization.
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
*     LNZ       (input/output) double precision array, dimension
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
*     TMPSIZ    (input) integer
*               Size of work array TEMP.
*
*     TEMP      (temporary) double precision array, dimension TMPSIZ
*               Work space for BLKLDL.
*
*     IWSIZE    (input) integer
*               Size of work array IWORK.
*
*     IWORK     (temporary) integer array, dimension 2*(N+NPANEL)
*               Hold various work arrays.
*
*     RWSIZE    (input) integer
*               Size of work array RWORK.
*
*     RWORK     (temporary) double precision array, dimension RWSIZE
*               Temporary work array; its size should be as big as
*               the largest panel in the matrix.
*
*     IFLAG     (output) integer
*               Error flag.
*
*               IFLAG =  0: successful factorization.
*               IFLAG = 31: insufficient work space in TEMP.
*               IFLAG = 32: insufficient work space in IWORK.
*               IFLAG = 33: insufficient work space in RWORK.
*
************************************************************************
*
      SUBROUTINE  BLKLDL  ( NPANEL, XPANEL, PNODE , XLINDX, LINDX , 
     &                      XLNZ  , LNZ   , DEFBLK, ASDEF , NDEF  ,
     &                      LBDEF ,
     &                      DEF   , TOL   , IPROW , IPCOL , TMPSIZ,
     &                      TEMP  , IWSIZE, IWORK , RWSIZE, RWORK ,
     &                      IFLAG                                   )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             DEFBLK, IFLAG , IWSIZE, LBDEF , NDEF  ,
     &                      NPANEL, RWSIZE, TMPSIZ, ASDEF
        DOUBLE PRECISION    TOL
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             XPANEL(*)
        INTEGER             PNODE(*)
        INTEGER             XLINDX(*), LINDX(*)
        INTEGER             XLNZ(*)
        DOUBLE PRECISION    LNZ(*)
        DOUBLE PRECISION    TEMP(*)
        INTEGER             IPROW(*), IPCOL(*)
        INTEGER             DEF(*)
        INTEGER             IWORK(*)
        DOUBLE PRECISION    RWORK(*)
*
**********************************************************************
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             N
*
*       ------------------------
*       External Subroutines ...
*       ------------------------
        EXTERNAL            BLKL2
*
**********************************************************************
*
        N = XPANEL(NPANEL+1) - 1
        IFLAG = 0
        IF  ( IWSIZE .LT. 2*(N+NPANEL) )  THEN
            IFLAG = 32
            RETURN
        END IF
        CALL  BLKL2 ( NPANEL, XPANEL, PNODE, XLINDX, LINDX, XLNZ,
     &                LNZ, DEFBLK, ASDEF, NDEF, LBDEF, DEF, TOL,
     &                IPROW,
     &                IPCOL,
     &                IWORK(1),
     &                IWORK(NPANEL+1),
     &                IWORK(2*NPANEL+1),
     &                IWORK(2*NPANEL+N+1),
     &                TMPSIZ, TEMP, RWSIZE, RWORK, IFLAG )
        RETURN
*
*       --------------
*       End of BLKLDL.
*       --------------
      END
