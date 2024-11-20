************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    Esmond G. Ng and Barry W. Peyton
*   Created:    05-02-1995 (Barry W. Peyton)
*   Modified:   09-23-1998 (Esmond G. Ng)
*
*   Computational Mathematics and Statistics Section
*   Oak Ridge National Laboratory
*
************************************************************************
************************************************************************
*****   ZASSMB ... indexed assembly operation (Complex data)       *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine performs an indexed assembly (i.e., scatter-add)
*     operation, assuming data structures used in some of our sparse
*     matrix factorization codes.
*
*   ----------
*   Arguments:
*   ----------
*
*     M         (input) integer
*               Number of rows in Y.
*
*     Q         (input) integer
*               Number of columns in Y.
*
*     Y         (input/output) complex*16 array, dimension M*Q
*               Block update to be incorporated into factor storage.
*
*     RELCOL    (input) integer array, dimension Q
*               List of relative column indices.
*
*     RELIND    (input) integer array, dimension M
*               List of relative indices for mapping the updates
*               onto the target columns.
*
*     XLNZ      (input) integer array, dimension *
*               Pointers to the start of each column in the target
*               matrix.
*
*     LNZ       (input/output) complex*16 array, dimension *
*               It contains columns modified by the update matrix.
*
*     LDLNZ     (input) integer
*               Leading dimension of LNZ.
*
************************************************************************
*
      SUBROUTINE  ZASSMB   ( M     , Q     , Y     , RELCOL, RELIND,
     &                      XLNZ  , LNZ   , LDLNZ                   )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             LDLNZ , M     , Q
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             XLNZ(*)
        COMPLEX*16          LNZ(*)
        INTEGER             RELCOL(*), RELIND(*)
        COMPLEX*16          Y(*)
*
************************************************************************
*
*       --------------
*       Parameters ...
*       --------------
        COMPLEX*16         ZERO
        PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             ICOL  , IL    , IR    , IY    , LBOT  ,
     &                      YCOL
*
************************************************************************
*
        IY = 0
        DO  ICOL = 1, Q
            YCOL = LDLNZ - RELCOL(ICOL)
            LBOT = XLNZ(YCOL+1) - 1
            DO  IR = 1, M
                IL      = LBOT - RELIND(IR)
                IY      = IY + 1
                LNZ(IL) = LNZ(IL) + Y(IY)
                Y(IY)   = ZERO
            END DO
        END DO
*
        RETURN
*
*       -------------
*       End of ZASSMB
*       -------------
      END
