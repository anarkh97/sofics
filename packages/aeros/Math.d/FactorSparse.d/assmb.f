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
*****   ASSMB ... indexed assembly operation                       *****
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
*     Y         (input/output) double precision array, dimension M*Q
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
*     LNZ       (input/output) double precision array, dimension *
*               It contains columns modified by the update matrix.
*
*     LDLNZ     (input) integer
*               Leading dimension of LNZ.
*
************************************************************************
*
      SUBROUTINE  ASSMB   ( M     , Q     , Y     , RELCOL, RELIND,
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
        DOUBLE PRECISION    LNZ(*)
        INTEGER             RELCOL(*), RELIND(*)
        DOUBLE PRECISION    Y(*)
*
************************************************************************
*
*       --------------
*       Parameters ...
*       --------------
        DOUBLE PRECISION    ZERO
        PARAMETER       (   ZERO = 0.0D0 )
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
*       End of ASSMB.
*       -------------
      END
