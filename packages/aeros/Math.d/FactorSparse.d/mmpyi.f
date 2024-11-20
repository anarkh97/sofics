************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    Esmond G. Ng and Barry W. Peyton
*   Created:    04-26-1996 (Barry W. Peyton)
*   Modified:   09-23-1998 (Esmond G. Ng)
*
*   Computational Mathematics and Statistics Section
*   Oak Ridge National Laboratory
*
************************************************************************
************************************************************************
*****   MMPYI ... matrix-matrix multiply                           *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine performs a matrix-matrix multiply, Z = Z + XY,
*     assuming data structures used in some of our sparse Cholesky
*     factorization codes.
*
*     Note:
*       Matrix X has only 1 column and matrix Y has only 1 row.
*
*   ----------
*   Arguments:
*   ----------
*
*     M         (input) integer
*               Number of rows in X and in Z.
*
*     Q         (input) integer
*               Number of columns in Y and Z.
*
*     ZINDXR    (input) integer array, dimension M
*               It gives the list of rows in Z.
*
*     ZINDXC    (input) integer array, dimension Q
*               It gives the list of columns in Z.
*
*     X         (input) double precision array, dimension M
*               It contains the rows of X.
*
*     Y         (input) double precision array, dimension Q
*               It contains the columns of Y.
*
*     IZ        (input) integer array, dimension Q+1
*               IZ(COL) POINTS TO THE BEGINNING OF COLUMN COL in Z.
*
*     Z         (input/output) double precision array, dimension M*Q
*               It is replaced by Z = Z + XY.
*
*     RELIND    (input) integer array, dimension N
*               Relative indices.
*
************************************************************************
*
      SUBROUTINE  MMPYI   ( M     , Q     , ZINDXR, ZINDXC, X     , 
     &                      Y     , IZ    , Z     , RELIND          )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             M     , Q
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             IZ(*)
        INTEGER             RELIND(*)
        INTEGER             ZINDXR(*), ZINDXC(*)
        DOUBLE PRECISION    X(*), Y(*), Z(*)
*
************************************************************************
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             COL   , I     , ISUB  , K     , ZLAST
        DOUBLE PRECISION    T
*
************************************************************************
*
        DO  K = 1, Q
            COL   = ZINDXC(K)
            ZLAST = IZ(COL+1) - 1
            T     = - Y(K)
            DO  I = 1, M
                ISUB    = ZINDXR(I)
                ISUB    = ZLAST - RELIND(ISUB)
                Z(ISUB) = Z(ISUB) + T*X(I)
            END DO
        END DO
        RETURN
*
*       -------------
*       End of MMPYI.
*       -------------
      END
