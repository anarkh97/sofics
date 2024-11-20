************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    Esmond G. Ng and Barry W. Peyton
*   Created:    12-27-1994 (Barry W. Peyton)
*   Modified:   09-23-1998 (Esmond G. Ng)
*
*   Computational Mathematics and Statistics Section
*   Oak Ridge National Laboratory
*
************************************************************************
************************************************************************
*****   IGATHR ... integer gather operation                        *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine performs a standard integer gather operation.
*
*   ----------
*   Arguments:
*   ----------
*
*     KLEN      (input) integer
*               Length of the list of global indices.
*
*     LINDX     (input) integer array, dimension KLEN
*               List of global indices.
*
*     INDMAP    (input) integer array, dimension N
*               It contains the required relative indices and is
*               indexed by the global indices.
*
*     RELIND    (output) integer array, dimension KLEN
*               List of relative indices.
*
************************************************************************
*
      SUBROUTINE  IGATHR  ( KLEN  , LINDX, INDMAP, RELIND           )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             KLEN
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             LINDX(*)
        INTEGER             INDMAP(*), RELIND(*)
*
************************************************************************
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             I
*
************************************************************************
*
        DO  I = 1, KLEN  
            RELIND(I) = INDMAP(LINDX(I))
        END DO
        RETURN
*
*       --------------
*       End of IGATHR.
*       --------------
      END
