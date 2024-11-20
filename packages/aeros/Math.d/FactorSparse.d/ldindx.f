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
*****   LDINDX ... load index vector                               *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine computes the second index vector used to
*     implement the doubly-indirect SAXPY-like loops that allow
*     us to accumulate update columns directly into factor storage.
*
*   ----------
*   Arguments:
*   ----------
*
*     JLEN      (input) integer
*               Length of the first column of the supernode,
*               including the diagonal entry.
*
*     LINDX     (input) integer array, dimension JLEN
*               The row indices of the nonzero entries of the first
*               column of the supernode.
*
*     INDMAP    (output) integer array, dimension N
*               This index vector maps every global row index of
*               nonzero entries in the first column of the supernode
*               to its position in the index list relative to the
*               last index in the list.  More precisely, it gives
*               the distance of each index from the last index in the
*               list.
*
************************************************************************
*
      SUBROUTINE  LDINDX  ( JLEN  , LINDX , INDMAP                  )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             JLEN
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             LINDX(*), INDMAP(*)
*
************************************************************************
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             CURLEN, J     , JSUB
*
************************************************************************
*
        CURLEN = JLEN
        DO  J = 1, JLEN
            JSUB         = LINDX(J)
            CURLEN       = CURLEN - 1
            INDMAP(JSUB) = CURLEN
        END DO
        RETURN
*
*       --------------
*       End of LDINDX.
*       --------------
      END
