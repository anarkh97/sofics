************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    Joseph W.H. Liu
*   Created:    07-17-1985 (Joseph W.H. Liu)
*   Modified:   12-27-1994 (Barry W. Peyton)
*   Modified:   09-23-1998 (Esmond G. Ng)
*
*   Computational Mathematics and Statistics Section
*   Oak Ridge National Laboratory
*
*   Note:
*     This subroutine is kindly provided by Joseph W.H. Liu.
*
************************************************************************
************************************************************************
*****   INVINV ... concatenation of two invp                       *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine performs the mapping of
*
*       ORIGINAL-INVP --> INTERMEDIATE-INVP --> NEW INVP
*
*     and the resulting ordering replaces INVP.  The new permutation
*     vector PERM is also computed.
*
*   ----------
*   Arguments:
*   ----------
*
*     N         (input) integer
*               Number of equations.
*
*     INVP      (input/output) integer array, dimension N
*               The inverse of the first permutation vector.  On
*               return, it contains the concatenation of INVP and
*               INVP2.
*
*     INVP2     (input) integer array, dimension N
*               The inverse of the second permutation vector.
*
*     PERM      (output) integer array, dimension N
*               The permutation vector.  (It can be the same as INVP2.)
*
************************************************************************
*
      SUBROUTINE  INVINV  ( N     , INVP  , INVP2 , PERM            )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             N
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             PERM(*), INVP(*), INVP2(*)
*
************************************************************************
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             I     , INTERM, NODE
*
************************************************************************
*
        DO  I = 1, N
            INTERM  = INVP(I)
            INVP(I) = INVP2(INTERM)
        END DO
*
        DO  I = 1, N
            NODE       = INVP(I)
            PERM(NODE) = I
        END DO
        RETURN
*
*       --------------
*       End of INVINV.
*       --------------
      END
