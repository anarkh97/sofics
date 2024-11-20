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
*****   ORDNAT ... natural ordering                                *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine returns the identity permutation in vectors
*     PERM and INVP.
*
*   ----------
*   Arguments:
*   ----------
*
*     N         (input) integer
*               Number of equations.
*
*     PERM      (output) integer array, dimension N
*               The permutation vector contains the "natural" ordering;
*               i.e., the initial ordering.
*
*     INVP      (output) integer array, dimension N
*               The inverse of the permutation vector.
*   
************************************************************************
*
      SUBROUTINE  ORDNAT  ( N     , PERM  , INVP                    )
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
        INTEGER             INVP(*), PERM(*)
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
        DO  I = 1, N
            PERM(I) = I
            INVP(I) = I
        END DO
        RETURN
*
*       --------------
*       End of ORDNAT.
*       --------------
      END
