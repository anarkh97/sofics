************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    Esmond G. Ng and Barry W. Peyton
*   Created:    01-18-1992 (Barry W. Peyton)
*   Modified:   11-22-1994 (Barry W. Peyton)
*   Modified:   12-27-1994 (Barry W. Peyton)
*   Modified:   09-23-1998 (Esmond G. Ng)
*
*   Computational Mathematics and Statistics Section
*   Oak Ridge National Laboratory
*
************************************************************************
************************************************************************
*****   FSUP2  ... find supernodes #2                              *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine is the second of two subroutines for finding
*     a maximal supernode partition.  It's sole purpose is to
*     construct the needed vector of length nsuper: XSUPER.  The
*     first subroutine FSUP1 computes the number of supernodes and
*     the supernode membership vector SNODE, which is of length N.
*
*     ------------
*     Assumptions:
*     ------------
*       This subroutine assumes a postordering of the elimination
*       tree.  It also assumes that the output from FSUP1 is
*       available.
*
*   ----------
*   Arguments:
*   ----------
*
*     N         (input) integer
*               Number of equations.
*
*     NSUPER    (input) integer
*               Number of supernodes.
*
*     SNODE     (input) integer array, dimension N
*               Membership of columns in the supernode partition.
*
*     XSUPER    (output) integer array, dimension NSUPER+1
*               The supernodal partition.
*
************************************************************************
*
      SUBROUTINE  FSUP2   ( N     , NSUPER, SNODE , XSUPER          )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             N     , NSUPER
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             SNODE(*), XSUPER(*)
*
************************************************************************
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             KCOL  , KSUP  , LSTSUP
*
************************************************************************
*
*       ----------------------------------------------
*       Compute the supernode partition vector XSUPER.
*       ----------------------------------------------
        LSTSUP = NSUPER + 1
        DO  KCOL = N, 1, -1
            KSUP = SNODE(KCOL)
            IF  ( KSUP .NE. LSTSUP )  THEN
                XSUPER(LSTSUP) = KCOL + 1
            END IF
            LSTSUP = KSUP
        END DO
        XSUPER(1) = 1
        RETURN
*
*       -------------
*       End of FSUP2.
*       -------------
      END
