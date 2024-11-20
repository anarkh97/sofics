************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    Esmond G. Ng and Barry W. Peyton
*   Created:    01-18-1992 (Barry W. Peyton)
*   Modified:   11-11-1994 (Barry W. Peyton)
*   Modified:   12-27-1994 (Barry W. Peyton)
*   Modified:   09-23-1998 (Esmond G. Ng)
*
*   Computational Mathematics and Statistics Section
*   Oak Ridge National Laboratory
*
************************************************************************
************************************************************************
*****   FSUP1 ... find supernodes #1                               *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine is the first of two subroutines for finding
*     a maximal supernode partition.  It returns only the number of
*     supernodes NSUPER and the supernode membership vector SNODE,
*     which is of length N.  The vector XSUPER of length NSUPER is
*     computed subsequently by the companion subroutine FSUP2.
*
*     -----------------------
*     Method and Assumptions:
*     -----------------------
*       This subroutine uses the elimination tree and the factor
*       column counts to compute the supernode partition; it also
*       assumes a postordering of the elimination tree.
*
*   ----------
*   Arguments:
*   ----------
*
*     MAXSUP    (input) integer
*               The maximum number of columns in each supernode.
*
*     DEFBLK    (input) integer
*               The size of the last supernode, which contains the
*               known deficiency.
*
*     N         (input) integer
*               Number of equations.
*
*     ETPAR     (input) integer array, dimension N
*               The elimination tree of the postordered matrix.
*
*     COLCNT    (input) integer array, dimension N
*               The number of nonzeros in each column of the Cholesky
*               factor, including the diagonal entry.
*
*     NOFSUB    (output) integer
*               Number of row indices in the compressed structure
*               representation.
*
*     NSUPER    (output) integer
*               Number of supernodes.
*
*     SNODE     (output) integer array, dimension N
*               Membership of columns in the supernode partition.
*
************************************************************************
*
      SUBROUTINE  FSUP1   ( MAXSIZ, DEFBLK, N     , ETPAR , COLCNT,
     &                      NOFSUB, NSUPER, SNODE                   )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             DEFBLK, MAXSIZ, N     , NOFSUB, NSUPER
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             ETPAR(*)
        INTEGER             COLCNT(*)
        INTEGER             SNODE(*)
*
************************************************************************
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             KCOL  , SUPSIZ , MAXDBK
*
************************************************************************
*
*       ------------------------------------------------
*       Compute the maximum size for the last supernode.
*       ------------------------------------------------
        IF  ( DEFBLK .NE. 0 )  THEN
            MAXDBK = N-1
            DO  KCOL = N, 2, -1
                IF  ( ETPAR(KCOL-1) .EQ. KCOL )  THEN
                    IF  ( COLCNT(KCOL-1) .EQ. COLCNT(KCOL)+1 )  THEN
                        GO TO 50
                    END IF
                END IF
                MAXDBK = N-KCOL+1
                EXIT
  50            CONTINUE
            END DO
            DEFBLK = MIN(MAXDBK,DEFBLK)
        ENDIF

*       --------------------------------------------
*       Compute the fundamental supernode partition.
*       --------------------------------------------
        NSUPER   = 1
        SNODE(1) = 1
        NOFSUB   = COLCNT(1)
        SUPSIZ   = 1
        DO  KCOL = 2, N-DEFBLK
            IF  ( ETPAR(KCOL-1) .EQ. KCOL )  THEN
                IF  ( COLCNT(KCOL-1) .EQ. COLCNT(KCOL)+1 )  THEN
                    IF  ( SUPSIZ .LT. MAXSIZ ) THEN
                        SUPSIZ      = SUPSIZ + 1
                        SNODE(KCOL) = NSUPER
                        GO TO 100
                    END IF
                END IF
            END IF
            NSUPER      = NSUPER + 1
            SUPSIZ      = 1
            SNODE(KCOL) = NSUPER
            NOFSUB      = NOFSUB + COLCNT(KCOL)
  100       CONTINUE
        END DO
        IF  ( DEFBLK .NE. 0 )  THEN
            NSUPER = NSUPER + 1
            DO  KCOL = N-DEFBLK+1, N
                SNODE(KCOL) = NSUPER
            END DO
            NOFSUB = NOFSUB + COLCNT(N-DEFBLK+1)
        END IF
        RETURN
*
*       -------------
*       End of FSUP1.
*       -------------
      END
