************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    Joseph W.H. Liu
*   Created:    07-17-1985 (Joseph w.H. Liu)
*   Modified:   12-27-1994 (Barry W. Peyton)
*   Modified:   09-23-1998 (Barry W. Peyton)
*
*   Computational Mathematics and Statistics Section
*   Oak Ridge National Laboratory
*
*   Note:
*     This subroutine is kindly provided by Joseph W.H. Liu.
*
************************************************************************
************************************************************************
*****   BETREE ... binary tree representation of elimination tree  *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine determines the binary tree representation of
*     the elimination tree given by the PARENT vector.  The returned
*     representation will be given by the first-son and brother
*     vectors.  the root of the binary tree is always N.
*
*   ----------
*   Arguments:
*   ----------
*
*     N         (input) integer
*               Number of equations.
*
*     PARENT    (input) integer array, dimension N
*               The parent vector of the elimination tree of the
*               symmetric matrix.  It is assumed that PARENT(I) > I
*               except of the roots.
*
*     FSON      (output) integer array, dimension N
*               The first son vector.
*
*     BROTHR    (output) integer array, dimension N
*               The brother vector.
*
************************************************************************
*
      SUBROUTINE  BETREE  ( N     , PARENT, FSON  , BROTHR          )
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
        INTEGER             PARENT(*)
        INTEGER             FSON(*), BROTHR(*)
*
************************************************************************
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             LROOT , NDPAR , NODE
*
************************************************************************
*
        IF  ( N .LE. 0 )  RETURN
*
*       ---------------
*       Initialization.
*       ---------------
        DO  NODE = 1, N
            FSON(NODE)   = 0
            BROTHR(NODE) = 0
        END DO
        LROOT = N
*       --------------------------------------------------------
*       For each NODE := N-1 step -1 downto 1, do the following.
*       --------------------------------------------------------
        IF  ( N .LE. 1 )  RETURN
        DO  NODE = N-1, 1, -1
            NDPAR = PARENT(NODE)
            IF  ( NDPAR .LE. 0  .OR.  NDPAR .EQ. NODE )  THEN
*               -------------------------------------------------
*               NODE has no parent.  Given structure is a forest.
*               Set node to be one of the roots of the trees.
*               -------------------------------------------------
                BROTHR(LROOT) = NODE
                LROOT         = NODE
            ELSE
*               -------------------------------------------
*               Otherwise, becomes first son of its parent.
*               -------------------------------------------
                BROTHR(NODE)  = FSON(NDPAR)
                FSON(NDPAR)   = NODE
            END IF
        END DO
        BROTHR(LROOT) = 0
*
        RETURN
*
*       --------------
*       End of BETREE.
*       --------------
      END
