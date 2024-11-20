************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    Esmond G. Ng and Barry W. Peyton
*   Created:    01-12-1995 (Barry W. Peyton)
*   Modified:   09-23-1998 (Esmond G. Ng)
*
*   Computational Mathematics and Statistics Section
*   Oak Ridge National Laboratory
*
************************************************************************
************************************************************************
*****   BTREE2 ... binary tree representation of elimination tree  *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine determines a binary tree representation of
*     the elimination tree, for which every "last child" has the
*     maximum possible column nonzero count in the factor.  The
*     returned representation will be given by the first-son and
*     brother vectors.  The root of the binary tree is always N.
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
*     COLCNT    (input) integer array, dimension N
*               The number of nonzeros in each column of the Cholesky
*               factor, including the diagonal entry.
*
*     FSON      (output) integer array, dimension N
*               The first son vector.
*
*     BROTHR    (output) integer array, dimension N
*               The brother vector.
*
*     LSON      (temporary) integer array, dimension N
*               The last son vector.
*
************************************************************************
*
      SUBROUTINE  BTREE2  ( N     , PARENT, COLCNT, FSON  , BROTHR,
     &                      LSON                                    )
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
        INTEGER             COLCNT(*)
        INTEGER             FSON(*), BROTHR(*)
        INTEGER             LSON(*)
*
************************************************************************
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             LROOT , NDLSON, NDPAR , NODE
*
************************************************************************
*
        IF  ( N .LE. 0 )  RETURN
*
        DO  NODE = 1, N
            FSON(NODE)   = 0
            BROTHR(NODE) = 0
            LSON(NODE)   = 0
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
*               Set NODE to be one of the roots of the trees.
*               -------------------------------------------------
                BROTHR(LROOT) = NODE
                LROOT         = NODE
            ELSE
*               -------------------------------------------
*               Otherwise, becomes first son of its parent.
*               -------------------------------------------
                NDLSON = LSON(NDPAR)
                IF  ( NDLSON .NE. 0 )  THEN
                    IF  ( COLCNT(NODE) .GE. COLCNT(NDLSON) )  THEN
                        BROTHR(NODE)   = FSON(NDPAR)
                        FSON(NDPAR)    = NODE
                    ELSE
                        BROTHR(NDLSON) = NODE
                        LSON(NDPAR)    = NODE
                    END IF
                ELSE
                    FSON(NDPAR) = NODE
                    LSON(NDPAR) = NODE
                END IF
            END IF
        END DO
        BROTHR(LROOT) = 0
*
        RETURN
*
*       --------------
*       End of BTREE2.
*       --------------
      END
