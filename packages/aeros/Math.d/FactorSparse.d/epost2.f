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
*****   EPOST2 ... elimination tree postordering #2                *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine determines a postordering from the binary
*     representation (first-son,brother) of the elimination tree.
*     The corresponding PARENT and COLCNT vectors are also modified
*     to reflect the reordering.
*
*   ----------
*   Arguments:
*   ----------
*
*     ROOT      (input) integer
*               Root of the elimination tree (usually it is N).
*
*     FSON      (input) integer array, dimension N
*               The first son vector.
*
*     BROTHR    (input) integer array, dimension N
*               The brother vector.
*
*     INVPOS    (output) integer array, dimension N
*               The inverse of the postordered permutation.
*
*     PARENT    (input/output) integer array, dimension N
*               The elimination tree associated with a postordering
*               of the symmetric matrix.
*
*     COLCNT    (input/output) integer array, dimension N
*               The number of nonzeros in each column of the Cholesky
*               factor, including the diagonal entry.
*
*     STACK     (temporary) integer array, dimension N
*               The stack for postorder traversal of the elimination
*               tree.
*
************************************************************************
*
      SUBROUTINE  EPOST2  ( ROOT  , FSON  , BROTHR, INVPOS, PARENT,
     &                      COLCNT, STACK                           )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             ROOT
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             PARENT(*)
        INTEGER             COLCNT(*)
        INTEGER             INVPOS(*)
        INTEGER             FSON(*), BROTHR(*)
        INTEGER             STACK(*)
*
************************************************************************
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             ITOP  , NDPAR , NODE  , NUM   , NUNODE
*
************************************************************************
*
        NUM  = 0
        ITOP = 0
        NODE = ROOT
*       -------------------------------------------------------------
*       Traverse along the first sons pointer and push the tree nodes
*       along the traversal into the stack.
*       -------------------------------------------------------------
  100   CONTINUE
            ITOP        = ITOP + 1
            STACK(ITOP) = NODE
            NODE        = FSON(NODE)
            IF  ( NODE .GT. 0 )  GO TO 100
*           ----------------------------------------------------------
*           If possible, pop a tree node from the stack and number it.
*           ----------------------------------------------------------
  200       CONTINUE
                IF  ( ITOP .LE. 0 )  GO TO 300
                NODE         = STACK(ITOP)
                ITOP         = ITOP - 1
                NUM          = NUM + 1
                INVPOS(NODE) = NUM
*               ----------------------------------------------------
*               Then, traverse to its younger brother if it has one.
*               ----------------------------------------------------
                NODE         = BROTHR(NODE)
                IF  ( NODE .LE. 0 )  GO TO 200
            GO TO 100
*
  300   CONTINUE
*       ------------------------------------------------------------
*       Determine the new PARENT vector of the postordering.  BROTHR
*       is used temporarily for the new PARENT vector.
*       ------------------------------------------------------------
        DO  NODE = 1, NUM
            NUNODE         = INVPOS(NODE)
            NDPAR          = PARENT(NODE)
            IF  ( NDPAR .GT. 0 )  NDPAR = INVPOS(NDPAR)
            BROTHR(NUNODE) = NDPAR
        END DO
*
        DO  NUNODE = 1, NUM
            PARENT(NUNODE) = BROTHR(NUNODE)
        END DO
*
*       ----------------------------------------------
*       Permute COLCNT(*) TO reflect the new ordering.
*       ----------------------------------------------
        DO  NODE = 1, NUM
            NUNODE        = INVPOS(NODE)
            STACK(NUNODE) = COLCNT(NODE)
        END DO
*
        DO  NODE = 1, NUM
            COLCNT(NODE) = STACK(NODE)
        END DO
*
        RETURN
*
*       --------------
*       End of EPOST2.
*       --------------
      END
