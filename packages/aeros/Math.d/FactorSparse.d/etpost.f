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
*       This subroutine is kindly provided by Joseph W.H. Liu.
*
************************************************************************
************************************************************************
*****   ETPOST ... elimination tree postordering                   *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*     This subroutine determines a postordering from the binary
*     representation (first-son,brother) of the elimination tree.
*     The corresponding PARENT vector is also modified to reflect
*     the reordering.
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
*     INVPOS    (input) integer array, dimension N
*               The inverse of the postordered permutation.
*
*     PARENT    (input/output) integer array, dimension N
*               The elimination tree associated with a postordering
*               of the symmetric matrix.
*
*     STACK     (temporary) integer array, dimension N
*               The stack for postorder traversal of the elimination
*               tree.
*
************************************************************************
*
      SUBROUTINE  ETPOST  ( ROOT  , FSON  , BROTHR, INVPOS, PARENT,
     &                      STACK                                   )
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
*       Determine the new parent vector of the postordering.  BROTHR
*       is used temporarily for the new parent vector.
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
        RETURN
*
*       --------------
*       End of ETPOST.
*       --------------
      END
