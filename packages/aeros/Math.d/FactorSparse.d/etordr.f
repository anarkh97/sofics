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
*****   ETORDR ... elimination tree reordering                     *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine determines an equivalent reordering based on
*     the structure of the elimination tree.  A postordering of the
*     given elimination tree is returned.
*
*   ----------
*   Arguments:
*   ----------
*
*     N         (input) integer
*               Number of equations.
*
*     XADJ      (input) integer array, dimension N+1
*               Pointers to the adjacency structure.
*
*     ADJNCY    (input) integer array, dimension XADJ(N+1)-1
*               The adjacency structure.
*
*     PERM      (input/output) integer array, dimension N
*               The permutation vector.  On output, PERM contains
*               an equivalent reordering, which is a postordering
*               of the elimination tree.
*
*     INVP      (input/output) integer array, dimension N
*               The inverse of the permutation vector.
*
*     PARENT    (output) integer array, dimension N
*               The elimination tree associated with the equivalent
*               reordering (which is a postordering) of the symmetric
*               matrix.
*
*     FSON      (temporary) integer array, dimension N
*               The first son vector.
*
*     BROTHR    (temporary) integer array, dimension N
*               The brother vector.
*
*     INVPOS    (temporary) integer array, dimension N
*               The inverse of the postordered permutation.
*
************************************************************************
*
      SUBROUTINE  ETORDR  ( N     , XADJ  , ADJNCY, PERM  , INVP  ,
     &                      PARENT, FSON  , BROTHR, INVPOS          )
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
        INTEGER             XADJ(*), ADJNCY(*)
        INTEGER             PERM(*), INVP(*)
        INTEGER             PARENT(*) 
        INTEGER             INVPOS(*)
        INTEGER             FSON(*), BROTHR(*)
*
*       ------------------------
*       External Subroutines ...
*       ------------------------
        EXTERNAL            BETREE, ETPOST, ETREE , INVINV
*
************************************************************************
*
*       -----------------------------
*       Compute the elimination tree.
*       -----------------------------
        CALL  ETREE ( N, XADJ, ADJNCY, PERM, INVP, PARENT, INVPOS )
*
*       --------------------------------------------------------
*       Compute a binary representation of the elimination tree.
*       --------------------------------------------------------
        CALL  BETREE ( N, PARENT, FSON, BROTHR )
*
*       -------------------------------
*       Postorder the elimination tree.
*       -------------------------------
        CALL  ETPOST ( N, FSON, BROTHR, INVPOS, PARENT, PERM )
*
*       --------------------------------------------------------
*       Combine the original ordering with the new postordering.
*       --------------------------------------------------------
        CALL  INVINV ( N, INVP, INVPOS, PERM )
*
        RETURN
*
*       --------------
*       End of ETORDR.
*       --------------
      END
