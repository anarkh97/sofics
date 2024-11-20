************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    Esmond G. Ng and Barry W. Peyton
*   Created:    12-27-1994 (Barry W. Peyton)
*   Modified:   09-23-1998 (Esmond G. Ng)
*
*   COMPUTATIONAL MATHEMATICS AND SATISTICS SECTION
*   Oak Ridge National Laboratory
*
************************************************************************
************************************************************************
*****   CHORDR ... child reordering                                *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine rearranges the children of each vertex so
*     that the last one maximizes (among the children) the number
*     of nonzeros in the corresponding column of the Cholesky
*     factor.  Also determine a new postordering based on the
*     structure of the modified elimination tree.
*
*   ----------
*   Arguments:
*   ----------
*
*     N         (input) integer
*               Number of equations.
*
*     PERM      (input/output) integer array, dimension N
*               The permutation vector.  On output, PERM contains
*               an equivalent reordering, which is a postordering
*               of the elimination tree.
*
*     INVP      (input/output) integer array, dimension N
*               The inverse of the permutation vector.
*
*     COLCNT    (input/output) integer array, dimension N
*               The number of nonzeros in each column of the Cholesky
*               factor, including the diagonal entry.
*
*     PARENT    (output) integer array, dimension N
*               The parent vector of the elimination tree associated
*               with the new postordering of the symmetric matrix.
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
      SUBROUTINE  CHORDR  ( N     , PERM  , INVP  , COLCNT, PARENT,
     &                      FSON  , BROTHR, INVPOS                  )
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
        INTEGER             PERM(*), INVP(*)
        INTEGER             INVPOS(*)
        INTEGER             PARENT(*)
        INTEGER             COLCNT(*)
        INTEGER             BROTHR(*), FSON(*)
*
*       ------------------------
*       External Subroutines ...
*       ------------------------
        EXTERNAL            BTREE2, EPOST2, INVINV
*
************************************************************************
*
*       ----------------------------------------------------------
*       Compute a binary representation of the elimination tree, 
*       so that each "last child" maximizes among its siblings the 
*       number of nonzero entries in the corresponding columns of
*       Cholesky factor.
*       ----------------------------------------------------------
        CALL  BTREE2 ( N, PARENT, COLCNT, FSON, BROTHR, INVPOS )
*
*       ----------------------------------------------------
*       Postorder the elimination tree (using the new binary  
*       representation.  
*       ----------------------------------------------------
        CALL  EPOST2 ( N, FSON, BROTHR, INVPOS, PARENT, COLCNT,
     &                 PERM ) 
*
*       --------------------------------------------------------
*       Combine the original ordering with the new postordering.
*       --------------------------------------------------------
        CALL  INVINV ( N, INVP, INVPOS, PERM )
*
        RETURN
*
*       --------------
*       End of CHORDR.
*       --------------
      END
