************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    Esmond G. Ng and Barry W. Peyton
*   Created:    11-14-1994 (Barry W. Peyton)
*   Modified:   01-12-1995 (Barry W. Peyton)
*   Modified:   09-23-1998 (Esmond G. Ng)
*
*   Computational Mathematics and Statistics Section
*   Oak Ridge National Laboratory
*
************************************************************************
************************************************************************
*****   SFINIT ... set up for symbolic factorization               *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*     This subroutine computes the storage requirements and sets up 
*     preliminary data structures for the symbolic factorization.
*
*     Note:
*       This version produces the maximal supernode partition (i.e.,
*       the one with the fewest possible supernodes).
*
*   ----------
*   Arguments:
*   ----------
*
*     N         (input) integer
*               Number of equations.
*
*     NADJ      (input) integer
*               Number of edges in the adjacency structure.
*
*     XADJ      (input) integer array, dimension N+1
*               Pointers to the adjacency structure.
*
*     ADJNCY    (input) integer array, dimension XADJ(N+1)-1
*               The adjacency structure.
*
*     PERM      (input/output) integer array, dimension N
*               The (postordered) permutation vector.
*
*     INVP      (input/output) integer array, dimension N
*               The inverse of the permutation vector.
*
*     MAXSUP    (input) integer
*               The maximum number of columns in each supernode.
*
*     DEFBLK    (input) integer
*               The size of the last supernode, which contains the
*               known deficiency.
*
*     COLCNT    (output) integer array, dimension N
*               The number of nonzeros in each column of the Cholesky
*               factor, including the diagonal entry.
*
*     NNZL      (output) integer
*               Number of nonzero entries in the Cholesky factor,
*               including the diagonal entries.
*
*     NSUB      (output) integer
*               Number of row indices in the compressed structure
*               representation.
*
*     NSUPER    (output) integer
*               Number of supernodes.
*
*     XSUPER    (output) integer array, dimension NSUPER+1
*               The supernodal partition.
*
*     SNODE     (output) integer array, dimension N
*               Membership of columns in the supernode partition.
*
*     IWSIZE    (input) integer
*               Size of work array IWORK.
*
*     IWORK     (temporary) integer array, dimension 7*N+3
*               Hold various work arrays.
*
*     IFLAG     (output) integer
*               Error flag.
*
*               IFLAG =  0: successful symbolic factorization
*                           initialization.
*               IFLAG = 21: insufficent work space in IWORK.
*
************************************************************************
*
      SUBROUTINE  SFINIT  ( N     , NADJ  , XADJ  , ADJNCY, PERM  ,
     &                      INVP  , MAXSUP, DEFBLK, COLCNT, NNZL  ,
     &                      NSUB  , NSUPER, XSUPER, SNODE , IWSIZE,
     &                      IWORK , IFLAG                           )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             DEFBLK, IFLAG , IWSIZE, MAXSUP, N     ,
     &                      NADJ  , NNZL  , NSUB  , NSUPER
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             XADJ(*), ADJNCY(*)
        INTEGER             INVP(*), PERM(*)
        INTEGER             XSUPER(*), SNODE(*)
        INTEGER             COLCNT(*)
        INTEGER             IWORK(*)
*
************************************************************************
*
*       ------------------------
*       External Subroutines ...
*       ------------------------
        EXTERNAL            CHORDR, ETORDR, FCNTHN, FSUP1 , FSUP2
*
************************************************************************
*
*       ----------------------------------------------------
*       Return if there is insufficient work space in IWORK.
*       ----------------------------------------------------
        IFLAG = 0
        IF  ( IWSIZE .LT. 7*N+3 )  THEN
            IFLAG = 1
            RETURN
        END IF
*
*       ------------------------------------------
*       Compute elimination tree and postordering.
*       ------------------------------------------
        CALL  ETORDR ( N, XADJ, ADJNCY, PERM, INVP,
     &                 IWORK(1),
     &                 IWORK(N+1),
     &                 IWORK(2*N+1),
     &                 IWORK(3*N+1) ) 
*
*       ---------------------------------------------------------
*       Compute row and column nonzero counts in Cholesky factor.
*       ---------------------------------------------------------
        CALL  FCNTHN ( N, NADJ, XADJ, ADJNCY, PERM, INVP,
     &                 IWORK(1), SNODE, COLCNT, NNZL,
     &                 IWORK(N+1),
     &                 IWORK(2*N+1),
     &                 XSUPER,
     &                 IWORK(3*N+1),
     &                 IWORK(4*N+2),
     &                 IWORK(5*N+3),
     &                 IWORK(6*N+4) )
*
*       ---------------------------------------------------------
*       Rearrange children so that the last child has the maximum 
*       number of nonzeros in its column of the Cholesky factor.
*       ---------------------------------------------------------
        CALL  CHORDR ( N, PERM, INVP, COLCNT, 
     &                 IWORK(1),
     &                 IWORK(N+1),
     &                 IWORK(2*N+1),
     &                 IWORK(3*N+1) )
*
*       ----------------
*       Find supernodes.
*       ----------------
        CALL  FSUP1 ( MAXSUP, DEFBLK, N, IWORK, COLCNT, NSUB,
     &                NSUPER, SNODE )
        CALL  FSUP2 ( N, NSUPER, SNODE, XSUPER )
        RETURN
*
*       --------------
*       End of SFINIT.
*       --------------
      END
