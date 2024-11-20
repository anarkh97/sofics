************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    Esmond G. Ng and Barry W. Peyton
*   Created:    02-13-1995 (Barry W. Peyton)
*   Modified:   09-23-1998 (Esmond G. Ng)
*
*   Computational Mathematics and Statistics Section
*   Oak Ridge National Laboratory
*
************************************************************************
************************************************************************
*****   SYMFCT ... symbolic Cholesky factorization                 *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine calls SYMFC2 which performs supernodal symbolic
*     factorization on a reordered symmetric matrix.
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
*     PERM      (input) integer array, dimension N
*               The (postordered) permutation vector.
*
*     INVP      (input) integer array, dimension N
*               The inverse of the permutation vector.
*
*     COLCNT    (input) integer array, dimension N
*               The number of nonzeros in each column of the Cholesky
*               factor, including the diagonal entry.
*
*     NSUPER    (input) integer
*               Number of supernodes.
*
*     XSUPER    (input) integer array, dimension NSUPER+1
*               The supernodal partition.
*
*     SNODE     (input) integer array, dimension N
*               Membership of columns in the supernode partition.
*
*     NOFSUB    (input) integer  
*               Number of row indices in the compressed structure
*               representation.
*
*     XLINDX    (output) integer array, dimension NSUPER+1
*               Pointers to column structure of the Cholesky factor.
*
*     LINDX     (output) integer array, dimension NOFSUB
*               Row indices of nonzero entries in the Cholesky factor,
*               stored by columns, in a compressed representation.
*
*     XLNZ      (output) integer array, dimension N+1
*               Pointers to nonzero entries in the Cholesky factor.
*
*     IWSIZE    (input) integer
*               Size of work array IWORK.
*
*     IWORK     (temporary) integer array, dimension NSUPER+2*N+1
*               Hold various work arrays.
*
*     IFLAG     (output) integer
*               Error flag.
*
*               IFLAG =  0: no error.
*               IFLAG = 22: insufficient work space in IWORK.
*               IFLAG = 23: inconsistancy in the input.
*
************************************************************************
*
      SUBROUTINE  SYMFCT  ( N     , NADJ  , XADJ  , ADJNCY, PERM  ,
     &                      INVP  , COLCNT, NSUPER, XSUPER, SNODE ,
     &                      NOFSUB, XLINDX, LINDX , XLNZ  , IWSIZ ,
     &                      IWORK , IFLAG                           )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             IFLAG , IWSIZ , N     , NADJ  , NOFSUB,
     &                      NSUPER
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             XADJ(*), ADJNCY(*)
        INTEGER             PERM(*), INVP(*)
        INTEGER             COLCNT(*)
        INTEGER             XSUPER(*)
        INTEGER             SNODE(*)
        INTEGER             XLINDX(*), LINDX(*)
        INTEGER             XLNZ(*)
        INTEGER             IWORK(*)
*
************************************************************************
*
*       ------------------------
*       External Subroutines ...
*       ------------------------
        EXTERNAL            SYMFC2
*
************************************************************************
*
        IFLAG = 0
        IF  ( IWSIZ .LT. NSUPER+2*N+1 )  THEN
            IFLAG = 22
            RETURN
        END IF
        CALL  SYMFC2 ( N, NADJ, XADJ, ADJNCY, PERM, INVP, COLCNT,
     &                 NSUPER, XSUPER, SNODE, NOFSUB, XLINDX, LINDX,
     &                 XLNZ,
     &                 IWORK(1),
     &                 IWORK(NSUPER+1),
     &                 IWORK(NSUPER+N+2),
     &                 IFLAG )
        RETURN
*
*       --------------
*       End of SYMFCT.
*       --------------
      END
