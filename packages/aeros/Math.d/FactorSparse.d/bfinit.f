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
*****   BFINIT ... initialization for block Cholesky factorization *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine computes items needed by the left-looking
*     block-to-block Cholesky factoritzation routine BLKLDL.
*
*   ----------
*   Arguments:
*   ----------
*
*     NSUPER    (input) integer
*               Number of supernodes.  (N = XSUPER(NSUPER+1)-1.)
*
*     XSUPER    (input) integer array, dimension NSUPER+1
*               The supernodal partition.
*
*     SNODE     (input) integer array, dimension N
*               Membership of columns in the supernode partition.
*
*     XLINDX    (input) integer array, dimension NSUPER+1
*               Pointers to column structure of the Cholesky factor.
*
*     LINDX     (input) integer array, dimension XSUPER(NSUPER+1)-1
*               Row indices of nonzero entries in the Cholesky factor,
*               stored by columns, in a compressed representation.
*
*     TMPSIZ    (output) integer
*               Size of work space in TEMP required by BLKLDL.
*
*     RWSIZE    (output) integer
*               Size of work space in RWORK required by BLKLDL.
*
************************************************************************
*
      SUBROUTINE  BFINIT  ( NSUPER, XSUPER, SNODE , XLINDX, LINDX ,
     &                      TMPSIZ, RWSIZE                          )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             NSUPER, RWSIZE, TMPSIZ
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             XLINDX(*), LINDX(*)
        INTEGER             XSUPER(*)
        INTEGER             SNODE(*)
*
*       ------------------------
*       External Subroutines ...
*       ------------------------
        EXTERNAL            FNTSIZ
*
************************************************************************
*
*       ------------------------------------------------
*       Determine floating point work space requirement.
*       ------------------------------------------------
        CALL  FNTSIZ ( NSUPER, XSUPER, SNODE, XLINDX, LINDX, TMPSIZ,
     &                 RWSIZE )
        RETURN
*
*       --------------
*       End of BFINIT.
*       --------------
      END
