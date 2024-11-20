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
*****   ORDMMD2 ... multiple minimum external degree                *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine calls Liu'S multiple minimum degree routine.
*
*   ----------
*   Arguments:
*   ----------
*
*     N         (input) integer
*               Number of equations.
*
*     XADJ      (input/output) integer array, dimension N+1
*               Pointers to the adjacency structure.
*
*     ADJNCY    (input/output) integer array, dimension XADJ(N+1)-1
*               The adjacency structure.
*
*     INVP      (output) integer array, dimension N
*               The inverse of the permutation vector.
*   
*     PERM      (output) integer array, dimension N
*               The permutation vector.
*
*     IWSIZE    (input) integer
*               Size of work array IWORK.
*
*     IWORK     (temporary) integer array, dimension 4*N
*               Hold various work arrays.
*
*     NOFSUB    (output) integer
*               An upper bound on the number of row indices for the
*               compressed storage scheme.
*
*     IFLAG     (output) integer
*               Error flag.
*
*               IFLAG =  0: successful ordering.
*               IFLAG = 11: insufficient work space in IWORK.
*
************************************************************************
*
      SUBROUTINE  ORDMMD2  ( N     , XADJ  , ADJNCY, INVP  , PERM  ,
     *                      IWSIZE, IWORK , NOFSUB, IFLAG           )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             IFLAG , IWSIZE, N     , NOFSUB
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             XADJ(*), ADJNCY(*)
        INTEGER             PERM(*), INVP(*)
        INTEGER             IWORK(*)
*
**********************************************************************
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             DELTA , MAXINT
*
*       ------------------------
*       External Subroutines ...
*       ------------------------
        EXTERNAL            GENMMD
*
**********************************************************************
*
        IFLAG = 0
        IF  ( IWSIZE .LT. 4*N )  THEN
            IFLAG = 11
            RETURN
        END IF
*
*       ------------------------------------------------------
*       DELTA   -   Tolerance value for multiple elimination.
*       MAXINT  -   Maximum machine representable integer
*                   (any smaller estimate will do) for marking
*                   nodes.
*       ------------------------------------------------------
*
        DELTA  = 0
        MAXINT = 2 000 000 000
        CALL GENMMD  (  N, XADJ, ADJNCY, INVP, PERM, DELTA, 
     &                  IWORK(1),
     &                  IWORK(N+1),
     &                  IWORK(2*N+1),
     &                  IWORK(3*N+1),
     &                  MAXINT, NOFSUB )
*
        RETURN
*
*       --------------
*       End of ORDMMD2.
*       --------------
      END
