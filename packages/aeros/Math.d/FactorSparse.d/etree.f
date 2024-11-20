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
*****   ETREE ... elimination tree                                 *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine determines the elimination tree from a given
*     ordering and the adjacency structure of a sparse symmetric
*     matrix.  The PARENT vector is returned.
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
*     PERM      (input) integer array, dimension N
*               The permutation vector.
*
*     INVP      (input) integer array, dimension N
*               The inverse of the permutation vector.
*
*     PARENT    (output) integer array, dimension N
*               The elimination tree of the symmetric matrix.
*
*     ANCSTR    (temporary) integer array, dimension N
*               The ancestor array keeps track what current roots of
*               various subtrees are.
*
************************************************************************
*
      SUBROUTINE  ETREE   ( N     , XADJ  , ADJNCY, PERM  , INVP  ,
     &                      PARENT, ANCSTR                          )
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
        INTEGER             ANCSTR(*)
*
************************************************************************
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             I     , J     , JBEG  , JEND  , NBR   ,
     &                      NEXT  , NODE
*
************************************************************************
*
        IF  ( N .LE. 0 )  RETURN
*
        DO  I = 1, N
            PARENT(I) = 0
            ANCSTR(I) = 0
            NODE      = PERM(I)
*
            JBEG = XADJ(NODE)
            JEND = XADJ(NODE+1) - 1
            IF  ( JBEG .LE. JEND )  THEN
                DO  J = JBEG, JEND
                    NBR = ADJNCY(J)
                    NBR = INVP(NBR)
                    IF  ( NBR .LT. I )  THEN
*                       -------------------------------------------
*                       For each NBR, find the root of its current
*                       elimination tree.  Perform path compression
*                       as the subtree is traversed.
*                       -------------------------------------------
  100                   CONTINUE
                            IF  ( ANCSTR(NBR) .EQ. I )  GO TO 200
                            IF  ( ANCSTR(NBR) .GT. 0 )  THEN
                                NEXT        = ANCSTR(NBR)
                                ANCSTR(NBR) = I
                                NBR         = NEXT
                                GO TO 100
                            END IF
*                       --------------------------------------------
*                       Now, NBR is the root of the subtree.  Make I
*                       the parent node of this root.
*                       --------------------------------------------
                        PARENT(NBR) = I
                        ANCSTR(NBR) = I
                    END IF
  200               CONTINUE
                END DO
            END IF
        END DO
*
        RETURN
*
*       -------------
*       End of ETREE.
*       -------------
      END
