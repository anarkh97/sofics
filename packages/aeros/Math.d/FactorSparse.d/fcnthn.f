************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    John R. Gilbert, Esmond G. Ng, and Barry W. Peyton
*   Created:    04-12-1990 (Barry W. Peyton)
*   Modified:   01-12-1995 (Barry W. Peyton)
*   Modified:   09-23-1998 (Esmond G. Ng)
*
*   Computational Mathematics and Statistics Section
*   Oak Ridge National Laboratory
*
************************************************************************
************************************************************************
*****   FCNTHN ... find nonzero counts                             *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine determines the row counts and column counts in
*     the Cholesky factor.  It uses a disjoint set union algorithm.
*
*     -----------
*     Techniques:
*     -----------
*       1) Supernode detection.
*       2) Path halving.
*       3) No union by rank.
*
*     ------
*     Notes:
*     ------
*       1) Assumes a postordering of the elimination tree.
*
*   ----------
*   Arguments:
*   ----------
*
*     N         (input) integer
*               Number of equations.
*
*     NADJ      (output) integer
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
*     ETPAR     (input) integer array, dimension N
*               The elimination tree of the postordered matrix.
*
*     ROWCNT    (output) integer array, dimension N
*               It containing the number of nonzero entries in each
*               row of the Cholesky factor, including the diagonal
*               entry.
*
*     COLCNT    (output) integer array, dimension N
*               It containing the number of nonzero entries in each
*               column of the Cholesky factor, including the diagonal
*               entry.
*
*     NNZL      (output) integer
*               Number of nonzero entries in the Cholesky factor,
*               including the diagonal entries.
*
*     SET       (temporary) integer array, dimension N
*               It maintains the disjoint sets (i.e., subtrees).
*
*     PRVLF     (temporary) integer array, dimension N
*               It records the previous leaf of each row subtree.
*
*     LEVEL     (temporary) integer array, dimension N+1
*               It contains the level numbers (distance from the
*               root of the elimination tree).
*
*     WEIGHT    (temporary) integer array, dimension N+1
*               It contains weights used to compute column counts.
*
*     FDESC     (temporary) integer array, dimension N+1
*               It contains the first (i.e., lowest-numbered)
*               descendant.
*
*     NCHILD    (temporary) integer array, dimension N+1
*               It contains the number of children.
*
*     PRVNBR    (temporary) integer array, dimension N
*               It records the previous "lower neighbor" of each node.
*
************************************************************************
*
      SUBROUTINE  FCNTHN  ( N     , NADJ  , XADJ  , ADJNCY, PERM  ,
     &                      INVP  , ETPAR , ROWCNT, COLCNT, NNZL  ,
     &                      SET   , PRVLF , LEVEL , WEIGHT, FDESC ,
     &                      NCHILD, PRVNBR                          )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             NADJ  , N     , NNZL
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             XADJ(N+1), ADJNCY(NADJ)
        INTEGER             PERM(N), INVP(N)
        INTEGER             ETPAR(N)
        INTEGER             ROWCNT(N), COLCNT(N)
        INTEGER             FDESC(0:N), LEVEL(0:N), NCHILD(0:N),
     &                      PRVLF(N), PRVNBR(N), SET(N), WEIGHT(0:N)
*
************************************************************************
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             HINBR , IFDESC, J     , JSTOP , JSTRT ,
     &                      K     , LAST1 , LAST2 , LCA   , LFLAG ,
     &                      LOWNBR, OLDNBR, PARENT, PLEAF , TEMP  ,
     &                      XSUP
*
************************************************************************
*
*       --------------------------------------------------
*       Compute LEVEL(*), FDESC(*), NCHILD(*).
*       Initialize XSUP, ROWCNT(*), COLCNT(*).
*       Initialize SET(*), PRVLF(*), WEIGHT(*), PRVNBR(*).
*       --------------------------------------------------
        XSUP     = 1
        LEVEL(0) = 0
        DO  K = N, 1, -1
            ROWCNT(K) = 1
            COLCNT(K) = 0
            SET(K)    = K
            PRVLF(K)  = 0
            LEVEL(K)  = LEVEL(ETPAR(K)) + 1
            WEIGHT(K) = 1
            FDESC(K)  = K
            NCHILD(K) = 0
            PRVNBR(K) = 0
        END DO
        NCHILD(0) = 0
        FDESC(0)  = 0
        DO  K = 1, N
            PARENT         = ETPAR(K)
            WEIGHT(PARENT) = 0
            NCHILD(PARENT) = NCHILD(PARENT) + 1
            IFDESC         = FDESC(K)
            IF  ( IFDESC .LT. FDESC(PARENT) )  THEN
                FDESC(PARENT) = IFDESC
            END IF
        END DO
*       ----------------------------------
*       For each "low neighbor" LOWNBR ...
*       ----------------------------------
        DO  LOWNBR = 1, N
            LFLAG  = 0
            IFDESC = FDESC(LOWNBR)
            OLDNBR = PERM(LOWNBR)
            JSTRT  = XADJ(OLDNBR)
            JSTOP  = XADJ(OLDNBR+1) - 1
*           --------------------------------------------
*           For each "high neighbor" HINBR of LOWNBR ...
*           --------------------------------------------
            DO  J = JSTRT, JSTOP
                HINBR = INVP(ADJNCY(J))
                IF  ( HINBR .GT. LOWNBR )  THEN
                    IF  ( IFDESC .GT. PRVNBR(HINBR) )  THEN
*                       -------------------------
*                       Increment WEIGHT(LOWNBR).
*                       -------------------------
                        WEIGHT(LOWNBR) = WEIGHT(LOWNBR) + 1
                        PLEAF          = PRVLF(HINBR)
*                       ---------------------------------------
*                       If HINBR has no previous "LOW NEIGHBOR"
*                       then ...
*                       ---------------------------------------
                        IF  ( PLEAF .EQ. 0 )  THEN
*                           -----------------------------------------
*                           ... accumulate LOWNBR-->HINBR path length
*                               in ROWCNT(HINBR).
*                           -----------------------------------------
                            ROWCNT(HINBR) = ROWCNT(HINBR) +
     &                              LEVEL(LOWNBR) - LEVEL(HINBR)
                        ELSE
*                           -----------------------------------------
*                           ... otherwise, LCA <-- FIND(PLEAF), which
*                               is the least common ancestor of PLEAF
*                               and LOWNBR.  (PATH HALVING.)
*                           -----------------------------------------
                            LAST1 = PLEAF
                            LAST2 = SET(LAST1)
                            LCA   = SET(LAST2)
  300                       CONTINUE
                                IF  ( LCA .NE. LAST2 )  THEN
                                    SET(LAST1) = LCA
                                    LAST1 = LCA
                                    LAST2 = SET(LAST1)
                                    LCA   = SET(LAST2)
                                    GO TO 300
                                END IF
*                           --------------------------------------
*                           Accumulate PLEAF-->LCA path length in
*                           ROWCNT(HINBR).  Decrement WEIGHT(LCA).
*                           --------------------------------------
                            ROWCNT(HINBR) = ROWCNT(HINBR) +
     &                              LEVEL(LOWNBR) - LEVEL(LCA)
                            WEIGHT(LCA)   = WEIGHT(LCA) - 1
                        END IF
*                       --------------------------------------------
*                       LOWNBR now becomes "previous leaf" of HINBR.
*                       --------------------------------------------
                        PRVLF(HINBR) = LOWNBR
                        LFLAG        = 1
                    END IF
*                   ------------------------------------------------
*                   LOWNBR now becomes "previous neighbor" of HINBR.
*                   ------------------------------------------------
                    PRVNBR(HINBR) = LOWNBR
                END IF
            END DO
*           ----------------------------------------------
*           Decrement WEIGHT(PARENT(LOWNBR)).
*           SET(P(LOWNBR)) <-- SET(P(LOWNBR)) + SET(XSUP).
*           ----------------------------------------------
            PARENT         = ETPAR(LOWNBR)
            WEIGHT(PARENT) = WEIGHT(PARENT) - 1
            IF  ( LFLAG .EQ. 1  .OR.  NCHILD(LOWNBR) .GE. 2 )  THEN
                XSUP = LOWNBR
            END IF
            SET(XSUP) = PARENT
        END DO
*       ---------------------------------------------------------
*       Use weights to compute column (and total) nonzero counts.
*       ---------------------------------------------------------
        NNZL = 0
        DO  K = 1, N
            TEMP      = COLCNT(K) + WEIGHT(K)
            COLCNT(K) = TEMP
            NNZL      = NNZL + TEMP
            PARENT    = ETPAR(K)
            IF  ( PARENT .NE. 0 )  THEN
                COLCNT(PARENT) = COLCNT(PARENT) + TEMP
            END IF
        END DO
*
        RETURN
*
*       --------------
*       End of FCNTHN.
*       --------------
      END
