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
*****   SYMFC2 ... symbolic Cholesky factorization                 *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose: 
*   --------
*
*     This subroutine performs supernodal symbolic factorization on
*     a reordered symmetric matrix.  It assumes access to the column
*     counts, supernode partition, and supernodal elimination tree
*     associated with the Cholesky factor matrix L.
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
*               Number of row indices for the compressed representation.
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
*     MRGLNK    (temporary) integer array, dimension NSUPER
*               It containing the children of each supernode as a
*               linked list.
*
*     RCHLNK    (temporary) integer array, dimension N+1
*               It containing the current linked list of merged indices
*               (the "REACH" set).
*
*     MARKER    (temporary) integer array, dimension N
*               It is used to mark indices as they are introduced
*               into each supernode's index set.
*
*     IFLAG     (output) integer
*               Error flag.
*
*               IFLAG =  0: no error.
*               IFLAG = 23: inconsistancy in the input.
*
************************************************************************
*
      SUBROUTINE  SYMFC2  ( N     , NADJ  , XADJ  , ADJNCY, PERM  , 
     &                      INVP  , COLCNT, NSUPER, XSUPER, SNODE ,
     &                      NOFSUB, XLINDX, LINDX , XLNZ  , MRGLNK,
     &                      RCHLNK, MARKER, IFLAG                   )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             IFLAG , N     , NADJ  , NOFSUB, NSUPER
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             XADJ(N+1), ADJNCY(NADJ)
        INTEGER             PERM(N), INVP(N)
        INTEGER             COLCNT(N)
        INTEGER             XSUPER(NSUPER+1)
        INTEGER             SNODE(N)
        INTEGER             XLINDX(NSUPER+1), LINDX(NOFSUB)
        INTEGER             XLNZ(N+1)
        INTEGER             MARKER(N)
        INTEGER             MRGLNK(NSUPER)
        INTEGER             RCHLNK(0:N)
*
************************************************************************
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             FSTCOL, HEAD  , I     , JCOL  , JNZBEG, 
     &                      JNZEND, JPTR  , JSUP  , JWIDTH, KNZ   , 
     &                      KNZBEG, KNZEND, KPTR  , KSUP  , LENGTH, 
     &                      LSTCOL, NEWI  , NEXTI , NODE  , NZBEG , 
     &                      NZEND , PCOL  , PSUP  , POINT , POINTL, 
     &                      TAIL  , WIDTH 
*
************************************************************************
*
        IFLAG = 0
        IF  ( N .LE. 0 )  RETURN
*
*       ---------------------------------------------------
*       Initializations ...
*           NZEND:   points to the last used slot in LINDX.
*           TAIL:    end of list indicator 
*                    (in RCHLNK(*), not MRGLNK(*)).
*           MRGLNK:  create empty lists.
*           MARKER:  "unmark" the indices.
*       ---------------------------------------------------
        NZEND = 0
        HEAD  = 0
        TAIL  = N + 1
        DO  I = 1, N
            MARKER(I) = 0
        END DO
        POINTL = 1
        DO  KSUP = 1, NSUPER
            FSTCOL = XSUPER(KSUP)
            LSTCOL = XSUPER(KSUP+1) - 1
            DO  JCOL = FSTCOL, LSTCOL
                XLNZ(JCOL) = POINTL
                POINTL     = POINTL + COLCNT(FSTCOL)
            END DO
        END DO
        XLNZ(N+1) = POINTL
        POINT     = 1
        DO  KSUP = 1, NSUPER
            MRGLNK(KSUP) = 0
            FSTCOL       = XSUPER(KSUP)
            XLINDX(KSUP) = POINT
            POINT        = POINT + COLCNT(FSTCOL)
        END DO
        XLINDX(NSUPER+1) = POINT
*
*       ---------------------------
*       For each supernode KSUP ... 
*       ---------------------------
        DO  KSUP = 1, NSUPER
*
*           --------------------------------------------------
*           Initializations ...
*               FSTCOL:  first column of supernode KSUP.
*               LSTCOL:  last column of supernode KSUP.
*               KNZ   :  will count the nonzero entries of L
*                        in column KCOL.
*               RCHLNK:  initialize empty index list for KCOL.
*           --------------------------------------------------
            FSTCOL       = XSUPER(KSUP)
            LSTCOL       = XSUPER(KSUP+1) - 1
            WIDTH        = LSTCOL - FSTCOL + 1
            LENGTH       = COLCNT(FSTCOL)
            KNZ          = 0
            RCHLNK(HEAD) = TAIL
            JSUP         = MRGLNK(KSUP)
*
*           --------------------------------------------------
*           If KSUP has children in the supernodal elimination
*           tree ...
*           --------------------------------------------------
            IF  ( JSUP .GT. 0 )  THEN
*               ---------------------------------------------
*               Copy the indices of the first child JSUP into 
*               the linked list, and mark each with the value 
*               KSUP.
*               ---------------------------------------------
                JWIDTH = XSUPER(JSUP+1) - XSUPER(JSUP)
                JNZBEG = XLINDX(JSUP) + JWIDTH
                JNZEND = XLINDX(JSUP+1) - 1
                DO  JPTR = JNZEND, JNZBEG, -1
                    NEWI         = LINDX(JPTR)
                    KNZ          = KNZ+1
                    MARKER(NEWI) = KSUP
                    RCHLNK(NEWI) = RCHLNK(HEAD)
                    RCHLNK(HEAD) = NEWI
                END DO
*               ------------------------------------------
*               For each subsequent child JSUP of KSUP ...
*               ------------------------------------------
                JSUP = MRGLNK(JSUP)
  100           CONTINUE
                IF  ( JSUP .NE. 0  .AND.  KNZ .LT. LENGTH )  THEN
*                   ----------------------------------------
*                   Merge the indices of JSUP into the list,
*                   and mark new indices with value KSUP.
*                   ----------------------------------------
                    JWIDTH = XSUPER(JSUP+1) - XSUPER(JSUP)
                    JNZBEG = XLINDX(JSUP) + JWIDTH
                    JNZEND = XLINDX(JSUP+1) - 1
                    NEXTI  = HEAD
                    DO  JPTR = JNZBEG, JNZEND
                        NEWI = LINDX(JPTR)
  200                   CONTINUE
                            I     = NEXTI
                            NEXTI = RCHLNK(I)
                            IF  ( NEWI .GT. NEXTI )  GO TO 200
                        IF  ( NEWI .LT. NEXTI )  THEN
                            KNZ          = KNZ+1
                            RCHLNK(I)    = NEWI
                            RCHLNK(NEWI) = NEXTI
                            MARKER(NEWI) = KSUP
                            NEXTI        = NEWI
                        END IF
                    END DO
                    JSUP = MRGLNK(JSUP)
                    GO TO 100
                END IF
            END IF
*           ----------------------------------------------------
*           Structure of A(*,FSTCOL) has not been examined yet.  
*           "Sort" its structure into the linked list, inserting
*           only those indices not already in the list.
*           ----------------------------------------------------
            IF  ( KNZ .LT. LENGTH )  THEN
                NODE   = PERM(FSTCOL)
                KNZBEG = XADJ(NODE)
                KNZEND = XADJ(NODE+1) - 1
                DO  KPTR = KNZBEG, KNZEND
                    NEWI = ADJNCY(KPTR)
                    NEWI = INVP(NEWI)
                    IF  ( NEWI .GT. FSTCOL  .AND.
     &                    MARKER(NEWI) .NE. KSUP )  THEN
*                       --------------------------------
*                       Position and insert NEWI in list
*                       and mark it with KCOL.
*                       --------------------------------
                        NEXTI = HEAD
  300                   CONTINUE
                            I     = NEXTI
                            NEXTI = RCHLNK(I)
                            IF  ( NEWI .GT. NEXTI )  GO TO 300
                        KNZ          = KNZ + 1
                        RCHLNK(I)    = NEWI
                        RCHLNK(NEWI) = NEXTI
                        MARKER(NEWI) = KSUP
                    END IF
                END DO
            END IF
*           -----------------------------------------------
*           If KSUP has no children, insert FSTCOL into the
*           linked list.
*           -----------------------------------------------
            IF  ( RCHLNK(HEAD) .NE. FSTCOL )  THEN
                RCHLNK(FSTCOL) = RCHLNK(HEAD)
                RCHLNK(HEAD)   = FSTCOL
                KNZ            = KNZ + 1
            END IF
*
*           --------------------------------------------
*           Copy indices from linked list into LINDX(*).
*           --------------------------------------------
            NZBEG = NZEND + 1
            NZEND = NZEND + KNZ
            IF  ( NZEND+1 .NE. XLINDX(KSUP+1) )  THEN
*               -----------------------------------------------
*               Inconsistency in data structure was discovered.
*               -----------------------------------------------
                IFLAG = 23
                RETURN
            END IF
            I = HEAD
            DO  KPTR = NZBEG, NZEND
                I           = RCHLNK(I)
                LINDX(KPTR) = I
            END DO
*
*           ---------------------------------------------------
*           If KSUP has a parent, insert KSUP into its parent's 
*           "merge" list.
*           ---------------------------------------------------
            IF  ( LENGTH .GT. WIDTH )  THEN
                PCOL         = LINDX ( XLINDX(KSUP) + WIDTH )
                PSUP         = SNODE(PCOL)
                MRGLNK(KSUP) = MRGLNK(PSUP)
                MRGLNK(PSUP) = KSUP
            END IF
*
        END DO
        RETURN
*
*       --------------
*       End of SYMFC2.
*       --------------
      END
*
