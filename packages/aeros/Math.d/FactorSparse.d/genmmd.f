************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    Joseph W.H. Liu
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
*--- SPARSPAK-A (ANSI FORTRAN) RELEASE III --- NAME = GENMMD
*  (C)  UNIVERSITY OF WATERLOO   JANUARY 1984
************************************************************************
************************************************************************
*****   GENMMD ... MULTIPLE MINIMUM EXTERNAL DEGREE                *****
************************************************************************
************************************************************************
*
*   PURPOSE:
*       THIS SUBROUTINE IMPLEMENTS THE MINIMUM DEGREE ALGORITHM.  IT
*       MAKES USE OF THE IMPLICIT REPRESENTATION OF ELIMINATION GRAPHS
*       BY QUOTIENT GRAPHS, AND THE NOTION OF INDISTINGUISHABLE NODES.
*       IT ALSO IMPLEMENTS THE MODIFICATIONS BY MULTIPLE ELIMINATION
*       AND MINIMUM EXTERNAL DEGREE.
*
*       --------------------------------------------------------
*       CAUTION - THE ADJACENCY VECTOR ADJNCY WILL BE DESTROYED.
*       --------------------------------------------------------
*
*   INPUT PARAMETERS:
*       NEQNS           -   NUMBER OF EQUATIONS.
*       (XADJ,ADJNCY)   -   THE ADJACENCY STRUCTURE.
*       DELTA           -   TOLERANCE VALUE FOR MULTIPLE ELIMINATION.
*       MAXINT          -   MAXIMUM MACHINE REPRESENTABLE INTEGER
*                           (ANY SMALLER ESTIMATE WILL DO) FOR MARKING
*                           NODES.
*
*   OUTPUT PARAMETERS:
*       PERM            -   THE MINIMUM DEGREE ORDERING.
*       INVP            -   THE INVERSE OF PERM.
*       NOFSUB          -   AN UPPER BOUND ON THE NUMBER OF NONZERO
*                           SUBSCRIPTS FOR THE COMPRESSED STORAGE
*                           SCHEME.
*
*   WORKING PARAMETERS:
*       DHEAD           -   VECTOR FOR HEAD OF DEGREE LISTS.
*       INVP            -   USED TEMPORARILY FOR DEGREE FORWARD LINK.
*       PERM            -   USED TEMPORARILY FOR DEGREE BACKWARD LINK.
*       QSIZE           -   VECTOR FOR SIZE OF SUPERNODES.
*       LLIST           -   VECTOR FOR TEMPORARY LINKED LISTS.
*       MARKER          -   A TEMPORARY MARKER VECTOR.
*
*   PROGRAM SUBROUTINES:
*       MMDELM, MMDINT, MMDNUM, MMDUPD.
*
************************************************************************
*
      SUBROUTINE  GENMMD  ( NEQNS , XADJ  , ADJNCY, INVP  , PERM  ,
     &                      DELTA , DHEAD , QSIZE , LLIST , MARKER,
     &                      MAXINT, NOFSUB                          )
*
************************************************************************
*
        INTEGER             ADJNCY(*), DHEAD(*) , INVP(*)  , LLIST(*) ,
     &                      MARKER(*), PERM(*)  , QSIZE(*)
        INTEGER             XADJ(*)
        INTEGER             DELTA , EHEAD , I     , MAXINT, MDEG  ,
     &                      MDLMT , MDNODE, NEQNS , NEXTMD, NOFSUB,
     &                      NUM, TAG
*
************************************************************************
*
        IF  ( NEQNS .LE. 0 )  RETURN
*
*       ------------------------------------------------
*       INITIALIZATION FOR THE MINIMUM DEGREE ALGORITHM.
*       ------------------------------------------------
        NOFSUB = 0
        CALL  MMDINT ( NEQNS, XADJ, ADJNCY, DHEAD, INVP, PERM,
     &                 QSIZE, LLIST, MARKER )
*
*       ----------------------------------------------
*       NUM COUNTS THE NUMBER OF ORDERED NODES PLUS 1.
*       ----------------------------------------------
        NUM = 1
*
*       -----------------------------
*       ELIMINATE ALL ISOLATED NODES.
*       -----------------------------
        NEXTMD = DHEAD(1)
  100   CONTINUE
            IF  ( NEXTMD .LE. 0 )  GO TO 200
                MDNODE = NEXTMD
                NEXTMD = INVP(MDNODE)
                MARKER(MDNODE) = MAXINT
                INVP(MDNODE) = - NUM
                NUM = NUM + 1
                GO TO 100
*
  200   CONTINUE
*       ----------------------------------------
*       SEARCH FOR NODE OF THE MINIMUM DEGREE.
*       MDEG IS THE CURRENT MINIMUM DEGREE;
*       TAG IS USED TO FACILITATE MARKING NODES.
*       ----------------------------------------
        IF  ( NUM .GT. NEQNS )  GO TO 1000
        TAG = 1
        DHEAD(1) = 0
        MDEG = 2
  300   CONTINUE
            IF  ( DHEAD(MDEG) .GT. 0 )  GO TO 400
                MDEG = MDEG + 1
                GO TO 300
  400       CONTINUE
*           -------------------------------------------------
*           USE VALUE OF DELTA TO SET UP MDLMT, WHICH GOVERNS
*           WHEN A DEGREE UPDATE IS TO BE PERFORMED.
*           -------------------------------------------------
            MDLMT = MDEG + DELTA
            EHEAD = 0
*
  500       CONTINUE
                MDNODE = DHEAD(MDEG)
                IF  ( MDNODE .GT. 0 )  GO TO 600
                    MDEG = MDEG + 1
                    IF  ( MDEG .GT. MDLMT )  GO TO 900
                        GO TO 500
  600           CONTINUE
*               ----------------------------------------
*               REMOVE MDNODE FROM THE DEGREE STRUCTURE.
*               ----------------------------------------
                NEXTMD = INVP(MDNODE)
                DHEAD(MDEG) = NEXTMD
                IF  ( NEXTMD .GT. 0 )  PERM(NEXTMD) = - MDEG
                INVP(MDNODE) = - NUM
                NOFSUB = NOFSUB + MDEG + QSIZE(MDNODE) - 2
                IF  ( NUM+QSIZE(MDNODE) .GT. NEQNS )  GO TO 1000
*               ----------------------------------------------
*               ELIMINATE MDNODE AND PERFORM QUOTIENT GRAPH
*               TRANSFORMATION.  RESET TAG VALUE IF NECESSARY.
*               ----------------------------------------------
                TAG = TAG + 1
                IF  ( TAG .LT. MAXINT )  GO TO 800
                    TAG = 1
                    DO  700  I = 1, NEQNS
                        IF  ( MARKER(I) .LT. MAXINT )  MARKER(I) = 0
  700               CONTINUE
  800           CONTINUE
                CALL  MMDELM ( MDNODE, XADJ, ADJNCY, DHEAD, INVP,
     &                         PERM, QSIZE, LLIST, MARKER, MAXINT,
     &                         TAG )
                NUM = NUM + QSIZE(MDNODE)
                LLIST(MDNODE) = EHEAD
                EHEAD = MDNODE
                IF  ( DELTA .GE. 0 )  GO TO 500
  900       CONTINUE
*           -------------------------------------------
*           UPDATE DEGREES OF THE NODES INVOLVED IN THE
*           MINIMUM DEGREE NODES ELIMINATION.
*           -------------------------------------------
            IF  ( NUM .GT. NEQNS )  GO TO 1000
            CALL  MMDUPD ( EHEAD, NEQNS, XADJ, ADJNCY, DELTA, MDEG,
     &                     DHEAD, INVP, PERM, QSIZE, LLIST, MARKER,
     &                     MAXINT, TAG )
            GO TO 300
*
 1000   CONTINUE
        CALL  MMDNUM ( NEQNS, PERM, INVP, QSIZE )
        RETURN
*
      END
