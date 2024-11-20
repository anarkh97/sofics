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
*--- SPARSPAK-A (ANSI FORTRAN) RELEASE III --- NAME = MMDELM
*  (C)  UNIVERSITY OF WATERLOO   JANUARY 1984
************************************************************************
************************************************************************
*****   MMDELM ... MULTIPLE MINIMUM DEGREE ELIMINATION             *****
************************************************************************
************************************************************************
*
*   PURPOSE:
*       THIS SUBROUTINE ELIMINATES THE NODE MDNODE OF MINIMUM DEGREE
*       FROM THE ADJACENCY STRUCTURE, WHICH IS STORED IN THE QUOTIENT
*       GRAPH FORMAT.  IT ALSO TRANSFORMS THE QUOTIENT GRAPH
*       REPRESENTATION OF THE ELIMINATION GRAPH.
*
*   INPUT PARAMETERS:
*       MDNODE          -   NODE OF MINIMUM DEGREE.
*       MAXINT          -   ESTIMATE OF MAXIMUM REPRESENTABLE INTEGER.
*       TAG             -   TAG VALUE.
*
*   UPDATED PARAMETERS:
*       (XADJ,ADJNCY)   -   UPDATED ADJACENCY STRUCTURE.
*       (DHEAD,DFORW,   -   DEGREE DOUBLY LINKED STRUCTURE.
*           DBAKW)
*       QSIZE           -   SIZE OF SUPERNODE.
*       MARKER          -   MARKER VECTOR.
*       LLIST           -   TEMPORARY LINKED LIST OF ELIMINATED NABORS.
*
************************************************************************
*
      SUBROUTINE  MMDELM  ( MDNODE, XADJ  , ADJNCY, DHEAD , DFORW ,
     &                      DBAKW , QSIZE , LLIST , MARKER, MAXINT,
     &                      TAG                                     )
*
************************************************************************
*
        INTEGER             ADJNCY(*), DBAKW(*) , DFORW(*) , DHEAD(*) ,
     &                      LLIST(*) , MARKER(*), QSIZE(*)
        INTEGER             XADJ(*)
        INTEGER             ELMNT , I     , ISTOP , ISTRT , J     ,
     &                      JSTOP , JSTRT , LINK  , MAXINT, MDNODE,
     &                      NABOR , NODE  , NPV   , NQNBRS, NXNODE,
     &                      PVNODE, RLMT  , RLOC  , RNODE , TAG   ,
     &                      XQNBR
*
************************************************************************
*
*       -----------------------------------------------
*       FIND REACHABLE SET AND PLACE IN DATA STRUCTURE.
*       -----------------------------------------------
        MARKER(MDNODE) = TAG
        ISTRT = XADJ(MDNODE)
        ISTOP = XADJ(MDNODE+1) - 1
*       -------------------------------------------------------
*       ELMNT POINTS TO THE BEGINNING OF THE LIST OF ELIMINATED
*       NABORS OF MDNODE, AND RLOC GIVES THE STORAGE LOCATION
*       FOR THE NEXT REACHABLE NODE.
*       -------------------------------------------------------
        ELMNT = 0
        RLOC = ISTRT
        RLMT = ISTOP
        DO  200  I = ISTRT, ISTOP
            NABOR = ADJNCY(I)
            IF  ( NABOR .EQ. 0 )  GO TO 300
                IF  ( MARKER(NABOR) .GE. TAG )  GO TO 200
                    MARKER(NABOR) = TAG
                    IF  ( DFORW(NABOR) .LT. 0 )  GO TO 100
                        ADJNCY(RLOC) = NABOR
                        RLOC = RLOC + 1
                        GO TO 200
  100               CONTINUE
                    LLIST(NABOR) = ELMNT
                    ELMNT = NABOR
  200   CONTINUE
  300   CONTINUE
*           -----------------------------------------------------
*           MERGE WITH REACHABLE NODES FROM GENERALIZED ELEMENTS.
*           -----------------------------------------------------
            IF  ( ELMNT .LE. 0 )  GO TO 1000
                ADJNCY(RLMT) = - ELMNT
                LINK = ELMNT
  400           CONTINUE
                    JSTRT = XADJ(LINK)
                    JSTOP = XADJ(LINK+1) - 1
                    DO  800  J = JSTRT, JSTOP
                        NODE = ADJNCY(J)
                        LINK = - NODE
                        IF  ( NODE )  400, 900, 500
  500                   CONTINUE
                        IF  ( MARKER(NODE) .GE. TAG  .OR.
     &                        DFORW(NODE) .LT. 0 )  GO TO 800
                            MARKER(NODE) = TAG
*                           ---------------------------------
*                           USE STORAGE FROM ELIMINATED NODES
*                           IF NECESSARY.
*                           ---------------------------------
  600                       CONTINUE
                                IF  ( RLOC .LT. RLMT )  GO TO 700
                                    LINK = - ADJNCY(RLMT)
                                    RLOC = XADJ(LINK)
                                    RLMT = XADJ(LINK+1) - 1
                                    GO TO 600
  700                       CONTINUE
                            ADJNCY(RLOC) = NODE
                            RLOC = RLOC + 1
  800               CONTINUE
  900           CONTINUE
                ELMNT = LLIST(ELMNT)
                GO TO 300
 1000   CONTINUE
        IF  ( RLOC .LE. RLMT )  ADJNCY(RLOC) = 0
*       --------------------------------------------------------
*       FOR EACH NODE IN THE REACHABLE SET, DO THE FOLLOWING ...
*       --------------------------------------------------------
        LINK = MDNODE
 1100   CONTINUE
            ISTRT = XADJ(LINK)
            ISTOP = XADJ(LINK+1) - 1
            DO  1700  I = ISTRT, ISTOP
                RNODE = ADJNCY(I)
                LINK = - RNODE
                IF  ( RNODE )  1100, 1800, 1200
 1200           CONTINUE
*               --------------------------------------------
*               IF RNODE IS IN THE DEGREE LIST STRUCTURE ...
*               --------------------------------------------
                PVNODE = DBAKW(RNODE)
                IF  ( PVNODE .EQ. 0  .OR.
     &                PVNODE .EQ. (-MAXINT) )  GO TO 1300
*                   -------------------------------------
*                   THEN REMOVE RNODE FROM THE STRUCTURE.
*                   -------------------------------------
                    NXNODE = DFORW(RNODE)
                    IF  ( NXNODE .GT. 0 )  DBAKW(NXNODE) = PVNODE
                    IF  ( PVNODE .GT. 0 )  DFORW(PVNODE) = NXNODE
                    NPV = - PVNODE
                    IF  ( PVNODE .LT. 0 )  DHEAD(NPV) = NXNODE
 1300           CONTINUE
*               ----------------------------------------
*               PURGE INACTIVE QUOTIENT NABORS OF RNODE.
*               ----------------------------------------
                JSTRT = XADJ(RNODE)
                JSTOP = XADJ(RNODE+1) - 1
                XQNBR = JSTRT
                DO  1400  J = JSTRT, JSTOP
                    NABOR = ADJNCY(J)
                    IF  ( NABOR .EQ. 0 )  GO TO 1500
                        IF  ( MARKER(NABOR) .GE. TAG )  GO TO 1400
                            ADJNCY(XQNBR) = NABOR
                            XQNBR = XQNBR + 1
 1400           CONTINUE
 1500           CONTINUE
*               ----------------------------------------
*               IF NO ACTIVE NABOR AFTER THE PURGING ...
*               ----------------------------------------
                NQNBRS = XQNBR - JSTRT
                IF  ( NQNBRS .GT. 0 )  GO TO 1600
*                   -----------------------------
*                   THEN MERGE RNODE WITH MDNODE.
*                   -----------------------------
                    QSIZE(MDNODE) = QSIZE(MDNODE) + QSIZE(RNODE)
                    QSIZE(RNODE) = 0
                    MARKER(RNODE) = MAXINT
                    DFORW(RNODE) = - MDNODE
                    DBAKW(RNODE) = - MAXINT
                    GO TO 1700
 1600           CONTINUE
*               --------------------------------------
*               ELSE FLAG RNODE FOR DEGREE UPDATE, AND
*               ADD MDNODE AS A NABOR OF RNODE.
*               --------------------------------------
                DFORW(RNODE) = NQNBRS + 1
                DBAKW(RNODE) = 0
                ADJNCY(XQNBR) = MDNODE
                XQNBR = XQNBR + 1
                IF  ( XQNBR .LE. JSTOP )  ADJNCY(XQNBR) = 0
*
 1700       CONTINUE
 1800   CONTINUE
        RETURN
*
      END
