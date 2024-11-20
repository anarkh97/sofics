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
*--- SPARSPAK-A (ANSI FORTRAN) RELEASE III --- NAME = MMDUPD
*  (C)  UNIVERSITY OF WATERLOO   JANUARY 1984
************************************************************************
************************************************************************
*****   MMDUPD ... MULTIPLE MINIMUM DEGREE UPDATE                  *****
************************************************************************
************************************************************************
*
*   PURPOSE:
*       THIS SUBROUTINE UPDATES THE DEGREES OF NODES AFTER A MULTIPLE
*       ELIMINATION STEP.
*
*   INPUT PARAMETERS:
*       EHEAD           -   THE BEGINNING OF THE LIST OF ELIMINATED
*                           NODES (I.E., NEWLY FORMED ELEMENTS).
*       NEQNS           -   NUMBER OF EQUATIONS.
*       (XADJ,ADJNCY)   -   ADJACENCY STRUCTURE.
*       DELTA           -   TOLERANCE VALUE FOR MULTIPLE ELIMINATION.
*       MAXINT          -   MAXIMUM MACHINE REPRESENTABLE INTEGER.
*
*   UPDATED PARAMETERS:
*       MDEG            -   NEW MINIMUM DEGREE AFTER DEGREE UPDATE.
*       (DHEAD,DFORW    -   DEGREE DOUBLY LINKED STRUCTURE.
*           DBAKW)
*       QSIZE           -   SIZE OF SUPERNODE.
*       LLIST           -   WORKING LINKED LIST.
*       MARKER          -   MARKER VECTOR FOR DEGREE UPDATE.
*       TAG             -   TAG VALUE.
*
************************************************************************
*
      SUBROUTINE  MMDUPD  ( EHEAD , NEQNS , XADJ  , ADJNCY, DELTA ,
     &                      MDEG  , DHEAD , DFORW , DBAKW , QSIZE ,
     &                      LLIST , MARKER, MAXINT, TAG             )
*
************************************************************************
*
        INTEGER             ADJNCY(*), DBAKW(*) , DFORW(*) , DHEAD(*) ,
     &                      LLIST(*) , MARKER(*), QSIZE(*)
        INTEGER             XADJ(*)
        INTEGER             DEG   , DEG0  , DELTA , EHEAD , ELMNT ,
     &                      ENODE , FNODE , I     , IQ2   , ISTOP ,
     &                      ISTRT , J     , JSTOP , JSTRT , LINK  ,
     &                      MAXINT, MDEG  , MDEG0 , MTAG  , NABOR ,
     &                      NEQNS , NODE  , Q2HEAD, QXHEAD, TAG
*
************************************************************************
*
        MDEG0 = MDEG + DELTA
        ELMNT = EHEAD
  100   CONTINUE
*           -------------------------------------------------------
*           FOR EACH OF THE NEWLY FORMED ELEMENT, DO THE FOLLOWING.
*           (RESET TAG VALUE IF NECESSARY.)
*           -------------------------------------------------------
            IF  ( ELMNT .LE. 0 )  RETURN
            MTAG = TAG + MDEG0
            IF  ( MTAG .LT. MAXINT )  GO TO 300
                TAG = 1
                DO  200  I = 1, NEQNS
                    IF  ( MARKER(I) .LT. MAXINT )  MARKER(I) = 0
  200           CONTINUE
                MTAG = TAG + MDEG0
  300       CONTINUE
*           ---------------------------------------------
*           CREATE TWO LINKED LISTS FROM NODES ASSOCIATED
*           WITH ELMNT: ONE WITH TWO NABORS (Q2HEAD) IN
*           ADJACENCY STRUCTURE, AND THE OTHER WITH MORE
*           THAN TWO NABORS (QXHEAD).  ALSO COMPUTE DEG0,
*           NUMBER OF NODES IN THIS ELEMENT.
*           ---------------------------------------------
            Q2HEAD = 0
            QXHEAD = 0
            DEG0 = 0
            LINK = ELMNT
  400       CONTINUE
                ISTRT = XADJ(LINK)
                ISTOP = XADJ(LINK+1) - 1
                DO  700  I = ISTRT, ISTOP
                    ENODE = ADJNCY(I)
                    LINK = - ENODE
                    IF  ( ENODE )  400, 800, 500
*
  500               CONTINUE
                    IF  ( QSIZE(ENODE) .EQ. 0 )  GO TO 700
                        DEG0 = DEG0 + QSIZE(ENODE)
                        MARKER(ENODE) = MTAG
*                       ----------------------------------
*                       IF ENODE REQUIRES A DEGREE UPDATE,
*                       THEN DO THE FOLLOWING.
*                       ----------------------------------
                        IF  ( DBAKW(ENODE) .NE. 0 )  GO TO 700
*                           ---------------------------------------
*                           PLACE EITHER IN QXHEAD OR Q2HEAD LISTS.
*                           ---------------------------------------
                            IF  ( DFORW(ENODE) .EQ. 2 )  GO TO 600
                                LLIST(ENODE) = QXHEAD
                                QXHEAD = ENODE
                                GO TO 700
  600                       CONTINUE
                            LLIST(ENODE) = Q2HEAD
                            Q2HEAD = ENODE
  700           CONTINUE
  800       CONTINUE
*           --------------------------------------------
*           FOR EACH ENODE IN Q2 LIST, DO THE FOLLOWING.
*           --------------------------------------------
            ENODE = Q2HEAD
            IQ2 = 1
  900       CONTINUE
                IF  ( ENODE .LE. 0 )  GO TO 1500
                IF  ( DBAKW(ENODE) .NE. 0 )  GO TO 2200
                    TAG = TAG + 1
                    DEG = DEG0
*                   ------------------------------------------
*                   IDENTIFY THE OTHER ADJACENT ELEMENT NABOR.
*                   ------------------------------------------
                    ISTRT = XADJ(ENODE)
                    NABOR = ADJNCY(ISTRT)
                    IF  ( NABOR .EQ. ELMNT )  NABOR = ADJNCY(ISTRT+1)
*                   ------------------------------------------------
*                   IF NABOR IS UNELIMINATED, INCREASE DEGREE COUNT.
*                   ------------------------------------------------
                    LINK = NABOR
                    IF  ( DFORW(NABOR) .LT. 0 )  GO TO 1000
                        DEG = DEG + QSIZE(NABOR)
                        GO TO 2100
 1000               CONTINUE
*                       --------------------------------------------
*                       OTHERWISE, FOR EACH NODE IN THE 2ND ELEMENT,
*                       DO THE FOLLOWING.
*                       --------------------------------------------
                        ISTRT = XADJ(LINK)
                        ISTOP = XADJ(LINK+1) - 1
                        DO  1400  I = ISTRT, ISTOP
                            NODE = ADJNCY(I)
                            LINK = - NODE
                            IF  ( NODE .EQ. ENODE )  GO TO 1400
                            IF  ( NODE )  1000, 2100, 1100
*
 1100                       CONTINUE
                            IF  ( QSIZE(NODE) .EQ. 0 )  GO TO 1400
                            IF  ( MARKER(NODE) .GE. TAG )  GO TO 1200
*                               -------------------------------------
*                               CASE WHEN NODE IS NOT YET CONSIDERED.
*                               -------------------------------------
                                MARKER(NODE) = TAG
                                DEG = DEG + QSIZE(NODE)
                                GO TO 1400
 1200                       CONTINUE
*                           ----------------------------------------
*                           CASE WHEN NODE IS INDISTINGUISHABLE FROM
*                           ENODE.  MERGE THEM INTO A NEW SUPERNODE.
*                           ----------------------------------------
                            IF  ( DBAKW(NODE) .NE. 0 )  GO TO 1400
                            IF  ( DFORW(NODE) .NE. 2 )  GO TO 1300
                                QSIZE(ENODE) = QSIZE(ENODE) +
     &                                         QSIZE(NODE)
                                QSIZE(NODE) = 0
                                MARKER(NODE) = MAXINT
                                DFORW(NODE) = - ENODE
                                DBAKW(NODE) = - MAXINT
                                GO TO 1400
 1300                       CONTINUE
*                           --------------------------------------
*                           CASE WHEN NODE IS OUTMATCHED BY ENODE.
*                           --------------------------------------
                            IF  ( DBAKW(NODE) .EQ.0 )
     &                            DBAKW(NODE) = - MAXINT
 1400                   CONTINUE
                        GO TO 2100
 1500           CONTINUE
*               ------------------------------------------------
*               FOR EACH ENODE IN THE QX LIST, DO THE FOLLOWING.
*               ------------------------------------------------
                ENODE = QXHEAD
                IQ2 = 0
 1600           CONTINUE
                    IF  ( ENODE .LE. 0 )  GO TO 2300
                    IF  ( DBAKW(ENODE) .NE. 0 )  GO TO 2200
                        TAG = TAG + 1
                        DEG = DEG0
*                       ---------------------------------
*                       FOR EACH UNMARKED NABOR OF ENODE,
*                       DO THE FOLLOWING.
*                       ---------------------------------
                        ISTRT = XADJ(ENODE)
                        ISTOP = XADJ(ENODE+1) - 1
                        DO  2000  I = ISTRT, ISTOP
                            NABOR = ADJNCY(I)
                            IF  ( NABOR .EQ. 0 )  GO TO 2100
                            IF  ( MARKER(NABOR) .GE. TAG )  GO TO 2000
                                MARKER(NABOR) = TAG
                                LINK = NABOR
*                               ------------------------------
*                               IF UNELIMINATED, INCLUDE IT IN
*                               DEG COUNT.
*                               ------------------------------
                                IF  ( DFORW(NABOR) .LT. 0 )  GO TO 1700
                                    DEG = DEG + QSIZE(NABOR)
                                    GO TO 2000
 1700                           CONTINUE
*                                   -------------------------------
*                                   IF ELIMINATED, INCLUDE UNMARKED
*                                   NODES IN THIS ELEMENT INTO THE
*                                   DEGREE COUNT.
*                                   -------------------------------
                                    JSTRT = XADJ(LINK)
                                    JSTOP = XADJ(LINK+1) - 1
                                    DO  1900  J = JSTRT, JSTOP
                                        NODE = ADJNCY(J)
                                        LINK = - NODE
                                        IF  ( NODE )  1700, 2000, 1800
*
 1800                                   CONTINUE
                                        IF  ( MARKER(NODE) .GE. TAG )
     &                                        GO TO 1900
                                            MARKER(NODE) = TAG
                                            DEG = DEG + QSIZE(NODE)
 1900                               CONTINUE
 2000                   CONTINUE
 2100               CONTINUE
*                   -------------------------------------------
*                   UPDATE EXTERNAL DEGREE OF ENODE IN DEGREE
*                   STRUCTURE, AND MDEG (MIN DEG) IF NECESSARY.
*                   -------------------------------------------
                    DEG = DEG - QSIZE(ENODE) + 1
                    FNODE = DHEAD(DEG)
                    DFORW(ENODE) = FNODE
                    DBAKW(ENODE) = - DEG
                    IF  ( FNODE .GT. 0 )  DBAKW(FNODE) = ENODE
                    DHEAD(DEG) = ENODE
                    IF  ( DEG .LT. MDEG )  MDEG = DEG
 2200               CONTINUE
*                   ----------------------------------
*                   GET NEXT ENODE IN CURRENT ELEMENT.
*                   ----------------------------------
                    ENODE = LLIST(ENODE)
                    IF  ( IQ2 .EQ. 1 )  GO TO 900
                        GO TO 1600
 2300       CONTINUE
*           -----------------------------
*           GET NEXT ELEMENT IN THE LIST.
*           -----------------------------
            TAG = MTAG
            ELMNT = LLIST(ELMNT)
            GO TO 100
*
      END
