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
*--- SPARSPAK-A (ANSI FORTRAN) RELEASE III --- NAME = MMDINT
*  (C)  UNIVERSITY OF WATERLOO   JANUARY 1984
************************************************************************
************************************************************************
*****   MMDINT ... MULT MINIMUM DEGREE INITIALIZATION              *****
************************************************************************
************************************************************************
*
*   PURPOSE:
*       THIS SUBROUTINE PERFORMS INITIALIZATION FOR THE MULTIPLE
*       ELIMINATION VERSION OF THE MINIMUM DEGREE ALGORITHM.
*
*   INPUT PARAMETERS:
*       NEQNS           -   NUMBER OF EQUATIONS.
*       (XADJ,ADJNCY)   -   ADJACENCY STRUCTURE.
*
*   OUTPUT PARAMETERS:
*       (DHEAD,DFORW,   -   DEGREE DOUBLY LINKED STRUCTURE.
*           DBAKW)
*       QSIZE           -   SIZE OF SUPERNODE (INITIALIZED TO ONE).
*       LLIST           -   LINKED LIST.
*       MARKER          -   MARKER VECTOR.
*
************************************************************************
*
      SUBROUTINE  MMDINT  ( NEQNS , XADJ  , ADJNCY, DHEAD , DFORW ,
     &                      DBAKW , QSIZE , LLIST , MARKER          )
*
************************************************************************
*
        INTEGER             ADJNCY(*), DBAKW(*) , DFORW(*) , DHEAD(*) ,
     &                      LLIST(*) , MARKER(*), QSIZE(*)
        INTEGER             XADJ(*)
        INTEGER             FNODE , NDEG  , NEQNS , NODE
*
************************************************************************
*
        DO  100  NODE = 1, NEQNS
            DHEAD(NODE) = 0
            QSIZE(NODE) = 1
            MARKER(NODE) = 0
            LLIST(NODE) = 0
  100   CONTINUE
*       ------------------------------------------
*       INITIALIZE THE DEGREE DOUBLY LINKED LISTS.
*       ------------------------------------------
        DO  200  NODE = 1, NEQNS
            NDEG = XADJ(NODE+1) - XADJ(NODE) + 1
            FNODE = DHEAD(NDEG)
            DFORW(NODE) = FNODE
            DHEAD(NDEG) = NODE
            IF  ( FNODE .GT. 0 )  DBAKW(FNODE) = NODE
            DBAKW(NODE) = - NDEG
  200   CONTINUE
        RETURN
*
      END
