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
*--- SPARSPAK-A (ANSI FORTRAN) RELEASE III --- NAME = MMDNUM
*  (C)  UNIVERSITY OF WATERLOO   JANUARY 1984
************************************************************************
************************************************************************
*****   MMDNUM ... MULTI MINIMUM DEGREE NUMBERING                  *****
************************************************************************
************************************************************************
*
*   PURPOSE:
*       THIS SUBROUTINE PERFORMS THE FINAL STEP IN PRODUCING THE
*       PERMUTATION AND INVERSE PERMUTATION VECTORS IN THE MULTIPLE
*       ELIMINATION VERSION OF THE MINIMUM DEGREE ORDERING ALGORITHM.
*
*   INPUT PARAMETERS:
*       NEQNS           -   NUMBER OF EQUATIONS.
*       QSIZE           -   SIZE OF SUPERNODES AT ELIMINATION.
*
*   UPDATED PARAMETERS:
*       INVP            -   INVERSE PERMUTATION VECTOR.  ON INPUT,
*                           IF QSIZE(NODE)=0, THEN NODE HAS BEEN
*                           MERGED INTO THE NODE -INVP(NODE);
*                           OTHERWISE, -INVP(NODE) IS ITS INVERSE
*                           LABELLING.
*
*   OUTPUT PARAMETERS:
*       PERM            -   THE PERMUTATION VECTOR.
*
************************************************************************
*
      SUBROUTINE  MMDNUM  ( NEQNS , PERM  , INVP  , QSIZE           )
*
************************************************************************
*
        INTEGER             INVP(*)  , PERM(*)  , QSIZE(*)
        INTEGER             FATHER, NEQNS , NEXTF , NODE  , NQSIZE,
     &                      NUM   , ROOT
*
************************************************************************
*
        DO  100  NODE = 1, NEQNS
            NQSIZE = QSIZE(NODE)
            IF  ( NQSIZE .LE. 0 )  PERM(NODE) = INVP(NODE)
            IF  ( NQSIZE .GT. 0 )  PERM(NODE) = - INVP(NODE)
  100   CONTINUE
*       ------------------------------------------------------
*       FOR EACH NODE WHICH HAS BEEN MERGED, DO THE FOLLOWING.
*       ------------------------------------------------------
        DO  500  NODE = 1, NEQNS
            IF  ( PERM(NODE) .GT. 0 )  GO TO 500
*               -----------------------------------------
*               TRACE THE MERGED TREE UNTIL ONE WHICH HAS
*               NOT BEEN MERGED, CALL IT ROOT.
*               -----------------------------------------
                FATHER = NODE
  200           CONTINUE
                    IF  ( PERM(FATHER) .GT. 0 )  GO TO 300
                        FATHER = - PERM(FATHER)
                        GO TO 200
  300           CONTINUE
*               -----------------------
*               NUMBER NODE AFTER ROOT.
*               -----------------------
                ROOT = FATHER
                NUM = PERM(ROOT) + 1
                INVP(NODE) = - NUM
                PERM(ROOT) = NUM
*               ------------------------
*               SHORTEN THE MERGED TREE.
*               ------------------------
                FATHER = NODE
  400           CONTINUE
                    NEXTF = - PERM(FATHER)
                    IF  ( NEXTF .LE. 0 )  GO TO 500
                        PERM(FATHER) = - ROOT
                        FATHER = NEXTF
                        GO TO 400
  500   CONTINUE
*       ----------------------
*       READY TO COMPUTE PERM.
*       ----------------------
        DO  600  NODE = 1, NEQNS
            NUM = - INVP(NODE)
            INVP(NODE) = NUM
            PERM(NUM) = NODE
  600   CONTINUE
        RETURN
*
      END
