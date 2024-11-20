        subroutine mstf11(elstif,sc1,sc2,sc3,u1,u2,u3,
     &                    v1,v2,v3,w1,w2,w3)
*****************************************************************
* THIS SUBROUTINE WILL BUILD THE ELEMENT STIFFNESS MATRIX FOR   *
* TORSIONAL SPRING ELEMENT.                                     *
*                                                               *
*****************************************************************
*                                                               *
*	AUTHOR  :  P.R. STERN
*	 DATE   :  AUGUST 1991
*	VERSION :  3.00 (Non-Homogeneous Mechanics Revision)
*
*****************************************************************
*
*	          GLOBAL VARIABLES
*
*
*	ELSTIF = ELEMENT STIFFNESS MATRIX
*	ELTDOF = LEADING DIMENSION OF ELSTIF
*          SC1 = SPRING CONSTANT FOR ROTATION ABOUT THE X
*          SC2 = SPRING CONSTANT FOR ROTATION ABOUT THE Y
*          SC3 = SPRING CONSTANT FOR ROTATION ABOUT THE Z 
*     X1,Y1,Z1 = VECTOR FOR 1 DIRECTION
*     X2,Y2,Z2 = VECTOR FOR 2 DIRECTION
*     X3,Y3,Z3 = VECTOR FOR 3 DIRECTION
*
******************************************************************
C
C.... DECLARE GLOBAL VARIABLES
C
C.... REAL CONSTANTS
C
        real*8 u1,u2,u3,v1,v2,v3,w1,w2,w3
        real*8 sc1,sc2,sc3
C
C.... REAL ARRAYS
C
        real*8 elstif(3,3)
C
C.... LOCAL VARIABLES
C
        double precision u(3,3),v(3,3),w(3,3)
C
        call outprd(u1,u2,u3,u,sc1)
        call outprd(v1,v2,v3,v,sc2)
        call outprd(w1,w2,w3,w,sc3)
C
        elstif(1,1) = u(1,1) + v(1,1) + w(1,1)
        elstif(1,2) = u(1,2) + v(1,2) + w(1,2)
        elstif(1,3) = u(1,3) + v(1,3) + w(1,3)
        elstif(2,1) = u(2,1) + v(2,1) + w(2,1)
        elstif(2,2) = u(2,2) + v(2,2) + w(2,2)
        elstif(2,3) = u(2,3) + v(2,3) + w(2,3)
        elstif(3,1) = u(3,1) + v(3,1) + w(3,1)
        elstif(3,2) = u(3,2) + v(3,2) + w(3,2)
        elstif(3,3) = u(3,3) + v(3,3) + w(3,3)
C
        return
        end
