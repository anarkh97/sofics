        subroutine mstf21(elstif,sc1,sc2,sc3,u1,u2,u3,
     &                    v1,v2,v3,w1,w2,w3)
*****************************************************************
* THIS SUBROUTINE WILL BUILD THE ELEMENT STIFFNESS MATRIX FOR   *
* TRANSLATIONAL SPRING-LINK ELEMENT                             * 
*                                                               * 
*****************************************************************
*
*	AUTHOR  :  P.R. STERN
*	 DATE   :  FEBURARY 1995
*	VERSION :  3.00 (Non-Homogeneous Mechanics Revision)
*
*****************************************************************
*
*	          GLOBAL VARIABLES
*
*	ELSTIF = ELEMENT STIFFNESS MATRIX
*	ELTDOF = LEADING DIMENSION OF ELSTIF
*          SC1 = SPRING CONSTANT FOR TRANSLATION ABOUT THE X
*          SC2 = SPRING CONSTANT FOR TRANSLATION ABOUT THE Y
*          SC3 = SPRING CONSTANT FOR TRANSLATION ABOUT THE Z 
*     U1,U1,U1 = VECTOR FOR 1 DIRECTION
*     V2,V2,V2 = VECTOR FOR 2 DIRECTION
*     W3,W3,W3 = VECTOR FOR 3 DIRECTION
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
        real*8 elstif(6,6)
C
C.... LOCAL VARIABLES
C
        integer i,j
        real*8 c(3,3)
C
        c(1,1) = (sc1*u1*u1)+(sc2*v1*v1)+(sc3*w1*w1)
        c(1,2) = (sc1*u1*u2)+(sc2*v1*v2)+(sc3*w1*w2)
        c(1,3) = (sc1*u1*u3)+(sc2*v1*v3)+(sc3*w1*w3)
        c(2,2) = (sc1*u2*u2)+(sc2*v2*v2)+(sc3*w2*w2)
        c(2,3) = (sc1*u2*u3)+(sc2*v2*v3)+(sc3*w2*w3)
        c(3,3) = (sc1*u3*u3)+(sc2*v3*v3)+(sc3*w3*w3)
        c(2,1) = c(1,2)
        c(3,1) = c(1,3)
        c(3,2) = c(2,3)
C
        do 100 i =1, 3
          do 200 j = 1, 3
            elstif(i  ,j  )   =  c(i,j)
            elstif(i+3,j+3)   =  c(i,j)
            elstif(i  ,j+3)   = -c(i,j)
            elstif(i+3,j  )   = -c(i,j)
200       continue
100     continue
C
        return
        end
