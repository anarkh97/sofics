        subroutine b3dstf(e,elstif,mxnseq,a,x,y,z,fevol)
*-----------------------------------------------------------------*
*THIS SUBROUTINE WILL FORM THE ELEMENT STIFFNESS MATRIX           *
*FOR THE 3-D TRUSS ELEMENT.                                       *
*******************************************************************
*       AUTHOR:  P.R. STERN                                       *
*       DATE:    JUNE 1989                                        *
*       VERSION: 1.0                                              *
*                                                                 *
*******************************************************************
*
*		VARIABLES
*
*	    A = CROSS SECTIONAL AREA OF THE ELEMENT
*           E = YOUNG'S MODULUS
*      KELSTF = STRUCTURAL ELEMENT STIFFNESS MATRIX
*      MXNSEQ = MAXIMUM DOF PER ELEMENT
*	    X = LOCAL X COORDINATES
*	    Y = LOCAL Y COORDINATES
*	    Z = LOCAL Z COORDINATES
*       FEVOL = VOLUME OF THE ELEMENT
*
******************************************************************
*
*	CALLED BY : STIFF
*
*       SUBROUTINES CALLED : ZRMTX
*
*----------------------------------------------------------------*
C
C.... DECLARE THE GLOBAL VARIABLES 
C
C.... INTEGER CONSTANTS
C
        integer mxnseq
C
C.... REAL CONSTANTS
C
        real*8 a,e,fevol
C
C.... REAL ARRAYS
C
        real*8 elstif(mxnseq,*),x(*),y(*),z(*)
C
C.... DECLARE LOCAL VARIABLES
C
        integer i,j 
        real*8 c11,c21,c31,dx,dy,dz,coeff,len
C
C.... ZERO THE ELMENT STIFFNESS MATRIX
C
CCCCC   call zrmtx(elstif,mxnseq,6,6)
        do 11 i=1,6
          do 21 j=1,6
            elstif(i,j) = 0.
 21       continue
 11     continue
C
C
C.... DETERMINE DISTANCE BETWEEN NODES
C
C       do 1 i =1 , 2
C          write(6,*)'X(',i,') =',x(i),'Y(',i,') =',y(i),
C     &              'Z(',i,') =',z(i)
C1	continue
        dx = x(2) - x(1)
        dy = y(2) - y(1)
        dz = z(2) - z(1)
C
        len = sqrt((dx*dx)+(dy*dy)+(dz*dz))
        c11 = dx/len
        c21 = dy/len
        c31 = dz/len
C
        coeff = E*A/len
C
C.... CALCULATE THE ELEMENT STIFFNESS MATRIX
C
        elstif(1,1) = c11*c11*coeff
        elstif(1,2) = c21*c11*coeff
        elstif(1,3) = c31*c11*coeff
C
        elstif(2,2) = c21*c21*coeff
        elstif(2,3) = c31*c21*coeff
C
        elstif(3,3) = c31*c31*coeff
C
C.... USE SYMMETRY TO BUILD THE REST OF THE ELEMENT STIFFNESS
C
        do 20 i=1,3
          do 10 j=1,i
            elstif(i,j) = elstif(j,i)
10        continue
20      continue
C
        do 40 i=1,3
          do 30 j=1,3
            elstif(i+3,j+3) = elstif(i,j)
            elstif(i+3,j) = -elstif(i,j)
            elstif(i,j+3) = -elstif(i,j)
30        continue
40      continue
C	do 41 i = 1, 6
C	  write(6,*)(elstif(i,j),j=1,6)
C41     continue
C	write(6,*)
C
C.... DETERMINE THE ELEMET'S VOLUME
C
        fevol = A*len

        return
        end
        
