        subroutine q4maslumpheat(x,y,p,elmass,m)
***************************************************************
* THIS SUBROUTINE COMPUTE THE ELEMENT MASS MATRIX FOR THE 4-  *
* NODE QUADRILATERAL IN A LUMPED FORM.                       *
*                                                             *
***************************************************************
*
*		AUTHOR : Hai
*	        DATE : JAN 1999
*
***************************************************************
*
*		VARIABLES
*
*	     P = PxP GAUSS QUADRATURE RULE
*	ELMASS = THE ELEMENT MASS MATRIX
*	     X = X COORDINATE ARRAY
*	     Y = Y COORDINATE ARRAY
*            M = NUMBER OF DEGREES OF FREEDOM
*
***************************************************************
*
*	SUBROUTINES CALLED : 	QGAUSS, Q4SHPE
*
*	CALLED BY : ThermQuadGal
*
***************************************************************
C
C..... DECLARE GLOBAL VARIABLES
C
C.... INTEGER CONSTANTS
C
        integer p, m
C
C.... REAL ARRAYS
C
        real*8 elmass(m,*),x(*),y(*)
C
C.... DECLARE LOCAL VARIABLES FOR Q4DMAS
C
        integer k,l
        real*8 a,c1,xi,eta,weight,det
        real*8 q(4),qx(4),qy(4)

C        WRITE(6,*) 'Using LUMPED mass'
C
C.... DETERMINE THE AREA OF THE QUADRILATERAL
C
        a = 0.0d0

        do 10 k=1, p
          do 20 l=1, p
            call qgauss(p,k,p,l,xi,eta,weight)
            call q4shpe(xi,eta,x,y,q,qx,qy,det)
            a = a + det
20	  continue
10	continue
C
C.... COMPUTE THE MATRIX COEFFICIENT
C
        c1 = a/4.0
C
C.... COMPUTE THE ELEMENT MASS MATRIX
C
        do 30 k=1, 4
          do 40 l=1, 4
            elmass(k,l) = 0.0d0
40        continue
30      continue

        do 50 k=1, 4
          elmass(k,k) = c1
50      continue
C
        return
        end
