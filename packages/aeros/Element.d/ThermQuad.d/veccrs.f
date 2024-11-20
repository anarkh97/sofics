        subroutine veccrs(a,b,c)

*********************************************************************
* THIS SUBROUTINE WILL TAKE THE VECTOR PRODUCT AxB = C.             *
*********************************************************************
*
*	AUTHOR : P.R. STERN
*	DATE : JULY 1990
*	VERSION : 1.00
*
**********************************************************************
*
*			GLOBAL VARIABLES
*
*	     A = 1ST VECTOR
*	     B = 2ND VECTOR
*	     C = RESULTANT VECTOR PRODUCT
*
**********************************************************************
C
C.... DECLARE GLOBAL VARIABLES
C
C.... REAL ARRAYS 
C
        real*8 a(3),b(3),c(3)
C
C.... LOCAL VARIABLES
C
        real*8 sr
C
C.... COMPUTE THE VECTOR PRODUCT
C
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
C
C.... COMPUTE NORM
C
        sr = dsqrt( c(1)*c(1) + c(2)*c(2) + c(3)*c(3) )
C
C.... NORMALIZE
C
        c(1) = c(1)/sr
        c(2) = c(2)/sr
        c(3) = c(3)/sr

        return
        end
