        subroutine e3dmas(rho,elmass,a,x,y,z,gamma,grvfor
     .,                                     grvflg,totmas,masflg)
*-----------------------------------------------------------------*
*THIS SUBROUTINE WILL FORM THE ELEMENT MASS MATRIX                *
*FOR THE 3-D EULER-BERNOULLI BEAM IN LUMPED FORM.                 *
*******************************************************************
*       AUTHOR:  P.R. STERN                                       *
*       DATE:    JANUARY  1991                                    *
*       VERSION: 1.0                                              *
*                                                                 *
*******************************************************************
*
*		VARIABLES
*
*	    A = CROSS SECTIONAL AREA OF THE ELEMENT
*         RHO = ELEMENT MASS DENSITY 
*      ELMASS = STRUCTURAL ELEMENT MASS MATRIX
*      MXNSEQ = MAXIMUM DOF PER ELEMENT
*	    X = LOCAL X COORDINATES
*	    Y = LOCAL Y COORDINATES
*	    Z = LOCAL Z COORDINATES
*       GAMMA = ACCELERATION DUE TO GRAVITY VECTOR
*      GRVFOR = FORCE DUE TO GRAVITY VECTOR
*      GRVFLG = LOGICAL FLAG FOR COMPUTING GRAVITY FORCE
*       MASFLG = FLAG FOR SUBDOMAIN MASS COMPUTATION
*       TOTMAS = ACCUMULATED SUBDOMAIN MASS
*
*
******************************************************************
*
*	CALLED BY : MASS6 
*
*----------------------------------------------------------------*
C
C.... DECLARE THE GLOBAL VARIABLES 
C
C.... REAL CONSTANTS
C
        real*8 a,rho,totmas
C
C.... REAL ARRAYS
C
        real*8 elmass(12,12),x(*),y(*),z(*)
        real*8 gamma(*),grvfor(*)
C
C.... LOGICAL CONSTANT
C
        logical grvflg,masflg
C
C.... DECLARE LOCAL VARIABLES
C
        integer i,j 
        real*8 dx,dy,dz,coeff1,coeff2,c1,len
C
C.... ZERO THE LUMPED ELMENT MASS MATRIX
C
        do 5 i=1,12
          do 15 j=1,12
            elmass(i,j) = 0.0
15        continue
5      continue
C
C.... DETERMINE DISTANCE BETWEEN NODES
C
        dx = x(2) - x(1)
        dy = y(2) - y(1)
        dz = z(2) - z(1)
C
        len = dsqrt((dx*dx)+(dy*dy)+(dz*dz))
C
        coeff1 = A*RHO*len/2.0
        coeff2 = A*RHO*(len**3)/24.0
C
C.... CALCULATE THE ELEMENT LUMPED MASS MATRIX
C
        elmass(1,1) = coeff1
        elmass(1,2) = 0.00 
        elmass(1,3) = 0.00 
        elmass(1,4) = 0.00 
        elmass(1,5) = 0.00 
        elmass(1,6) = 0.00 
C
        elmass(2,2) = coeff1
        elmass(2,3) = 0.00
        elmass(2,4) = 0.00
        elmass(2,5) = 0.00
        elmass(2,6) = 0.00
C
        elmass(3,3) = coeff1
        elmass(3,4) = 0.00
        elmass(3,5) = 0.00
        elmass(3,6) = 0.00
C
        elmass(4,4) = coeff2
        elmass(4,5) = 0.00
        elmass(4,6) = 0.00
C
        elmass(5,5) = coeff2
        elmass(5,6) = 0.00
C
        elmass(6,6) = coeff2
C
C.... USE SYMMETRY TO BUILD THE REST OF THE ELEMENT MASS MATRIX 
C
        do 20 i=1,6
          do 10 j=1,i
            elmass(i,j) = elmass(j,i)
10        continue
20      continue
C
        do 40 i=1,6
          do 30 j=1,6
            elmass(i+6,j+6) = elmass(i,j)
            elmass(i+6,j)   = 0.00
            elmass(i,j+6)   = 0.00
30        continue
40       continue
C
C..... COMPUTE THE BODY FORCE DUE TO GRAVITY IF NEEDED
C
        if (grvflg) then
          c1 = 2.0d0*coeff1
          grvfor(1) = c1*gamma(1)
          grvfor(2) = c1*gamma(2)
          grvfor(3) = c1*gamma(3)
        endif
C
C.... COMPUTE THE SUBDOMAIN MASS
C
        if (masflg) then
          totmas = totmas + 2.0d0*coeff1
        endif
C
        return
        end
        
