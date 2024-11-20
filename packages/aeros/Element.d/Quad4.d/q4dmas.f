        subroutine q4dmas(p,rho,elmass,mxnseq,h,x,y,gamma,grvfor
     .,                                     grvflg,totmas,masflg)
***************************************************************
* THIS SUBROUTINE COMPUTE THE ELEMENT MASS MATRIX FOR THE 4-  *
* NODE QUADRILATERAL IN A LUMPPED FORM.                       *
*                                                             *
***************************************************************
*
*		AUTHOR : P.R. STERN
*	        DATE : APRIL 1990
*	        VER : 1.00
*
***************************************************************
*
*		VARIABLES
*
*	     P = PxP GAUSS QUADRATURE RULE
*	   RHO = ELEMENT DENSITY
*	ELMASS = THE ELEMENT MASS MATRIX
*	MXNSEQ = LEADING DIMENSION OF ELMASS
*	     H = ARRAY CONTAINING NODAL THICKNESS
*	     X = X COORDINATE ARRAY
*	     Y = Y COORDINATE ARRAY
*        GAMMA = ACCELERATIONS DUE TO GRAVITY VECTOR
*       GRVFOR = BODY FORCES DUE TO GRAVITY
*       GRVFLG = LOGICAL FLAG TO INDICATE GRAVITY FORCES NEEDED
*       MASFLG = FLAG FOR SUBDOMAIN MASS COMPUTATION
*       TOTMAS = ACCUMULATED SUBDOMAIN MASS
*
***************************************************************
*
*	SUBROUTINES CALLED : 	QGAUSS, Q4SHPE
*
*	CALLED BY : MASS2 
*
***************************************************************
C
C..... DECLARE GLOBAL VARIABLES
C
C.... INTEGER CONSTANTS
C
        integer p,mxnseq
C
C.... REAL CONSTANTS
C
        real*8 rho,totmas
C
C.... REAL ARRAYS
C
        real*8 elmass(mxnseq,*),h(*),x(*),y(*)
        real*8 gamma(*),grvfor(*)
C
C.... LOGICAL CONSTANTS
C
        logical grvflg,masflg
C
C.... DECLARE LOCAL VARIABLES FOR Q4DMAS
C
        integer k,l
        real*8 a,c1,xi,eta,weight,det
        real*8 q(4),qx(4),qy(4)
C
C.... DETERMINE THE AREA OF THE QUADRILATERAL
C
        a = 0.0d0
        do 10 k=1, p
          do 20 l=1, p
            call qgauss(p,k,p,l,xi,eta,weight)
            call q4shpe(xi,eta,x,y,q,qx,qy,det)
            a = a + det
20        continue
10      continue
C
C.... COMPUTE THE MATRIX COEFFICIENT
C
        c1 = rho*h(1)*a/4.0
C
C.... COMPUTE THE ELEMENT MASS MATRIX
C
        elmass(1,1) = c1
        elmass(1,2) = 0.0
        elmass(1,3) = 0.0
        elmass(1,4) = 0.0
        elmass(2,2) = c1
        elmass(2,3) = 0.0
        elmass(2,4) = 0.0
        elmass(3,3) = c1
        elmass(3,4) = 0.0
        elmass(4,4) = c1
C
C.... USE SYMMETRY TO BUILD THE REST OF THE MATRIX
C
        do 30 k=1, 4
          do 40 l=1, k
            elmass(k,l) = elmass(l,k)
40        continue
30      continue
C
        do 50 k=1, 4
          do 60 l=1, 4
            elmass(k+4,l+4) = elmass(k,l)
            elmass(k+4,l) = 0.0d0
            elmass(k,l+4) = 0.0d0
60        continue
50      continue

C
C..... COMPUTE THE BODY FORCE DUE TO GRAVITY IF NEEDED
C
        if (grvflg) then
          grvfor(1) = 4.0d0*c1*gamma(1)
          grvfor(2) = 4.0d0*c1*gamma(2)
        endif
C
C.... COMPUTE THE SUBDOMAIN TOTAL MASS
C
        if (masflg) then
          totmas = totmas + 4.0d0*c1
        endif
C
        return
        end
