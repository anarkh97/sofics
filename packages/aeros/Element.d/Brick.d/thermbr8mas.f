        subroutine thermbr8mas(p,rho,elmass,mxnseq,x,y,z,totmas)

***************************************************************
* THIS SUBROUTINE COMPUTE THE ELEMENT MASS MATRIX FOR THE 8-  *
* NODE BRICK IN A LUMPPED FORM.                               *
*                                                             *
***************************************************************
*                                                             *
*		AUTHOR : P.R. STERN                           *
*	        DATE : APRIL 1994                             *
*	        VER : 1.00                                    *
*							      *
*		MODIFIED : FEBRUARY 1997 TO RETURN TOTAL MASS *
***************************************************************
*                                                             *
*		VARIABLES 				      *
*                                                             *
*	     P = PxPxP GAUSS QUADRATURE RULE                  *
*	   RHO = ELEMENT DENSITY                              *
*	ELMASS = THE ELEMENT MASS MATRIX                      *
*	MXNSEQ = LEADING DIMENSION OF ELMASS                  *
*	     X = X COORDINATE ARRAY                           *
*	     Y = Y COORDINATE ARRAY                           *
*	     Z = Z COORDINATE ARRAY                           *
*       GAMMA = ACCELERATION DUE TO GRAVITY VECTOR            *
*      GRVFOR = FORCE DUE TO GRAVITY VECTOR                   *
*      GRVFLG = LOGICAL FLAG FOR COMPUTING GRAVITY FORCE      *
*
***************************************************************
*                                                             *
*	SUBROUTINES CALLED : 	HXGAUS, H8SHPE 		      *
*                                                             *
*	CALLED BY : MASS8 				      *
*                                                             *
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
        real*8 rho, totmas
C
C.... REAL ARRAYS
C
        real*8 elmass(mxnseq,*),z(*),x(*),y(*)
C
C.... DECLARE LOCAL VARIABLES FOR Q4DMAS
C
        integer k,l,jj
        real*8 v,c1,xi,eta,emu,weight,det
        real*8 q(8),qx(8),qy(8),qz(8)
C
C.... DETERMINE THE VOLUME OF THE BRICK 
C
        v = 0.0d0
        do 10 k=1, p
          do 20 l=1, p
            do 30 jj = 1, p
              call hxgaus(p,k,p,l,p,jj,xi,eta,emu,weight)
              call h8shpe(xi,eta,emu,x,y,z,q,qx,qy,qz,det)
              v = v + det
30          continue
20        continue
10      continue

C
C.... COMPUTE THE TOTAL ELEMENT MASS
C
         totmas = rho*v
C
C.... COMPUTE THE MATRIX COEFFICIENT
C
        c1 = totmas/8.0
C
C.... COMPUTE THE ELEMENT MASS MATRIX
C
        do 40 k = 1, 8
          do 50 l = 1, 8
            elmass(k,l) = 0.0d0
50        continue
40      continue

        do 60 k=1, 8
          elmass(k,k) = c1
60      continue
C
        return
        end
