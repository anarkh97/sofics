        subroutine br8mas(p,rho,elmass,x,y,z,gamma,grvfor
     .,                   grvflg,totmas)
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
        integer p
C
C.... REAL CONSTANTS
C
        real*8 rho, totmas, coef
C
C.... REAL ARRAYS
C
        real*8 elmass(24,24),z(*),x(*),y(*)
        real*8 gamma(*),grvfor(*)
C
C.... LOGICAL CONSTANT
C
        logical grvflg
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
C ... KHP
C ... Modified 5-6-98
C
        if(v .le. 0.0) then
          write(*,*) 'Negative/zero volume'
C          v = -v
        endif
C
C.... COMPUTE THE TOTAL ELEMENT MASS
C
         totmas = rho*v
C
C.... COMPUTE THE MATRIX COEFFICIENT
C
C       c1 = totmas/8.0
C
        c1 = 0.125*totmas
C
C.... COMPUTE BRICK ELEMENT MASS MATRIX
C
        do 40 k = 1, 24
          do 50 l = 1, 24
            elmass(k,l) = 0.0d0
50        continue
40      continue

        do 60 k=1, 24
          elmass(k,k) = c1
60      continue
C
C..... COMPUTE THE BODY FORCE DUE TO GRAVITY IF NEEDED
C
        if (grvflg) then
          coef = 8.0d0*c1
          grvfor(1) = coef*gamma(1)
          grvfor(2) = coef*gamma(2)
          grvfor(3) = coef*gamma(3)
        endif
C
        return
        end
