***************************************************************
* THIS SUBROUTINE COMPUTE THE ELEMENT MASS MATRIX FOR THE 20- *
* NODE BRICK IN A LUMPPED FORM.                               *
*                                                             *
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
        subroutine br20mas(p,rho,elmass,x,y,z,gamma,grvfor,
     *                     grvflg,totmas)
C
C..... DECLARE GLOBAL VARIABLES
C
C.... INTEGER CONSTANTS
C        implicit none
C
        integer p
C
C.... REAL CONSTANTS
C
        real*8 rho, totmas, coef
C
C.... REAL ARRAYS
C
        real*8 elmass(60,60),z(*),x(*),y(*)
        real*8 gamma(*),grvfor(*)
C
C.... LOGICAL CONSTANT
C
        logical grvflg
C
C.... DECLARE LOCAL VARIABLES FOR Q4DMAS
C
        integer k,l,jj
        real*8 v,c1,xi,eta,mu,weight,det
        real*8 q(20),qx(20),qy(20),qz(20)
C
C.... DETERMINE THE VOLUME OF THE BRICK 
C
        v = 0.0d0
        do 10 k=1, p
          do 20 l=1, p
            do 30 jj = 1, p
              call hxgaus20(p,k,p,l,p,jj,xi,eta,mu,weight)
              call h20shpe(xi,eta,mu,x,y,z,q,qx,qy,qz,det)
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
        endif
C
C.... COMPUTE THE TOTAL ELEMENT MASS
C
         totmas = rho*v
C
C.... COMPUTE THE MATRIX COEFFICIENT
C
cPJSA        c1 = 0.125*totmas
        c1 = 0.05*totmas
C
C.... COMPUTE BRICK ELEMENT MASS MATRIX
C
        do 40 k = 1, 60
          do 50 l = 1, 60
            elmass(k,l) = 0.0d0
50        continue
40      continue

        do 60 k=1, 60
          elmass(k,k) = c1
60      continue
C
C..... COMPUTE THE BODY FORCE DUE TO GRAVITY IF NEEDED
C
        if (grvflg) then
cPJSA          coef = 8.0d0*c1
          coef = 20.0d0*c1
          grvfor(1) = coef*gamma(1)
          grvfor(2) = coef*gamma(2)
          grvfor(3) = coef*gamma(3)
        endif
C
        return
        end
