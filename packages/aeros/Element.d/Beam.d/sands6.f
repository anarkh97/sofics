        subroutine sands6(area,e,elm,stress,maxsze,maxgus,maxstr,eframe,
     &                    ix,iy,iz,nu,x,y,z,ug,alpha,tref,Temp)
*********************************************************************
*	THIS SUBROUTINE WILL COMPUTE THE INTERNAL FORCES FOR THE    *
* EULER-BERNOULLI BEAM ELEMENT AND STORE THEM IN THE STRESS ARRAY   *
*                                                                   *
*********************************************************************
*
*	AUTHOR  :  P.R. STERN
*	DATE    :  SEPTEMBER 1991
*	VERSION :  3.0 Multielement Revison
*
*********************************************************************
*
*		DEFINE THE GLOBAL VARIABLES
*
*	  AREA = CROSS-SECTIONAL AREA OF THE BEAM
*	     E = YOUNGS MODULUS FORTHE BEAM
*        ALPHA = DILATATION COEFFICIENT
*	   ELM = CURRENT ELEMENT NUMBER
*	STRESS = INTERNAL FORCE ARRAY
*       MAXSZE = LEADING DIMENSION OF STRESS AND STRAIN ARRAYS
*	MAXGUS = SECOND DIMENSION OF STRESS
*	MAXSTR = THIRD DIMENSION OF STRESS
*	EFRAME = ELEMENT REFERENCE FRAMES
*     IX,IY,IZ = BEAM MOMENTS OF INERTIA
*	    NU = POISSON'S RATIO
*        X,Y,Z = COORDINATES FOR THE BEAM
*           UG = GLOBAL DISPLACEMENT VECTOR FOR ELEMENT #ELM
*         TEMP = NODAL TEMPERATURE
*         TREF = REFERENCE TEMPERATURE
*      TSTRESS = THERMAL STRESS
*
*********************************************************************
*
*		CALLED BY : MDERIV
*
*********************************************************************
C
C.... DECLARE ALL GLOBAL VARIABLES
C
C.... INTEGER CONSTANTS
C
        integer maxsze,maxstr,maxgus,elm
C
C.... REAL CONSTANTS
C
        real*8 area,e,nu,ix,iy,iz,eiy,eiz,length2,length3
        real*8 alpha, tref
C
C.... REAL ARRAYS
C
        real*8 x(2),y(2),z(2),eframe(*),ug(*)
        real*8 stress(maxsze,maxstr,maxgus)
        real*8 temp(2), tl(2)
C
C.... LOCAL VARIABLES
C
        integer ii,jj,kk,ic
C
        real*8 pi,G,J,dx,dy,dz,length
        real*8 stress_0dxy,stress_0dyz,stress_0dxz
        real*8 stress_1dxy,stress_1dyz,stress_1dxz
        real*8 vec(10),ke(12,12),tug(12),res(12)
        real*8 u(3),v(3),w(3)
        real*8 tran(12,12),t33(3,3)
C
C.... INITIALIZE LOCAL CONSTANTS
C
        pi = 4.0d00*atan(1.0)
        G  = e / (2.0d00*(1.0d00+nu))
**      J  = 0.5d00*area*area / pi 
        J  = Ix
C
C.... COMPUTE THE LENGTH OF THE BEAM
C
        dx = x(2) - x(1)
        dy = y(2) - y(1)
        dz = z(2) - z(1)
        length = dsqrt(dx*dx + dy*dy + dz*dz)
C
C.... COMPUTE THE CONSTITUTIVE RELATIONS FOR THE BEAM
C
        eiy     = e*iy
        eiz     = e*iz
        length2 = length*length
        length3 = length2*length
        vec(1)  =  E * area / length
        vec(2)  = 12.0d00 * eiy / length3
        vec(3)  = 12.0d00 * eiz / length3
        vec(4)  =  6.0d00 * eiy / length2      
        vec(5)  =  6.0d00 * eiz / length2 
        vec(6)  =  G * J / length              
        vec(7)  =  4.0d00 * eiy / length      
        vec(8)  =  4.0d00 * eiz / length     
        vec(9)  =  2.0d00 * eiy / length    
        vec(10) =  2.0d00 * eiz / length   

C
C.... INITIALIZE MATRICES AND VECTORS  
C
        do 110 jj = 1, 12
          tug(jj) = 0.0d00
          res(jj) = 0.0d00
          do 110 ii = 1, 12
            ke(ii,jj) = 0.0d00
            tran(ii,jj) = 0.0d00
 110	continue
C
C.... BUILD THE ELEMENT STIFFNESS MATRIX
C
        ke(1,1) = vec(1)
        ke(2,2) = vec(3)
        ke(3,3) = vec(2)
        ke(4,4) = vec(6)
        ke(5,5) = vec(7)
        ke(6,6) = vec(8)
        ke(7,7) = vec(1)
        ke(8,8) = vec(3)
        ke(9,9) = vec(2)
        ke(10,10) = vec(6)
        ke(11,11) = vec(7)
        ke(12,12) = vec(8)
C
        ke(5,3) = -vec(4)
        ke(3,5) = ke(5,3)
        ke(6,2) = vec(5)
        ke(2,6) = ke(6,2)
        ke(7,1) = -vec(1)
        ke(1,7) = ke(7,1)
        ke(8,2) = -vec(3)
        ke(2,8) = ke(8,2)
        ke(8,6) = -vec(5)
        ke(6,8) = ke(8,6)
        ke(9,3) = -vec(2)
        ke(3,9) = ke(9,3)
        ke(9,5) = vec(4)
        ke(5,9) = ke(9,5)
        ke(10,4) = -vec(6)
        ke(4,10) = ke(10,4)
        ke(11,3) = -vec(4)
        ke(3,11) = ke(11,3)
        ke(11,5) = vec(9)
        ke(5,11) = ke(11,5)
        ke(11,9) = vec(4)
        ke(9,11) = ke(11,9)
        ke(12,2) = vec(5)
        ke(2,12) = ke(12,2)
        ke(12,6) = vec(10)
        ke(6,12) = ke(12,6)
        ke(12,8) = -vec(5)
        ke(8,12) = ke(12,8)
C
        do 120 ii = 1, 12
          do 120 jj = 1,12
            ke(jj,ii) = ke(ii,jj)
120     continue
C
C.... COMPUTE THE TRANFORMATION MATRIX  
C
        call frame6(eframe,u,v,w)
C
        t33(1,1) = u(1)
        t33(1,2) = u(2)
        t33(1,3) = u(3)
        t33(2,1) = v(1)
        t33(2,2) = v(2)
        t33(2,3) = v(3)
        t33(3,1) = w(1)
        t33(3,2) = w(2)
        t33(3,3) = w(3)
C
        do 130 kk = 1, 4
          ic = 3*(kk-1)
          do 140 ii = 1, 3
            do 140 jj = 1, 3
              tran(ic+ii,ic+jj) = t33(ii,jj)
140       continue
130     continue
C
C.... ROTATE THE GLOBAL DISPLACEMENT VECTOR BACK TO LOCAL FRAME
C
        call dgemv('N',12,12,1.0d00,tran,12,ug,1,0.0d00,tug,1)
C
C.... COMPUTE THE INTERNAL FORCE RESULTANTS
C
        call dgemv('N',12,12,1.0d00,ke,12,tug,1,0.0d00,res,1)

C
C.... COMPUTE THERMAL STRESS

       do i=1,2
        tl(i) = temp(i) - tref
       enddo

       tstress = e*alpha*area*(tl(1)+tl(2))/2.
C
C.... WRITE THE RESULTANTS INTO THE STRESS ARRAY WHERE
C                 
C        FORX = STRESSXX      FORY=STRESSYY     FORZ=STRESSZZ
C        MOMX = STRESSXY      MOMY=STRESSXZ     MOMZ=STRESSYZ
C
        stress(elm,1,1) = res(1) + tstress
        stress(elm,2,1) = res(2)
        stress(elm,3,1) = res(3)
        stress(elm,4,1) = res(4)
        stress(elm,5,1) = res(5)
        stress(elm,6,1) = res(6)
        stress(elm,1,2) = res(7) - tstress 
        stress(elm,2,2) = res(8)
        stress(elm,3,2) = res(9)
        stress(elm,4,2) = res(10)
        stress(elm,5,2) = res(11)
        stress(elm,6,2) = res(12)
C
C.... WRITE THE VON MISES RESULTANT INTO THE STRESS ARRAY
C
        stress_0dxy = stress(elm,1,1) - stress(elm,2,1)
        stress_0dyz = stress(elm,2,1) - stress(elm,3,1)
        stress_0dxz = stress(elm,1,1) - stress(elm,3,1)
        stress_1dxy = stress(elm,1,2) - stress(elm,2,2)
        stress_1dyz = stress(elm,2,2) - stress(elm,3,2)
        stress_1dxz = stress(elm,1,2) - stress(elm,3,2)
C
        stress(elm,7,1) = dsqrt( 0.5d00*(stress_0dxy*stress_0dxy + 
     +        stress_0dyz*stress_0dyz + stress_0dxz*stress_0dxz) + 
     +                   3.0d00*(stress(elm,4,1)*stress(elm,4,1) + 
     +                           stress(elm,5,1)*stress(elm,5,1) + 
     +                           stress(elm,6,1)*stress(elm,6,1)) )
        stress(elm,7,2) = dsqrt( 0.5d00*(stress_1dxy*stress_1dxy + 
     +        stress_1dyz*stress_1dyz + stress_1dxz*stress_1dxz) + 
     +                   3.0d00*(stress(elm,4,2)*stress(elm,4,2) + 
     +                           stress(elm,5,2)*stress(elm,5,2) + 
     +                           stress(elm,6,2)*stress(elm,6,2)) )
C
c       write(6,*) 'Element Number :',elm
c       write(6,*) 'Force_x :' , res(1),res(7), stress(elm,1,1)
c       write(6,*) 'Moment_x :', res(4),res(10)
c       write(6,*) 'Moment_z :', res(6),res(12)
C
        return
        end
