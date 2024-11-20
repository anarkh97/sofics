        subroutine mass6(rho,Ix,Iy,Iz,elmass,a,x,y,z,gamma,grvfor
     .,                  grvflg,totmas,masflg,eframe)
*-----------------------------------------------------------------*
*THIS SUBROUTINE WILL FORM THE ELEMENT MASS MATRIX                *
*FOR THE 3-D EULER-BERNOULLI BEAM IN CONSISTENT FORM.             *
*******************************************************************
*                                                                 *
*  The form of the Mass Matrix is taken from                      *
*     Theory of Matrix Structural Analysis                        *
*     J.S. Przemieniecki                                          *
*                                                                 *
*******************************************************************
*
*		VARIABLES
*
*	    A = CROSS SECTIONAL AREA OF THE ELEMENT
*         RHO = ELEMENT MASS DENSITY 
*    Ix,Iy,Iz = Moments of Inertia 
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
*       eframe = Local Frame of Element
*
*
******************************************************************
*
*----------------------------------------------------------------*
C
C.... DECLARE THE GLOBAL VARIABLES 
C
C.... REAL CONSTANTS
C
        real*8 a,rho,totmas,Ix,Iy,Iz
C	real*8 Iztemp
C
C.... REAL ARRAYS
C
        real*8 elmass(12,12),x(*),y(*),z(*)
        real*8 gamma(*),grvfor(*)
        real*8 eframe(9)
C
C.... LOGICAL CONSTANT
C
        logical grvflg,masflg
C
C.... DECLARE LOCAL VARIABLES
C
        integer i,j 
        real*8 dx,dy,dz,len
        real*8 me(12,12), mass
        real*8 u(3),v(3),w(3)
        real*8 tran(12,12),trant(12,12),ttme(12,12),gmass(12,12)
C
C.... ZERO THE  ELMENT MASS MATRIX
CXX        write(*,*) 'MASS6.F: Consistent Mass for Beam'
CXX        call flush(6)
C
        do 5 i=1,12
          do 15 j=1,12
            elmass(i,j) = 0.d0
            me(i,j) = 0.d0
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
        mass = A*RHO*len
C
C.... CALCULATE THE LOCAL MASS MATRIX
C
C       write(6,*) "ZEROING IZ in EULERBEAM FOR TEMPORARY USE"
C       Iztemp = Iz
C       Iz = 0.d0

        me(1,1)   = 1.d0/3.d0
        me(1,7)   = 1.d0/6.d0
C
        me(2,2)   = (13.d0)/(35.d0) + (6.d0*Iz)/(5.d0*A*len*len)
        me(2,6)   = (11.d0*len)/(210.d0) + (Iz)/(10.d0*A*len)
        me(2,8)   = (9.d0)/(70.d0) - (6.d0*Iz)/(5.d0*A*len*len)
        me(2,12)  = (-13.d0*len)/(420.d0) + (Iz)/(10.d0*A*len)
C
        me(3,3)   = (13.d0)/(35.d0) + (6.d0*Iy)/(5.d0*A*len*len)
        me(3,5)   = (-11.d0*len)/(210.d0) - (Iy)/(10.d0*A*len)
        me(3,9)   = (9.d0)/(70.d0) - (6.d0*Iy)/(5.d0*A*len*len)
        me(3,11)  = (13.d0*len)/(420.d0) - (Iy)/(10.d0*A*len)
C
        me(4,4)   = Ix/(3.d0*A)
        me(4,10)  = Ix/(6.d0*A)
C
        me(5,5)   = (len*len)/(105.d0) + (2*Iy)/(15.d0*A)
        me(5,9)   = -1.d0*me(3,11) 
        me(5,11)  = (-1.d0*len*len)/(140.d0) - (Iy)/(30.d0*A)
C
        me(6,6)   = (len*len)/(105.d0) + (2*Iz)/(15.d0*A)
        me(6,8)   = -1.d0*me(2,12) 
        me(6,12)  = (-1.d0*len*len)/(140.d0) - (Iz)/(30.d0*A)
C
        me(7,7)   = me(1,1) 
C
        me(8,8)   = me(2,2) 
        me(8,12)  = -1.d0*me(2,6) 
C
        me(9,9)   = me(3,3) 
        me(9,11)  = -1.d0*me(3,5) 
C
        me(10,10) = me(4,4) 
C
        me(11,11) = me(5,5) 
C
        me(12,12) = me(6,6) 

C       Iz= Iztemp 
C
C.... USE SYMMETRY TO BUILD THE REST OF THE LOCAL MASS MATRIX 
C
        do 20 i=1,12
          do 10 j=1,i
            me(i,j) = me(j,i)
10        continue
20      continue
C
C
C Get Transformation Matrices
        call frame6(eframe,u,v,w)
        call transf(u,v,w,tran,trant)
C
C Transform Mass Matrix into Global Frame
C dgemm is a BLAS routine
C
      call dgemm('N','N',12,12,12,1.0d00,
     &           trant,12,me,12,0.0d00,ttme,12)
      call dgemm('N','N',12,12,12,1.0d00,
     &           ttme,12,tran,12,0.0d00,gmass,12)

C Copy into output and multiply by rho*area*lengh
        do 40 i=1,12
          do 30 j=1,12
            elmass(i,j) = mass*gmass(i,j)
30        continue
40      continue
C
C..... COMPUTE THE BODY FORCE DUE TO GRAVITY IF NEEDED
C
        if (grvflg) then
          grvfor(1) = mass*gamma(1)
          grvfor(2) = mass*gamma(2)
          grvfor(3) = mass*gamma(3)
        endif
C
C.... COMPUTE THE SUBDOMAIN MASS
C
        if (masflg) then
          totmas = totmas + mass
        endif
C
        return
        end
        
