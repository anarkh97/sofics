        subroutine vmelmv(stress,maxgus,maxstr,msize,elm,nno)
*
***************************************************************
*                                                             *
*	THIS SUBROUTINE COMPUTE THE VON MISES STRESS AT THE   *
*       NODES OF THE FOLLOWING ELEMENTS                       *
*                                                             *
*       Quad-4, Brick-8, Tetra-4, Tetra-10, Penta-6           *
*                                                             *
***************************************************************
*
*	AUTHOR	:	C. FARHAT
*	DATE	:	SEPT. 96
*
*       modified by:   
*
*	AUTHOR	:	K. MAUTE
*	DATE	:	NOV. 2000
*
***************************************************************
*
*	STRESS = 3-D ARRRAY HOLDING ELEMENTAL STRESSES
*	MAXGUS = THIRD DIMENSION OF STRESS
*	MAXSTR = SECOND DIMENSION OF STRESS ARRAY
*	 MSIZE = LEADING DIMENSION OF STRESS ARRAY
*	   ELM = CURRENT ELEMENT NUMBER
*
***************************************************************
*
*	CALLED BY :  MDERIV.F
*
***************************************************************
C
C.... DECLARE ALL GLOBAL VARIABLES
C
C.... INTEGER CONSTANTS
C
        integer maxgus,maxstr,msize,elm,nno
C
C.... REAL ARRAYS
C
        real*8 stress(msize,maxstr,maxgus)

C
C.... DECLARE ALL LOCAL VARIABLES
C
        integer n
        real*8  sxx,syy,szz,sxy,sxz,syz
        real*8  dsxx,dsyy,dszz,dsxy,dsxz,dsyz
        real*8  j2,comp,vms
c
        vms = 0.0d0
C
C.... LOOP OVER ALL NODES
C
        do 10 n = 1, nno
c
          sxx = stress(elm,1,n)
          syy = stress(elm,2,n)
          szz = stress(elm,3,n)
          sxy = stress(elm,4,n)
          syz = stress(elm,5,n)
          sxz = stress(elm,6,n)
C
C.... COMPUTE THE FIRST DEVEATORIC STRESSES
C
          comp = (sxx + syy + szz)/3.0d0
          dsxx = sxx - comp
          dsyy = syy - comp
          dszz = szz - comp
          dsxy = sxy
          dsyz = syz
          dsxz = sxz
C
C.... COMPUTE THE SECOND DEVEATORIC STRESS
C
          j2 = ((dsxx*dsxx)+(dsyy*dsyy)+(dszz*dszz))/2.0d0 +
     &          (dsxy*dsxy)+(dsyz*dsyz)+(dsxz*dsxz)
C
C.... COMPUTE THE VON MISES STRESS
C
C          vms = vms +  dsqrt(3.0d0 * j2)
CC
C   10   continue
CC
C        vms = vms/nno
CC
CC.... DISTRIBUTE OUT TO THE NODES 
CC
C      do 11 n = 1, nno
C        stress(elm,7,n) = vms
C   11 continue
CC
C PJSA 10/7/2010 return the unaveraged nodal von mises effective stresses
C If averaging is required it will be done elsewhere
          stress(elm,7,n) = dsqrt(3.0d0 * j2)
C
   10   continue
        return
        end

