        subroutine vmelmvc(stress,maxgus,maxstr,msize,elm,nno)
*
***************************************************************
*                                                             *
*	THIS SUBROUTINE COMPUTE THE COMPLEX VON MISES STRESS  *
*       AT THE NODES OF THE FOLLOWING ELEMENTS                *
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
        double complex stress(msize,maxstr,maxgus)

C
C.... DECLARE ALL LOCAL VARIABLES
C
        integer n
        double complex  sxx,syy,szz,sxy,sxz,syz
        double complex  dsxx,dsyy,dszz,dsxy,dsxz,dsyz
        double complex  j2,comp,vms
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
          vms = vms +  sqrt(3.0d0 * j2)
C
   10   continue
C
        vms = vms/nno
C
C.... DISTRIBUTE OUT TO THE NODES 
C
      do 11 n = 1, nno
        stress(elm,7,n) = vms
   11 continue
C
        return
        end

