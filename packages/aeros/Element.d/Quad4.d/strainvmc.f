        subroutine strainvmc(strain,maxgus,maxstr,msize,numnod)
***************************************************************
*       THIS ROUTINE CALCULATES EQUIVLAENT COMPLEX VON MISES  *
*       STRAIN                                                *
***************************************************************
*                                                             *
*	AUTHOR	:	K.H. PIERSON 			      *
*	DATE	:	MARCH 1997			      *
*	VERSION	:       FEM-C++ 1.00    		      *
*
*       modified by:
*
*       AUTHOR  :       G. BROWN
*       DATE    :       DEC. 2000
*                                                             *
***************************************************************
*                                                             *
*	STRAIN = 3-D ARRRAY HOLDING ELEMENTAL STRAINS         *
*	MAXGUS = THIRD DIMENSION OF STRAIN                    *
*	MAXSTR = SECOND DIMENSION OF STRAIN ARRAY             *
*	 MSIZE = LEADING DIMENSION OF STRAIN ARRAY            *
*       NUMNOD = # OF NODES IN THE ELEMENT                    *
*                                                             *
***************************************************************
*                                                             *
*	CALLED BY :  SANDS19.F                                *
*                                                             *
***************************************************************
C
C.... INTEGER CONSTANTS
C
        integer maxgus,maxstr,msize,numnod
C
C.... REAL CONSTANTS
C
        double precision nno
C
C.... REAL ARRAYS
C
        double complex strain(msize,maxstr,maxgus)
C
C.... DECLARE ALL LOCAL VARIABLES
C
        integer n
        double complex exx,eyy,ezz,exy,exz,eyz
        double complex dexx,deyy,dezz,dexy,dexz,deyz
        double complex j2,comp,vms
C
        vms = 0.0d0
        nno = int(numnod)
C
C.... LOOP OVER ALL NODES
C
        do 10 n = 1, numnod
C
          exx = strain(1,1,n)
          eyy = strain(1,2,n)
          ezz = strain(1,3,n)
C convert to strain tensor from engineering strain
          exy = strain(1,4,n)/2.0D0
          eyz = strain(1,5,n)/2.0D0
          exz = strain(1,6,n)/2.0D0

C
C ... COMPUTE THE MEAN HYDROSTATIC STRAIN
C
          comp = (exx + eyy + ezz)/3.0d0
C
C.... COMPUTE THE FIRST DEVIATORIC STRAINS
C
          dexx = exx - comp
          deyy = eyy - comp
          dezz = ezz - comp
          dexy = exy
          deyz = eyz
          dexz = exz
C
C.... COMPUTE THE SECOND DEVIATORIC STRAIN
C
          j2 = ((dexx*dexx)+(deyy*deyy)+(dezz*dezz))/2.0d0 +
     &          (dexy*dexy)+(deyz*deyz)+(dexz*dexz)
C
C.... COMPUTE THE VON MISES STRAIN
C
          vms = vms + sqrt(3.0d0 * j2)
C
 10     continue
C
        vms = vms/nno
C
C.... DISTRIBUTE OUT TO THE NODES
C
        do 20 n = 1, numnod 
          strain(1,7,n) = vms
20      continue
C
        return
        end
