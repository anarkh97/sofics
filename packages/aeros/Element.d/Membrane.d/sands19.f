        subroutine sands19(xl,yl,zl,e,nu,h,v,stress,msize
     .,                      maxstr,maxgus,elm,strainFlg)
*********************************************************************
*
* THIS SUBROUTINE EVALUATES THE VON MISES STRESS AT THE CENTROID
* OF THE TRIANGLE. IT COMPUTES THE VALUES AT THE TOP AND BOTTOM
* AND PICKS THE MAXIMUN AS OUTPUT
* 
*-------------------------------------------------------------------*
*  CALLED BY : Membrane.C 
*
*  SUBROUTINES CALLED:
*                      ROTATION
*                      TRIROT
*                      MEMBRA
*                      VONMISES
*-------------------------------------------------------------------* 
C
C.... DECLARE ALL GLOBAL VARIABLES
C
C.... INTEGER CONSTANTS
C 
        integer msize,maxstr,maxgus,elm,strainFlg
C
C.... REAL CONSTANTS
C
        double precision e,nu
C
C.... REAL ARRAYS
C
        real*8 xl(*),yl(*),zl(*),h(*),v(*)
        real*8 stress(msize,maxstr,maxgus)
C
C.... DECALRE ALL LOCAL VARIABLES
C
        integer nd,i,j
        integer le(9) 
        double precision dm(3,3),r1(3,3)
        double precision xp(3),yp(3),zp(3),xlp(3),ylp(3),zlp(3)
        double precision dll(18),v1n(3),v2n(3),v3n(3)
        double precision cb,x21,y21,z21,x32,y32,z32,x13,y13,z13
        double precision rx,ry,rz,bx,by,bz,rlr,rlb,bpr,area,ylr,zlr
        double precision ycg,xcg,zcg,xlcg,ylcg,zlcg,t
        double precision trM(6,6)
        double precision rmem(18,3)
        double precision rmx,rmy,rmxy,rnx,rny,rnxy,ebar
        double precision rmmx,rmmy,rmmxy,rnnx,rnny,rnnxy,sbf
        double precision xg(3),yg(3),zg(3),str(6)
        character*10 status
C
C.... INIALIZE DATA
C
        data le/1,2,7,8,13,14,6,12,18/
        nd=18
        t  = h(1)
   	cb=e*t/(1.0-nu*nu)

         do 5 i = 1, 18
           do 5 j = 1, 3
           rmem(i,j) = 0.0d0
5        continue

C
C.... INITIALIZE  ELASTIC MATRIX
C
         do 10 i=1,3
           do 10 j=1,3
             dm(i,j)=0.0
10       continue

         dm(1,1) =    cb
         dm(1,2) = nu*cb
         dm(2,1) = nu*cb
         dm(2,2) =    cb
         dm(3,3) = (1.0-nu)/2.0*cb
C
C.... DIMENSION VARIABLES
C
         x21 = xl(2) - xl(1)
         y21 = yl(2) - yl(1)
         z21 = zl(2) - zl(1)
         x32 = xl(3) - xl(2)
         y32 = yl(3) - yl(2)
         z32 = zl(3) - zl(2)
         x13 = xl(1) - xl(3)
         y13 = yl(1) - yl(3)
         z13 = zl(1) - zl(3)
C
C....  TRIANGLE IN SPACE : WE COMPUTE THE LENGTH 
C....  OF ONE SIDE AND THE DISTANCE OF THE
C....  OPPOSING NODE TO THAT SIDE TO COMPUTE THE AREA
C
      rx   = x21
      ry   = y21
      rz   = z21
      bx   = x32
      by   = y32
      bz   = z32
      rlr  = dsqrt( rx*rx + ry*ry + rz*rz )
      rlb  = dsqrt( bx*bx + by*by + bz*bz )
      bpr  = dsqrt((rx * bx + ry * by + rz *bz )**2)/rlr
      area = rlr*(dsqrt(rlb**2-bpr**2))/2.0d+00
C
C.... DIRECTION COSINES OF THE LOCAL SYSTEM . 
C.... X' IS DIRCTED PARALLEL TO THE 2-1 SIDE
C.... Z' IS THE EXTERNAL NORMAL (COUNTERCLOCKWISE). 
C.... Y' COMPUTED AS Z'x X'
C
      xp(1) = x21/rlr
      xp(2) = y21/rlr
      xp(3) = z21/rlr
C
C.... Z'
C
      zp(1) = y21 * z32 - z21 * y32
      zp(2) = z21 * x32 - x21 * z32
      zp(3) = x21 * y32 - y21 * x32
      zlr   = dsqrt( zp(1)**2 + zp(2)**2 + zp(3)**2 )
      zp(1) = zp(1)/zlr
      zp(2) = zp(2)/zlr
      zp(3) = zp(3)/zlr
C
C.... Y'
C
      yp(1) = zp(2) * xp(3) - zp(3) * xp(2)
      yp(2) = zp(3) * xp(1) - zp(1) * xp(3)
      yp(3) = zp(1) * xp(2) - zp(2) * xp(1)
      ylr   = dsqrt( yp(1)**2 + yp(2)**2 + yp(3)**2 )
      yp(1) = yp(1)/ylr
      yp(2) = yp(2)/ylr
      yp(3) = yp(3)/ylr
C
C.... CENTER OF GRAVITY
C
      xcg = (xl(1) + xl(2) + xl(3))/3.0d+00
      ycg = (yl(1) + yl(2) + yl(3))/3.0d+00
      zcg = (zl(1) + zl(2) + zl(3))/3.0d+00
C
C.... COMPUTING LOCAL COORDINATES 
C
      do 43 i=1,3
        xlcg   = xl(i)-xcg
        ylcg   = yl(i)-ycg
        zlcg   = zl(i)-zcg
        xlp(i) = xp(1) * xlcg + xp(2) * ylcg + xp(3) * zlcg
        ylp(i) = yp(1) * xlcg + yp(2) * ylcg + yp(3) * zlcg
        zlp(i) = zp(1) * xlcg + zp(2) * ylcg + zp(3) * zlcg
43    continue
C
C.... COMPUTING NODAL ROTATION MATRICES
C
        v1n(1)= 1.0d+00
        v1n(2)= 0.0d+00
        v1n(3)= 0.0d+00
        v2n(1)= 0.0d+00
        v2n(2)= 1.0d+00
        v2n(3)= 0.0d+00
        v3n(1)= 0.0d+00
        v3n(2)= 0.0d+00
        v3n(3)= 1.0d+00

        call rotation(xp,yp,zp,v1n,v2n,v3n,r1)
C
C.... COMPUTE THE VON MISES STRESS
C
C.... ROTATE NODAL DISPLACEMENTS TO LOCAL SYSTEM
C
        do 270 i=1,18
          dll(i)= 0.0d+00
270     continue

        do 280 i=1,3
          do 280 j=1,3
            dll(i)   =dll(i)   +r1(j,i)*v(j)
            dll(i+3) =dll(i+3) +r1(j,i)*v(j+3)
            dll(i+6) =dll(i+6) +r1(j,i)*v(j+6)
            dll(i+9) =dll(i+9) +r1(j,i)*v(j+9)
            dll(i+12)=dll(i+12)+r1(j,i)*v(j+12)
            dll(i+15)=dll(i+15)+r1(j,i)*v(j+15)
280     continue
C
C.... COMPUTE CENTROIDAL MEMBRANE STRAINS
C
        call membra(xlp,ylp,1.5d+00,le,rmem,status)
C
        rmx  = 0.0
        rmy  = 0.0
        rmxy = 0.0
C 
        rnx  = 0.0
        rny  = 0.0
        rnxy = 0.0
C
        do 290 j=1,18
          rnx  = rnx  + rmem(j,1)*dll(j)
          rny  = rny  + rmem(j,2)*dll(j)
          rnxy = rnxy + rmem(j,3)*dll(j)
290    continue

C
C ...  COMPUTE EQUIVALENT STRAINS
C
       xg(1) = 1.0
       xg(2) = 0.0
       xg(3) = 0.0
       yg(1) = 0.0
       yg(2) = 1.0
       yg(3) = 0.0
       zg(1) = 0.0
       zg(2) = 0.0
       zg(3) = 1.0
       if(strainFlg .eq. 1) then
         str(1) = rnx
         str(2) = rny
         str(3) = 0.0
         str(4) = 0.5*rnxy
         str(5) = -(nu/(1.0-nu))*(rnx + rny)
         str(6) = 0.0
         call transform(xp,yp,zp,xg,yg,zg,str,trM)
         stress(elm,1,1) = str(1)
         stress(elm,1,2) = str(1)
         stress(elm,1,3) = str(1)
         stress(elm,2,1) = str(2)
         stress(elm,2,2) = str(2)
         stress(elm,2,3) = str(2)
         stress(elm,3,1) = str(3)
         stress(elm,3,2) = str(3)
         stress(elm,3,3) = str(3)
         stress(elm,4,1) = 2.0D0*str(4)
         stress(elm,4,2) = 2.0D0*str(4)
         stress(elm,4,3) = 2.0D0*str(4)
         stress(elm,5,1) = 2.0D0*str(5)
         stress(elm,5,2) = 2.0D0*str(5)
         stress(elm,5,3) = 2.0D0*str(5)
         stress(elm,6,1) = 2.0D0*str(6)
         stress(elm,6,2) = 2.0D0*str(6)
         stress(elm,6,3) = 2.0D0*str(6)
C
         ebar= ((rnx-rny)*(rnx-rny)+rnx*rnx+rny*rny)/6.0d0
         ebar = ebar + 0.25*(rnxy*rnxy)
         ebar = dsqrt(3.0d0*ebar)
C
C        ebar= dsqrt(2.0d0/3.0d0)*dsqrt((rnx*rnx + rny*rny + rnxy*rnxy))
C
         stress(elm,7,1) = ebar
         stress(elm,7,2) = ebar
         stress(elm,7,3) = ebar
         return
       end if

C
C....  COMPUTE CENTROIDAL STRESS RESULTANTS
C
       rmmx  = 0.0d0 
       rmmy  = 0.0d0 
       rmmxy = 0.0d0 

       rnnx  = dm(1,1)*rnx+dm(1,2)*rny+dm(1,3)*rnxy
       rnny  = dm(2,1)*rnx+dm(2,2)*rny+dm(2,3)*rnxy
       rnnxy = dm(3,1)*rnx+dm(3,2)*rny+dm(3,3)*rnxy 
C
C ...  COMPUTE VON MISES STRESS RESULTANT
C
       call vonmis(rmmx,rmmy,rmmxy,rnnx,rnny,rnnxy,t,sbf,1)
C
C ... COMPUTE SIGMAXX SIGMAYY AND SIGMAXY
C
      str(1) = rnnx/t
      str(2) = rnny/t
      str(3) = 0.0
      str(4) = rnnxy/t
      str(5) = 0.0
      str(6) = 0.0

      call transform(xp,yp,zp,xg,yg,zg,str,trM)

      stress(elm,1,1) = str(1)
      stress(elm,1,2) = str(1)
      stress(elm,1,3) = str(1)

      stress(elm,2,1) = str(2)
      stress(elm,2,2) = str(2)
      stress(elm,2,3) = str(2)

      stress(elm,3,1) = str(3)
      stress(elm,3,2) = str(3)
      stress(elm,3,3) = str(3)

      stress(elm,4,1) = str(4)
      stress(elm,4,2) = str(4)
      stress(elm,4,3) = str(4)

      stress(elm,5,1) = str(5)
      stress(elm,5,2) = str(5)
      stress(elm,5,3) = str(5)

      stress(elm,6,1) = str(6)
      stress(elm,6,2) = str(6)
      stress(elm,6,3) = str(6)

      stress(elm,7,1) = sbf
      stress(elm,7,2) = sbf
      stress(elm,7,3) = sbf

      return
      end
