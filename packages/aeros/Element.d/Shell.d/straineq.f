        subroutine straineq(rmx,rmy,rmxy,rnx,rny,rnxy,nu,surface,ebar)
***************************************************************
*       THIS ROUTINE CALCULATES THE EQUIVALENT VON MISES      *
*       STRAIN FOR THE 3 NODE SHELL ELEMENT                   *
***************************************************************
*                                                             *
*       AUTHOR  :       K.H. PIERSON                          *
*       DATE    :       MARCH 1997                            *
*       VERSION :       FEM-C++ 1.00                          *
*                                                             *
***************************************************************
*                                                             *
*	rmx  = centroidal bending curvature (kxxc)            *
*	rmy  = centroidal bending curvature (kyyc)            *
*	rmxy = centroidal bending curvature (kxyc)            *
*	rnx  = centroidal membrane strain   (exxc)            *
*	rny  = centroidal membrane strain   (eyyc)            *
*	rnxy = centroidal membrane strain   (exyc)            *
*	nu   = Poisson's ratio                                *
*	ebar = equivalent strain                              *
*                                                             *
***************************************************************
*                                                             *
*       CALLED BY :  SANDS8.F                                 *
*                                                             *
***************************************************************

C ... ARGUMENTS
        integer surface
        double precision rmx,rmy,rmxy,rnx,rny,rnxy,nu,ebar

C ... LOCAL VARIABLES
        double precision ex,ey,ez,exy,etop,ebot,emid

C ... CALCULATE STRAINS AT MEDIAN SURFACE,
C ... COMPUTE EQUIVALENT STRAIN AT MEDIAN SURFACE,
C ... AND RETURN IF MEDIAN SURFACE IS REQUESTED
        if(surface .eq. 2) then
          ex  = rnx
          ey  = rny
          ez  = -nu/(1.0-nu)*(ex + ey)
          exy = 0.5*rnxy
          call equiv(ex,ey,ez,exy,emid)
          ebar = emid
          return
        endif

C ... CALCULATE STRAINS AT TOP SURFACE
        ex  = rnx + rmx
        ey  = rny + rmy
        ez  = -nu/(1.0-nu)*(ex + ey)
        exy = 0.5*(rnxy + rmxy)

C ... COMPUTE EQUIVALENT STRAIN AT TOP SURFACE
        call equiv(ex,ey,ez,exy,etop)

C ... RETURN IF TOP SURFACE VALUE IS REQUESTED
        if(surface .eq. 1) then
          ebar = etop
          return
        endif

C ... CALCULATE STRAINS AT BOTTOM SURFACE
        ex   = rnx - rmx
        ey   = rny - rmy
        ez   = -nu/(1.0-nu)*(ex + ey)
        exy  = 0.5*(rnxy - rmxy)

C ... COMPUTE EQUIVALENT STRAIN AT BOTTOM SURFACE
        call equiv(ex,ey,ez,exy,ebot)

C ... RETURN IF BOTTOM SURFACE VALUE IS REQUESTED
        if(surface .eq. 3) then
          ebar = ebot
          return
        endif

C ... RETURN THE MAXIMUM EQUIVALENT STRAIN
        ebar = max(etop,ebot)

        return
        end

C
C ... SUBROUTINE TO CALCULATE EQUIVALENT STRAIN
C
        subroutine equiv(ex,ey,ez,exy,eq)

C ... ARGUMENTS
        double precision ex,ey,ez,exy,eq

C ... LOCAL VARIABLES
        double precision e0,dex,dey,dez
        
C ... COMPUTE MEAN HYDROSTATIC STRAIN
        e0 = (ex + ey + ez)/3.0d0

C ... COMPUTE DEVIATORIC STRAINS
        dex = ex - e0
        dey = ey - e0
        dez = ez - e0

C ... COMPUTE EQUIVALENT STRAIN
        eq = ((dex*dex + dey*dey + dez*dez)/2.0d0) + (exy*exy)
        eq = dsqrt(3.0d0 * eq)

        return
        end
