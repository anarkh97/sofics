       subroutine sands8(xl,yl,zl,e,nu,h,v,stress,
     &            strainFlg,maxsze,maxstr,maxgus,elm,surface,thrmStr)
C-----------------------------------------------------------------
C This subroutine evaluates the von mises stress at the centroid
C of the triangle. it computes the values at the top and bottom
C and picks the maximum as output
C for the spatial 18 d.o.f three node triangle obtained as
C a combination of the aqr bending triangle plus the membrane
C with driling d.o.f. developed by Felippa et al.
C
C MODIFIED: 10-07-97 	K. H. Pierson
C Added stress calculations for: sigmaxx, sigmayy, sigmaxy,
C epsilonxx, epsilonyy, epsilonzz, epsilonxy and an equivalent
C strain similar to the vonmises stress. These stresses and strains
C can be calculated at the top, median or bottom surfaces.
C Stresses and strains are computed locally and then transformed
C to the global coordinate system.
C 
*-------------------------------------------------------------------*
*  CALLED BY : ThreeNodeShell.C 
*
*  SUBROUTINES CALLED:
*                      ROTATION
*                      TRIROT
*                      MEMBRA
*                      MOMEN
*			  TRANSFORMATION
*                      VONMISES
*-------------------------------------------------------------------* 
C
C t = thickness of triangle
C
       integer maxsze,maxstr,maxgus
       integer elm,surface
       double precision  e,nu,thrmStr
       double precision  epsxx,epsyy,epszz,epsxy
       double precision  rnxt,rnyt,rnxyt,rmxt,rmyt,rmxyt
       double precision  t2,sx,sy,sxy
       double precision  xl(*),yl(*),zl(*),h(*),v(*)
       double precision  stress(maxsze,maxstr,maxgus)
       double precision  db(3,3),dm(3,3)
       double precision  xp(3),yp(3),zp(3),xlp(3),ylp(3),zlp(3)
       double precision  r1(3,3),dll(18)
       double precision  xg(3),yg(3),zg(3),str(6)
       double precision  cb,x21,y21,z21,x32,y32,z32,x13,y13,z13
       double precision  rx,ry,rz,bx,by,bz,rlr,rlb,bpr,area,ylr,zlr
       double precision  ycg,xcg,zcg,xlcg,ylcg,zlcg,f,t
       double precision  rmom(18,3),rmem(18,3)
       double precision  rmx,rmy,rmxy,rnx,rny,rnxy,clr,cqr
       double precision  rmmx,rmmy,rmmxy,rnnx,rnny,rnnxy,sbf
       double precision  ebar
       double precision  factor
       integer lb(9),le(9), strainFlg
       character*10 status
       integer i,j
       data lb/3,4,5,9,10,11,15,16,17/
       data le/1,2,7,8,13,14,6,12,18/
C
       do 5 i = 1, 18
         rmom(i,1) = 0.0d0
         rmom(i,2) = 0.0d0
         rmom(i,3) = 0.0d0
         rmem(i,1) = 0.0d0
         rmem(i,2) = 0.0d0
         rmem(i,3) = 0.0d0
5      continue
C
        f   = 1.0d+00
        clr = 0.0d+00
        cqr = 1.0d+00
        t = h(1)
C
C set the bending constitutive matrix
C
   	cb=e*(t*t*t)/12.0/(1.0-nu*nu)
        db(1,1) =     cb
        db(1,2) =  nu*cb
        db(1,3) = 0.0
        db(2,1) =  db(1,2)
        db(2,2) =     cb
        db(2,3) = 0.0
        db(3,1) = 0.0
        db(3,2) = 0.0
        db(3,3) = 0.5*(1.0-nu)*cb
C
C set the membrane constitutive matrix
C
   	cb=e*t/(1.0-nu*nu)
        dm(1,1) =     cb
        dm(1,2) =  nu*cb
        dm(1,3) = 0.0
        dm(2,1) =  dm(1,2)
        dm(2,2) =     cb
        dm(2,3) = 0.0
        dm(3,1) = 0.0
        dm(3,2) = 0.0
        dm(3,3) = 0.5*(1.0-nu)*cb
C
C triangular dimension variables
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
C  triangle in space : we compute the length of one side and the distance of the
C  opposing node to that side to compute the area
       rx   = x21
       ry   = y21
       rz   = z21
       bx   = x32
       by   = y32
       bz   = z32
       rlr  = dsqrt( rx*rx + ry*ry + rz*rz )
       rlb  = dsqrt( bx*bx + by*by + bz*bz )
       bpr  = dsqrt((rx * bx + ry * by + rz *bz )**2)/rlr
       area = 0.5*rlr*dsqrt(rlb*rlb-bpr*bpr)
C
C Direction cosines of the local system . X' is directed parallel to the 
C 2-1 side Z' is the external normal (counterclockwise). 
C Y' computed as Z' x X'
C
       xp(1) = x21/rlr
       xp(2) = y21/rlr
       xp(3) = z21/rlr
C
C Z' local axis
C
       zp(1) = y21 * z32 - z21 * y32
       zp(2) = z21 * x32 - x21 * z32
       zp(3) = x21 * y32 - y21 * x32
       zlr   = dsqrt( zp(1)*zp(1) + zp(2)*zp(2) + zp(3)*zp(3) )
       zp(1) = zp(1)/zlr
       zp(2) = zp(2)/zlr
       zp(3) = zp(3)/zlr
C
C Y' local axis
C
       yp(1) = zp(2) * xp(3) - zp(3) * xp(2)
       yp(2) = zp(3) * xp(1) - zp(1) * xp(3)
       yp(3) = zp(1) * xp(2) - zp(2) * xp(1)
       ylr   = dsqrt( yp(1)*yp(1) + yp(2)*yp(2) + yp(3)*yp(3) )
       yp(1) = yp(1)/ylr
       yp(2) = yp(2)/ylr
       yp(3) = yp(3)/ylr
C
C center of gravity
C
       xcg = (xl(1) + xl(2) + xl(3))/3.0d+00
       ycg = (yl(1) + yl(2) + yl(3))/3.0d+00
       zcg = (zl(1) + zl(2) + zl(3))/3.0d+00
C
C computing local coordinates
C
       do 43 i=1,3
         xlcg   = xl(i) - xcg
         ylcg   = yl(i) - ycg
         zlcg   = zl(i) - zcg
         xlp(i) = xp(1) * xlcg + xp(2) * ylcg + xp(3) * zlcg
         ylp(i) = yp(1) * xlcg + yp(2) * ylcg + yp(3) * zlcg
         zlp(i) = zp(1) * xlcg + zp(2) * ylcg + zp(3) * zlcg
 43    continue
C
C Set Global axes
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
C
C computing nodal rotation matrix
C
       call rotation(xp,yp,zp,xg,yg,zg,r1)
C
C compute the von mises stress
C rotate nodal displacements to local system
C
       do 270 i=1,18
          dll(i)= 0.0d+00
270    continue
       do 280 i=1,3
         do 280 j=1,3
           dll(i)    = dll(i)    + r1(j,i)*v(j)
           dll(i+3)  = dll(i+3)  + r1(j,i)*v(j+3)
           dll(i+6)  = dll(i+6)  + r1(j,i)*v(j+6)
           dll(i+9)  = dll(i+9)  + r1(j,i)*v(j+9)
           dll(i+12) = dll(i+12) + r1(j,i)*v(j+12)
           dll(i+15) = dll(i+15) + r1(j,i)*v(j+15)
280    continue

c      compute centroidal membrane strains
       call membra(xlp,ylp,1.5d+00,le,rmem,status)

c      compute centroidal bending strains (curvatures (1/radius))
       call momen (xlp,ylp,lb,rmom,status)

c      pjsa 8/11/2011: the matrix rmom returned by momen function is
c      off by a factor of -1/sqrt(area) 
       factor = -1.0d0 / dsqrt( area )

       rmx  = 0.0d0
       rmy  = 0.0d0
       rmxy = 0.0d0
       rnx  = 0.0d0
       rny  = 0.0d0
       rnxy = 0.0d0
       do 290 j=1,18
          rmx  =  rmx + factor*rmom(j,1)*dll(j)
          rmy  =  rmy + factor*rmom(j,2)*dll(j)
          rmxy = rmxy + factor*rmom(j,3)*dll(j)
          rnx  =  rnx + rmem(j,1)*dll(j)
          rny  =  rny + rmem(j,2)*dll(j)
          rnxy = rnxy + rmem(j,3)*dll(j)
290    continue
C
C      COMPUTE VON MISES STRAIN RESULTANT
C
        if(strainFlg .eq. 1) then

          t2 = 0.5*t

          rmx  = t2*rmx
          rmy  = t2*rmy
          rmxy = t2*rmxy

C ... COMPUTE STRAINS AT MEDIAN SURFACE
          if(surface .eq. 2) then
            epsxx = rnx
            epsyy = rny
            epsxy = 0.5*rnxy
C ... COMPUTE STRAINS AT BOTTOM SURFACE
          else if(surface .eq. 3) then
            epsxx =  rnx -  rmx
            epsyy =  rny -  rmy
            epsxy = 0.5*(rnxy - rmxy)
C ... COMPUTE STRAINS AT TOP SURFACE
          else
            epsxx =  rnx +  rmx
            epsyy =  rny +  rmy
            epsxy = 0.5*(rnxy + rmxy)
          endif

          epszz  = -nu/(1.0 - nu)*(epsxx + epsyy)

          str(1) = epsxx
          str(2) = epsyy
          str(3) = epszz
          str(4) = epsxy
          str(5) = 0.0
          str(6) = 0.0

          call transform(xp,yp,zp,xg,yg,zg,str)

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

          call straineq(rmx,rmy,rmxy,rnx,rny,rnxy,nu,surface,ebar)

          stress(elm,7,1) = ebar 
          stress(elm,7,2) = ebar
          stress(elm,7,3) = ebar 

          return
        endif
C
C     subtract off thermal strain portions
C     pjsa 6/11/2014: this should be not affect the return values when strainFlg.eq.1
      rnx = rnx - thrmStr
      rny = rny - thrmStr
C
C     compute centroidal stress resultants
C
C     bending resultants
      rmmx  = db(1,1)*rmx + db(1,2)*rmy + db(1,3)*rmxy
      rmmy  = db(2,1)*rmx + db(2,2)*rmy + db(2,3)*rmxy
      rmmxy = db(3,1)*rmx + db(3,2)*rmy + db(3,3)*rmxy

C     membrane resultants
      rnnx  = dm(1,1)*rnx + dm(1,2)*rny + dm(1,3)*rnxy
      rnny  = dm(2,1)*rnx + dm(2,2)*rny + dm(2,3)*rnxy
      rnnxy = dm(3,1)*rnx + dm(3,2)*rny + dm(3,3)*rnxy 

      t = abs(t)
      t2 = t*t

      rnxt  =  rnnx/t
      rnyt  =  rnny/t
      rnxyt = rnnxy/t

      rmxt  = 6.0* rmmx / t2
      rmyt  = 6.0* rmmy / t2
      rmxyt = 6.0*rmmxy / t2

C COMPUTE SIGMAXX, SIGMAYY AND SIGMAXY AT TOP SURFACE (DEFAULT)
      sx  =  rnxt +  rmxt
      sy  =  rnyt +  rmyt
      sxy = rnxyt + rmxyt

C COMPUTE SIGMAXX, SIGMAYY AND SIGMAXY AT MEDIAN SURFACE
      if(surface .eq. 2) then
        sx  = rnxt
        sy  = rnyt
        sxy = rnxyt
      endif

C COMPUTE SIGMAXX, SIGMAYY AND SIGMAXY AT BOTTOM SURFACE
      if(surface .eq. 3) then
        sx  =  rnxt -  rmxt
        sy  =  rnyt -  rmyt
        sxy = rnxyt - rmxyt
      endif

      str(1) = sx
      str(2) = sy
      str(3) = 0.0 
      str(4) = sxy
      str(5) = 0.0 
      str(6) = 0.0

      call transform(xp,yp,zp,xg,yg,zg,str)

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

C     compute von mises stress resultant

      call vonmis(rmmx,rmmy,rmmxy,rnnx,rnny,rnnxy,t,sbf,surface)

      stress(elm,7,1) = sbf
      stress(elm,7,2) = sbf
      stress(elm,7,3) = sbf

      return
      end
