C -------------------------------------------------------------- C
C subroutine to compute the equivalent von mises stress          C
C in a plate or shell in the top surface (z=+t/2)                C
C and in the bottom surface (z=-t/2)                             C
C returns maximum value of the two stresses                      C
C                                                                C
C MODIFIED: March 3, 1997                                        C 
C BY: K. H. PIERSON                                              C
C WHY: the z component deviatoric stress was ignored previously  C
C                                                                C
C -------------------------------------------------------------- C
      subroutine vonmis(rmx,rmy,rmxy,rnx,rny,rnxy,t,sv,surface)

C.... GLOBAL VARIABLES
      integer surface
      double precision rmx,rmy,rmxy,rnx,rny,rnxy,t,sv

C.... LOCAL VARIABLES
C st = von mises stress in top surface
C sm = von mises stress in median surface
C sb = von mises stress in bottom surface

      double precision sx,sy,sxy,st,sb,sm,t2,sq3
      double precision rnxt,rnyt,rnxyt,rmxt,rmyt,rmxyt

      t    = abs(t)
      t2   = t*t

      sq3  = dsqrt(3.0d0)

      rnxt  =  rnx/t
      rnyt  =  rny/t
      rnxyt = rnxy/t

      rmxt  = 6.0d0 * rmx / t2
      rmyt  = 6.0d0 * rmy / t2
      rmxyt = 6.0d0 *rmxy / t2

C ... COMPUTE VON MISES STRESS IN BOTTOM SURFACE

      sx  =  rnxt -  rmxt
      sy  =  rnyt -  rmyt
      sxy = rnxyt - rmxyt

      call compj2(sx,sy,sxy,sb)

      sb = sq3 * sb

      if(surface .eq. 3) then 
        sv = sb
        return
      endif

C ... COMPUTE VON MISES STRESS IN MEDIAN SURFACE

      if(surface .eq. 2) then

        sx  = rnxt
        sy  = rnyt
        sxy = rnxyt

        call compj2(sx,sy,sxy,sm)

        sv = sq3 * sm

        return

      endif

C ... COMPUTE VON MISES STRESS IN TOP SURFACE

      sx  =  rnxt +  rmxt
      sy  =  rnyt +  rmyt
      sxy = rnxyt + rmxyt

      call compj2(sx,sy,sxy,st)

      st = sq3 *st

      if(surface .eq. 1) then 
        sv = st
        return
      endif

C ... VON MISES STRESS = max(sb,st)

      sv = max(sb,st)

      return
      end
C
C ... SUBROUTINE TO CALCULATE J2
C
      subroutine compj2(sx,sy,sxy,svm)

C ... GLOBAL VARIABLES
      double precision sx,sy,sxy,svm

C ... LOCAL VARIABLES
      double precision sz,s0,dsx,dsy,dsz,j2

C ... SET sz = 0 TO REMIND USER OF THIS ASSUMPTION
      sz = 0.0d0

C ... COMPUTE AVERAGE HYDROSTATIC STRESS
      s0 = (sx + sy + sz)/3.0d0

C ... COMPUTE DEVIATORIC STRESSES
      dsx  = sx - s0
      dsy  = sy - s0
      dsz  = sz - s0 
C
C ... RCFEM IGNORED THE DEVIATORIC Z COMPONENT: 
C     dsz  = 0.0

C ... COMPUTE J2

      j2 = 0.5*((dsx*dsx) + (dsy*dsy) + (dsz*dsz)) + (sxy*sxy)

C ... COMPUTE VON MISES STRESS
      svm = dsqrt(j2)

      return
      end
