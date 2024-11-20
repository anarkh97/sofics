CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C This subroutine computes the mass matrix for the
C AQR element, shell version. The main idea is to distri-
C bute the complete mass of the element over three
C fictitious beams allocated over the sides. The big
C problem is that we must construct the consistent mass
C matrix in the global system and then lump it.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C
        subroutine mass8(xl,yl,zl,h,dens,fv,eltdof,gamma,grvfor
     .,                                    grvflg,totmas,masflg)
        integer i,j,i1,i2,i3,eltdof,nd
        real*8 xl(3),yl(3),zl(3),h(3)
        real*8 fv(eltdof,eltdof),rl(3)
        real*8 x21,y21,z21,x32,y32,z32,x13,y13,z13
        real*8 rx,ry,rz,bx,by,bz,rlr,rlb,bpr,area,rmas
        real*8 dens,esp,totmas
        real*8 gamma(*),grvfor(*)
        logical grvflg,masflg
        nd=18
C initialize mass matrix to zero
          esp=h(1)
          do 10 i = 1, eltdof
            do 10 j = 1, eltdof
              fv(i,j) = 0.0d0
10	  continue
C dimension variables
          x21= xl(2)-xl(1)
          y21= yl(2)-yl(1)
          z21= zl(2)-zl(1)
          x32= xl(3)-xl(2)
          y32= yl(3)-yl(2)
          z32= zl(3)-zl(2)
          x13= xl(1)-xl(3)
          y13= yl(1)-yl(3)
          z13= zl(1)-zl(3)
          rl(1)=dsqrt(x21**2+y21**2+z21**2)
          rl(2)=dsqrt(x32**2+y32**2+z32**2)
          rl(3)=dsqrt(x13**2+y13**2+z13**2)
C  triangle in space : we compute the length of one side and the distance of the
C  opposing node to that side to compute the area
      rx= x21
      ry= y21
      rz= z21
      bx= x32
      by= y32
      bz= z32
      rlr= dsqrt( rx*rx + ry*ry + rz*rz )
      rlb= dsqrt( bx*bx + by*by + bz*bz )
      bpr= dsqrt((rx * bx + ry * by + rz *bz )**2)/rlr
      area= rlr*(dsqrt(rlb**2-bpr**2))/2.0d+00
      rmas= dens*area*esp/3.0
C     
C constructing mass matrix
      do 20 i=1,3
        i2=6+i
        i3=12+i
        fv(i,i)=rmas
        fv(i2,i2)=rmas
        fv(i3,i3)=rmas
  20  continue
      do 30 i=1,3
        i1=i+3
        i2=i+9
        i3=i+15
        fv(i1,i1)=(rl(1)**2+rl(3)**2)/420.*rmas
        fv(i2,i2)=(rl(2)**2+rl(1)**2)/420.*rmas
        fv(i3,i3)=(rl(3)**2+rl(2)**2)/420.*rmas
  30  continue
C
C..... COMPUTE THE BODY FORCE DUE TO GRAVITY IF NEEDED
C
      if (grvflg) then
        coef = 3.0d0*rmas
        grvfor(1) = coef*gamma(1)
        grvfor(2) = coef*gamma(2)
        grvfor(3) = coef*gamma(3)
      endif

C
C.... ACCUMULATE THE SUBDOMAIN MASS
C
        if (masflg) then
          totmas = totmas + 3.0d0*rmas
        endif
C
      return
      end
