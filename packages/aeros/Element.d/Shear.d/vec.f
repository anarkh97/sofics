c
c some simple vector subroutines used to create the shear element
c
c***********************************************************************
      subroutine cross(u,v,w)
c     w = u X v
c
      real*8     u(3),v(3),w(3)
c
      w(1) = (u(2)*v(3)) - (u(3)*v(2))
      w(2) = (u(3)*v(1)) - (u(1)*v(3))
      w(3) = (u(1)*v(2)) - (u(2)*v(1))
c
      return
      end
c
c***********************************************************************
      subroutine unitv(u,v,w)
c     unit vector from point u to v
c
      real*8     u(3),v(3),w(3)
c
      do 10 i=1,3
        w(i) = v(i)-u(i)
 10   continue
      call normve(w)
c
      return
      end
c
c***********************************************************************
      subroutine normve(u)
c     make vector into unit vector
c
      real*8     u(3)
      real*8     len,zero
      data    zero  /0.000000D+00/
c
      len = zero
      do 10 i=1,3
        len  = len + u(i)*u(i)
 10   continue
      if (len.eq.zero) then
        write(*,*) 'Vector is of length zero'
        stop
      else
        len = abs(sqrt(len))
      endif
      do 20 i=1,3
        u(i) = u(i)/len
 20   continue
c
      return
      end
c
c***********************************************************************
      subroutine length(u,len)
c     length of u
c
      real*8     u(3),len
      real*8     zero
      data    zero  /0.000000D+00/
c
      len = zero
      do 10 i=1,3
        len  = len + u(i)*u(i)
 10   continue
      if (len.eq.zero) then
        write(*,*) 'Vector is of length zero'
        stop
      else
        len = abs(sqrt(len))
      endif
c
      return
      end
c
