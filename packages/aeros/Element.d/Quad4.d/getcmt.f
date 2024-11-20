      subroutine    GETCMT (rip, e, nu, c)
c
c     rip = 0 : plane stress
c     rip = 1 : plane strain
c
      implicit none

      double precision  rip, e, nu, c(3,3)
      double precision  omn,om2n,opn

      if (rip.eq.0) then
         c(1,1) =  e / (1.0d0-nu*nu)
         c(2,2) =  c(1,1)
         c(3,3) =  0.5d0*c(1,1)*(1.0d0-nu)
         c(1,2) =  c(1,1)*nu
         c(2,1) =  c(1,2)
         c(1,3) =  0.0d0
         c(2,3) =  0.0d0
         c(3,1) =  0.0d0
         c(3,2) =  0.0d0
         return
      endif

      if (rip.eq.1.0) then

         omn     = 1.0d0-nu
         om2n    = 1.0d0-2.0d0*nu
         opn     = 1.0d0+nu
      
         c(1,1) =  e*omn/opn/om2n
         c(2,2) =  c(1,1)
         c(3,3) =  e/opn/2.0d0
         c(1,2) =  e*nu /opn/om2n
         c(2,1) =  c(1,2)
         c(1,3) =  0.0d0
         c(2,3) =  0.0d0
         c(3,1) =  0.0d0
         c(3,2) =  0.0d0

         return
      endif

      write(*,*) " *** ERROR: wrong option in getcmt"
      stop

      end
