        subroutine vol20(p,x,y,z,volume)
C
C       implicit none
        integer p
        real*8  volume
        real*8  z(*),x(*),y(*)

        integer k,l,jj
        real*8 v,xi,eta,emu,weight,det
        real*8 q(20),qx(20),qy(20),qz(20)
C
C.... DETERMINE THE VOLUME OF THE BRICK 
C
        v = 0.0d0
        do 10 k=1, p
          do 20 l=1, p
            do 30 jj = 1, p
              call hxgaus20(p,k,p,l,p,jj,xi,eta,emu,weight)
              call h20shpe(xi,eta,emu,x,y,z,q,qx,qy,qz,det)
              v = v + det
30          continue
20        continue
10      continue

        volume = v

        return
        end
