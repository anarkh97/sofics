        subroutine vol17(p,x,y,z,volume)
C
        integer p
        real*8  volume
        real*8  z(*),x(*),y(*)

        integer k,l,jj
        real*8 v,xi,eta,emu,weight,det
        real*8 q(8),qx(8),qy(8),qz(8)
C
C.... DETERMINE THE VOLUME OF THE BRICK 
C
        v = 0.0d0
        do 10 k=1, p
          do 20 l=1, p
            do 30 jj = 1, p
              call hxgaus(p,k,p,l,p,jj,xi,eta,emu,weight)
              call h8shpe(xi,eta,emu,x,y,z,q,qx,qy,qz,det)
              v = v + det
30          continue
20        continue
10      continue

        volume = v

        return
        end
