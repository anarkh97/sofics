C     HXGAUS returns the hexahedron coordinates of
C     sample points and weights for a Gauss-product integration
C     rule over an isoparametric hexahedron.
C
C     Input arguments:
C
C       P1        Number of Gauss points in the XI direction
C       I1        Index of sample point in the XI direction
C       P2        Number of Gauss points in the ETA direction
C       I2        Index of sample point in the ETA direction
C       P3        Number of Gauss points in the EMU direction
C       I3        Index of sample points in the EMU direction
C        
C
C     Outputs arguments:
C
C       XI, ETA, EMU   Hexahedral coordinates of sample point
C       WEIGHT    Weight factor
C
      subroutine    hxgaus(p1,i1,p2,i2,p3,i3,xi,eta,emu,weight)
C
C                   A R G U M E N T S
C
      integer           p1, i1, p2, i2, p3, i3
      real*8  xi, eta, emu, weight
C
C                   L O C A L   V A R I A B L E S
C
      real*8  w1, w2, w3
C
C                   L O G I C
C
      call lgauss(p1, i1, xi,  w1)
      call lgauss(p2, i2, eta, w2)
      call lgauss(p3, i3, emu, w3)
      weight =  w1 * w2 * w3
      return
      end
