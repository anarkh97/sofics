C--------------------------------------------------------------------
C     HXGAUS20 returns the hexahedron coordinates of
C     sample points and weights for a Gauss-product integration
C     rule over an isoparametric hexahedron.
C
C     Input arguments:
C
C       P1        Number of Gauss points in the XI direction
C       I1        Index of sample point in the XI direction
C       P2        Number of Gauss points in the ETA direction
C       I2        Index of sample point in the ETA direction
C       P3        Number of Gauss points in the MU direction
C       I3        Index of sample points in the MU direction
C        
C     Outputs arguments:
C
C       XI, ETA, MU   Hexahedral coordinates of sample point
C       WEIGHT    Weight factor
C--------------------------------------------------------------------
      subroutine    hxgaus20(p1,i1,p2,i2,p3,i3,xi,eta,mu,weight)
C
C                   A R G U M E N T S
C
      integer           p1, i1, p2, i2, p3, i3
      double precision  xi, eta, mu, weight
C
C                   L O C A L   V A R I A B L E S
C
      double precision  w1, w2, w3
C
C                   L O G I C
C
      call LINEGAUSSQ (p1, i1, xi,  w1)
      call LINEGAUSSQ (p2, i2, eta, w2)
      call LINEGAUSSQ (p3, i3, mu,  w3)

      weight =  w1 * w2 * w3

      return
      end
C
C-------------------------------------------------------------------
C
C     LINEGAUSSQ  returns the absissae and weight factors of
C     the p-th Gauss-Legendre integration rule (p=1,2,3,4) for use
C     over the interval xi:(-1,+1).
C
C     Input arguments:
C
C       P         Number of points in the integration rule (1 to 4)
C                 If le 1 assume 1; if gt 4 assume 4.
C       I         Index of sample point (1 to P)
C
C     Outputs arguments:
C
C       XI        Absissa of sample point (zero of Legendre polynomial)
C       WEIGHT    Weight factor
C--------------------------------------------------------------------
      subroutine    LINEGAUSSQ(p, i, xi, weight)
C
C                   A R G U M E N T S
C
      integer           p, i
      double precision  xi, weight
C
      if (p .le. 1)            then
        xi =       0.0d0
        weight =   2.0d0
        return
      else if (p .eq. 2)       then
        if (i .eq. 1)          then
          xi =     -1.0d0/sqrt(3.d0)
        else
          xi =      1.0d0/sqrt(3.d0)
        end if
        weight =   1.0d0
        return
      else if (p .eq. 3)       then
        if (i .eq. 1)          then
          xi =     -sqrt(0.6d0)
          weight =  5.d0/9.0d0
        else if (i .eq. 2)     then
          xi =      0.0d0
          weight =  8.d0/9.0d0
        else
          xi =      sqrt(0.6d0)
          weight =  5.d0/9.0d0
        end if
        return
      else
        if (i .eq. 1)          then
          xi =     -0.861136311594053d0
          weight =  0.347854845137454d0
        else if (i .eq. 2)     then
          xi =     -0.339981043584856d0
          weight =  0.652145154862546d0
        else if (i .eq. 3)     then
          xi =      0.339981043584856d0
          weight =  0.652145154862546d0
        else
          xi =      0.861136311594053d0
          weight =  0.347854845137454d0
        end if
      end if
      return
      end
