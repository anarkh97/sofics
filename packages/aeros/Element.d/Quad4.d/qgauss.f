C
C     QUADGAUSSQ returns the quadriliteral coordinates of
C     sample points and weights for a Gauss-product integration
C     rule over an isoparametric quadrilateral.
C
C
C     The calling sequence is
C
C       CALL      QGAUSS  (P1, I1, P2, I2, XI, ETA, WEIGHT)
C
C     Input arguments:
C
C       P1        Number of Gauss points in the XI direction
C       I1        Index of sample point in the XI direction
C       P2        Number of Gauss points in the ETA direction
C       I2        Index of sample point in the ETA direction
C
C     Outputs arguments:
C
C       XI, ETA   Quadrilateral coordinates of sample point
C       WEIGHT    Weight factor
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine    qgauss(p1, i1, p2, i2, xi, eta, weight)
C
C                   A R G U M E N T S
C
      integer           p1, i1, p2, i2
      double precision  xi, eta, weight
C
C                   L O C A L   V A R I A B L E S
C
      double precision  w1, w2
C
C                   L O G I C
C
      call      lgauss (p1, i1, xi, w1)
      call      lgauss (p2, i2, eta, w2)
      weight =  w1 * w2
      return
      end
C=END FORTRAN
 
