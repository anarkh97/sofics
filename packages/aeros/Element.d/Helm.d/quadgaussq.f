C=DECK QUADGAUSSQ
C=PURPOSE Get s-point abscissas and weight for quad product Gauss rule
C=AUTHOR C. A. Felippa, April 1969
C=VERSION May 1982 (Fortran 77)
C=EQUIPMENT Machine independent
C=KEYWORDS quadrilateral Gauss integration rule absissae weigth
C=BLOCK ABSTRACT
C
C     QUADGAUSSQ returns the quadriliteral coordinates of
C     sample points and weights for a Gauss-product integration
C     rule over an isoparametric quadrilateral.
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       call      QUADGAUSSQ  (P1, I1, P2, I2, XI, ETA, WEIGHT)
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
      subroutine    QUADGAUSSQ
     $             (p1, i1, p2, i2, xi, eta, weight)
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
      call      LINEGAUSSQ (p1, i1, xi, w1)
      call      LINEGAUSSQ (p2, i2, eta, w2)
      weight =  w1 * w2
      return
      end
C=END FORTRAN
