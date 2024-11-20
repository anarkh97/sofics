C
C     TGAUSS resturns the absissae  and weight factors of
C     the p-th Gauss-Legendre integration rule over isoparametric triangle
C
C
C     The calling sequence is
C
C       CALL      lgauss  (p, i, xi, weight)
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
C       ETA        Absissa of sample point (zero of Legendre polynomial)
C       WEIGHT    Weight factor
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine    tgauss(p, i, xi, eta, weight)
C
C                   A R G U M E N T S
C
      integer           p, i
      double precision  xi, eta, weight
C
C                   L O C A L   V A R I A B L E S
C
C
C                   L O G I C
C
      if (p .le. 3) .or. (p.gt.3)           then
          p = 3
      endif

        if (i .eq. 1)          then
          xi  =     0.5D0
          eta =     0.0D0
          weight =  1.D0/3.
        else if (i .eq. 2)     then
          xi  =      0.0D0
          eta =      0.5D0
          weight =  1.D0/3.
        else
          xi  =      0.5D0
          eta =      0.5D0
          weight =  1.D0/3.
        end if

      return
      end
