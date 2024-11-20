C
C     LINEGAUSSQ  resturns the absissae  and weight factors of
C     the p-th Gauss-Legendre integration rule over the
C     segment (-1,+1).
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
C       WEIGHT    Weight factor
C
C=END USAGE
C=BLOCK FORTRAN
      subroutine    lgauss(p, i, xi, weight)
C
C                   A R G U M E N T S
C
      integer           p, i
      double precision  xi, weight
C
C                   L O C A L   V A R I A B L E S
C
C
C                   L O G I C
C
      if (p .le. 1)            then
        xi =      0.0
        weight =  2.0
      else if (p .eq. 2)       then
        if (i .eq. 1)          then
          xi =     -1.0/dsqrt(3.D0)
        else
          xi =      1.0/dsqrt(3.D0)
        end if
        weight =   1.0
      else if (p .eq. 3)       then
        if (i .eq. 1)          then
          xi =     -dsqrt(0.6D0)
          weight =  5.D0/9.
        else if (i .eq. 2)     then
          xi =      0.0
          weight =  8.D0/9.
        else
          xi =      dsqrt(0.6D0)
          weight =  5.D0/9.
        end if
      else
        if (i .eq. 1)          then
          xi =     -0.861136311594053D0
          weight =  0.347854845137454D0
        else if (i .eq. 2)     then
          xi =     -0.339981043584856D0
          weight =  0.652145154862546D0
        else if (i .eq. 3)     then
          xi =      0.339981043584856D0
          weight =  0.652145154862546D0
        else
          xi =      0.861136311594053D0
          weight =  0.347854845137454D0
        end if
      end if
      return
      end
