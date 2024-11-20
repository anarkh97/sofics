C=DECK LINEGAUSSQ
C=PURPOSE Get s-point abscissas and weight for quad product Gauss rule
C=AUTHOR C. A. Felippa, April 1969
C=VERSION May 1982 (Fortran 77)
C=EQUIPMENT Machine independent
C=KEYWORDS line Gauss integration rule absissae weigth
C=BLOCK ABSTRACT
C
C     LINEGAUSSQ  returns the absissae and weight factors
C     of the p-th Gauss-Legendre integration rule (=1,2,3,4) over
C     the interval xi:(-1,+1).
C
C=END ABSTRACT
C=BLOCK USAGE
C
C     The calling sequence is
C
C       call      LINEGAUSSQ  (P, I, XI, WEIGHT)
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
      subroutine    LINEGAUSSQ
     $             (p, i, xi, weight)
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
          xi =     -1.0/sqrt(3.D0)
        else
          xi =      1.0/sqrt(3.D0)
        end if
        weight =   1.0
      else if (p .eq. 3)       then
        if (i .eq. 1)          then
          xi =     -sqrt(0.6D0)
          weight =  5.D0/9.
        else if (i .eq. 2)     then
          xi =      0.0
          weight =  8.D0/9.
        else
          xi =      sqrt(0.6D0)
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
C=END FORTRAN
