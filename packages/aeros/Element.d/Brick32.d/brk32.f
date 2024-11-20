C***********************************************************************
      SUBROUTINE BRK32(FUN, IFUN, DER, IDER, JDER, XI, ETA, ZETA,
     *     ITEST)
C-----------------------------------------------------------------------
C PURPOSE
C      CALCULATES THE VALUES OF THE SHAPE FUNCTIONS AND THEIR
C      DERIVATIVES AT A POINT FOR A 32-NODED BRICK ELEMENT.
C      FUNCTION CONTINUOUS ACROSS ELEMENT BOUNDARIES
C
C HISTORY
C
C      COPYRIGHT (C) 1997 : CLRC, RUTHERFORD APPLETON LABORATORY
C                           CHILTON, DIDCOT, OXFORDSHIRE OX11 0QX
C
C      RELEASE 1.1  29 OCT 1979 (CG)
C      COMMENTED    10 OCT 1980 (KR)
C
C ARGUMENTS IN
C      IFUN    DIMENSION OF VECTOR FUN (.GE.32)
C      IDER    FIRST DIMENSION OF ARRAY DER (.GE.3)
C      JDER    SECOND DIMENSION OF ARRAY DER (.GE.32)
C      XI      FIRST LOCAL COORDINATE
C      ETA     SECOND LOCAL COORDINATE
C      ZETA    THIRD LOCAL COORDINATE
C      ITEST   ERROR CHECKING OPTION
C
C ARGUMENTS OUT
C      FUN     VECTOR CONTAINING SHAPE FUNCTIONS.  FUN(I)
C              CONTAINS THE VALUE OF THE I'TH SHAPE FUNCTION AT
C              (XI,ETA,ZETA)
C      DER     ARRAY CONTAINING THE DERIVATIVES OF THE SHAPE
C              FUNCTIONS.  DER(I,J) CONTAINS THE VALUE OF THE
C              DERIVATIVE OF THE J'TH SHAPE FUNCTION WITH
C              RESPECT TO THE I'TH COORDINATE AT (XI,ETA,ZETA)
C
C ROUTINES CALLED
C      ERRMES
C
C
C     SUBROUTINE BRK32(FUN, IFUN, DER, IDER, JDER, XI, ETA, ZETA,
C    *     ITEST)
C***********************************************************************
C
      INTEGER IDER, IFUN, ITEST, JDER
C      INTEGER ERRMES, IDER, IERROR, IFUN, ITEST, JDER
      DOUBLE PRECISION DER, DNCX, DNCY, DNCZ, DNSXX, DNSXY,
     *     DNSXZ, DNSYX, DNSYY, DNSYZ, DNSZX, DNSZY, DNSZZ,
     *     ETA, FUN, M1, MTHRD, NC, NSX, NSY, NSZ, P1,
     *     PTHRD, X, XI, Y, Z, ZETA
C     *     PTHRD, SRNAME, VAL, VEPS, X, XI, Y, Z, ZETA, DUMMY
      DIMENSION DER(IDER,JDER), FUN(IFUN)
C      DATA SRNAME /8H BRK32  /
C
C     STATEMENT FUNCTIONS
C
      NC(X,Y,Z) = 1.D0/64.D0*(1.D0+XI*X)*(1.D0+ETA*Y)*(1.D0+ZETA*Z)
     *     *(9.D0*(XI*XI+ETA*ETA+ZETA*ZETA)-19.D0)
      NSX(X,Y,Z) = 9.D0/64.D0*(1.D0-XI*XI)*(1.D0+9.D0*XI*X)*(1.D0+ETA*
     *     Y)*(1.D0+ZETA*Z)
      NSY(X,Y,Z) = 9.D0/64.D0*(1.D0-ETA*ETA)*(1.D0+9.D0*ETA*Y)*
     *     (1.D0+XI*X)*(1.D0+ZETA*Z)
      NSZ(X,Y,Z) = 9.D0/64.D0*(1.D0-ZETA*ZETA)*(1.D0+9.D0*ZETA*Z)*
     *     (1.D0+XI*X)*(1.D0+ETA*Y)
      DNCX(X,Y,Z) = 1.D0/64.D0*(1.D0+ETA*Y)*(1.D0+ZETA*Z)*(X*(9.D0*
     *     (XI*XI+ETA*ETA+ZETA*ZETA)-19.D0)+18.D0*XI*(1.D0+XI*X))
      DNCY(X,Y,Z) = 1.D0/64.D0*(1.D0+XI*X)*(1.D0+ZETA*Z)*(Y*(9.D0*
     *     (XI*XI+ETA*ETA+ZETA*ZETA)-19.D0)+18.D0*ETA*(1.D0+ETA*Y))
      DNCZ(X,Y,Z) = 1.D0/64.D0*(1.D0+XI*X)*(1.D0+ETA*Y)*(Z*(9.D0*
     *     (XI*XI+ETA*ETA+ZETA*ZETA)-19.D0)+18.D0*ZETA*(1.D0+ZETA*
     *     Z))
      DNSXX(X,Y,Z) = 9.D0/64.D0*(1.D0+ETA*Y)*(1.D0+ZETA*Z)*(9.D0*X*
     *     (1.D0-XI*XI)-2.D0*XI*(1.D0+9.D0*XI*X))
      DNSXY(X,Y,Z) = 9.D0/64.D0*Y*(1.D0-XI*XI)*(1.D0+9.D0*XI*X)*
     *     (1.D0+ZETA*Z)
      DNSXZ(X,Y,Z) = 9.D0/64.D0*Z*(1.D0-XI*XI)*(1.D0+9.D0*XI*X)*
     *     (1.D0+ETA*Y)
      DNSYX(X,Y,Z) = 9.D0/64.D0*X*(1.D0-ETA*ETA)*(1.D0+9.D0*ETA*Y)*
     *     (1.D0+ZETA*Z)
      DNSYY(X,Y,Z) = 9.D0/64.D0*(1.D0+XI*X)*(1.D0+ZETA*Z)*(9.D0*Y*
     *     (1.D0-ETA*ETA)-2.D0*ETA*(1.D0+9.D0*ETA*Y))
      DNSYZ(X,Y,Z) = 9.D0/64.D0*Z*(1.D0-ETA*ETA)*(1.D0+9.D0*ETA*Y)*
     *     (1.D0+XI*X)
      DNSZX(X,Y,Z) = 9.D0/64.D0*X*(1.D0-ZETA*ZETA)*(1.D0+9.D0*ZETA*
     *     Z)*(1.D0+ETA*Y)
      DNSZY(X,Y,Z) = 9.D0/64.D0*Y*(1.D0-ZETA*ZETA)*(1.D0+9.D0*ZETA*
     *     Z)*(1.D0+XI*X)
      DNSZZ(X,Y,Z) = 9.D0/64.D0*(1.D0+XI*X)*(1.D0+ETA*Y)*(9.D0*Z*
     *     (1.D0-ZETA*ZETA)-2.D0*ZETA*(1.D0+9.D0*ZETA*Z))
C
C     PARAMETER CHECKING
C
C      IF (ITEST.EQ.-1) GO TO 1010
C      IERROR = 0
C      IF (IFUN.LT.32) IERROR = 1
C      IF (IDER.LT.3 .OR. JDER.LT.32) IERROR = 2
C      VAL = 1.0D0 + VEPS(DUMMY)
C      IF (DABS(XI).GT.VAL .OR. DABS(ETA).GT.VAL .OR.
C     *     DABS(ZETA).GT.VAL) IERROR = 3
C      ITEST = ERRMES(ITEST,IERROR,SRNAME)
C      IF (ITEST.NE.0) RETURN
C
C     BODY OF CODE
C
 1010 P1 = 1.D0
      M1 = -1.D0
      PTHRD = 1.D0/3.D0
      MTHRD = -1.D0/3.D0
C
C     SET SHAPE FUNCTIONS
C
      FUN(1) = NC(M1,M1,M1)
      FUN(2) = NSY(M1,MTHRD,M1)
      FUN(3) = NSY(M1,PTHRD,M1)
      FUN(4) = NC(M1,P1,M1)
      FUN(5) = NSX(MTHRD,P1,M1)
      FUN(6) = NSX(PTHRD,P1,M1)
      FUN(7) = NC(P1,P1,M1)
      FUN(8) = NSY(P1,PTHRD,M1)
      FUN(9) = NSY(P1,MTHRD,M1)
      FUN(10) = NC(P1,M1,M1)
      FUN(11) = NSX(PTHRD,M1,M1)
      FUN(12) = NSX(MTHRD,M1,M1)
      FUN(13) = NSZ(M1,M1,MTHRD)
      FUN(14) = NSZ(M1,P1,MTHRD)
      FUN(15) = NSZ(P1,P1,MTHRD)
      FUN(16) = NSZ(P1,M1,MTHRD)
      FUN(17) = NSZ(M1,M1,PTHRD)
      FUN(18) = NSZ(M1,P1,PTHRD)
      FUN(19) = NSZ(P1,P1,PTHRD)
      FUN(20) = NSZ(P1,M1,PTHRD)
      FUN(21) = NC(M1,M1,P1)
      FUN(22) = NSY(M1,MTHRD,P1)
      FUN(23) = NSY(M1,PTHRD,P1)
      FUN(24) = NC(M1,P1,P1)
      FUN(25) = NSX(MTHRD,P1,P1)
      FUN(26) = NSX(PTHRD,P1,P1)
      FUN(27) = NC(P1,P1,P1)
      FUN(28) = NSY(P1,PTHRD,P1)
      FUN(29) = NSY(P1,MTHRD,P1)
      FUN(30) = NC(P1,M1,P1)
      FUN(31) = NSX(PTHRD,M1,P1)
      FUN(32) = NSX(MTHRD,M1,P1)
C
C     SET DERIVATIVES
C
      DER(1,1) = DNCX(M1,M1,M1)
      DER(2,1) = DNCY(M1,M1,M1)
      DER(3,1) = DNCZ(M1,M1,M1)
      DER(1,2) = DNSYX(M1,MTHRD,M1)
      DER(2,2) = DNSYY(M1,MTHRD,M1)
      DER(3,2) = DNSYZ(M1,MTHRD,M1)
      DER(1,3) = DNSYX(M1,PTHRD,M1)
      DER(2,3) = DNSYY(M1,PTHRD,M1)
      DER(3,3) = DNSYZ(M1,PTHRD,M1)
      DER(1,4) = DNCX(M1,P1,M1)
      DER(2,4) = DNCY(M1,P1,M1)
      DER(3,4) = DNCZ(M1,P1,M1)
      DER(1,5) = DNSXX(MTHRD,P1,M1)
      DER(2,5) = DNSXY(MTHRD,P1,M1)
      DER(3,5) = DNSXZ(MTHRD,P1,M1)
      DER(1,6) = DNSXX(PTHRD,P1,M1)
      DER(2,6) = DNSXY(PTHRD,P1,M1)
      DER(3,6) = DNSXZ(PTHRD,P1,M1)
      DER(1,7) = DNCX(P1,P1,M1)
      DER(2,7) = DNCY(P1,P1,M1)
      DER(3,7) = DNCZ(P1,P1,M1)
      DER(1,8) = DNSYX(P1,PTHRD,M1)
      DER(2,8) = DNSYY(P1,PTHRD,M1)
      DER(3,8) = DNSYZ(P1,PTHRD,M1)
      DER(1,9) = DNSYX(P1,MTHRD,M1)
      DER(2,9) = DNSYY(P1,MTHRD,M1)
      DER(3,9) = DNSYZ(P1,MTHRD,M1)
      DER(1,10) = DNCX(P1,M1,M1)
      DER(2,10) = DNCY(P1,M1,M1)
      DER(3,10) = DNCZ(P1,M1,M1)
      DER(1,11) = DNSXX(PTHRD,M1,M1)
      DER(2,11) = DNSXY(PTHRD,M1,M1)
      DER(3,11) = DNSXZ(PTHRD,M1,M1)
      DER(1,12) = DNSXX(MTHRD,M1,M1)
      DER(2,12) = DNSXY(MTHRD,M1,M1)
      DER(3,12) = DNSXZ(MTHRD,M1,M1)
      DER(1,13) = DNSZX(M1,M1,MTHRD)
      DER(2,13) = DNSZY(M1,M1,MTHRD)
      DER(3,13) = DNSZZ(M1,M1,MTHRD)
      DER(1,14) = DNSZX(M1,P1,MTHRD)
      DER(2,14) = DNSZY(M1,P1,MTHRD)
      DER(3,14) = DNSZZ(M1,P1,MTHRD)
      DER(1,15) = DNSZX(P1,P1,MTHRD)
      DER(2,15) = DNSZY(P1,P1,MTHRD)
      DER(3,15) = DNSZZ(P1,P1,MTHRD)
      DER(1,16) = DNSZX(P1,M1,MTHRD)
      DER(2,16) = DNSZY(P1,M1,MTHRD)
      DER(3,16) = DNSZZ(P1,M1,MTHRD)
      DER(1,17) = DNSZX(M1,M1,PTHRD)
      DER(2,17) = DNSZY(M1,M1,PTHRD)
      DER(3,17) = DNSZZ(M1,M1,PTHRD)
      DER(1,18) = DNSZX(M1,P1,PTHRD)
      DER(2,18) = DNSZY(M1,P1,PTHRD)
      DER(3,18) = DNSZZ(M1,P1,PTHRD)
      DER(1,19) = DNSZX(P1,P1,PTHRD)
      DER(2,19) = DNSZY(P1,P1,PTHRD)
      DER(3,19) = DNSZZ(P1,P1,PTHRD)
      DER(1,20) = DNSZX(P1,M1,PTHRD)
      DER(2,20) = DNSZY(P1,M1,PTHRD)
      DER(3,20) = DNSZZ(P1,M1,PTHRD)
      DER(1,21) = DNCX(M1,M1,P1)
      DER(2,21) = DNCY(M1,M1,P1)
      DER(3,21) = DNCZ(M1,M1,P1)
      DER(1,22) = DNSYX(M1,MTHRD,P1)
      DER(2,22) = DNSYY(M1,MTHRD,P1)
      DER(3,22) = DNSYZ(M1,MTHRD,P1)
      DER(1,23) = DNSYX(M1,PTHRD,P1)
      DER(2,23) = DNSYY(M1,PTHRD,P1)
      DER(3,23) = DNSYZ(M1,PTHRD,P1)
      DER(1,24) = DNCX(M1,P1,P1)
      DER(2,24) = DNCY(M1,P1,P1)
      DER(3,24) = DNCZ(M1,P1,P1)
      DER(1,25) = DNSXX(MTHRD,P1,P1)
      DER(2,25) = DNSXY(MTHRD,P1,P1)
      DER(3,25) = DNSXZ(MTHRD,P1,P1)
      DER(1,26) = DNSXX(PTHRD,P1,P1)
      DER(2,26) = DNSXY(PTHRD,P1,P1)
      DER(3,26) = DNSXZ(PTHRD,P1,P1)
      DER(1,27) = DNCX(P1,P1,P1)
      DER(2,27) = DNCY(P1,P1,P1)
      DER(3,27) = DNCZ(P1,P1,P1)
      DER(1,28) = DNSYX(P1,PTHRD,P1)
      DER(2,28) = DNSYY(P1,PTHRD,P1)
      DER(3,28) = DNSYZ(P1,PTHRD,P1)
      DER(1,29) = DNSYX(P1,MTHRD,P1)
      DER(2,29) = DNSYY(P1,MTHRD,P1)
      DER(3,29) = DNSYZ(P1,MTHRD,P1)
      DER(1,30) = DNCX(P1,M1,P1)
      DER(2,30) = DNCY(P1,M1,P1)
      DER(3,30) = DNCZ(P1,M1,P1)
      DER(1,31) = DNSXX(PTHRD,M1,P1)
      DER(2,31) = DNSXY(PTHRD,M1,P1)
      DER(3,31) = DNSXZ(PTHRD,M1,P1)
      DER(1,32) = DNSXX(MTHRD,M1,P1)
      DER(2,32) = DNSXY(MTHRD,M1,P1)
      DER(3,32) = DNSXZ(MTHRD,M1,P1)
C
      RETURN
      END
C***********************************************************************