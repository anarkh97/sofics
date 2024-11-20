      SUBROUTINE Eg3C2C(NNO,NPO,X,Y,Z,NPI,IJT,POIDS,VP1,
     +                  VDPQ2,VDPQ1,ro,
     +                  DELTA,DF,DFINV,XINT,YINT,ZINT,A2,
     +                  GAMMA,GRVFOR)
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  COMPUTE GRAVITY FORCE VECTOR
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  PARAMETRES D'ENTREE :
C  -------------------
C  NNO         : NOMBRE DE NOEUDS DE L ELEMENT
C  NPO         : NOMBRE DE POINTS DE L ELEMENT
C  X,Y,Z       : COORDONNEES DES NOEUDS DE L ELEMENT
C  NPI         : NOMBRE DE POINTS  D INTEGRATION  SUR LE VOLUME
C  IJT         : PERMUTATION
C  POIDS       : POIDS DE LA FORMULE D INTEGRATION SUR LE CUBE
C  VP1         : VALEURS DES  DES POL DE BASE
C  NBAR        : NUMERO DU BARYCENTRE
C  VDPQ2       : VALEURS DES DERIVES DES POL DE BASE DE Q'2(R3)
C                AUX POINTS D INTEGRATION
C  VDPQ1       : VALEURS DES DERIVES DES POLYNOMES DE BASE DE Q1(R3)
C                AUX POINTS D INTEGRATION
C  ELAS        : MATRICE DE L ELASTICITE
C  DELT
C  DF
C  DFINV
C  XINT,YINT,ZINT
C  A2
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     PROGRAMMEUR  : MARINA VIDRASCU INRIA 1989
C  ..................................................................
      real*8 gamma(*),grvfor(*)
      real*8 XINT(NPI),YINT(NPI),ZINT(NPI)
      REAL*8 X(NPO),Y(NPO),Z(NPO)
      real*8 POIDS(NPI),VP1(NNO,NPI),
     +           VDPQ1(3,NPO,NPI),VDPQ2(3,NNO,NPI),DELTA(NPI),
     +           A2(1),DFINV(3,3),DF(3,3) ,zero,ro,lforce(NNO)
c     logical grvflg,masflg
      parameter (zero = 0.)
      INTEGER IJT( *),npi,nno,npo,i,indice,l
c
      INDICE = 0
      CALL FOBASE(3,3,NNO,NPO,NPI,VP1,VP1,VDPQ2,VDPQ1,X,Y,Z,
     +            A2,XINT,YINT,ZINT,DELTA,DFINV,DF,INDICE)


        DO 3 I = 1,NNO
          lforce(i) = 0.0
 3      CONTINUE

        DO 1 L  = 1,NPI
          DO 2 I = 1,NNO
            lforce(i) = lforce(i) + POIDS(L)*delta(l)*VP1(i,l)
 2        CONTINUE
 1      CONTINUE

C..... COMPUTE THE BODY FORCE DUE TO GRAVITY IF NEEDED
C
        DO 4 I = 1,NNO
          grvfor((3*i)-2) = lforce(i)*ro*gamma(1)
          grvfor((3*i)-1) = lforce(i)*ro*gamma(2)
          grvfor((3*i)  ) = lforce(i)*ro*gamma(3)
 4      CONTINUE
C

       END
