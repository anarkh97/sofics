      SUBROUTINE Em3C2C(NNO,NPO,X,Y,Z,NPI,IJT,POIDS,VP1,
     +                  VDPQ2,VDPQ1,ro,AE,
     +                  DELTA,DF,DFINV,XINT,YINT,ZINT,A2,
     +                  GAMMA,GRVFOR,GRVFLG,TOTMAS,MASFLG)
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  BUT :   CALCUL DE LA MATRICE ELEMENTAIRE DE MASSE
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
C  PARAMETRES RESULTATS :
C  --------------------
C  AE          : MATRICE ELEMENTAIRE DE MASSE 
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     PROGRAMMEUR  : MARINA VIDRASCU INRIA 1989
C  ..................................................................
      real*8 gamma(*),grvfor(*),totmas
      real*8 XINT(NPI),YINT(NPI),ZINT(NPI)
      REAL*8 X(NPO),Y(NPO),Z(NPO)
      real*8 POIDS(NPI),AE(*),VP1(NNO,NPI),
     +           VDPQ1(3,NPO,NPI),VDPQ2(3,NNO,NPI),DELTA(NPI),
     +           A2(1),DFINV(3,3),DF(3,3) ,zero,ro,volume
      logical grvflg,masflg
      parameter (zero = 0.)
      INTEGER IJT( *),npi,nno,npo,i,indice,l
c
      INDICE = 0
      CALL FOBASE(3,3,NNO,NPO,NPI,VP1,VP1,VDPQ2,VDPQ1,X,Y,Z,
     +            A2,XINT,YINT,ZINT,DELTA,DFINV,DF,INDICE)



C       calcul du volume
C       -------------------
        volume = 0.0
        DO 1 L  = 1,NPI
              volume = volume+POIDS(L)*delta(l)
 1      CONTINUE
c
        DO 9 I = 1,3*NNO
         ae(i) = volume *RO/NNO
 9     CONTINUE
C
C..... COMPUTE THE BODY FORCE DUE TO GRAVITY IF NEEDED
C
        if (grvflg) then
          grvfor(1) = volume*ro*gamma(1)
          grvfor(2) = volume*ro*gamma(2)
          grvfor(3) = volume*ro*gamma(3)
        endif
C
C.... ACCUMULATE THE TOTAL MASS
C
        if (masflg) then
          totmas = totmas + volume*ro
        endif

       END
