      SUBROUTINE ER3C2C(NNO,NPO,X,Y,Z,NPI,IJT,POIDS,VP1,
     +                  VDPQ2,VDPQ1,e,xnu,ELAS,AE,
     +                  DELTA,DF,DFINV,XINT,YINT,ZINT,A2,status)
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  BUT :   CALCUL DE LA MATRICE ELEMENTAIRE DE RAIDEUR
C  ---     D UN ELEMENT  HEXA3Q2D OU HEXA3Q2C ISOTROPE OU ORTHOTROPE
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
C  AE          : MATRICE ELEMENTAIRE DE RAIDEUR STOCKEE SYMETRIQUE
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     PROGRAMMEUR  : MARINA VIDRASCU INRIA 1989
C  ..................................................................
      real*8 XINT(NPI),YINT(NPI),ZINT(NPI)
      real*8 X(NPO),Y(NPO),Z(NPO)
      real*8 POIDS(NPI),AE(*),VP1(NNO,NPI),
     +           VDPQ1(3,NPO,NPI),VDPQ2(3,NNO,NPI),DELTA(NPI),
     +           A2(3,NNO,NPI),ELAS(9),
     +           DFINV(3,3),DF(3,3) ,zero,e,xnu,xmu,xlam
      parameter (zero = 0.)
      integer nno,npi,npo,i,j,ii,jj,l,k1,k2,k3,nm,indice
      INTEGER IJT( *), status
c
      status = 0
      NM=NNO*3*(NNO*3+1)/2
      call dcopy(nm,zero,0,ae,1)
      INDICE = 1
      CALL FOBASE(3,3,NNO,NPO,NPI,VP1,VP1,VDPQ2,VDPQ1,X,Y,Z,
     +            A2,XINT,YINT,ZINT,DELTA,DFINV,DF,INDICE)
      do i=1,NPI
        if(delta(i).eq.0) then
           write(6,*) 'Nul jacobian in ER3C2C '
           status = -1
           return
        endif
        if(delta(i).le.0) then
           write(6,*) 'Negative jacobian in ER3C2C '
           status = -1
           return
        endif
      enddo


C           --------------
               XMU   = E/(2.*(1.+XNU))
               XLAM  = E*XNU/((1.+XNU)* (1.-2.*XNU))
               ELAS(1) = XLAM + 2. * XMU
               ELAS(3) = XLAM + 2. * XMU
               ELAS(6) = XLAM + 2. * XMU
               ELAS(2) = XLAM
               ELAS(4) = XLAM
               ELAS(5) = XLAM
               ELAS(7) = XMU
               ELAS(8) = XMU
               ELAS(9) = XMU
C
c       pour traiter le nonisotrope donner directement elas!
C       LES BLOCS DIAGONAUX
C       -------------------
        DO 9 J = 1,NNO
        DO 9 I = 1,J
C
C          BLOC 11
C          -------
           II = MIN(IJT(I),IJT(J))
           JJ = MAX(IJT(I),IJT(J))
           K1 = JJ*(JJ-1)/2+II
C
C          BLOC 22
C          -------
           II = MIN(IJT(NNO+I),IJT(NNO+J))
           JJ = MAX(IJT(NNO+I),IJT(NNO+J))
           K2  = JJ*(JJ-1)/2+II
C
C          BLOC 33
C          -------
           II = MIN(IJT(2*NNO+I),IJT(2*NNO+J))
           JJ = MAX(IJT(2*NNO+I),IJT(2*NNO+J))
           K3  = JJ*(JJ-1)/2+II
           DO 1 L  = 1,NPI
              AE(K1) = AE(K1)+POIDS(L)/DELTA(L)*(
     +                    ELAS(1) * A2(1,I,L) * A2(1,J,L) +
     +                    ELAS(7) * A2(2,I,L) * A2(2,J,L) +
     +                    ELAS(9) * A2(3,I,L) * A2(3,J,L) )
              AE(K2) = AE(K2)+POIDS(L)/DELTA(L)*(
     +                    ELAS(7) * A2(1,I,L) * A2(1,J,L) +
     +                    ELAS(3) * A2(2,I,L) * A2(2,J,L) +
     +                    ELAS(8) * A2(3,I,L) * A2(3,J,L) )
              AE(K3) = AE(K3)+POIDS(L)/DELTA(L)*(
     +                    ELAS(9) * A2(1,I,L) * A2(1,J,L) +
     +                    ELAS(8) * A2(2,I,L) * A2(2,J,L) +
     +                    ELAS(6) * A2(3,I,L) * A2(3,J,L) )
 1         CONTINUE
 9     CONTINUE
C
C       LES BLOCS EXTRA-DIAGONAUX
C       -------------------------
        DO 10 J = 1,NNO
        DO 10 I = 1,NNO
C
C          BLOC 12
C          -------
           II = MIN(IJT(I),IJT(NNO+J))
           JJ = MAX(IJT(I),IJT(NNO+J))
           K1 = JJ*(JJ-1)/2+II
C
C          BLOC 13
C          -------
           II = MIN(IJT(I),IJT(2*NNO+J))
           JJ = MAX(IJT(I),IJT(2*NNO+J))
           K2  = JJ*(JJ-1)/2+II
C
C          BLOC 23
C          -------
           II = MIN(IJT(NNO+I),IJT(2*NNO+J))
           JJ = MAX(IJT(NNO+I),IJT(2*NNO+J))
           K3  = JJ*(JJ-1)/2+II
           DO 4 L  = 1,NPI
              AE(K1) = AE(K1)+POIDS(L)/DELTA(L)*(
     +                    ELAS(7) * A2(2,I,L) * A2(1,J,L) +
     +                    ELAS(2) * A2(1,I,L) * A2(2,J,L) )
              AE(K2) = AE(K2)+POIDS(L)/DELTA(L)*(
     +                    ELAS(9) * A2(3,I,L) * A2(1,J,L) +
     +                    ELAS(4) * A2(1,I,L) * A2(3,J,L) )
              AE(K3) = AE(K3)+POIDS(L)/DELTA(L)*(
     +                    ELAS(8) * A2(3,I,L) * A2(2,J,L) +
     +                    ELAS(5) * A2(2,I,L) * A2(3,J,L) )
  4        CONTINUE
 10     CONTINUE
       END
