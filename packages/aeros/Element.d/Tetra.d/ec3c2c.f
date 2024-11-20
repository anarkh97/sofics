       SUBROUTINE EC3C2C(elm,NNO,NPO,X,Y,Z,NPI,IJT,VP1,
     +                   VDPQ2,VDPQ1,e,xnu,ELAS,
     +                   stress,SIGMAE,strain, straine,
     +                   DELTA,DF,DFINV,XINT,YINT,ZINT,A2,u,
     +                   maxgus,maxstr,msize,outerr,vmflg,
     +                   strainFlg)
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  BUT : CALCUL DES TABLEAUX SIGMAE ET SIGMAE
C  ---        CONTRAINTES ELEMENTAIRES D UN ELEMENT VOLUMIQUE
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C   PARAMETRES D'ENTREE :
C  -------------------
C  IJT          : PERMUTATION
C  NPI          : NOMBRE DE POINTS DE LA FORMULE D INTEGRATION
C                 SUR LE CUBE DE REFERENCE
C  X,Y,Z        : COORDONNEES DES NOEUDS DE L ELEMENT
C  ELAS         : MATRICE D ELASTICITE DE L ELEMENT
C  VDPQ2        : VALEURS DES DERIVES DES POL DE BASE DE Q'2(R3)
C  VDPQ1        : VALEURS DES DERIVES DES POLYNOMES DE BASE DE Q1(R3)
C  POLQ2        : VALEURS DES POLYNOMES DE BASE DE Q'2(R3)
C  POLQ1        : VALEURS DES POLYNOMES DE BASE DE Q1(R3)
C  POUR CES POL. OU DP. CES VALEURS  SONT AUX POINTS D INTEGRATION
C  PARAMETRES RESULTATS :
C  --------------------
C  SIGMAE       : VALEURS DES CONTRAINTES ELEMENTAIRES
C  SIGMA        : contraintes
C  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  PROGRAMMEURS  : Marina Vidrascu 1996
C  ..................................................................
      INTEGER maxgus,maxstr,msize,outerr,elm
      INTEGER indice,i,j,j1,l,NNO,NPI,npo
      logical vmflg,strainFlg
      real*8  stress(msize,maxstr,maxgus),u(*)
      real*8  strain(msize,maxstr,maxgus)
      real*8  straine(6,3*NNO,NPI)
      REAL*8 XINT(NPI),YINT(NPI),ZINT(NPI)
      REAL*8 X(NPO),Y(NPO),Z(NPO)
      REAL*8 VP1(NNO,NPI),
     +           VDPQ1(3,NPO,NPI),VDPQ2(3,NNO,NPI),DELTA(NPI),
     +           A2(3,NNO,NPI),ELAS(9),
     +           DFINV(3,3),DF(3,3),
     +           SIGMAE(6,3*NNO,NPI),
     +           e,xnu,xlam,xmu
      INTEGER IJT(*)
C
C   -----  CALCUL DE DELTA  -----
C
      INDICE = 1
      CALL FOBASE(3,3,NNO,NPO,NPI,VP1,VP1,VDPQ2,VDPQ1,X,Y,Z,
     +            A2,XINT,YINT,ZINT,DELTA,DFINV,DF,INDICE)
c
c
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
      DO 10 L=1,NPI
C
C        B = DELTA * (DF)-1 * VDPQ2(L)
C
      DO 11 J=1,NNO
C
C   --- BLOC1
C
        J1=IJT(J)
            straine(1,J1,L) = A2(1,J,L) / DELTA(L)
            straine(2,J1,L) = 0.0
            straine(3,J1,L) = 0.0
C           straine(2,J1,L) = A2(1,J,L) / DELTA(L)
C           straine(3,J1,L) = A2(1,J,L) / DELTA(L)
            straine(4,J1,L) = A2(2,J,L) / DELTA(L)
            straine(5,J1,L) = 0.0
            straine(6,J1,L) = A2(3,J,L) / DELTA(L)

            SIGMAE(1,J1,L) = ELAS(1)*A2(1,J,L) / DELTA(L)
            SIGMAE(2,J1,L) = ELAS(2)*A2(1,J,L) / DELTA(L)
            SIGMAE(3,J1,L) = ELAS(4)*A2(1,J,L) / DELTA(L)
            SIGMAE(4,J1,L) = ELAS(7)*A2(2,J,L) / DELTA(L)
            SIGMAE(5,J1,L) = 0
            SIGMAE(6,J1,L) = ELAS(9)*A2(3,J,L) / DELTA(L)
c
C   --- BLOC2
C
        J1=IJT(J+NNO)
            straine(1,J1,L) = 0.0
C           straine(1,J1,L) = A2(2,J,L) / DELTA(L)
            straine(2,J1,L) = A2(2,J,L) / DELTA(L)
            straine(3,J1,L) = 0.0
C           straine(3,J1,L) = A2(2,J,L) / DELTA(L)
            straine(4,J1,L) = A2(1,J,L) / DELTA(L)
            straine(5,J1,L) = A2(3,J,L) / DELTA(L)
            straine(6,J1,L) = 0.0

            SIGMAE(1,J1,L) = ELAS(2)*A2(2,J,L) / DELTA(L)
            SIGMAE(2,J1,L) = ELAS(3)*A2(2,J,L) / DELTA(L)
            SIGMAE(3,J1,L) = ELAS(5)*A2(2,J,L) / DELTA(L)
            SIGMAE(4,J1,L) = ELAS(7)*A2(1,J,L) / DELTA(L)
            SIGMAE(5,J1,L) = ELAS(8)*A2(3,J,L) / DELTA(L)
            SIGMAE(6,J1,L) = 0
C
C   --- BLOC3
C
        J1=IJT(J+NNO*2)
            straine(1,J1,L) = 0.0
            straine(2,J1,L) = 0.0
C           straine(1,J1,L) = A2(3,J,L) / DELTA(L)
C           straine(2,J1,L) = A2(3,J,L) / DELTA(L)
            straine(3,J1,L) = A2(3,J,L) / DELTA(L)
            straine(4,J1,L) = 0.0
            straine(5,J1,L) = A2(2,J,L) / DELTA(L)
            straine(6,J1,L) = A2(1,J,L) / DELTA(L)

            SIGMAE(1,J1,L) = ELAS(4)*A2(3,J,L) / DELTA(L)
            SIGMAE(2,J1,L) = ELAS(5)*A2(3,J,L) / DELTA(L)
            SIGMAE(3,J1,L) = ELAS(6)*A2(3,J,L) / DELTA(L)
            SIGMAE(4,J1,L) = 0
            SIGMAE(5,J1,L) = ELAS(8)*A2(2,J,L) / DELTA(L)
            SIGMAE(6,J1,L) = ELAS(9)*A2(1,J,L) / DELTA(L)
 11    CONTINUE
c
c       multiplier la contraine elementaire par le deplacement
c
        do 2 i = 1,6
           stress(elm,i,l) = 0.0
           strain(elm,i,l) = 0.0
           do 2 j = 1, 3*NNO
              stress(elm,i,l) = stress(elm,i,l) +  sigmae(i,j,l)*u(j)
              strain(elm,i,l) = strain(elm,i,l) + straine(i,j,l)*u(j)
 2      continue

   10 CONTINUE
C
C  les coordonnees des points de calcul sont xint,yint,zint
C

C
C.... COMPUTE THE VON MISES STRESS
C
        if (vmflg) then
          call vmelmv(stress,maxgus,maxstr,msize,elm,nno)
        endif

C
C.... COMPUTE THE VON MISES STRAIN
C
        if(strainFlg) then
          call strainvm(strain,maxgus,maxstr,msize,nno)
        endif

      END
