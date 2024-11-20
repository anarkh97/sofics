      SUBROUTINE FOBASE(NDIM,NDIMS,NNO,NPO,NPI,VP,VP1,VDPQ2,VDPQ1,X,Y,Z,
     +                 VDPELC,XINT,YINT,ZINT,DELTA,DFINV,DF,INDIC)
C .....................................................................
C BUT : CALCULER  SUR L'ELEMENT COURANT, AU POINTS D'INTEGRATION
C ---      -  LES VALEURS DES GRADIENTS DES FONCTIONS DE BASE  (VDPELC)
C          -  LES VALEURS DES JACOBIENS                        (DELTA)
C          -  LES COORDONNEES DES POINTS D'INTEGRATION
C .....................................................................
C PROGRAMMEUR : MARINA VIDRASCU  INRIA 1989
C .....................................................................
C PARAMETRES D'ENTREE :
C -------------------
C  NDIM           :  DIMENSION
C  NDIMS          :  3 VOLUME
C                    2 SURFACE
C  NNO            :  NOMBRE DE NOEUDS
C  NPO            :  NOMBRE DE POINTS
C  NPI            :  NOMBRE DE POINTS D'INTEGRATION
C  VP             :  VALEUR DES POLYNOMES DE BASE (FONCTIONS)
C  VP1            :  VALEUR DES POLYNOMES DE BASE (GEOMETRIE - FT)
C  VDPQ2          :  DERIVEES DES POLYNOMES DE BASE (FONCTIONS)
C  VDPQ1          :  DERIVEES DES DES POLYNOMES DE BASE (GEOMETRIE - FT)
C  X,Y,Z          :  COORDONNEES DES POINTS
C  VDPELC         :  DERIVEES DES POLS DE BASE SUR ELEMENT  COURANT
C PARAMETRES DE SORTIE :
C --------------------
C  XINT,YINT,ZINT :  COORDONNEES DES POINTS D'INTEGRATION
C  DELTA          :  JACOBIEN
C  INDIC          :  1 CALCUL DES DERIVEES DES POLYNOMES
C                    2 CALCUL DES COORDONNEES DES POINTS D'INTEGRATION
C                    3 CALCUL DE DELTA
C ---------------------------------------------------------------------
      real*8     VP(NNO,NPI),VP1(NPO,NPI),VDPQ2(NDIM,NNO,NPI),
     +           VDPQ1(NDIMS,NPO,NPI),X(NPO),Y(NPO),
     +           Z(NPO),VDPELC(NDIM,NNO,NPI),XINT(NPI)
     +           ,YINT(NPI),ZINT(NPI),DELTA(NPI),
     +           DF(NDIMS,3),DFINV(3,3)
      integer    ndim,ndims,nno,npo,npi,indic,i,j,n,l
      DO 25 L=1,NPI
C
         DO 21 I = 1,NDIMS
         DO 21 J = 1,NDIM
            Df(I,J) = 0
  21     CONTINUE
C
         DO 1 I=1,NDIMS
         DO 1 N=1,NPO
            DF(I,1) = DF(I,1) + VDPQ1(I,N,L) * X(N)
            DF(I,2) = DF(I,2) + VDPQ1(I,N,L) * Y(N)
            DF(I,3) = DF(I,3) + VDPQ1(I,N,L) * Z(N)
   1     CONTINUE
         IF (NDIMS .EQ. 3) THEN
            DFINV(1,1) = DF(2,2)*DF(3,3)-DF(2,3)*DF(3,2)
            DFINV(2,2) = DF(1,1)*DF(3,3)-DF(3,1)*DF(1,3)
            DFINV(3,3) = DF(1,1)*DF(2,2)-DF(2,1)*DF(1,2)
            DFINV(1,2) = DF(1,3)*DF(3,2)-DF(1,2)*DF(3,3)
            DFINV(2,1) = DF(2,3)*DF(3,1)-DF(2,1)*DF(3,3)
            DFINV(1,3) = DF(1,2)*DF(2,3)-DF(1,3)*DF(2,2)
            DFINV(3,1) = DF(2,1)*DF(3,2)-DF(2,2)*DF(3,1)
            DFINV(2,3) = DF(2,1)*DF(1,3)-DF(2,3)*DF(1,1)
            DFINV(3,2) = DF(3,1)*DF(1,2)-DF(3,2)*DF(1,1)
C
            DELTA(L) = DF(1,1) * DFINV(1,1) +
     +                    DF(2,1) * DFINV(1,2) +
     +                    DF(3,1) * DFINV(1,3)
         ELSE
C
            DELTA(L) = SQRT(
     +                 (DF(1,1)*DF(2,2)-DF(1,2)*DF(2,1))**2
     +               + (DF(1,2)*DF(2,3)-DF(1,3)*DF(2,2))**2
     +               + (DF(1,3)*DF(2,1)-DF(1,1)*DF(2,3))**2)
         END IF
         IF (INDIC .GT. 0 .AND. INDIC .LE. 2) THEN
C
          DO 2 I = 1,3
          DO 2 J = 1,NNO
              VDPELC(I,J,L) =  DFINV(I,1) * VDPQ2(1,J,L)  +
     +                            DFINV(I,2) * VDPQ2(2,J,L)  +
     +                            DFINV(I,3) * VDPQ2(3,J,L)
    2     CONTINUE
        END IF
C
        IF (INDIC .GE. 2) THEN
              XINT(L) = 0
              YINT(L) = 0
              ZINT(L) = 0
           DO 3 I=1,NPO
              XINT(L)  = XINT(L) + VP1(I,L) * X(I)
              YINT(L)  = YINT(L) + VP1(I,L) * Y(I)
              ZINT(L)  = ZINT(L) + VP1(I,L) * Z(I)
   3       CONTINUE
       END IF
 25   CONTINUE
      END
