      SUBROUTINE sands25(ELM,X,Y,Z,e,xnu,u,stress,strain,maxgus
     &                  ,maxstr,msize,outerr,vmflg,strainFlg)
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  BUT   : CALCUL DES Contraintes
C  ---     TETRAEDRE P2 COURBE
C          formule d'integration a 15 points, les sommets
C
C ************************************************************
C
C     STRESSES ARE COMPUTED AT THE NODES
C
C ************************************************************
C
C******************************************************************
C
C  PARAMETRES D ENTREE  :
C  -------------------
C   X,Y,Z   : TABLEAUX DES COORDONNEES DES POINTS DE L ELEMENT
c   ijt     : permutation pour oasser de la numerotation par inconnues
c             a celle par noeuds
c   nno     : nombre de noeuds de l'element
c   npo     : nombre de points
c   npi     : nombre de points d'integration
c   dp      : valeur des derivees des polynomes de base aux points
c             d'integration sur l'element de reference
c   vp1     : valeur des polynomes de base aux points
c             d'integration sur l'element de reference
c
c  tableaux de travail :
c  -------------------
c   delta   : jacobien aux points d'integration
c   (x y z)int : coordonnees des points d'integration sur 
c              l'element courant
C
C  PARAMETRE DE SORTIE  :
C  --------------------
C   sigma   : contraintes
C           
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  PROGRAMMEURS  : Charbel Farhat  1997
C ...................................................................
C
      integer elm,ndim,npo,npi,nno,maxgus,maxstr,msize,outerr
      logical vmflg,strainFlg
C     parameter (ndim=3, nno=10, npo=10, npi=15)
      parameter (ndim=3, nno=10, npo=10, npi=10)
      integer IJT(ndim*nno)
      real*8  stress(msize,maxstr,maxgus)
      real*8  strain(msize,maxstr,maxgus)
      real*8  straine(6,3*nno,npi)
      real*8 sigmae(6,nno*ndim,npi)
      real*8 u(nno*ndim)
      real*8 X(npo),Y(npo),Z(npo),xint(npi),yint(npi),zint(npi)
      real*8 elas(9),delta(npi),df(ndim,ndim),
     &                 a2(ndim,nno,npi),vp1(nno,npi),dfinv(ndim,ndim)
     &                  ,DP(ndim,nno,npi),e,xnu
      DATA IJT/ 1, 4, 7,10,13,16,19,22,25,28,
     &          2, 5, 8,11,14,17,20,23,26,29,
     &          3, 6, 9,12,15,18,21,24,27,30/
C
      DATA VP1 /  1.0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,
     &             .0, 1.0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,
     &             .0,  .0, 1.0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,
     &             .0,  .0,  .0, 1.0,  .0,  .0,  .0,  .0,  .0,  .0, 
     &             .0,  .0,  .0,  .0, 1.0,  .0,  .0,  .0,  .0,  .0, 
     &             .0,  .0,  .0,  .0,  .0, 1.0,  .0,  .0,  .0,  .0,
     &             .0,  .0,  .0,  .0,  .0,  .0, 1.0,  .0,  .0,  .0,
     &             .0,  .0,  .0,  .0,  .0,  .0,  .0, 1.0,  .0,  .0, 
     &             .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0, 1.0,  .0, 
     &             .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0, 1.0 /
C
      DATA DP/-3.0,-3.0,-3.0,-1.0,  .0,  .0,  .0,-1.0,  .0,  .0,  .0,
     &-1.0, 4.0,  .0,  .0,  .0,  .0,  .0,  .0, 4.0,  .0,  .0,  .0,
     & 4.0,  .0,  .0,  .0,  .0,  .0,  .0, 1.0, 1.0, 1.0, 3.0,  .0,  .0,
     &  .0,-1.0,  .0,  .0,  .0,-1.0,-4.0,-4.0,-4.0,  .0, 4.0,  .0,  .0,
     &  .0,  .0,  .0,  .0,  .0,  .0,  .0, 4.0,  .0,  .0,  .0, 1.0, 1.0,
     & 1.0,-1.0,  .0,  .0,  .0, 3.0,  .0,  .0,  .0,-1.0,  .0,  .0,  .0,
     & 4.0,  .0,  .0,-4.0,-4.0,-4.0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,
     &  .0, 4.0, 1.0, 1.0, 1.0,-1.0,  .0,  .0,  .0,-1.0,  .0,  .0,  .0,
     & 3.0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,  .0,-4.0,-4.0,-4.0,
     & 4.0,  .0,  .0,  .0, 4.0,  .0,-1.0,-1.0,-1.0, 1.0,  .0,  .0,  .0,
     &-1.0,  .0,  .0,  .0,-1.0,  .0,-2.0,-2.0,  .0, 2.0,  .0,  .0, 2.0,
     &  .0,  .0,  .0, 2.0,  .0,  .0, 2.0,  .0,  .0,  .0, 1.0, 1.0, 1.0,
     & 1.0,  .0,  .0,  .0, 1.0,  .0,  .0,  .0,-1.0,-2.0,-2.0,-2.0, 2.0,
     & 2.0,  .0,-2.0,-2.0,-2.0,  .0,  .0,  .0,  .0,  .0, 2.0,  .0,  .0,
     & 2.0,-1.0,-1.0,-1.0,-1.0,  .0,  .0,  .0, 1.0,  .0,  .0,  .0,-1.0,
     & 2.0,  .0,  .0, 2.0,  .0,  .0,-2.0,  .0,-2.0,  .0,  .0, 2.0,  .0,
     &  .0,  .0,  .0,  .0, 2.0,-1.0,-1.0,-1.0,-1.0,  .0,  .0,  .0,-1.0,
     &  .0,  .0,  .0, 1.0, 2.0,  .0,  .0,  .0,  .0,  .0,  .0, 2.0,  .0,
     &-2.0,-2.0,  .0, 2.0,  .0,  .0,  .0, 2.0,  .0, 1.0, 1.0, 1.0, 1.0,
     &  .0,  .0,  .0,-1.0,  .0,  .0,  .0, 1.0,-2.0,-2.0,-2.0,  .0, 2.0,
     &  .0,  .0,  .0,  .0,-2.0,-2.0,-2.0, 2.0,  .0, 2.0,  .0, 2.0,  .0,
     & 1.0, 1.0, 1.0,-1.0,  .0,  .0,  .0, 1.0,  .0,  .0,  .0, 1.0,  .0,
     &  .0,  .0, 2.0,  .0,  .0,-2.0,-2.0,-2.0,-2.0,-2.0,-2.0, 2.0,  .0,
     &  .0,  .0, 2.0, 2.0/
C
C     CAS ISOTROPE
C     -------------------------------------------
      CALL EC3C2C(ELM,    NNO,    NNO,       X,     Y,
     +              Z,    NPI,    IJT,     vp1,
     +             dp,     dp,      e,     xnu,  elas,
     +         stress, sigmae, strain, straine, delta,
     +             df,  DFINV,   xint,    yint,  zint,
     +             a2,      u, maxgus,  maxstr, msize,
     +         outerr,  vmflg, strainFlg)
c
      END

