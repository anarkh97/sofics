      SUBROUTINE sands24(ELM,X,Y,Z,e,xnu,u,stress,strain,
     &           maxgus,maxstr,msize,outerr,vmflg,strainFlg)
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C ************************************************************
C
C     STRESSES ARE COMPUTED AT THE NODES
C
C ************************************************************
C
C  BUT   : CALCUL DEs contraintes
C  ---     pentaedre r1 DROIT
c          formule d'integration a 6 points (gauss) 3*2
C *******************************************************************
C 6 nodes PENTAEDRA 6 gauss integration points
C *******************************************************************
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
C   stress   : contraintes
C           
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  PROGRAMMEURS  : Marina Vidrascu 1996
C ...................................................................
C
      integer elm,ndim,npo,npi,nno,maxgus,maxstr,msize,outerr
      logical vmflg,strainFlg
      parameter (ndim=3, nno=6, npo=6, npi=6)
      integer ijt(3*nno)
      real*8 stress(msize,maxstr,maxgus)
      real*8 strain(msize,maxstr,maxgus)
      real*8 sigmae(6,nno*ndim,npi),u(nno*ndim)
      real*8 straine(6,nno*ndim,npi)
      real*8 X(npo),Y(npo),Z(npo),xint(npi),yint(npi),zint(npi)
      real*8 elas(9),delta(npi),df(ndim,ndim),
     &                 a2(ndim,nno,npi),vp1(nno,npi),dfinv(ndim,ndim)
     &                  ,DP(ndim,nno,npi),e,xnu
c
      DATA IJT/1,4,7,10,13,16,
     +         2,5,8,11,14,17,
     +         3,6,9,12,15,18/
c
      DATA DP / -1.0, -1.0, -1.0, 1.0, .0, .0, .0, 1.0, .0,
     &          .0, .0, 1.0, .0, .0, .0, .0, .0, .0, -1.0,
     &          -1.0, .0, 1.0, .0, -1.0, .0, 1.0, .0, .0,
     &          .0, .0, .0, .0, 1.0, .0, .0, .0, -1.0, -1.0,
     &          .0, 1.0, .0, .0, .0, 1.0, -1.0, .0, .0, .0,
     &          .0, .0, .0, .0, .0, 1.0, .0, .0, -1.0, .0, .0,
     &          .0, .0, .0, .0, -1.0, -1.0, 1.0, 1.0, .0, .0, .0,
     &          1.0, .0, .0, .0, .0, .0, .0,-1.0, .0, .0, .0, -1.0,
     &         -1.0, .0, 1.0, .0, 1.0, .0, 1.0, .0, .0, .0, .0, .0,
     &         .0, .0, .0, .0, -1.0, -1.0, -1.0, .0, 1.0, .0, .0, .0,
     &         1.0, 1.0/
C
      DATA VP1 / 1.0, .0, .0, .0, .0, .0, .0, 1.0, .0,
     &           .0, .0, .0, .0, .0, 1.0, .0, .0, .0, .0,
     &           .0, .0, 1.0, .0, .0, .0, .0, .0, .0,
     &           1.0, .0, .0, .0, .0, .0, .0, 1.0/

C
C     CAS ISOTROPE
C     -------------------------------------------
      CALL EC3C2C (ELM,NNO,NNO,X,Y,Z,NPI,IJT,vp1,
     +                   dp,dp,e,xnu,
     +                   elas,stress,sigmae,strain,
     +                   straine,delta,df,DFINV,
     +                   xint,yint,zint,a2,u,
     +                   maxgus,maxstr,msize,outerr,vmflg,
     +			 strainFlg)
c
       END
