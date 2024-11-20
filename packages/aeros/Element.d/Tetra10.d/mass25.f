      SUBROUTINE MASS25(X,Y,Z,ro,A,nmax,gamma,grvfor,grvflg,
     &                  totmas,masflg,dp,vp1)
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  BUT   : CALCUL DE LA MATRICE ELEMENTAIRE ELASTIQUE DE MASSE
C  ---     TETRAEDRE P2 COURBE
c          formule d'integration a 15 points, les sommets
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
c   poids   : poids de la formule d'integration
c
c  tableaux de travail :
c  -------------------
c   delta   : jacobien aux points d'integration
c   (x y z)int : coordonnees des points d'integration sur 
c              l'element courant
C
C  PARAMETRE DE SORTIE  :
C  --------------------
C   A      : MATRICE DE MASSE. LUMPED.
C            STockee dans a(nmax,nmax)
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  PROGRAMMEURS  : Charbel Farhat 1997
C ...................................................................
C
      logical rdf25
      common /E25/rdf25
      integer ndim,nno,nmax,i,j,npo,npi
      parameter (ndim=3, nno=10, npo=10, npi=15)
      real*8 A(nmax,nmax),ae(NNO*3)
      real*8 X(npo),Y(npo),Z(npo),xint(npi),yint(npi),zint(npi)
      integer IJT(ndim*nno)
      real*8  dp(3,10,15), vp1(10,15)
      real*8 poids(npi),delta(npi),df(ndim,ndim),
     &                 a2(ndim,nno,npi),dfinv(ndim,ndim),ro
c
*     common /E25DATA/DP,VP1
c
      logical grvflg,masflg
      real*8  totmas,gamma(*),grvfor(*)
c
      DATA IJT/ 1, 4, 7,10,13,16,19,22,25,28,
     &          2, 5, 8,11,14,17,20,23,26,29,
     &          3, 6, 9,12,15,18,21,24,27,30/
C
      DATA POIDS / 1.975308731198311E-02, 1.198951396316977E-02,
     &             1.198951396316977E-02, 1.198951396316977E-02,
     &             1.198951396316977E-02, 1.151136787104540E-02,
     &             1.151136787104540E-02, 1.151136787104540E-02,
     &             1.151136787104540E-02, 8.818342350423336E-03,
     &             8.818342350423336E-03, 8.818342350423336E-03,
     &             8.818342350423336E-03, 8.818342350423336E-03,
     &             8.818342350423336E-03/
C
C*****************************************************************
C***** READ DP and VP1 FROM INPUT FILES
C*****************************************************************
C
*      if(rdf25) then
*      open(unit=1,file="DPE25")
*      read(1,*) (((DP(i,j,k),i=1,ndim),j=1,nno),k=1,npi)
*      close(1)
*      open(unit=1,file="VP1E25")
*      read(1,*) ((VP1(i,j),i=1,nno),j=1,npi)
*      close(1)
       rdf25 = .false.
*      endif
C*****
C*****************************************************************
C
C     CAS ISOTROPE
C     -------------------------------------------
      CALL Em3C2C (NNO,NNO,X,Y,Z,NPI,IJT,poids,vp1,
     +                   dp,dp,ro,
     +                   AE,delta,df,DFINV,
     +                   xint,yint,zint,a2,
     +                   gamma,grvfor,grvflg,totmas,masflg)
c
c     symmetriser a
c
      do 1 i = 1,3*nno
      a(i,i) = ae(i)
      do 1 j = i+1,3*nno
         a(i,j) = 0
         a(j,i) = a(i,j)
 1    continue
      END
