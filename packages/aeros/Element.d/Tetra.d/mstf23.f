      SUBROUTINE MSTF23(X,Y,Z,e,xnu,A,nmax,status)
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  BUT   : CALCUL DE LA MATRICE ELEMENTAIRE ELASTIQUE DE RIGIDITE
C  ---     TETRAEDRE P1 DROIT
c          formule d'integration a 4 points, les sommets
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
C   AE     : MATRICE DE RIGIDITE . SYMETRIQUE
C            STockee dans A(nmax,nmax)
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  PROGRAMMEURS  : Marina Vidrascu 1996
C ...................................................................
C
      integer ndim,nno,npo,npi,nmax,i,j
      parameter (ndim=3, nno=4, npo=4, npi=4)
      real*8  A(nmax,nmax),ae(NNO*3*(NNO*3+1)/2)
      real*8  X(npo),Y(npo),Z(npo),xint(npi),yint(npi),zint(npi)
      integer IJT(ndim*nno), status
      real*8    poids(npi),elas(9),delta(npi),df(ndim,ndim),
     &          a2(ndim,nno,npi),vp1(nno,npi),dfinv(ndim,ndim),
     &          DP(ndim,nno,npi),e,xnu
c
      DATA DP/-1.D0 ,-1.D0 ,-1.D0 , 1.D0,3*0.D0,1.D0,3*0.D0,1.D0,
     &        -1.D0 ,-1.D0 ,-1.D0 , 1.D0,3*0.D0,1.D0,3*0.D0,1.D0,   
     &        -1.D0 ,-1.D0 ,-1.D0 , 1.D0,3*0.D0,1.D0,3*0.D0,1.D0,   
     &        -1.D0 ,-1.D0 ,-1.D0 , 1.D0,3*0.D0,1.D0,3*0.D0,1.D0/
      DATA IJT / 1, 4, 7, 10,  2, 5, 8, 11,  3, 6, 9, 12/
      data vp1 /1.,0.,0.,0.,  0.,1.,0.,0.,  0.,0.,1.,0.,
     &          0.,0.,0.,1./
C
      data poids / 4*  4.166666666666666E-02/
c     poids = .25/6.
C
C     CAS ISOTROPE
C     -------------------------------------------
      CALL ER3C2C (NNO,NNO,X,Y,Z,NPI,IJT,poids,vp1,
     +                   dp,dp,e,xnu,
     +                   elas,AE,delta,df,DFINV,
     +                   xint,yint,zint,a2,status)
c
c     symmetriser a
c
      do 1 i = 1,3*nno
      do 1 j = i,3*nno
         a(i,j) = ae(j*(j-1)/2+i)
         a(j,i) = a(i,j)
 1    continue
      END
*******************************
