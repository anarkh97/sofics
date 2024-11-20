C=====================================================================C
C        1         2         3         4         5         6         7C
C234567890123456789012345678901234567890123456789012345678901234567890C
C=====================================================================C
      subroutine compms( x       , y       , z     , h      , rho     ,
     $                   emass   , medof   , nttly , ncmpfr , elm     ,
     $                   idlay   , mtlay   , cmpfr , iatt   , ctyp    ,
     $                   catt    , cfrm    , gamma , grvfor , grvflg  ,
     $                   totmas  , area  , masflg                     )
C=====================================================================C
C                                                                     C
C     Performs =   This subroutine will form the elemental mass       C
C     ----------   matrix of the 3D 3-node ANDES composite shell.     C
C                  Lumping is assumed here.                           C
C                                                                     C
C                                                                     C
C     Inputs/Outputs =                                                C
C     ----------------                                                C
C     X        <input>   nodal coordinates in the X-direction         C
C     Y        <input>   nodal coordinates in the Y-direction         C
C     Z        <input>   nodal coordinates in the Z-direction         C
C     H        <input>   element thicknesses (assumed constant)       C
C     RHO      <input>   density                                      C
C     EMASS    <output>  element mass matrix                          C
C     MEDOF    <input>   maximum number of DOFs per FE                C
C     NTTLY    <input>   total number of composite layers             C
C     NCMPFR   <input>   number of composite frames                   C
C     ELM      <input>   finite element number                        C
C     IDLAY    <input>   identificators for the composite layers      C
C     MTLAY    <input>   material properties of the composite layers  C
C     IATT     <input>   attribute number of the composite shell      C
C     CTYP     <input>   type of composite law (type 0, 1, 2, or 3)   C
C     CATT     <input>   composite attribute number of the element    C
C     CFRM     <input>   frame number of the composite shell element  C
C                                                                     C
C                                                                     C
C     Computations =                                                  C
C     ---------------                                                 C
C                                                                     C
C     The lumped mass matrix [M] is equal to:                         C
C                                                                     C
C           [ mt 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ]  C
C           [ 0  mt 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ]  C
C           [ 0  0  mt 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 ]  C
C           [ 0  0  0  m1 0  0  0  0  0  0  0  0  0  0  0  0  0  0 ]  C
C           [ 0  0  0  0  m1 0  0  0  0  0  0  0  0  0  0  0  0  0 ]  C
C           [ 0  0  0  0  0  m1 0  0  0  0  0  0  0  0  0  0  0  0 ]  C
C           [ 0  0  0  0  0  0  mt 0  0  0  0  0  0  0  0  0  0  0 ]  C
C           [ 0  0  0  0  0  0  0  mt 0  0  0  0  0  0  0  0  0  0 ]  C
C           [ 0  0  0  0  0  0  0  0  mt 0  0  0  0  0  0  0  0  0 ]  C
C     [M] = [ 0  0  0  0  0  0  0  0  0  m2 0  0  0  0  0  0  0  0 ]  C
C           [ 0  0  0  0  0  0  0  0  0  0  m2 0  0  0  0  0  0  0 ]  C
C           [ 0  0  0  0  0  0  0  0  0  0  0  m2 0  0  0  0  0  0 ]  C
C           [ 0  0  0  0  0  0  0  0  0  0  0  0  mt 0  0  0  0  0 ]  C
C           [ 0  0  0  0  0  0  0  0  0  0  0  0  0  mt 0  0  0  0 ]  C
C           [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  mt 0  0  0 ]  C
C           [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  m3 0  0 ]  C
C           [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  m3 0 ]  C
C           [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  m3]  C
C                                                                     C
C     with the following ordering of local degrees of freedom:        C
C                                                                     C
C                                [     U_x1 ]                         C
C                                [     U_y1 ]                         C
C                                [     U_z1 ]                         C
C                                [ theta_x1 ]                         C
C                                [ theta_y1 ]                         C
C                                [ theta_z1 ]                         C
C                                [     U_x2 ]                         C
C                                [     U_y2 ]                         C
C                                [     U_z2 ]                         C
C     Ordering_of_Local_DOFs  =  [ theta_x2 ]                         C
C                                [ theta_y2 ]                         C
C                                [ theta_z2 ]                         C
C                                [     U_x3 ]                         C
C                                [     U_y3 ]                         C
C                                [     U_z3 ]                         C
C                                [ theta_x3 ]                         C
C                                [ theta_y3 ]                         C
C                                [ theta_z3 ]                         C
C                                                                     C
C     No rotation of local-to-global basis is implemented since the   C
C     mass matrix [M] is formed of 3 by 3 blocks proportional to the  C
C     identity. The lumping factors are equal to:                     C
C                                                                     C
C     [mt] = [rho] * [A] * [h]        /    3.0                        C
C     [m1] = [rho] * [A] * [h] * [Ix] / 1260.0                        C
C     [m2] = [rho] * [A] * [h] * [Iy] / 1260.0                        C
C     [m3] = [rho] * [A] * [h] * [Iz] / 1260.0                        C
C                                                                     C
C     where:                                                          C
C                                                                     C
C     [rho]                 equivalent density of the composite       C
C     [h]                   equivalent thickness of the composite     C
C     [A]                   area                                      C
C     [Ix], [Iy] and [Iz]   equivalent pseudo-moments of inertia      C
C                                                                     C
C     The computation of these various quantities is performed        C
C     according to the type of composite law prescribed by the        C
C     variable [CTYP]. Four types of composite laws are available:    C
C                                                                     C
C     type-0: isotropic element                                       C
C                                                                     C
C     type-1: constitutive coefficients are given                     C
C                                                                     C
C     type-2: properties of each layer are given and no coupling      C
C             is assumed between bending and membrane                 C
C                                                                     C
C     type-3: properties of each layer are given and coupling         C
C             between bending and membrane is assumed                 C
C                                                                     C
C     The meaning of each input argument is specified in the          C
C     following as a function of the type of constitutive law.        C
C                                                                     C
C     1. Constitutive Law of Type-0                                   C
C     - - - - - - - - - - - - - - -                                   C
C                                                                     C
C     [rho]                 density of the isotropic material         C
C     [h]                   thickness of the isotropic shell          C
C     [A]                   area of the shell                         C
C     [Ix], [Iy] and [Iz]   pseudo-moments of inertia of the shell    C
C                                                                     C
C     Quantities [A], [Ix], [Iy] and [Iz] are obtained via            C
C     numerical integration.                                          C
C                                                                     C
C     2. Constitutive Law of Type-1                                   C
C     - - - - - - - - - - - - - - -                                   C
C                                                                     C
C     The density parameter must be initialized in the input file as  C
C     the density per unit surface of the composite element. In that  C
C     case, the thickness [h] is not required whereas quantities [A], C
C     [Ix], [Iy] and [Iz] are obtained via numerical integration once C
C     again. The lumping factors are equal to:                        C
C                                                                     C
C     [mt] = [rho_surface] * [A]        /    3.0                      C
C     [m1] = [rho_surface] * [A] * [Ix] / 1260.0                      C
C     [m2] = [rho_surface] * [A] * [Iy] / 1260.0                      C
C     [m3] = [rho_surface] * [A] * [Iz] / 1260.0                      C
C                                                                     C
C     where the density per unit surface [rho_surface] is stored in   C
C     the same variable as before (type-0), that is, [rho].           C
C                                                                     C
C     3. Constitutive Law of Type-2                                   C
C     - - - - - - - - - - - - - - -                                   C
C                                                                     C
C     The density parameter must be initialized in the input file as  C
C     the non-structural mass density per unit surface of the         C
C     composite element.                                              C
C     The mass coefficients are computed by adding together the       C
C     contributions of each layer of the composite shell, and         C
C     also the non-structural mass.                                   C
C                                                                     C
C     [mt] = ( nsm + sum{ [rho_k] * [h_k] } ) * [A]        /    3.0   C
C     [m1] = ( nsm + sum{ [rho_k] * [h_k] } ) * [A] * [Ix] / 1260.0   C
C     [m2] = ( nsm + sum{ [rho_k] * [h_k] } ) * [A] * [Iy] / 1260.0   C
C     [m3] = ( nsm + sum{ [rho_k] * [h_k] } ) * [A] * [Iz] / 1260.0   C
C                                                                     C
C     where:                                                          C
C                                                                     C
C     [rho_k]               density of the layer number [k]           C
C     [h_k]                 thickness of the layer number [k]         C
C     [A]                   area of the shell                         C
C     [Ix], [Iy] and [Iz]   pseudo-moments of inertia of the shell    C
C                                                                     C
C     where the non-structural mass density per unit surface [nsm] is C
C     stored in the same variable as before (type-0), that is, [rho]. C
C                                                                     C
C     Quantities [A], [Ix], [Iy] and [Iz] are obtained via            C
C     numerical integration. Densities (per unit volume) [rho_k] and  C
C     thicknesses [h_k] for each layer [k] of the composite are       C
C     given in the input file and retrieved from arrays [IDLAY] and   C
C     [MTLAY] using the information stored in the column that         C
C     corresponds to the particular layer number [k]:                 C
C                                                                     C
C     [IDLAY] Row 1: attribute number as read in the input file       C
C     [IDLAY] Row 2: number of layers of the composite shell          C
C     [IDLAY] Row 3: layer number (that is, number [k])               C
C     [IDLAY] Row 4: type of constitutive law (2 or 3)                C
C     [IDLAY] Row 5: frame number for the fibers reference vector     C
C                                                                     C
C     [MTLAY] Row 1: orthotropic Young modulus E1                     C
C     [MTLAY] Row 2: orthotropic Young modulus E2                     C
C     [MTLAY] Row 3: orthotropic Poisson's ratio nu12                 C
C     [MTLAY] Row 4: orthotropic shear mogulus G12                    C
C     [MTLAY] Row 5: first coefficient of mutual influence mu1,12     C
C     [MTLAY] Row 6: second coefficient of mutual influence mu2,12    C
C     [MTLAY] Row 7: density of the layer number [k]                  C
C     [MTLAY] Row 8: thickness of the layer number [k]                C
C     [MTLAY] Row 9: angle (degrees) between fibers and reference     C
C                                                                     C
C     4. Constitutive Law of Type-3                                   C
C     - - - - - - - - - - - - - - -                                   C
C                                                                     C
C     Same as Type-2.                                                 C
C                                                                     C
C                                                                     C
C     Caution =                                                       C
C     ---------                                                       C
C     The finite element is assumed to have constant thickness so     C
C     that no numerical interpolation is required. The array storing  C
C     the thicknesses at each nodal points [h] is reduced to a scalar C
C     in this routine. Also, The outputed element mass matrix is a    C
C     (diagonal) 18 by 18 block stored in the left upper corner of    C
C     the output storage [EMASS].                                     C
C                                                                     C
C                                                                     C
C     Outputs =   no output.                                          C
C     ---------                                                       C
C                                                                     C
C=====================================================================C
C=Author  = Francois M. Hemez                                         C
C=Date    = June 9, 1995                                              C
C=Version = 2.0                                                       C
C=Comment =                                                           C
C=====================================================================C
C
C     -----------
C     DECLARATION
C     -----------
C
C.....Globlal Variables
C
      integer    nttly , ctyp , catt , cfrm , elm
      integer    idlay(5,nttly) , iatt , ncmpfr , medof
C
      real*8     rho , emass(medof,medof) , mtlay(12,nttly)
      real*8     h(3) , cmpfr(9,ncmpfr) , x(3) , y(3) , z(3)
      real*8     totmas , area , gamma(*) , grvfor(*)

      logical    grvflg , masflg
C
C.....Local Maximum Number of Layers per Composite Element
C
      integer    maxlayer
      parameter( maxlayer = 1000 )
C
C.....Local Variables
C
      integer    i , j , i1 , i2 , i3 , ilayer , nlayer
C
      real*8     zero , thick , mass0 , mass1 , mass2 , mass3
      real*8     x13 , y13 , z13 , hlayer(maxlayer)
      real*8     x32 , y32 , z32 , rholayer(maxlayer)
      real*8     x21 , y21 , z21 , rhoh
      real*8     dist(3) , rlr , rlb , bpr , twicearea2
      real*8     Ix , Iy , Iz
C
C     ----
C     DATA
C     ----
C
      data zero  /0.000000D+00/
C
C     -----
C     LOGIC
C     -----
C
C.....CLEAR THE MASS LUMPING FACTORS
C
      mass0 = zero
      mass1 = zero
      mass2 = zero
      mass3 = zero
C
C.....CLEAR THE LOCAL STORAGES FOR DENSITY AND THICKNESS PER LAYER
C
      do 1001 ilayer=1,maxlayer
         rholayer(ilayer) = zero
 1001 continue
C
      do 1002 ilayer=1,maxlayer
         hlayer(ilayer) = zero
 1002 continue
C
C.....CLEAR THE NUMBER OF LAYERS OF THE ELEMENT
C
      nlayer = 0
C
C.....CLEAR THE OUTPUT MASS MATRIX
C
      do 1003 j=1,medof
         do 1004 i=1,medof
            emass(i,j) = zero
 1004    continue
 1003 continue
C
C     --------------------------
C     CHECKS AND INITIALIZATIONS
C     --------------------------
C
C.....CHECK THE TYPE OF CONSTITUTIVE LAW FOR THIS ELEMENT
C
      if (      (ctyp.ne.0)
     $     .and.(ctyp.ne.1)
     $     .and.(ctyp.ne.2)
     $     .and.(ctyp.ne.3) ) go to 100
C
C.....CHECK THE ADDRESSING IN ARRAYS [IDLAY] AND [MTLAY]
C
      if ( (ctyp.eq.2).or.(ctyp.eq.3) ) then
         if ( (catt.lt.1).or.(catt.gt.nttly) ) go to 200
      endif
C
C.....EXTRACT THE LAYER INFORMATION FOR TYPES 2 AND 3 LAWS
C
      if ( (ctyp.eq.2).or.(ctyp.eq.3) ) then
         nlayer = idlay(2,catt)
         if ( nlayer.gt.maxlayer ) go to 300
         do 2001 i=catt,(catt+nlayer-1)
            ilayer           = idlay(3,i)
            if ( ilayer.gt.nlayer ) go to 400
            rholayer(ilayer) = mtlay(7,i)
            hlayer(ilayer)   = mtlay(8,i)
 2001    continue
      endif
C
C.....COMPUTE THE DISTANCE BETWEEN X-, Y- AND Z- NODAL COORDINATES
C
      x21 = x(2) - x(1)
      y21 = y(2) - y(1)
      z21 = z(2) - z(1)
C
      x32 = x(3) - x(2)
      y32 = y(3) - y(2)
      z32 = z(3) - z(2)
C
      x13 = x(1) - x(3)
      y13 = y(1) - y(3)
      z13 = z(1) - z(3)
C
C.....COMPUTE THE DISTANCE BETWEEN NODES 1-2, 2-3 AND 3-1
C
      dist(1) = sqrt( x21*x21 + y21*y21 + z21*z21 )
      dist(2) = sqrt( x32*x32 + y32*y32 + z32*z32 )
      dist(3) = sqrt( x13*x13 + y13*y13 + z13*z13 )
C
C.....COMPUTE THE LENGTH OF SIDE 1-2
C
      rlr = sqrt( x21*x21 + y21*y21 + z21*z21 )
C
C.....CHECK FOR ZERO-SIDE LENGTH
C
      if ( rlr.eq.zero ) go to 500
C
C.....COMPUTE THE DISTANCE OF THE OPPOSING NODE (3) TO THAT SIDE (1-2)
C
      rlb = sqrt( x32*x32 + y32*y32 + z32*z32 )
      bpr =  abs( x21*x32 + y21*y32 + z21*z32 )/rlr
C
C.....COMPUTE THE SQUARE OF TWICE THE TRIANGLE'S AREA
C
      twicearea2 = rlb*rlb - bpr*bpr
C
C.....CHECK IF THE TRIANGLE'S AREA IS POSITIVE
C
      if ( twicearea2.le.zero ) go to 600
C
C.....COMPUTE THE AREA OF THE TRIANGLE
C
      area = 0.50D+00*rlr*sqrt(twicearea2)
C
C.....COMPUTE THE THREE PSEUDO MOMENTS OF INERTIA
C
      Ix = dist(1)*dist(1) + dist(3)*dist(3)
      Iy = dist(1)*dist(1) + dist(2)*dist(2)
      Iz = dist(2)*dist(2) + dist(3)*dist(3)
C
C     -----------------------
C     TYPE-0 CONSTITUTIVE LAW
C     -----------------------
C
      if ( ctyp.eq.0 ) then
C
C.....INITIALIZE THE ELEMENT'S CONSTANT THICKNESS
C
      thick = h(1)
C
C.....FORM THE MASS COEFFICIENTS PER DEGREE OF FREEDOM
C
      mass0 = rho*area*thick/3.00D+00
      mass1 = rho*area*thick*Ix/1260.00D+00
      mass2 = rho*area*thick*Iy/1260.00D+00
      mass3 = rho*area*thick*Iz/1260.00D+00
C
C.....END OF TREATMENT FOR A TYPE-0 CONSTITUTIVE LAW
C
      endif
C
C     -----------------------
C     TYPE-1 CONSTITUTIVE LAW
C     -----------------------
C
      if ( ctyp.eq.1 ) then
C
C.....FORM THE MASS COEFFICIENTS PER DEGREE OF FREEDOM
C.....WARNING: THE "DENSITY" PARAMETER PRESCRIBED IN THE
C.....ATTRIBUTE SECTION OF THE INPUT FILE AND USED HERE
C.....MUST BE A DENSITY PER UNIT SURFACE
C
      mass0 = rho*area/3.00D+00
      mass1 = rho*area*Ix/1260.00D+00
      mass2 = rho*area*Iy/1260.00D+00
      mass3 = rho*area*Iz/1260.00D+00
C
C.....INITIALIZE THE ELEMENT'S CONSTANT THICKNESS
C
      thick = h(1)
C
C.....END OF TREATMENT FOR A TYPE-1 CONSTITUTIVE LAW
C
      endif
C
C     -----------------------------------
C     TYPE-2 AND TYPE-3 CONSTITUTIVE LAWS
C     -----------------------------------
C
      if ( (ctyp.eq.2).or.(ctyp.eq.3) ) then
C
C.....ACCUMULATE THE PRODUCT DENSITY BY THICKNESS PER LAYER
C.....ALSO ACCOUNTING FOR NON-STRUCTURAL MASS
C
      rhoh = rho
C
      do 3001 ilayer=1,nlayer
         rhoh = rhoh + rholayer(ilayer)*hlayer(ilayer)
 3001 continue

C
C.....FORM THE MASS COEFFICIENTS PER DEGREE OF FREEDOM
C
      mass0 = rhoh*area/3.00D+00
      mass1 = rhoh*area*Ix/1260.00D+00
      mass2 = rhoh*area*Iy/1260.00D+00
      mass3 = rhoh*area*Iz/1260.00D+00
C
C.....END OF TREATMENT FOR TYPE-2 AND TYPE-3 CONSTITUTIVE LAWS
C
      endif
C
C     -------------------------------------
C     ASSEMBLY OF THE ELEMENTAL MASS MATRIX
C     -------------------------------------
C
C.....FORM THE LUMPED ELEMENT MASS MATRIX
C
      do 4001 i=1,3
         i2           = i+6
         i3           = i+12
         emass(i ,i ) = mass0
         emass(i2,i2) = mass0
         emass(i3,i3) = mass0
 4001 continue
C
      do 4002 i=1,3
         i1           = i+3
         i2           = i+9
         i3           = i+15
         emass(i1,i1) = mass1
         emass(i2,i2) = mass2
         emass(i3,i3) = mass3
 4002 continue
C
C..... COMPUTE THE BODY FORCE DUE TO GRAVITY IF NEEDED
C
      if (grvflg) then
        grvfor(1) = 3.0d0*mass0*gamma(1)
        grvfor(2) = 3.0d0*mass0*gamma(2)
        grvfor(3) = 3.0d0*mass0*gamma(3)
      endif

C
C.... ACCUMULATE THE SUBDOMAIN MASS
C
        if (masflg) then
          totmas = totmas + 3.0d0*mass0
        endif
C
C     ------
C     RETURN
C     ------
C
      return
C
C     ---------------
C     ERROR TREATMENT
C     ---------------
C
C.....ERROR-MESSAGE IF UNKNOWN TYPE OF CONSTITUTIVE LAW
C
  100 continue
      write(*,*) "*** FATAL ERROR in routine COMPMS. ***"
      write(*,*) "*** The constitutive law type      ***"
      write(*,*) "*** must be either 0, 1, 2, or 3.  ***"
      write(*,*) "*** Execution terminated.          ***"
      stop
C
C.....ERROR-MESSAGE IF THE ADDRESSING IS NOT CORRECT
C
  200 continue
      write(*,*) "*** FATAL ERROR in Routine COMPMS ***"
      write(*,*) "*** The Address in Arrays [IDLAY] ***"
      write(*,*) "*** and [MTLAY] is Out-of-Bounds  ***"
      write(*,*) "*** Execution Terminated Here     ***"
      stop
C
C.....ERROR-MESSAGE IF THE MAXIMUM NUMBER OF LAYERS IS EXCEEDED
C
  300 continue
      write(*,*) "*** FATAL ERROR in Routine COMPMS     ***"
      write(*,*) "*** Maximum Number of Layers Exceeded ***"
      write(*,*) "*** Boost Local Parameter [MAXLAYER]  ***"
      write(*,*) "*** Execution Terminated Here         ***"
      stop
C
C.....ERROR-MESSAGE IF THE TOTAL NUMBER OF LAYERS IS NOT CONSISTENT
C
  400 continue
      write(*,*) "*** FATAL ERROR in Routine COMPMS      ***"
      write(*,*) "*** A Layer Number Exceeds the Total   ***"
      write(*,*) "*** Number of Layers Stored in [IDLAY] ***"
      write(*,*) "*** Execution Terminated Here          ***"
      stop
C
C.....ERROR-MESSAGE IF A SIDE HAS ZERO-LENGTH
C
  500 continue
      write(*,*) "*** FATAL ERROR in routine COMPMS     ***"
      write(*,*) "*** The Side 1-2 has Zero Length      ***"
      write(*,*) "*** Check Coordinates and FE Topology ***"
      write(*,*) "*** Execution Terminated Here         ***"
      stop
C
C.....ERROR-MESSAGE IF THE AREA IS NEGATIVE OR ZERO
C
  600 continue
      write(*,*) "*** FATAL ERROR in routine COMPMS     ***"
      write(*,*) "*** The Area is Negative or Zero      ***"
      write(*,*) "*** Check Coordinates and FE Topology ***"
      write(*,*) "*** ... Counterclock Nodal Numbering? ***"
      write(*,*) "*** Execution Terminated Here         ***"
      stop
C
C     ---
C     END
C     ---
C
      end
C=end of routine "COMPMS"
C=======================C
