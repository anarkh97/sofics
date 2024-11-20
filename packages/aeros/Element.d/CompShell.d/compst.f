C===================================================================C
      subroutine compst(  e     ,  elm   , h      , estiff , medof  ,
     $                    nu    ,  X     , Y      , Z      , nttco  ,
     $                    nttly , ncmpfr , cmpco  , idlay  , mtlay  ,
     $                    cmpfr , iatt   , ctyp   , catt   , cfrm,
     $                    flag   )
C===================================================================C
C                                                                   C
C     Perform =    This subroutine will form the element stiffness  C
C     ---------    matrix for the 3D 3-node composite shell element C
C                  derived from the ANDES formulation (ANS shell    C
C                  element with Assumed Quadratic Rotations).       C
C                                                                   C
C                                                                   C
C     Inputs/Outputs =                                              C
C     ----------------                                              C
C     E       <input>  Young modulus                                C
C     ELM     <input>  finite element number                        C
C     H       <input>  element thickness (assumed constant)         C
C     ESTIFF  <output> element stiffness matrix                     C
C     MEDOF   <input>  maximum number of DOFs per FE                C
C     NU      <input>  Poisson's ratio                              C
C     X       <input>  nodal coordinates in the X-direction         C
C     Y       <input>  nodal coordinates in the Y-direction         C
C     Z       <input>  nodal coordinates in the Z-direction         C
C     NTTCO   <input>  number of attributes for laws of type-1      C
C     NTTLY   <input>  total number of composite layers             C
C     NCMPFR  <input>  number of frames for composite shells        C
C     CMPCO   <input>  constitutive coefficients of the attributes  C
C     IDLAY   <input>  identificators for the composite layers      C
C     MTLAY   <input>  material properties of the composite layers  C
C     CMPFR   <input>  storage of frames for composite shells       C
C     IATT    <input>  attribute number of the element              C
C     CTYP    <input>  composite attribute number                   C
C     CATT    <input>  starting address for the constitutive law    C
C     CFRM    <input>  frame number for definition of the shell     C
C     FLAG    <input>  integer specifying whether to return         C
C                      transformed element stiffness matrix or      C
C                      global element stiffness matrix              C
C                                                                   C
C                                                                   C
C     Computations =                                                C
C     --------------                                                C
C     This subroutine evaluates the stiffness matrix for the 18     C
C     degrees of freedom 3-node composite triangle obtained as a    C
C     combination of the Assumed Quadratic Rotations bending        C
C     triangle plus the membrane with driling degrees of freedom    C
C     developed by Militello, Felippa et al. For documentation,     C
C     see Carmelo Militello's doctoral dissertation, pp112-113.     C
C                                                                   C
C     The original version of the ANS shell element has been        C
C     generalized here to the case of a composite shell element     C
C     with complete bending-membrane coupling. Four types of        C
C     composite laws are available via the input:                   C
C        type-0: isotropic element                                  C
C        type-1: constitutive coefficients are given                C
C        type-2: properties of each layer are given and no coupling C
C                   is assumed between bending and membrane         C
C        type-3: properties of each layer are given and coupling    C
C                   between bending and membrane is assumed         C
C                                                                   C
C     The stiffness matrix [K] is assembled as the combination of   C
C     the basic stiffness and the higher order stiffness matrices.  C
C     There are two of such matrices for each physical effect.      C
C     In the most general case, the stiffness [K] is given as:      C
C        [K] = [K_basic       _Bending _Bending ]                   C
C            + [K_basic       _Membrane_Membrane]                   C
C            + [K_basic       _Bending _Membrane]                   C
C            + [K_basic       _Membrane_Bending ]                   C
C            + [K_higher_order_Bending _Bending ]                   C
C            + [K_higher_order_Membrane_Membrane]                   C
C            + [K_higher_order_Bending _Membrane]                   C
C            + [K_higher_order_Membrane_Bending ]                   C
C     In general, a particular stiffness matrix may be expressed as C
C     the product of the following quantities:                      C
C         [K] = sum_{i=1,2,3} [B1_i]^T [C] [B2_i]                   C
C     where [i] represents the numerical integration index; [B1_i]  C
C     and [B2_i] are either the moment-curvature or force-strain    C
C     matrices (for bending or membrane effect, respectively); and  C
C     [C] represents the constitutive law. See the following        C
C     routines for details:                                         C
C     "compbBB" assembly of [K_basic       _Bending _Bending ]      C
C     "compbMM" assembly of [K_basic       _Membrane_Membrane]      C
C     "compbBM" assembly of [K_basic       _Bending _Membrane]      C
C     "compbMB" assembly of [K_basic       _Membrane_Bending ]      C
C     "comphBB" assembly of [K_higher_order_Bending _Bending ]      C
C     "comphMM" assembly of [K_higher_order_Membrane_Membrane]      C
C     "comphBM" assembly of [K_higher_order_Bending _Membrane]      C
C     "comphMB" assembly of [K_higher_order_Membrane_Bending ]      C
C                                                                   C
C                                                                   C
C     Caution =                                                     C
C     ---------                                                     C
C     The finite element is assumed to have a constant thickness    C
C     so that no numerical interpolation is required.  The array    C
C     storing thicknesses at each nodal points, [H], has a size     C
C     equal to 3 here. Also, The outputed element stiffness matrix  C
C     is a 18 by 18 block stored in the left upper corner of the    C
C     output storage [ESTIFF]. Finally, the maximum number of       C
C     composite layers per finite element is assumed to be equal to C
C     one thousand (upgrade local parameter [MAXLAYER] if not large C
C     enough).                                                      C
C                                                                   C
C                                                                   C
C     Output =   no output.                                         C
C     --------                                                      C
C                                                                   C
C===================================================================C
C=Author  = Francois M. Hemez                                       C
C=Date    = June 9, 1995                                            C
C=Version = 2.0                                                     C
C=Comment =                                                         C
C===================================================================C
C
C     -----------
C     DECLARATION
C     -----------
C
C.....GLOBAL VARIABLES
C
      integer    medof , elm , idlay(5,nttly)
      integer    nttco , nttly , ncmpfr, flag
      integer    iatt , ctyp , catt , cfrm
C
      real*8     nu , h(3) , e , estiff(medof,medof)
      real*8     X(3) , Y(3) , Z(3) , cmpco(36,nttco)
      real*8     mtlay(12,nttly) , cmpfr(9,ncmpfr)
C
C.....LOCAL DIMENSION OF THE STIFFNESS MATRIX
C.....(3 NODES AND 6 DOFS PER NODE = 18 DOFS TOTAL)
C
      integer    ndof
      parameter( ndof = 18 )
C
C.....LOCAL MAXIMUM NUMBER OF LAYERS PER COMPOSITE ELEMENT
C.....(SET TO 1,000 HERE)
C
      integer    maxlayer
      parameter( maxlayer = 1000 )
C
C.....LOCAL VARIABLES
C
      integer    i , j , nlayer , ilayer
      integer    rowb(9) , colb(9) , rowm(9) , colm(9)
      integer    idcmp23(5,maxlayer)
C
      real*8     rot(6,6) , xlp(3) , ylp(3) , zlp(3)
      real*8     zero , one , thick , cstcoef(36)
      real*8     cstbb(3,3) , cstmm(3,3) , onehalf
      real*8     cstbm(3,3) , cstmb(3,3) , point32
      real*8     kbBB(ndof,ndof) , kbMM(ndof,ndof)
      real*8     kbBM(ndof,ndof) , kbMB(ndof,ndof)
      real*8     khBB(ndof,ndof) , khMM(ndof,ndof)
      real*8     khBM(ndof,ndof) , khMB(ndof,ndof)
      real*8     Lb(9,3) , Pb(9,3) , mtcmp23(8,maxlayer)
      real*8     Lh1(9,3) , Lh2(9,3) , Lh3(9,3)
      real*8     Ph1(9,3) , Ph2(9,3) , Ph3(9,3)
      real*8     aframe(3,3) , eframe(3,3)
      logical    fastcal
C
C     ----
C     DATA
C     ----
C
      data zero    /0.000000D+00/
      data one     /1.000000D+00/
      data onehalf /1.500000D+00/
      data point32 /0.320000D+00/
C
C.....INITIALIZE THE LOGICAL FOR ENFORCING FASTER COMPUTATIONS
C
      data fastcal /.true./
C
C     -----
C     LOGIC
C     -----
C
C.....CHECK DIMENSION OF THE ELEMENTAL STIFFNESS MATRIX
C
      if ( medof.ne.ndof ) go to 100
C
C.....CLEAR THE OUTPUT ELEMENT STIFFNESS MATRIX
C
      do 1001 j=1,medof
         do 1002 i=1,medof
            estiff(i,j) = zero
 1002    continue
 1001 continue
C
C.....CLEAR THE LOCAL CONSTITUTIVE MATRICES
C
      do 1003 j=1,3
         do 1004 i=1,3
            cstbb(i,j) = zero
 1004    continue
         do 1005 i=1,3
            cstmm(i,j) = zero
 1005    continue
         do 1006 i=1,3
            cstbm(i,j) = zero
 1006    continue
         do 1007 i=1,3
            cstmb(i,j) = zero
 1007    continue
 1003 continue
C
C.....CLEAR THE LOCAL CONTRIBUTIONS TO THE STIFFNESS MATRIX
C
      do 1008 j=1,ndof
C         do 1009 i=1,ndof
            kbBB(i,j) = zero
C 1009    continue
C         do 1010 i=1,ndof
            kbMM(i,j) = zero
C 1010    continue
C         do 1011 i=1,ndof
            kbBM(i,j) = zero
C 1011    continue
C         do 1012 i=1,ndof
            kbMB(i,j) = zero
C 1012    continue
C         do 1013 i=1,ndof
            khBB(i,j) = zero
C 1013    continue
C         do 1014 i=1,ndof
            khMM(i,j) = zero
C 1014    continue
C         do 1015 i=1,ndof
            khBM(i,j) = zero
C 1015    continue
C         do 1016 i=1,ndof
            khMB(i,j) = zero
C 1016    continue
 1008 continue
C
C.....CLEAR THE LOCAL CONSTITUTIVE COEFFICIENT ARRAY
C
      do 1017 i=1,36
         cstcoef(i) = zero
 1017 continue
C
C.....CLEAR THE TRIANGULAR COORDINATES
C
      do 1018 i=1,3
         xlp(i) = zero
 1018 continue
C
      do 1019 i=1,3
         ylp(i) = zero
 1019 continue
C
      do 1020 i=1,3
         zlp(i) = zero
 1020 continue
C
C.....CLEAR THE DEGREE OF FREEDOM POINTERS
C
C      do 1023 i=1,9
C         rowb(i) = 0
C 1023 continue
C
C      do 1024 i=1,9
C         colb(i) = 0
C 1024 continue
C
C      do 1025 i=1,9
C         rowm(i) = 0
C 1025 continue
C
C      do 1026 i=1,9
C         colm(i) = 0
C 1026 continue
C
C.....CLEAR THE LOCAL INTEGRATED STRAIN-TO-DISPLACEMENT AND
C.....CURVATURE-TO-DISPLACEMENT MATRICES
C
      do 1027 j=1,3
         do 1028 i=1,9
            Lb(i,j)  = zero
 1028    continue
         do 1029 i=1,9
            Pb(i,j)  = zero
 1029    continue
         do 1030 i=1,9
            Lh1(i,j) = zero
 1030    continue
         do 1031 i=1,9
            Lh2(i,j) = zero
 1031    continue
         do 1032 i=1,9
            Lh3(i,j) = zero
 1032    continue
         do 1033 i=1,9
            Ph1(i,j) = zero
 1033    continue
         do 1034 i=1,9
            Ph2(i,j) = zero
 1034    continue
         do 1035 i=1,9
            Ph3(i,j) = zero
 1035    continue
 1027 continue
C
C.....CLEAR THE NUMBER OF LAYERS OF THE COMPOSITE ELEMENT
C
      nlayer = 0
C
C.....CLEAR THE LOCAL STORAGES FOR LAYER PROPERTIES
C
      do 1036 j=1,maxlayer
         do 1037 i=1,5
            idcmp23(i,j) = 0
 1037    continue
         do 1038 i=1,8
            mtcmp23(i,j) = zero
 1038    continue
 1036 continue
C
C.....CLEAR THE ARBITRARY FRAME OF THE COMPOSITE ELEMENT
C
      do 1039 j=1,3
         do 1040 i=1,3
            aframe(i,j) = zero
 1040    continue
 1039 continue
C
C.....CLEAR THE ELEMENT LEVEL FRAME
C
      do 1041 j=1,3
         do 1042 i=1,3
            eframe(i,j) = zero
 1042    continue
 1041 continue
C
C.....INITIALIZE THE ELEMENT'S CONSTANT THICKNESS
C
      thick = h(1)
C
C.....CHECK THE TYPE OF CONSTITUTIVE LAW
C
      if (      (ctyp.ne.0)
     $     .and.(ctyp.ne.1)
     $     .and.(ctyp.ne.2)
     $     .and.(ctyp.ne.3) ) go to 200
C
C.....CHECK THE ADDRESSING IN STORAGE [CMPCO]
C
      if ( ctyp.eq.1 ) then
         if ( (catt.lt.1).or.(catt.gt.nttco) ) go to 300
      endif
C
C.....CHECK THE ADDRESSING IN ARRAYS [IDLAY] AND [MTLAY]
C
      if ( (ctyp.eq.2).or.(ctyp.eq.3) ) then
         if ( (catt.lt.1).or.(catt.gt.nttly) ) go to 400
      endif
C
C.....CHECK THE ADDRESSING IN ARRAY [CMPFR]
C
      if ( (ctyp.eq.1).or.(ctyp.eq.2).or.(ctyp.eq.3) ) then
         if ( (cfrm.lt.0).or.(cfrm.gt.ncmpfr) ) go to 500
      endif
C
C.....INITIALIZE THE CONSTITUTIVE COEFFICIENTS IN CASE OF A TYPE-1 LAW
C
      if ( ctyp.eq.1 ) then
         do 2001 i=1,36
            cstcoef(i) = cmpco(i,catt)
 2001    continue
      endif
C
C.....INITIALIZE THE LAYER PROPERTIES FOR TYPE-2 AND TYPE-3 LAWS
C
      if ( (ctyp.eq.2).or.(ctyp.eq.3) ) then
         nlayer = idlay(2,catt)
         if ( nlayer.gt.maxlayer ) go to 600
         do 2002 i=catt,(catt+nlayer-1)
            ilayer = idlay(3,i)
            if ( ilayer.gt.nlayer ) go to 700
            idcmp23(1,ilayer) = idlay(1,i)
            idcmp23(2,ilayer) = idlay(2,i)
            idcmp23(3,ilayer) = idlay(3,i)
            idcmp23(4,ilayer) = idlay(4,i)
            idcmp23(5,ilayer) = idlay(5,i)
            mtcmp23(1,ilayer) = mtlay(1,i)
            mtcmp23(2,ilayer) = mtlay(2,i)
            mtcmp23(3,ilayer) = mtlay(3,i)
            mtcmp23(4,ilayer) = mtlay(4,i)
            mtcmp23(5,ilayer) = mtlay(5,i)
            mtcmp23(6,ilayer) = mtlay(6,i)
            mtcmp23(7,ilayer) = mtlay(8,i)
            mtcmp23(8,ilayer) = mtlay(9,i)
 2002    continue
      endif
C
C.....CHECK CONSISTENCY OF BASIC PARAMETERS FOR TYPES 2 AND 3
C
      if ( (ctyp.eq.2).or.(ctyp.eq.3) ) then
         do 2003 i=1,nlayer
            if ( idcmp23(1,i).ne.idcmp23(1,1) ) go to 800
            if ( idcmp23(2,i).ne.nlayer       ) go to 800
 2003    continue
         do 2004 i=(nlayer+1),maxlayer
            if ( idcmp23(1,i).ne.0 ) go to 800
            if ( idcmp23(2,i).ne.0 ) go to 800
 2004    continue
      endif
C
C.....GET THE ARBITRARY FRAME FOR DEFINITION OF THE CONSTITUTIVE LAW
C
      if ( cfrm .eq. 0 ) then
C
C.....INITIALIZE WITH THE IDENTITY IF THE FRAME NUMBER IS ZERO
C
         aframe(1,1) = one
         aframe(2,1) = zero
         aframe(3,1) = zero
         aframe(1,2) = zero
         aframe(2,2) = one
         aframe(3,2) = zero
         aframe(1,3) = zero
         aframe(2,3) = zero
         aframe(3,3) = one
C
      else
C
         aframe(1,1) = cmpfr(1,cfrm)
         aframe(2,1) = cmpfr(2,cfrm)
         aframe(3,1) = cmpfr(3,cfrm)
         aframe(1,2) = cmpfr(4,cfrm)
         aframe(2,2) = cmpfr(5,cfrm)
         aframe(3,2) = cmpfr(6,cfrm)
         aframe(1,3) = cmpfr(7,cfrm)
         aframe(2,3) = cmpfr(8,cfrm)
         aframe(3,3) = cmpfr(9,cfrm)
C
      endif
C
C.....GET THE ELEMENT TRIANGULAR COORDINATES
C.....GET THE ROTATION MATRIX
C.....GET THE DEGREE OF FREEDOM POINTERS
C
      call compcrd( elm  , ctyp , X    , Y    , Z    , rot  ,
     $              xlp  , ylp  , zlp  , rowb , colb , rowm ,
     $              colm                                    )
C
C.....GET THE ELEMENT LEVEL FRAME FROM THE 6 x 6 ROTATION MATRIX
C
      eframe(1,1) = rot(1,1)
      eframe(2,1) = rot(2,1)
      eframe(3,1) = rot(3,1)
      eframe(1,2) = rot(1,2)
      eframe(2,2) = rot(2,2)
      eframe(3,2) = rot(3,2)
      eframe(1,3) = rot(1,3)
      eframe(2,3) = rot(2,3)
      eframe(3,3) = rot(3,3)
C
C.....GET THE CONSTITUTIVE MATRIX FOR PURE BENDING
C
      call compcst( e       , thick   , nu     , cstcoef , nlayer ,
     $              idcmp23 , mtcmp23 , xlp    , ylp     , zlp    ,
     $              cstbb   , ctyp    , eframe , aframe  , "BB"   )
C
C.....GET THE CONSTITUTIVE MATRIX FOR PURE MEMBRANE
C
      call compcst( e       , thick   , nu     , cstcoef , nlayer ,
     $              idcmp23 , mtcmp23 , xlp    , ylp     , zlp    ,
     $              cstmm   , ctyp    , eframe , aframe  , "MM"   )
C
C.....GET THE CONSTITUTIVE MATRIX FOR BENDING-MEMBRANE COUPLING
C
      call compcst( e       , thick   , nu     , cstcoef , nlayer ,
     $              idcmp23 , mtcmp23 , xlp    , ylp     , zlp    ,
     $              cstbm   , ctyp    , eframe , aframe  , "BM"   )
C
C.....GET THE CONSTITUTIVE MATRIX FOR MEMBRANE-BENDING COUPLING
C
      if ( fastcal ) then
C
      cstmb(1,1) = cstbm(1,1)
      cstmb(2,1) = cstbm(1,2)
      cstmb(3,1) = cstbm(1,3)
      cstmb(1,2) = cstbm(2,1)
      cstmb(2,2) = cstbm(2,2)
      cstmb(3,2) = cstbm(2,3)
      cstmb(1,3) = cstbm(3,1)
      cstmb(2,3) = cstbm(3,2)
      cstmb(3,3) = cstbm(3,3)
C
      else
C
      call compcst( e       , thick   , nu     , cstcoef , nlayer ,
     $              idcmp23 , mtcmp23 , xlp    , ylp     , zlp    ,
     $              cstmb   , ctyp    , eframe , aframe  , "MB"   )
C
      endif
C
C.....CHECK THE CONSTITUTIVE MATRIX
C.....(COMMENTED OUT HERE - USED FOR DEBUGGING ONLY)
C
C     call compchk( elm , cstbb , cstmm , cstbm , cstmb )
C
C.....FORM THE LOCAL BASIC STIFFNESS FOR PURE BENDING
C
      call compbBB( elm , ctyp , xlp  , ylp  , cstbb ,
     $              one , zero , one  , rowb , colb  ,
     $              rot , Lb   , kbBB                )
C
C.....FORM THE LOCAL BASIC STIFFNESS FOR PURE MEMBRANE
C
      call compbMM( elm  , ctyp , xlp  , ylp  , cstmm ,
     $              onehalf     , one  , rowm , colm  ,
     $              rot  , Pb   , kbMM                )
C
C.....FORM THE LOCAL BASIC STIFFNESS FOR BENDING-MEMBRANE COUPLING
C
      call compbBM( elm , ctyp , xlp  , ylp  , cstbm ,
     $              one , zero , one  , onehalf      ,
     $              one , rowb , colm , rot  , Lb    ,
     $              Pb  , fastcal     , kbBM         )
C
C.....FORM THE LOCAL BASIC STIFFNESS FOR MEMBRANE-BENDING COUPLING
C
      call compbMB( elm , ctyp , xlp  , ylp  , cstmb ,
     $              one , zero , one  , onehalf      ,
     $              one , rowm , colb , rot  , Lb    ,
     $              Pb  , fastcal     , kbMB         )
C
C.....FORM THE LOCAL HIGHER ORDER STIFFNESS FOR PURE BENDING
C
      call comphBB( elm , ctyp , xlp  , ylp , cstbb ,
     $              one , rowb , colb , rot , Lh1   ,
     $              Lh2 , Lh3  , khBB               )
C
C.....FORM THE LOCAL HIGHER ORDER STIFFNESS FOR PURE MEMBRANE
C
      call comphMM( elm     , ctyp , xlp  , ylp , cstmm ,
     $              point32 , rowm , colm , rot , Ph1   ,
     $              Ph2     , Ph3  , khMM               )
C
C.....FORM THE LOCAL HIGHER ORDER STIFFNESS FOR BENDING-MEMBRANE COUPLING

C
      call comphBM( elm , ctyp    , xlp  , ylp  , cstbm ,
     $              one , point32 , rowb , colm , rot   ,
     $              Lh1 , Lh2     , Lh3  , Ph1  , Ph2   ,
     $              Ph3 , fastcal , khBM                )
C
C.....FORM THE LOCAL HIGHER ORDER STIFFNESS FOR MEMBRANE-BENDING COUPLING

C
      call comphMB( elm , ctyp    , xlp  , ylp  , cstmb ,
     $              one , point32 , rowm , colb , rot   ,
     $              Lh1 , Lh2     , Lh3  , Ph1  , Ph2   ,
     $              Ph3 , fastcal , khMB                )
C
C.....ADD ALL LOCAL STIFFNESS MATRICES
C
      do 3001 j=1,ndof
         do 3002 i=1,ndof
            estiff(i,j) = kbBB(i,j) + kbMM(i,j)
     $                  + kbBM(i,j) + kbMB(i,j)
     $                  + khBB(i,j) + khMM(i,j)
     $                  + khBM(i,j) + khMB(i,j)
 3002    continue
 3001 continue
C
C.....ROTATE ELEMENT STIFFNESS TO GLOBAL COORDINATES
C
C      call compmrot(estiff, rot, rot, rot)
C     
C     if flag equals 1, transform from local to global frame
C     else return local element stiffness matrix.
C
      if(flag .eq. 1) then
        call trirotation( estiff, eframe ) 
      endif
C
C.....CHECK THE POSITIVITY OF THE OUTPUT STIFFNESS MATRIX
C.....(COMMENTED OUT HERE - USED FOR DEBUGGING ONLY)
C
C     call compchk2( elm , estiff )
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
C.....ERROR-MESSAGE IF DIMENSIONS DO NOT AGREE
C
  100 continue
      write(*,*) "*** FATAL ERROR in Routine COMPST     ***"
      write(*,*) "*** The Elemental Stiffness Matrix    ***"
      write(*,*) "*** Must Be a Square 18 by 18 Matrix: ***"
      write(*,*) "*** Dimensions Do Not Check Out...    ***"
      write(*,*) "*** Execution Terminated Here         ***"
      stop
C
C.....ERROR-MESSAGE IF UNKNOWN TYPE OF CONSTITUTIVE LAW
C
  200 continue
      write(*,*) "*** FATAL ERROR in Routine COMPST ***"
      write(*,*) "*** The Type of Constitutive Law  ***"
      write(*,*) "*** Must Either be 0, 1, 2, or 3. ***"
      write(*,*) "*** Execution Terminated Here     ***"
      stop
C
C.....ERROR-MESSAGE IF ADDRESSING IN [CMPCO] IS NOT CORRECT
C
  300 continue
      write(*,*) "*** FATAL ERROR in Routine COMPST ***"
      write(*,*) "*** The Address in Array [CMPCO]  ***"
      write(*,*) "*** is Out-of-Bounds.             ***"
      write(*,*) "*** Execution Terminated Here     ***"
      stop
C
C.....ERROR-MESSAGE IF ADDRESSING IN [IDLAY]/[MTLAY] IS NOT CORRECT
C
  400 continue
      write(*,*) "*** FATAL ERROR in Routine COMPST ***"
      write(*,*) "*** The Address in Arrays [IDLAY] ***"
      write(*,*) "*** and [MTLAY] is Out-of-Bounds. ***"
      write(*,*) "*** Execution Terminated Here     ***"
      stop
C
C.....ERROR-MESSAGE IF ADDRESSING IN [CMPFR] IS NOT CORRECT
C
  500 continue
      write(*,*) "*** FATAL ERROR in Routine COMPST ***"
      write(*,*) "*** The Address in Array [CMPFR]  ***"
      write(*,*) "*** is Out-of-Bounds.             ***"
      write(*,*) "*** Execution Terminated Here     ***"
      stop
C
C.....ERROR-MESSAGE IF THE MAXIMUM NUMBER OF LAYERS HAS BEEN EXCEEDED
C
  600 continue
      write(*,*) "*** FATAL ERROR in Routine COMPST         ***"
      write(*,*) "*** The Maximum Number of Layers Has Been ***"
      write(*,*) "*** Exceeded: Boost Parameter [MAXLAYER]. ***"
      write(*,*) "*** Execution Terminated Here             ***"
      stop
C
C.....ERROR-MESSAGE IF THE TOTAL NUMBER OF LAYERS IS NOT CONSISTENT
C
  700 continue
      write(*,*) "*** FATAL ERROR in Routine COMPST       ***"
      write(*,*) "*** A Layer Number Exceeds the Total    ***"
      write(*,*) "*** Number of Layers Stored in [IDLAY]. ***"
      write(*,*) "*** Execution Terminated Here           ***"
      stop
C
C.....ERROR-MESSAGE IF BASIC PARAMETERS FOR TYPES 2 & 3 ARE INCONSISTENT
C
  800 continue
      write(*,*) "*** FATAL ERROR in Routine COMPST         ***"
      write(*,*) "*** The Element Numbers and Total Number  ***"
      write(*,*) "*** of Composite Layers Do Not Check Out. ***"
      write(*,*) "*** Execution Terminated Here             ***"
      stop
C
C     ---
C     END
C     ---
C
      end
C=end of routine "COMPST"
C=======================C
