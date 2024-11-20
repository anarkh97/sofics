C=====================================================================C
C        1         2         3         4         5         6         7C
C234567890123456789012345678901234567890123456789012345678901234567890C
C=====================================================================C
      subroutine compcst( e       , thick   , nu     , coef   , nlayer,
     $                    idlayer , mtlayer , x      , y      , z     ,
     $                    d       , type    , eframe , aframe , effect)
C=====================================================================C
C                                                                     C
C     Perform =   Assembles the 3 by 3 Constitutive Matrix According  C
C     ---------   to the Type of Constitutive Law Requested.          C
C                                                                     C
C                                                                     C
C     Input/Output =                                                  C
C     --------------                                                  C
C     E       <input>  Young modulus                                  C
C     THICK   <input>  thickness (assumed constant over the element)  C
C     NU      <input>  Poisson's ratio                                C
C     COEF    <input>  coefficients of the constitutive law           C
C     NLAYER  <input>  number of layers of the composite element      C
C     IDLAYER <input>  identificators for each layer                  C
C     MTLAYER <input>  material properties of each layer              C
C     X       <input>  triangular coordinates in local x-direction    C
C     Y       <input>  triangular coordinates in local y-direction    C
C     Z       <input>  triangular coordinates in local z-direction    C
C     D       <output> 3 by 3 constitutive matrix                     C
C     TYPE    <input>  type of constitutive law                       C
C     EFRAME  <input>  element level 3x3 frame                        C
C     AFRAME  <input>  arbitrary 3x3 frame of the constitutive law    C
C     EFFECT  <input>  type of matrix [D] requested                   C
C                                                                     C
C                                                                     C
C     Computations =                                                  C
C     --------------                                                  C
C     Four different matrices [D] are assembled according to whether  C
C     the flag [effect] is equal to "BB", "MM", "BM", or "MB":        C
C                                                                     C
C                             [ [D_mm]  [D_mb] ]                      C
C     [Constitutive_Matrix] = [                ]                      C
C            6 by 6           [ [D_bm]  [D_bb] ]                      C
C                                                                     C
C     where "b" and "m" stand for bending and membrane, respectively: C
C                                                                     C
C     [effect] = "BB"  =>  Assemble [D_bb] in [D]                     C
C     [effect] = "MM"  =>  Assemble [D_mm] in [D]                     C
C     [effect] = "BM"  =>  Assemble [D_bm] in [D]                     C
C     [effect] = "MB"  =>  Assemble [D_mb] in [D]                     C
C                                                                     C
C     The constitutive matrix [D_bb] relates the element's curvatures C
C     to the bending moments and the constitutive matrix [D_mm]       C
C     relates the element's normal efforts to the strains. Similarly, C
C     the constitutive matrices [D_bm] and [D_mb] couple the bending  C
C     and membrane effects:                                           C
C                                                                     C
C     [ M_x  ]            [ k_x  ]                                    C
C     [ M_y  ] = [D_bb] * [ k_y  ]                                    C
C     [ M_xy ]            [ k_xy ]                                    C
C                                                                     C
C     [ N_x  ]            [ e_x  ]                                    C
C     [ N_y  ] = [D_mm] * [ e_y  ]                                    C
C     [ N_xy ]            [ e_xy ]                                    C
C                                                                     C
C     [ M_x  ]            [ e_x  ]                                    C
C     [ M_y  ] = [D_bm] * [ e_y  ]                                    C
C     [ M_xy ]            [ e_xy ]                                    C
C                                                                     C
C     [ N_x  ]            [ k_x  ]                                    C
C     [ N_y  ] = [D_mb] * [ k_y  ]                                    C
C     [ N_xy ]            [ k_xy ]                                    C
C                                                                     C
C     The symmetry constraints are:                                   C
C                                                                     C
C     [D_bb]^T = [D_bb]                                               C
C     [D_mm]^T = [D_mm]                                               C
C     [D_bm]^T = [D_mb]                                               C
C                                                                     C
C     The assembly of these matrices is performed according to the    C
C     type of constitutive law available and prescribed by [type].    C
C     In the following, the assembly of each one of the four matrices C
C     is briefly summarized.                                          C
C                                                                     C
C     1.  Constitutive Law of type-0:                                 C
C     - - - - - - - - - - - - - - - -                                 C
C     The material is assumed isotropic and known via the Young       C
C     modulus [E], the Poisson's ratio [nu] and thickness [thick]:    C
C                                                                     C
C     1.1  Pure Bending ("BB"):                                       C
C     - - - - - - - - - - - - -                                       C
C                                                                     C
C              [ d_11  d_12   0   ]                                   C
C     [D_bb] = [ d_12  d_11   0   ]                                   C
C              [  0     0    d_33 ]                                   C
C                                                                     C
C     with:                                                           C
C                                                                     C
C              [E]*[thick]^3                                          C
C     [d_11] = -------------                                          C
C              12*(1-[nu]^2)                                          C
C                                                                     C
C              [nu]*[E]*[thick]^3                                     C
C     [d_12] = ------------------                                     C
C                12*(1-[nu]^2)                                        C
C                                                                     C
C              [E]*[thick]^3                                          C
C     [d_33] = -------------                                          C
C               24*(1+[nu])                                           C
C                                                                     C
C     1.2 Pure Membrane ("MM"):                                       C
C     - - - - - - - - - - - - -                                       C
C                                                                     C
C              [ d_44  d_45   0   ]                                   C
C     [D_mm] = [ d_45  d_55   0   ]                                   C
C              [  0     0    d_66 ]                                   C
C                                                                     C
C     with:                                                           C
C                                                                     C
C              [E]*[thick]                                            C
C     [d_44] = -----------                                            C
C              (1-[nu]^2)                                             C
C                                                                     C
C              [nu]*[E]*[thick]                                       C
C     [d_45] = ----------------                                       C
C                 (1-[nu]^2)                                          C
C                                                                     C
C              [E]*[thick]                                            C
C     [d_66] = ------------                                           C
C               2*(1+[nu])                                            C
C                                                                     C
C     1.3 Coupling Bending-Membrane ("BM"):                           C
C     - - - - - - - - - - - - - - - - - - -                           C
C     [D_bm] = zero (no coupling for isotropic material)              C
C                                                                     C
C     1.4 Coupling Membrane-Bending ("MB"):                           C
C     - - - - - - - - - - - - - - - - - - -                           C
C     [D_mb] = zero (no coupling for isotropic material)              C
C                                                                     C
C     2.  Constitutive Law of type-1:                                 C
C     - - - - - - - - - - - - - - - -                                 C
C     The coefficients [d_ij] for i=1 to 6 and j=1 to 6 are given in  C
C     the output. They are stored in a vector of length 36 according  C
C     to the following convention:                                    C
C                                                                     C
C     [  d_11  d_12  d_13  d_14  d_15  d_16  ]                        C
C     [  d_12  d_22  d_23  d_24  d_25  d_26  ]                        C
C     [  d_13  d_23  d_33  d_34  d_35  d_36  ]                        C
C     [  d_14  d_24  d_33  d_44  d_45  d_46  ]                        C
C     [  d_15  d_25  d_33  d_44  d_55  d_56  ]                        C
C     [  d_16  d_26  d_33  d_44  d_55  d_66  ]                        C
C                                                                     C
C     is stored in the vector [coef] of size 36 at the following      C
C     location:                                                       C
C                                                                     C
C     [   01    02    03    04    05    06   ]                        C
C     [   07    08    09    10    11    12   ]                        C
C     [   13    14    15    16    17    18   ]                        C
C     [   19    20    21    22    23    24   ]                        C
C     [   25    26    27    28    29    30   ]                        C
C     [   31    32    33    34    35    36   ]                        C
C                                                                     C
C     Therefore, the entry on the ith row and jth column is stored at C
C     position number 6*(i-1)+j in the vector [coef].                 C
C                                                                     C
C     2.1  Pure Bending ("BB"):                                       C
C     - - - - - - - - - - - - -                                       C
C                                                                     C
C              [ d_11 d_12 d_13 ]                    [ 22  23  24 ]   C
C     [D_bb] = [ d_12 d_22 d_23 ] found at positions [ 28  29  30 ]   C
C              [ d_13 d_23 d_33 ]                    [ 34  35  36 ]   C
C                                                                     C
C     2.2 Pure Membrane ("MM"):                                       C
C     - - - - - - - - - - - - -                                       C
C                                                                     C
C              [ d_44 d_42 d_43 ]                    [ 01  02  03 ]   C
C     [D_mm] = [ d_45 d_55 d_56 ] found at positions [ 07  08  09 ]   C
C              [ d_43 d_56 d_66 ]                    [ 13  14  15 ]   C
C                                                                     C
C     2.3 Coupling Bending-Membrane ("BM"):                           C
C     - - - - - - - - - - - - - - - - - - -                           C
C                                                                     C
C              [ d_14 d_15 d_16 ]                    [ 19  20  21 ]   C
C     [D_bm] = [ d_24 d_25 d_26 ] found at positions [ 25  26  27 ]   C
C              [ d_34 d_35 d_36 ]                    [ 31  32  33 ]   C
C                                                                     C
C     2.4 Coupling Membrane-Bending ("MB"):                           C
C     - - - - - - - - - - - - - - - - - - -                           C
C                                                                     C
C              [ d_41 d_42 d_43 ]                    [ 04  05  06 ]   C
C     [D_mb] = [ d_51 d_52 d_53 ] found at positions [ 10  11  12 ]   C
C              [ d_61 d_62 d_63 ]                    [ 16  17  18 ]   C
C                                                                     C
C     3.  Constitutive Law of type-2:                                 C
C     - - - - - - - - - - - - - - - -                                 C
C     The material properties of each layer are known and integration C
C     through the thickness of the composite material is performed.   C
C     It is assumed that there is NO coupling between bending and     C
C     membrane effects (even though these terms are found non-zero    C
C     when numerical integration through the thickness is performed). C
C                                                                     C
C     3.1  Pure Bending ("BB"):                                       C
C     - - - - - - - - - - - - -                                       C
C                                                                     C
C     3.2 Pure Membrane ("MM"):                                       C
C     - - - - - - - - - - - - -                                       C
C                                                                     C
C     3.3 Coupling Bending-Membrane ("BM"):                           C
C     - - - - - - - - - - - - - - - - - - -                           C
C     [D_bm] = zero (assumed)                                         C
C                                                                     C
C     3.4 Coupling Membrane-Bending ("MB"):                           C
C     - - - - - - - - - - - - - - - - - - -                           C
C     [D_mb] = zero (assumed)                                         C
C                                                                     C
C     4.  Constitutive Law of type-3:                                 C
C     - - - - - - - - - - - - - - - -                                 C
C     The material properties of each layer are known and integration C
C     through the thickness of the composite material is performed.   C
C                                                                     C
C     4.1  Pure Bending ("BB"):                                       C
C     - - - - - - - - - - - - -                                       C
C                                                                     C
C     4.2 Pure Membrane ("MM"):                                       C
C     - - - - - - - - - - - - -                                       C
C                                                                     C
C     4.3 Coupling Bending-Membrane ("BM"):                           C
C     - - - - - - - - - - - - - - - - - - -                           C
C                                                                     C
C     4.4 Coupling Membrane-Bending ("MB"):                           C
C     - - - - - - - - - - - - - - - - - - -                           C
C                                                                     C
C                                                                     C
C     Caution =   It is assumed that the element has a constant       C
C     ---------   thickness so that no numerical interpolation is     C
C                 required. It is also assumed that the symmetry of   C
C                 the [D] matrix has been checked when its 36 entries C
C                 are inputed ([type]=1). (See routines "precmp.f"    C
C                 and "reacmp.f" in directory "Input".)               C
C                                                                     C
C=====================================================================C
C=Author  = Francois M. Hemez                                         C
C=Date    = June 9th, 1995                                            C
C=Version = 2.0                                                       C
C=Comment =                                                           C
C=====================================================================C
C
C     ------------
C     DECLARATIONS
C     ------------
C
C.....Global Variables
C
      integer     type , nlayer , idlayer(5,nlayer)
      real*8      e , nu , thick , coef(36) , d(3,3)
      real*8      x(3) , y(3) , z(3) , mtlayer(8,nlayer)
      real*8      eframe(3,3) , aframe(3,3)
      character   effect*2
C
C.....Local Variables
C
      integer     i , j , k , ilayer , layernumber , irot
      real*8      zero , one , pi , twopi , intthick
      real*8      E1 , E2 , nu12 , G12 , mu1 , mu2
      real*8      zsup , zinf , thetaF , thetaD , theta
      real*8      S11 , S12 , S13 , S22 , S23 , S33 , detS
      real*8      Q(3,3) , Qbar(3,3) , T(3,3) , R(3)
      real*8      invT(3,3) , costheta , sintheta
      real*8      z0 , qt1 , qt2 , qt3
      real*8      rotate(3,3) , RotateD(3,3) , refvec(3)
      real*8      norm1 , norm2 , normref , proj1 , proj2
      real*8      cosine1 , cosine2 , orifiber(3)
      logical     pureben , puremem , cbenmem , cmemben
c PJSA 5-11-2007 tensor rotation matrix and inverse
      real*8      rtrmat(3,3), tmatinv(3,3)
      real*8      c,s
C
C     ----
C     DATA
C     ----
C
      data zero   /0.000000D+00/
      data one    /1.000000D+00/
C
C.....INITIALIZE THE MAIN DIAGONAL OF REUTER'S MATRIX
C
      data R      /1.000000D+00, 1.000000D+00, 2.000000D+00/
C
C     -----
C     LOGIC
C     -----
C
c pjsa start new method for type 1 with generalized tensor transformation based on types 2&3
      if( type.eq.1 ) then 
         call compcst1(e, thick, nu, coef, nlayer,
     $                 idlayer, mtlayer, x, y, z,
     $                 d, type, eframe, aframe, effect)
         return
      endif
c pjsa end 
C
      pi    = acos(-one)
      twopi = 2.000000D+00*pi
C
C.....CHECK IF THE ELEMENT IS A COMPOSITE SHELL
C
      if ( type.eq.-1 ) go to 100
C
C.....CHECK IF TYPE OF CONTITUTIVE LAW IS CORRECT
C
      if (      (type.ne.0)
     $     .and.(type.ne.1)
     $     .and.(type.ne.2)
     $     .and.(type.ne.3) ) go to 200
C
C.....CHECK THE PHYSICAL EFFECT
C
      if (      (effect.ne."BB")
     $     .and.(effect.ne."MM")
     $     .and.(effect.ne."BM")
     $     .and.(effect.ne."MB") ) go to 300
C
C.....SET THE TYPE OF PHYSICAL EFFECT
C
      pureben = ( effect.eq."BB" )
      puremem = ( effect.eq."MM" )
      cbenmem = ( effect.eq."BM" )
      cmemben = ( effect.eq."MB" )
C
C.....CLEAR THE 3 BY 3 CONSTITUTIVE MATRIX
C
      do 1001 j=1,3
         do 1002 i=1,3
            d(i,j) = zero
 1002    continue
 1001 continue
C
C     ------------------------------------------------
C     ISOTROPIC (ESSENTIALLY PLANE) MATERIAL
C     NO COUPLING BETWEEN BENDING AND MEMBRANE EFFECTS
C     ------------------------------------------------
C
      if ( type.eq.0 ) then
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE BENDING
C
      if ( pureben ) then
C
         d(1,1) = e*(thick*thick*thick)/(12.00D+00*(one-nu*nu))
         d(1,2) = nu*e*(thick*thick*thick)/(12.00D+00*(one-nu*nu))
         d(1,3) = zero
         d(2,1) = nu*e*(thick*thick*thick)/(12.00D+00*(one-nu*nu))
         d(2,2) = e*(thick*thick*thick)/(12.00D+00*(one-nu*nu))
         d(2,3) = zero
         d(3,1) = zero
         d(3,2) = zero
         d(3,3) = e*(thick*thick*thick)/(24.00D+00*(one+nu))
C
      endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE MEMBRANE
C
      if ( puremem ) then
C
         d(1,1) = e*thick/(one-nu*nu)
         d(1,2) = nu*e*thick/(one-nu*nu)
         d(1,3) = zero
         d(2,1) = nu*e*thick/(one-nu*nu)
         d(2,2) = e*thick/(one-nu*nu)
         d(2,3) = zero
         d(3,1) = zero
         d(3,2) = zero
         d(3,3) = e*thick/(2.00D+00*(one+nu))
C
      endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING BENDING-MEMBRANE
C
      if ( cbenmem ) then
C
         d(1,1) = zero
         d(1,2) = zero
         d(1,3) = zero
         d(2,1) = zero
         d(2,2) = zero
         d(2,3) = zero
         d(3,1) = zero
         d(3,2) = zero
         d(3,3) = zero
C
      endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING MEMBRANE-BENDING
C
      if ( cmemben ) then
C
         d(1,1) = zero
         d(1,2) = zero
         d(1,3) = zero
         d(2,1) = zero
         d(2,2) = zero
         d(2,3) = zero
         d(3,1) = zero
         d(3,2) = zero
         d(3,3) = zero
C
      endif
C
C.....END OF TREATMENT FOR ISOTROPIC MATERIAL
C
      endif
C
C     -------------------------------------------------
C     COMPOSITE MATERIAL WITH KNOWN CONSTITUTIVE MATRIX
C     -------------------------------------------------
C
      if ( type.eq.1 ) then
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE BENDING
C
      if ( pureben ) then
C
         d(1,1) = coef(22)
         d(1,2) = coef(23)
         d(1,3) = coef(24)
         d(2,1) = coef(28)
         d(2,2) = coef(29)
         d(2,3) = coef(30)
         d(3,1) = coef(34)
         d(3,2) = coef(35)
         d(3,3) = coef(36)
C
      endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE MEMBRANE
C
      if ( puremem ) then
C
         d(1,1) = coef( 1)
         d(1,2) = coef( 2)
         d(1,3) = coef( 3)
         d(2,1) = coef( 7)
         d(2,2) = coef( 8)
         d(2,3) = coef( 9)
         d(3,1) = coef(13)
         d(3,2) = coef(14)
         d(3,3) = coef(15)
C
      endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING BENDING-MEMBRANE
C
      if ( cbenmem ) then
C
         d(1,1) = coef(19)
         d(1,2) = coef(20)
         d(1,3) = coef(21)
         d(2,1) = coef(25)
         d(2,2) = coef(26)
         d(2,3) = coef(27)
         d(3,1) = coef(31)
         d(3,2) = coef(32)
         d(3,3) = coef(33)
C
      endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING MEMBRANE-BENDING
C
      if ( cmemben ) then
C
         d(1,1) = coef( 4)
         d(1,2) = coef( 5)
         d(1,3) = coef( 6)
         d(2,1) = coef(10)
         d(2,2) = coef(11)
         d(2,3) = coef(12)
         d(3,1) = coef(16)
         d(3,2) = coef(17)
         d(3,3) = coef(18)
C
      endif
C
C.....COMPUTE THE ROTATION MATRIX FROM THE ARBITRARY FRAME SYSTEM
C.....TO THE ELEMENT LEVEL FRAME SYSTEM ASSUMING THAT EACH ONE OF
C.....THESE FRAMES ARE ORTHONORMAL, RIGHT-HANDED SYSTEMS (THE I-TH
C.....VECTOR, I=1,2 OR 3, OF THE ARBITRARY FRAME SYSTEM IS ROTATED
C.....ONTO THE I-TH VECTOR OF THE ELEMENTAL FRAME SYSTEM)
C
c pjsa debug
c      write(*,*) 'aframe = ', aframe
c      write(*,*) 'eframe = ', eframe
c      write(*,*) 'd before = ', d
      do 1101 j=1,3
         do 1102 i=1,3
            rotate(i,j) = zero
 1102    continue
 1101 continue
C
      do 1103 j=1,3
         do 1104 i=1,3
            do 1105 k=1,3
               rotate(i,j) = rotate(i,j) + eframe(i,k)*aframe(j,k)
 1105       continue
 1104    continue
 1103 continue
c PJSA 5-11-2007: need tensor rotation matrix, including conversion from engineering strains 
c      write(*,*) 'rotate = ',rotate
      c = rotate(1,1)
      s = rotate(1,2)
      rtrmat(1,1) = c*c
      rtrmat(1,2) = s*s
      rtrmat(1,3) = c*s
      rtrmat(2,1) = s*s
      rtrmat(2,2) = c*c
      rtrmat(2,3) = -c*s
      rtrmat(3,1) = -2*c*s
      rtrmat(3,2) = 2*c*s
      rtrmat(3,3) = c*c-s*s
c      write(*,*) 'rtrmat = ',rtrmat
      tmatinv(1,1) = c*c
      tmatinv(1,2) = s*s
      tmatinv(1,3) = -2*c*s
      tmatinv(2,1) = s*s
      tmatinv(2,2) = c*c
      tmatinv(2,3) = 2*c*s
      tmatinv(3,1) = c*s
      tmatinv(3,2) = -c*s
      tmatinv(3,3) = c*c-s*s
c      write(*,*) 'tmatinv = ',tmatinv
C
C.....ROTATE THE CONSTITUTIVE MATRIX IN THE ELEMENT LEVEL FRAME
C
      do 1106 j=1,3
         do 1107 i=1,3
            RotateD(i,j) = zero
 1107    continue
 1106 continue
C
      do 1108 j=1,3
         do 1109 i=1,3
            do 1110 k=1,3
c PJSA 5-11-2007 use tensor rotation matrix
c               RotateD(i,j) = RotateD(i,j) + d(i,k)*rotate(k,j)
               RotateD(i,j) = RotateD(i,j) + d(i,k)*rtrmat(k,j)
 1110       continue
 1109    continue
 1108 continue
C
      do 1111 j=1,3
         do 1112 i=1,3
            d(i,j) = zero
 1112    continue
 1111 continue
C
      do 1113 j=1,3
         do 1114 i=1,3
            do 1115 k=1,3
c PJSA 5-11-2007 use tensor rotation matrix
c               d(i,j) = d(i,j) + rotate(k,i)*RotateD(k,j)
               d(i,j) = d(i,j) + tmatinv(i,k)*RotateD(k,j)
 1115       continue
 1114    continue
 1113 continue
c pjsa debug
c      write(*,*) 'd after = ', d
C
C.....END OF TREATMENT FOR TYPE-1 CONSTITUTIVE LAW
C
      endif
C
C     ----------------------------------------------
C     COMPOSITE MATERIAL WITH KNOWN LAYER PROPERTIES
C         (COMPUTATIONS COMMUN TO TYPES 2 AND 3)
C     ----------------------------------------------
C
      if ( (type.eq.2).or.(type.eq.3) ) then
C
C.....EXIT (SPEED UP COMPUTATIONS) IF COUPLING IS REQUESTED FOR TYPE-2
C.....BECAUSE TYPE-2 IS PRECISELY THE NO COUPLING CASE
C
      if ( type.eq.2 ) then
         if ( cbenmem ) return
         if ( cmemben ) return
      endif
C
C.....CALCULATE THE TOTAL HALF-HEIGHT OF THE LAYER
C
      z0 = zero
C
      do 2001 ilayer=1,nlayer
         z0 = z0 + mtlayer(7,ilayer)
 2001 continue
C
      z0 = -0.50D+00*z0
C
C.....LOOP ON LAYERS OF THE COMPOSITE SHELL ELEMENT
C
      do 2002 ilayer=1,nlayer
C
C.....INITIALIZE THE LAYER MATERIAL PROPERTIES
C
      layernumber = idlayer(3,ilayer)
C
      E1          = mtlayer(1,ilayer)
      E2          = mtlayer(2,ilayer)
      nu12        = mtlayer(3,ilayer)
      G12         = mtlayer(4,ilayer)
      mu1         = mtlayer(5,ilayer)
      mu2         = mtlayer(6,ilayer)
      thetaF      = mtlayer(8,ilayer)
C
C.....CHECK FOR OBVIOUS ERRORS IN THE INPUT DATA
C
      if ( (layernumber.lt.1).or.(layernumber.gt.nlayer) ) go to 400
C
      if (   E1.le.zero ) go to 500
      if (   E2.le.zero ) go to 500
C     if ( nu12.le.zero ) go to 500
      if (  G12.le.zero ) go to 500
      if (  mu1.lt.zero ) go to 500
      if (  mu2.lt.zero ) go to 500
C
C.....TRANSFORM ANGLE IN THE RANGE BETWEEN 0-360      
C
      if ( (thetaF.lt.zero).or.(thetaF.gt.360.00D+00) ) then 
        irot = thetaF / 360.0D0
        thetaF = thetaF - real(irot) * 360.0D0
        if (thetaF.lt.0.0D0) thetaF = 360.0D0 + thetaF
      endif
C
C.....TRANSFORM FROM DEGREE TO RADIAN THE ANGLE BETWEEN THE
C.....REFERENCE ORIENTATION VECTOR AND THE ORIENTATION OF THE FIBERS
C
      thetaF = (pi*thetaF)/180.00D+00
C
C.....SET THE REFERENCE VECTOR FOR ORIENTING THE FIBERS OF THE LAYER
C.....(ALWAYS TAKE THE FIRST VECTOR OF THE FRAME)
C
      refvec(1) = aframe(1,1)
      refvec(2) = aframe(2,1)
      refvec(3) = aframe(3,1)
C
C.....INITIALIZE THE UPPER AND LOWER [z] COORDINATES FOR THE LAYER
C
      zinf = zero
      zsup = zero
C
      do 2003 i=1,nlayer
         if ( idlayer(3,i).lt.layernumber ) then
            zinf = zinf + mtlayer(7,i)
         endif
         if ( idlayer(3,i).le.layernumber ) then
            zsup = zsup + mtlayer(7,i)
         endif
 2003 continue
C
      zinf = z0 + zinf
      zsup = z0 + zsup
C
C.....PROJECT THE REFERENCE VECTOR INTO THE PLANE OF
C.....THE ELEMENT TO GET THE FIBER ORIENTATION VECTOR
C
      orifiber(1) = zero
      orifiber(2) = zero
      orifiber(3) = zero
C
      norm1   = zero
      norm2   = zero
      normref = zero
      proj1   = zero
      proj2   = zero
C
      do 2004 i=1,3
         norm1   =   norm1 + eframe(i,1)*eframe(i,1)
         norm2   =   norm2 + eframe(i,2)*eframe(i,2)
         normref = normref +   refvec(i)*refvec(i)
         proj1   =   proj1 + eframe(i,1)*refvec(i)
         proj2   =   proj2 + eframe(i,2)*refvec(i)
 2004 continue
C
      norm1   = sqrt(norm1)
      norm2   = sqrt(norm2)
      normref = sqrt(normref)
C
      if ( normref.eq.zero ) then
         cosine1 = one
         cosine2 = zero
      else
         cosine1 = proj1/(norm1*normref)
         cosine2 = proj2/(norm2*normref)
      endif
C
      do 2005 i=1,3
         orifiber(i) = cosine1*eframe(i,1) + cosine2*eframe(i,2)
 2005 continue
C
C.....CALCULATE THE ANGLE FROM THE [x] FRAME OF THE LOCAL
C.....COORDINATE SYSTEM TO THE REFERENCE ORIENTATION VECTOR
C
      proj1  = zero
      proj2  = zero
      thetaD = zero
C
      do 2006 i=1,3
         proj1 = proj1 + eframe(i,1)*orifiber(i)
         proj2 = proj2 + eframe(i,2)*orifiber(i)
 2006 continue
C
      if ( proj1.eq.zero ) then
         if ( proj2.eq.zero ) go to 700
         if ( proj2.gt.zero ) thetaD = 0.50D+00*pi
         if ( proj2.lt.zero ) thetaD = 1.50D+00*pi
      else
         if ( proj2.eq.zero ) then
            if ( proj1.eq.zero ) go to 700
            if ( proj1.gt.zero ) thetaD = zero
            if ( proj1.lt.zero ) thetaD = pi
         else
            thetaD = atan(proj2/proj1)
         endif
      endif
C
      if ( thetaD.lt.zero ) then
         thetaD = thetaD + twopi
      endif
C
      if ( (thetaD.lt.zero).or.(thetaD.gt.twopi) ) go to 800
C
C.....CALCULATE THE ANGLE FROM THE [x] FRAME OF THE LOCAL
C.....COORDINATE SYSTEM TO THE DIRECTION OF THE FIBER
C
      theta = thetaD + thetaF
C
      if ( theta.gt.twopi ) then
         theta = theta - twopi
      endif
C
C.....CALCULATE THE COMPLIANCE MATRIX [S] WHICH RELATES THE STRESSES
C.....[s1], [s2] AND [s12] TO THE STRAINS [e1], [e2] AND [e12] IN
C.....THE COORDINATE SYSTEM {1;2} ASSOCIATED WITH THE FIBER ORIENTATION
C
      S11 = one/E1
      S12 = -nu12/E1
      S13 = mu1/G12
      S22 = one/E2
      S23 = mu2/G12
      S33 = 1/G12
C
C.....CALCULATE THE DETERMINANT OF THE COMPLIANCE MATRIX
C
      detS = S33*(S11*S22-S12*S12) - S11*S23*S23 - S22*S13*S13
      detS = detS + 2.00D+00*S12*S13*S23
C
      if ( detS.eq.zero ) go to 900
C
C.....CALCULATE THE INVERSE OF THE COMPLIANCE MATRIX WHICH RELATES
C.....THE STRAINS [e1], [e2] AND [e12] TO THE STRESSES [s1], [s2]
C.....AND [s12] IN THE COORDINATE SYSTEM {1;2} OF THE FIBER ORIENTATION
C
      do 2202 j=1,3
         do 2203 i=1,3
            Q(i,j) = zero
 2203    continue
 2202 continue
C
      Q(1,1) = S22*S33 - S23*S23
      Q(1,2) = S13*S23 - S12*S33
      Q(1,3) = S12*S23 - S13*S22
C
      Q(2,1) = Q(1,2)
      Q(2,2) = S11*S33 - S13*S13
      Q(2,3) = S12*S13 - S11*S23
C
      Q(3,1) = Q(1,3)
      Q(3,2) = Q(2,3)
      Q(3,3) = S11*S22 - S12*S12
C
      do 2204 j=1,3
         do 2205 i=1,3
            Q(i,j) = Q(i,j)/detS
 2205    continue
 2204 continue
C
C.....INITIALIZE THE ROTATION MATRIX FROM THE FIBER COORDINATE
C.....SYSTEM {1;2} TO THE ELEMENT TRIANGULAR SYSTEM {x;y}
C
      do 2206 j=1,3
         do 2207 i=1,3
            T(i,j) = zero
 2207    continue
 2206 continue
C
      costheta = cos(theta)
      sintheta = sin(theta)
C
      T(1,1)   =  costheta*costheta
      T(1,2)   =  sintheta*sintheta
      T(1,3)   =  2.00D+00*costheta*sintheta
C
      T(2,1)   =  sintheta*sintheta
      T(2,2)   =  costheta*costheta
      T(2,3)   = -2.00D+00*costheta*sintheta
C
      T(3,1)   = -costheta*sintheta
      T(3,2)   =  costheta*sintheta
      T(3,3)   =  costheta*costheta - sintheta*sintheta
C
C.....COMPUTE THE INVERSE OF [T]:
C.....[invT] = inverse(diag[R]) * [T]^t * diag[R]
C
      do 2208 j=1,3
         do 2209 i=1,3
            invT(i,j) = zero
 2209    continue
 2208 continue
C
      do 2210 j=1,3
         do 2211 i=1,3
            invT(i,j) = T(j,i)*(R(j)/R(i))
 2211    continue
 2210 continue
C
C.....CALCULATE THE INVERSE OF THE COMPLIANCE MATRIX WHICH RELATES
C.....THE STRAINS [ex], [ey] AND [exy] TO THE STRESSES [sx], [sy]
C.....AND [sxy] IN THE TRIANGULAR COORDINATE SYSTEM {x;y}:
C.....[Qbar] = [invT] * [Q] * [invT]^t
C
      do 2212 j=1,3
         do 2213 i=1,3
            Qbar(i,j) = zero
 2213    continue
 2212 continue
C
      do 2214 j=1,3
C
         qt1 = Q(1,1)*invT(j,1) + Q(1,2)*invT(j,2) + Q(1,3)*invT(j,3)
         qt2 = Q(2,1)*invT(j,1) + Q(2,2)*invT(j,2) + Q(2,3)*invT(j,3)
         qt3 = Q(3,1)*invT(j,1) + Q(3,2)*invT(j,2) + Q(3,3)*invT(j,3)
C
         do 2215 i=1,3
            Qbar(i,j) = qt1*invT(i,1) + qt2*invT(i,2) + qt3*invT(i,3)
 2215    continue
C
 2214 continue
C
C     ------------------------------------------------
C      COMPOSITE MATERIAL WITH KNOWN LAYER PROPERTIES
C     NO COUPLING BETWEEN BENDING AND MEMBRANE EFFECTS
C      (NUMERICAL INTEGRATION THROUGH THE THICKNESS)
C     ------------------------------------------------
C
      if ( type.eq.2 ) then
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE BENDING
C
      if ( pureben ) then
C
         intthick = (zsup*zsup*zsup - zinf*zinf*zinf)/3.00D+00
C
         do 3001 j=1,3
            do 3002 i=1,3
               d(i,j) = d(i,j) + Qbar(i,j)*intthick
 3002       continue
 3001    continue
C
      endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE MEMBRANE
C
      if ( puremem ) then
C
         intthick = zsup - zinf
C
         do 3003 j=1,3
            do 3004 i=1,3
               d(i,j) = d(i,j) + Qbar(i,j)*intthick
 3004       continue
 3003    continue
C
      endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING BENDING-MEMBRANE
C
      if ( cbenmem ) then
C
         do 3005 j=1,3
            do 3006 i=1,3
               d(i,j) = zero
 3006       continue
 3005    continue
C
      endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING MEMBRANE-BENDING
C
      if ( cmemben ) then
C
         do 3007 j=1,3
            do 3008 i=1,3
               d(i,j) = zero
 3008       continue
 3007    continue
C
      endif
C
C.....END OF TREATMENT FOR TYPE-2 CONSTITUTIVE LAW
C
      endif
C
C     --------------------------------------------------
C       COMPOSITE MATERIAL WITH KNOWN LAYER PROPERTIES
C     WITH COUPLING BETWEEN BENDING AND MEMBRANE EFFECTS
C       (NUMERICAL INTEGRATION THROUGH THE THICKNESS)
C     --------------------------------------------------
C
      if ( type.eq.3 ) then
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE BENDING
C
      if ( pureben ) then
C
         intthick = (zsup*zsup*zsup - zinf*zinf*zinf)/3.00D+00
C
         do 4001 j=1,3
            do 4002 i=1,3
               d(i,j) = d(i,j) + Qbar(i,j)*intthick
 4002       continue
 4001    continue
C
      endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR PURE MEMBRANE
C
      if ( puremem ) then
C
         intthick = zsup - zinf
C
         do 4003 j=1,3
            do 4004 i=1,3
               d(i,j) = d(i,j) + Qbar(i,j)*intthick
 4004       continue
 4003    continue
C
      endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING BENDING-MEMBRANE
C.....(THE TRANSPOSED OF MATRIX [Qbar] IS TAKEN TO GET BEND-MEMB COUPLING)
C
      if ( cbenmem ) then
C
         intthick = 0.50D+00*(zsup*zsup - zinf*zinf)
C
         do 4005 j=1,3
            do 4006 i=1,3
               d(i,j) = d(i,j) + Qbar(j,i)*intthick
 4006       continue
 4005    continue
C
      endif
C
C.....ASSEMBLE THE CONSTITUTIVE MATRIX FOR COUPLING MEMBRANE-BENDING
C
      if ( cmemben ) then
C
         intthick = 0.50D+00*(zsup*zsup - zinf*zinf)
C
         do 4007 j=1,3
            do 4008 i=1,3
               d(i,j) = d(i,j) + Qbar(i,j)*intthick
 4008       continue
 4007    continue
C
      endif
C
C.....END OF TREATMENT FOR TYPE-3 CONSTITUTIVE LAW
C
      endif
C
C.....END OF TREATMENT FOR TYPES 2 AND 3 CONSTITUTIVE LAW
C
 2002 continue
C
      endif
C
C     ------
C     RETURN
C     ------
C
      return
C
C     ------
C     FORMAT
C     ------
C
   91 format("*** Type Used is: ",I10,12x," ***")
   92 format("*** Type Used is: ",5x,A2,9x," ***")
C
C     ---------------
C     ERROR TREATMENT
C     ---------------
C
C.....ERROR-MESSAGE IF THE ELEMENT IS NOT A COMPOSITE SHELL
C
  100 continue
      write(*,*) "*** FATAL ERROR in Routine COMPCST        ***"
      write(*,*) "*** The Finite Element is not a Composite ***"
      write(*,*) "*** Shell Element: Inconsistency Detected ***"
      write(*,*) "*** STOP ALL TREATMENTS RIGHT HERE        ***"
      stop
C
C.....ERROR-MESSAGE IF THE CONSTITUTIVE LAW IS INCORRECT
C
  200 continue
      write(*,*) "*** FATAL ERROR in Routine COMPCST       ***"
      write(*,*) "*** Wrong Type of Constitutive Law       ***"
      write(*,91) type
      write(*,*) "*** Types Allowed are:                   ***"
      write(*,*) "*** 0 = isotropic material (default)     ***"
      write(*,*) "***     (no coupling bending/membrane)   ***"
      write(*,*) "*** 1 = given constitutive law           ***"
      write(*,*) "*** 2 = given layers properties          ***"
      write(*,*) "***     (no coupling bending/membrane)   ***"
      write(*,*) "*** 3 = given layers properties          ***"
      write(*,*) "***     (with coupling bending/membrane) ***"
      write(*,*) "*** STOP ALL TREATMENTS RIGHT HERE       ***"
      stop
C
C.....ERROR-MESSAGE IF THE PHYSICAL EFFECT IS INCORRECT
C
  300 continue
      write(*,*) "*** FATAL ERROR in Routine COMPCST ***"
      write(*,*) "*** Wrong Type of Physical Effect  ***"
      write(*,92) effect
      write(*,*) "*** Types Allowed are:             ***"
      write(*,*) "*** BB = pure bending              ***"
      write(*,*) "*** MM = pure membrane             ***"
      write(*,*) "*** BM = coupling bending-membrane ***"
      write(*,*) "*** MB = coupling membrane-bending ***"
      write(*,*) "*** STOP ALL TREATMENTS RIGHT HERE ***"
      stop
C
C.....ERROR-MESSAGE IF A LAYER IDENTIFICATION NUMBER IS NOT CORRECT
C
  400 continue
      write(*,*) "*** FATAL ERROR in Routine COMPCST  ***"
      write(*,*) "*** The Local Layer Number is Not   ***"
      write(*,*) "*** Correct or Out of Bounds        ***"
      write(*,*) "*** STOP ALL TREATMENTS RIGHT HERE  ***"
      stop
C
C.....ERROR-MESSAGE IF A LAYER PROPERTY IS NOT CORRECT
C
  500 continue
      write(*,*) "*** FATAL ERROR in Routine COMPCST     ***"
      write(*,*) "*** One of the Layer Material Property ***"
      write(*,*) "*** is Negative or Zero!               ***"
      write(*,*) "*** STOP ALL TREATMENTS RIGHT HERE     ***"
      stop
C
C.....ERROR-MESSAGE IF THE FIBER ANGLE IS OUT-OF-BOUNDS
C
  600 continue
      write(*,*) "*** FATAL ERROR in Routine COMPCST    ***"
      write(*,*) "*** The Angle From the Reference      ***"
      write(*,*) "*** Direction to the Direction of     ***"
      write(*,*) "*** Fibers is Out-of-Bounds: it Must  ***"
      write(*,*) "*** be Within the Range 0-360 Degrees ***"
      write(*,*) "*** STOP ALL TREATMENTS RIGHT HERE    ***"
      stop
C
C.....ERROR-MESSAGE IF THE REFERENCE ORIENTATION IS BUGGY
C
  700 continue
      write(*,*) "*** FATAL ERROR in Routine COMPCST   ***"
      write(*,*) "*** The Reference Orientation Vector ***"
      write(*,*) "*** is Parallel to the Two Inplane   ***"
      write(*,*) "*** and Orthogonal Local Frames!     ***"
      write(*,*) "*** STOP ALL TREATMENTS RIGHT HERE   ***"
      stop
C
C.....ERROR-MESSAGE IF THE ORIENTATION ANGLE IS OUT-OF-BOUNDS
C
  800 continue
      write(*,*) "*** FATAL ERROR in Routine COMPCST    ***"
      write(*,*) "*** The Angle From the Local [x]      ***"
      write(*,*) "*** Axis of the Triangular Coordinate ***"
      write(*,*) "*** System to the Reference Direction ***"
      write(*,*) "*** is Out-of-Bounds: it Must be      ***"
      write(*,*) "*** Within the Range 0-2pi Radians    ***"
      write(*,*) "*** STOP ALL TREATMENTS RIGHT HERE    ***"
      stop
C
C.....ERROR-MESSAGE IF THE COMPLIANCE MATRIX IS SINGULAR
C
  900 continue
      write(*,*) "*** FATAL ERROR in Routine COMPCST    ***"
      write(*,*) "*** The Compliance Matrix is Singular ***"
      write(*,*) "*** ... Check Material Properties ... ***"
      write(*,*) "*** STOP ALL TREATMENTS RIGHT HERE    ***"
      stop
C
C     ---
C     END
C     ---
C
      end
C=end of routine "COMPCST"
C========================C
