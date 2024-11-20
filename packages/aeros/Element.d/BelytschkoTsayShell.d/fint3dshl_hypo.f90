! ==================================
! internal force: bt shell / hypo-elastic
! ==================================
!      type                  name                              argument
!      ----                  ----                              --------
! 1.  subroutine        updstrsbt              (delt,ematpro,gqpoin,gqweigt,locbvec,ecurn,evelo,
!                                               effpstrn,hardvar,sigvoitloc,strnvoitloc, effstrs)
! 2.  subroutine        updstrn2bt             (delt,ematpro,zeta,ecurnloc,eveloloc, effpstrn,hardvar,ipstrn,ipstrs, effstrs)
! 3.  subroutine        gqfintbt               (delt,ematpro,gqpoin,gqweigt,area,sigvoitloc,bmat1pt,bcmat1pt,bsmat1pt, gqfint)
! 4.  subroutine        getsighypobt1          (optdmg,delt,ematpro,ehleng,zeta,ecurnloc,eveloloc, damage,ipstrn, ipstrsdot)
! 5.  subroutine        getsighypobt2          (ematpro,ecurnloc,eveloloc, tsstrsdot)
! 6.  subroutine        getstrndotbt1          (thick,zeta,ecurnloc,eveloloc, ltens2d,ipstrndot)
! 7.  subroutine        getstrndotbt2          (ecurnloc,eveloloc, tsstrndot)
!
! =========================================================================================================



subroutine updstrsbt(optctv,optdmg,delt,ematpro,area,ipstrndot,tsstrndot, &
                     effpstrn,hardvar,sigvoitloc,strnvoitloc)
  !=======================================================================
  !  updstrsbt = update j2 (or hypoelas) cauchy stress of belytschko tsay element 
  !
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  delt : time increment
  ! 
  !  ematpro(*) : material property
  !
  !  inoutput:
  !  --------
  !  effpstrn : effective plastic strain (j2) or effective strain (hypoelas)
  !
  !  hardvar : hardening variable (j2) or damage (hypoelas)
  !
  !  sigvoitloc(6,1) : co-rotational cauchy stress
  !
  !  strnvoitloc(6,1) : co-rotational strain
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optctv, optdmg
  real(8), intent(in) :: delt
  real(8), dimension(*), intent(in) :: ematpro
  real(8), intent(in) :: area
  real(8), dimension(3,1), intent(in) :: ipstrndot
  real(8), dimension(2,1), intent(in) :: tsstrndot

  ! ------------------------------------
  real(8), intent(inout) :: effpstrn,hardvar
  real(8), dimension(6,1), intent(inout) :: sigvoitloc, strnvoitloc
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(3,1) :: ipstrnloc, ipstrsloc
  real(8), dimension(2,1) :: tsstrslocdot
  real(8) :: ehleng
  ! ====================================

  ! -------------------------------------------------------------------------------
  ! update in plane stress 
  ! -----------------------------------------------

  ! extract in plane strain components
  ipstrnloc(1,1)= strnvoitloc(1,1)
  ipstrnloc(2,1)= strnvoitloc(2,1)
  ipstrnloc(3,1)= strnvoitloc(6,1)

  ! extract in plane stress components
  ipstrsloc(1,1)= sigvoitloc(1,1)
  ipstrsloc(2,1)= sigvoitloc(2,1)
  ipstrsloc(3,1)= sigvoitloc(6,1)

  ! update hypoelastic in-plane stress of belytschko tsay shell element
  if(optctv .eq. 1) then
    ehleng= dsqrt(area)
    call getsighypobt1(optdmg,delt,ematpro,ehleng,ipstrndot,ipstrnloc, &
                       hardvar,ipstrsloc)
     ! input : optdmg,delt,ematpro,ehleng,ipstrnloc
     ! inoutput : damage,ipstrsloc
  end if

  ! set updated in plane stress components
  sigvoitloc(1,1)= ipstrsloc(1,1) ! sig_x
  sigvoitloc(2,1)= ipstrsloc(2,1) ! sig_y
  sigvoitloc(6,1)= ipstrsloc(3,1) ! sig_xy

  ! -------------------------------------------------------------------------------
  ! compute tranverse shear stress rate and update
  ! ----------------------------------------------
  ! compute hypoelastic transverse shear stress rate of belytschko tsay shell element
  call getsighypobt2(ematpro,tsstrndot, tsstrslocdot)
     ! input : ematpro,tsstrndot
     ! output : tsstrslocdot

  ! set updated transverse stress components
  sigvoitloc(3,1)= 0.0d0 ! sig_z (plane stress)
  sigvoitloc(4,1)= sigvoitloc(4,1) + tsstrslocdot(2,1) * delt ! sig_yz
  sigvoitloc(5,1)= sigvoitloc(5,1) + tsstrslocdot(1,1) * delt ! sig_xz
  ! -------------------------------------------------------------------------------
  
  return
end subroutine updstrsbt


subroutine updstrnbt(optcor,delt,ematpro,gqpoin,eveloloc,bmat1pt,bcmat1pt,bsmat1pt, &
                     strnvoitloc, &
                     ipstrndot,tsstrndot) 
  !=======================================================================
  !  updstrnbt = update strain of belytschko tsay element 
  !
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  delt : time increment
  ! 
  !  ematpro(*) : material property
  !
  !  gqpoin, gqweigt : through thickness gq point and weight
  !
  !  eveloloc(5,4) : local nodal velocity: v_x, v_y, v_z, theta_x, theta_y
  !
  !  inoutput:
  !  --------
  !  strnvoitloc(6,1) : co-rotational strain
  !
  !  output:
  !  ------
  !  ipstrndot(3,1), tsstrndot(2,1) : in plane and transverse shear rate of deformation
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, dimension(2), intent(in) :: optcor
  real(8), intent(in) :: delt
  real(8), dimension(*), intent(in) :: ematpro
  real(8), intent(in) :: gqpoin
  real(8), dimension(6,4), intent(in) :: eveloloc
  real(8), dimension(2,4), intent(in) :: bmat1pt
  real(8), dimension(2,4), intent(in) :: bcmat1pt
  real(8), dimension(2,3,4), intent(in) :: bsmat1pt

  ! ------------------------------------
  real(8), dimension(6,1), intent(inout) :: strnvoitloc
  ! ------------------------------------

  real(8), dimension(3,1), intent(out) :: ipstrndot
  real(8), dimension(2,1), intent(out) :: tsstrndot
  ! ====================================

  ! ----------------------------------------------------------------
  ! compute in plane rate of deformation
  ! ---------------------------
  ! ipstrndot= [d_x, d_y, 2d_xy]
  call getstrndotbt1(optcor,ematpro(20),gqpoin,eveloloc,bmat1pt,bcmat1pt, ipstrndot)
     ! input : optcor,thick,zeta,eveloloc,bmat1pt,bcmat1pt
     ! output : ipstrndot

  ! -------------------------------------------------------------------------------
  ! update in plane strain
  ! -----------------------------------------------
  strnvoitloc(1,1)= strnvoitloc(1,1) + delt*ipstrndot(1,1) ! eps_x
  strnvoitloc(2,1)= strnvoitloc(2,1) + delt*ipstrndot(2,1) ! eps_y
  strnvoitloc(6,1)= strnvoitloc(6,1) + delt*ipstrndot(3,1) ! eps_xy

  ! ----------------------------------------------------------------
  ! compute transverse shear rate of deformation
  ! ---------------------------
  ! tsstrndot= [2d_xz, 2dyz]
  call getstrndotbt2(optcor,eveloloc,bmat1pt,bsmat1pt, tsstrndot)
     ! input : optcor,eveloloc,bmat1pt,bsmat1pt
     ! output : tsstrndot

  ! -------------------------------------------------------------------------------
  ! update transverse shear strain
  ! -----------------------------------------------
  ! note: previously this wasn't ever being computed TODO check
  strnvoitloc(3,1)= 0.d0 ! eps_z
  strnvoitloc(4,1)= strnvoitloc(4,1) + delt*tsstrndot(1,1) ! eps_xz
  strnvoitloc(5,1)= strnvoitloc(5,1) + delt*tsstrndot(2,1) ! eps_yz
  
  return
end subroutine updstrnbt



subroutine gqfintbt(optcor,delt,ematpro,gqpoin,gqweigt,area,sigvoitloc, &
                    bmat1pt,bcmat1pt,bsmat1pt, efintloc)
  !=======================================================================
  !  qfintbt = add local internal force of belytschko tsay shell element
  !            based on qfintbtopt but with some optimizations 
  !
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  delt : time increment
  ! 
  !  ematpro(*) : material property
  !
  !  gqpoin, gqweigt : through thickness gq point and weight
  !
  !  ecurnloc(3,4) : local current nodal coordinate
  !
  !  sigvoitloc(6,1) : voight form of local cauchy stress
  !
  !  output:
  !  ------
  !  efint(24,1) : element nodal force
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, dimension(2), intent(in) :: optcor
  real(8), intent(in) :: delt
  real(8), dimension(*), intent(in) :: ematpro
  real(8), intent(in) :: gqpoin, gqweigt
  real(8), intent(in) :: area
  real(8), dimension(6,1), intent(in) :: sigvoitloc
  real(8), dimension(2,4), intent(in) :: bmat1pt
  real(8), dimension(2,4), intent(in) :: bcmat1pt
  real(8), dimension(2,3,4), intent(in) :: bsmat1pt

  real(8), dimension(24,1), intent(inout) :: efintloc
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: rk, thick
  real(8) :: psibar, dtdpsi
  real(8) :: fx, fy, fxy, fxz, fyz
  real(8) :: mx, my, mxy

  integer :: iloc

  ! loop index
  integer :: inode
  ! ====================================

  ! ------------------------------------------------
  ! material properties
  rk= ematpro(19)     ! shear correction factor: 0.840d0 
  thick= ematpro(20)

  ! -----------------------------------------------------------
  ! pseudo-thickness : eq (5)
  psibar= (gqpoin*thick)/2.0d0

  ! gqweight * d thickness / d psi
  dtdpsi= (gqweigt*thick)/2.0d0

  ! in plane stress
  fx=  sigvoitloc(1,1) *dtdpsi ! f_x
  fy=  sigvoitloc(2,1) *dtdpsi ! f_y
  fxy= sigvoitloc(6,1) *dtdpsi ! f_xy

  ! tranverse shear stress
  fxz= sigvoitloc(5,1) *dtdpsi ! f_xz
  fyz= sigvoitloc(4,1) *dtdpsi ! f_yz

  ! in plane moment
  mx=  psibar * fx ! m_x
  my=  psibar * fy ! m_y
  mxy= psibar * fxy ! m_xy

  ! -----------------------------------------------------------
  ! compute local internal nodal force
  ! ----------------------------------

  do inode=1, 4

    ! address
    iloc= 6*(inode-1)

    ! compute internal nodal force

    if((optcor(1) > 0).and.(optcor(2) > 0)) then
      ! --------------------------------------------------------------------------------------
      ! the original BLT element with warping and shear corrections
      ! -----------------------------------------------------------
      efintloc(iloc+1,1)= efintloc(iloc+1,1) + area*( bmat1pt(1,inode)*fx + bmat1pt(2,inode)*fxy &
                          + bcmat1pt(1,inode)*mx + bcmat1pt(2,inode)*mxy )  ! f_xi

      efintloc(iloc+2,1)= efintloc(iloc+2,1) + area*( bmat1pt(2,inode)*fy + bmat1pt(1,inode)*fxy &
                          + bcmat1pt(2,inode)*my + bcmat1pt(1,inode)*mxy )  ! f_yi

      efintloc(iloc+3,1)= efintloc(iloc+3,1) + area*rk*( bsmat1pt(1,1,inode)*fxz + bsmat1pt(2,1,inode)*fyz ) ! f_zi

      efintloc(iloc+4,1)= efintloc(iloc+4,1) + area*( rk*bsmat1pt(1,2,inode)*fxz + rk*bsmat1pt(2,2,inode)*fyz &
                          - bmat1pt(2,inode)*my - bmat1pt(1,inode)*mxy)  ! m_xi

      efintloc(iloc+5,1)= efintloc(iloc+5,1) + area*( rk*bsmat1pt(1,3,inode)*fxz + rk*bsmat1pt(2,3,inode)*fyz &
                          + bmat1pt(1,inode)*mx + bmat1pt(2,inode)*mxy)  ! m_yi

      efintloc(iloc+6,1)= efintloc(iloc+6,1) + 0.0d0  ! m_zi
 
    else if(optcor(1) > 0) then
      ! --------------------------------------------------------------------------------------
      ! the original BLT element with warping correction
      ! -----------------------------------------------------------
      efintloc(iloc+1,1)= efintloc(iloc+1,1) + area*( bmat1pt(1,inode)*fx + bmat1pt(2,inode)*fxy &
                          + bcmat1pt(1,inode)*mx + bcmat1pt(2,inode)*mxy )  ! f_xi

      efintloc(iloc+2,1)= efintloc(iloc+2,1) + area*( bmat1pt(2,inode)*fy + bmat1pt(1,inode)*fxy &
                          + bcmat1pt(2,inode)*my + bcmat1pt(1,inode)*mxy )  ! f_yi

      efintloc(iloc+3,1)= efintloc(iloc+3,1) + area*rk*( bmat1pt(1,inode)*fxz + bmat1pt(2,inode)*fyz ) ! f_zi

      efintloc(iloc+4,1)= efintloc(iloc+4,1) + area*( - bmat1pt(2,inode)*my - bmat1pt(1,inode)*mxy - 0.250d0*rk*fyz ) ! m_xi

      efintloc(iloc+5,1)= efintloc(iloc+5,1) + area*( bmat1pt(1,inode)*mx + bmat1pt(2,inode)*mxy + 0.250d0*rk*fxz )  ! m_yi

      efintloc(iloc+6,1)= efintloc(iloc+6,1) + 0.0d0  ! m_zi

    else
      ! --------------------------------------------------------------------------------------
      ! the original BLT element
      ! ---------------------------
      efintloc(iloc+1,1)= efintloc(iloc+1,1) + area*( bmat1pt(1,inode)*fx + bmat1pt(2,inode)*fxy )  ! f_xi

      efintloc(iloc+2,1)= efintloc(iloc+2,1) + area*( bmat1pt(2,inode)*fy + bmat1pt(1,inode)*fxy )  ! f_yi

      efintloc(iloc+3,1)= efintloc(iloc+3,1) + area*rk*( bmat1pt(1,inode)*fxz + bmat1pt(2,inode)*fyz ) ! f_zi

      efintloc(iloc+4,1)= efintloc(iloc+4,1) + area*( - bmat1pt(2,inode)*my - bmat1pt(1,inode)*mxy - 0.250d0*rk*fyz ) ! m_xi

      efintloc(iloc+5,1)= efintloc(iloc+5,1) + area*( bmat1pt(1,inode)*mx + bmat1pt(2,inode)*mxy + 0.250d0*rk*fxz )  ! m_yi

      efintloc(iloc+6,1)= efintloc(iloc+6,1) + 0.0d0  ! m_zi

    end if

  end do

  return
end subroutine gqfintbt



subroutine getsighypobt1(optdmg,delt,ematpro,ehleng,ipstrndot,ipstrn, &
                         damage,ipstrs)
  !=======================================================================
  !  getsighypobt1 = compute rate of belytschko-tsay local element cauchy stress 
  !
  !                  note:
  !                  ----
  !                  belytschko, wong and chiang, CMAME, 1992, vol. 96, pp. 93-107
  !                  advances in one point quadrature shell elements
  !                  (see, apendix a.3 )
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optdmg : damage model option handler
  !
  !  delt : integration time step
  !
  !  ematpro(*) : material property
  !
  !  ehleng : characteristic element length
  !
  !  inoutput:
  !  --------
  !  damage : material damage parameter
  !
  !  ipstrn(3,1) : in plane strain
  !                ipstrn(3,1)= strn_x
  !                ipstrn(3,2)= strn_y
  !                ipstrn(3,3)= strn_xy
  !
  !  output:
  !  ------
  !  ipstrsdot(3,1) : in plane stress rate
  !                   ipstrsdot(1,1)= strs_x dot
  !                   ipstrsdot(2,1)= strs_y dot
  !                   ipstrsdot(3,1)= strs_xy dot
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optdmg
  real(8), intent(in) :: delt
  real(8), dimension(*), intent(in) :: ematpro
  real(8), intent(in) :: ehleng
  real(8), dimension(3,1) :: ipstrndot

  ! ------------------------------------
!  real(8), intent(out) :: effstrn
  real(8), intent(inout) :: damage
  real(8), dimension(3,1), intent(in) :: ipstrn
  real(8), dimension(3,1), intent(inout) :: ipstrs
  ! ------------------------------------

  ! ====================================
  ! local variable
  ! ==============
  real(8) :: young, poiss, thick

  real(8), dimension(2,2) :: etens2d

  real(8), dimension(3,3) :: elsvoit2d, csntvoit2d, ctanvoit2d
  real(8), dimension(2,2,2,2) :: ctantens2d

  real(8), dimension(3,1) :: ipstrsdot
  ! ====================================

  ! ------------------------------------------------
  ! get material properties
  young= ematpro(1)
  poiss= ematpro(2)
  thick= ematpro(20)
  ! ------------------------------------------------

  ! ----------------------------------------------------------------
  ! compute material damage
  ! -----------------------
  select case(optdmg)
  case(0)
     ! no damage
     damage= 0.0d0

     ! compute voight form elastic moduli matrix
     call getelsvoit2d(1,young,poiss, elsvoit2d)
       ! input : 1(optpty=plane-strs),young,poiss
       ! output : elsvoit2d

  case(1) ! lematire damage model
     etens2d(1,1)= ipstrn(1,1)
     etens2d(2,2)= ipstrn(2,1)
     etens2d(1,2)= 0.50d0*ipstrn(3,1)
     etens2d(2,1)= 0.50d0*ipstrn(3,1)

     ! compute tangent elastic moduli based on damage
     call getlemdmg2d0(1,ematpro,etens2d, damage, csntvoit2d,ctanvoit2d,ctantens2d)
        ! input : 1(optpty:p-strs),ematpro,etens2d
        ! inoutput : damage
        ! output : csntvoit2d,ctanvoit2d,ctantens2d

     ! compute voight form elastic moduli matrix
     elsvoit2d(1:3,1:3)= ctanvoit2d(1:3,1:3)

  case(2) ! linear softening with scaling
     etens2d(1,1)= ipstrn(1,1)
     etens2d(2,2)= ipstrn(2,1)
     etens2d(1,2)= 0.50d0*ipstrn(3,1)
     etens2d(2,1)= 0.50d0*ipstrn(3,1)

     ! compute tangent elastic moduli based on damage
     call getlineardmg2d0(1,ematpro,ehleng,etens2d, damage, csntvoit2d,ctanvoit2d,ctantens2d)
        ! input : 1(optpty:p-strs),ematpro,ehleng,etens2d
        ! inoutput : damage
        ! output : csntvoit2d,ctanvoit2d

     ! compute voight form elastic moduli matrix
     elsvoit2d(1:3,1:3)= ctanvoit2d(1:3,1:3)

  case default
     write(*,*) "not implemented damage model: gqfinthypo2d0"
     write(nout5,*) "not implemented damage model: gqfinthypo2d0"
     stop

  end select


  ! ----------------------------------------------------------------
  ! compute stress rate
  ! -------------------
  ! compute stress rate: plane stress law is used
  ! in-plane stress rate: [strs_x, strs_y. strs_xy]
  call dgemv('n',3,3,1.0d0,elsvoit2d,3,ipstrndot,1,0.0d0,ipstrsdot,1)
     ! input : 3,3,0, 3,1,0, 3,1, elsvoit2d,ipstrndot
     ! output : ipstrsdot

  ! update in-plane cauchy stress
  ipstrs(1,1)= ipstrs(1,1) + ipstrsdot(1,1) * delt ! sig_x
  ipstrs(2,1)= ipstrs(2,1) + ipstrsdot(2,1) * delt ! sig_y
  ipstrs(3,1)= ipstrs(3,1) + ipstrsdot(3,1) * delt ! sig_xy

  return
end subroutine getsighypobt1



subroutine getsighypobt2(ematpro,tsstrndot, tsstrsdot)
  !=======================================================================
  !  getsighypobt2 = compute rate of belytschko-tsay local element cauchy stress 
  !
  !                  note:
  !                  ----
  !                  belytschko, wong and chiang, CMAME, 1992, vol. 96, pp. 93-107
  !                  advances in one point quadrature shell elements
  !                  (see, apendix a.3 )
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ematpro(*) : material property
  !
  !  tsstrndot(2,1) : tranverse shear strain rate [2dxz,2dyz]
  !
  !  output:
  !  ------
  !  tsstrsdot(2,1) : transverse shear stress rate
  !                   tsstrsdot(1,1)= strs_xz dot
  !                   tsstrsdot(2,1)= strs_yz dot
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(*), intent(in) :: ematpro
  real(8), dimension(2,1), intent(in) :: tsstrndot

  real(8), dimension(2,1), intent(out) :: tsstrsdot
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: young, poiss, rk
  real(8) :: mu

  ! ====================================

  ! ------------------------------------------------
  ! get material properties
  young= ematpro(1)
  poiss= ematpro(2)
  rk= ematpro(19)     ! shear correction factor: 0.840d0

  ! ------------------------------------------------

  ! 2nd lame constant: shear modulus 
  mu= young / ( 2.0d0 * (1.0d0+poiss) )

  ! ----------------------------------------------------------------
  ! compute stress rate
  ! -------------------
  ! tranverse stress rate: [strs_xz, strs_yz]
  tsstrsdot(1,1)= mu * tsstrndot(1,1)
  tsstrsdot(2,1)= mu * tsstrndot(2,1)


  return
end subroutine getsighypobt2




subroutine getstrndotbt1(optcor,thick,zeta,eveloloc,bmat1pt,bcmat1pt, ipstrndot)
  !=======================================================================
  !  getstrndotbt1 = compute in plane local rate of deformation (velocity strain) of tb shell
  !                  with one point integration and local z method
  !                  for computational efficiency, explicit components form is used
  !
  !                 note:
  !                 ----
  !                 belytschko, wong and chiang, CMAME, 1992, vol. 96, pp. 93-107
  !                 advances in one point quadrature shell elements
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  thick : thickness of shell
  !
  !  zeta : [-1, +1] : pseudo-thickness parameter
  !
  !  ecurnloc(3,4) : local current nodal coordinate
  !
  !  eveloloc(5,4) : local nodal velocity: v_x, v_y, v_z, theta_x, theta_y
  !
  !  output:
  !  ------
  !  ltens2d(2,2) : velocity gradient tensor in 2d plane
  !  ipstrndot(3,1) : in plane strain rate [dx,dy,2dxy]
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, dimension(2), intent(in) :: optcor
  real(8), intent(in) :: thick, zeta
  real(8), dimension(5,4), intent(in) :: eveloloc
  real(8), dimension(2,4), intent(in) :: bmat1pt
  real(8), dimension(2,4), intent(in) :: bcmat1pt

  real(8), dimension(3,1), intent(out) :: ipstrndot
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: psibar
  real(8) :: dvxdx, dvxdy, dvydy, dvydx

  real(8) :: dx, dy, dxy

  ! loop index
  integer :: inode
  ! ====================================

  ! ----------------------------------------------------------------
  ! compute in plane local velocity strain: eq (18) and see, research note
  ! --------------------------------------
  ! pseudo-thickness : eq (5)
  psibar= ( zeta * thick ) / 2.0d0

  dvxdx= 0.0d0 ! initialize
  dvydy= 0.0d0
  dvxdy =0.0d0
  dvydx =0.0d0

  dx= 0.0d0
  dy= 0.0d0
  dxy= 0.0d0

  do inode=1, 4

    if(optcor(1) > 0) then
      ! --------------------------------------------------------------------------------------
      ! the original BLT element with warping correction
      ! ------------------------------------------------
      ! vx,x= b_xi v_xi + psibar*( b^c _xi v_xi + b_xi theta_yi )
      dvxdx= dvxdx + bmat1pt(1,inode)*eveloloc(1,inode) &
             + psibar*( bcmat1pt(1,inode)*eveloloc(1,inode) + bmat1pt(1,inode)*eveloloc(5,inode) )

      ! vy,y= b_yi v_yi + psibar*( b^c _yi v_yi - b_yi theta_xi )
      dvydy= dvydy + bmat1pt(2,inode)*eveloloc(2,inode) &
             + psibar*( bcmat1pt(2,inode)*eveloloc(2,inode) - bmat1pt(2,inode)*eveloloc(4,inode) )

      ! vx,y= b_yi v_xi + psibar*( b^c _yi v_xi + b_yi theta_yi)
      dvxdy= dvxdy + bmat1pt(2,inode)*eveloloc(1,inode) &
             + psibar*( bcmat1pt(2,inode)*eveloloc(1,inode) + bmat1pt(2,inode)*eveloloc(5,inode) )

      ! vy,x= b_xi v_yi + psibar*( b^c _xi v_yi - b_xi theta_xi )
      dvydx= dvydx + bmat1pt(1,inode)*eveloloc(2,inode)  &
             + psibar*( bcmat1pt(1,inode)*eveloloc(2,inode) - bmat1pt(1,inode)*eveloloc(4,inode) )

    else
      ! --------------------------------------------------------------------------------------
      ! the original BLT element
      ! -------------------------
      ! dx=
      dx= dx + bmat1pt(1,inode)*eveloloc(1,inode) + psibar*( bmat1pt(1,inode)*eveloloc(5,inode) )

      ! dy=
      dy= dy + bmat1pt(2,inode)*eveloloc(2,inode) - psibar*( bmat1pt(2,inode)*eveloloc(4,inode) )

      ! dxy=
      dxy= dxy + 0.50d0 * ( bmat1pt(2,inode)*eveloloc(1,inode) + bmat1pt(1,inode)*eveloloc(2,inode) &
           + psibar*( bmat1pt(2,inode)*eveloloc(5,inode) - bmat1pt(1,inode)*eveloloc(4,inode) ) )

    end if

  end do

  if(optcor(1) > 0) then
    ! --------------------------------------------------------------------------------------
    ! the original BLT element with warping correction
    ! ------------------------------------------------

    ! set in plane velocity strain rate
    ipstrndot(1,1)= dvxdx ! d_x
    ipstrndot(2,1)= dvydy ! d_y
    ipstrndot(3,1)= dvxdy + dvydx ! 2.0 d_xy

  else
    ! --------------------------------------------------------------------------------------
    ! the original BLT element
    ! ------------------------

    ! set in plane velocity strain rate
    ipstrndot(1,1)= dx
    ipstrndot(2,1)= dy
    ipstrndot(3,1)= 2.0d0 * dxy

  end if

  return
end subroutine getstrndotbt1




subroutine getstrndotbt2(optcor,eveloloc,bmat1pt,bsmat1pt, tsstrndot)
  !=======================================================================
  !  getstrndotbt2 = compute transverse shear local rate of deformation (velocity strain)
  !                  of tb shell with one point integration and local z method
  !                  for computational efficiency, explicit components form is used
  !
  !                  note:
  !                  ----
  !                  belytschko, wong and chiang, CMAME, 1992, vol. 96, pp. 93-107
  !                  advances in one point quadrature shell elements
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  eveloloc(5,4) : local nodal velocity: v_x, v_y, v_z, theta_x, theta_y
  !
  !  output:
  !  ------
  !  tsstrndot(2,1) : tranverse shear strain rate [2dxz,2dyz]
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, dimension(2), intent(in) :: optcor
  real(8), dimension(5,4), intent(in) :: eveloloc
  real(8), dimension(2,4), intent(in) :: bmat1pt
  real(8), dimension(2,3,4), intent(in) :: bsmat1pt

  real(8), dimension(2,1), intent(out) :: tsstrndot
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: dxz, dyz

  ! loop index
  integer :: inode
  ! ====================================

  ! ----------------------------------------------------------------
  ! compute transverse shear local velocity strain: eq (34) and see, research note
  ! ----------------------------------------------
  dxz= 0.0d0 ! initialize
  dyz= 0.0d0

  do inode=1, 4

    if(optcor(2) > 0) then
      ! --------------------------------------------------------------------------------------
      ! the original BLT element with shear correction
      ! ----------------------------------------------
      ! dxz= 0.50d0 * ( b^s _x1i v_zi + b^s _x2i theta_xi + b^s _x3i theta_yi )
      dxz= dxz + 0.50d0*( bsmat1pt(1,1,inode)*eveloloc(3,inode) &
                        + bsmat1pt(1,2,inode)*eveloloc(4,inode) &
                        + bsmat1pt(1,3,inode)*eveloloc(5,inode) )

      ! dyz= 0.50d0 * ( b^s _y1i v_zi + b^s _y2i theta_xi + b^s _y3i theta_yi )
      dyz= dyz + 0.50d0*( bsmat1pt(2,1,inode)*eveloloc(3,inode) &
                        + bsmat1pt(2,2,inode)*eveloloc(4,inode) &
                        + bsmat1pt(2,3,inode)*eveloloc(5,inode) )
 
    else

      ! --------------------------------------------------------------------------------------
      ! the original BLT element
      ! ------------------------
      ! dxz=
      dxz= dxz + 0.50d0*( bmat1pt(1,inode)*eveloloc(3,inode) + 0.250d0*eveloloc(5,inode) )

      ! dyz=
      dyz= dyz + 0.50d0*( bmat1pt(2,inode)*eveloloc(3,inode) - 0.250d0*eveloloc(4,inode) )

    end if

  end do

  ! ----------------------------------------------------------------
  ! set tranverse shear velocity strain
  ! -----------------------------------
  tsstrndot(1,1)= 2.0d0 * dxz
  tsstrndot(2,1)= 2.0d0 * dyz
  ! ----------------------------------------------------------------

  return
end subroutine getstrndotbt2

