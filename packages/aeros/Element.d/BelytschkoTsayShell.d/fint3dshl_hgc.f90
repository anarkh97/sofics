! ==================================
! hourglass control force: bt shell
! ==================================
!      type                  name                              argument
!      ----                  ----                              --------
! 2.  subroutine         updhgcstrsbt         (prmhgc,delt,ematpro,locbvec,ecurn,evelo, hgcvoitloc)
! 3.  subroutine         gethgcstrsdotbt      (prmhgc,ematpro,ecurnloc,eveloloc, hgcbstrsdot,hgcmstrsdot)
! 4.  subroutine         gethgcstrndotbt      (ecurnloc,eveloloc, hgcbstrndot,hgcmstrndot)
! 5.  subroutine         gethgconstbt         (prmhgc,ematpro,ecurnloc, c1,c2,c3)
! 6.  subroutine         gqfhgcbt             (locbvec,ecurn,evelo,hgcvoitloc, gqfhgc)
!
! =========================================================================================================



subroutine updhgcstrsbt(prmhgc,delt,ematpro,eveloloc,area,bmat1pt, &
                        gamma,zgamma, hgcvoitloc)
  !=======================================================================
  !  updhgcstrsbt = update hourglass control stress
  !
  !               note:
  !               ----
  !               chiang, m.s. thesis, northwestern univ., 1992
  !               advances in one point quadrature shell element
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !
  !  inoutput:
  !  --------
  !  hgcvoitloc(6,1) : voightform hourglass control stress
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(*), intent(in) :: prmhgc
  real(8), intent(in) :: delt
  real(8), dimension(*), intent(in) :: ematpro
  real(8), dimension(5,4), intent(in) :: eveloloc
  real(8), intent(in) :: area
  real(8), dimension(2,4), intent(in) :: bmat1pt
  real(8), dimension(4,1), intent(in) :: gamma
  real(8), intent(in) :: zgamma
  
  ! ------------------------------------

  real(8), dimension(6,1), intent(inout) :: hgcvoitloc
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(3,1) :: hgcbstrsdot
  real(8), dimension(2,1) :: hgcmstrsdot
  ! ====================================

  ! initialize: do not initialize

  ! -------------------------------------------------------------------------------
  ! compute hourglass stress rate
  ! -----------------------------
  call gethgcstrsdotbt(prmhgc,ematpro,eveloloc,area,bmat1pt, &
                       gamma,zgamma, hgcbstrsdot,hgcmstrsdot)
     ! input : prmhgc,ematpro,eveloloc,area,bmat1pt,gamma,zgamma
     ! output : hgcbstrsdot,hgcmstrsdot


  ! -------------------------------------------------------------------------------
  ! update hourglass control stress
  ! -------------------------------
  hgcvoitloc(1,1)= hgcvoitloc(1,1) + hgcbstrsdot(1,1) * delt ! sig_x
  hgcvoitloc(2,1)= hgcvoitloc(2,1) + hgcbstrsdot(2,1) * delt ! sig_y
  hgcvoitloc(3,1)= 0.0d0 ! sig_z (plane stress)

  hgcvoitloc(4,1)= hgcvoitloc(4,1) + hgcmstrsdot(2,1) * delt ! sig_yz
  hgcvoitloc(5,1)= hgcvoitloc(5,1) + hgcmstrsdot(1,1) * delt ! sig_xz
  hgcvoitloc(6,1)= hgcvoitloc(6,1) + hgcbstrsdot(3,1) * delt ! sig_xy
  ! -------------------------------------------------------------

  return
end subroutine updhgcstrsbt


subroutine gethgcstrsdotbt(prmhgc,ematpro,eveloloc,area,bmat1pt, &
                           gamma,zgamma, hgcbstrsdot,hgcmstrsdot)
  !=======================================================================
  !  gethgcstrsdotbt = compute hourglass stress rate
  !
  !                 note:
  !                 ----
  !                 chiang, m.s. thesis, northwestern univ., 1992
  !                 advances in one point quadrature shell element
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !
  !  output:
  !  ------
  !  hgcbstrsdot(3,1) : bending hourglass control stress rate
  !
  !  hgcmstrsdot(2,1) : membrane hourglass control stress rate
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(*), intent(in) :: prmhgc
  real(8), dimension(*), intent(in) :: ematpro
  real(8), dimension(5,4), intent(in) :: eveloloc
  real(8), intent(in) :: area
  real(8), dimension(2,4), intent(in) :: bmat1pt
  real(8), dimension(4,1), intent(in) :: gamma
  real(8), intent(in) :: zgamma

  real(8), dimension(3,1), intent(out) :: hgcbstrsdot
  real(8), dimension(2,1), intent(out) :: hgcmstrsdot
  ! ====================================
  ! local variable
  ! ==============
  real(8):: c1, c2, c3

  real(8), dimension(3,1) :: hgcbstrndot
  real(8), dimension(2,1) :: hgcmstrndot
  ! ====================================

  ! compute hourglass mode control constant
  call gethgconstbt(prmhgc,ematpro,area,bmat1pt, c1,c2,c3)
     ! input : prmhgc,ematpro
     ! output : c1,c2,c3


  ! compute hourglass strain rate
  call gethgcstrndotbt(eveloloc,gamma,zgamma, hgcbstrndot,hgcmstrndot)
     ! input : eveloloc,gamma,zgamma
     ! output : hgbstrndot,hgmstrndot 


  ! -------------------------------------------------------------
  ! compute hourglass stress rate
  ! bending
  hgcbstrsdot(1,1)= c1*hgcbstrndot(1,1)
  hgcbstrsdot(2,1)= c1*hgcbstrndot(2,1)
  hgcbstrsdot(3,1)= c2*hgcbstrndot(3,1)

  ! membrane
  hgcmstrsdot(1,1)= c3*hgcmstrndot(1,1)
  hgcmstrsdot(2,1)= c3*hgcmstrndot(2,1)
  ! -------------------------------------------------------------


  return
end subroutine gethgcstrsdotbt



subroutine gethgcstrndotbt(eveloloc,gamma,zgamma, hgcbstrndot,hgcmstrndot)
  !=======================================================================
  !  gethgcstrndotbt = compute hourglass strain rate for memebrane and bending
  !                  in membrane strain rate, warping correction is considered
  !
  !                  note:
  !                  ----
  !                  chiang, m.s. thesis, northwestern univ., 1992
  !                  advances in one point quadrature shell element
  !                  (see, eq (47), (48))
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !
  !  output:
  !  ------
  !  hgcbstrndot(3,1) : bending hourglass strain rate
  !
  !  hgcmstrndot(2,1) : membrane hourglass strain rate
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(5,4), intent(in) :: eveloloc
  real(8), dimension(4,1), intent(in) :: gamma
  real(8), intent(in) :: zgamma

  real(8), dimension(3,1), intent(out) :: hgcbstrndot
  real(8), dimension(2,1), intent(out) :: hgcmstrndot
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(3,1) :: gammatv
  real(8), dimension(2,1) :: gammatrot
  real(8), dimension(2,1) :: strot

  ! loop index
  integer :: inode, idime
  ! ====================================

  ! compute gamma^t v
  gammatv(:,:)= 0.0d0 ! initialize
  do inode=1, 4
     do idime=1, 3
        gammatv(idime,1)= gammatv(idime,1) + gamma(inode,1)*eveloloc(idime,inode)
     end do
  end do


  ! compute gamma^t theta
  gammatrot(:,:)= 0.0d0 ! initialize
  do inode=1, 4
     do idime=1, 2 ! only have x and y dof in rotation
        gammatrot(idime,1)= gammatrot(idime,1) + gamma(inode,1)*eveloloc(3+idime,inode)
     end do
  end do


  ! compute s ^t theta: s= [1 1 1 1] for rigid body motion
  strot(:,:)= 0.0d0 ! initialize
  do inode=1, 4
     do idime=1, 2 ! only have x and y dof in rotation
        strot(idime,1)= strot(idime,1) + eveloloc(3+idime,inode)
     end do
  end do


  ! -------------------------------------------------------------
  ! compute strain rate
  ! -------------------
  ! bending hourglass strain rate
  hgcbstrndot(1,1)= gammatrot(1,1) 
  hgcbstrndot(2,1)= gammatrot(2,1)
  hgcbstrndot(3,1)= gammatv(3,1)


  ! membrane hourglass strain rate
  hgcmstrndot(1,1)= gammatv(1,1) - (1.0d0/4.0d0)*zgamma*strot(2,1)
  hgcmstrndot(2,1)= gammatv(2,1) + (1.0d0/4.0d0)*zgamma*strot(1,1)
  ! -------------------------------------------------------------


  return
end subroutine gethgcstrndotbt





subroutine gethgconstbt(prmhgc,ematpro,area,bmat1pt, c1,c2,c3)
  !=======================================================================
  !  gethgconstbt = compute constant for hourglass mode control
  !                 warping correction is considered
  !
  !                 ref.:
  !                 ----
  !                 chiang, m.s. thesis, northwestern univ., 1992
  !                 advances in one point quadrature shell element
  !                 (see, eq (49))
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !
  !  output:
  !  ------
  !  c1, c2, c3 : hourglass mode control constant
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(*), intent(in) :: prmhgc
  real(8), dimension(*), intent(in) :: ematpro
  real(8), intent(in) :: area
  real(8), dimension(2,4), intent(in) :: bmat1pt

  real(8), intent(out) :: c1,c2,c3
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: rt, rw, rm, rk

  real(8) :: young, poiss, thick

  real(8) :: mu

  real(8) :: btb

  real(8) :: const1, const2, const3

  ! loop index
  integer :: inode, idime
  ! ====================================

  ! ------------------------------------------------
  ! hourglass mode control parameter
  rt= prmhgc(1) ! 0.010d0 
  rw= prmhgc(2) ! 0.010d0 
  rm= prmhgc(3) ! 0.010d0
  ! ------------------------------------------------
  ! get material properties
  young= ematpro(1)
  poiss= ematpro(2)
  thick= ematpro(20)

  ! shear correction factor
  rk= ematpro(19) ! 0.840d0 
  ! ------------------------------------------------

  ! compute shear modulus
  mu= young /( 2.0d0*(1.0d0+poiss) )

  ! compute b_xi b_xi + b_yi b_yi
  btb= 0.0d0 ! initialize
  do inode=1, 4
    do idime=1, 2
        btb= btb + bmat1pt(idime,inode) * bmat1pt(idime,inode)
    end do
  end do


  ! -------------------------------------------------------------
  ! compute c1
  const1= rt/192.0d0
  const2= young * thick**3 * area 
  const3= 1.0d0 + (2.0d0*rk*area)/(3.0d0*thick**2)
  c1= const1 * const2 * const3 * btb

  ! compute c2
  const1= rw/12.0d0
  const2= rk * mu * thick**3
  c2= const1 * const2 * btb

  ! compute c3
  const1= rm/8.0d0
  const2= young * thick * area
  c3= const1 * const2 * btb
  ! -------------------------------------------------------------


  return
end subroutine gethgconstbt


subroutine gqfhgcbt(gamma,zgamma,hgcvoitloc, gqfintloc)
  !=======================================================================
  !  elefhgclocbt = add local hourglass control force for bt shell
  !
  !                 note:
  !                 ----
  !                 chiang, m.s. thesis, northwestern univ., 1992
  !                 advances in one point quadrature shell element
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !
  !  in/output:
  !  ------
  !  efloc(24,1) : element hourglass control force in local coordinate
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(4,1), intent(in) :: gamma
  real(8), intent(in) :: zgamma

  ! -----------------------------------

  real(8), dimension(6,1), intent(inout) :: hgcvoitloc
  real(8), dimension(24,1), intent(inout) :: gqfintloc
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: const

  integer :: iloc

  ! loop index
  integer :: inode
  ! ====================================

  const = (1.0d0/4.0d0)*zgamma

  ! -------------------------------------------------------------------------------
  ! set local hourglass stabilization force
  ! ---------------------------------------
  do inode=1, 4

     ! address
     iloc= 6*(inode-1)

     ! add local hourglass control force
     gqfintloc(iloc+1,1)= gqfintloc(iloc+1,1) + gamma(inode,1)*hgcvoitloc(5,1)  ! f^hgc_x
     gqfintloc(iloc+2,1)= gqfintloc(iloc+2,1) + gamma(inode,1)*hgcvoitloc(4,1)  ! f^hgc_y
     gqfintloc(iloc+3,1)= gqfintloc(iloc+3,1) + gamma(inode,1)*hgcvoitloc(6,1)  ! f^hgc_z

     ! add local hourglass moment
     gqfintloc(iloc+4,1)= gqfintloc(iloc+4,1) + gamma(inode,1)*hgcvoitloc(1,1) + &
                          const*hgcvoitloc(4,1)  ! m^hgc_x
     gqfintloc(iloc+5,1)= gqfintloc(iloc+5,1) + gamma(inode,1)*hgcvoitloc(2,1) - &
                          const*zgamma*hgcvoitloc(5,1)  ! m^hgc_y
     !gqfintloc(iloc+6,1)= gqfintloc(iloc+6,1) + 0.0d0 ! m^hgc_z

  end do

  return
end subroutine gqfhgcbt
