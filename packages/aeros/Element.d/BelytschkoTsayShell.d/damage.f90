! ==================================
! continuum based damage model
! ==================================
!      type                  name                              argument
!      ----                  ----                              --------
! 1d model
! --------
! 1.  subroutine          getlemdmg1d             (ematpro,strn, damage)
! 2.  subroutine          getlineardmg1d          (ematpro,strn, damage)
!
! 2d model
! --------
! 3.  subroutine          getlemdmg2d0            (optpty,ematpro,strntens, damage, csntvoit,ctanvoit,ctantens)
! 5.  subroutine          getlineardmg2d0         (optpty,ematpro,ehleng,strntens, damage, csntvoit,ctanvoit,ctantens)
! 6.  subroutine          scalstrcurv1            (young,hele,ft,gf, nflag,eps0,eps1)  ! scalling linear softening curve
!
! =========================================================================================================



subroutine getlemdmg1d(ematpro,strn, damage)
  !=======================================================================
  !  getlemdmg2d = 1d lemaitre material damage model
  !                lemaitre, mechanics of solid materials, p. 449
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ematpro(*) : material properties
  !
  !  strn : 1d strain value
  !
  !  inoutput:
  !  --------
  !  damage : damage at last step
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(*), intent(in) :: ematpro
  real(8), intent(in) :: strn

  real(8), intent(inout) :: damage
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: cofa, cofb, epscr
  real(8) :: young
  real(8) :: damage0
  real(8) :: epsstar, f
  ! ====================================

  ! initialize: do not initialize damage


  ! store history value
  damage0= damage ! save the last damage

  ! ------------------------------------
  ! get material properties
  young= ematpro(1)

  ! set damagemodel parameter
  cofa= ematpro(4) ! softening parameter: cofa=1 : complete softening w/o residual stress
  cofb= ematpro(5) ! softening parameter: softening slope
  epscr= ematpro(6) ! critical strain
  ! ------------------------------------

  ! calculate the equivalent strain and corresponding increment
  epsstar= max(strn,0.0d0)**2 ! consider tension only
  epsstar= dsqrt(epsstar) ! the current equivalent strain

  ! always damage criterion should be less than 0.0d0: f <= 0.0d0
  f= max(epsstar - epscr, 0.0d0) ! damage criterion

  ! calculate damage parameter
  damage=1.0d0 - epscr*(1.0d0 - cofa)/max(epsstar, 0.10e-20) - cofa/dexp(cofb*f)

  ! damage can not recover
  if(damage < damage0) damage= damage0 ! same as previous damage

  ! maximum value of damage parameter is greater than 0.0d0 and less than 1.0d0
  damage= max(0.0d0, damage)
  damage= min(1.0d0, damage)



  return
end subroutine getlemdmg1d





subroutine getlineardmg1d(ematpro,strn, damage)
  !=======================================================================
  !  getlineardmg1d = 1d isotropic linear softening damage model
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ematpro(*) : material properties
  !
  !  strn : 1d strain value
  !
  !  inoutput:
  !  --------
  !  damage : damage at last step
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(*), intent(in) :: ematpro
  real(8), intent(in) :: strn

  real(8), intent(inout) :: damage
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: ft, gf, eps0, eps1
  real(8) :: young
  real(8) :: const
  real(8) :: damage0
  real(8) :: epsstar
  ! ====================================

  ! initialize: do not initialize damage


!     strs
!       ^
!    ft |       *
!       |      *|  *
!       |     * |     *
!       |    *  |        *
!       |   *   |           *
!       |  *    |              *
!       | *     |                 *
!       |*      |                    *
!       *-------|-----------------------*--------> strn
!              eps0                    eps1
!
!       fracture energy: gf
!       gf = (1/2) * ft * eps1= (1/2) * young * eps0 * eps1



  ! ------------------------------------
  ! get material properties
  young= ematpro(1)

  ! set damagemodel parameter
  ft= ematpro(4) ! maximum tensile strength
  gf= ematpro(5) ! fracture energy

  ! ------------------------------------

  ! compute eps0
  eps0= ft / young

  ! compute eps1: eps1= gf *  ( 2 / (young * eps0 ) )
  eps1= 2.0d0 * gf / ( young * eps0 )

  ! check negetive tangential softening
  if ( eps1 < eps0 ) then
    write(*,*) "assigned fracture enenrgy is too small: getlineardmg1d"
    write(nout5,*) "assigned fracture enenrgy is too small: getlineardmg1d"
    stop
  end if
  ! ------------------------------------

  ! store history value
  damage0= damage ! save the last damage

  ! calculate the equivalent strain and corresponding increment
  epsstar= max( strn,0.0d0 )**2 ! tension only
  epsstar= dsqrt(epsstar) ! the current equivalent strain

  ! compute damage parameter for liner softening
   if ( epsstar <= eps0 ) then
      damage= 0.0d0

  else if ( eps0 < epsstar .and. epsstar < eps1 ) then
     const= eps0 / ( eps1 - eps0 )
     damage= 1.0d0 - const * ( eps1 / epsstar - 1.0d0 )
  
  else if ( eps1 <= epsstar ) then
      damage= 1.0d0

  end if

  ! damage can not recover
  if(damage < damage0) damage= damage0 ! same as previous damage

  ! maximum value of damage parameter is greater than 0.0d0 and less than 1.0d0
  damage= max(0.0d0, damage)
  damage= min(1.0d0, damage)



  return
end subroutine getlineardmg1d





subroutine getlemdmg2d0(optpty,ematpro,strntens, damage, csntvoit,ctanvoit,ctantens)
  !=======================================================================
  !  getlemdmg2d0 = 2d lemaitre material damage model
  !                lemaitre, mechanics of solid materials, p. 449
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optptyp : plane stress/strain option
  !
  !  ematpro(*) : material properties
  !
  !  strntens(2,2) : 2d strain tensor
  !
  !  inoutput:
  !  --------
  !  damage : damage at last step
  !
  !  output:
  !  ------
  !  csntvoit(3,3) : voighr form secant modulus
  !
  !  ctanvoit(3,3) : voight form tangent modulus
  !
  !  ctantens(2,2,2) : tangent modulus tensor
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optpty
  real(8), dimension(*), intent(in) :: ematpro
  real(8), dimension(2,2), intent(in) :: strntens

  real(8), intent(inout) :: damage

  real(8), dimension(3,3), intent(out) :: csntvoit, ctanvoit
  real(8), dimension(2,2,2,2), intent(out) :: ctantens
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: cofa, cofb, epscr
  real(8) :: young, poiss
  real(8), dimension(3,3) :: c0
  real(8) :: damage0

  real(8) :: pval1, pval2
  real(8) :: epsstar, f

  real(8) :: c1111, c1122, c2211, c2222, c1212
  real(8) :: e11, e22, e12
  real(8) :: delta
  real(8) :: ee11, ee22
  real(8) :: heavi
  real(8) :: pestarpee11, pestarpee22
  real(8) :: pee11pe11, pee11pe22, pee22pe11, pee22pe22, pee11pe12, pee22pe12
  real(8) :: pestarpe11, pestarpe22, pestarpe12
  real(8) :: pdpestar
  real(8) :: pdpe11, pdpe22, pdpe12
  real(8) :: s11, s22, s12

  ! loop index
  integer :: idime
  ! ====================================

  ! initialize: do not initialize damage
  csntvoit(1:3,1:3)= 0.0d0
  ctanvoit(1:3,1:3)= 0.0d0
  ctantens(1:2,1:2,1:2,1:2)= 0.0d0


  ! store history value
  damage0= damage ! save the last damage

  ! ------------------------------------
  ! get material properties
  young= ematpro(1)
  poiss= ematpro(2)

  ! set damagemodel parameter
  cofa= ematpro(4) ! softening parameter: cofa=1 : complete softening w/o residual stress
  cofb= ematpro(5) ! softening parameter: softening slope
  epscr= ematpro(6) ! strain thresholds
  ! ------------------------------------

  ! get vergin elastic moduli
  call getelsvoit2d(optpty,young,poiss, c0)
     ! input : optpty,young,poiss
     ! output : elsvoit2d

  ! compute principal tensile strain
  e11= strntens(1,1)
  e22= strntens(2,2)
  e12= strntens(1,2)
  pval1= (e11+e22)/2.0d0 + dsqrt( ((e11-e22)/2.0d0)**2 + e12**2 )
  pval2= (e11+e22)/2.0d0 - dsqrt( ((e11-e22)/2.0d0)**2 + e12**2 )

  ! calculate the equivalent strain and corresponding increment
  epsstar= max( pval1,0.0d0 )**2 + max( pval2,0.0d0 )**2 ! consider tension only
  epsstar= dsqrt( epsstar ) ! the current equivalent strain

  ! always damage criterion should be less than 0.0d0: f <= 0.0d0
  f= max(epsstar - epscr, 0.0d0) ! damage criterion

  ! calculate damage parameter
  damage=1.0d0 - epscr*(1.0d0 - cofa)/max(epsstar, 0.10e-20) - cofa/dexp(cofb*f)

  ! damage can not recover
  if(damage < damage0) damage= damage0 ! same as previous damage

  ! maximum value of damage parameter is greater than 0.0d0 and less than 1.0d0
  damage= max(0.0d0, damage)
  damage= min(1.0d0, damage)

  ! -------------------------------------------------------
  ! set result: voight form secant modulus
  csntvoit(1:3,1:3)= (1.0d0 - damage) * c0(1:3,1:3)
  ! -------------------------------------------------------

  ! calculate the tangent modulus
  ! set elastic moduli components
  c1111=c0(1,1)
  c1122=c0(1,2)
  c2211=c0(2,1)
  c2222=c0(2,2)
  c1212=c0(3,3)

  ! set strain components
  e11= strntens(1,1)
  e22= strntens(2,2)
  e12= strntens(1,2)

  delta= dsqrt( (e11+e22)**2 - 4.0d0 * ( e11*e22 - e12**2 ) )

  if ( delta == 0.0d0 ) delta= 1.0e-20

  ee11= 0.50d0*((e11+e22)+delta)
  ee22= 0.50d0*((e11+e22)-delta)

  if ( epsstar /= 0.0d0 ) then
     pestarpee11= ee11/epsstar * heavi(ee11)
     pestarpee22= ee22/epsstar * heavi(ee22)

  else
     pestarpee11= 0.0d0
     pestarpee22= 0.0d0

  end if
   
  pee11pe11= 0.50d0 * ( 1.0d0 + (e11-e22)/delta )
  pee11pe22= 0.50d0 * ( 1.0d0 - (e11-e22)/delta )
  pee22pe11= 0.50d0 * ( 1.0d0 - (e11-e22)/delta )
  pee22pe22= 0.50d0 * ( 1.0d0 + (e11-e22)/delta )
  pee11pe12= 0.50d0 * ( 4.0d0 * e12 / delta )
  pee22pe12= 0.50d0 * ( -4.0d0 * e12 / delta )

  pestarpe11= pestarpee11 * pee11pe11 + pestarpee22 * pee22pe11
  pestarpe22= pestarpee11 * pee11pe22 + pestarpee22 * pee22pe22
  pestarpe12= pestarpee11 * pee11pe12 + pestarpee22 * pee22pe12

  ! a damage / a eps^star
  if ( damage > damage0 ) then
     pdpestar= epscr * ( 1.0d0 - cofa ) / epsstar**2 + cofa*cofb / dexp( cofb * f )
  else
     pdpestar= 0.0d0
  end if

  pdpe11= pdpestar * pestarpe11
  pdpe22= pdpestar * pestarpe22
  pdpe12= pdpestar * pestarpe12

  ! compute fictious stress
  s11= c1111 * e11 + c1122 * e22
  s22= c2211 * e11 + c2222 * e22
  s12= 2.0d0 * c1212 * e12

  ! compute components of tangent modulus
  ctantens(1,1,1,1)= (1.0d0-damage)*c1111 - pdpe11*s11
  ctantens(1,1,2,2)= (1.0d0-damage)*c1122 - pdpe22*s11
  ctantens(2,2,1,1)= (1.0d0-damage)*c2211 - pdpe11*s22
  ctantens(2,2,2,2)= (1.0d0-damage)*c2222 - pdpe22*s22

  ctantens(1,1,1,2)= -pdpe12 * s11
  ctantens(2,2,1,2)= -pdpe12 * s22
  ctantens(1,2,1,1)= -pdpe11 * s12
  ctantens(1,2,2,2)= -pDpe22 * s12

  ctantens(1,1,2,1)= ctantens(1,1,1,2)
  ctantens(2,2,2,1)= ctantens(2,2,1,2)
  ctantens(2,1,1,1)= ctantens(1,2,1,1)
  ctantens(2,1,2,2)= ctantens(1,2,2,2)

  ctantens(1,2,1,2)= (1.0d0-damage)*c1212 - pdpe12*s12
  ctantens(1,2,2,1)= ctantens(1,2,1,2)
  ctantens(2,1,1,2)= ctantens(1,2,1,2)
  ctantens(2,1,2,1)= ctantens(1,2,1,2)

  ! -------------------------------------------------------
  ! set result: voight form tangent modulus

  ctanvoit(1,1)= ctantens(1,1,1,1)
  ctanvoit(1,2)= ctantens(1,1,2,2)
  ctanvoit(1,3)= ctantens(1,1,1,2)

  ctanvoit(2,1)= ctantens(2,2,1,1)
  ctanvoit(2,2)= ctantens(2,2,2,2)
  ctanvoit(2,3)= ctantens(2,2,1,2)

  ctanvoit(3,1)= ctantens(1,2,1,1)
  ctanvoit(3,2)= ctantens(1,2,2,2)
  ctanvoit(3,3)= ctantens(1,2,1,2)

  ! -------------------------------------------------------



  return
end subroutine getlemdmg2d0





subroutine getlineardmg2d0(optpty,ematpro,ehleng,strntens, damage, csntvoit,ctanvoit,ctantens)
  !=======================================================================
  !  getlineardmg2d0 = 2d isotropic linear softening model
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optptyp : plane stress/strain option
  !
  !  ematpro(*) : material properties
  !
  !  ehleng : characteristic element size
  !
  !  strntens(2,2) : 2d strain tensor
  !
  !  inoutput:
  !  --------
  !  damage : damage at last step
  !
  !  output:
  !  ------
  !  csntvoit(3,3) : voighr form secant modulus
  !
  !  ctanvoit(3,3) : voight form tangent modulus
  !
  !  ctantens(2,2,2) : tangent modulus tensor
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optpty
  real(8), dimension(*), intent(in) :: ematpro
  real(8), intent(in) :: ehleng
  real(8), dimension(2,2), intent(in) :: strntens

  real(8), intent(inout) :: damage

  real(8), dimension(3,3), intent(out) :: csntvoit, ctanvoit
  real(8), dimension(2,2,2,2), intent(out) :: ctantens
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: ft, gf, scalftcr
  real(8) :: hele
  real(8) :: eps0, eps1
  integer :: nflag
  real(8) :: young, poiss
  real(8), dimension(3,3) :: c0
  real(8) :: damage0

  real(8) :: pval1, pval2
  real(8) :: epsstar
  real(8) :: const

  real(8) :: c1111, c1122, c2211, c2222, c1212
  real(8) :: e11, e22, e12
  real(8) :: delta
  real(8) :: ee11, ee22
  real(8) :: heavi
  real(8) :: pestarpee11, pestarpee22
  real(8) :: pee11pe11, pee11pe22, pee22pe11, pee22pe22, pee11pe12, pee22pe12
  real(8) :: pestarpe11, pestarpe22, pestarpe12
  real(8) :: pdpestar
  real(8) :: pdpe11, pdpe22, pdpe12
  real(8) :: s11, s22, s12

  ! loop index
  integer :: idime
  ! ====================================

  ! initialize: do not initialize damage
  csntvoit(1:3,1:3)= 0.0d0
  ctanvoit(1:3,1:3)= 0.0d0
  ctantens(1:2,1:2,1:2,1:2)= 0.0d0


!     strs
!       ^
!    ft |       *
!       |      *|  *
!       |     * |     *
!       |    *  |        *
!       |   *   |           *
!       |  *    |              *
!       | *     |                 *
!       |*      |                    *
!       *-------|-----------------------*--------> strn
!              eps0                    eps1
!
!       fracture energy: gf
!       gf / hele = (1/2) * ft * eps1= (1/2) * young * eps0 * eps1



  ! ------------------------------------
  ! get material properties
  young= ematpro(1)
  poiss= ematpro(2)

  ! set damagemodel parameter
  ft= ematpro(4) ! maximum tensile strength
  gf= ematpro(5) ! fracture energy
  scalftcr= ematpro(6) ! element size scaling factor
                       ! note: scalftcr= 1.0d0 means no scalling                      
  ! ------------------------------------

  ! store history value
  damage0= damage ! save the last damage

  ! set element size for fracture energy scalling
  if ( scalftcr == 1.0d0 ) then
     hele= 1.0d0

  else
     hele= ehleng

  end if
 
  ! scaling softening slop according to element size
  call scalstrcurv1(young,hele,ft,gf, nflag,eps0,eps1)
     ! input : young,hele,ft,gf
     ! output : nflag,eps0,eps1

  ! get vergin elastic moduli
  call getelsvoit2d(optpty,young,poiss, c0)
     ! input : optpty,young,poiss
     ! output : elsvoit2d

  ! compute principal tensile strain
  e11= strntens(1,1)
  e22= strntens(2,2)
  e12= strntens(1,2)
  pval1= (e11+e22)/2.0d0 + dsqrt( ((e11-e22)/2.0d0)**2 + e12**2 )
  pval2= (e11+e22)/2.0d0 - dsqrt( ((e11-e22)/2.0d0)**2 + e12**2 )

  ! calculate the equivalent strain and corresponding increment
  epsstar= max( pval1,0.0d0 )**2 + max( pval2,0.0d0 )**2 ! tension only
  epsstar= dsqrt( epsstar ) ! the current equivalent strain

  ! compute damage parameter for liner softening
   if ( epsstar <= eps0 ) then
      damage= 0.0d0

  else if ( eps0 < epsstar .and. epsstar < eps1 ) then
     const= eps0 / ( eps1 - eps0 )
     damage= 1.0d0 - const * ( eps1 / epsstar - 1.0d0 )
  
  else if ( eps1 <= epsstar ) then
      damage= 1.0d0

  end if

  ! damage can not recover
  if(damage < damage0) damage= damage0 ! same as previous damage

  ! maximum value of damage parameter is greater than 0.0d0 and less than 1.0d0
  damage= max(0.0d0, damage)
  damage= min(1.0d0, damage)

  ! -------------------------------------------------------
  ! set result: voight form secant modulus
  csntvoit(1:3,1:3)= (1.0d0 - damage) * c0(1:3,1:3)
  ! -------------------------------------------------------

  ! calculate the tangent modulus
  ! set elastic moduli components
  c1111=c0(1,1)
  c1122=c0(1,2)
  c2211=c0(2,1)
  c2222=c0(2,2)
  c1212=c0(3,3)

  ! set strain components
  e11= strntens(1,1)
  e22= strntens(2,2)
  e12= strntens(1,2)

  delta= dsqrt( (e11+e22)**2 - 4.0d0 * ( e11*e22 - e12**2 ) )

  if ( delta == 0.0d0 ) delta= 1.0e-20

  ee11= 0.50d0*((e11+e22)+delta)
  ee22= 0.50d0*((e11+e22)-delta)

  if ( epsstar /= 0.0d0 ) then
     pestarpee11= ee11/epsstar * heavi(ee11)
     pestarpee22= ee22/epsstar * heavi(ee22)

  else
     pestarpee11= 0.0d0
     pestarpee22= 0.0d0

  end if
   
  pee11pe11= 0.50d0 * ( 1.0d0 + (e11-e22)/delta )
  pee11pe22= 0.50d0 * ( 1.0d0 - (e11-e22)/delta )
  pee22pe11= 0.50d0 * ( 1.0d0 - (e11-e22)/delta )
  pee22pe22= 0.50d0 * ( 1.0d0 + (e11-e22)/delta )
  pee11pe12= 0.50d0 * ( 4.0d0 * e12 / delta )
  pee22pe12= 0.50d0 * ( -4.0d0 * e12 / delta )


  pestarpe11= pestarpee11 * pee11pe11 + pestarpee22 * pee22pe11
  pestarpe22= pestarpee11 * pee11pe22 + pestarpee22 * pee22pe22
  pestarpe12= pestarpee11 * pee11pe12 + pestarpee22 * pee22pe12

  ! a damage / a eps^star
  if ( damage > damage0 ) then
     const= eps0 / ( eps1 - eps0 )
     pdpestar= const * eps1 / epsstar**2

  else
     pdpestar= 0.0d0
  end if

  pdpe11= pdpestar * pestarpe11
  pdpe22= pdpestar * pestarpe22
  pdpe12= pdpestar * pestarpe12

  ! compute fictious stress
  s11= c1111 * e11 + c1122 * e22
  s22= c2211 * e11 + c2222 * e22
  s12= 2.0d0 * c1212 * e12

  ! compute components of tangent modulus
  ctantens(1,1,1,1)= (1.0d0-damage)*c1111 - pdpe11*s11
  ctantens(1,1,2,2)= (1.0d0-damage)*c1122 - pdpe22*s11
  ctantens(2,2,1,1)= (1.0d0-damage)*c2211 - pdpe11*s22
  ctantens(2,2,2,2)= (1.0d0-damage)*c2222 - pdpe22*s22

  ctantens(1,1,1,2)= -pdpe12 * s11
  ctantens(2,2,1,2)= -pdpe12 * s22
  ctantens(1,2,1,1)= -pdpe11 * s12
  ctantens(1,2,2,2)= -pdpe22 * s12

  ctantens(1,1,2,1)= ctantens(1,1,1,2)
  ctantens(2,2,2,1)= ctantens(2,2,1,2)
  ctantens(2,1,1,1)= ctantens(1,2,1,1)
  ctantens(2,1,2,2)= ctantens(1,2,2,2)

  ctantens(1,2,1,2)= ( 1.0d0 - damage ) * c1212 - pdpe12*s12
  ctantens(1,2,2,1)= ctantens(1,2,1,2)
  ctantens(2,1,1,2)= ctantens(1,2,1,2)
  ctantens(2,1,2,1)= ctantens(1,2,1,2)

  ! -------------------------------------------------------
  ! set result: voight form tangent modulus
  ctanvoit(1,1)= ctantens(1,1,1,1)
  ctanvoit(1,2)= ctantens(1,1,2,2)
  ctanvoit(1,3)= ctantens(1,1,1,2)

  ctanvoit(2,1)= ctantens(2,2,1,1)
  ctanvoit(2,2)= ctantens(2,2,2,2)
  ctanvoit(2,3)= ctantens(2,2,1,2)

  ctanvoit(3,1)= ctantens(1,2,1,1)
  ctanvoit(3,2)= ctantens(1,2,2,2)
  ctanvoit(3,3)= ctantens(1,2,1,2)



  return
end subroutine getlineardmg2d0





subroutine scalstrcurv1(young,ehleng,ft,gf, nflag,eps0,eps1)
  !=======================================================================
  !  scalstrcurv1 = scaling stress-strain curve which has linear softening
  !
  !                 note:
  !                 ----
  !                 we need to scale stress-strain curve
  !                 to remove mesh dependency in softening material.
  !
  !                 refer fracture ans size effect by prof. bazant, page 251
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  young : young's modulus
  !
  !  ehleng : characteristic element size
  !
  !  ft : maximum tensile strength
  !
  !  gf : fracture energy
  !
  !  output:
  !  ------
  !  nflag : result flag
  !          nflag= 0 : scalling is sucess
  !          nflag= 1 : negative tangential softening is happened
  !
  !  eps0 : strain at peak
  !
  !  eps1 : strain at vanish
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), intent(in) :: young, ehleng
  real(8), intent(in) :: ft, gf
  
  integer, intent(out) :: nflag
  real(8), intent(out) :: eps0, eps1
  ! ====================================
  ! local variable
  ! ==============

  ! loop index
  ! ====================================

  ! initilaize
  nflag= 0
  eps0= 0.0d0
  eps1= 0.0d0


!     strs
!       ^
!    ft |       *
!       |      *|  *
!       |     * |     *
!       |    *  |        *
!       |   *   |           *
!       |  *    |              *
!       | *     |                 *
!       |*      |                    *
!       *-------|-----------------------*--------> strn
!              eps0                    eps1
!
!       fracture energy: gf
!       gf / hele = (1/2) * ft * eps1= (1/2) * young * eps0 * eps1


  ! compute eps0
  eps0= ft / young

  ! compute eps1: eps1= ( gf / hele ) *  ( 2 / (young * eps0 ) )
  eps1= 2.0d0 * gf / ( ehleng * young * eps0 )
!  eps1= 2.0d0 * gf / ( ehleng * young * eps0 ) + eps0 ! exclude elastic regim from fracture energy

  ! check negetive tangential softening
  if ( eps1 < eps0 ) then
    nflag= 1
    write(*,*) "element size is too big: scalstrcurv1"
    write(nout5,*) "element size is too big: scalstrcurv1"
    stop
  end if



  return
end subroutine scalstrcurv1
