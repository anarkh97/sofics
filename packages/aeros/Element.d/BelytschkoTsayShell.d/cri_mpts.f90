! ================================
! crack growth criterion
! ================================
!      type                  name                              argument
!      ----                  ----                              --------
! 1.  subroutine           maxcri2d1             (prmcri,mgqpt1,crival,criang, nflagcr,crvalue,crangle)
! 2.  subroutine           chkmptcri2d           (optstr,svoit2d, crival,criang)
! 3.  subroutine           chkmptscri2d          (optstr,svoit2d, crival,criang)
! 4.  subroutine           getptstr2d            (optstr,svoit2d, pval1,pval2,pang1,pang2)
! 5.  subroutine           getpsstr2d            (optstr,svoit2d, pval1,pang1)
!
! =========================================================================================================



subroutine getptstr2d(optstr,svoit2d, pval1,pval2,pang1,pang2)
  !=======================================================================
  !  getptstr2d = compute principal tensile stress/strain
  !
  !               note:
  !               ----
  !               minimum principal direction angle: -90 < pang(2) < 90
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optstr : stress/strain option handler
  !           opt=0 : stress
  !           opt=1 : strain
  !
  !  svoit2d(3,1) : 2d voight form stress or strain 
  ! 
  !  output:
  !  ------
  !  pval1= maximum principle tensile stress
  !
  !  pval2= minimum principle tensile stress
  !
  !  pang1= maximum principle tensile stress angle
  !
  !  pang2= minimum principle tensile stress angle
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optstr
  real(8), dimension(3,1), intent(in) :: svoit2d

  real(8), intent(out) :: pval1, pval2
  real(8), intent(out) :: pang1, pang2
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: s11, s12, s22
  ! ====================================

  ! initialize
  pval1= 0.0d0
  pval2= 0.0d0
  pang1= 0.0d0
  pang2= 0.0d0


  ! set component
  s11= svoit2d(1,1)
  s22= svoit2d(2,1)

  select case(optstr)
  case(0) ! stress
     s12= svoit2d(3,1)

  case(1) ! strain
     s12= 0.50d0 * svoit2d(3,1)

  case default
     write(*,*) "wrong opt: getptstr2d"
     write(nout5,*) "wrong opt: getptstr2d"
     stop

  end select

  ! compute principal stress/strain
  pval1= ( s11 + s22 ) / 2.0d0 + dsqrt( ( ( s11 - s22 ) / 2.0d0 )**2 + s12**2 ) ! maximum
  pval2= ( s11 + s22 ) / 2.0d0 - dsqrt( ( ( s11 - s22 ) / 2.0d0 )**2 + s12**2 ) ! minimum

  ! -----------------------------------
  ! note: description of generic ATAN2
  ! ----
  !
  ! The result type is the same as x. The value lies in the range -pi < ATAN2 (y, x) <= pi.
  ! If x /= zero, the result is approximately equal to the value of arctan (y/x). 
  ! -----------------------------------
  ! If y > zero, the result is positive. 
  ! If y < zero, the result is negative. 
  ! If y = zero, the result is zero (if x > zero) or pi (if x < zero). 
  ! If x = zero, the absolute value of the result is pi/2
  ! -----------------------------------

  ! maximum principal direction
  pang1= 0.50d0 * datan2( 2.0d0*s12,s11-s22 )

  ! minimum principal direction: -90 < pang2 <= 90
  ! pang1= 0 or (+/-)180
  if ( pang1 == 0.0d0 .or. pang1 >= 180.0d0*d2r .or. pang1 <= -180.0d0*d2r ) then
      pang2= 90.0d0 * d2r

  ! pang1=(+/-)90
  else if ( pang1 == 90.0d0*d2r .or. pang1 == -90.0d0*d2r ) then
      pang2= 0.0d0

  ! 0 < pang1 < 90
  else if ( 0.0d0 < pang1 .and. pang1 < 90.0d0*d2r ) then
      pang2= pang1 - 90.0d0 * d2r

  ! 90 < pang1 < 180
  else if ( 90.0d0*d2r < pang1 .and. pang1 < 180.0d0*d2r ) then
      pang2= pang1 - 90.0d0 * d2r

  ! -90 < pang1 < 0
  else if ( -90.0d0*d2r < pang1 .and. pang1< 0.0d0 ) then
      pang2= pang1 + 90.0d0 * d2r

  ! -180 < pang1 < -90
  else if ( -180.0d0*d2r < pang1 .and. pang1< -90.0d0*d2r ) then
      pang2= pang1 + 90.0d0 * d2r

  end if


  
  return
end subroutine getptstr2d


