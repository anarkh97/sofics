! ================================
! crack growth criterion: effective strain / mtps
! ================================
!      type                  name                              argument
!      ----                  ----                              --------
! 1.  subroutine           chkeffstrcri2d        (optstr,effpstr,svoit2d, crival,criang)
!
! =========================================================================================================



subroutine chkeffstrcri2d(optstr,effpstr,svoit2d, crival,criang)
  !=======================================================================
  !  chkeffstrcri2d= 2d effective stress/strain criterion in 2d
  !
  !                  note:
  !                  ----
  !                  von-mises criterion
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optstr : stress/strain option handler
  !           opt=0 : stress
  !           opt=1 : strain
  !
  !  effpstr : effective stress/strain
  !
  !  svoit2d(3,1) : 2d voight form stress or strain 
  !
  !  output:
  !  ------
  !  crival : criterion value
  !
  !  criang : crack kink angle form criterion
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optstr
  real(8), intent(in) :: effpstr
  real(8), dimension(3,1), intent(in) :: svoit2d

  real(8), intent(out) :: crival
  real(8), intent(out) :: criang
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: pval1, pval2
  real(8) :: pang1, pang2
  ! ====================================

  ! initialize
  crival= 0.0d0
  criang= 0.0d0


  ! compute principal tensile stress/strain
  call getptstr2d(optstr,svoit2d, pval1,pval2,pang1,pang2)
     ! input : optstr,svoit2d
     ! output : pval1,pval2,pang1,pang2

  ! set results
  crival= effpstr


! ### crack direction is computed from both of tension and compression
!  if ( abs(pval1) > abs(pval2) ) then

     ! minimum principal stress angle
!     criang= pang2

!  else if ( abs(pval2) > abs(pval1) ) then

     ! minimum principal stress angle
!     criang= pang1

!  end if

  ! perpendicular to maximum principal direction
  criang= pang2



  return
end subroutine chkeffstrcri2d
