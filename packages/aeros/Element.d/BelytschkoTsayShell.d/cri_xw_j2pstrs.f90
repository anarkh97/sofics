! ================================
! crack growth criterion : bai-wierzbicki criterion for von-mises material
! ================================
!      type                  name                              argument
!      ----                  ----                              --------
! 1.  subroutine           xwj2pstrscri2d        (prmcri,mgqpt1,crival,criang,etatri, nflagcr,crvalue,crangle)
! 2.  subroutine           chkxwj2pstrscri2d     (optstr,effpstr,svoit2d, crival,criang,etatri)
! 3.  subroutine           xwpstrnfj2pstrs       (c1,c2,a,n,etatri, pstrnf)
!
! =========================================================================================================



subroutine xwj2pstrscri2d(prmcri,mgqpt1,crival,criang,etatri, nflagcr,crvalue,crangle)
  !=======================================================================
  !  maxcri2d1= 
  !
  !            note:
  !            ----
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  prmcri(*) : criterion parameter
  !
  !  mgqpt1 : number of gq point
  !
  !  crival(mgqpt1) : pre computed crack criterion value
  !
  !  criang(mgqpt1) : pre computed crack propagation angle
  !
  !  output:
  !  ------
  !  nflagcr : crack initiation flag
  !            nflagcr=0 : crack criterion is not satisfied
  !            nflagcr=1 : crack criterion is satisfied
  !
  !  crvalue : crack propagation value
  !
  !  crangle : crack propagation angle
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(*), intent(in) :: prmcri
  integer, intent(in) :: mgqpt1
  real(8), dimension(mgqpt1), intent(in) :: crival, criang, etatri

  integer, intent(out) :: nflagcr
  real(8), intent(out) :: crvalue
  real(8), intent(out) :: crangle
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: c1, c2, a, n
  real(8) :: pstrnf

  ! loop index
  integer :: igaus
  ! ====================================

  ! initialize
  nflagcr= 0
  crvalue= 0.0d0
  crangle= 0.0d0


  ! -------------------------------------------
  ! get criterion parameters
  c1= prmcri(5) ! critical value
  c2= prmcri(5) ! critical value
  a= prmcri(5) ! critical value
  n= prmcri(5) ! critical value
  ! -------------------------------------------

  ! compute

  ! loop over gauss point
  do igaus=1, mgqpt1

     ! compute fracture plastic strain
     call xwpstrnfj2pstrs(c1,c2,a,n,etatri(igaus), pstrnf)
        ! inpute : c1,c2,a,n,etatri(igaus)
        ! output : pstrnf

     if ( crival(igaus) >= pstrnf ) then

        ! increase counter
        nflagcr= nflagcr + 1

        ! crack propagation
        crvalue= crvalue + crival(igaus)
        
        ! crack propagation angle
        crangle= crangle + criang(igaus)

     end if

  end do ! do igaus

  ! ------------------------------------------------------------------
  ! set flag and average angle
  if ( nflagcr /= 0 ) then

     crvalue= crvalue / real(nflagcr)

     crangle= crangle / real(nflagcr)

     nflagcr= 1

  end if
  ! ------------------------------------------------------------------



  return
end subroutine xwj2pstrscri2d





subroutine chkxwj2pstrscri2d(optstr,effpstr,svoit2d, crival,criang,etatri)
  !=======================================================================
  !  chkxwj2pstrscri2d= 
  !
  !               note:
  !               ----
  !               2d plane stress case
  !               von-mises material
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
  !  etatri : stress triaxiality parameter for plane stress case
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
  real(8), intent(out) :: etatri
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: pval1, pval2
  real(8) :: pang1, pang2
  
  real(8) :: alpha
  real(8) :: term1, term2
  ! ====================================

  ! initialize
  crival= 0.0d0
  criang= 0.0d0
  etatri= 0.0d0


  ! compute principal tensile stress/strain
  call getptstr2d(optstr,svoit2d, pval1,pval2,pang1,pang2)
     ! input : optstr,svoit2d
     ! output : pval1,pval2,pang1,pang2

  ! set criterion value
  crival= effpstr

  ! perpendicular to maximum principal direction
  criang= pang2

  ! compute triaxial stress parameter for plane stress case
  if ( pval1 /= 0.0d0 ) then
     alpha= pval2 / pval1

     term1= 1.0d0 + alpha
     term2= dsqrt( 3.0d0 * (1.0d0 + alpha + alpha**2 ) )  
     etatri= term1 / term2

  end if



  return
end subroutine chkxwj2pstrscri2d





subroutine xwpstrnfj2pstrs(c1,c2,a,n,etatri, pstrnf)
  !=======================================================================
  !  xwpstrnfj2pstrs= compute effective plastic strain to fracture by xue-wierzbicki criterion
  !
  !         note:
  !         ----
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !
  !  output:
  !  ------
  !  pstrnf : effective plastic fracture strain
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), intent(in) :: c1, c2
  real(8), intent(in) :: a
  real(8), intent(in) :: n
  real(8), intent(in) :: etatri
  
  real(8), intent(out) :: pstrnf
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: term1, term2, term3
  real(8) :: f1, f2

  ! ====================================

  ! initialize
  pstrnf= 0.0d0


  ! compute f1 and f2
  term1= -27.0d0 / 2.0d0 * etatri * ( etatri**2 - 1.0d0 / 3.0d0 )
  term2= 1.0d0 / 3.0d0 * asin( term1 )
  f1= cos(term2)
  f2= sin(term2)
  
  ! compute effective plastic fracture strain
  term1= a / c2
  term2= dsqrt( (1.0d0 + c1**2) / 3.0d0 ) * f1 
  term3= c1 * ( etatri + f2 / 3.0d0 )
  
  pstrnf = (term1 * ( term2 + term3 ) )**(-1.0d0 / n)
  


  return
end subroutine xwpstrnfj2pstrs
