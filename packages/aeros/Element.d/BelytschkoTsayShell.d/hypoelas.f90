! ==================================
! hypo elasticity  
! ==================================
!      type                  name                              argument
!      ----                  ----                              --------
! 1.  subroutine          getlameconst                     (young,poiss, lamd,mu)
! 2.  subroutine          getbulkconst                     (young,poiss, bulk)
! 3.  subroutine          getelsvoit2d                     (optpty,young,poiss, elsvoit2d)
! 4.  subroutine          getelsvoit3d                     (young,poiss, elsvoit3d)
! 5.  subroutine          getelstens                       (nndex,lamd,mu, elstens)
! 6.  real(8) function    elscomp                          (lamd,mu, i,j,k,l)
!
! =========================================================================================================



subroutine getlameconst(young,poiss, lamd,mu)
  !=======================================================================
  !  getlameconst = compute lame material constant and bulk modulus
  !
  !                note:
  !                ----
  !                malvern, pp. 293
  !                introduction to the mechanics of a continuous medium
  !
  !                TB's book, pp.  228
  !                nonlinear finite elements for continua and structures
  !
  !                lamd= ( e * nu ) /( (1+nu)*(1-2*nu) )
  !                
  !                mu= e / (2*(1+nu))
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  young : young's modulus
  !
  !  poiss : poisson's ratio
  !  
  !  output:
  !  ------
  !  lamd : 1st lame constant lambda
  !
  !  mu : 2nd lame constant mu (= shear modulus G)
  !                            
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), intent(in) :: young, poiss

  real(8), intent(out) :: lamd, mu
  ! ====================================

  ! initialize
  lamd= 0.0d0
  mu= 0.0d0


  ! 1st lame constant
  lamd= (poiss * young) / ( (1.0d0 + poiss) * (1.0d0 - 2.0d0*poiss) )
  ! lamd= ( e * nu ) /( 1 - nu**2 ) <- plane stress

  ! 2nd lame constant: shear modulus 
  mu= young / ( 2.0d0 * (1.0d0+poiss) )



  return
end subroutine getlameconst





subroutine getbulkconst(young,poiss, bulk)
  !=======================================================================
  !  getbulkconst = compute bulk modulus
  !
  !                note:
  !                ----
  !                malvern, pp. 293
  !                introduction to the mechanics of a continuous medium
  !
  !                TB's book, pp.  228
  !                nonlinear finite elements for continua and structures
  !
  !                bulk= lamd + (2/3)*mu
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  young : young's modulus
  !
  !  poiss : poisson's ratio
  !  
  !  output:
  !  ------
  !  bulk : bulk modulus k
  !                            
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), intent(in) :: young, poiss

  real(8), intent(out) :: bulk
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: lamd, mu

  ! ====================================

  ! initialize
  bulk= 0.0d0


  ! bulk modulus
  bulk= young / ( 3.0d0 * ( 1.0d0 - 2.0d0 * poiss) )



  return
end subroutine getbulkconst





subroutine getelsvoit2d(optpty,young,poiss, elsvoit2d)
  !=======================================================================
  !  getelsvoit2d= compute voight form 2d isotropic elsatic moduli for 
  !
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optpty : 2d problem type
  !           optpty=1 : plane-stress
  !           optpty=2 : plane-strain
  ! 
  !  young : young's modulus
  !
  !  poiss : poisson's ratio
  !  
  !  output:
  !  ------
  !  elsvoit(3,3) : constitutive law voight form matrix
  !                            
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optpty
  real(8), intent(in) :: young, poiss

  real(8), dimension(3,3),intent(out) :: elsvoit2d
  ! ====================================
  ! local variable
  ! ==============
  character(len=10) :: handler

  real(8) :: const
  ! ====================================

  ! initialize
  elsvoit2d(:,:)= 0.0d0
  

  select case(optpty)
  case(1) ! plane stress case
     const=young/(1.0d0-poiss**2)

     elsvoit2d(1,1)=const
     elsvoit2d(1,2)=const*poiss
     elsvoit2d(1,3)=0.0d0

     elsvoit2d(2,1)=const*poiss
     elsvoit2d(2,2)=const
     elsvoit2d(2,3)=0.0d0

     elsvoit2d(3,1)=0.0d0
     elsvoit2d(3,2)=0.0d0
     elsvoit2d(3,3)=(1.0d0-poiss)*const/2.0d0

  case(2) ! plane strain case
     const=young/((1.0d0+poiss)*(1.0d0-2.0d0*poiss))

     elsvoit2d(1,1)=const*(1.0d0-poiss)
     elsvoit2d(1,2)=const*poiss
     elsvoit2d(1,3)=0.0d0

     elsvoit2d(2,1)=const*poiss
     elsvoit2d(2,2)=const*(1.0d0-poiss)
     elsvoit2d(2,3)=0.0d0

     elsvoit2d(3,1)=0.0d0
     elsvoit2d(3,2)=0.0d0
     elsvoit2d(3,3)=const*(1.0d0-2.0d0*poiss)/2.0d0

  end select



  return
end subroutine getelsvoit2d





subroutine getelsvoit3d(young,poiss, elsvoit3d)
  !=======================================================================
  !  getelsvoit3d= compute voight form 3d isotropic elsatic moduli for 
  !
  !                note:
  !                ----
  !                malvern, pp. 290
  !                introduction to the mechanics of a continuous medium
  !
  !                TB's book, pp.  228
  !                nonlinear finite elements for continua and structures
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  young : young's modulus
  !
  !  poiss : poisson's ratio
  !  
  !  output:
  !  ------
  !  elsvoit3d(6,6) : constitutive law voight form matrix
  !                            
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), intent(in) :: young, poiss

  real(8), dimension(6,6),intent(out) :: elsvoit3d
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: lamd, mu

  ! ====================================

  ! initialize
  elsvoit3d(:,:)= 0.0d0
  

  ! get lame constant
  call getlameconst(young,poiss, lamd,mu)
     ! input : e,nu
     ! output : lamd,mu

  ! ------------------------------------
  elsvoit3d(1,1)= lamd + 2.0d0 * mu

  elsvoit3d(2,2)= elsvoit3d(1,1)
  elsvoit3d(3,3)= elsvoit3d(1,1)

  ! ------------------------------------
  elsvoit3d(1,2)= lamd

  elsvoit3d(1,3)= elsvoit3d(1,2)
  elsvoit3d(2,1)= elsvoit3d(1,2)
  elsvoit3d(2,3)= elsvoit3d(1,2)
  elsvoit3d(3,1)= elsvoit3d(1,2)
  elsvoit3d(3,2)= elsvoit3d(1,2)

  ! ------------------------------------
  elsvoit3d(4,4)= mu

  elsvoit3d(5,5)= elsvoit3d(4,4)
  elsvoit3d(6,6)= elsvoit3d(4,4)
  ! ------------------------------------



  return
end subroutine getelsvoit3d





subroutine getelstens(nndex,lamd,mu, elstens)
  !=======================================================================
  !  getelstens = compute 4th order isotropic linear elastic moduli tensor C_ijkl
  !
  !               C_ijkl= lamd*( d_ij*d_kl ) + mu*( d_ik*d_jk + d_il*d_jk )
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nndex : dimension of 4th order tensor
  !
  !  lamd, mu : lame constant
  !
  !  output:
  !  ------
  !  elstens(nndex,nndex,nndex,nndex) : 4th order elastic moduli tensor C_ijkl
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nndex
  real(8), intent(in) :: lamd, mu 

  real(8), dimension(nndex,nndex,nndex,nndex), intent(out) :: elstens
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: cijkl
  real(8) :: elscomp

  integer :: i, j, k, l ! loop index
  ! ====================================

  ! initialize
  elstens(:,:,:,:)= 0.0d0


  ! compute tensorial form 4th order elastic moduli:
  ! -----------------------------------------------
  ! c_ijkl= lamd*( d_ij*d_kl ) + mu*( d_ik*d_jk + d_il*d_jk )
  do i=1, nndex
     do j=1, nndex
        do k=1, nndex
           do l=1, nndex

              elstens(i,j,k,l)= elscomp(lamd,mu, i,j,k,l)

           end do
        end do
     end do
  end do



  return
end subroutine getelstens





real(8) function elscomp(lamd,mu, i,j,k,l)
  !=======================================================================
  !  elscomp = compute component of 4th order isotropic elastic moduli tensor
  !
  !            isotropic material:
  !            ------------------
  !            c_ijkl= lamd*( d_ij*d_kl ) + mu*( d_ik*d_jk + d_il*d_jk )
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  lamd, mu : lame constant
  !
  !  i, j, k, l : index of c_ijkl
  !
  !  output:
  !  ------
  !  elscomp : isotropic linear elastic tensor component
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), intent(in) :: lamd, mu

  integer, intent(in) :: i, j, k, l
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: krodelt
  real(8) :: dij, dkl, dik, djl, dil, djk
  ! ====================================

  ! initialize
  elscomp= 0.0d0


  ! get kronelker delta value
  dij= krodelt(i,j) ! d_ij
  dkl= krodelt(k,l) ! d_kl
  dik= krodelt(i,k) ! d_ik
  djl= krodelt(j,l) ! d_jl
  dil= krodelt(i,l) ! d_il
  djk= krodelt(j,k) ! d_jk

  ! compute elastic moduli component: c_ijkl= lamd*( d_ij*d_kl ) + mu*( d_ik*d_jk + d_il*d_jk )
  ! --------------------------------
  elscomp= lamd*( dij*dkl ) + mu*( dik*djl + dil*djk )



  return
end function elscomp
