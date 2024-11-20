! ======================================
! mechanical tensors
! ======================================
!      type                  name                              argument
!      ----                  ----                              --------
! 1.  real(8) function       krodelt         (i,j)
! 2.  real(8) function       sprdts          (opt,nndex,atens,btens)
! 3.  subroutine             sprdtt          (opt,nndex,atens,btens, ctens)
! 4.  subroutine             dyadt           (nndex,atens,btens, ctens)
! 5.  subroutine             symtens         (nndex,atens, btens)
! 6.  subroutine             antisymtens     (nndex,atens, btens)
! 7.  real(8) function       tr              (nndex, atens)
! 8.  real(8) function       invrnt3d        (opt,strstens3d)
! 9.  integer function       nvoit           (ndime)
! 10. subroutine             ind2voit        (ntype,ndime,tensor, voigt)
! 11. subroutine             voit2ind        (ntype,ndime,voigt, tensor)
! 12. subroutine             kir2caucy       (nndex,ftens,kirtens, sigtens)        ! kir -> cauchy
! 13. subroutine             caucy2kir       (nndex,ftens,sigtens, kirtens)        ! cauchy -> kir
! 14. subroutine             caucy2pk2       (nndex,ftens,sigtens, pk2tens)        ! cauchy -> pk2
! 15. subroutine             pk22caucy       (nndex,ftens,pk2tens, sigtens)        ! pk2 -> cauchy
! 16. subroutine             pk22nomi        (nndex,ftens,pk2tens, ptens)          ! pk2 -> nominal
! 17. subroutine             kir2pk2         (nndex,ftens,kirtens, pk2tens)        ! kirchoff -> pk2
! 18. subroutine             getftens        (ndime,nnode,edisp,cartd, ftens)      ! deformation gradient tensor F
! 19. subroutine             getetens        (ndime,nnode,edisp,cartd, etens)      ! green strain tensor E
! 20. subroutine             getltens        (ndime,nnode,evel,cartd,ftens, ltens) ! velocity gradient tensor L
! 21. subroutine             getrotens2d     (ndime,theta, rotens2d)
! 22. subroutine             getdevtens      (nndex,tens, devtens)                 ! deviatoric tensor
! 23. subroutine             geteffstrs      (nndex,devstrstens, effstrs)          ! effective stress
! 24. subroutine             getj2dirtens    (nndex,devstrstens,effstrs, pdirtens) ! j2 flow direction
! 25. subroutine             updjgeo         (nndex,midltens,delt, kirtens)        ! jaumann rate geo update
! 26. subroutine             updhwgeo        (nndex,midltens,delt, strstens)       ! hughes-winget geo update
! 27. subroutine             voit3d2tens2d4  (voit3d, tens2d)                      ! vonvert 4th order tensor
! 28. subroutine             getcstrtens3d4  (sigtens3d, cstrtens3d)
! 29. subroutine             voit2ind3d4     (voigt3d4, tensor3d4)
! 30. subroutine             ctantenscj2ct3d (ctantenscj3d,sigtens3d, ctantensct3d)
! 31. subroutine             getatens3d4     (ctantens3d,sigtens3d, atens3d)
! 32. subroutine             curn2inivec     (ndime,ftens,veccurn, vecini)
! 33. subroutine             curn2iniang2d   (ftens,angcurn, angini)
! 34. real(8) function       trsprs          (nsprs,krow,kcol,sprsmat)             ! trace of a sparse matrix
! 35. subroutine             getbtens        (ndime,ftens, btens)                  ! left cauchy-green deformation tensor B
! 36. subroutine             getbbartens     (ndime,ftens, bbartens)               ! volume-preserving part of left cauchy-green deformation tensor B
! 37. subroutine             getctens        (ndime,ftens, ctens)                  ! right cauchy-green deformation tensor c
!
! =========================================================================================================



real(8) function krodelt(i,j)
  !=======================================================================
  !  krodelt = kronecker delta function value
  !
  !            del_ij = 1   if i==j
  !                   = 0   if i/=j
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  i, j : i and j index of kronecker delta
  !
  !  output:
  !  ------
  !  krodelt
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: i, j
  ! ====================================

  ! initialize
  krodelt= 0.0d0

  ! kronecker delta value
  if ( i == j ) then

     krodelt= 1.0d0

  end if



  return
end function krodelt





real(8) function sprdts(opt,nndex,atens,btens)
  !=======================================================================
  !  sprdts = compute scalar product of two given 2nd order tensor
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  opt : option handler
  !        opt=0: sprdts= A:B = A_ij * B_ij
  !
  !        opt=1: sprdts= A..B= A_ij * B_ji
  !
  !  nndex : dimension of a and b tensor
  !
  !  atens(nndex, nndex) : A tensor
  !
  !  btens(nndex, nndex) : B tensor
  !
  !  output:
  !  ------
  !  sprdts : result of scalar product of two tensor
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: opt
  integer, intent(in) :: nndex
  real(8), dimension(nndex,nndex), intent(in) :: atens
  real(8), dimension(nndex,nndex), intent(in) :: btens
  ! ====================================
  ! local variable
  ! ==============

  ! loop index
  integer :: i, j
  ! ====================================

  ! initialize
  sprdts= 0.0d0


  select case(opt)
  case(0)

     ! compute sclar product: sprdts= a:b= a_ij * b_ij
     ! ---------------------
     do i=1, nndex
        do j=1, nndex
           sprdts= sprdts + atens(i,j)*btens(i,j)
        end do
     end do

  case(1)
     ! compute sclar product: sprdts= a..b= a_ij * b_ji
     ! ---------------------
     do i=1, nndex
        do j=1, nndex
           sprdts= sprdts + atens(i,j)*btens(j,i)
        end do
     end do

  end select



  return
end function sprdts





subroutine sprdtt(opt,nndex,atens,btens, ctens)
  !=======================================================================
  !  sprdtt = compute scalar product between 4th order and 2nd order tensor
  !
  !              note:
  !              ----
  !              result is 2nd order tensor
  !
  !              opt=0: C= A:B
  !                     C_ij= A_ijkl * B_kl
  !
  !              opt=1: C= A..B
  !                     C_ij= A_ijkl * B_lk
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  opt : option handler
  !        opt=0: C= A:B
  !               C_ij= A_ijkl * B_kl
  !
  !        opt=1: C= A..B
  !               C_ij= A_ijkl * B_lk
  !
  !  nndex : dimension of tensor
  !
  !  atens(nndex,nndex,nndex,nndex) : 4th order tensor
  !
  !  btens(nndex,nndex) : 2nd order tensor
  !
  !  output:
  !  ------
  !  ctens(nndex,nndex) : 2nd order tensor
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: opt
  integer, intent(in) :: nndex
  real(8), dimension(nndex,nndex,nndex,nndex), intent(in) :: atens
  real(8), dimension(nndex,nndex), intent(in) :: btens

  real(8), dimension(nndex,nndex), intent(out) :: ctens
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: sum

  integer :: i, j, k, l ! loop index
  ! ====================================

  ! initialize
  ctens(:,:)=0.0d0


  select case(opt)
  case(0)
     ! c_ij= a_ijkl * b_kl
     ! -------------------
     do i=1, nndex
        do j=1, nndex

           sum= 0.0d0 ! initialize
           do k=1, nndex
              do l=1, nndex
                 sum= sum + atens(i,j,k,l)*btens(k,l)
              end do
           end do        
           ctens(i,j)= sum

        end do
     end do


  case (1)
     ! c_ij= a_ijkl * b_lk
     ! -------------------
     do i=1, nndex
        do j=1, nndex

           sum= 0.0d0 ! initialize
           do k=1, nndex
              do l=1, nndex
                 sum= sum + atens(i,j,k,l)*btens(l,k)
              end do
           end do        
           ctens(i,j)= sum

        end do
     end do 

  end select



  return
end subroutine sprdtt





subroutine dyadt(nndex,atens,btens, ctens)
  !=======================================================================
  !  dyadt = compute dyadic of two given 2nd order tensor
  !
  !          note:
  !          ----
  !          result is 4th order tensor
  !
  !          c= a (x) b
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nndex : dimension of tensor
  !
  !  atens(nndex,nndex) : 2nd order tensor
  !
  !  btens(nndex,nndex) : 2nd order tensor
  !
  !  output:
  !  ------
  !  ctens(nndex,nndex,nndex,nndex) : 4th order tensor
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nndex
  real(8), dimension(nndex,nndex), intent(in) :: atens
  real(8), dimension(nndex,nndex), intent(in) :: btens

  real(8), dimension(nndex,nndex,nndex,nndex), intent(out) :: ctens
  ! ====================================
  ! local variable
  ! ==============
  integer :: i, j, k, l ! loop index
  ! ====================================

  ! initialize
  ctens(:,:,:,:)= 0.0d0 

  ! compute dyad: C= A (x) B
  ! ------------
  do i=1, nndex
     do j=1, nndex
        do k=1, nndex
           do l=1, nndex
              ctens(i,j,k,l)= atens(i,j)*btens(k,l)
           end do
        end do        
     end do
  end do



  return
end subroutine dyadt





subroutine symtens(nndex,atens, btens)
  !=======================================================================
  !  symtens = compute symmetric part of given tensor
  !
  !          b= 1/2 * ( a + a^T)
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nndex : dimension of tensor
  !
  !  atens(nndex,nndex) : 2nd order tensor
  !
  !  output:
  !  ------
  !  btens(nndex,nndex) : symmetric part of A tensor
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nndex
  real(8), dimension(nndex,nndex), intent(in) :: atens

  real(8), dimension(nndex,nndex), intent(out) :: btens
  ! ====================================
  ! local variable
  ! ==============
  integer :: i, j ! loop index
  ! ====================================

  ! initialize
  btens(:,:)= 0.0d0 


  ! compute symmetric part of atens: b= 1/2 * ( a + a^T )
  ! -------------------------------
  do i=1, nndex
     do j=1, nndex
        btens(i,j)= 0.50d0 * ( atens(i,j) + atens(j,i) )
     end do
  end do



  return
end subroutine symtens





subroutine antisymtens(nndex,atens, btens)
  !=======================================================================
  !  antisymtens = compute anti-symmetric part of given tensor
  !
  !          b= 1/2 * ( a - a^T)
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nndex : dimension of tensor
  !
  !  atens(nndex,nndex) : 2nd order tensor
  !
  !  output:
  !  ------
  !  btens(nndex,nndex) : anti-symmetric part of a tensor
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nndex
  real(8), dimension(nndex,nndex), intent(in) :: atens

  real(8), dimension(nndex,nndex), intent(out) :: btens
  ! ====================================
  ! local variable
  ! ==============
  integer :: i, j ! loop index
  ! ====================================

  ! initialize
  btens(:,:)= 0.0d0 


  ! compute anti-symmetric part of atens: b= 1/2 * ( a - a^T )
  ! ------------------------------------
  do i=1, nndex
     do j=1, nndex
        btens(i,j)= 0.50d0 * ( atens(i,j) - atens(j,i) )
     end do
  end do

  return
end subroutine antisymtens





real(8) function tr(nndex, atens)
  !=======================================================================
  !  tr = compute trace of given tensor
  !
  !           tr= Tr(a) = a_kk
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nndex : dimension of a and b tensor
  !
  !  atens(nndex, nndex) : A tensor
  !
  !  output:
  !  ------
  !  tr : trace of tensor
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nndex
  real(8), dimension(nndex,nndex), intent(in) :: atens
  ! ====================================
  ! local variable
  ! ==============
  
  ! loop index
  integer :: i 
  ! ====================================

  ! initialize
  tr= 0.0d0


  ! compute trace of given tensor
  do i=1, nndex

     tr= tr + atens(i,i)

  end do



  return
end function tr





real(8) function invrnt3d(opt,strstens3d)
  !=======================================================================
  !  invrnt3d = compute three invariants of given 2nd order 3x3 tensor
  !
  !              note:
  !              ----
  !              malvern, pp. 89
  !              introduction to the mechanics of a continuous medium
  !              
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  opt : invariant option handler
  !
  !  strstens3d(3,3) : 2nd order 3d tensor
  !
  !  output:
  !  ------
  !  invrnt3d : tensor invariant
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: opt
  real(8), dimension(3,3), intent(in) :: strstens3d
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: tr, sprdts, getdet

  ! ====================================

  ! initialize
  invrnt3d= 0.0d0


  ! compute tensor invariant
  select case(opt)
  case(1) ! 1st invariant : j1= tr(tens)
     invrnt3d= tr(3,strstens3d)

  case(2) ! 2nd invariant : j2= 0.50d0 * ( tens:tens - tr(tens)**2 )
     invrnt3d= 0.50d0 * ( sprdts(0,3,strstens3d,strstens3d) - tr(3,strstens3d)**2 )

  case(3) ! 3rd invariant : j3= det(tens)
     invrnt3d= getdet(3,strstens3d)

  case default
     write(*,*) "wrong option : invrnt3d"
     write(nout5,*) "wrong option : invrnt3d"
     stop

  end select



  return
end function invrnt3d





integer function nvoit(ndime)
  !=======================================================================
  !  nvoit= get the total number of voight form components
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : dimension of tensor
  !
  !  output:
  !  ------
  !  nvoit
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime
  ! ====================================

  nvoit= ndime * ( ndime + 1 ) / 2



  return
end function nvoit



subroutine ind2voit(ntype,ndime,tensor, voigt)
  !=======================================================================
  !  ind2voit= convert indicial form tensor to voit form matrix
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ntype : define the type of considered quantity
  !          ntype=0 : kinetic (ex. stress)
  !          ntype=1 : kinematic (ex. strain)
  !
  !  ndime : dimension of inditial matrix
  !
  !  tensor(ndime,ndime) : inditial format tensor : ndime by ndime
  !
  !  output:
  !  ------
  !  voigt(ndime*(ndime+1)/2,1) : voigt format matrix : nvoit by 1
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ntype
  integer, intent(in) :: ndime
  real(8), dimension(ndime,ndime), intent(in) :: tensor

  real(8), dimension(ndime*(ndime+1)/2,1), intent(out) :: voigt
  ! ====================================
  ! local variable
  ! ==============
  character(len=20) :: handler
  ! ====================================

  ! initialize
  voigt(:,:)=0.0d0


  ! set the handler
  if(ntype==0.and.ndime==2) handler='2d_kinetic'
  if(ntype==0.and.ndime==3) handler='3d_kinetic'
  if(ntype==1.and.ndime==2) handler='2d_kinematic'
  if(ntype==1.and.ndime==3) handler='3d_kinematic'

  select case(handler)
  case('2d_kinetic') ! 2d kinetic
     voigt(1,1)= tensor(1,1)
     voigt(2,1)= tensor(2,2)
     voigt(3,1)= tensor(1,2)

  case('3d_kinetic') ! 3d kinetic
     voigt(1,1)= tensor(1,1)
     voigt(2,1)= tensor(2,2)
     voigt(3,1)= tensor(3,3)
     voigt(4,1)= tensor(2,3)
     voigt(5,1)= tensor(1,3)
     voigt(6,1)= tensor(1,2)

  case('2d_kinematic') ! 2d kinematic
     voigt(1,1)= tensor(1,1)
     voigt(2,1)= tensor(2,2)
     voigt(3,1)= 2.0d0 * tensor(1,2)

  case('3d_kinematic') ! 3d kinematic
     voigt(1,1)= tensor(1,1)
     voigt(2,1)= tensor(2,2)
     voigt(3,1)= tensor(3,3)
     voigt(4,1)= 2.0d0 * tensor(2,3)
     voigt(5,1)= 2.0d0 * tensor(1,3)
     voigt(6,1)= 2.0d0 * tensor(1,2)

  case default
     write(*,*) "invalid option: ind2voit"
     write(nout5,*) "invalid option: ind2voit"
     stop

  end select



  return
end subroutine ind2voit





subroutine voit2ind(ntype,ndime,voigt, tensor)
  !=======================================================================
  !  voit2ind= convert voit form matrix indicial form tensor
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ntype : define the type of considered quantity
  !          ntype=0 : kinetic (ex. stress)
  !          ntype=1 : kinematic (ex. strain)
  !
  !  ndime : dimension of inditial matrix
  !
  !  voigt(ndime*(ndime+1)/2,1) : voigt format matrix : nvoit by 1
  !
  !  output:
  !  ------
  !  tensor(ndime,ndime) : inditial format tensor : ndime by ndime
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ntype
  integer, intent(in) :: ndime
  real(8), dimension(ndime*(ndime+1)/2,1), intent(in) :: voigt

  real(8), dimension(ndime,ndime), intent(out) :: tensor
  ! ====================================
  ! local variable
  ! ==============
  character(len=20) :: handler
  ! ====================================

  ! initialize
  tensor(:,:)=0.0d0


  ! set the handler
  if(ntype==0.and.ndime==2) handler='2d_kinetic'
  if(ntype==0.and.ndime==3) handler='3d_kinetic'
  if(ntype==1.and.ndime==2) handler='2d_kinematic'
  if(ntype==1.and.ndime==3) handler='3d_kinematic'

  select case(handler)
  case('2d_kinetic') ! 2d kinetic
     tensor(1,1)= voigt(1,1)
     tensor(1,2)= voigt(3,1)
     tensor(2,2)= voigt(2,1)
     tensor(2,1)= tensor(1,2) ! symmetric parts

  case('3d_kinetic') ! 3d kinetic
     tensor(1,1)= voigt(1,1)
     tensor(2,2)= voigt(2,1)
     tensor(3,3)= voigt(3,1)
     tensor(2,3)= voigt(4,1)
     tensor(1,3)= voigt(5,1)
     tensor(1,2)= voigt(6,1)
     tensor(3,2)= tensor(2,3) ! symmetric parts
     tensor(3,1)= tensor(1,3)
     tensor(2,1)= tensor(1,2)

  case('2d_kinematic') ! 2d kinematic
     tensor(1,1)= voigt(1,1)
     tensor(2,2)= voigt(2,1)
     tensor(1,2)= 0.50d0 * voigt(3,1)
     tensor(2,1)= tensor(1,2) ! symmetric part

  case('3d_kinematic') ! 3d kinematic
     tensor(1,1)= voigt(1,1)
     tensor(2,2)= voigt(2,1)
     tensor(3,3)= voigt(3,1)
     tensor(2,3)= 0.50d0 * voigt(4,1)
     tensor(1,3)= 0.50d0 * voigt(5,1)
     tensor(1,2)= 0.50d0 * voigt(6,1)
     tensor(3,2)= tensor(2,3) ! symmetric parts
     tensor(3,1)= tensor(1,3)
     tensor(2,1)= tensor(1,2)

  case default
     write(*,*) "invalid option: voit2ind"
     write(nout5,*) "invalid option: voit2ind"
     stop

  end select



  return
end subroutine voit2ind





subroutine caucy2kir(nndex,ftens,sigtens, kirtens)
  !=======================================================================
  !  caucy2kir= convert cauchy stress tensor to kirchhoff stress tensor
  !
  !              kir= det(f) * sig (see, TB's book page 103)
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nndex : dimension of input tensors
  !
  !  ftens(nndex,nndex) : deformation gradient tensor
  !
  !  sigtens(nndex,nndex) : cauchy stress tensor
  ! 
  !  output:
  !  ------
  !  kirtens(nndex,nndex) : kirchoff tensor
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nndex
  real(8), dimension(nndex,nndex), intent(in) :: ftens
  real(8), dimension(nndex,nndex), intent(in) :: sigtens

  real(8), dimension(nndex,nndex), intent(out) :: kirtens 
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: getdet, detf
  ! ====================================

  ! initialize
  kirtens(:,:)= 0.0d0


  ! compute determinant of F tensor: det(f)
  detf= getdet(nndex,ftens)

  ! kir= det(f) * sig
  kirtens(1:nndex,1:nndex)= detf * sigtens(1:nndex,1:nndex)



  return
end subroutine caucy2kir





subroutine caucy2pk2(nndex,ftens,sigtens, pk2tens)
  !=======================================================================
  !  caucy2pk2= convert cauchy stress tensor to 2nd piola-kirchhoff stress tensor
  !
  !              pk2= det(f) f^-1.sig.(f^t)^-1 (see, TB's book page 103)
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nndex : dimension of input tensors
  !
  !  ftens(nndex,nndex) : deformation gradient tensor
  !
  !  sigtens(nndex,nndex) : cauchy stress tensor
  !
  !  output:
  !  ------
  !  pk2tens(nndex,nndex) : 2nd piola-kirchoff tensor
  ! 
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nndex
  real(8), dimension(nndex,nndex), intent(in) :: ftens
  real(8), dimension(nndex,nndex), intent(in) :: sigtens

  real(8), dimension(nndex,nndex), intent(out) :: pk2tens
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(nndex,nndex) :: fttens
  real(8), dimension(nndex,nndex) :: invftens
  real(8), dimension(nndex,nndex) :: invfttens

  real(8), dimension(nndex,nndex) :: temp, temp1 ! tempolary work space

  real(8) :: getdet, detf

  ! loop index
  integer :: index, jndex
  ! ====================================

  ! initialize
  pk2tens(:,:)= 0.0d0


  ! pk2= det(f) f^-1.sig.(f^t)^-1
  ! -----------------------------

  ! compute transpose of deformation gradient tensor: ftens^T
  do index=1, nndex
     do jndex=1, nndex
        fttens(index,jndex)= ftens(jndex,index)
     end do
  end do

  ! compute inverse of deformation gradient tensor: inv(ftens)
  call getinv(nndex,ftens, invftens)
     ! input : nndex,ftens
     ! output : invftens

  ! compute inverse of transposed deformation gradient tensor: inv(ftens^T)
  call getinv(nndex,fttens, invfttens)
     ! input : nndex,fttens
     ! output : invfttens

  ! compute determinant of f tensor: det(f)
  detf= getdet(nndex,ftens)

  ! compute inv(f).sig
  call matprd(nndex,nndex,0, nndex,nndex,0, nndex,nndex, invftens,sigtens, temp)
    ! input : nndex,nndex,0, nndex,nndex,0, nndex,nndex, invftens,sigtens
    ! output : temp

  ! compute (inv(f).sig).inv(f^t)
  call matprd(nndex,nndex,0, nndex,nndex,0, nndex,nndex, temp,invfttens, temp1)
    ! input : nndex,nndex,0, nndex,nndex,0, nndex,nndex, temp,invfttens
    ! output : temp1

  ! multiply det(f)
  if( detf /= 0.0d0 ) then
    pk2tens(1:nndex,1:nndex)= detf * temp1(1:nndex,1:nndex)

  else
    pk2tens(1:nndex,1:nndex)= 0.0d0

  end if



  return
end subroutine caucy2pk2





subroutine kir2caucy(nndex,ftens,kirtens, sigtens)
  !=======================================================================
  !  kir2caucy= convert kirchhoff stress tensor to cauchy stress tensor
  !
  !              sig= det(f)^-1 kir (see, TB's book page 103)
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nndex : dimension of input tensors
  !
  !  ftens(nndex,nndex) : deformation gradient tensor
  !
  !  kirtens(nndex,nndex) : kirchoff tensor
  ! 
  !  output:
  !  ------
  !  sigtens(nndex,nndex) : cauchy stress tensor
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nndex
  real(8), dimension(nndex,nndex), intent(in) :: ftens
  real(8), dimension(nndex,nndex), intent(in) :: kirtens

  real(8), dimension(nndex,nndex), intent(out) :: sigtens
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: getdet, detf
  ! ====================================

  ! initialize
  sigtens(:,:)=0.0d0


  ! sig= det(f)^-1 kir
  ! ----------------
  ! compute determinant of F tensor: det(f)
  detf= getdet(nndex,ftens)

  ! multiply 1/det(f)
  if( detf /= 0.0d0 ) then
    sigtens(1:nndex,1:nndex)= (1.0d0/detf) * kirtens(1:nndex,1:nndex)

  else
    sigtens(1:nndex,1:nndex)= 0.0d0

  end if



  return
end subroutine kir2caucy





subroutine pk22caucy(nndex,ftens,pk2tens, sigtens)
  !=======================================================================
  !  pk22caucy= convert 2nd piola-kirchhoff stress tensor to cauchy stress tensor
  !
  !              sig= det(f)^-1 f.pk2.f^t (see, TB's book page 103)
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nndex : dimension of input tensors
  !
  !  ftens(nndex,nndex) : deformation gradient tensor
  !
  !  pk2tens(nndex,nndex) : 2nd piola-kirchoff tensor
  ! 
  !  output:
  !  ------
  !  sigtens(nndex,nndex) : cauchy stress tensor
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nndex
  real(8), dimension(nndex,nndex), intent(in) :: ftens
  real(8), dimension(nndex,nndex), intent(in) :: pk2tens

  real(8), dimension(nndex,nndex), intent(out) :: sigtens
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(nndex,nndex) :: temp, temp1 ! tempolary work space

  real(8) :: getdet, detf
  ! ====================================

  ! initialize
  sigtens(:,:)= 0.0d0


  ! sig= det(f)^-1 f.pk2.f^t
  ! ----------------
  ! compute f.pk2
  call matprd(nndex,nndex,0, nndex,nndex,0, nndex,nndex, ftens,pk2tens, temp)
    ! input : nndex,nndex,0, nndex,nndex,0, nndex,nndex, ftens,pk2tens
    ! output : temp

  ! compute (f.pk2).f^t
  call matprd(nndex,nndex,0, nndex,nndex,1, nndex,nndex, temp,ftens, temp1)
    ! input : nndex,nndex,0, nndex,nndex,1, nndex,nndex, temp,ftens
    ! output : temp1

  ! compute determinant of F tensor: det(f)
  detf= getdet(nndex,ftens)

  ! multiply 1/det(f)
  if( detf /= 0.0d0 ) then
    sigtens(1:nndex,1:nndex)= (1.0d0/detf) * temp1(1:nndex,1:nndex)

  else
    sigtens(1:nndex,1:nndex)= 0.0d0

  end if



  return
end subroutine pk22caucy





subroutine pk22nomi(nndex,ftens,pk2tens, ptens)
  !=======================================================================
  !  pk22nomi= convert 2nd piola-kirchhoff (pk2) stress to nominal stress tensor (p)
  !
  !              p=pk2.f^t (see, TB's book page 103)
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nndex : dimension of input tensors
  !
  !  ftens(nndex,nndex) : deformation gradient tensor
  !
  !  pk2tens(nndex,nndex) : 2nd piola-kirchoff tensor
  ! 
  !  output:
  !  ------
  !  ptens(nndex,nndex) : nominal stress tensor
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nndex
  real(8), dimension(nndex,nndex), intent(in) :: ftens
  real(8), dimension(nndex,nndex), intent(in) :: pk2tens

  real(8), dimension(nndex,nndex), intent(out) :: ptens
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(nndex,nndex) :: temp ! tempolary work space
  ! ====================================

  ! initialize
  ptens(:,:)=0.0d0


  ! p= pk2.f^t
  ! ---------
  ! compute pk2.f^t
  call matprd(nndex,nndex,0, nndex,nndex,1, nndex,nndex, pk2tens,ftens, temp)
    ! input : nndex,nndex,0, nndex,nndex,1, nndex,nndex, pk2tens,ftens
    ! output : temp

  ! copy result to ptens
  ptens(1:nndex,1:nndex)= temp(1:nndex,1:nndex)



  return
end subroutine pk22nomi





subroutine kir2pk2(nndex,ftens,kirtens, pk2tens)
  !=======================================================================
  !  kir2pk2= convert kirchoff stress tensor to 2nd piola-kirchhoff stress tensor
  !
  !              pk2= f^-1.kir.f^-t (see, TB's book page 103)
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nndex : dimension of input tensors
  !
  !  ftens(nndex,nndex) : deformation gradient tensor
  !
  !  kirtens(nndex,nndex) : kirchoff stress tensor
  ! 
  !  output:
  !  ------
  !  pk2tens(nndex,nndex) : 2nd piola-kirchoff tensor
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nndex
  real(8), dimension(nndex,nndex), intent(in) :: ftens
  real(8), dimension(nndex,nndex), intent(in) :: kirtens

  real(8), dimension(nndex,nndex), intent(out) :: pk2tens
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(nndex,nndex) :: invftens, fttens, invfttens 
  real(8), dimension(nndex,nndex) :: temp, temp1 ! tempolary work space

  ! loop index
  integer :: i, j
  ! ====================================

  ! initialize
  pk2tens(:,:)= 0.0d0


  ! compute inverse of ftens
  call getinv(nndex,ftens, invftens)
     ! input : nndex,ftens
     ! output : invftens

  ! compute transpose of ftens
  do i=1, nndex
     do j=1, nndex
        fttens(i,j)= ftens(j,i)
     end do
  end do

  ! compute inverse of transpose of ftens
  call getinv(nndex,fttens, invfttens)
     ! input : nndex,fttens
     ! output : invfttens

  ! pk2= f^-1.kir.f^-t
  ! ------------------
  ! compute f^-1.kir
  call matprd(nndex,nndex,0, nndex,nndex,0, nndex,nndex, invftens,kirtens, temp)
     ! input : nndex,nndex,0, nndex,nndex,0, nndex,nndex, invftens,kirtens
     ! output : temp

  ! compute (f^-1.kir).f^-t
  call matprd(nndex,nndex,0, nndex,nndex,0, nndex,nndex, temp,invfttens, temp1)
     ! input : nndex,nndex,0, nndex,nndex,0, nndex,nndex, temp,invfttens
     ! output : temp1

  ! copy result to pk2tens
  pk2tens(1:nndex,1:nndex)= temp1(1:nndex,1:nndex)



  return
end subroutine kir2pk2





subroutine getftens(ndime,nnode,edisp,cartd, ftens)
  !=======================================================================
  !  getftens= compute deformation gradient tensor
  !
  !            f_ij= a x_i / a X_j
  !
  !            note: to prevent rundoff error we use h tensor (see TB's book page 201) 
  !            ----
  !            f_ij= a x_i / a X_j     <- let, x= u+X
  !                = a (u_i + X_i)/ a X_j
  !                = a u_i / a X_j + del_ij
  !                = h_ij + del_ij
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nnode : total number of node per element
  !
  !  nndex : dimension
  !
  !  edisp(ndime,nnode) : element nodal displacement
  !           edisp(1,i) : u_ix
  !           edisp(2,i) : u_iy
  !           edisp(3,i) : u_iz
  !
  !  cartd(ndime,nnode) : cartesian derivatives of sahpe function
  !           cartd(1,i) : dn_i/dX
  !           cartd(2,i) : dn_i/dY
  !           cartd(3,i) : dn_i/dZ
  !
  !  output:
  !  ------
  !  ftens(ndime,ndime) : deformation gradient tensor
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime, nnode
  real(8), dimension(ndime,nnode), intent(in) :: edisp
  real(8), dimension(ndime,nnode), intent(in) :: cartd

  real(8), dimension(ndime,ndime), intent(out) :: ftens
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(ndime,ndime) :: htens ! gradient of displacement w.r.t X-Y
  real(8) :: krodelt

  ! loop index
  integer :: idime, jdime
  ! ====================================

  ! initialize
  ftens(:,:)= 0.0d0


  ! compute h tensor: h_ij= a u_i / a X_j   
  ! ----------------
  call matprd(ndime,nnode,0, ndime,nnode,1, ndime,ndime, edisp,cartd, htens)
    ! input : ndime,nnode,0, ndime,nnode,1, ndime,ndime, edisp,cartd
    ! output : htens


  ! compute deformation gradient tensor: f_ij= h_ij + del_ij
  ! -----------------------------------
  do idime=1, ndime
     do jdime=1, ndime

        ftens(idime,jdime)= htens(idime,jdime) + krodelt(idime,jdime)

     end do
  end do



  return
end subroutine getftens





subroutine getetens(ndime,nnode,edisp,cartd, etens)
  !=======================================================================
  !  getftens= green (green-lagrange) strain tensor
  !
  !            e_ij= 1/2*( a u_i / a X_j + a u_j / a X_i + (a u_k / a X_i)*(a u_k / a X_j) )
  !
  !            note: to prevent rundoff error we use h tensor    (see TB's book page 93, 201) 
  !            ----
  !            e= 1/2*(h + h^T + h^t.h )
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nnode : total number of node per element
  !
  !  nndex : dimension
  !
  !  edisp(ndime,nnode) : element nodal displacement
  !           edisp(1,i) : u_ix
  !           edisp(2,i) : u_iy
  !           edisp(3,i) : u_iz
  !
  !  cartd(ndime,nnode) : cartesian derivatives of sahpe function
  !           cartd(1,i) : dn_i/dX
  !           cartd(2,i) : dn_i/dY
  !           cartd(3,i) : dn_i/dZ
  !
  !  output:
  !  ------
  !  etens(ndime,ndime) : green (green-lagrange) strain tensor
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime, nnode 
  real(8), dimension(ndime,nnode), intent(in) :: edisp
  real(8), dimension(ndime,nnode), intent(in) :: cartd

  real(8), dimension(ndime,ndime), intent(out) :: etens
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(ndime,ndime) :: htens ! gradient of displacement w.r.t X-Y
  real(8), dimension(ndime,ndime) :: hthtens ! h^t.h tensor

  ! loop index
  integer :: idime, jdime
  ! ====================================

  ! initialize
  etens(:,:)= 0.0d0


  ! compute h tensor: h_ij= a u_i / a X_j   
  ! ----------------
  call matprd(ndime,nnode,0, ndime,nnode,1, ndime,ndime, edisp,cartd, htens)
    ! input : ndime,nnode,0, ndime,nnode,1, ndime,ndime, edisp,cartd
    ! output : htens

  ! compute h^t.h
  call matprd(ndime,ndime,1, ndime,ndime,0, ndime,ndime, htens,htens, hthtens)
    ! input : ndime,ndime,1, ndime,ndime,0, ndime,ndime, htens,htens
    ! output : hthtens

  ! compute green strain tensor: e= 1/2 * ( h + h^T + h^t.h )
  ! ---------------------------
  do idime=1, ndime
     do jdime=1, ndime

        etens(idime,jdime)= 0.50d0 * ( htens(idime,jdime) + htens(jdime,idime) + hthtens(idime,jdime) )

     end do
  end do



  return
end subroutine getetens





subroutine getltens(ndime,nnode,evel,cartd,ftens, ltens)
  !=======================================================================
  !  getltens= compute velocity gradient tensor
  !
  !            l_ij= a v_i / a x_j=  (a v_i / a X_k) * ( a X_k / a x_j)
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nnode : total number of node per element
  !
  !  nndex : dimension
  !
  !  evel(ndime,nnode) : element nodal velocity
  !           evel(1,i) : v_ix
  !           evel(2,i) : v_iy
  !           evel(3,i) : v_iz
  !
  !  cartd(ndime,nnode) : point value of shape function's derivation in X-Y domain
  !           cartd(1,i) : dn_i/dX
  !           cartd(2,i) : dn_i/dY
  !           cartd(3,i) : dn_i/dZ
  !
  !  ftens(ndime,nnode) : deformation gradient tensor
  !
  !  output:
  !  ------
  !  ltens(ndime,ndime) : velocity gradient tensor
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime, nnode
  real(8), dimension(ndime,nnode), intent(in) :: evel
  real(8), dimension(ndime,nnode), intent(in) :: cartd
  real(8), dimension(ndime,ndime), intent(in) :: ftens

  real(8), dimension(ndime,ndime), intent(out) :: ltens
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(ndime,ndime) :: temp, temp1
  real(8), dimension(ndime,ndime) :: invftens
  ! ====================================

  ! initialize
  ltens(:,:)= 0.0d0


  ! compute h tensor: temp_ij= a v_i / a X_j   
  ! ----------------
  call matprd(ndime,nnode,0, ndime,nnode,1, ndime,ndime, evel,cartd, temp)
    ! input : ndime,nnode,0, ndime,nnode,1, ndime,ndime, evel,cartd
    ! output : temp

  ! compute inverse of f tensor: inv(f)=a X_j / a x_i
  call getinv(ndime,ftens, invftens)
    ! input : ndime,ftens
    ! output : invftens

  ! compute velocity gradient tensor: 
  call matprd(ndime,ndime,0, ndime,ndime,0, ndime,ndime, temp,invftens, ltens)
    ! input : ndime,ndime,0, ndime,ndime,0, ndime,ndime, temp,invftens
    ! output : ltens



  return
end subroutine getltens





subroutine getrotens2d(ndime,theta, rotens2d)
  !=======================================================================
  !  getrotens2d= 2d rotation tensor
  !
  !                v'_i = a^t _i ^j v_j
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : dimension
  !
  !  theta : rotation angle between original and new "coordinate"
  !
  !  output:
  !  ------
  !  rotens2d(ndime,ndime) : 2 dimensional rotation tensor for "coordinate"
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime
  real(8), intent(in) :: theta

  real(8), dimension(ndime,ndime), intent(out) :: rotens2d
  ! ====================================

  ! initialize
  rotens2d(:,:)= 0.0d0

  ! error check
  if(ndime /= 2) then
     write(*,*) "this subroutine is only for 2d tensor: getrotens2d"
     write(nout5,*) "this subroutine is only for 2d tensor: getrotens2d"
     stop
  end if

  ! set rotation tensor
  rotens2d(1,1)= dcos(theta)
  rotens2d(1,2)= -dsin(theta) ! cos(90+theta)

  rotens2d(2,1)= dsin(theta) ! cos(90-theta)
  rotens2d(2,2)= dcos(theta)



  return
end subroutine getrotens2d



subroutine getdevtens(nndex,tens, devtens)
  !=======================================================================
  !  getdevtens = compute deviatoric tensor
  !
  !               dev_ij= strs_ij - (1/nndex)*Tr(strs)*del_ij
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nndex : dimension of tensor
  !
  !  tens(nndex,nndex) : tensor
  !
  !  output:
  !  ------
  !  devtens(nndex,nndex) : deviatoric tensor
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nndex
  real(8), dimension(nndex,nndex), intent(in) :: tens

  real(8), dimension(nndex,nndex), intent(out) :: devtens
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: trsum, tr
  real(8) :: krodelt

  integer :: i, j ! loop index
  ! ====================================

  ! initialize
  !devtens(:,:)= 0.0d0


  ! compute trace of tensor: strs_kk= tr(strs)
  trsum= tr(nndex,tens)
    ! input : nndex,tens

  ! compute deviatoric tensor
  do i=1, nndex
     do j=1, nndex

        !devtens(i,j)= tens(i,j) - ( 1.0d0 / 3.0d0 ) * trsum * krodelt(i,j)
        devtens(i,j)= tens(i,j) - onethird * trsum * krodelt(i,j)

     end do
  end do



  return
end subroutine getdevtens





subroutine geteffstrn(nndex,devstrntens, effstrn)
  !=======================================================================
  !  geteffstrn = compute effectice stress
  !
  !               effstrs= sqrt( 2.0d0 / 3.0d0 * (dev : dev) )
  !                      = sqrt( 2.0d0 / 3.0d0 * (dev_ij * dev_ij) )
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nndex : dimension of tensor
  !
  !  devstrntens(nndex,nndex) : deviatoric strain tensor
  !
  !  output:
  !  ------
  !  effstrn : effective strain value
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nndex
  real(8), dimension(nndex,nndex), intent(in) :: devstrntens

  real(8), intent(out) :: effstrn
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: sprdts, devdev
  ! ====================================

  ! initialize
  effstrn= 0.0d0


  ! compute scalar product between deviatoric stress tensor
  devdev= sprdts(0,nndex,devstrntens,devstrntens)
    ! input : 0(opt= :),nndex,devstrntens,devstrntens

  ! compute effetive strain
  effstrn= dsqrt( 2.0d0 / 3.0d0 * devdev )



  return
end subroutine geteffstrn





subroutine geteffstrs(nndex,devstrstens, effstrs)
  !=======================================================================
  !  geteffstrs = compute effectice stress
  !
  !               effstrs= sqrt( 3.0d0 / 2.0d0 * (dev : dev) )
  !                      = sqrt( 3.0d0 / 2.0d0 * (dev_ij * dev_ij) )
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nndex : dimension of tensor
  !
  !  devstrstens(nndex,nndex) : deviatoric stress tensor
  !
  !  output:
  !  ------
  !  effstrs : effective stress value
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nndex
  real(8), dimension(nndex,nndex), intent(in) :: devstrstens

  real(8), intent(out) :: effstrs
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: sprdts, devdev

  ! ====================================

  ! compute scalar product between deviatoric stress tensor
  devdev= sprdts(0,nndex,devstrstens,devstrstens)
    ! input : 0(opt= :),nndex,devstrstens,devstrstens

  ! compute effetive stress
  effstrs= dsqrt( 3.0d0 / 2.0d0 * devdev )



  return
end subroutine geteffstrs





subroutine getj2dirtens(nndex,devstrstens,effstrs, j2dirtens)
  !=======================================================================
  !  getj2dirtens = compute plastic flow direction tensor of j2 material
  !
  !                j2dirtens= (3/2) * (devstrstens/effstrs)
  !
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nndex : dimension of tensor
  !
  !  devstrstens(nndex,nndex) : deviatoric stress tensor
  !
  !  output:
  !  ------
  !  j2dirtens(nndex,nndex) : plastic flow direction tensor of j2 material
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nndex
  real(8), dimension(nndex,nndex), intent(in) :: devstrstens
  real(8), intent(in) :: effstrs

  real(8), dimension(nndex,nndex), intent(out) :: j2dirtens
  ! ====================================
  ! local variable
  ! ==============
  integer :: indx, jndx ! loop index
  ! ====================================

  ! initialize
  j2dirtens(:,:)= 0.0d0


  ! compute plastic flow direction
  if ( effstrs == 0.0d0 ) then
     j2dirtens(:,:)= 0.0d0

  else
     ! j2 dir tensor: ptens= (3/2) * (devstrstens/effstrs)
     do indx=1, nndex
        do jndx=1, nndex
           j2dirtens(indx,jndx)= (3.0d0/2.0d0) * devstrstens(indx,jndx)/effstrs
        end do
     end do
  end if



  return
end subroutine getj2dirtens





subroutine updjgeo(nndex,midltens,delt, strsjtens) 
  !=======================================================================
  !  updjgeo = geometry update for jaumann rate 
  !
  !             strs_n+1= strs_n + (w_n+1/2.strs_n + strs_n.w_n+1/2^t)*delt
  !                                 ---------------------------------
  !                                          geo update
  !
  !                            + (strs_n+1/2 ^jauman)*delt
  !                               ------------------
  !                                 material update
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nndex : dimension of tensor
  !
  !  midltens(nndex,nndex) : velocity gradient tensor at [n+1/2] : predicted value
  !
  !  delt : time increment
  !
  !  inoutput:
  !  --------
  !  strsjtens(nndex,nndex) : stress tensor: geometry update [n] -> [n+1]
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nndex
  real(8), dimension(nndex,nndex), intent(in) :: midltens
  real(8), intent(in) :: delt

  real(8), dimension(nndex,nndex), intent(inout) :: strsjtens
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(nndex,nndex) :: midspintens

  real(8), dimension(nndex,nndex) :: temp1, temp2
  ! ====================================
  ! initialize: do not initialize strsjtens


  ! compute spin tensor: mid step [n+1/2]
  call antisymtens(nndex,midltens, midspintens)
     ! input : nndex,ltens
     ! output : spintens

  ! compute w.kir tensor
  call matprd(nndex,nndex,0, nndex,nndex,0, nndex,nndex, midspintens,strsjtens,temp1)
    ! input : nndex,nndex,0, nndex,nndex,0, nndex,nndex, midspintens,strsjtens
    ! output : temp1

  ! compute kir.w^t tensor
  call matprd(nndex,nndex,0, nndex,nndex,1, nndex,nndex, strsjtens,midspintens,temp2)
    ! input : nndex,nndex,0, nndex,nndex,1, nndex,nndex, strsjtens,midspintens
    ! output : temp2


  ! update rotation of geometry
  strsjtens(1:nndex,1:nndex)= strsjtens(1:nndex,1:nndex) + delt * ( temp1(1:nndex,1:nndex) + temp2(1:nndex,1:nndex) )



  return
end subroutine updjgeo





subroutine updhwgeo(nndex,midltens,delt, strstens) 
  !=======================================================================
  !  updhwgeo = hughes-winget geometry update 
  !
  !             strs_n+1= Q_n+1/2.strs_n.Q_n+1/2^T + (strs_n+1/2 ^jauman)*delt
  !                       ------------------------   -------------------------
  !                           geo update                  material update
  !                                                                                  
  !               note:
  !               ----
  !               hughes and winget, IJNME, 1980, vol.15 , pp. 1862-1867
  !               finite rotation effects in numerical integration
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nndex : dimension of tensor
  !
  !  midltens(nndex,nndex) : velocity gradient tensor at [n+1/2] : predicted value
  !
  !  delt : time increment
  !
  !  inoutput:
  !  --------
  !  strstens(nndex,nndex) : stress tensor: geometry update [n] -> [n+1]
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nndex
  real(8), dimension(nndex,nndex), intent(in) :: midltens
  real(8), intent(in) :: delt

  real(8), dimension(nndex,nndex), intent(inout) :: strstens
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(nndex,nndex) :: midspintens, midwtens, qtens
  real(8) :: krodelt

  real(8), dimension(nndex,nndex) :: invtemp1
  real(8), dimension(nndex,nndex) :: temp1, temp2, temp3, temp4

  ! loop index
  integer :: index, jndex
  ! ====================================

  ! initialize: do not initialize strsjtens


  ! compute spin tensor: mid step [n+1/2]
  call antisymtens(nndex,midltens, midspintens)
     ! input : nndex,midltens
     ! output : midspintens

  ! time integrals over the step of the spin tensor
  midwtens(1:nndex,1:nndex)= delt * midspintens(1:nndex,1:nndex)
 
  ! compute: temp_ij= del_ij - alpha * midwtens_ij: see, eq(13)
  do index=1, nndex
     do jndex=1, nndex
        temp1(index,jndex)= krodelt(index,jndex) - 0.50d0 * midwtens(index,jndex)
     end do
  end do

  ! compute: inv(temp_ij)
  call getinv(nndex,temp1, invtemp1)
     ! input : nndex,temp1
     ! output : invtemp1

  ! compute: inv(temp_ij).midwtens
  call matprd(nndex,nndex,0, nndex,nndex,0, nndex,nndex, invtemp1,midwtens, temp2)
    ! input : nndex,nndex,0, nndex,nndex,0, nndex,nndex, invtemp1,midwtens
    ! output : temp2

  ! compute: Q_ij= del_ij + temp2_ij
  do index=1, nndex
     do jndex=1, nndex
        qtens(index,jndex)= krodelt(index,jndex) + temp2(index,jndex)
     end do
  end do

  ! compute: temp3=Q.strstens
  call matprd(nndex,nndex,0, nndex,nndex,0, nndex,nndex, qtens,strstens, temp3)
    ! input : nndex,nndex,0, nndex,nndex,0, nndex,nndex, qtens,strstens
    ! output : temp3

  ! compute: temp4= temp3.Q^T
  call matprd(nndex,nndex,0, nndex,nndex,1, nndex,nndex, temp3,qtens, temp4)
    ! input : nndex,nndex,0, nndex,nndex,1, nndex,nndex, temp3,qtens
    ! output : temp4

  ! update rotation of geometry
  strstens(1:nndex,1:nndex)= temp4(1:nndex,1:nndex)



  return
end subroutine updhwgeo





subroutine voit3d2tens2d4(voit3d, tens2d)
  !=======================================================================
  !  voit3d2tens2d4 = convert 3d voight form to 2d tensor form
  !                                                                                  
  !                   note:
  !                   ----
  !                   4th order tnesor
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  voit3d(6,6) : 3d voight form of 4th order thensor
  !
  !  output:
  !  ------
  !  tens2d(2,2,2,2) : 2d 4th order tensor
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(6,6), intent(in) :: voit3d

  real(8), dimension(2,2,2,2), intent(out) :: tens2d
  ! ====================================
  ! local variable
  ! ==============

  ! ====================================

  ! initialize
  tens2d(:,:,:,:)= 0.0d0


  ! set major components
  tens2d(1,1,1,1)= voit3d(1,1)
  tens2d(1,1,2,2)= voit3d(1,2) 
  tens2d(2,2,1,1)= voit3d(2,1)
  tens2d(2,2,2,2)= voit3d(2,2)

  tens2d(1,1,1,2)= voit3d(1,6)
  tens2d(2,2,1,2)= voit3d(2,6)
  tens2d(1,2,1,1)= voit3d(6,1)
  tens2d(1,2,2,2)= voit3d(6,2)

  tens2d(1,2,1,2)= voit3d(6,6)

  ! using symmetry condition
  tens2d(1,1,2,1)= tens2d(1,1,1,2)
  tens2d(2,2,2,1)= tens2d(2,2,1,2)
  tens2d(2,1,1,1)= tens2d(1,2,1,1)
  tens2d(2,1,2,2)= tens2d(1,2,2,2)

  tens2d(1,2,2,1)= tens2d(1,2,1,2)
  tens2d(2,1,1,2)= tens2d(1,2,1,2)
  tens2d(2,1,2,1)= tens2d(1,2,1,2)



  return
end subroutine voit3d2tens2d4





subroutine getcstrtens3d4(sigtens3d, cstrtens3d)
  !=======================================================================
  !  getcstrtens3d4 = compute 4th order 3d c^* tensor
  !
  !                   note:
  !                   ----
  !                   see, TB's book page 231, eq. (5.4.28) and page 232, box 5.1
  !
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  sigtens3d(3,3) : 3d cauchy stress tensor
  !
  !  output:
  !  ------
  !  cstrtens3d(3,3,3,3) : 3d 4th order c^* tensor
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(3,3), intent(in) :: sigtens3d

  real(8), dimension(3,3,3,3), intent(out) :: cstrtens3d
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: term1, term2, term3, term4, term5
  real(8) :: krodelt

  ! loop index
  integer :: i, j, k, l
  ! ====================================

  ! initialize
  cstrtens3d(:,:,:,:)= 0.0d0


  do i=1, 3
     do j=1, 3
        do k=1, 3
           do l=1, 3

              term1= krodelt(i,k) * sigtens3d(j,l)
              term2= krodelt(i,l) * sigtens3d(j,k)
              term3= krodelt(j,k) * sigtens3d(i,l)
              term4= krodelt(j,l) * sigtens3d(i,k)
              term5= sigtens3d(i,j) * krodelt(k,l)

              cstrtens3d(i,j,k,l)= 0.50d0 * ( term1 + term2 + term3 + term4 ) - term5

           end do
        end do
     end do
  end do



  return
end subroutine getcstrtens3d4





subroutine voit2ind3d4(voigt3d4, tensor3d4)
  !=======================================================================
  !  voit2ind= convert 3d 4th order voit form matrix to indicial form tensor
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  voigt3d4(6,6) : voigt format matrix
  !
  !  output:
  !  ------
  !  tensor3d4(3,3,3,3) : inditial format tensor
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(6,6), intent(in) :: voigt3d4

  real(8), dimension(3,3,3,3), intent(out) :: tensor3d4
  ! ====================================
  ! local variable
  ! ==============
  integer, dimension(3,3) :: ij2a
  integer :: a, b

  ! loop index
  integer :: i, j, k, l
  ! ====================================

  ! initialize
  tensor3d4(:,:,:,:)= 0.0d0


  ! indicial to voight
  data ij2a /1,6,5, 6,2,4, 5,4,3/

  ! convert voight form to tensor form
  do i=1, 3
     do j=1, 3
        do k=1,3
           do l=1, 3

              a= ij2a(i,j)
              b= ij2a(k,l)

              tensor3d4(i,j,k,l)= voigt3d4(a,b)

           end do
        end do
     end do
  end do



  return
end subroutine voit2ind3d4





subroutine ctantenscj2ct3d(ctantenscj3d,sigtens3d, ctantensct3d)
  !=======================================================================
  !  ctantenscj2ct3d= convert 
  !
  !                   note:
  !                   ----
  !                   tb's book page 232, box 5.1
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ctantenscj3d(3,3,3,3) : tangent moduli of cauchy stress jauman rate
  !
  !  sigtens3d(3,3) : cauchy stress tensor
  !
  !  output:
  !  ------
  !  ctantensct3d(3,3,3,3) : tangent moduli of cauchy stress truesdell rate
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(3,3,3,3), intent(in) :: ctantenscj3d
  real(8), dimension(3,3), intent(in) :: sigtens3d

  real(8), dimension(3,3,3,3), intent(out) ::ctantensct3d
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(3,3,3,3) :: cstrtens3d

  ! ====================================

  ! initialize
  ctantensct3d(:,:,:,:)= 0.0d0


  ! compute c^* tensor: tb's book page 231, eq.(5.4.28)
  call getcstrtens3d4(sigtens3d, cstrtens3d)
     ! input : sigtens3d
     ! output : cstrtens3d

  ! convert tangent moduli of caychy stress truesdell rate
  ctantensct3d(1:3,1:3,1:3,1:3)= ctantenscj3d(1:3,1:3,1:3,1:3) - cstrtens3d(1:3,1:3,1:3,1:3)



  return
end subroutine ctantenscj2ct3d





subroutine getatens3d4(ctantens3d,sigtens3d, atens3d)
  !=======================================================================
  !  getatens3d4= compute perturbed tangent moduli tensor
  !
  !               note:
  !               ----
  !               for the checking material stability in the current configuration
  !               tb's book, page 387, eq. (6.7.18)
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ctantens3d(3,3,3,3) : tangent moduli of cauchy stress truesedell rate
  !
  !  sigtens3d(3,3) : 3d cauchy stress tensor
  !
  !  output:
  !  ------
  !  atens3d(3,3,3,3) : perturbed tangent moduli tensor
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(3,3,3,3), intent(in) :: ctantens3d
  real(8), dimension(3,3), intent(in) :: sigtens3d

  real(8), dimension(3,3,3,3), intent(out) :: atens3d
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: krodelt

  ! loop index
  integer :: i, j, k, l
  ! ====================================

  ! initialize
  atens3d(:,:,:,:)= 0.0d0


  ! compute perturbed tangnet stiffness: a_ijkl= c_ijkl + sig_il*kdelt_jk
  do i=1, 3
     do j=1, 3
        do k=1, 3
           do l=1, 3

              atens3d(i,j,k,l)= ctantens3d(i,j,k,l) + sigtens3d(i,l) * krodelt(j,k)

           end do
        end do
     end do
  end do



  return
end subroutine getatens3d4





subroutine curn2inivec(ndime,ftens,veccurn, vecini)
  !=======================================================================
  !  curn2inivec= convert a vector in the current domain to initial dommain
  !
  !               note:
  !               ----
  !               
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !
  !  output:
  !  ------
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime
  real(8), dimension(ndime,ndime), intent(in) :: ftens
  real(8), dimension(ndime,1), intent(in) :: veccurn

  real(8), dimension(ndime,1), intent(out) :: vecini
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: getdet
  real(8), dimension(ndime,ndime) :: invftens

  ! loop index
  integer :: i, j, k, l
  ! ====================================

  ! initialize
  vecini(:,:)= 0.0d0


  ! check determinant of f tensor
  if( getdet(ndime,ftens) == 0.0d0 ) then ! no deformation

     vecini(1:ndime,1)= veccurn(1:ndime,1)

  else if ( getdet(ndime,ftens) > 0.0d0 ) then

     ! compute inverse of f tensor: inv(f)=a X_j / a x_i
     call getinv(ndime,ftens, invftens)
       ! input : ndime,ftens
       ! output : invftens

     call matprd(ndime,ndime,0, ndime,1,0, ndime,1, invftens,veccurn, vecini)
       ! input : ndime,ndime,0, ndime,1,0, ndime,1, invftens,veccurn
       ! output : vecini

  else
     write(*,*) "determinant of f tensor is negative: curn2inivec"
     write(nout5,*) "determinant of f tensor is negative: curn2inivec"
     stop

  end if



  return
end subroutine curn2inivec





subroutine curn2iniang2d(ftens,angcurn, angini)
  !=======================================================================
  !  curn2iniang2d= convert a direction angle in the current domain to initial dommain
  !
  !               note:
  !               ----
  !               
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !
  !  output:
  !  ------
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(2,2), intent(in) :: ftens
  real(8), intent(in) :: angcurn

  real(8), intent(out) :: angini
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(2,1) :: veccurn
  real(8), dimension(2,1) :: vecini

  ! loop index
  integer :: i, j, k, l
  ! ====================================

  ! initialize
  angini= 0.0d0


  ! set current direction vector
  veccurn(1,1)= dcos(angcurn)
  veccurn(2,1)= dsin(angcurn)

  ! convert the direction vector in current domain to initial domain
  call curn2inivec(2,ftens,veccurn, vecini)
     ! input : 2(ndime),ftens,veccurn
     ! output : vecini

  ! maximum principal direction
  angini= datan2( vecini(2,1),vecini(1,1) )

  ! angini== +/- 180
  if ( angini >= 180.0d0*d2r .or. angini <= -180.0d0*d2r ) then
      angini= 0.0d0

  ! 90 < angini < 180
  else if ( 90.0d0*d2r < angini .and. angini < 180.0d0*d2r ) then
      angini= angini - 180.0d0 * d2r

  ! -180 < angini < -90
  else if ( -180.0d0*d2r < angini .and. angini< -90.0d0*d2r ) then
      angini= angini + 180.0d0 * d2r

  end if



  return
end subroutine curn2iniang2d




  

real(8) function trsprs(nsprs,krow,kcol,sprsmat)
  !=======================================================================
  !  trsprs = compute trace of given sparse matrix
  !
  !           tr= Tr(a) = a_kk
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nsprs : total number of sparse components
  !
  !  krow(*), kcol(*) : row and column index
  !
  !  sprsmat(*) : sparse matrix
  !
  !  output:
  !  ------
  !  tr : trace of tensor
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nsprs
  integer, dimension(*), intent(in) :: krow, kcol
  real(8), dimension(*), intent(in) :: sprsmat
  ! ====================================
  ! local variable
  ! ==============
  
  ! loop index
  integer :: isprs
  ! ====================================

  ! initialize
  trsprs= 0.0d0


  ! compute trace of given tensor
  do isprs=1, nsprs

     if ( krow(isprs) == kcol(isprs) ) then

        trsprs= trsprs + sprsmat(isprs)

     end if

  end do



  return
end function trsprs





subroutine getbtens(ndime,ftens, btens)
  !=======================================================================
  !  getbtens= compute left cauchy-green deformation tensor
  !
  !            b = f.f^t
  !
  !            note:
  !            ----
  !            rotation-independent deformation measures
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : dimension of tensor
  !
  !  ftens(ndime,ndime) : deformation gradient tensor
  !  
  !  output:
  !  ------
  !  btens(ndime,ndime) : left cauchy-green deformation tensor
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime
  real(8), dimension(ndime,ndime), intent(in) :: ftens

  real(8), dimension(ndime,ndime), intent(out) :: btens
  ! ====================================
  ! local variable
  ! ==============
  
  ! ====================================

  ! initialize
  btens(:,:)= 0.0d0


  ! compute left cauchy-green deformation tensor: b = f.f^t
  call matprd(ndime,ndime,0, ndime,ndime,1, ndime,ndime, ftens,ftens, btens)
     ! input : ndime,ndime,0, ndime,ndime,1, ndime,ndime, ftens,ftens
     ! output : btens



  return
end subroutine getbtens





subroutine getbbartens(ndime,ftens, bbartens)
  !=======================================================================
  !  getbtens= compute volume preserving part of left cauchy-green deformation tensor
  !
  !            b = fbar.fbar^t
  !
  !            fbar= |f|^(-1/3) f
  !
  !            note:
  !            ----
  !            rotation-independent deformation measures
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : dimension of tensor
  !
  !  ftens(ndime,ndime) : deformation gradient tensor
  !  
  !  output:
  !  ------
  !  bbartens(ndime,ndime) : volume preserving part of left cauchy-green deformation tensor
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime
  real(8), dimension(ndime,ndime), intent(in) :: ftens

  real(8), dimension(ndime,ndime), intent(out) :: bbartens
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: jdet, getdet
  real(8), dimension(ndime,ndime) :: fbartens
  
  ! ====================================

  ! initialize
  bbartens(:,:)= 0.0d0


  ! compute determinant of ftens tensor: jdet
  jdet= getdet(ndime,ftens)

  ! volume-preserving part of deformation gradient tensor
  fbartens(1:ndime,1:ndime)= jdet**( -1.0d0 / 3.0d0 ) * ftens(1:ndime,1:ndime)

  ! compute volume-preserving part of left cauchy-green deformation tensor: bbar = fbar.fbar^t
  call matprd(ndime,ndime,0, ndime,ndime,1, ndime,ndime, fbartens,fbartens, bbartens)
     ! input : ndime,ndime,0, ndime,ndime,1, ndime,ndime, fbartens,fbartens
     ! output : bbartens



  return
end subroutine getbbartens





subroutine getctens(ndime,ftens, ctens)
  !=======================================================================
  !  getbtens= compute right cauchy-green deformation tensor
  !
  !            c = f^t.f
  !
  !            note:
  !            ----
  !            rotation-independent deformation measures
  !          
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : dimension of tensor
  !
  !  ftens(ndime,ndime) : deformation gradient tensor
  !  
  !  output:
  !  ------
  !  ctens(ndime,ndime) : right cauchy-green deformation tensor
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime
  real(8), dimension(ndime,ndime), intent(in) :: ftens

  real(8), dimension(ndime,ndime), intent(out) :: ctens
  ! ====================================
  ! local variable
  ! ==============
  
  ! ====================================

  ! initialize
  ctens(:,:)= 0.0d0


  ! compute right cauchy-green deformation tensor: b = f^t.f
  call matprd(ndime,ndime,1, ndime,ndime,0, ndime,ndime, ftens,ftens, ctens)
     ! input : ndime,ndime,1, ndime,ndime,0, ndime,ndime, ftens,ftens
     ! output : ctens



  return
end subroutine getctens
