! ==================================
! interpret reading input data 
! ==================================
!      type                  name                              argument
!      ----                  ----                              --------
! 4.  subroutine          getgqsize                 (optmhd,optele,ngqpt, mgaus,mgqpt)
!
! =========================================================================================================



subroutine getgqsize(optmhd,optele,ngqpt, mgaus,mgqpt)
  !=======================================================================
  !  getgqsize = get general parameter
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optmhd : numerical method option
  !
  !  optele : element type option
  !
  !  ngqpt(3) : gq rule
  !             ngqpt(1) : gq rule for regular element 
  !             ngqpt(2) : gq rule for enriched element 
  !             ngqpt(3) : gq rule through thickness
  !
  !  output:
  !  ------
  !  mgaus(3) : the number of quadrature point
  !             mgaus(1)= regular element
  !             mgaus(2)= enriched element
  !             mgaus(3)= through thickness
  !
  !  mgqpt(2) : the total number of quadrature point
  !             mgqpt(1)= mgaus(1) * mgaus(3) : regular element
  !             mgqpt(2)= mgaus(2) * mgaus(3) : enriched element
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optmhd, optele
  integer, dimension(3), intent(in) :: ngqpt
  
  integer, dimension(3), intent(out) :: mgaus
  integer, dimension(2), intent(out) :: mgqpt
  ! ====================================

  ! initialize
  mgaus(:)= 0
  mgqpt(:)= 0


  ! --------------------------------------------------------------------
  ! set mgaus(1:3)
  ! --------------

  ! in plane regular element: ngqpt(1)
  ! ------------------------
  if( optele==0 ) then ! 1d line
     mgaus(1)= ngqpt(1)

  else if( optele==1 ) then ! 2d tri
     mgaus(1)= ngqpt(1)

  else if( optele==2 ) then ! 2d quad
     mgaus(1)= ngqpt(1)**2
  
  else if( optele==3 ) then ! 3d tb shell
     mgaus(1)= ngqpt(1)**2

  else if( optele==4 ) then ! 3d brick
     mgaus(1)= ngqpt(1)**3

  end if
 
  ! in plane enriched element: ngqpt(2)
  ! -------------------------
  select case(optmhd)
  case(0:1) ! conventional fem or element delection
     mgaus(2)= 0

  case(2) ! discrete phantom nodes crack method
     mgaus(2)= 1

  case(3:4) ! phantom nodes method
     if( optele==1 ) then ! 2d tri
        mgaus(2)= ngqpt(1)

     else if( optele==2 ) then ! 2d quad
        mgaus(2)= ngqpt(1)**2

     else if( optele==3 ) then ! 3d tb shell
        mgaus(2)= ngqpt(1)**2

     else if( optele==4 ) then ! 3d brick
        mgaus(2)= ngqpt(1)**3

     end if

  case(5:6) ! xfem method
     if( optele==1 ) then ! 2d tri
        mgaus(2)= ngqpt(2)

     else if( optele==2 ) then ! 2d quad
        mgaus(2)= ngqpt(2)**2

     else if( optele==3 ) then ! 3d tb shell
        mgaus(2)= ngqpt(2)**2

     else if( optele==4 ) then ! 3d brick
        mgaus(2)= ngqpt(2)**3

     end if

  end select

  ! through thickness: ngqpt(3)
  ! -----------------
  if( optele==0 ) then ! 1d line
     mgaus(3)= 1

  else if( optele==1 ) then ! 2d tri
     mgaus(3)= 1

  else if( optele==2 ) then ! 2d quad
     mgaus(3)= 1

  else if( optele==3 ) then ! 3d tb shell
     mgaus(3)= ngqpt(3)

  else if( optele==4 ) then ! 3d brick
     mgaus(3)= 1

  end if

  ! --------------------------------------------------------------------
  ! set mgqpt(1:2)
  ! --------------

  ! total gauss point: total= inplane
  ! -----------------
  mgqpt(1)= mgaus(1) * mgaus(3) ! regular element
  mgqpt(2)= mgaus(2) * mgaus(3) ! enriched element

  ! --------------------------------------------------------------------



  return
end subroutine getgqsize

