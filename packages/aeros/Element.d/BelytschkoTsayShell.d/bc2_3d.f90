! ==========================
! apply boundary condition2 : element surface
! ==========================
!      type                  name                              argument
!      ----                  ----                              --------
! 1.  subroutine          elefbc3dbrkshl2      (opttrc,ecord3d,edisp3d,trac, efbc)
! 2.  subroutine          ele3dsurf2d4nod0     (ecord3d,edisp3d, area,nvec)
!
! =========================================================================================================



subroutine elefbc3dbrkshl2(opttrc,ecord3d,edisp3d,trac, efbc)
  !=======================================================================
  !  elefbc3dbrkshl2 = compute 3d bt shell or birck element force due to
  !                    element surface stress or pressure
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  opttrc : traction type option handler
  !           opttrc= 0 : pressure
  !           opttrc= 1 : traction
  !
  !  output:
  !  ------
  !  efbc(12,1) : element nodal force vector ( translational dof only)
  !                          
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: opttrc
  real(8), dimension(3,4), intent(in) :: ecord3d
  real(8), dimension(3,4), intent(in) :: edisp3d
  real(8), dimension(3,1), intent(in) :: trac

  real(8), dimension(12,1), intent(out) :: efbc
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: area
  real(8), dimension(3,1) :: nvec
  real(8), dimension(3,1) :: tracvec

  ! ====================================

  ! initialize
  efbc(:,:)= 0.0d0


  ! ---------------------------------------------------------------------
  ! compute surface area and normal vector
  ! --------------------------------------
  call ele3dsurf2d4nod0(ecord3d,edisp3d, area,nvec)
     ! input : ecord3d,edisp3d
     ! output : area,nvec

  ! ---------------------------------------------------------------------
  ! get global traction vector
  ! --------------------------
  select case(opttrc)
  case(0) ! pressure
     ! set pressure vector
     tracvec(1:3,1)= trac(1,1) * nvec(1:3,1)

  case(1) ! traction
     tracvec(1:3,1)= trac(1:3,1)

  case default
     write(*,*) "wrong opttrc: elefbc3dbrkshl2"
     write(nout5,*) "wrong opttrc: elefbc3dbrkshl2"
     stop

  end select

  ! ---------------------------------------------------------------------
  ! set element nodal force
  ! -----------------------
  efbc(1:3,1)= area *  tracvec(1:3,1) / 4.0d0    ! f_x, f_y and f_z at node 1
  efbc(4:6,1)= area *  tracvec(1:3,1) / 4.0d0    ! f_x, f_y and f_z at node 2
  efbc(7:9,1)= area *  tracvec(1:3,1) / 4.0d0    ! f_x, f_y and f_z at node 3
  efbc(10:12,1)= area *  tracvec(1:3,1) / 4.0d0  ! f_x, f_y and f_z at node 4
  ! ---------------------------------------------------------------------



  return
end subroutine elefbc3dbrkshl2





subroutine ele3dsurf2d4nod0(ecord3d,edisp3d, area,nvec)
  !=======================================================================
  !  ele3dsurf2d4nod0 = compute 3d bt shell or birck element
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
  real(8), dimension(3,4), intent(in) :: ecord3d
  real(8), dimension(3,4), intent(in) :: edisp3d

  real(8), intent(out) :: area
  real(8), dimension(3,1), intent(out) :: nvec
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(3,4) :: ecurn
  real(8), dimension(3,3) :: locbvec
  real(8), dimension(3,4) :: ecurnloc

  real(8), dimension(3,1) :: norlvecloc

  ! loop index
  ! ====================================

  ! initialize
  area= 0.0d0
  nvec(:,:)= 0.0d0


  ! -------------------------------------------------------------------------------
  ! construct corotational coordinate
  ! ---------------------------------
  ! set element current coordinate
  ecurn(1:3,1:4)= ecord3d(1:3,1:4) + edisp3d(1:3,1:4)

  ! compute co rotational local base vector: locbvec
  call getlocbvecbt(ecurn, locbvec)
     ! input : ecurn
     ! output : locbvec

  ! -------------------------------------------------------------------------------
  ! compute surface area
  ! --------------------
  ! get local nodal coordinate: glb -> loc
  call glb2locnodv(3,4,3,locbvec,ecurn, ecurnloc)
     ! input : 3(ndime),4(nnode),3(ntrndof),locbvec,ecurn
     ! output : ecurnloc

  ! compute area
  area= 0.50d0*( (ecurnloc(1,3)-ecurnloc(1,1))*(ecurnloc(2,4)-ecurnloc(2,2)) &
                +(ecurnloc(1,2)-ecurnloc(1,4))*(ecurnloc(2,3)-ecurnloc(2,1)) )

  ! -------------------------------------------------------------------------------
  ! compute surface unit normal vector
  ! ---------------------------------
  ! set normal to corotational element surface
  norlvecloc(1,1)= 0.0d0  ! x^
  norlvecloc(2,1)= 0.0d0  ! y^
  norlvecloc(3,1)= -1.0d0 ! z^

  ! convert local to global
  call loc2glbv(3,3,locbvec,norlvecloc, nvec)
     ! input : 3(ndime),3(ntrndof),locbvec,norlvecloc
     ! output : nvec



  return
end subroutine ele3dsurf2d4nod0


subroutine elefbc3dtri(opttrc,optele,ecord3d,edisp3d,trac, efbc)
  !=======================================================================
  !  elefbc3dtri = compute 3d tri shell element force due to
  !                element surface stress or pressure
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  opttrc : traction type option handler
  !           opttrc= 0 : pressure
  !           opttrc= 1 : traction
  !
  !  output:
  !  ------
  !  efbc(9,1) : element nodal force vector ( translational dof only)
  !                          
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: opttrc, optele
  real(8), dimension(3,3), intent(in) :: ecord3d
  real(8), dimension(3,3), intent(in) :: edisp3d
  real(8), dimension(3,1), intent(in) :: trac

  real(8), dimension(9,1), intent(out) :: efbc
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: area
  real(8), dimension(3,1) :: nvec
  real(8), dimension(3,1) :: tracvec

  ! ====================================

  ! initialize
  efbc(:,:)= 0.0d0


  ! ---------------------------------------------------------------------
  ! compute surface area and normal vector
  ! --------------------------------------
  call ele3dsurf2d3nod0(ecord3d,edisp3d, area,nvec)
     ! input : ecord3d,edisp3d
     ! output : area,nvec

  ! ---------------------------------------------------------------------
  ! get global traction vector
  ! --------------------------
  select case(opttrc)
  case(0) ! pressure
     ! set pressure vector
     tracvec(1:3,1)= trac(1,1) * nvec(1:3,1)

  case(1) ! traction
     tracvec(1:3,1)= trac(1:3,1)

  case default
     write(*,*) "wrong opttrc: elefbc3dbrkshl2"
     write(nout5,*) "wrong opttrc: elefbc3dbrkshl2"
     stop

  end select


  ! ---------------------------------------------------------------------
  ! set element nodal force
  ! -----------------------
  efbc(1:3,1)= area *  tracvec(1:3,1) / 3.0d0    ! f_x, f_y and f_z at node 1
  efbc(4:6,1)= area *  tracvec(1:3,1) / 3.0d0    ! f_x, f_y and f_z at node 2
  efbc(7:9,1)= area *  tracvec(1:3,1) / 3.0d0    ! f_x, f_y and f_z at node 3
  ! ---------------------------------------------------------------------



  return
end subroutine elefbc3dtri





subroutine ele3dsurf2d3nod0(ecord3d,edisp3d, area,nvec)
  !=======================================================================
  !  ele3dsurf2d3nod0 = compute 3d tri shell 
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
  real(8), dimension(3,3), intent(in) :: ecord3d
  real(8), dimension(3,3), intent(in) :: edisp3d

  real(8), intent(out) :: area
  real(8), dimension(3,1), intent(out) :: nvec
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(3,3) :: ecurn
  real(8), dimension(3,3) :: locbvec
  real(8), dimension(3,3) :: ecurnloc

  real(8), dimension(3,1) :: norlvecloc

  ! loop index
  ! ====================================

  ! initialize
  area= 0.0d0
  nvec(:,:)= 0.0d0


  ! -------------------------------------------------------------------------------
  ! construct corotational coordinate
  ! ---------------------------------
  ! set element current coordinate
  ecurn(1:3,1:3)= ecord3d(1:3,1:3) + edisp3d(1:3,1:3)

  ! compute co rotational local base vector: locbvec
  call getlocbvectri(ecurn, locbvec)
     ! input : ecurn
     ! output : locbvec

  ! -------------------------------------------------------------------------------
  ! compute surface area
  ! --------------------
  ! get local nodal coordinate: glb -> loc
  call glb2locnodv(3,3,3,locbvec,ecurn, ecurnloc)
     ! input : 3(ndime),3(nnode),3(ntrndof),locbvec,ecurn
     ! output : ecurnloc

  ! compute area
  area= 0.50d0*( (ecurnloc(1,2)-ecurnloc(1,1))*(ecurnloc(2,3)-ecurnloc(2,1)) &
                +(ecurnloc(1,1)-ecurnloc(1,3))*(ecurnloc(2,2)-ecurnloc(2,1)) )

  ! -------------------------------------------------------------------------------
  ! compute surface unit normal vector
  ! ---------------------------------
  ! set normal to corotational element surface
  norlvecloc(1,1)= 0.0d0  ! x^
  norlvecloc(2,1)= 0.0d0  ! y^
  norlvecloc(3,1)= -1.0d0 ! z^

  ! convert local to global
  call loc2glbv(3,3,locbvec,norlvecloc, nvec)
     ! input : 3(ndime),3(ntrndof),locbvec,norlvecloc
     ! output : nvec



  return
end subroutine ele3dsurf2d3nod0


subroutine elefbc3dbrkshl2opt(area,trac,tmftval, efint)
  !=======================================================================
  !  elefbc3dbrkshl2 = subtract 3d bt shell or brick element force due to
  !                    element surface stress or pressure
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !
  !  output:
  !  ------
  !  efbc(24,1) : local element nodal force vector
  !                          
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), intent(in) :: area
  real(8), dimension(3,1), intent(in) :: trac
  real(8), intent(in) :: tmftval

  real(8), dimension(24,1), intent(inout) :: efint
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(3,1) :: efbc

  ! ====================================

  ! ---------------------------------------------------------------------
  ! update element nodal force
  ! -----------------------
  efbc(1:3,1)= tmftval * area * trac(1:3,1) / 4.0d0 ! f_x, f_y and f_z at each node

  efint(1:3,1)= efint(1:3,1) - efbc(1:3,1)
  efint(7:9,1)= efint(7:9,1) - efbc(1:3,1)
  efint(13:15,1)= efint(13:15,1) - efbc(1:3,1)
  efint(19:21,1)= efint(19:21,1) - efbc(1:3,1)
  ! ---------------------------------------------------------------------

  return
end subroutine elefbc3dbrkshl2opt
