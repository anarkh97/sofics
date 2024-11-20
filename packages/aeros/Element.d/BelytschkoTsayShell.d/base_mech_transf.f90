! ==================================
! coordinate transformation
! ==================================
!      type                  name                              argument
!      ----                  ----                              --------
! 1.  subroutine             gettransmat                (ndime,ntrndof,locbvec, transmat)
! 2.  subroutine             loc2glbtens                (ndime,ntrndof,locbvec,loctens, glbtens)
! 3.  subroutine             glb2loctens                (ndime,ntrndof,locbvec,glbtens, loctens)
! 4.  subroutine             glb2locnodv                (ndime,nnode,ntrndof,locbvec,eglbnodv, elocnodv)
! 5.  subroutine             loc2glbnodv                (ndime,nnode,ntrndof,locbvec,elocnodv, eglbnodv)
! 6.  subroutine             glb2locv                   (ndime,ntrndof,locbvec,glbvec, locvec)
! 7.  subroutine             loc2glbv                   (ndime,ntrndof,locbvec,locvec, glbvec)
! 8.  subroutine             loc2loc0v                  (ndime,ntrndof,locbvec,locbvec0,locvec, locvec0)
!
! bt shell
! --------
! 9.  subroutine             getlocbvecbt               (ecurn, locbvec)
!
! =========================================================================================================



subroutine gettransmat(ndime,ntrndof,locbvec, transmat)
  !=======================================================================
  !  gettransmat = get transformation matrix
  !                transmat is the matrix of direction cosine
  !                between local(x,y,z) and global(1,2,3)
  !
  !                transformation matrix:
  !                ---------------------
  !                a= [ e_x1, e_x2, e_x3 ]
  !                   [ e_y1, e_y2, e_y3 ]
  !                   [ e_z1, e_z2, e_z3 ]
  !
  !                note:
  !                ----
  !                loc_vec = a   * glb_vec
  !                glb_vec = a^t * loc_vec
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : dimension of base vector
  !
  !  ntrndof : number of dof to transform
  !
  !  locbvec(ndime,ndime) : local base vector
  !                         locbvec(i,j)= [ e_x1, e_y1, e_z1 ]
  !                                       [ e_x2, e_y2, e_z2 ]
  !                                       [ e_x3, e_y3, e_z3 ]
  !                     
  !  output:
  !  ------
  !  transmat(ndime,1) : local vector
  !                            
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime, ntrndof
  real(8), dimension(ndime,ndime), intent(in) :: locbvec

  real(8), dimension(ntrndof,ntrndof), intent(out) :: transmat
  ! ====================================
  ! local variable
  ! ==============

  ! loop index
  integer :: idime, jdime
  ! ====================================

  ! initialize
  transmat(:,:)= 0.0d0


  if ( ntrndof == ndime ) then

     do idime=1, ndime
        do jdime=1, ndime
           transmat(idime,jdime)= locbvec(jdime,idime)
        end do
     end do

  else if ( ndime == 3 .and. ntrndof == 5 ) then
     !    a= [ e_x1, e_x2, e_x3,    0,    0 ]
     !       [ e_y1, e_y2, e_y3,    0,    0 ]
     !       [ e_z1, e_z2, e_z3,    0,    0 ]
     !       [    0,    0,    0, e_x1, e_x2 ]
     !       [    0,    0,    0, e_y1, e_y2 ]

     do idime=1, 3
        do jdime=1, 3
           transmat(idime,jdime)= locbvec(jdime,idime)
        end do
     end do

     transmat(4:5,4:5)= transmat(1:2,1:2)

  else if ( ndime == 3 .and. ntrndof == 6 ) then
     !    a= [ e_x1, e_x2, e_x3,    0,    0,    0 ]
     !       [ e_y1, e_y2, e_y3,    0,    0,    0 ]
     !       [ e_z1, e_z2, e_z3,    0,    0,    0 ]
     !       [    0,    0,    0, e_x1, e_x2, e_x3 ]
     !       [    0,    0,    0, e_y1, e_y2, e_y3 ]
     !       [    0,    0,    0, e_z1, e_z2, e_z3 ]

     do idime=1, 3
        do jdime=1, 3
           transmat(idime,jdime)= locbvec(jdime,idime)
        end do
     end do

     transmat(4:6,4:6)= transmat(1:3,1:3)

   else
     write(*,*) "not available yet: glb2locv"
     write(nout5,*) "not available yet: glb2locv"
     stop

   end if



  return
end subroutine gettransmat






subroutine loc2glbtens(ndime,ntrndof,locbvec,loctens, glbtens)
  !=======================================================================
  !  loc2glbtens = transform local tensor to global tensor
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : dimension of base vector
  !
  !  ntrndof : number of dof to transform
  !
  !  locbvec(ndime,ndime) : local base vector
  !
  !  loctens(ntrndof,ntrndof) : local tensor
  !
  !  output:
  !  ------
  !  glbtens(ntrndof,ntrndof) : global tensor
  !                            
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime, ntrndof
  real(8), dimension(ndime,ndime), intent(in) :: locbvec
  real(8), dimension(ntrndof,ntrndof), intent(in) :: loctens

  real(8), dimension(ntrndof,ntrndof), intent(out) :: glbtens
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(ntrndof,ntrndof) :: transmat
  real(8), dimension(ntrndof,ntrndof) :: tmptens

  ! loop index
  integer :: idime, jdime
  ! ====================================

  ! initialize
  glbtens(:,:)= 0.0d0

  ! get transformation matrix
  call gettransmat(ndime,ntrndof,locbvec, transmat)
     ! input : ndime,ntrndof,locbvec
     ! output : transmat

  ! transform vector: local -> global : glb_tens = a^t * loc_tens * a
  ! ----------------
  call matprd(ntrndof,ntrndof,1, ntrndof,ntrndof,0, ntrndof,ntrndof, transmat,loctens, tmptens)
     ! input : ntrndof,ntrndof,1, ntrndof,ntrndof,0, ntrndof,ntrndof, transmat,loctens
     ! output : tmptens

  call matprd(ntrndof,ntrndof,0, ntrndof,ntrndof,0, ntrndof,ntrndof, tmptens,transmat, glbtens)
     ! input : ntrndof,ntrndof,0, ntrndof,ntrndof,0, ntrndof,ntrndof, tmptens,transmat
     ! output : glbtens



  return
end subroutine loc2glbtens





subroutine glb2loctens(ndime,ntrndof,locbvec,glbtens, loctens)
  !=======================================================================
  !  glb2loctens = transform global tensor to local tensor
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : dimension of base vector
  !
  !  ntrndof : number of dof to transform
  !
  !  locbvec(ndime,ndime) : local base vector
  !
  !  glbtens(ntrndof,ntrndof) : global tensor
  !
  !  output:
  !  ------
  !  loctens(ntrndof,ntrndof) : local tensor
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime, ntrndof
  real(8), dimension(ndime,ndime), intent(in) :: locbvec
  real(8), dimension(ntrndof,ntrndof), intent(in) :: glbtens

  real(8), dimension(ntrndof,ntrndof), intent(out) :: loctens
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(ntrndof,ntrndof) :: transmat
  real(8), dimension(ntrndof,ntrndof) :: tmptens

  ! loop index
  integer :: idime, jdime
  ! ====================================

  ! initialize
  loctens(:,:)= 0.0d0

  ! get transformation matrix
  call gettransmat(ndime,ntrndof,locbvec, transmat)
     ! input : ndime,ntrndof,locbvec
     ! output : transmat

  ! transform vector: global -> local : loc_tens = a * glb_tens * a^t
  ! ----------------
  call matprd(ntrndof,ntrndof,0, ntrndof,ntrndof,0, ntrndof,ntrndof, transmat,glbtens, tmptens)
     ! input : ntrndof,ntrndof,0, ntrndof,ntrndof,0, ntrndof,ntrndof, transmat,glbtens
     ! output : tmptens

  call matprd(ntrndof,ntrndof,0, ntrndof,ntrndof,1, ntrndof,ntrndof, tmptens,transmat, loctens)
     ! input : ntrndof,ntrndof,0, ntrndof,ntrndof,1, ntrndof,ntrndof, tmptens,transmat
     ! output : loctens



  return
end subroutine glb2loctens





subroutine glb2locnodv(ndime,nnode,ntrndof,locbvec,eglbnodv, elocnodv)
  !=======================================================================
  !  glb2locnod= convert global nodal vector to local nodal vector
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  npoin : 
  !  
  !
  !  output:
  !  ------
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime, nnode, ntrndof
  real(8), dimension(ndime,ndime), intent(in) :: locbvec
  real(8), dimension(ntrndof,nnode), intent(in) :: eglbnodv

  real(8), dimension(ntrndof,nnode), intent(out) :: elocnodv
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(ntrndof,1) :: glbvec, locvec
  real(8), dimension(ntrndof,ntrndof) :: transmat !pjsa

  ! loop index
  integer :: inode
  ! ====================================

  ! initialize
  !elocnodv(:,:)= 0.0d0
  call gettransmat(ndime,ntrndof,locbvec, transmat) !pjsa

  ! convert global to local
  ! -----------------------
  do inode=1, nnode

     ! get global nodal vector
     glbvec(1:ntrndof,1)= eglbnodv(1:ntrndof,inode)

     ! convert global to local
     !call glb2locv(ndime,ntrndof,locbvec,glbvec, locvec)
     call dgemv('n',ntrndof,ntrndof,1.0d0,transmat,ntrndof,glbvec,1,0.0d0,locvec,1) !pjsa
        ! input : ndime,ntrndof,locbvec,glbvec
        ! output : locvec

     ! set local nodal vector
     elocnodv(1:ntrndof,inode)= locvec(1:ntrndof,1)     

  end do


  return
end subroutine glb2locnodv



subroutine loc2glbnodv(ndime,nnode,ntrndof,locbvec,elocnodv, eglbnodv)
  !=======================================================================
  !  loc2glbnodv= convert local nodal vector to global nodal vector
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  npoin : 
  !  
  !
  !  output:
  !  ------
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime, nnode, ntrndof
  real(8), dimension(ndime,ndime), intent(in) :: locbvec
  real(8), dimension(ntrndof,nnode), intent(in) :: elocnodv

  real(8), dimension(ntrndof,nnode), intent(out) :: eglbnodv
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(ntrndof,1) :: glbvec, locvec
  real(8), dimension(ntrndof,ntrndof) :: transmat

  ! loop index
  integer :: inode
  ! ====================================

  ! initialize
  !eglbnodv(:,:)= 0.0d0
  call gettransmat(ndime,ntrndof,locbvec, transmat) !pjsa

  ! convert global to local
  ! -----------------------
  do inode=1, nnode

     ! get global nodal vector
     locvec(1:ntrndof,1)= elocnodv(1:ntrndof,inode)

     ! convert local to global
     !call loc2glbv(ndime,ntrndof,locbvec,locvec, glbvec)
     call dgemv('t',ntrndof,ntrndof,1.0d0,transmat,ntrndof,locvec,1,0.0d0,glbvec,1) !pjsa
        ! input : ndime,ntrndof,locbvec,locvec
        ! output : glbvec

     ! set local nodal vector
     eglbnodv(1:ntrndof,inode)= glbvec(1:ntrndof,1)     

  end do


  return
end subroutine loc2glbnodv



subroutine glb2locv(ndime,ntrndof,locbvec,glbvec, locvec)
  !=======================================================================
  !  glb2locv = transform global vector to local vector
  !
  !             note:
  !             ----
  !             belytschko, lin and tsay, CMAME, 1984, vol. 42, pp. 225-251
  !             explicit algorithms for the nonlinear dynamics of sheels
  !             (see, apendix a)
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : dimension of base vector
  !
  !  ntrndof : number of dof to transform
  !
  !  locbvec(ndime,ndime) : local base vector
  !                     
  !  glbvec(ntrndof,1) : global vector
  !
  !  output:
  !  ------
  !  locvec(ntrndof,1) : local vector
  !                            
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime, ntrndof
  real(8), dimension(ndime,ndime), intent(in) :: locbvec
  real(8), dimension(ntrndof,1), intent(in) :: glbvec

  real(8), dimension(ntrndof,1), intent(out) :: locvec
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(ntrndof,ntrndof) :: transmat

  ! ====================================

  ! initialize
  !locvec(:,:)= 0.0d0

  ! get transformation matrix
  call gettransmat(ndime,ntrndof,locbvec, transmat)
     ! input : ndime,ntrndof,locbvec
     ! output : transmat


  ! transform vector: global -> local : loc_vec = a * glb_vec
  ! ----------------
  !call matprd(ntrndof,ntrndof,0, ntrndof,1,0, ntrndof,1, transmat,glbvec, locvec)
  call dgemv('n',ntrndof,ntrndof,1.0d0,transmat,ntrndof,glbvec,1,0.0d0,locvec,1)
     ! input : ntrndof,ntrndof,0, ntrndof,1,0, ntrndof,1, transmat,glbvec
     ! output : locvec



  return
end subroutine glb2locv



subroutine loc2glbv(ndime,ntrndof,locbvec,locvec, glbvec)
  !=======================================================================
  !  loc2glbv = transform local vector to global vector
  !
  !             note:
  !             ----
  !             belytschko, lin and tsay, CMAME, 1984, vol. 42, pp. 225-251
  !             explicit algorithms for the nonlinear dynamics of sheels
  !             (see, apendix a)
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : dimension of base vector
  !
  !  ntrndof : number of dof to transform
  !
  !  locbvec(ndime,ndime) : local base vector
  !
  !  locvec(ndime,1) : local vector
  !
  !  output:
  !  ------
  !  glbvec(ndime,1) : global vector
  !                            
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime, ntrndof
  real(8), dimension(ndime,ndime), intent(in) :: locbvec
  real(8), dimension(ntrndof,1), intent(in) :: locvec

  real(8), dimension(ntrndof,1), intent(out) :: glbvec
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(ntrndof,ntrndof) :: transmat

  ! ====================================

  ! initialize
  !glbvec(:,:)= 0.0d0

  ! get transformation matrix
  call gettransmat(ndime,ntrndof,locbvec, transmat)
     ! input : ndime,ntrndof,locbvec
     ! output : transmat

  ! transform vector: local -> global : glb_vec = a^t * loc_vec
  ! ----------------
  !call matprd(ntrndof,ntrndof,1, ntrndof,1,0, ntrndof,1, transmat,locvec, glbvec)
  call dgemv('t',ntrndof,ntrndof,1.0d0,transmat,ntrndof,locvec,1,0.0d0,glbvec,1)
     ! input : ntrndof,ntrndof,1, ntrndof,1,0, ntrndof,1, transmat,locvec
     ! output : glbvec


  return
end subroutine loc2glbv



subroutine loc2loc0v(ndime,ntrndof,locbvec,locbvec0,locvec, locvec0)
  !=======================================================================
  !  loc2loc0v = transform current local vector to initial local vector
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : dimension of base vector
  !
  !  ntrndof : number of dof to transform
  !
  !  locbvec(ndime,ndime) : current local base vector
  !
  !  locbvec0(ndime,ndime) : initial local base vector
  !
  !  locvec(ndime,1) : local vector
  !
  !  output:
  !  ------
  !  locvec0(ndime,1) : initial local vector
  !                            
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime, ntrndof
  real(8), dimension(ndime,ndime), intent(in) :: locbvec
  real(8), dimension(ndime,ndime), intent(in) :: locbvec0
  real(8), dimension(ntrndof,1), intent(in) :: locvec

  real(8), dimension(ntrndof,1), intent(out) :: locvec0
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(ntrndof,ntrndof) :: transmat
  real(8), dimension(ntrndof,ntrndof) :: transmat0
  real(8), dimension(ntrndof,1) :: glbvec


  ! loop index
  integer :: idime, jdime
  ! ====================================

  ! initialize
  locvec0(:,:)= 0.0d0


  ! get transformation matrix
  call gettransmat(ndime,ntrndof,locbvec, transmat)
     ! input : ndime,ntrndof,locbvec
     ! output : transmat

  ! transform vector: local -> global : glb_vec = a^t * loc_vec
  ! ----------------
  call matprd(ntrndof,ntrndof,1, ntrndof,1,0, ntrndof,1, transmat,locvec, glbvec)
     ! input : ntrndof,ntrndof,1, ntrndof,1,0, ntrndof,1, transmat,locvec
     ! output : glbvec


  ! get transformation matrix
  call gettransmat(ndime,ntrndof,locbvec0, transmat0)
     ! input : ndime,ntrndof,locbvec0
     ! output : transmat0


  ! transform vector: global -> local0 : loc_vec0 = a0 * glb_vec
  ! ----------------
  call matprd(ntrndof,ntrndof,0, ntrndof,1,0, ntrndof,1, transmat0,glbvec, locvec0)
     ! input : ntrndof,ntrndof,0, ntrndof,1,0, ntrndof,1, transmat0,glbvec
     ! output : locvec0



  return
end subroutine loc2loc0v




subroutine getlocbvecbt(ecurn, locbvec)
  !=======================================================================
  !  getlocbvecbt = compute co-rotational local base vector of BT shell
  !
  !                 note:
  !                 ----
  !                 belytschko, wong and chiang, CMAME, 1992, vol. 96, pp. 93-107
  !                 advances in one point quadrature shell elements
  !                 (direction of local e1, e2 vector is redefinend)
  !
  !                 belytschko and leviathan, CMAME, 1994, vol. 113, pp. 321-350
  !                 physical stabilization of the 4 node shell element with one point quadrature
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ecurn(3,4) : current element nodal coordinate data
  !               x= u + X
  !
  !  output:
  !  ------
  !  locbvec(3,3) : local base vector
  !                            
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(3,4), intent(in) :: ecurn

  real(8), dimension(3,3), intent(out) :: locbvec
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(3,1) :: mid1, mid2, mid3, mid4
  real(8), dimension(3,1) :: r13vec, r42vec

  real(8), dimension(3,1) :: e1vec, e2vec, e3vec

  ! loop index
  integer :: inode, idime
  ! ====================================

  ! initialize
  locbvec(:,:)= 0.0d0

 
  !                 x2_loc
  !                 ^
  !                 |
  !                 |
  !                 |
  !               (mid3)
  !        4--------o--------3
  !        |    r13 ^        |
  !        |        |        |
  !        |        |  r42vec|
  ! (mid4) o --------------> o (mid2)  ---> x1_loc
  !        |        |        |
  !        |        |        |
  !        |        |        |
  !        1--------o--------2
  !               (mid1)
  !

  ! mid side point
  mid1(1:3,1)= 0.50d0*( ecurn(1:3,1)+ecurn(1:3,2) )
  mid2(1:3,1)= 0.50d0*( ecurn(1:3,2)+ecurn(1:3,3) )
  mid3(1:3,1)= 0.50d0*( ecurn(1:3,3)+ecurn(1:3,4) )
  mid4(1:3,1)= 0.50d0*( ecurn(1:3,4)+ecurn(1:3,1) )

  ! mid surface tangent vector: r13, r42
  r13vec(1:3,1)= mid3(1:3,1)-mid1(1:3,1)
  r42vec(1:3,1)= mid2(1:3,1)-mid4(1:3,1)


  ! local 3 axis dir base vector: e3= r42 x r13 / ||r42 x r13|| 
  ! ----------------------------
  call crsprdt3d(1,r42vec,r13vec, e3vec)
    ! input : 1(opt:normalize),r42vec,r13vec
    ! output : e3vec

  ! local 1 axis dir base vector: e1 
  ! ----------------------------
  call unitvec2(3,r42vec, e1vec)
    ! input : 3(ndime),r42vec
    ! output : e1vec

  ! local 2 axis dir base vector: e2 
  ! ----------------------------
  call crsprdt3d(1,e3vec,e1vec, e2vec)
    ! input : 1(opt:normalize),e3vec,e1vec
    ! output : e2vec


  ! set return result
  locbvec(1:3,1)= e1vec(1:3,1)
  locbvec(1:3,2)= e2vec(1:3,1)
  locbvec(1:3,3)= e3vec(1:3,1)



  return
end subroutine getlocbvecbt


subroutine getlocbvectri(ecurn, locbvec)
  !=======================================================================
  !  getlocbvectri = compute co-rotational local base vector of tri shell
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ecurn(3,3) : current element nodal coordinate data
  !               x= u + X
  !
  !  output:
  !  ------
  !  locbvec(3,3) : local base vector
  !                            
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(3,3), intent(in) :: ecurn

  real(8), dimension(3,3), intent(out) :: locbvec
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(3,1) :: r12vec, r13vec

  real(8), dimension(3,1) :: e1vec, e2vec, e3vec

  ! loop index
  integer :: inode, idime
  ! ====================================

  ! initialize
  locbvec(:,:)= 0.0d0

 
  ! edges
  r12vec(1:3,1)= ecurn(1:3,2)-ecurn(1:3,1)
  r13vec(1:3,1)= ecurn(1:3,3)-ecurn(1:3,1)


  ! local 3 axis dir base vector: e3= r12 x r13 / ||r12 x r13|| 
  ! ----------------------------
  call crsprdt3d(1,r12vec,r13vec, e3vec)
    ! input : 1(opt:normalize),r12vec,r13vec
    ! output : e3vec

  ! local 1 axis dir base vector: e1 
  ! ----------------------------
  call unitvec2(3,r12vec, e1vec)
    ! input : 3(ndime),r42vec
    ! output : e1vec

  ! local 2 axis dir base vector: e2 
  ! ----------------------------
  call crsprdt3d(1,e3vec,e1vec, e2vec)
    ! input : 1(opt:normalize),e3vec,e1vec
    ! output : e2vec


  ! set return result
  locbvec(1:3,1)= e1vec(1:3,1)
  locbvec(1:3,2)= e2vec(1:3,1)
  locbvec(1:3,3)= e3vec(1:3,1)



  return
end subroutine getlocbvectri
