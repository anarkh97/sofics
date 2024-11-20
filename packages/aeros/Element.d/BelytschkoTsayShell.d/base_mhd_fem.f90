! ==================================
! finite element method 
! ==================================
!      type                  name                              argument
!      ----                  ----                              --------
! 1.  subroutine          intrplfe2d                     (optele,nnode,psi,eta,nodval, intplval)
! 2.  subroutine          intrplpt                       (optele,ndime,nnode,enodval,psi,eta, pt)
! 3.  subroutine          getparcord2d                   (optele,nnode,nodloc, parnod)
! 4.  subroutine          getshape1d                     (nnode,psi, shap,deriv)
! 5.  subroutine          getshape2d                     (optele,nnode,psi,eta, shap,deriv)
! 6.  subroutine          getshape3d                     (optele,nnode,psi,eta,zet, shap,deriv)
! 7.  subroutine          jacob0                         (ndime,nnode,deriv,ecord, djacob)
! 8.  subroutine          jacob1                         (ndime,nnode,deriv,ecord, djacob,cartd)
! 9.  subroutine          jacob2                         (ndime,nnode,deriv,ecord, djacob,xjaci,cartd)
! 10. subroutine          getgqele                       (optele,ngqdim,ngqpt1,mgaus1, gqpoin,gqweigt)
! 11. subroutine          getgq3dhexa                    (ngqpt1,mgaus1, gqpoin,gqweigt)
! 12. subroutine          getgq2dquad                    (ngqpt1,mgaus1, gqpoin,gqweigt)
! 13. subroutine          getgq2dtri                     (ngaus, posgp,weigp)
! 14. subroutine          getgq1d                        (ngaus, posgp,weigp)
! 15. subroutine          getgqcod                       (ndime,nnode,shap,ecord, gqcod)
! 16. subroutine          getb0mat2d                     (nnode,cartd,ftens, b0mat2d) ! for pk2 stress tensor
! 17. subroutine          getbmat2d                      (nnode,cartd, bmat2d)
! 18. subroutine          getb0mat3d                     (nnode,cartd,ftens, b0mat3d) ! for pk2 stress tensor 
! ---------------------------------------------------------------------------------------------  ! bt shell
! 19. subroutine          getbmat1pt                     (ecordloc, bmat1pt)
! 20. subroutine          getbcmat1pt                    (ecordloc, bcmat1pt)
! 21. subroutine          getbsmat1pt                    (ecordloc, bsmat1pt)
! 22. subroutine          rotprojbt1                     (ecord,edisp,evelo, edisp0,evelo0)
! 23. subroutine          rotprojbt2                     (ecord,edisp,evelo,eaccl, edisp0,evelo0,eaccl0)
! 24. subroutine          rotprojbt3                     (ecord,edisp,efvec0, efvec)
! 25. subroutine          getrotpmatbt                   (inode,ecurn, rotpmat)
! ---------------------------------------------------------------------------------------------
! 26. subroutine          getenergy                      (msize,nsize,delt,velo,fvec,fvecold, wenrg)
! ---------------------------------------------------------------------------------------------
! 27. subroutine          getsq1d                        (nsimp, possp,weisp)
! 28. subroutine          getsq2dquad                    (nsqpt,msimp, sqpoin,sqweigt)
! 29. subroutine          getsqele                       (optele,nsqdim,nsqpt,msimp, sqpoin,sqweigt)
!
! =========================================================================================================



subroutine intrplfe2d(optele,nnode,psi,eta,nodval, intplval)
  !=======================================================================
  !  intrplfe2d= 2d fem interpolation  
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optele : element type option handler
  !
  !  nnode : total number of node per element
  !  
  !  psi,eta : parent coordinate
  !
  !  nodval(1,nnode) : nodal value
  !
  !  output:
  !  ------
  !  intplval : interpolated value
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optele, nnode
  real(8), intent(in) :: psi,eta
  real(8), dimension(1,nnode), intent(in) :: nodval

  real(8), intent(out) :: intplval
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(nnode) :: shap
  real(8), dimension(2,nnode) :: deriv

  ! loop index
  integer :: inode
  ! ====================================
  ! initialize
  intplval= 0.d0

  ! get shape function value at given pasi and eta
  call getshape2d(optele,nnode,psi,eta, shap,deriv)
     ! input : optele,nnode,psi,eta
     ! output : shap,deriv

  ! compute gq position in xy domain
  intplval= 0.0d0 ! initialize
  do inode=1, nnode
     intplval= intplval + shap(inode)*nodval(1,inode)
  end do



  return
end subroutine intrplfe2d





subroutine intrplpt(optele,ndime,nnode,enodval,psi,eta, pt)
  !=======================================================================
  !  intrplpt= interpolate at given point in element
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optele, ndime, nnode : problem dimensions
  !
  !  enodval(ndime,nnode) : element nodal value
  !
  !  psi, eta : psi and eta coordinate
  !
  !  output: 
  !  ------
  !  pt(ndime,1) : interpolated point
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optele
  integer, intent(in) :: ndime, nnode
  real(8), dimension(ndime,nnode), intent(in) :: enodval
  real(8), intent(in) :: psi, eta

  real(8), dimension(ndime,1), intent(out) :: pt
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(1,nnode) :: nodval
  real(8) :: intplval

  ! loop index
  integer :: idime
  ! ====================================
  ! initialize
  pt(:,:)= 0.0d0

  do idime=1, ndime

     ! set nodal value
     nodval(1,1:nnode)= enodval(idime,1:nnode)

     ! interpolate origin in current coordinate    
     call intrplfe2d(optele,nnode,psi,eta,nodval, intplval)
        ! input : optele,nnode,psi,eta,nodval
        ! output : intplval

     ! set result
     pt(idime,1)= intplval

  end do ! end do idime



  return
end subroutine intrplpt





subroutine getparcord2d(optele,nnode,nodloc, parnod)
  !=======================================================================
  !  getparcord2d= get 2d parent nodal coordinate
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optele : element type option
  !
  !  nnode : the total number of node per element
  !
  !  nnodeloc : local node number
  !
  !  output:
  !  ------
  !  parnod : the total number of edge
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optele,nnode
  integer, intent(in) :: nodloc

  real(8), dimension(2,1), intent(out) :: parnod
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(2,4) :: quad4nod
  real(8), dimension(2,3) :: tri3nod
  real(8), dimension(2,6) :: tri6nod

  ! loop index
  ! ====================================

  ! initialize
  parnod(:,:)= 0.0d0

  data quad4nod / -1.0d0,-1.0d0, 1.0d0,-1.0d0, 1.0d0,1.0d0, -1.0d0,1.0d0 / 
  data tri3nod / 0.0d0,0.0d0, 1.0d0,0.0d0, 0.0d0,1.0d0 / 
  data tri6nod / 0.0d0,0.0d0, 1.0d0,0.0d0, 0.0d0,1.0d0, 0.50d0,0.0d0, 0.50d0,0.50d0, 0.0d0,0.50d0 / 

  select case(optele)
  case(1) ! 2d tri
     if ( nnode==3 ) then ! 3 node tri.
        parnod(1:2,1)= tri3nod(1:2,nodloc)
     else if ( nnode==6 ) then ! 6 node tri.
        parnod(1:2,1)= tri6nod(1:2,nodloc)
     end if
    
  case(2) ! 2d quad
     if ( nnode==4 ) then ! 4 node quad.
        parnod(1:2,1)= quad4nod(1:2,nodloc)

     end if
  
  case(3) ! 3d 4 node bt shell
     parnod(1:2,1)= quad4nod(1:2,nodloc)

  end select 



  return
end subroutine getparcord2d 





subroutine getshape1d(nnode,psi, shap,deriv)
  !=======================================================================
  !  getshape1d = compute shape functions and their derivatives
  !               for 1d element
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  psi : psi-coordinate value
  !
  !  optput:
  !  ------
  !  shap(nnode) : point value of shape function at (psi) 
  !  deriv(1,nnode) : point value of shape function's derivation
  !                       in psi-eta domain
  !           deriv(1,i) :  dn_i/d(psi) value of ith node
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nnode
  real(8), intent(in) :: psi

  real(8), dimension(nnode), intent(out) :: shap
  real(8), dimension(1,nnode), intent(out) :: deriv
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: s
  ! ====================================

  ! initialize
  shap(:)= 0.0d0
  deriv(:,:)= 0.0d0

  ! set initial value for psi and eta
  !---------------
  s= psi
  !---------------


  select case(nnode)

  ! 2d line element
  case (2)
     ! set shape function's value
     shap(1)= (1.0d0-s) / 2.0d0
     shap(2)= (1.0d0+s) / 2.0d0

     ! set derivation of shape function's value
     deriv(1,1)= -0.50d0 ! dn/d(psi)
     deriv(1,2)= 0.50d0
    
  case default
     write(*,*) "not avaliable shape function: getshape1d"
     write(nout5,*) "not avaliable shape function: getshape1d"
     stop

  end select


  
  return
end subroutine getshape1d





subroutine getshape2d(optele,nnode,psi,eta, shap,deriv)
  !=======================================================================
  !  getshape2d = compute shape functions and their derivatives
  !               for isoparametric 2-d elements
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optele : element type
  !
  !  nnode : the total number of node per element
  !
  !  psi : psi-coordinate value
  !  eta : eta-coordinate value
  !
  !  optput:
  !  ------
  !  shap(nnode) : point value of shape function at (psi,eta) 
  !  deriv(2,nnode) : point value of shape function's derivation
  !                       in psi-eta domain
  !           deriv(1,i) :  dn_i/d(psi) value of ith node
  !           deriv(2,i) :  dn_i/d(eta) value of ith node
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optele, nnode
  real(8), intent(in) :: psi,eta

  real(8), dimension(nnode), intent(out) :: shap
  real(8), dimension(2,nnode), intent(out) :: deriv
  ! ====================================
  ! local variable
  ! ==============
  character(len=10) :: handler
  real(8) :: s,t,st,ss,tt,s2,t2,sst,stt,st2
  ! ====================================

  ! initialize variable
  shap(:)=0.0d0
  deriv(:,:)=0.0d0

  ! set initial value for psi and eta
  !---------------
  s= psi
  t= eta
  !---------------
  st= psi * eta
  ss= psi**2
  tt= eta**2
  !---------------
  s2= psi * 2.0d0
  t2= eta * 2.0d0
  st2= psi * eta * 2.0d0
  !---------------
  sst= psi * psi * eta
  stt= psi * eta * eta
  !---------------

  ! ------------------------------------------  
  ! 2d element
  if(optele==1 .and. nnode==3) handler='2d_t3'
  if(optele==1 .and. nnode==6) handler='2d_t6'

  if(optele==2 .and. nnode==4) handler='2d_q4'
  if(optele==2 .and. nnode==8) handler='2d_q8'

  ! 3d bt shell element
  if(optele==3 .and. nnode==4) handler='2d_q4'
  ! ------------------------------------------

  select case(handler)

  ! three-nodes triangle element : t3
  ! ----------------------------
  !    eta
  !     ^
  !     |
  ! (3) 0 
  !     | * 
  !     |   * 
  !     |     *
  !     |       * 
  !     |         *
  !     |           *
  !     0-------------0----> psi
  !  (1)              (2)
  !
  case ('2d_t3')
     ! set shape function's value
     shap(1)=1.0d0-s-t
     shap(2)=s
     shap(3)=t

     ! set derivation of shape function's value
     deriv(1,1)=-1.0d0 ! dn/d(psi)
     deriv(1,2)=1.0d0
     deriv(1,3)=0.0d0
     deriv(2,1)=-1.0d0 ! dn/d(eta)
     deriv(2,2)=0.0d0
     deriv(2,3)=1.0d0


  ! six-nodes triangle element : t6
  ! --------------------------
  !    eta
  !     ^
  !     |
  ! (3) 0 
  !     |  * 
  !     |    * 
  !     |      *
  ! [6] 0        0  [5]
  !     |          * 
  !     |            *
  !     |              *
  !     0--------0-------0----> psi
  !  (1)        [4]       (2)
  !

  case('2d_t6')
     ! set shape function's value
     shap(1)=1.0d0-3.0d0*(s+t)+4.0d0*st+2.0d0*(ss+tt)
     shap(2)=s*(2.0d0*s-1.0d0)
     shap(3)=t*(2.0d0*t-1.0d0)
     shap(4)=4.0d0*s*(1.0d0-s-t) ! mid node
     shap(5)=4.0d0*s*t           ! mid node
     shap(6)=4.0d0*t*(1.0d0-s-t) ! mid node

     ! set derivation of shape function's value
     deriv(1,1)=4.0d0*(s+t)-3.0d0 ! dn/d(psi)
     deriv(1,2)=4.0d0*s-1.0d0 
     deriv(1,3)=0.0d0
     deriv(1,4)=4.0d0*(1.0d0-t-2.0d0*s) ! mid node
     deriv(1,5)=4.0d0*t                 ! mid node
     deriv(1,6)=-4.0d0*t                ! mid node

     deriv(2,1)=4.0d0*(s+t)-3.0d0 ! dn/d(eta)
     deriv(2,2)=0.0d0
     deriv(2,3)=4.0d0*t-1.0d0
     deriv(2,4)=-4.0d0*s                ! mid node
     deriv(2,5)=4.0d0*s                 ! mid node
     deriv(2,6)=4.0d0*(1.0d0-s-2.0d0*t) ! mid node



  ! four-nodes quadrialteral element : q4
  ! --------------------------------
  !             eta
  !              ^
  !              |
  !  (4)         |         (3) 
  !     0--------|--------0 
  !     |        |        |
  !     |        |        |
  !     |        |        |
  !     |        ---------|---------> psi
  !     |                 |
  !     |                 |
  !     |                 |
  !     0-----------------0
  !  (1)                   (2)
  !

  case('2d_q4')
     ! shape functions for 4 noded element
     shap(1)=(1.d0-t-s+st)/4.0d0
     shap(2)=(1.d0-t+s-st)/4.0d0
     shap(3)=(1.d0+t+s+st)/4.0d0
     shap(4)=(1.d0+t-s-st)/4.0d0

     ! shape function derivatives for 4 noded element
     deriv(1,1)=(-1.d0+t)/4.0d0   ! dn/d(psi)
     deriv(1,2)=( 1.d0-t)/4.0d0
     deriv(1,3)=( 1.d0+t)/4.0d0
     deriv(1,4)=(-1.d0-t)/4.0d0
     deriv(2,1)=(-1.d0+s)/4.0d0   ! dn/d(eta)
     deriv(2,2)=(-1.d0-s)/4.0d0
     deriv(2,3)=( 1.d0+s)/4.0d0
     deriv(2,4)=( 1.d0-s)/4.0d0


  ! eight-nodes quadrialteral element : q8
  ! ---------------------------------
  !             eta
  !              ^
  !              |
  !              |
  !  (4)        [7]        (3) 
  !     0--------0--------0 
  !     |        |        |
  !     |        |        |
  !     |        |        |
  ! [8] 0        ---------0---------> psi
  !     |                 | [6]
  !     |                 |
  !     |                 |
  !     0--------0--------0
  !  (1)        [5]        (2)
  !

  case('2d_q8')
     shap(1)=(-1.0d0+st+ss+tt-sst-stt)/4.0d0
     shap(2)=(-1.0d0-st+ss+tt-sst+stt)/4.0d0
     shap(3)=(-1.0d0+st+ss+tt+sst+stt)/4.0d0
     shap(4)=(-1.0d0-st+ss+tt+sst-stt)/4.0d0
     shap(5)=(1.0d0-t-ss+sst)/2.0d0 ! mid node
     shap(6)=(1.0d0+s-tt-stt)/2.0d0 ! mid node
     shap(7)=(1.0d0+t-ss-sst)/2.0d0 ! mid node
     shap(8)=(1.0d0-s-tt+stt)/2.0d0 ! mid node

     ! shape function derivatives
     deriv(1,1)=(t+s2-st2-tt)/4.0d0 ! dn/d(psi)
     deriv(1,2)=(-t+s2-st2+tt)/4.0d0
     deriv(1,3)=(t+s2+st2+tt)/4.0d0
     deriv(1,4)=(-t+s2+st2-tt)/4.0d0
     deriv(1,5)=-s+st             ! mid node
     deriv(1,6)=(1.0d0-tt)/2.0d0  ! mid node
     deriv(1,7)=-s-st             ! mid node
     deriv(1,8)=(-1.0d0+tt)/2.0d0 ! mid node
     deriv(2,1)=(s+t2-ss-st2)/4.0d0 ! dn/d(eta)
     deriv(2,2)=(-s+t2-ss+st2)/4.0d0
     deriv(2,3)=(s+t2+ss+st2)/4.0d0
     deriv(2,4)=(-s+t2+ss-st2)/4.0d0
     deriv(2,5)=(-1.0d0+ss)/2.0d0 ! mid node
     deriv(2,6)=-t-st             ! mid node
     deriv(2,7)=(1.0d0-ss)/2.0d0  ! mid node
     deriv(2,8)=-t+st             ! mid node
  
  case default
     write(*,*) "not avaliable shape function: getshape2d"
     write(nout5,*) "not avaliable shape function: getshape2d"
     stop

  end select


  
  return
end subroutine getshape2d





subroutine getshape3d(optele,nnode,psi,eta,zet, shap,deriv)
  !=======================================================================
  !  getshape2d = compute shape functions and their derivatives
  !               for isoparametric 3-d elements
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optele : element type
  !
  !  nnode : the total number of node per element
  !
  !  psi : psi-coordinate value
  !  eta : eta-coordinate value
  !  zet : zeta-coordinate value
  !
  !  optput:
  !  ------
  !  shap(nnode) : point value of shape function at (psi,eta,zet) 
  !  deriv(3,nnode) : point value of shape function's derivation
  !                   in psi-eta-zet domain
  !           deriv(1,i) :  dn_i/d(psi) value of ith node
  !           deriv(2,i) :  dn_i/d(eta) value of ith node
  !           deriv(3,i) :  dn_i/d(zet) value of ith node
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optele, nnode
  real(8), intent(in) :: psi, eta, zet

  real(8), dimension(nnode), intent(out) :: shap
  real(8), dimension(3,nnode), intent(out) :: deriv
  ! ====================================
  ! local variable
  ! ==============
  character(len=10) :: handler

  real(8), dimension(3,8) :: stzq8i
  real(8) :: s, t, z
  real(8) :: si, ti, zi


  ! loop index
  integer :: inode
  ! ====================================

  ! initialize variable
  shap(:)= 0.0d0
  deriv(:,:)= 0.0d0

  ! set initial value for psi, eta and zet
  !---------------
  s= psi
  t= eta
  z= zet
  !---------------

  ! set nodal position in psi,eta-zet domain
  ! ----------------------------------------
  data stzq8i / -1.0d0,-1.0d0,-1.0d0, 1.0d0,-1.0d0,-1.0d0, 1.0d0,1.0d0,-1.0d0, -1.0d0,1.0d0,-1.0d0, &
                -1.0d0,-1.0d0,1.0d0, 1.0d0,-1.0d0,1.0d0,  1.0d0,1.0d0,1.0d0, -1.0d0,1.0d0,1.0d0 /

  ! ------------------------------------------  
  if( optele==4 .and. nnode==8 ) handler='3d_q8'
  ! ------------------------------------------

  select case(handler)
  case('3d_q8')
     do inode=1, 8

        ! set nodal psotion
        si= stzq8i(1,inode)
        ti= stzq8i(2,inode) 
        zi= stzq8i(3,inode)

        ! get nodal shape function
        shap(inode)= 1.0d0 / 8.0d0 * ( 1.0d0 + si * s ) * ( 1.0d0 + ti * t ) * ( 1.0d0 + zi * z )

        ! shape function derivatives
        deriv(1,inode)= 1.0d0 / 8.0d0 * si * ( 1.0d0 + ti * t ) * ( 1.0d0 + zi * z ) ! dn/d(psi)
        deriv(2,inode)= 1.0d0 / 8.0d0 * ti * ( 1.0d0 + si * s ) * ( 1.0d0 + zi * z ) ! dn/d(eta)
        deriv(3,inode)= 1.0d0 / 8.0d0 * zi * ( 1.0d0 + si * s ) * ( 1.0d0 + ti * t ) ! dn/d(zet)

     end do

  case default
     write(*,*) "not avaliable shape function: getshape3d"
     write(nout5,*) "not avaliable shape function: getshape3d"
     stop

  end select


  
  return
end subroutine getshape3d





subroutine jacob0(ndime,nnode,deriv,ecord, djacob)
  !=======================================================================
  ! jacob0 = compute determinant of jacobian matrix
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : problem dimension
  !
  !  nnode : the total number of node per element
  !
  !  deriv(ndime, nnode) : point value of shape function's derivation in psi-eta domain
  !           deriv(1,i) :  dn_i/d(psi) value of ith node
  !           deriv(2,i) :  dn_i/d(eta) value of ith node
  !           deriv(3,i) :  dn_i/d(zeta) value of ith node
  !
  !  ecord(ndime,nnode) : nodal coordinates of current element
  !           ecord(1,i) : ith node's x-coordinate
  !           ecord(2,i) : ith node's y-coordinate
  !           ecord(3,i) : ith node's z-coordinate
  !
  !  output:
  !  ------
  !  djacob : determinant of jacobian matrix     
  !
  ! ======================================================================
 
  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime, nnode
  real(8), dimension(ndime,nnode), intent(in) :: deriv
  real(8), dimension(ndime,nnode), intent(in) :: ecord

  real(8), intent(out) :: djacob
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(ndime,ndime) :: xjacm ! jacobian matrix
  real(8) :: getdet
  ! ====================================

  ! initialize
  djacob= 0.0d0

  ! calculate jacobian matrix wrt refernce configuration : xjacm
  ! xjacm(ndime,ndime)=deriv(ndime,nnode).ecord(ndime,nnode)^t
  call matprd(ndime,nnode,0, ndime, nnode,1, ndime,ndime, deriv,ecord, xjacm)
     ! input : ndime,nnode,0, ndime, nnode,1, ndime,ndime, deriv,ecord
     ! output : xjacm

  ! check determinant of jacobian matrix : djacb
  djacob= getdet(ndime,xjacm)


  return
end subroutine jacob0





subroutine jacob1(ndime,nnode,deriv,ecord, djacob,cartd)
  !=======================================================================
  ! jacob1 = compute jacobian matrix and its determinant
  !          also, compute derivation of shape function in cartesian domain
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : problem dimension
  !
  !  nnode : the total number of node per element
  !
  !  deriv(ndime, nnode) : point value of shape function's derivation in psi-eta domain
  !           deriv(1,i) :  dn_i/d(psi) value of ith node
  !           deriv(2,i) :  dn_i/d(eta) value of ith node
  !           deriv(3,i) :  dn_i/d(zeta) value of ith node
  !
  !  ecord(ndime,nnode) : nodal coordinates of current element
  !           ecord(1,i) : ith node's x-coordinate
  !           ecord(2,i) : ith node's y-coordinate
  !           ecord(3,i) : ith node's z-coordinate
  !
  !  output:
  !  ------
  !  djacob : determinant of jacobian matrix     
  !
  !  cartd(ndime,nnode) : point value of shape function's derivation in x-y domain
  !           cartd(1,i) : dn_i/dx
  !           cartd(2,i) : dn_i/dy
  !           cartd(3,i) : dn_i/dz
  !                            
  ! ======================================================================
 
  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime, nnode
  real(8), dimension(ndime,nnode), intent(in) :: deriv
  real(8), dimension(ndime,nnode), intent(in) :: ecord

  real(8), intent(out) :: djacob
  real(8), dimension(ndime,nnode), intent(out) :: cartd
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(ndime,ndime) :: xjacm ! jacobian matrix
  real(8), dimension(ndime,ndime) :: xjaci ! inverse of jacobian matrix
  real(8) :: getdet
  ! ====================================

  ! initialize
  djacob= 0.0d0
  cartd(:,:)= 0.0d0

  ! calculate jacobian matrix wrt refernce configuration : xjacm
  ! xjacm(ndime,ndime)=deriv(ndime,nnode).ecord(ndime,nnode)^t
  call matprd(ndime,nnode,0, ndime,nnode,1, ndime,ndime, deriv,ecord, xjacm)
     ! input : ndime,nnode,0, ndime,nnode,1, ndime,ndime, deriv,ecord
     ! output : xjacm

  ! check determinant of jacobian matrix : djacb
  djacob= getdet(ndime,xjacm)

  ! calculate inverse of jacobian matrix : xjaci
  call getinv(ndime,xjacm, xjaci)
    ! input : ndime,xjacm
    ! output : xjaci

  ! determine derivation value at sampling point in physical coordinate
  ! cartd(ndime,nnode)=xjaci(ndime,ndime).deriv(ndime,nnode)
  call matprd(ndime,ndime,0, ndime,nnode,0, ndime,nnode, xjaci,deriv, cartd)
    ! input : ndime,ndime,0, ndime, nnode,0, ndime,nnode, xjaci,deriv
    ! output : cartd

  return
end subroutine jacob1




subroutine jacob2(ndime,nnode,deriv,ecord, djacob,xjaci,cartd)
  !=======================================================================
  ! jacob2 = evaluates the jacobian matrix and its determinant
  !          also, compute derivation of shape function in cartesian domain
  !          additionally return inverse of jacobian matrix
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : problem dimension
  !
  !  nnode : the total number of node per element
  !
  !  deriv(ndime, nnode) : point value of shape function's derivation in psi-eta domain
  !           deriv(1,i) :  dn_i/d(psi) value of ith node
  !           deriv(2,i) :  dn_i/d(eta) value of ith node
  !           deriv(3,i) :  dn_i/d(zeta) value of ith node
  !
  !  ecord(ndime,nnode) : nodal coordinates of current element
  !           ecord(1,i) : ith node's x-coordinate
  !           ecord(2,i) : ith node's y-coordinate
  !           ecord(3,i) : ith node's z-coordinate
  !
  !  output:
  !  ------
  !  djacob : determinant of jacobian matrix
  !
  !  xjaci(ndime,ndime) : inverse of jacobian matrix
  !
  !  cartd(ndime,nnode) : point value of shape function's derivation in x-y domain
  !           cartd(1,i) : dn_i/dx
  !           cartd(2,i) : dn_i/dy
  !           cartd(3,i) : dn_i/dz
  !                            
  ! ======================================================================
 
  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime, nnode
  real(8), dimension(ndime,nnode), intent(in) :: deriv
  real(8), dimension(ndime,nnode), intent(in) :: ecord

  real(8), intent(out) :: djacob
  real(8), dimension(ndime,ndime), intent(out) :: xjaci
  real(8), dimension(ndime,nnode), intent(out) :: cartd
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(ndime,ndime) :: xjacm ! jacobian matrix
  real(8) :: getdet
  ! ====================================

  ! initialize
  djacob= 0.0d0
  cartd(:,:)= 0.0d0

  ! calculate jacobian matrix wrt refernce configuration : xjacm
  ! xjacm(ndime,ndime)=deriv(ndime,nnode).ecord(ndime,nnode)T
  call matprd(ndime,nnode,0, ndime, nnode,1, ndime,ndime, deriv,ecord, xjacm)
     ! input : ndime,nnode,0, ndime, nnode,1, ndime,ndime, deriv,ecord
     ! output : xjacm

  ! check determinant of jacobian matrix : djacb
  djacob= getdet(ndime,xjacm)

  ! calculate inverse of jacobian matrix : xjaci
  call getinv(ndime,xjacm, xjaci)
    ! input : ndime,xjacm
    ! output : xjaci

  ! determine derivation value at sampling point in x-y coordinate
  ! cartd(ndime,nnode)=xjaci(ndime,ndime).deriv(ndime,nnode)
  call matprd(ndime,ndime,0, ndime,nnode,0, ndime,nnode, xjaci,deriv, cartd)
    ! input : ndime,ndime,0, ndime, nnode,0, ndime,nnode, xjaci,deriv
    ! output : cartd

  return
end subroutine jacob2





subroutine getgqele(optele,ngqdim,ngqpt1,mgaus1, gqpoin,gqweigt)
  !=======================================================================
  !  getgqele = get rule for given gq dimension and element
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optele : element type option handler
  !
  !  ngqdim : dimension of gq
  !
  !  ngqpt1 : 1d gq rule
  !
  !  mgaus1 : total number of gq
  !           for tri. ele. : mgaus1= ngqpt1
  !           for quad. ele. : mgaus1= ngqpt1**2 :psi and eta direction
  !
  !  output:
  !  ------
  !  gqpoin(2,mgaus1) : gq position
  !
  !  gqweigt(mgaus1) : gq weight
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optele
  integer, intent(in) :: ngqdim
  integer, intent(in) :: ngqpt1, mgaus1

  real(8), dimension(ngqdim,mgaus1), intent(out) :: gqpoin 
  real(8), dimension(mgaus1), intent(out) :: gqweigt 
  ! ====================================
  ! local variable
  ! ==============

  ! loop index
  ! ====================================

  ! initialize
  gqpoin(:,:)= 0.0d0
  gqweigt(:)= 0.0d0
  

  if ( ngqdim== 2 .and. optele == 1 ) then ! 2d tri
     ! get gq rule for 2d triangular domain
     call getgq2dtri(ngqpt1, gqpoin,gqweigt)
        ! input : ngaus1
        ! output : gqpoin,gqweigt

  else if ( ngqdim== 2 .and. optele == 2 ) then ! 2d quad
     ! get gq rule for 2d quad
     call getgq2dquad(ngqpt1,mgaus1, gqpoin,gqweigt)
        ! iput : ngqpt1,mgaus1
        ! output : gqpoin,gqweigt

  else if ( ngqdim== 2 .and. optele == 3 ) then ! 3d bt shell (in plane)
     ! get gq rule for 2d quad
     call getgq2dquad(ngqpt1,mgaus1, gqpoin,gqweigt)
        ! iput : ngqpt1,mgaus1
        ! output : gqpoin,gqweigt

  else if ( ngqdim== 3 .and. optele == 4 ) then ! 3d brick
     ! get gq rule for 3d brick
     call getgq3dhexa(ngqpt1,mgaus1, gqpoin,gqweigt)
        ! iput : ngqpt1,mgaus1
        ! output : gqpoin,gqweigt

  else
     write(*,*) "not available: getgqele"
     write(nout5,*) "not available: getgqele"
     stop

  end if



  return
end subroutine getgqele





subroutine getgq3dhexa(ngqpt1,mgaus1, gqpoin,gqweigt)
  !=======================================================================
  !  getgq3dhexa = get 3d gq rule for hexahedron
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ngqpt1 : 1d gq rule
  !
  !  mgaus1 : total number of gq
  !           mgaus1 = ngqpt1**3 : psi, eta and zeta direction
  !
  !  output:
  !  ------
  !  gqpoin(3,mgaus1) : gq position
  !
  !  gqweigt(mgaus1) : gq weight
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ngqpt1, mgaus1

  real(8), dimension(3,mgaus1), intent(out) :: gqpoin 
  real(8), dimension(mgaus1), intent(out) :: gqweigt 
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(ngqpt1) :: posgp
  real(8), dimension(ngqpt1) :: weigp

  integer :: mgaus

  ! loop index
  integer :: igaus, jgaus, kgaus
  ! ====================================

  ! initialize
  gqpoin(:,:)= 0.0d0
  gqweigt(:)= 0.0d0
  

  ! get gq rule for 1d line domain
  call getgq1d(ngqpt1, posgp,weigp)
     ! iput : ngqpt1
     ! output : posgp,weigp

  ! set total gq point and it's weight
  mgaus=0 ! initialize
  do igaus=1, ngqpt1
     do jgaus=1, ngqpt1
        do kgaus=1, ngqpt1

           ! increase counter
           mgaus= mgaus+1

           ! set 3d position
           gqpoin(1,mgaus)= posgp(igaus)
           gqpoin(2,mgaus)= posgp(jgaus)
           gqpoin(3,mgaus)= posgp(kgaus)

           ! set weight
           gqweigt(mgaus)= weigp(igaus) * weigp(jgaus) * weigp(kgaus)

        end do ! igaus
     end do ! jgaus
  end do ! kgaus

  ! error check
  if( mgaus /= mgaus1) then
     write(*,*) "inconsistency in gq rule: getgq3dhexa"
     write(nout5,*) "inconsistency in gq rule: getgq3dhexa"
     stop

  end if



  return
end subroutine getgq3dhexa





subroutine getgq2dquad(ngqpt1,mgaus1, gqpoin,gqweigt)
  !=======================================================================
  !  getgq2dquad = get 2d gq rule for quad
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ngqpt1 : 1d gq rule
  !
  !  mgaus1 : total number of gq
  !           mgaus1 = ngqpt1**2 : psi and eta direction
  !
  !  output:
  !  ------
  !  gqpoin(2,mgaus1) : gq position
  !
  !  gqweigt(mgaus1) : gq weight
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ngqpt1, mgaus1

  real(8), dimension(2,mgaus1), intent(out) :: gqpoin 
  real(8), dimension(mgaus1), intent(out) :: gqweigt 
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(ngqpt1) :: posgp
  real(8), dimension(ngqpt1) :: weigp

  integer :: mgaus

  ! loop index
  integer :: igaus, jgaus
  ! ====================================

  ! initialize
  gqpoin(:,:)= 0.0d0
  gqweigt(:)= 0.0d0
  

  ! get gq rule for 1d line domain
  call getgq1d(ngqpt1, posgp,weigp)
     ! iput : ngqpt1
     ! output : posgp,weigp

  ! set total gq point and it's weight
  mgaus=0 ! initialize
  do igaus=1, ngqpt1
     do jgaus=1, ngqpt1

        ! increase counter
        mgaus= mgaus+1

        ! set 2d position
        gqpoin(1,mgaus)= posgp(igaus)
        gqpoin(2,mgaus)= posgp(jgaus)

        ! set weight
        gqweigt(mgaus)= weigp(igaus) * weigp(jgaus)

     end do ! jgaus
  end do ! igaus

  ! error check
  if( mgaus /= mgaus1) then
     write(*,*) "inconsistency in gq rule: getgq2dquad"
     write(nout5,*) "inconsistency in gq rule: getgq2dquad"
     stop

  end if



  return
end subroutine getgq2dquad





subroutine getgq2dtri(ngaus, posgp,weigp)
  !=======================================================================
  !  getgq2dtri = sets up the gauss-legendre integration constants
  !               for 2d triangular domain
  !
  !  input :
  !  -------
  !  ngaus : the number of g.q 
  !
  !  output :
  !  --------
  !  posgp(2,ngaus) : g.q position for 2d triangular element in psi-eta domain
  !
  !  weigp(ngasu) : g.q weight values
  !
  !=======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ngaus

  real(8), dimension(2,ngaus), intent(out) :: posgp 
  real(8), dimension(ngaus), intent(out) :: weigp 
  ! ====================================

  ! initialize
  posgp(:,:)=0.0d0
  weigp(:)=0.0d0 


  ! set gauss quadrature point ! 1st order
  if(ngaus==1) then
     posgp(1,1) = 0.33333333333330d0
     posgp(2,1) = 0.33333333333330d0

     weigp(1)   = 1.0d0

  elseif(ngaus==3) then ! 2nd order
     posgp(1,1) = 0.16666666666670d0
     posgp(2,1) = 0.16666666666670d0
     posgp(1,2) = 0.66666666666670d0
     posgp(2,2) = 0.16666666666670d0
     posgp(1,3) = 0.16666666666670d0
     posgp(2,3) = 0.66666666666670d0

     weigp(1)   = 0.33333333333330d0 
     weigp(2)   = 0.33333333333330d0 
     weigp(3)   = 0.33333333333330d0   

  elseif(ngaus==7) then ! 5th order
     posgp(1,1) = 0.10128650732350d0
     posgp(2,1) = 0.10128650732350d0
     posgp(1,2) = 0.79742698535310d0
     posgp(2,2) = 0.10128650732350d0
     posgp(1,3) = 0.10128650732350d0
     posgp(2,3) = 0.79742698535310d0
     posgp(1,4) = 0.47014206410510d0
     posgp(2,4) = 0.05971587178980d0
     posgp(1,5) = 0.47014206410510d0
     posgp(2,5) = 0.47014206410510d0
     posgp(1,6) = 0.05971587178980d0
     posgp(2,6) = 0.47014206410510d0 
     posgp(1,7) = 0.33333333333330d0 
     posgp(2,7) = 0.33333333333330d0

     weigp(1)   = 0.12593918054480d0 
     weigp(2)   = 0.12593918054480d0 
     weigp(3)   = 0.12593918054480d0 
     weigp(4)   = 0.13239415278850d0
     weigp(5)   = 0.13239415278850d0
     weigp(6)   = 0.13239415278850d0
     weigp(7)   = 0.2250d0

  else
     write(*,*) "not available: getgq2dtri"
     write(nout5,*) "not available: getgq2dtri"
     stop

  endif

  ! area fraction
  weigp(:)= weigp(:) * 0.50d0



  return
end subroutine getgq2dtri




subroutine getgq1d(ngaus, posgp,weigp)
  !=======================================================================
  !  getgq1d = sets up the gauss-legendre integration constants for 1d line
  !
  !  input :
  !  -------
  !  ngaus : the number of g.q 
  !
  !  output :
  !  --------
  !  posgp(ngaus) : g.q position in psi-eta domain
  !
  !  weigp(ngasu) : g.q weight values
  !
  !=======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ngaus

  real(8), dimension(ngaus), intent(out) :: posgp
  real(8), dimension(ngaus), intent(out) :: weigp
  ! ====================================
  ! local variable
  ! ==============

  ! ====================================

  ! initialize
  posgp(:)=0.d0
  weigp(:)=0.d0 


  ! set 1d mension gauss point and it's weight value
  if( ngaus==1 )then
     posgp(1) =  0.0d0

     weigp(1) =  2.0d0

  else if( ngaus==2 )then
     posgp(1) = -0.5773502691896260d0
     posgp(2) =  0.5773502691896260d0

     weigp(1) =  1.0d0
     weigp(2) =  1.0d0

  else if( ngaus==3 )then
     posgp(1) = -0.7745966692414830d0
     posgp(2) =  0.0d0
     posgp(3) =  0.7745966692414830d0

     weigp(1) =  0.55555555555555560d0
     weigp(2) =  0.88888888888888890d0
     weigp(3) =  0.55555555555555560d0

  else if( ngaus==4 )then
     posgp(1) = -0.8611363115940530d0
     posgp(2) = -0.3399810435848560d0
     posgp(3) =  0.3399810435848560d0
     posgp(4) =  0.8611363115940530d0

     weigp(1) =  0.3478548451374540d0
     weigp(2) =  0.6521451548625460d0
     weigp(3) =  0.6521451548625460d0
     weigp(4) =  0.3478548451374540d0

  else if( ngaus==5 )then
     posgp(1) = -0.9061798459386640d0
     posgp(2) = -0.5384693101056830d0
     posgp(3) =  0.0d0
     posgp(4) =  0.5384693101056830d0
     posgp(5) =  0.9061798459386640d0

     weigp(1) =  0.2369268850561890d0
     weigp(2) =  0.4786286704993660d0
     weigp(3) =  0.5688888888988890d0
     weigp(4) =  0.4786286704993660d0
     weigp(5) =  0.2369268850561890d0

  else if( ngaus==6 )then
     posgp(1) = -0.9324695142031520d0
     posgp(2) = -0.6612093864662650d0
     posgp(3) = -0.2386191860831970d0
     posgp(4) =  0.2386191860831970d0
     posgp(5) =  0.6612093864662650d0
     posgp(6) =  0.9324695142031520d0

     weigp(1) =  0.1713244923791700d0
     weigp(2) =  0.3607615730481390d0
     weigp(3) =  0.4679139345726910d0
     weigp(4) =  0.4679139345726910d0
     weigp(5) =  0.3607615730481390d0
     weigp(6) =  0.1713244923791700d0

  else if( ngaus==7 )then
     posgp(1) = -0.949107912300000d0
     posgp(2) = -0.741531185600000d0
     posgp(3) = -0.405845151400000d0
     posgp(4) =  0.d0
     posgp(5) =  0.405845151400000d0
     posgp(6) =  0.741531185600000d0
     posgp(7) =  0.949107912300000d0

     weigp(1) =  0.129484966200000d0
     weigp(2) =  0.279705391500000d0
     weigp(3) =  0.381830050500000d0
     weigp(4) =  0.417959183700000d0
     weigp(5) =  0.381830050500000d0
     weigp(6) =  0.279705391500000d0
     weigp(7) =  0.129484966200000d0

  else if( ngaus==8 )then    
     posgp(1) = -0.9602898564975360d0
     posgp(2) = -0.7966664774136270d0
     posgp(3) = -0.5255324099163290d0
     posgp(4) = -0.1834346424956500d0
     posgp(5) =  0.1834346424956500d0
     posgp(6) =  0.5255324099163290d0
     posgp(7) =  0.7966664774136270d0
     posgp(8) =  0.9602898564975360d0

     weigp(1) =  0.1012285362903760d0
     weigp(2) =  0.2223810344533740d0
     weigp(3) =  0.3137066458778870d0
     weigp(4) =  0.3626837833783620d0
     weigp(5) =  0.3626837833783620d0
     weigp(6) =  0.3137066458778870d0
     weigp(7) =  0.2223810344533740d0
     weigp(8) =  0.1012285362903760d0

  else if( ngaus==9 )then
     posgp(1) = -0.9681602395076260d0
     posgp(2) = -0.8360311073266360d0
     posgp(3) = -0.6133714327005900d0
     posgp(4) = -0.324253423403809d0
     posgp(5) =  0.d0
     posgp(6) =  0.324253423403809d0
     posgp(7) =  0.613371432700590d0
     posgp(8) =  0.836031107326636d0
     posgp(9) =  0.968160239507626d0

     weigp(1) =  0.081274388361574d0
     weigp(2) =  0.180648160694857d0
     weigp(3) =  0.260610696402935d0
     weigp(4) =  0.312347077040003d0
     weigp(5) =  0.330239355001260d0
     weigp(6) =  0.312347077040003d0
     weigp(7) =  0.260610696402935d0
     weigp(8) =  0.180648160694857d0
     weigp(9) =  0.081274388361574d0

  else if( ngaus==10 )then
     posgp(1) = -0.973906528517172d0
     posgp(2) = -0.865063366688984d0
     posgp(3) = -0.679409568299024d0
     posgp(4) = -0.433395394129247d0
     posgp(5) = -0.148874338981631d0
     posgp(6) =  0.148874338981631d0
     posgp(7) =  0.433395394129247d0
     posgp(8) =  0.679409568299024d0
     posgp(9) =  0.865063366688984d0
     posgp(10)=  0.973906528517172d0

     weigp(1) =  0.066671344308682d0
     weigp(2) =  0.149451349150580d0
     weigp(3) =  0.219086362515982d0
     weigp(4) =  0.269266719309991d0
     weigp(5) =  0.295524224714752d0
     weigp(6) =  0.295524224714752d0
     weigp(7) =  0.269266719309991d0
     weigp(8) =  0.219086362515982d0
     weigp(9) =  0.149451349150580d0
     weigp(10)=  0.066671344308682d0
  else if( ngaus==11 )then
     posgp(1) = -0.978228658146100d0
     posgp(2) = -0.887062599768100d0
     posgp(3) = -0.730152005574000d0
     posgp(4) = -0.519096129206800d0
     posgp(5) = -0.269543155952300d0
     posgp(6) =  0.000000000000000d0
     posgp(7) =  0.269543155952300d0
     posgp(8) =  0.519096129206800d0
     posgp(9) =  0.730152005574000d0
     posgp(10)=  0.887062599768100d0
     posgp(11)=  0.978228658146100d0

     weigp(1) =  0.055668567116170d0
     weigp(2) =  0.125580369464900d0
     weigp(3) =  0.186290210927700d0
     weigp(4) =  0.233193764591990d0
     weigp(5) =  0.262804544510200d0
     weigp(6) =  0.272925086777900d0
     weigp(7) =  0.262804544510200d0
     weigp(8) =  0.233193764591990d0
     weigp(9) =  0.186290210927700d0
     weigp(10)=  0.125580369464900d0
     weigp(11)=  0.055668567116170d0

  else if( ngaus==12 )then           
     posgp(1) = -0.981560634246719d0
     posgp(2) = -0.904117256370475d0
     posgp(3) = -0.769902674194305d0
     posgp(4) = -0.587317954286617d0
     posgp(5) = -0.367831498918180d0
     posgp(6) = -0.125333408511469d0
     posgp(7) =  0.125333408511469d0
     posgp(8) =  0.367831498918180d0
     posgp(9) =  0.587317954286617d0
     posgp(10)=  0.769902674194305d0
     posgp(11)=  0.904117256370475d0
     posgp(12)=  0.981560634246719d0

     weigp(1) =  0.047175336386512d0
     weigp(2) =  0.106939325995318d0
     weigp(3) =  0.160078328543346d0
     weigp(4) =  0.203167426723066d0
     weigp(5) =  0.233492536538355d0
     weigp(6) =  0.249147045813403d0
     weigp(7) =  0.249147045813403d0
     weigp(8) =  0.233492536538355d0
     weigp(9) =  0.203167426723066d0
     weigp(10)=  0.160078328543346d0
     weigp(11)=  0.106939325995318d0
     weigp(12)=  0.047175336386512d0

  else if( ngaus==14 )then           
     posgp(1) = -0.986283808696812d0
     posgp(2) = -0.928434883663574d0
     posgp(3) = -0.827201315069765d0
     posgp(4) = -0.687292904811685d0
     posgp(5) = -0.515248636358154d0
     posgp(6) = -0.319112368927890d0
     posgp(7) = -0.108054948707344d0
     posgp(8) =  0.108054948707344d0
     posgp(9) =  0.319112368927890d0
     posgp(10)=  0.515248636358154d0
     posgp(11)=  0.687292904811685d0
     posgp(12)=  0.827201315069765d0
     posgp(13)=  0.928434883663574d0
     posgp(14)=  0.986283808696812d0

     weigp(1) =  0.035119460331752d0
     weigp(2) =  0.080158087159760d0
     weigp(3) =  0.121518570687903d0
     weigp(4) =  0.157203167158194d0
     weigp(5) =  0.185538397477938d0
     weigp(6) =  0.205198463721296d0
     weigp(7) =  0.215263853463158d0
     weigp(8) =  0.215263853463158d0
     weigp(9) =  0.205198463721296d0
     weigp(10)=  0.185538397477938d0
     weigp(11)=  0.157203167158194d0
     weigp(12)=  0.121518570687903d0
     weigp(13)=  0.080158087159760d0
     weigp(14)=  0.035119460331752d0

  else if( ngaus==16 )then           
     posgp(1) = -0.989400934991650d0
     posgp(2) = -0.944575023073233d0
     posgp(3) = -0.865631202387832d0
     posgp(4) = -0.755404408355003d0
     posgp(5) = -0.617876244402644d0
     posgp(6) = -0.458016777657227d0
     posgp(7) = -0.281603550779259d0
     posgp(8) = -0.095012509837637d0
     posgp(9) =  0.095012509837637d0
     posgp(10)=  0.281603550779259d0
     posgp(11)=  0.458016777657227d0
     posgp(12)=  0.617876244402644d0
     posgp(13)=  0.755404408355003d0
     posgp(14)=  0.865631202387832d0
     posgp(15)=  0.944575023073233d0
     posgp(16)=  0.989400934991650d0

     weigp(1) =  0.027152459411754d0
     weigp(2) =  0.062253523938648d0
     weigp(3) =  0.095158511682493d0
     weigp(4) =  0.124628971255534d0
     weigp(5) =  0.149595988816577d0
     weigp(6) =  0.169156519395003d0
     weigp(7) =  0.182603415044924d0
     weigp(8) =  0.189450610455069d0
     weigp(9) =  0.189450610455069d0
     weigp(10)=  0.182603415044924d0
     weigp(11)=  0.169156519395003d0
     weigp(12)=  0.149595988816577d0
     weigp(13)=  0.124628971255534d0
     weigp(14)=  0.095158511682493d0
     weigp(15)=  0.062253523938648d0
     weigp(16)=  0.027152459411754d0

  else if( ngaus==24 )then           
     posgp(1) = -0.99518721999702d0
     posgp(2) = -0.97472855597131d0
     posgp(3) = -0.93827455200273d0
     posgp(4) = -0.88641552700440d0
     posgp(5) = -0.82000198597390d0
     posgp(6) = -0.74012419157855d0
     posgp(7) = -0.64809365193697d0
     posgp(8) = -0.54542147138884d0
     posgp(9) = -0.43379350762604d0
     posgp(10)= -0.31504267969616d0
     posgp(11)= -0.19111886747361d0
     posgp(12)= -0.06405689286260d0
     posgp(13)=  0.06405689286605d0
     posgp(14)=  0.19111886747361d0
     posgp(15)=  0.31504267969616d0
     posgp(16)=  0.43379350762604d0
     posgp(17)=  0.54542147138883d0
     posgp(18)=  0.64809365193697d0
     posgp(19)=  0.74012419157855d0
     posgp(20)=  0.82000198597390d0
     posgp(21)=  0.88641552700440d0
     posgp(22)=  0.93827455200273d0
     posgp(23)=  0.97472855597130d0
     posgp(24)=  0.99518721999702d0

     weigp(1) = 0.01234122979998d0 
     weigp(2) = 0.02853138862893d0 
     weigp(3) = 0.04427743881741d0 
     weigp(4) = 0.05929858491543d0 
     weigp(5) = 0.07334648141107d0 
     weigp(6) = 0.08619016153195d0 
     weigp(7) = 0.09761865210411d0 
     weigp(8) = 0.10744427011596d0 
     weigp(9) = 0.11550566805372d0 
     weigp(10)= 0.12167047292780d0 
     weigp(11)= 0.12583745634682d0 
     weigp(12)= 0.12793819534675d0 
     weigp(13)= 0.12793819534675d0 
     weigp(14)= 0.12583745634682d0 
     weigp(15)= 0.12167047292780d0 
     weigp(16)= 0.11550566805372d0 
     weigp(17)= 0.10744427011596d0 
     weigp(18)= 0.09761865210411d0 
     weigp(19)= 0.08619016153195d0 
     weigp(20)= 0.07334648141107d0 
     weigp(21)= 0.05929858491543d0 
     weigp(22)= 0.04427743881741d0 
     weigp(23)= 0.02853138862893d0 
     weigp(24)= 0.01234122979998d0 

  endif



  return
end subroutine getgq1d




subroutine getnmat(nnode,nndof,shap, nmat)
  !=======================================================================
  !  getnmat = construct n matrix which is
  !            corresponding to dimension of current problem
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nnode : the total number of node per element
  !
  !  nndof : dof per node
  !
  !  shap(nnode) : shape function value
  !
  !  output:
  !  ------
  !  nmat(nndof,nndof*nnode) : shape function matrix
  !                             [ n1  0   n2  0   ...]  
  !                             [ 0   n1  0   n2  ...]
  !                             [ ...                ]
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nnode, nndof
  real(8), dimension(nnode), intent(in) :: shap

  real(8), dimension(nndof,nndof*nnode), intent(out) :: nmat
  ! ====================================
  ! local variable
  ! ==============
  integer :: icol

  ! loop index
  integer :: inode, idof
  ! ====================================
  ! initialize
  nmat(:,:)= 0.0d0

  ! set n matrix
  do inode=1, nnode
     do idof=1, nndof
        icol=nndof*(inode-1)+idof
        nmat(idof,icol)= shap(inode)
     end do
  end do



  return
end subroutine getnmat




subroutine getgqcod(ndime,nnode,shap,ecord, gqcod)
  !=======================================================================
  !  getgqcod= get gq coordinate in physical domain
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : dimension of problem
  !
  !  nnode : total number of node per element
  !  
  !  shap(1,nnode) : shape function value
  !
  !  ecord(ndime,nnode) : nodal coordinates
  !
  !  output:
  !  ------
  !  gqcod(ndime,1) : xy coordinate of gq point
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime, nnode
  real(8), dimension(nnode), intent(in) :: shap
  real(8), dimension(ndime,nnode), intent(in) :: ecord

  real(8), dimension(ndime,1), intent(out) :: gqcod
   ! ====================================
  ! local variable
  ! ==============
  
  ! loop index
  integer :: inode, idime
  ! ====================================

  ! initialize
  gqcod(:,:)= 0.0d0


  ! compute gq position in xy domain
  do inode=1, nnode
     do idime=1, ndime
        gqcod(idime,1)= gqcod(idime,1) + shap(inode) * ecord(idime,inode)
     end do
  end do



  return
end subroutine getgqcod





subroutine getb0mat2d(nnode,cartd,ftens, b0mat2d)
  !=======================================================================
  !  getb0mat2d = 2d voight form of B matrix for {S} (see TB's book page 202)
  !
  !               b_0= sym_{i,k} ( a n_I / a X_i * F_jk)
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nnode : the total number of node per element
  !
  !  cartd(ndime,nnode) : cartesian derivation
  !
  !  ftens(ndime,ndime) : deformation tensor
  !
  !  output:
  !  ------
  !  b0mat(ndime,nnode*ndime) : b_0= sym_{i,k} ( a N_I / a X_i * F_jk)
  !                            
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nnode
  real(8), dimension(2,nnode), intent(in) :: cartd
  real(8), dimension(2,2), intent(in) :: ftens

  real(8), dimension(3,2*nnode), intent(out) :: b0mat2d
  ! ====================================
  ! local variable
  ! ==============
  integer :: icol

  ! loop index
  integer :: inode
  ! ====================================

  ! initialize
  b0mat2d(:,:)= 0.0d0


  ! set b0 matrix for pk2 stress tensor
  do inode=1, nnode
     icol= 2 * inode - 1

     b0mat2d(1,icol)= cartd(1,inode) * ftens(1,1)
     b0mat2d(1,icol+1)= cartd(1,inode) * ftens(2,1)

     b0mat2d(2,icol)= cartd(2,inode) * ftens(1,2)
     b0mat2d(2,icol+1)= cartd(2,inode) * ftens(2,2)

     b0mat2d(3,icol)= cartd(1,inode) * ftens(1,2) + cartd(2,inode) * ftens(1,1)
     b0mat2d(3,icol+1)= cartd(1,inode) * ftens(2,2) + cartd(2,inode) * ftens(2,1)
  end do



  return
end subroutine getb0mat2d




subroutine getbmat2d(nnode,cartd, bmat2d)
  !=======================================================================
  !  getbmat2d = 2d voight form of b matrix
  !
  !              note:
  !              ----
  !              conventional b matrix for static
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nnode : the total number of node per element
  !
  !  cartd(ndime,nnode) : cartesian derivation
  !
  !  output:
  !  ------
  !  bmat2d(3,2*nnode) : conventional fe b matrix
  !                            
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nnode
  real(8), dimension(2,nnode), intent(in) :: cartd

  real(8), dimension(3,2*nnode), intent(out) :: bmat2d
  ! ====================================
  ! local variable
  ! ==============
  integer :: icol

  ! loop index
  integer :: inode
  ! ====================================

  ! initialize
  bmat2d(:,:)= 0.0d0


  ! set conventional fe b matrix
  do inode=1, nnode
     icol= 2 * inode - 1

     bmat2d(1,icol)= cartd(1,inode)
     bmat2d(1,icol+1)= 0.0d0

     bmat2d(2,icol)= 0.0d0
     bmat2d(2,icol+1)= cartd(2,inode)

     bmat2d(3,icol)= cartd(2,inode)
     bmat2d(3,icol+1)= cartd(1,inode)
  end do



  return
end subroutine getbmat2d





subroutine getb0mat3d(nnode,cartd,ftens, b0mat3d)
  !=======================================================================
  !  getb0mat3d = 3d voight form of B matrix for {S} (see TB's book page 203)
  !
  !               b_0= sym_{i,k} ( a n_I / a X_i * F_jk)
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nnode : the total number of node per element
  !
  !  cartd(ndime,nnode) : cartesian derivation
  !
  !  ftens(ndime,ndime) : deformation tensor
  !
  !  output:
  !  ------
  !  b0mat3d(ndime,nnode*ndime) : b_0= sym_{i,k} ( a N_I / a X_i * F_jk)
  !                            
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nnode
  real(8), dimension(3,nnode), intent(in) :: cartd
  real(8), dimension(3,3), intent(in) :: ftens

  real(8), dimension(6,3*nnode), intent(out) :: b0mat3d
  ! ====================================
  ! local variable
  ! ==============
  integer :: icol

  ! loop index
  integer :: inode
  ! ====================================

  ! initialize
  b0mat3d(:,:)= 0.0d0

  ! set b0 matrix for pk2 stress tensor
  do inode=1, nnode
     icol= 3 * inode - 2

     b0mat3d(1,icol)= cartd(1,inode) * ftens(1,1)
     b0mat3d(1,icol+1)= cartd(1,inode) * ftens(2,1)
     b0mat3d(1,icol+2)= cartd(1,inode) * ftens(3,1)

     b0mat3d(2,icol)= cartd(2,inode) * ftens(1,2)
     b0mat3d(2,icol+1)= cartd(2,inode) * ftens(2,2)
     b0mat3d(2,icol+2)= cartd(2,inode) * ftens(3,2)

     b0mat3d(3,icol)= cartd(3,inode) * ftens(1,3)
     b0mat3d(3,icol+1)= cartd(3,inode) * ftens(2,3)
     b0mat3d(3,icol+2)= cartd(3,inode) * ftens(3,3)

     b0mat3d(4,icol)= cartd(2,inode) * ftens(1,3) + cartd(3,inode) * ftens(1,2)
     b0mat3d(4,icol+1)= cartd(2,inode) * ftens(2,3) + cartd(3,inode) * ftens(2,2)
     b0mat3d(4,icol+2)= cartd(2,inode) * ftens(3,3) + cartd(3,inode) * ftens(3,2)

     b0mat3d(5,icol)= cartd(1,inode) * ftens(1,3) + cartd(3,inode) * ftens(1,1)
     b0mat3d(5,icol+1)= cartd(1,inode) * ftens(2,3) + cartd(3,inode) * ftens(2,1)
     b0mat3d(5,icol+2)= cartd(1,inode) * ftens(3,3) + cartd(3,inode) * ftens(3,1)

     b0mat3d(6,icol)= cartd(1,inode) * ftens(1,2) + cartd(2,inode) * ftens(1,1)
     b0mat3d(6,icol+1)= cartd(1,inode) * ftens(2,2) + cartd(2,inode) * ftens(2,1)
     b0mat3d(6,icol+2)= cartd(1,inode) * ftens(3,2) + cartd(2,inode) * ftens(3,1)
  end do



  return
end subroutine getb0mat3d




subroutine rotprojbt1(ecord,edisp,evelo, edisp0,evelo0)
  !=======================================================================
  !  rotprojbt1 = rotation projection: project 6 dof to 5 dof: nodal solution
  !
  !  arguments description
  !  ---------------------
  !  input
  !  -----
  !  ecord(3,4) : local element nodal coordinate
  !
  !  edisp(6,4) : d_x, d_y, d_z, r_x, r_y, r_z
  !
  !  evelo(6,4) : v_x, v_y, v_z, rv_x, rv_y, rv_z
  !
  !  output:
  !  ------
  !  edisp0(5,4) : d_x, d_y, d_z, r_x, r_y
  !
  !  evelo0(5,4) : v_x, v_y, v_z, rv_x, rv_y
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(3,4), intent(in) :: ecord
  real(8), dimension(6,4), intent(in) :: edisp
  real(8), dimension(6,4), intent(in) :: evelo

  real(8), dimension(5,4), intent(out) :: edisp0
  real(8), dimension(5,4), intent(out) :: evelo0
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(3,4) :: ecurn
  real(8), dimension(2,3) :: rotpmat

  real(8), dimension(3,1) :: erot, erotdot
  real(8), dimension(2,1) :: erot0, erotdot0

  ! loop index
  integer :: inode
  ! ====================================

  ! initialize
  edisp0(:,:)= 0.0d0
  evelo0(:,:)= 0.0d0

  ! get current nodal coordinate
  ecurn(1:3,1:4)= ecord(1:3,1:4) + edisp(1:3,1:4)

  ! loop over nodes
  do inode=1, 4

     ! -----------------------------------------------------
     ! set translational dof
     edisp0(1:3,inode)= edisp(1:3,inode)
     evelo0(1:3,inode)= evelo(1:3,inode)

     ! -----------------------------------------------------
     ! get rotation projection matrix
     call getrotpmatbt(inode,ecurn, rotpmat)
        ! input : inode,ecurn
        ! output : rotpmat

     ! -----------------------------------------------------
     ! extract rotation
     erot(1:3,1)= edisp(4:6,inode)

     ! rotation projection: rot0= rotpmat.erot
     call matprd(2,3,0, 3,1,0, 2,1, rotpmat,erot, erot0)
        ! input : 2,3,0, 3,1,0, 2,1, rotpmat,erot
        ! output : erot0

     ! set results
     edisp0(4:5,inode)= erot0(1:2,1)

     ! -----------------------------------------------------
     ! extract angular velocity
     erotdot(1:3,1)= evelo(4:6,inode)

     ! rotation projection: rot0= rotpmat.erot
     call matprd(2,3,0, 3,1,0, 2,1, rotpmat,erotdot, erotdot0)
        ! input : 2,3,0, 3,1,0, 2,1, rotpmat,erotdot
        ! output : erotdot0

     ! set results
     evelo0(4:5,inode)= erotdot0(1:2,1)

     ! -----------------------------------------------------

  end do




  return
end subroutine rotprojbt1





subroutine rotprojbt2(ecord,edisp,evelo,eaccl, edisp0,evelo0,eaccl0)
  !=======================================================================
  !  rotprojbt2 = rotation projection: project 6 dof to 5 dof: nodal solution
  !
  !  arguments description
  !  ---------------------
  !  input
  !  -----
  !  ecord(3,4) : local element nodal coordinate
  !
  !  edisp(6,4) : d_x, d_y, d_z, r_x, r_y, r_z
  !
  !  evelo(6,4) : v_x, v_y, v_z, rv_x, rv_y, rv_z
  !
  !  eaccl(6,4) : a_x, a_y, a_z, ra_x, ra_y, ra_z
  !
  !  output:
  !  ------
  !  edisp0(5,4) : d_x, d_y, d_z, r_x, r_y
  !
  !  evelo0(5,4) : v_x, v_y, v_z, rv_x, rv_y
  !
  !  eaccl0(5,4) : a_x, a_y, a_z, ra_x, ra_y
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(3,4), intent(in) :: ecord
  real(8), dimension(6,4), intent(in) :: edisp
  real(8), dimension(6,4), intent(in) :: evelo
  real(8), dimension(6,4), intent(in) :: eaccl

  real(8), dimension(5,4), intent(out) :: edisp0
  real(8), dimension(5,4), intent(out) :: evelo0
  real(8), dimension(5,4), intent(out) :: eaccl0
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(3,4) :: ecurn
  real(8), dimension(2,3) :: rotpmat

  real(8), dimension(3,1) :: erot, erotdot
  real(8), dimension(2,1) :: erot0, erotdot0

  ! loop index
  integer :: inode
  ! ====================================

  ! initialize
  edisp0(:,:)= 0.0d0
  evelo0(:,:)= 0.0d0
  eaccl0(:,:)= 0.0d0

  ! get current nodal coordinate
  ecurn(1:3,1:4)= ecord(1:3,1:4) + edisp(1:3,1:4)

  ! loop over nodes
  do inode=1, 4

     ! -----------------------------------------------------
     ! set translational dof
     edisp0(1:3,inode)= edisp(1:3,inode)
     evelo0(1:3,inode)= evelo(1:3,inode)
     eaccl0(1:3,inode)= eaccl(1:3,inode)

     ! -----------------------------------------------------
     ! get rotation projection matrix
     call getrotpmatbt(inode,ecurn, rotpmat)
        ! input : inode,ecurn
        ! output : rotpmat

     ! -----------------------------------------------------
     ! extract rotation
     erot(1:3,1)= edisp(4:6,inode)

     ! rotation projection: rot0= rotpmat.erot
     call matprd(2,3,0, 3,1,0, 2,1, rotpmat,erot, erot0)
        ! input : 2,3,0, 3,1,0, 2,1, rotpmat,erot
        ! output : erot0

     ! set results
     edisp0(4:5,inode)= erot0(1:2,1)

     ! -----------------------------------------------------
     ! extract angular velocity
     erotdot(1:3,1)= evelo(4:6,inode)

     ! rotation projection: rot0= rotpmat.erot
     call matprd(2,3,0, 3,1,0, 2,1, rotpmat,erotdot, erotdot0)
        ! input : 2,3,0, 3,1,0, 2,1, rotpmat,erotdot
        ! output : erotdot0

     ! set results
     evelo0(4:5,inode)= erotdot0(1:2,1)

     ! -----------------------------------------------------
     ! extract angular acceleration
     erotdot(1:3,1)= eaccl(4:6,inode)

     ! rotation projection: rot0= rotpmat.erot
     call matprd(2,3,0, 3,1,0, 2,1, rotpmat,erotdot, erotdot0)
        ! input : 2,3,0, 3,1,0, 2,1, rotpmat,erotdot
        ! output : erotdot0

     ! set results
     eaccl0(4:5,inode)= erotdot0(1:2,1)

     ! -----------------------------------------------------

  end do



  return
end subroutine rotprojbt2





subroutine rotprojbt3(ecord,edisp,efvec0, efvec)
  !=======================================================================
  !  rotprojbt3 = rotation projection: project 5 dof to 6 dof : force vector
  !
  !
  !  arguments description
  !  ---------------------
  !  input
  !  -----
  !  ecord(3,4) : local element nodal coordinate
  !
  !  edisp(6,4) : d_x, d_y, d_z, r_x, r_y, r_z
  !
  !  efvec0(20,1) : force vector f_x, f_y, f_z, m_x, m_y
  !
  !  output:
  !  ------
  !  efvec(24,1) : force vector f_x, f_y, f_z, m_x, m_y, m_z
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(3,4), intent(in) :: ecord
  real(8), dimension(6,4), intent(in) :: edisp
  real(8), dimension(20,1), intent(in) :: efvec0

  real(8), dimension(24,1), intent(out) :: efvec
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(3,4) :: ecurn
  real(8), dimension(2,3) :: rotpmat

  real(8), dimension(3,1) :: emomnt
  real(8), dimension(2,1) :: emomnt0

  integer :: iloc0, iloc
  ! loop index
  integer :: inode
  ! ====================================

  ! initialize
  efvec(:,:)= 0.0d0


  ! get current nodal coordinate
  ecurn(1:3,1:4)= ecord(1:3,1:4) + edisp(1:3,1:4)


  ! loop over nodes
  do inode=1, 4

     ! address
     iloc0= 5*(inode-1)
     iloc= 6*(inode-1)

     ! -----------------------------------------------------
     ! force
     efvec(iloc+1,1)= efvec0(iloc0+1,1) ! f_x
     efvec(iloc+2,1)= efvec0(iloc0+2,1) ! f_y
     efvec(iloc+3,1)= efvec0(iloc0+3,1) ! f_z

     ! -----------------------------------------------------
     ! moment
     emomnt0(1,1)= efvec0(iloc0+4,1) ! m_x
     emomnt0(2,1)= efvec0(iloc0+5,1) ! m_y

     ! get rotation projection matrix
     call getrotpmatbt(inode,ecurn, rotpmat)
        ! input : inode,ecurn
        ! output : rotpmat

     ! rotation projection: rot0= rotpmat.erot
     call matprd(2,3,1, 2,1,0, 3,1, rotpmat,emomnt0, emomnt)
        ! input : 2,3,1, 2,1,0, 3,1, rotpmat,emomnt0
        ! output : emomnt

     ! set results
     efvec(iloc+4,1)= emomnt(1,1) ! m_x
     efvec(iloc+5,1)= emomnt(2,1) ! m_y
     efvec(iloc+6,1)= emomnt(3,1) ! m_z
     ! -----------------------------------------------------

  end do



  return
end subroutine rotprojbt3





subroutine getrotpmatbt(inode,ecurn, rotpmat)
  !=======================================================================
  !  getrotpmatbt = compute rotation projection matrix
  !
  !                note:
  !                ----
  !                theta^bar = theta - ( theta . e3 ) e3
  !
  !  arguments description
  !  ---------------------
  !  input
  !  -----
  !  inode : current node number
  !
  !  ecurn(3,4) : current element nodal coordinate
  !
  !  output:
  !  ------
  !  rotpmat(2,3) : rotation projection matrix
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: inode
  real(8), dimension(3,4), intent(in) :: ecurn

  real(8), dimension(2,3), intent(out) :: rotpmat
  ! ====================================
  ! local variable
  ! ==============
  integer, dimension(3,4) :: ijkndx

  real(8), dimension(3,1) :: nodi, nodj, nodk
  real(8), dimension(3,1) ::e1vec, e2vec, e3vec
  ! ====================================

  ! initialize
  rotpmat(:,:)= 0.d0

  data ijkndx /1,2,4, 2,3,1, 3,4,2, 4,1,3/

  !   k
  !    o
  !    |
  !    |
  !    |
  ! e2 ^
  !    |
  !    |
  !    |
  !    o-------->--------o
  !   i        e1         j
  
  ! define node
  nodi(1:3,1)= ecurn(1:3,ijkndx(1,inode))
  nodj(1:3,1)= ecurn(1:3,ijkndx(2,inode))
  nodk(1:3,1)= ecurn(1:3,ijkndx(3,inode))

  ! element edge vector: e1, e2
  e1vec(1:3,1)= nodj(1:3,1) - nodi(1:3,1)
  e2vec(1:3,1)= nodk(1:3,1) - nodi(1:3,1)

  ! normalize 
  call unitvec1(3, e1vec)
     ! input : 3(ndime)
     ! inoutput : e1vec

  call unitvec1(3, e2vec)
     ! input : 3(ndime)
     ! inoutput : e2vec

  ! e3 vector: e3= e1 x e2 / ||e1 x e2|| 
  ! ---------
  call crsprdt3d(1,e1vec,e2vec, e3vec)
    ! input : 1(opt:normalize),e1vec,e2vec
    ! output : e3vec


  ! set rotation projection matrix
  rotpmat(1,1)= 1.0d0 - e3vec(1,1)*e3vec(1,1)
  rotpmat(1,2)= - e3vec(1,1)*e3vec(2,1)
  rotpmat(1,3)= - e3vec(1,1)*e3vec(3,1)

  rotpmat(2,1)= - e3vec(2,1)*e3vec(1,1)
  rotpmat(2,2)= 1.0d0 - e3vec(2,1)*e3vec(2,1)
  rotpmat(2,3)= - e3vec(2,1)*e3vec(3,1)



  return
end subroutine getrotpmatbt





subroutine getenergy(msize,nsize,delt,velo,fvec,fvecold, wenrg)
  !=======================================================================
  !  getenergy = compute energy
  !
  !              note:
  !              ----
  !              w^n+1 = w^n + delt*velo*0.5*(fvec+fvecold)
  !
  !  arguments description
  !  ---------------------
  !  input
  !  -----
  !
  !  inoutput:
  !  --------
  !  wenrg : energy
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: msize, nsize
  real(8), intent(in) :: delt
  real(8), dimension(msize,1), intent(in) :: velo
  real(8), dimension(msize,1), intent(in) :: fvec
  real(8), dimension(msize,1), intent(in) :: fvecold
  ! ------------------------------------

  real(8), intent(inout) :: wenrg
  ! ====================================
  ! local variable
  ! ==============

  ! loop index
  integer :: indx
  ! ====================================

  ! initialize: do not initialzie
  
  
  do indx=1, nsize

     wenrg= wenrg + delt* velo(indx,1) * 0.50d0 * ( fvec(indx,1) + fvecold(indx,1) )

  end do



  return
end subroutine getenergy





subroutine getsq1d(nsimp, possp,weisp)
  !=======================================================================
  !  getsq1d = sets up the simpson integration constants for 1d line
  !
  !  input :
  !  -------
  !  nsimp : the number of simpson quadarture
  !
  !  output :
  !  --------
  !  possp(nsimp) : simpson quadrature position in psi domain
  !
  !  weisp(nsimp) : simpson quadrature weight values
  !
  !=======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nsimp

  real(8), dimension(nsimp), intent(out) :: possp
  real(8), dimension(nsimp), intent(out) :: weisp
  ! ====================================
  ! local variable
  ! ==============

  ! loop index
  integer :: isimp
  ! ====================================

  ! initialize
  possp(:)=0.d0
  weisp(:)=0.d0 


  ! set 1d mension simpson quadrature and it's weight value
  do isimp=1, nsimp
     
     possp(isimp)= -1.0d0 + 1.0d0 / real(nsimp) + 2.0d0 / real(nsimp) * real( isimp - 1 )
     
     weisp(isimp)= 2.0d0 / real(nsimp)

  end do


  return
end subroutine getsq1d





subroutine getsq2dquad(nsqpt,msimp, sqpoin,sqweigt)
  !=======================================================================
  !  getsq2dquad = get 2d sq rule for quad
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nsqpt : 1dsq rule
  !
  !  msimp : total number of sq
  !           msimp = nsqpt**2 : psi and eta direction
  !
  !  output:
  !  ------
  !  gqpoin(2,msimp) : gq position
  !
  !  gqweigt(msimp) : gq weight
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nsqpt, msimp

  real(8), dimension(2,msimp), intent(out) :: sqpoin 
  real(8), dimension(msimp), intent(out) :: sqweigt 
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(nsqpt) :: possp
  real(8), dimension(nsqpt) :: weisp

  integer :: nsimp

  ! loop index
  integer :: isimp, jsimp
  ! ====================================

  ! initialize
  sqpoin(:,:)= 0.0d0
  sqweigt(:)= 0.0d0
  

  ! get simpson qquadrature rule for 1d line domain
  call getsq1d(nsqpt, possp,weisp)
     ! iput : nsqpt
     ! output : possp,weisp

  ! set total gq point and it's weight
  nsimp=0 ! initialize
  do isimp=1, nsqpt
     do jsimp=1, nsqpt

        ! increase counter
        nsimp= nsimp+1

        ! set 2d position
        sqpoin(1,nsimp)= possp(isimp)
        sqpoin(2,nsimp)= possp(jsimp)

        ! set weight
       sqweigt(nsimp)= weisp(isimp) * weisp(jsimp)

     end do ! jgaus
  end do ! igaus

  ! error check
  if( nsimp /= msimp) then
     write(*,*) "inconsistency in sq rule: getsq2dquad"
     write(nout5,*) "inconsistency in sq rule: getsq2dquad"
     stop

  end if



  return
end subroutine getsq2dquad





subroutine getsqele(optele,nsqdim,nsqpt,msimp, sqpoin,sqweigt)
  !=======================================================================
  !  getsqele = get rule for given sq dimension and element
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optele : element type option handler
  !
  !  nsqdim : dimension of sq
  !
  !  nsqpt : 1d sq rule
  !
  !  msimp : total number of sq
  !           for quad. ele. : msimp= nsqpt**2 :psi and eta direction
  !
  !  output:
  !  ------
  !  sqpoin(2,msimp) : sq position
  !
  !  sqweigt(msimp) : sq weight
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optele
  integer, intent(in) :: nsqdim
  integer, intent(in) :: nsqpt, msimp

  real(8), dimension(nsqdim,msimp), intent(out) :: sqpoin 
  real(8), dimension(msimp), intent(out) :: sqweigt 
  ! ====================================
  ! local variable
  ! ==============

  ! loop index
  ! ====================================

  ! initialize
  sqpoin(:,:)= 0.0d0
  sqweigt(:)= 0.0d0
  

  if ( nsqdim== 2 .and. optele == 2 ) then ! 2d quad

     ! get sq rule for 2d quad
     call getsq2dquad(nsqpt,msimp, sqpoin,sqweigt)
        ! iput : nsqpt,msimp
        ! output : sqpoin,sqweigt

  else
     write(*,*) "not available: getsqele"
     write(nout5,*) "not available: getsqele"
     stop

  end if



  return
end subroutine getsqele





subroutine getbmat1pt(ecordloc, area, bmat1pt)
  !=======================================================================
  !  getbmat1pt = compute b matrix of 4 node quad integrated with 1 point gq
  !
  !               note:
  !               ----
  !               belytschko, wong and chiang, CMAME, 1992, vol. 96, pp. 93-107
  !               advances in one point quadrature shell elements
  !               (see, eq (19))
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ecordloc(3,4) : local element nodal coordinate
  !
  !  output:
  !  ------
  !  bmat1pt(2,4) : b matrix of 4 node quad integrated with 1 point gq
  !                            
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(3,4), intent(in) :: ecordloc
  real(8), intent(in) :: area

  real(8), dimension(2,4), intent(out) :: bmat1pt
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: const
  real(8) :: x1,x2,x3,x4
  real(8) :: y1,y2,y3,y4

  ! ====================================

  ! define components
  x1= ecordloc(1,1)
  x2= ecordloc(1,2)
  x3= ecordloc(1,3)
  x4= ecordloc(1,4)

  y1= ecordloc(2,1)
  y2= ecordloc(2,2)
  y3= ecordloc(2,3)
  y4= ecordloc(2,4)

  ! compute constant
  const= 1.0d0 /( 2.0d0*area)

  ! set b_x components
  bmat1pt(1,1)= const*( y2-y4 ) 
  bmat1pt(1,2)= const*( y3-y1 ) 
  bmat1pt(1,3)= const*( y4-y2 ) 
  bmat1pt(1,4)= const*( y1-y3 ) 

  ! set b_y components
  bmat1pt(2,1)= const*( x4-x2 ) 
  bmat1pt(2,2)= const*( x1-x3 ) 
  bmat1pt(2,3)= const*( x2-x4 ) 
  bmat1pt(2,4)= const*( x3-x1 ) 



  return
end subroutine getbmat1pt



subroutine getbcmat1pt(ecordloc, area, gamma, zgamma, bcmat1pt)
  !=======================================================================
  !  getbcmat1pt = compute b^c matrix for warping correction
  !                with one point integration and local z method
  !                for computational efficiency, explicit components form is used
  !
  !                note:
  !                ----
  !                belytschko, wong and chiang, CMAME, 1992, vol. 96, pp. 93-107
  !                advances in one point quadrature shell elements
  !                (see, eq (30))
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ecordloc(3,4) : local element nodal coordinate
  !
  !  output:
  !  ------
  !  bcmat1pt(2,4) : bmatrix for warping correction
  !                            
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(3,4), intent(in) :: ecordloc
  real(8), intent(in) :: area
  real(8), dimension(4,1), intent(in) :: gamma
  real(8), intent(in) :: zgamma

  real(8), dimension(2,4), intent(out) :: bcmat1pt
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: const
  real(8) :: x1,x2,x3,x4
  real(8) :: y1,y2,y3,y4

  ! loop index
  integer :: inode
  ! ====================================

  ! define components
  x1= ecordloc(1,1)
  x2= ecordloc(1,2)
  x3= ecordloc(1,3)
  x4= ecordloc(1,4)

  y1= ecordloc(2,1)
  y2= ecordloc(2,2)
  y3= ecordloc(2,3)
  y4= ecordloc(2,4)


  ! compute constant
  const= (2.0d0*zgamma)/area**2

  ! set b_x ^c components
  bcmat1pt(1,1)= const*( x1-x3 ) 
  bcmat1pt(1,2)= const*( x4-x2 ) 
  bcmat1pt(1,3)= const*( x3-x1 ) 
  bcmat1pt(1,4)= const*( x2-x4 ) 

  ! set b_y ^c components
  bcmat1pt(2,1)= const*( y1-y3 ) 
  bcmat1pt(2,2)= const*( y4-y2 ) 
  bcmat1pt(2,3)= const*( y3-y1 ) 
  bcmat1pt(2,4)= const*( y2-y4 ) 



  return
end subroutine getbcmat1pt




subroutine getbsmat1pt(ecordloc, bsmat1pt)
  !=======================================================================
  !  getbsmat1pt = compute b^s matrix for shear projection
  !                with one point integration
  !                for computational efficiency, explicit components form is used
  !
  !                note:
  !                ----
  !                belytschko, wong and chiang, CMAME, 1992, vol. 96, pp. 93-107
  !                advances in one point quadrature shell elements
  !                (see, eq(35))
  !
  !  arguments description
  !  ---------------------
  !  input
  !  -----
  !  ecordloc(3,4) : local element nodal coordinate
  !
  !  output:
  !  ------
  !  bsmat1pt(2,3,4) : b^s matrix for shear projection
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(3,4), intent(in) :: ecordloc

  real(8), dimension(2,3,4), intent(out) :: bsmat1pt
  ! ====================================
  ! local variable
  ! ==============
  integer, dimension(3,4) :: ijkndx

  real(8) :: xi, xj, xk
  real(8) :: yi, yj, yk
  real(8) :: xji, yji, xik, yik
  real(8) :: lji, lik

  ! loop index
  integer :: inode
  ! ====================================

  ! initialize
  data ijkndx /1,2,4, 2,3,1, 3,4,2, 4,1,3/

  
  do inode=1, 4

     ! extract current and neighbor nodal coordinate
     xi= ecordloc(1,ijkndx(1,inode))
     xj= ecordloc(1,ijkndx(2,inode))
     xk= ecordloc(1,ijkndx(3,inode))

     yi= ecordloc(2,ijkndx(1,inode))
     yj= ecordloc(2,ijkndx(2,inode))
     yk= ecordloc(2,ijkndx(3,inode))

     ! compute relative nodal coordinate difference
     xji= xj-xi
     yji= yj-yi

     xik= xi-xk
     yik= yi-yk

     ! compute edge length: l_ji, l_ik
     lji= dsqrt( xji**2 + yji**2 )
     lik= dsqrt( xik**2 + yik**2 )


     ! set b_x ^s components
     bsmat1pt(1,1,inode)= (1.0d0/2.0d0)*( xji/(lji**2) - xik/(lik**2) )
     bsmat1pt(1,2,inode)= (1.0d0/4.0d0)*( xji*yji/(lji**2) + xik*yik/(lik**2) ) 
     bsmat1pt(1,3,inode)= -(1.0d0/4.0d0)*( (xji/lji)**2 + (xik/lik)**2 ) 

     ! set b_y ^s components
     bsmat1pt(2,1,inode)= (1.0d0/2.0d0)*( yji/(lji**2) - yik/(lik**2) )
     bsmat1pt(2,2,inode)= (1.0d0/4.0d0)*( (yji/lji)**2 + (yik/lik)**2 ) 
     bsmat1pt(2,3,inode)= -(1.0d0/4.0d0)*( xji*yji/(lji**2) + xik*yik/(lik**2) ) 

  end do



  return
end subroutine getbsmat1pt
