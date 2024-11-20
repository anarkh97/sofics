! =======================================
! basic useful functions
! =======================================
!      type                  name                              argument
!      ----                  ----                              --------
! 2.  subroutine             matprd          (nar,nac,itrna,nbr,nbc,itrnb,ncr,ncc,mata,matb,matc)
! 3.  real(8) function       area            (nvrtx,vrtx)
! 4.  logical function       chkcross        (ndime,line1,line2)
! 5.  logical function       ptinarea        (nvrtx,vrtx,point)
! 6.  real(8) function       permu           (indx,jndx,kndx)
! 7.  real(8) function       distance        (ndime, line, point)
! 8.  real(8) function       hsign           (ndime, line, point)
! 9.  real(8) function       heavi           (xval)
! 10. real(8) function       angle           (ndime,point)
! 12. subroutine             getparent2d     (optele,nnode,ecord,poinx,poiny, s,t)
! 13. subroutine             weightff        (opttyp,optpol,rmax,dist, weight)
! 14. subroutine             angfilter       (ndime,crpt, angle)
! 15. subroutine             ludcmp          (n, a, indx,d)
! 16. subroutine             lubksb          (n,a,indx, b)
! 17. real(8) function       getdet          (ndime,amat)
! 18. subroutine             getinv          (ndime,amat, bmat)
! 19. subroutine             getpertub       (nptb,minval,maxval, ptbval)
! 20. subroutine             gettvec         (ndime,line, tvec)
! 21. subroutine             hpsort          (n, ra, indx)
! 22. subroutine             rotpt2d         (ndime,theta, poin)
! 23. subroutine             getlineseg      (ndime,seglen,segang,segorg, segline)
! 24. subroutine             crsprdt3d       (opt,avec,bvec, cvec)
! 25. real(8) function       dotprdt         (ndime,avec,bvec)
! 26. subroutine             unitvec1        (ndime, avec)
! 27. subroutine             unitvec2        (ndime,avec, bvec)
! 28. subroutine             mlsfit2dqd      (npt,rmax,ptcord,ptdist,ptval,chkpt, dvaldx,dvaldy,fitcoef)
! 29. subroutine             mlsfit2dln      (npt,rmax,ptcord,ptdist,ptval, dvaldx,dvaldy,fitcoef)
! 30. subroutine             mlsrhsmat       (npt,npol,pmat,wmat,ptval, rhsmat)
! 31. subroutine             mlsmmat         (npt,npol,pmat,wmat, mmat)
! 32. subroutine             elearea         (optele,ndime,ecord, earea)
! 33. subroutine             elehleng0       (optele,ndime,ecord, ehleng)
!
! =========================================================================================================



subroutine matprd(nar,nac,itrna,nbr,nbc,itrnb,ncr,ncc,mata,matb,matc)
  !=======================================================================
  !  matprd= matrix multiplication
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nar, nac : size of A matrix : nar by nac
  !  itrna : option handler to check transpose
  !          itrna=0: w/o transpose
  !          itrna=1: consider transposed A matrix
  !
  !  nbr, nbc : size of B matrix : nbr by nbc
  !  itrnb : option handler to check transpose
  !          itrnb=0: w/o transpose
  !          itrnb=1: consider transposed B matrix
  !
  !  ncr,ncc : size of result C matrix
  !
  !  mata, matb : A and B matrix
  !
  !  output:
  !  ------
  !  matc : result C matrix
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nar,nac ! size of a matrix
  integer, intent(in) :: itrna
  integer, intent(in) :: nbr,nbc ! size of b matrix
  integer, intent(in) :: itrnb
  integer, intent(in) :: ncr,ncc ! size of c matrix
  real(8), dimension(nar,nac), intent(in) :: mata
  real(8), dimension(nbr,nbc), intent(in) :: matb

  real(8), dimension(ncr,ncc), intent(out) :: matc
  ! ====================================
  ! local variable
  ! ==============
  integer :: index, jndex, kndex
  ! ====================================

  ! initialize
! PJSA
!  matc(:,:)=0.0d0

  ! c=a.b
  ! -----
  if((itrna.eq.0).and.(itrnb.eq.0)) then
!     ! check compatibility
!     if(nac/=nbr) then
!        write(*,*) "incompatible matrix size: matprd"
!        write(nout5,*) "incompatible matrix size: matprd"
!        stop
!     end if
!
!     ! matrix mulitiplication: c_ij=a_ik.b_kj
!     do index=1, nar
!        do jndex=1, nbc
!           do kndex=1, nac
!              matc(index,jndex)=matc(index,jndex)+mata(index,kndex)*matb(kndex,jndex)
!           end do
!        end do
!     end do
!
! PJSA
    if(nbc.eq.1) then
      call dgemv('n',nar,nac,1.0d0,mata,nar,matb,nbc,0.0d0,matc,ncc)
    else
      call dgemm('n','n',nar,nbc,nac,1.0d0,mata,nar,matb,nbr,0.0d0,matc,ncr)
    end if

  end if

  ! c=aT.b
  ! -------
  if((itrna.eq.1).and.(itrnb.eq.0)) then
!     ! check compatibility
!     if(nar/=nbr) then
!        write(*,*) "incompatible matrix size: matprd"
!        write(nout5,*) "incompatible matrix size: matprd"
!        stop
!     end if
!
!     ! matrix mulitiplication: c_ij=a_ki.b_kj
!     do index=1, nac
!        do jndex=1, nbc
!           do kndex=1, nar
!              matc(index,jndex)=matc(index,jndex)+mata(kndex,index)*matb(kndex,jndex)
!           end do
!        end do
!     end do
! PJSA
    if(nbc.eq.1) then
      call dgemv('t',nar,nac,1.0d0,mata,nar,matb,nbc,0.0d0,matc,ncc)
    else
      call dgemm('t','n',nar,nbc,nac,1.0d0,mata,nar,matb,nbr,0.0d0,matc,ncr)
    endif

  end if

  ! c=a.bT
  ! -------
  if((itrna.eq.0).and.(itrnb.eq.1)) then
!     ! check compatibility
!     if(nac/=nbc) then
!        write(*,*) "incompatible matrix size: matprd"
!        write(nout5,*) "incompatible matrix size: matprd"
!        stop
!     end if
!
!     ! matrix mulitiplication: c_ij=a_ik.b_jk
!     do index=1, nar
!        do jndex=1, nbr
!           do kndex=1, nac
!              matc(index,jndex)=matc(index,jndex)+mata(index,kndex)*matb(jndex,kndex)
!           end do
!        end do
!     end do
! PJSA
    call dgemm('n','t',nar,nbc,nac,1.0d0,mata,nar,matb,nbr,0.0d0,matc,ncr)

  end if



  return
end subroutine matprd



real(8) function area(nvrtx,vrtx)
  !=======================================================================
  !  area= compute 2d area of given polygon
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nvrtx : the total number of vertex point
  !
  !  vrtx(2,nvrtx) : vertex coordinate of polygon
  !
  !  output:
  !  ------
  !  area : 2d area of given polygon
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nvrtx
  real(8), dimension(2,nvrtx), intent(in) :: vrtx
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: centerx, centery
  real(8) :: x1,x2,x3,y1,y2,y3

  ! loop index
  integer :: ivrtx, jvrtx
  ! ====================================

  ! initialize
  area= 0.0d0

  ! find center point which is located inside of polygon
  centerx= 0.0d0 ! initialize
  centery= 0.0d0

  do ivrtx=1, nvrtx
     centerx= centerx + vrtx(1,ivrtx) / real(nvrtx)
     centery= centery + vrtx(2,ivrtx) / real(nvrtx)
  end do

  ! compute area of triangle(vx_i,vy_i)-(vx_j,vy_j)-(centerx,centery)
  do ivrtx=1, nvrtx

     if(ivrtx==nvrtx) then
        jvrtx=1
     else
        jvrtx=ivrtx+1
     end if

     x1=vrtx(1,ivrtx) ! current vertex point
     y1=vrtx(2,ivrtx)
     x2=vrtx(1,jvrtx) ! next vertex point
     y2=vrtx(2,jvrtx)
     x3=centerx ! center of polygon
     y3=centery
     area= area + 0.50d0*abs(x1*y2+y1*x3+y3*x2-y2*x3-y1*x2-x1*y3)
  end do



  return
end function area





logical function chkcross(ndime,line1,line2)
  !=======================================================================
  !  chkcross= checking whether two finite length lines are cross or not
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  line1(ndime,2) : start and end point of finite length line
  !
  !  line2(ndime,2) : start and end point of finite length line
  !
  !  output:
  !  ------
  !  chkcross : logical value
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime
  real(8), dimension(ndime,2), intent(in) :: line1, line2
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: ax,ay,bx,by,cx,cy,dx,dy
  real(8) :: det1,det2

  ! loop index
  ! ====================================

  ! initialize
  chkcross=.false.

  select case(ndime)
  case(2)
     ! set useful symbolic value, see idea note
     ax=line1(1,1)
     ay=line1(2,1)

     bx=line1(1,2)
     by=line1(2,2)

     cx=line2(1,1)
     cy=line2(2,1)

     dx=line2(1,2)
     dy=line2(2,2)

     det1= (ax-cx)*(by-cy) - (ay-cy)*(bx-cx)
     det2= (ax-dx)*(by-dy) - (ay-dy)*(bx-dx)

     
     if(det1*det2 < 0.0d0) then
        chkcross=.true.

     else if(det1*det2 >= 0.0d0) then ! note: det1*det2=0 means any one line tip is on the other line
        chkcross=.false.

     end if

  end select


  return
end function chkcross





logical function ptinarea(nvrtx,vrtx,point)
  !=======================================================================
  !  ptinarea = check location of given point in 2d
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nvrtx : the total number of vertex point
  !
  !  vrtx(2,nvrtx) : vertex coordinate of polygon
  !
  !  point(2,1) : given point
  !
  !  inout:
  !  ------
  !  ptinarea : logical flag
  !             ptinarea= .true. : given point located in polygon
  !             ptinarea= .false. : given point located out of polygon
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nvrtx
  real(8), dimension(2,nvrtx), intent(in) :: vrtx
  real(8), dimension(2,1), intent(in) :: point
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: area, area0, area1
  real(8) :: x1,x2,x3,y1,y2,y3

  ! loop index
  integer :: ivrtx, jvrtx
  ! ====================================

  ! initialize
  ptinarea= .false.

  ! compute area of given polyogon
  area0= area(nvrtx,vrtx)

  ! compute area of triangle(vx_i,vy_i)-(vx_j,vy_j)-(point(1,1),point(2,1))
  area1= 0.0d0
  do ivrtx=1, nvrtx

     if(ivrtx==nvrtx) then
        jvrtx=1
     else
        jvrtx=ivrtx+1
     end if

     x1=vrtx(1,ivrtx) ! current vertex point
     y1=vrtx(2,ivrtx)
     x2=vrtx(1,jvrtx) ! next vertex point
     y2=vrtx(2,jvrtx)
     x3=point(1,1) ! given point
     y3=point(2,1)
     area1=area1+0.50d0*abs(x1*y2+y1*x3+y3*x2-y2*x3-y1*x2-x1*y3)
  end do

  ! make decesion
  if( abs( (area0-area1)/area0 ) <= toler(1) ) then
       ptinarea= .true.
  else
       ptinarea= .false.
  end if
  


  return
end function ptinarea





real(8) function permu(indx,jndx,kndx)
  !=======================================================================
  !  permu = compute the value of permutation symbol: e_ijk
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  indx, jndx, kndx : index
  !  
  !  output:
  !  ------
  !  permu : the permutation symbol value
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: indx,jndx,kndx

  ! ====================================
  ! local variable
  ! ==============
  integer :: ncount

  ! loop index
  integer :: index
  ! ====================================

  ! initialize
  permu= 0.0d0

  if ( ( indx==1 .and. jndx==2 .and. kndx==3 ) .or. & 
       ( indx==2 .and. jndx==3 .and. kndx==1 ) .or. &
       ( indx==3 .and. jndx==1 .and. kndx==2 ) ) then
     permu= 1.0d0

  else if ( ( indx==1 .and. jndx==3 .and. kndx==2 ) .or. &
            ( indx==2 .and. jndx==1 .and. kndx==3 ) .or. &
            ( indx==3 .and. jndx==2 .and. kndx==1 ) ) then
     permu= -1.0d0

  else
     permu= 0.0d0

  end if


  return
end function permu





real(8) function distance(ndime, line, point)
  !=======================================================================
  !  distance = compute the distance between a point and line
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : dimension
  !
  !  line(ndime,2) : start and end points of line
  !
  !  point(ndime,1) : coordinate of point
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
  real(8), dimension(ndime,2), intent(in) :: line
  real(8), dimension(ndime,1), intent(in) :: point

  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(ndime) :: vline, vpoin
  real(8) :: vllength, vplength, length

  ! loop index
  integer :: idime
  ! ====================================
  ! initialize
  distance=0.0d0

  ! compute line vector
  vllength=0.0d0 ! initialize
  do idime=1, ndime
     vline(idime)=line(idime,2)-line(idime,1)
     vllength=vllength+vline(idime)**2
  end do
  ! compute length of line vector 
  vllength=dsqrt(vllength)
  ! compute unit vector
  do idime=1, ndime
     vline(idime)=vline(idime)/vllength
  end do

  ! compute point vector
  vplength=0.0d0 ! initialize
  do idime=1, ndime
     vpoin(idime)=point(idime,1)-line(idime,1)
     vplength=vplength+vpoin(idime)**2
  end do
  ! compute length of line vector 
  vplength=dsqrt(vplength)


  ! compute projected length
  length=0.0d0
  do idime=1, ndime
     length=length+vline(idime)*vpoin(idime)
  end do

  ! compute
  distance=dsqrt(vplength**2-length**2)



  return
end function distance





real(8) function hsign(ndime, line, point)
  !=======================================================================
  !  hsign = compute the sign value
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : dimension
  !
  !  line(ndime,2) : start and end points of line
  !
  !  point(ndime,1) : coordinate of point
  ! 
  !  output:
  !  ------
  !  hsign
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime
  real(8), dimension(ndime,2), intent(in) :: line
  real(8), dimension(ndime,1), intent(in) :: point
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: c11, c12 ,c21, c22
  ! ====================================

  ! initialize
  hsign= 0.0d0

  if(ndime /= 2) stop

  c11=line(1,1)-point(1,1)
  c12=line(2,1)-point(2,1)

  c21=line(1,2)-point(1,1)
  c22=line(2,2)-point(2,1)

  if(c11*c22-c12*c21 <= 0.0d0 ) then
     hsign=-1.0d0
  else
     hsign=1.0d0
  end if



  return
end function hsign





real(8) function heavi(xval)
  !=======================================================================
  !  heavi = heaviside function
  !
  !          note:
  !          ----
  !          heavi= 1.0    for xval  > 0.0d0
  !          heavi= 0.0    for xval <= 0.0d0
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  xval : value
  !
  !  output:
  !  ------
  !  heavi : heaviside function value
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), intent(in) :: xval

  ! ====================================

  ! initialize
  heavi= 0.0d0

  if ( xval > 0.0d0 ) then
     heavi= 1.0d0
  else
     heavi= 0.0d0
  end if



  return
end function heavi





real(8) function angle(ndime,point)
  !=======================================================================
  !  angle = compute angle between given line and global x-coordinate
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : dimension
  !
  !  point(ndime,2) : start and end point of line
  ! 
  !  output:
  !  ------
  !  angle : angle between given line and global x-coordinate
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime
  real(8), dimension(ndime,2), intent(in) :: point
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(ndime) :: midpt, vec1, vec2
  real(8) :: length1, length2


  ! loop index
  integer :: idime
  ! ====================================
  ! initialize
  angle=0.d0


  ! find mid point of given line
  do idime=1, ndime
     midpt(idime)=(point(idime,1)+point(idime,2))/2.0d0
  end do

  select case(ndime)

     !
     ! pt1<-----------midpt----------->pt2
     !         vec1           vec2
     !


  case(2) ! 2d case
     ! compute line vector which has the origin at mid point
     length1=0.0d0
     do idime=1, ndime
        vec1(idime)=point(idime,1)-midpt(idime)
        length1=length1+vec1(idime)**2
     end do
     length1=dsqrt(length1)
     ! normalize
     vec1(1:ndime)=vec1(1:ndime)/length1 ! unit vector


     length2=0.0d0
     do idime=1, ndime
        vec2(idime)=point(idime,2)-midpt(idime)
        length2=length2+vec2(idime)**2
     end do
     length2=dsqrt(length2)
     ! normalize
     vec2(1:ndime)=vec2(1:ndime)/length2 ! unit vector

     ! compute angle
     if(vec1(2) > 0.0d0 .and. vec2(2) < 0.0d0) then ! y component
        angle=dacos(vec1(1)) ! (1,0).(vec1x, vec1y)=cos angle

     else if(vec1(2) < 0.0d0 .and. vec2(2) < 0.0d0) then ! y component
        angle=dacos(vec2(1)) ! (1,0).(vec2x, vec2y)=cos angle

     else if(vec1(2) == 0.0d0 .and. vec2(2) == 0.0d0) then
        angle=180.0d0*d2r

     else
        write(*,*) "can not compute angle: angle"
        write(nout5,*) "can not compute angle: angle"
        stop
   
     end if

  end select

  return
end function angle





subroutine getparent2d(optele,nnode,ecord,poinx,poiny, s,t)
  ! =======================================================================
  !  getparent = compute parent coordinate corresponding to current physical coordinate
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optele : element type
  !
  !  ecord(2,nnode) : element nodal coordinate
  !
  !  poinx, poiny : x and y coordinate in physical domain
  !
  !  output:
  !  ------
  !  s,t : corresponding psi and eta coordinate in parent domain
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optele, nnode
  real(8), dimension(2,nnode), intent(in) :: ecord
  real(8), intent(in) :: poinx, poiny

  real(8), intent(out) :: s,t
  ! ====================================
  ! local variable
  ! ==============
  integer :: ncount

  real(8), dimension(nnode) :: shap
  real(8), dimension(2,nnode) :: deriv

  real(8) :: rx, ry ! residual
  real(8) :: drxds, drxdt, dryds, drydt ! components of tangent matrix
  real(8) :: det ! determinant of tangent matrix
  real(8) :: dels, delt ! increment of psi and eta

  ! loop index
  integer :: inode
  ! ====================================
  ! initialize: initial guess
  select case (optele)
  case(1) ! tri
     s=0.50d0 ! psi
     t=0.50d0 ! eta

  case(2:3) ! quad, bt shell(quad)
     s=0.0d0 ! psi
     t=0.0d0 ! eta

  end select


  ! iteration loop to find crosspoints
  ncount=0 ! initialize
  do

     ! iteration counter
     ncount=ncount+1

     ! compute shape function value at parent domain
     call getshape2d(optele,nnode,s,t, shap,deriv)
        ! input : optele,nnode,s,t
        ! output : shap,deriv

     ! compute residual
     rx=0.0d0 ! initialize
     ry=0.0d0
     do inode=1, nnode
        rx=rx+shap(inode)*ecord(1,inode)
        ry=ry+shap(inode)*ecord(2,inode)
     end do
     rx=rx-poinx
     ry=ry-poiny

     ! compute components of tangent matrix and determinant
     drxds=0.0d0 ! initialize
     drxdt=0.0d0
     dryds=0.0d0
     drydt=0.0d0
     do inode=1, nnode
        drxds=drxds+deriv(1,inode)*ecord(1,inode)
        drxdt=drxdt+deriv(2,inode)*ecord(1,inode)
        dryds=dryds+deriv(1,inode)*ecord(2,inode)
        drydt=drydt+deriv(2,inode)*ecord(2,inode)
     end do
     det=drxds*drydt-drxdt*dryds ! determinant

     ! compute increment : del_psi, del_eta
     dels=-( drydt*rx-drxdt*ry)/det
     delt=-(-dryds*rx+drxds*ry)/det

     ! update psi and eta
     s=s+dels
     t=t+delt

     ! check criterion
     if(dsqrt(dels**2+delt**2) <= toler(1)) then
        exit ! quit iteration loop

     else if(ncount >= 50 ) then
        write(*,*) "can not find parent coordinate: getparent2d"
        write(nout5,*) "can not find parent coordinate: getparent2d"
        stop

     end if

  end do

  select case (optele)
  case(1) ! tri
     if ( abs(0.0d0 - s) < toler(1) ) s= 0.0d0
     if ( abs(1.0d0 - s) < toler(1) ) s= 1.0d0

     if ( abs(0.0d0 - t) < toler(1) ) t= 0.0d0
     if ( abs(1.0d0 - t) < toler(1) ) t= 1.0d0


  case(2:3) ! quad, bt shell(quad)
     if ( abs(1.0d0 - s) < toler(1) ) s= 1.0d0
     if ( abs(-1.0d0 - s) < toler(1) ) s= -1.0d0

     if ( abs(1.0d0 - t) < toler(1) ) t= 1.0d0
     if ( abs(-1.0d0 - t) < toler(1) ) t= -1.0d0

  end select


  return
end subroutine getparent2d




subroutine weightff(opttyp,optpol,rmax,dist, weight)
  !=======================================================================
  !  weightff = weight function value
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  opttyp : weight function type
  !           opttyp= 1 : spline function
  !
  !  optpol : polynimial order
  !
  !  rmax : size of influance domain
  !
  !  dist : distance
  !
  !  output:
  !  ------
  !  weight : corresponding weight value
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: opttyp,optpol
  real(8), intent(in) :: rmax,dist

  real(8), intent(out) :: weight
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: normd

  ! ====================================

  ! initialize
  weight=0.0d0


  ! normalized distance
  normd= dist/rmax

  select case(opttyp) ! weight function type
  case(1) ! spline function
     if ( optpol == 3 ) then ! cubic spline
        ! compute weight function value: cubic spline
        if(0.0d0 <= normd .and. normd <= 0.50d0) then
           weight= 4.0d0*( normd - 1.0d0 )*normd**2 + (2.0d0/3.0d0)

        else if( 0.50d0 < normd .and. normd <= 1.0d0 ) then
           weight= (4.0d0/3.0d0) * ( 1.0d0 - normd )**3

        else
           weight=0.0d0

        endif
     else
        write(*,*) "not implemented yet: weightff"
        write(nout5,*) "not implemented yet: weightff"
        stop

     end if

  case default
     write(*,*) "not implemented yet: weightff"
     write(nout5,*) "not implemented yet: weightff"
     stop

  end select


  return
end subroutine weightff





subroutine angfilter(ndime,crpt, angle)
  !=======================================================================
  !  angfilter = filtering crack propagation angle by comparing initial crack angle
  !              with new propagation angle
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  crpt(ndime,2) : crack coordinates
  !
  !  inoutput:
  !  --------
  !  angle : crack propagation angle
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime
  real(8), dimension(ndime,2), intent(in) :: crpt

  real(8), intent(inout) :: angle
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(ndime,1) :: vect1, vect2
  real(8) :: length

  real(8) :: inprd ! inner product

  ! loop index
  integer :: idime   
  ! ====================================

  !                             x
  !                           .   
  !                         .
  !                       .
  !  o---------------->o ---------> [ vect1 ]
  ! tail              tip
  ! crpt1             crpt2

  ! compute crack_tail->crack_tip vector
  do idime=1, ndime
     vect1(idime,1) = crpt(idime,2) - crpt(idime,1)
  end do

  ! compute vector length
  length= 0.0d0
  do idime=1, ndime
     length= length + vect1(idime,1)**2
  end do
  length= dsqrt(length)

  ! get unit vector
  do idime=1, ndime
    vect1(idime,1)= vect1(idime,1) / length
  end do


  select case (ndime)
  case (2)

     ! compute virtual crack vector
     vect2(1,1)= dcos(angle)*1.0d0
     vect2(2,1)= dsin(angle)*1.0d0

     ! compute length
     length= dsqrt(vect2(1,1)**2 + vect2(2,1)**2)

     ! get unit vector
     do idime=1, ndime
        vect2(idime,1)= vect2(idime,1) / length
     end do

  end select

  ! compute inner product
  inprd=0.0d0
  do idime=1, ndime
     inprd= inprd + vect1(idime,1)*vect2(idime,1)
  end do

  ! consider proper opposite direction
  if( inprd < 0.0d0 ) then
     angle= angle+180.0d0*d2r
  end if

  return
end subroutine angfilter





subroutine ludcmp(n, a, indx,d)
  !=======================================================================
  !  ludcmp = decompose given a matrix to LU matrix
  !
  !           ref:
  !           ---
  !           numerical recipes in fortran 77, pp. 36
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  n : the total number of the unknown
  !
  !  inoutput:
  !  --------
  !  a(n,n) : LU decomposed a matrix
  !
  !  output:
  !  ------
  !  indx(n) : records the row permutation effected by the partial pivoting
  !
  !  d : +/-1 depending on whether the row number interchanges was even or odd
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: n

  real(8), dimension(n,n), intent(inout) :: a

  integer, dimension(n), intent(out) :: indx
  real(8), intent(out) :: d
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(n) :: vv

  real(8) :: sum, dum
  real(8) :: aamax

  integer :: imax

  integer :: i, j, k ! loop index
  ! ====================================

  ! initialize
  indx(:)=0

  d=1.0d0
  ! =========================================
  do i=1, n
     aamax=0.0d0

    ! -------------------------------
     do j=1, n
        if ( abs(a(i,j)) > aamax ) then
		   aamax=abs(a(i,j))
        end if
     end do

     if (aamax == 0.0d0) then
        write(*,*) 'singular matrix in ludcmp'
        write(nout5,*) 'singular matrix in ludcmp'
        stop
     else
        vv(i)= 1.0d0/aamax
     end if
 
  end do
  ! =========================================
  do j=1, n

     ! -------------------------------
     do i=1, j-1
        sum= a(i,j)

        do k=1, i-1
           sum= sum - a(i,k)*a(k,j)
        end do

        a(i,j)=sum
     end do
     ! -------------------------------
     aamax=0.d0
     do i=j, n
        sum=a(i,j)
       
        do k=1, j-1
           sum= sum - a(i,k)*a(k,j)
        end do

        a(i,j)=sum
        dum= vv(i)*abs(sum)

        if (dum >= aamax) then
           imax=i
           aamax=dum
        endif

     end do
     ! -------------------------------
     if (j /= imax) then
        do k=1, n
           dum= a(imax,k)
           a(imax,k)= a(j,k)
           a(j,k)= dum
        end do

        d= -d
        vv(imax)=vv(j)
     end if
     ! -------------------------------
     indx(j)= imax
     if( a(j,j) == 0.0d0 ) then 
        a(j,j)= 1.0e-10
     end if

     if( j /= n ) then
        dum= 1.0d0 / a(j,j)
        do i=j+1, n
           a(i,j)= a(i,j)*dum
        end do
     endif
     ! -------------------------------

  end do ! do j=1, n

  return
end subroutine ludcmp





subroutine lubksb(n,a,indx, b)
  !=======================================================================
  !  lubksb = solve lu decomposed matrix by back subsitution
  !           a.x=b -> x=a-1.b
  !
  !           ref:
  !           ---
  !           numerical recipes in fortran 77, pp. 36
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  n : the total number of the unknown
  !
  !  a(n,n) : a matrix
  !
  !  indx(n) : original index matrix
  !
  !  inoutput:
  !  --------
  !  b(n,1) : input -> b matrix of a.x=b
  !           output-> x matrix
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: n
  real(8), dimension(n,n), intent(in) :: a
  integer, dimension(n), intent(in) :: indx

  real(8), dimension(n,1), intent(inout) :: b
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: sum

  integer :: i, j, ii
  integer :: ll
  ! ====================================

  ii=0 ! initialize

  ! --------------------------------
  do i=1, n
     ll = indx(i)
     sum = b(indx(i),1)
     b(ll,1) = b(i,1)

     ! -----------------------------
     if ( ii /= 0 ) then
        do j=ii, i-1
           sum = sum - a(i,j)*b(j,1)
        end do
     else if ( sum /= 0.0d0 ) then
        ii=i
     endif
     ! -----------------------------

     b(i,1)=sum

  end do
  ! --------------------------------
  do i=n, 1, -1

     sum=b(i,1)

     ! -----------------------------
     do j=i+1, n
        sum = sum - a(i,j)*b(j,1)
     end do
     ! -----------------------------

     b(i,1)= sum / a(i,i)

  end do
  ! --------------------------------

  return
end subroutine lubksb





real(8) function getdet(ndime,amat)
  !=======================================================================
  !  getdet = compute determinant of given matrix
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : dimension of tensor
  ! 
  !  amat : given matrix
  !
  !  output:
  !  --------
  !  getdet : determinant value
  !
  !======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime
  real(8), dimension(ndime,ndime), intent(in) :: amat
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(ndime,ndime) :: lumat
  integer, dimension(ndime) :: indx

  ! loop index
  integer :: idime
  ! ====================================

  ! initialize
  getdet= 0.0d0


  select case(ndime)
  case(1) ! 1 by 1 matrix 
     getdet= amat(1,1)

  case(2) ! 2 by 2 matrix 
     getdet= amat(1,1)*amat(2,2) - amat(1,2)*amat(2,1)

  case(3) ! 3 by 3 matrix 
     getdet= -amat(1,3)*amat(2,2)*amat(3,1) + amat(1,2)*amat(2,3)*amat(3,1) &
             +amat(1,3)*amat(2,1)*amat(3,2) - amat(1,1)*amat(2,3)*amat(3,2) &
             -amat(1,2)*amat(2,1)*amat(3,3) + amat(1,1)*amat(2,2)*amat(3,3)

  case(4:) ! general n by n matrix
     ! copy original matrix
     lumat(1:ndime,1:ndime)= amat(1:ndime,1:ndime)

     ! lu decomposition
     call ludcmp(ndime, lumat, indx,getdet)
        ! input : ndime
        ! inoutput : lumat
        ! output : indx,getdet

     ! compute determinant
     do idime=1, ndime
        getdet= getdet * lumat(idime,idime)
     end do

  end select



  return
end function getdet





subroutine getinv(ndime,amat, bmat)
  !=======================================================================
  !  getinv = compute inverse matrix
  !
  !           note:
  !           ----
  !           for 2d and 3d matrix, we use explicit expression
  !           for high order matrix, we compute inverse matrix by
  !           using lu decomposition
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : dimension of matrix
  !
  !  amat(ndime,ndime) : given matrix
  !
  !  output:
  !  ------
  !  bmat(ndime,ndime) : inverse of amat
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime
  real(8), dimension(ndime,ndime), intent(in) :: amat

  real(8), dimension(ndime,ndime), intent(out) :: bmat
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: det, getdet

  real(8), dimension(ndime,ndime) :: lumat
  real(8), dimension(ndime,1) :: temp

  integer, dimension(ndime) :: indx

  real(8) :: d

  ! loop index
  integer :: idime
  ! ====================================

  ! initialize
  bmat(1:ndime,1:ndime)= 0.0d0


  select case(ndime)
  case(1)
    det= amat(1,1)
    if ( det == 0.0d0 ) then
      write(*,*) "determinant of the matrix is zero: getinv"
      write(nout5,*) "determinant of the matrix is zero: getinv"
      stop
    end if

    ! set component
    bmat(1,1)= 1.0d0 / amat(1,1)

  case(2)
    ! compute determinant of amat
    det= getdet(ndime, amat)
    if ( det == 0.0d0 ) then
      write(*,*) "determinant of the matrix is zero: getinv"
      write(nout5,*) "determinant of the matrix is zero: getinv"
      stop
    end if

    ! set component
    bmat(1,1)= amat(2,2)
    bmat(1,2)=-amat(1,2)

    bmat(2,1)=-amat(2,1)
    bmat(2,2)= amat(1,1)

    ! multiply 1/det[a]
    bmat(1:ndime,1:ndime)= (1.0d0/det) * bmat(1:ndime,1:ndime)

  case(3)
    ! compute determinant of amat
    det= getdet(ndime, amat)
    if ( det == 0.0d0 ) then
      write(*,*) 'determinant of the matrix is zero: getinv'
      write(nout5,*) 'determinant of the matrix is zero: getinv'
      stop
    end if

    ! set component
    bmat(1,1)=-amat(2,3) * amat(3,2) + amat(2,2) * amat(3,3)
    bmat(1,2)= amat(1,3) * amat(3,2) - amat(1,2) * amat(3,3)
    bmat(1,3)=-amat(1,3) * amat(2,2) + amat(1,2) * amat(2,3)

    bmat(2,1)= amat(2,3) * amat(3,1) - amat(2,1) * amat(3,3)
    bmat(2,2)=-amat(1,3) * amat(3,1) + amat(1,1) * amat(3,3)
    bmat(2,3)= amat(1,3) * amat(2,1) - amat(1,1) * amat(2,3)

    bmat(3,1)=-amat(2,2) * amat(3,1) + amat(2,1) * amat(3,2)
    bmat(3,2)= amat(1,2) * amat(3,1) - amat(1,1) * amat(3,2)
    bmat(3,3)=-amat(1,2) * amat(2,1) + amat(1,1) * amat(2,2)

    ! multiply 1/det[a]
    bmat(1:ndime,1:ndime)= (1.0d0/det) * bmat(1:ndime,1:ndime)


  case(4:) ! using LU dicomposition and back substitution

     ! copy original matrix
     lumat(1:ndime,1:ndime)= amat(1:ndime,1:ndime)

     ! lu decomposition
     call ludcmp(ndime, lumat, indx,d)
        ! input : ndime
        ! inoutput : lumat
        ! output : indx,d

     ! back substitutin with identity matrix
     do idime=1, ndime

        ! set each column vector of ndime by ndime identity matrix
        temp(1:ndime,1)= 0.0d0 
        temp(idime,1)= 1.0d0

        ! lu back subsitution
        call lubksb(ndime,lumat,indx, temp)
           ! input : ndime,lumat,indx
           ! output : temp

        ! set result to output matrix
        bmat(1:ndime,idime)= temp(1:ndime,1)

     end do

  end select



  return
end subroutine getinv



subroutine getpertub(nptb,minval,maxval, ptbval)
  !=======================================================================
  !  getpertub = get randomly perturbed value within the range of [ minval - maxval ]
  !
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  nptb : the total number of perturbed value
  !
  !  minval, maxval : range of perturbed value
  !
  !  output:
  !  ------
  !  ptbval(nptb) : perturbed value within the rage of [ minval - maxval ]
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: nptb
  real(8), intent(in) :: minval, maxval

  real(8), dimension(nptb), intent(out) :: ptbval
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: rvalue

  integer :: indx ! loop index
  ! ====================================

  ! initialize
  ptbval(:)= 0.0d0


  ! check minval and maxval
  if (minval > maxval) then
     write(*,*) "minval > maxval: getpertub"
     write(nout5,*) "minval > maxval: getpertub"
     stop
  end if

  ! make seed
  call random_seed()

  do indx=1, nptb
     ! get random number
     call random_number(rvalue) ! rvalue is within the range of [0,1]

     ! set perturbation range
     ptbval(indx)= minval + (maxval - minval) * rvalue
  end do


  return
end subroutine getpertub





subroutine gettvec(ndime,line, tvec)
  !=======================================================================
  !  gettvec = get tangential vector
  !
  !            note:
  !            ----
  !              o-------------------->o
  !        line(ndime,1)        line(ndime,2)
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : dimension
  !
  !  line(ndime,2) : start and end point of line
  !
  !  output:
  !  ------
  !  tvec(ndime,1) : unit tangential vector
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime
  real(8), dimension(ndime,2), intent(in) :: line

  real(8), dimension(ndime,1), intent(out) :: tvec
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: length

  integer :: idime ! loop index
  ! ====================================

  ! initialize
  tvec(:,:)= 0.0d0


  !      o-------------------->o
  ! line(ndime,1)        line(ndime,2)

  ! compute tangential vector
  length=0.0d0 ! initialize
  do idime=1, ndime
     tvec(idime,1)= line(idime,2) - line(idime,1)
     length= length + tvec(idime,1)**2
  end do

  ! compute length of line vector 
  length=dsqrt(length)

  ! compute unit vector
  tvec(:,1)= tvec(:,1) / length

  return
end subroutine gettvec





subroutine hpsort(n, ra, indx)
  !=======================================================================
  !  hpsort = sorts an array ra(1:n) into ascending numerical order
  !           using the heap sort algorithm.
  !           also, return the original array index of sorted data
  !
  !           note:
  !           ----
  !           data will be sorted with ascending order
  !
  !           ref:
  !           ----
  !           numerical recipes in fortran 77, pp. 327
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  n : the total number of input data
  !
  !  inoutput:
  !  --------
  !  value(*) : replaced on output by its sorted rearrangement
  !
  !  output:
  !  ------
  !  indx(*) : the original array index of the sorted data
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) ::  n

  real(8), dimension(*), intent(inout) :: ra

  integer, dimension(*), intent(out) :: indx
  ! ====================================
  ! local variable
  ! ==============
  integer :: i,ir,j,l
  real(8) :: rra ! copy ra
  integer :: rndx ! copy indx
  ! ====================================

  ! initialize: do not initialize ra array
  do i=1, n
     indx(i)= i
  end do

  ! need at leat 1 component
  if (n <= 1) return

  l=n/2+1
  ir=n

10 continue

  if(l > 1)then ! still in hiring phase
     l= l-1
     rra= ra(l)
     ! -------------
     rndx= indx(l)
     ! -------------

  else ! in retirement and promotion phase
     rra= ra(ir) ! clear a space at end of array
     ra(ir)= ra(1) ! retire the top of the heap into it
     ! -------------
     rndx= indx(ir)  ! for index
     indx(ir)=indx(1)
     ! -------------

     ir=ir-1 ! decrease the size of the corporation

     if(ir == 1)then ! done with the last promotion
        ra(1)= rra ! the least competent worker of all
        ! -------------
        indx(1)= rndx
        ! -------------
        return
     end if

  end if

  i=l   ! whether in the hiring phase or promotion phase,
  j=l+l ! we here set up to sift down element rra to its proper level

20 continue

  if (j <= ir) then ! do while j <= ir
 
     if (j < ir) then
        if (ra(j) < ra(j+1)) j=j+1 ! compare to the better underling
     end if

     if (rra < ra(j)) then ! demote rra
        ra(i)= ra(j)
        ! -------------
        indx(i)= indx(j)
        ! -------------
        i= j
        j= j+j
     else ! this is rra's level. Set j to terminate the sift-down
        j= ir+1
     end if

     goto 20
  endif

  ra(i)= rra ! put rra into its slot
  ! -------------
  indx(i)= rndx
  ! -------------
  goto 10


  return
end subroutine hpsort





subroutine rotpt2d(ndime,theta, poin)
  !=======================================================================
  !  rotpt2d= compute rotated point coordinate w.r.t original x coordinate
  !
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : dimension
  !
  !  theta : rotation angle
  !
  !  output:
  !  ------
  !  poin(ndime,1) : rotated point
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime
  real(8), intent(in) :: theta
 
  real(8), dimension(ndime,1), intent(out) :: poin
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(ndime,1) :: xdirvec

  real(8), dimension(ndime,ndime) :: rotens2d
  ! ====================================

  ! initialize
  poin(:,:)= 0.0d0

  ! error check
  if(ndime /= 2) then
     write(*,*) "this subroutine is only for 2d: rotpt2d"
     write(nout5,*) "this subroutine is only for 2d: rotpt2d"
     stop
  end if

  ! set x direction base vector
  xdirvec(1,1)= 1.0d0
  xdirvec(2,1)= 0.0d0

  ! get rotation tensor
  call getrotens2d(ndime,theta, rotens2d)
     ! input : ndime,theta
     ! output : rotens2d

  ! compute 
  call matprd(ndime,ndime,0, ndime,1,0, ndime,1, rotens2d,xdirvec, poin)
     ! input : ndime,ndime,0, ndime,1,0, ndime,1, rotens2d,xdirvec
     ! output : poin

  return
end subroutine  rotpt2d





subroutine getlineseg(ndime,seglen,segang,segorg, segline)
  !=======================================================================
  !  getcrseg= create line segment
  !
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : dimension
  !
  !  seglen : segment length
  !
  !  segang : segmeent angle
  !
  !  segorg(ndime,1) : origin of line segment
  !
  !  output:
  !  ------
  !  segline(ndime,2) : start and end point of line segmenet
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime
  real(8), intent(in) :: seglen, segang
  real(8), dimension(ndime,1), intent(in) :: segorg

  real(8), dimension(ndime,2), intent(out) :: segline
  ! ====================================
  ! local variable
  ! ============== 
  real(8), dimension(ndime,1) :: segtip
  ! ====================================

  ! initialize
  segline(:,:)= 0.0d0

  select case(ndime)
  case(2)
     ! get the segment tip unit vector
     call rotpt2d(ndime,segang, segtip)
        ! input : ndime,segang
        ! output : segtip

     ! set the size of the segmenet tip vector
     segtip(:,:)= segtip(:,:)*seglen

     ! transform origin of segmenet vector and set start and end point
     segline(1:ndime,1)= segorg(1:ndime,1)
     segline(1:ndime,2)= segorg(1:ndime,1) + segtip(1:ndime,1)

  case default
     write(*,*) "not implemented yet: getlineseg"
     write(nout5,*) "not implemented yet: getlineseg"
     stop

  end select



  return
end subroutine getlineseg





subroutine crsprdt3d(opt,avec,bvec, cvec)
  !=======================================================================
  !  crsprdt3d = cross product of given vector
  !
  !            c= a x b
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  opt : normalize option
  !        opt=0 : c= a x b
  !        opt=1 : c= a x b / ||a x b||
  ! 
  !  avec(3,1) : given vector a
  !                     
  !  bvec(3,1) : given vector b
  !
  !  output:
  !  ------
  !  cvec(3,1) : cross product resultant vector
  !                            
  ! ======================================================================
  
  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: opt
  real(8), dimension(3,1), intent(in) :: avec, bvec

  real(8), dimension(3,1), intent(out) :: cvec
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: length
  ! ====================================

  ! initialize
  cvec(:,:)= 0.0d0


  ! set cross product result to c vector: c= a x b
  cvec(1,1)=  avec(2,1)*bvec(3,1) - avec(3,1)*bvec(2,1)
  cvec(2,1)= -avec(1,1)*bvec(3,1) + avec(3,1)*bvec(1,1) 
  cvec(3,1)=  avec(1,1)*bvec(2,1) - avec(2,1)*bvec(1,1)

  ! check normalize option
  if ( opt == 1 ) then
     ! normalize result vector
     call unitvec1(3,cvec)
        ! input : ndime
        ! inoutput : cvec
  end if
  

  return
end subroutine crsprdt3d



real(8) function dotprdt(ndime,avec,bvec)
  !=======================================================================
  !  dotprdt = dot product
  !            c= a^t . b
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ndime : dimension of given vector
  !
  !  avec(ndime,1) : given vector a
  !                     
  !  bvec(ndime,1) : given vector b
  !
  !  output:
  !  ------
  !  dotprdt : dot product result
  !                            
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime
  real(8), dimension(ndime,1), intent(in) :: avec, bvec
  ! ====================================
  ! local variable
  ! ==============

  ! loop index
  integer :: idime
  ! ====================================

  ! initialize
  dotprdt= 0.0d0

  ! compute dot product
  do idime=1, ndime
     dotprdt= dotprdt + avec(idime,1)*bvec(idime,1)
  end do



  return
end function dotprdt



subroutine unitvec1(ndime, avec)
  !=======================================================================
  !  unitvec1 = normalize given vector
  !            a= a / ||a||
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  opt : normalize option
  ! 
  !  ndime : dimension of given vector
  !
  !  inoutput:
  !  --------
  !  avec(ndime,1) : unit vector
  !                            
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime

  real(8), dimension(ndime,1), intent(inout) :: avec
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: length

  ! loop index
  integer :: idime
  ! ====================================

  ! initialize : do not initialize avec

  ! compute length of given vector
  length= 0.0d0 ! initialize
  do idime=1, ndime
     length= length + avec(idime,1)**2
  end do
  length= dsqrt(length)

  if( length == 0.0d0 ) then
     write(*,*) "given vector is zero vector: unitvec1"
     write(nout5,*) "given vector is zero vector: unitvec1"

  else
     ! normalize given vector
     avec(1:ndime,1)= avec(1:ndime,1) / length

  end if



  return
end subroutine unitvec1



subroutine unitvec2(ndime,avec, bvec)
  !=======================================================================
  !  unitvec2 = normalize given vector
  !            b= a / ||a||
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  opt : normalize option
  ! 
  !  ndime : dimension of given vector
  !
  !  avec(ndime,1) : given vector
  !
  !  output:
  !  --------
  !  bvec(ndime,1) : unit vector
  !                            
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: ndime
  real(8), dimension(ndime,1), intent(in) :: avec

  real(8), dimension(ndime,1), intent(out) :: bvec
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: length

  ! loop index
  integer :: idime
  ! ====================================

  ! initialize
  bvec(:,:)= 0.0d0

  ! compute length of given vector
  length= 0.0d0 ! initialize
  do idime=1, ndime
     length= length + avec(idime,1)**2
  end do
  length= dsqrt(length)

  if( length == 0.0d0 ) then
     write(*,*) "given vector is zero vector: unitvec2"
     write(nout5,*) "given vector is zero vector: unitvec2"

  else
     ! normalize given vector
     bvec(1:ndime,1)= avec(1:ndime,1)/length

  end if



  return
end subroutine unitvec2




subroutine mlsfit2dqd(npt,rmax,ptcord,ptdist,ptval,chkpt, dvaldx,dvaldy,fitcoef)
  !=======================================================================
  !  mlsfit2dln = linear moving least square fit at given point
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  npt : the number of point
  !
  !  rmax : maximum distance (radius of influance domain)
  !
  !  ptcord(2,npt) : coordinate of point
  !
  !  ptdist(npt,1) : distance
  !
  !  ptval(npt,1) : value at point
  !
  !  chkpt(2,1) : point at calculating derivatives
  !
  !  inoutput:
  !  --------
  !  dvaldx : derivative wrt x
  !
  !  dvaldy : derivative wrt y
  !
  !  fitcoef(3,1) : mls fit coefficient at given point
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: npt
  real(8), intent(in) :: rmax
  
  real(8), dimension(2,npt), intent(in) :: ptcord
  real(8), dimension(npt,1), intent(in) :: ptdist
  real(8), dimension(npt,1), intent(in) :: ptval
  real(8), dimension(2,1), intent(in) :: chkpt
  
  real(8), intent(out) :: dvaldx,dvaldy
  real(8), dimension(6,1), intent(out) :: fitcoef
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: dist, weight

  real(8), dimension(npt,npt) :: wmat
  real(8), dimension(npt,6) :: pmat
  real(8), dimension(6,6) :: mmat, mmatinv
  real(8), dimension(6,1) :: rhsmat

  ! loop index
  integer :: ipt
  ! ====================================

  ! initialize
  dvaldx= 0.0d0
  dvaldy= 0.0d0
  fitcoef(:,:)= 0.0d0


  ! check the number of data point
  if ( npt < 6 ) then
     write(*,*) "we have few number of data points: mlsfit2dqd"
     write(nout5,*) "we have few number of data points: mlsfit2dqd"
     stop

  end if

  ! mls fitting
  wmat(1:npt,1:npt)= 0.0d0 ! initialize weight matrix
  pmat(1:npt,1:6)= 0.0d0 ! initialize polynomial matrix
  do ipt=1, npt

     ! get distance
     dist= ptdist(ipt,1)

     ! ----------------------------------------------
     ! compute w matrix
     ! ----------------
     ! compute weight function
     call weightff(1,3,rmax,dist, weight)
        ! input : 1(opttyp: cubic spline),3(optpol: polynomial order),rmax,dist
        ! output : weight

     ! set weight matrix
     wmat(ipt,ipt)= weight

     ! ----------------------------------------------
     ! compute p matrix: linear polynomial
     ! ----------------
     pmat(ipt,1)= 1.0d0
     pmat(ipt,2)= ptcord(1,ipt) ! x
     pmat(ipt,3)= ptcord(2,ipt) ! y
     pmat(ipt,4)= ptcord(1,ipt) * ptcord(2,ipt) ! xy
     pmat(ipt,5)= ptcord(1,ipt)**2 ! x**2
     pmat(ipt,6)= ptcord(2,ipt)**2 ! y**2

     ! ----------------------------------------------

  end do

  ! compute m matrix
  call mlsmmat(npt,6,pmat,wmat, mmat)
     ! input : npt,6(npol),pmat,wmat
     ! output : mmat

  ! compute inverse of mmat: mmatinv
  call getinv(6,mmat, mmatinv)
     ! input : 6(ndime),mmat
     ! output : mmatinv

  ! compute rhs matrix
  call mlsrhsmat(npt,6,pmat,wmat,ptval, rhsmat)
     ! input : npt,6(npol),pmat,wmat,ptval
     ! output : rhsmat


  ! compute [mmat]-1.rhsmat
  call matprd(6,6,0, 6,1,0, 6,1, mmatinv,rhsmat, fitcoef)
     ! input : 6,6,0, 6,1,0, 6,1, mmatinv,rhsmat
     ! output : fitcoef

  ! derivative wrt x
  dvaldx= fitcoef(2,1) + fitcoef(4,1)*chkpt(2,1) + 2.0d0*fitcoef(5,1)*chkpt(1,1)

  ! derivative wrt y
  dvaldy= fitcoef(3,1) + fitcoef(4,1)*chkpt(1,1) + 2.0d0*fitcoef(6,1)*chkpt(2,1)



  return
end subroutine mlsfit2dqd





subroutine mlsfit2dln(npt,rmax,ptcord,ptdist,ptval, dvaldx,dvaldy,fitcoef)
  !=======================================================================
  !  mlsfit2dln = linear moving least square fit at given point
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  npt : the number of point
  !
  !  rmax : maximum distance (radius of influance domain)
  !
  !  ptcord(2,npt) : coordinate of point
  !
  !  ptdist(npt,1) : distance
  !
  !  ptval(npt,1) : value at point
  !
  !  inoutput:
  !  --------
  !  dvaldx : derivative wrt x
  !
  !  dvaldy : derivative wrt y
  !
  !  fitcoef(3,1) : mls fit coefficient at given point
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: npt
  real(8), intent(in) :: rmax
  
  real(8), dimension(2,npt), intent(in) :: ptcord
  real(8), dimension(npt,1), intent(in) :: ptdist
  real(8), dimension(npt,1), intent(in) :: ptval
  
  real(8), intent(out) :: dvaldx,dvaldy
  real(8), dimension(3,1), intent(out) :: fitcoef
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: dist, weight

  real(8), dimension(npt,npt) :: wmat
  real(8), dimension(npt,3) :: pmat
  real(8), dimension(3,3) :: mmat, mmatinv
  real(8), dimension(3,1) :: rhsmat

  ! loop index
  integer :: ipt
  ! ====================================

  ! initialize
  dvaldx= 0.0d0
  dvaldy= 0.0d0
  fitcoef(:,:)= 0.0d0


  ! check the number of data point
  if ( npt < 3 ) then
     write(*,*) "we have few number of data points: mlsfit2dln"
     write(nout5,*) "we have few number of data points: mlsfit2dln"
     stop

  end if

  ! mls fitting
  wmat(1:npt,1:npt)= 0.0d0 ! initialize weight matrix
  pmat(1:npt,1:3)= 0.0d0 ! initialize polynomial matrix
  do ipt=1, npt

     ! get distance
     dist= ptdist(ipt,1)

     ! ----------------------------------------------
     ! compute w matrix
     ! ----------------
     ! compute weight function
     call weightff(1,3,rmax,dist, weight)
        ! input : 1(opttyp: cubic spline),3(optpol: polynomial order),rmax,dist
        ! output : weight

     ! set weight matrix
     wmat(ipt,ipt)= weight

     ! ----------------------------------------------
     ! compute p matrix: linear polynomial
     ! ----------------
     pmat(ipt,1)= 1.0d0
     pmat(ipt,2)= ptcord(1,ipt) ! x_ipt
     pmat(ipt,3)= ptcord(2,ipt) ! y_ipt
     ! ----------------------------------------------

  end do

  ! compute m matrix
  call mlsmmat(npt,3,pmat,wmat, mmat)
     ! input : npt,3(npol),pmat,wmat
     ! output : mmat

  ! compute inverse of mmat: mmatinv
  call getinv(3,mmat, mmatinv)
     ! input : 3(ndime),mmat
     ! output : mmatinv

  ! compute rhs matrix
  call mlsrhsmat(npt,3,pmat,wmat,ptval, rhsmat)
     ! input : npt,3(npol),pmat,wmat,ptval
     ! output : rhsmat


  ! compute [mmat]-1.rhsmat
  call matprd(3,3,0, 3,1,0, 3,1, mmatinv,rhsmat, fitcoef)
     ! input : 3,3,0, 3,1,0, 3,1, mmatinv,rhsmat
     ! output : fitcoef

  ! derivative wrt x
  dvaldx= fitcoef(2,1)

  ! derivative wrt y
  dvaldy= fitcoef(3,1)



  return
end subroutine mlsfit2dln





subroutine mlsrhsmat(npt,npol,pmat,wmat,ptval, rhsmat)
  !=======================================================================
  !  mlsrhsmat = m matrix for moving least square fit
  !
  !              note:
  !              ----
  !              rhsmat= pmat^t.wmat.ptval
  !              [npol by 1]= [npt by npol]^t.[npt by npt].[npt by 1] 
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  npt : the total number of points
  !
  !  npol : the number of componenets of polynomial set
  !
  !  pmat(npt,npol) : polynomial set matrix
  !
  !  wmat(npt,npt) : weight matrix
  !
  !  ptval(npt,1) : value at neighbour point
  !
  !  output:
  !  ------
  !  rhsmat(npol,1) : rhs matrix
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: npt, npol
  real(8), dimension(npt,npol), intent(in) :: pmat
  real(8), dimension(npt,npt), intent(in) :: wmat
  real(8), dimension(npt,1), intent(in) :: ptval
  
  real(8), dimension(npol,1), intent(out) :: rhsmat
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(npol,npt) :: temp

  ! ====================================

  ! initialize
  rhsmat(:,:)= 0.0d0


  ! compute pmat^t.wmat
  call matprd(npt,npol,1, npt,npt,0, npol,npt, pmat,wmat, temp)
     ! input : npt,npol,1, npt,npt,0, npol,npt, pmat,wmat
     ! output : temp

  ! compute [pmat^t.wmat].ptval
  call matprd(npol,npt,0, npt,1,0, npol,1, temp,ptval, rhsmat)
     ! input : npol,npt,0, npt,1,0, npol,1, temp,ptval
     ! output : rhsmat



  return
end subroutine mlsrhsmat





subroutine mlsmmat(npt,npol,pmat,wmat, mmat)
  !=======================================================================
  !  mlsmmat = m matrix for moving least square fit
  !
  !           note:
  !           ----
  !           mmat= pmat^t.wmat.pmat
  !           [npol by npol]= [npt by npol]^t.[npt by npt].[npt by npol] 
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  npt : the total number of points
  !
  !  npol : the number of componenets of polynomial set
  !
  !  pmat(npt,npol) : polynomial set matrix
  !
  !  wmat(npt,npt) : weight matrix
  !
  !  output:
  !  ------
  !  mmat(npol,npol) : m matrix
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: npt, npol
  real(8), dimension(npt,npol), intent(in) :: pmat
  real(8), dimension(npt,npt), intent(in) :: wmat
  
  real(8), dimension(npol,npol), intent(out) :: mmat
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(npol,npt) :: temp

  ! ====================================

  ! initialize
  mmat(:,:)= 0.0d0


  ! compute pmat^t.wmat
  call matprd(npt,npol,1, npt,npt,0, npol,npt, pmat,wmat, temp)
     ! input : npt,npol,1, npt,npt,0, npol,npt, pmat,wmat
     ! output : temp

  ! compute [pmat^t.wmat].pmat
  call matprd(npol,npt,0, npt,npol,0, npol,npol, temp,pmat, mmat)
     ! input : npol,npt,0, npt,npol,0, npol,npol, temp,pmat
     ! output : mmat



  return
end subroutine mlsmmat





subroutine elearea(optele,ndime,ecord, earea)
  !=======================================================================
  !  earea= compute 2d area of given element
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optele : element type option handler
  !
  !  ndime : problem dimension
  !
  !  ecord(2,*) : element nodal coordinates
  !
  !  output:
  !  ------
  !  earea : 2d area of given element
  !
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optele, ndime
  real(8), dimension(ndime,*), intent(in) :: ecord

  real(8), intent(out) :: earea
  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(3,3) :: locbvec
  real(8), dimension(3,4) :: ecordloc

  real(8), dimension(2,3) :: vrtx3
  real(8), dimension(2,4) :: vrtx4

  real(8) :: area
  ! ====================================

  ! initialize
  earea= 0.0d0


  ! compute element area
  select case(optele)
  case(1) ! 2d tri
     vrtx3(1:2,1:3)= ecord(1:2,1:3)
     earea= area(3,vrtx3)

  case(2) ! 2d quad
     vrtx4(1:2,1:4)= ecord(1:2,1:4) 
     earea= area(4,vrtx4)

  case(3) ! 3d tb shell
     ! compute co rotational local base vector: locbvec
     call getlocbvecbt(ecord, locbvec)
        ! input : ecord
        ! output : locbvec

     ! get local nodal coordinate: glb -> loc
     call glb2locnodv(3,4,3,locbvec,ecord, ecordloc)
        ! input : 3(ndime),4(nnode),3(ntrndof),locbvec,ecord
        ! output : ecordloc

     vrtx4(1:2,1:4)= ecordloc(1:2,1:4) 
     earea= area(4,vrtx4)

  case default
     write(*,*) "not available: elearea"
     write(nout5,*) "not available: elearea"
     stop

  end select



  return
end subroutine elearea





subroutine elehleng0(optele,ndime,ecord, ehleng)
  !=======================================================================
  !  elehleng0 = compute average element size h
  !
  !              2d triangular, quadrilateral
  !              3d shell
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  optele : element type option
  !
  !  ndime : element basic dimension
  !
  !  ecord(ndime,*) : nodal coordinate
  !
  !  output:
  !  ------
  !  ehleng : characteristic element h length
  ! 
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  integer, intent(in) :: optele
  integer, intent(in) :: ndime
  real(8), dimension(ndime,*), intent(in) :: ecord

  real(8), intent(out) :: ehleng
  ! ====================================
  ! local variable
  ! ==============
  integer :: nvrtx
  real(8), dimension(3,3) :: locbvec
  real(8), dimension(3,4) :: ecordloc

  real(8), dimension(2,3) :: vrtx3
  real(8), dimension(2,4) :: vrtx4
  real(8) :: area

  real(8) :: vol3d8nod
  ! ====================================

  ! initialize
  ehleng= 0.0d0


  select case(optele)
  case(0) ! 1d line element
     ehleng= dsqrt( ( ecord(1,1) - ecord(1,2) )**2 )

  case(1) ! 2d tri
     nvrtx=3
     vrtx3(1:2,1:3)= ecord(1:2,1:3)
     ehleng= 1.50d0 * dsqrt( area(nvrtx,vrtx3) )

  case(2) ! 2d quad
     nvrtx=4
     vrtx4(1:2,1:4)= ecord(1:2,1:4) 
     ehleng= dsqrt( area(nvrtx,vrtx4) )

  case(3) ! 3d tb shell
     ! compute co rotational local base vector: locbvec
     call getlocbvecbt(ecord, locbvec)
        ! input : ecord
        ! output : locbvec

     ! get local nodal coordinate: glb -> loc
     call glb2locnodv(3,4,3,locbvec,ecord, ecordloc)
        ! input : 3(ndime),4(nnode),3(ntrndof),locbvec,ecord
        ! output : ecordloc

     nvrtx=4
     vrtx4(1:2,1:4)= ecordloc(1:2,1:4) 
     ehleng= dsqrt( area(nvrtx,vrtx4) )

  case(4) ! 3d brick
     ehleng= ( vol3d8nod(ecord) ) ** (1.0d0/3.0d0)

  case default
     write(*,*) "not available: gethleng"
     write(nout5,*) "not available: gethleng"
     stop

  end select



  return
end subroutine elehleng0

