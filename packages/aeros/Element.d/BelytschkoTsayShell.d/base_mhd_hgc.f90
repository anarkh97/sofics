! ==================================
! hourglass control scheme
! ==================================
!      type                  name                              argument
!      ----                  ----                              --------
! 1.  subroutine          getgamma4nod                  (ecordloc, gamma)
! 2.  real(8) function    vol3d8nod                     (ecord)
! 3.  subroutine          getcijk3d8nod                 (cijk)
!
! =========================================================================================================




subroutine getgamma4nod(ecordloc, area, gamma, zgamma)
  !=======================================================================
  !  getgamma4nod = compute gamma projection operator of 4 node quad. element
  !
  !                 note:
  !                 ----
  !                 flanagan and belytschko, IJNME, 1981, vol. 17, pp. 679-706
  !                 a uniform strain hexahedron and quadrilateral
  !                 with orthogonal hourglass control
  !                 (see, eq (72))
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ecordloc(2,4) : local element nodal coordinate
  !
  !  output:
  !  ------
  !  gamma(4,1) : gamma projection operator
  !                            
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(3,4), intent(in) :: ecordloc
  real(8), intent(in) :: area

  real(8), dimension(4,1), intent(out) :: gamma
  real(8), intent(out) :: zgamma
  ! ====================================
  ! local variable
  ! ==============
  real(8) :: const
  real(8) :: x1,x2,x3,x4
  real(8) :: y1,y2,y3,y4
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
  const= 1.0d0 / ( 4.0d0 * area )

  ! set gamm projection operator for one point integration
  gamma(1,1)= const*( x2*(y3-y4) + x3*(y4-y2) + x4*(y2-y3) )
  gamma(2,1)= const*( x3*(y1-y4) + x4*(y3-y1) + x1*(y4-y3) )
  gamma(3,1)= const*( x4*(y1-y2) + x1*(y2-y4) + x2*(y4-y1) )
  gamma(4,1)= const*( x1*(y3-y2) + x2*(y1-y3) + x3*(y2-y1) )

  ! compute gamma_i z_i
  zgamma= 0.0d0 ! initialize
  do inode=1, 4
     zgamma= zgamma + ecordloc(3,inode)*gamma(inode,1)
  end do

  return
end subroutine getgamma4nod




real(8) function vol3d8nod(ecord)
  !=======================================================================
  !  vol3d8nod = compute appriximated the volume of 3d 8 node hexahedron element
  !
  !                  note:
  !                  ----
  !                  flanagan and belytschko, IJNME, 1981, vol.17 , pp. 683, eq.(19)
  !                  a univorm strain hexahedron and quadrilateral with orthogonal hourglas control
  !
  !  arguments description
  !  ---------------------
  !  input:
  !  -----
  !  ecord(3,8) : element nodal coordinate data
  !
  !  output:
  !  ------
  !  vol3d8nod : approximated volume
  !                            
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(3,8), intent(in) :: ecord

  ! ====================================
  ! local variable
  ! ==============
  real(8), dimension(8,8,8) :: cijk

  ! loop index
  integer :: inode, jnode, knode
  ! ====================================

  ! initialize
  vol3d8nod= 0.0d0

  ! set c_ijk value: eq.(20)
  call getcijk3d8nod(cijk)
     ! output : cijk

  ! loop over node
  vol3d8nod= 0.0d0 ! initialize
  do inode=1, 8
     do jnode=1, 8
        do knode=1, 8

           vol3d8nod= vol3d8nod + ecord(1,inode)*ecord(2,jnode)*ecord(3,knode)*cijk(inode,jnode,knode)

        end do
     end do
  end do



  return
end function vol3d8nod





subroutine getcijk3d8nod(cijk)
  !=======================================================================
  !  getcijk3d8nod = compute c_ijk for 3d 8node hexahedron
  !
  !                  note:
  !                  ----
  !                  flanagan and belytschko, NJNME, 1981, vol.17 , pp. 683, eq.(20)
  !                  a univorm strain hexahedron and quadrilateral with orthogonal hourglas control
  !
  !  arguments description
  !  ---------------------
  !  output:
  !  ------
  !  cijk(8,8,8): pre computed c_ijk value
  !                            
  ! ======================================================================

  include 'preset.fi'
  ! ====================================
  ! subroutine argument
  ! ===================
  real(8), dimension(8,8,8), intent(out) :: cijk
  ! ====================================

  ! initialize
  cijk(1:8,1:8,1:8)= 0.0d0


  ! set pre computed nonzero components only
  cijk(1,2,3)= -1.0d0 / 12.0d0
  cijk(1,2,4)= -1.0d0 / 12.0d0
  cijk(1,2,5)=  1.0d0 / 12.0d0
  cijk(1,2,6)=  1.0d0 / 12.0d0
  cijk(1,3,2)=  1.0d0 / 12.0d0
  cijk(1,3,4)= -1.0d0 / 12.0d0
  cijk(1,4,2)=  1.0d0 / 12.0d0
  cijk(1,4,3)=  1.0d0 / 12.0d0
  cijk(1,4,5)= -1.0d0 / 12.0d0
  cijk(1,4,8)= -1.0d0 / 12.0d0
  cijk(1,5,2)= -1.0d0 / 12.0d0
  cijk(1,5,4)=  1.0d0 / 12.0d0
  cijk(1,5,6)= -1.0d0 / 12.0d0
  cijk(1,5,8)=  1.0d0 / 12.0d0
  cijk(1,6,2)= -1.0d0 / 12.0d0
  cijk(1,6,5)=  1.0d0 / 12.0d0
  cijk(1,8,4)=  1.0d0 / 12.0d0
  cijk(1,8,5)= -1.0d0 / 12.0d0
  cijk(2,1,3)=  1.0d0 / 12.0d0
  cijk(2,1,4)=  1.0d0 / 12.0d0
  cijk(2,1,5)= -1.0d0 / 12.0d0
  cijk(2,1,6)= -1.0d0 / 12.0d0
  cijk(2,3,1)= -1.0d0 / 12.0d0
  cijk(2,3,4)= -1.0d0 / 12.0d0
  cijk(2,3,6)=  1.0d0 / 12.0d0
  cijk(2,3,7)=  1.0d0 / 12.0d0
  cijk(2,4,1)= -1.0d0 / 12.0d0
  cijk(2,4,3)=  1.0d0 / 12.0d0
  cijk(2,5,1)=  1.0d0 / 12.0d0
  cijk(2,5,6)= -1.0d0 / 12.0d0
  cijk(2,6,1)=  1.0d0 / 12.0d0
  cijk(2,6,3)= -1.0d0 / 12.0d0
  cijk(2,6,5)=  1.0d0 / 12.0d0
  cijk(2,6,7)= -1.0d0 / 12.0d0
  cijk(2,7,3)= -1.0d0 / 12.0d0
  cijk(2,7,6)=  1.0d0 / 12.0d0
  cijk(3,1,2)= -1.0d0 / 12.0d0
  cijk(3,1,4)=  1.0d0 / 12.0d0
  cijk(3,2,1)=  1.0d0 / 12.0d0
  cijk(3,2,4)=  1.0d0 / 12.0d0
  cijk(3,2,6)= -1.0d0 / 12.0d0
  cijk(3,2,7)= -1.0d0 / 12.0d0
  cijk(3,4,1)= -1.0d0 / 12.0d0
  cijk(3,4,2)= -1.0d0 / 12.0d0
  cijk(3,4,7)=  1.0d0 / 12.0d0
  cijk(3,4,8)=  1.0d0 / 12.0d0
  cijk(3,6,2)=  1.0d0 / 12.0d0
  cijk(3,6,7)= -1.0d0 / 12.0d0
  cijk(3,7,2)=  1.0d0 / 12.0d0
  cijk(3,7,4)= -1.0d0 / 12.0d0
  cijk(3,7,6)=  1.0d0 / 12.0d0
  cijk(3,7,8)= -1.0d0 / 12.0d0
  cijk(3,8,4)= -1.0d0 / 12.0d0
  cijk(3,8,7)=  1.0d0 / 12.0d0
  cijk(4,1,2)= -1.0d0 / 12.0d0
  cijk(4,1,3)= -1.0d0 / 12.0d0
  cijk(4,1,5)=  1.0d0 / 12.0d0
  cijk(4,1,8)=  1.0d0 / 12.0d0
  cijk(4,2,1)=  1.0d0 / 12.0d0
  cijk(4,2,3)= -1.0d0 / 12.0d0
  cijk(4,3,1)=  1.0d0 / 12.0d0
  cijk(4,3,2)=  1.0d0 / 12.0d0
  cijk(4,3,7)= -1.0d0 / 12.0d0
  cijk(4,3,8)= -1.0d0 / 12.0d0
  cijk(4,5,1)= -1.0d0 / 12.0d0
  cijk(4,5,8)=  1.0d0 / 12.0d0
  cijk(4,7,3)=  1.0d0 / 12.0d0
  cijk(4,7,8)= -1.0d0 / 12.0d0
  cijk(4,8,1)= -1.0d0 / 12.0d0
  cijk(4,8,3)=  1.0d0 / 12.0d0
  cijk(4,8,5)= -1.0d0 / 12.0d0
  cijk(4,8,7)=  1.0d0 / 12.0d0
  cijk(5,1,2)=  1.0d0 / 12.0d0
  cijk(5,1,4)= -1.0d0 / 12.0d0
  cijk(5,1,6)=  1.0d0 / 12.0d0
  cijk(5,1,8)= -1.0d0 / 12.0d0
  cijk(5,2,1)= -1.0d0 / 12.0d0
  cijk(5,2,6)=  1.0d0 / 12.0d0
  cijk(5,4,1)=  1.0d0 / 12.0d0
  cijk(5,4,8)= -1.0d0 / 12.0d0
  cijk(5,6,1)= -1.0d0 / 12.0d0
  cijk(5,6,2)= -1.0d0 / 12.0d0
  cijk(5,6,7)=  1.0d0 / 12.0d0
  cijk(5,6,8)=  1.0d0 / 12.0d0
  cijk(5,7,6)= -1.0d0 / 12.0d0
  cijk(5,7,8)=  1.0d0 / 12.0d0
  cijk(5,8,1)=  1.0d0 / 12.0d0
  cijk(5,8,4)=  1.0d0 / 12.0d0
  cijk(5,8,6)= -1.0d0 / 12.0d0
  cijk(5,8,7)= -1.0d0 / 12.0d0
  cijk(6,1,2)=  1.0d0 / 12.0d0
  cijk(6,1,5)= -1.0d0 / 12.0d0
  cijk(6,2,1)= -1.0d0 / 12.0d0
  cijk(6,2,3)=  1.0d0 / 12.0d0
  cijk(6,2,5)= -1.0d0 / 12.0d0
  cijk(6,2,7)=  1.0d0 / 12.0d0
  cijk(6,3,2)= -1.0d0 / 12.0d0
  cijk(6,3,7)=  1.0d0 / 12.0d0
  cijk(6,5,1)=  1.0d0 / 12.0d0
  cijk(6,5,2)=  1.0d0 / 12.0d0
  cijk(6,5,7)= -1.0d0 / 12.0d0
  cijk(6,5,8)= -1.0d0 / 12.0d0
  cijk(6,7,2)= -1.0d0 / 12.0d0
  cijk(6,7,3)= -1.0d0 / 12.0d0
  cijk(6,7,5)=  1.0d0 / 12.0d0
  cijk(6,7,8)=  1.0d0 / 12.0d0
  cijk(6,8,5)=  1.0d0 / 12.0d0
  cijk(6,8,7)= -1.0d0 / 12.0d0
  cijk(7,2,3)=  1.0d0 / 12.0d0
  cijk(7,2,6)= -1.0d0 / 12.0d0
  cijk(7,3,2)= -1.0d0 / 12.0d0
  cijk(7,3,4)=  1.0d0 / 12.0d0
  cijk(7,3,6)= -1.0d0 / 12.0d0
  cijk(7,3,8)=  1.0d0 / 12.0d0
  cijk(7,4,3)= -1.0d0 / 12.0d0
  cijk(7,4,8)=  1.0d0 / 12.0d0
  cijk(7,5,6)=  1.0d0 / 12.0d0
  cijk(7,5,8)= -1.0d0 / 12.0d0
  cijk(7,6,2)=  1.0d0 / 12.0d0
  cijk(7,6,3)=  1.0d0 / 12.0d0
  cijk(7,6,5)= -1.0d0 / 12.0d0
  cijk(7,6,8)= -1.0d0 / 12.0d0
  cijk(7,8,3)= -1.0d0 / 12.0d0
  cijk(7,8,4)= -1.0d0 / 12.0d0
  cijk(7,8,5)=  1.0d0 / 12.0d0
  cijk(7,8,6)=  1.0d0 / 12.0d0
  cijk(8,1,4)= -1.0d0 / 12.0d0
  cijk(8,1,5)=  1.0d0 / 12.0d0
  cijk(8,3,4)=  1.0d0 / 12.0d0
  cijk(8,3,7)= -1.0d0 / 12.0d0
  cijk(8,4,1)=  1.0d0 / 12.0d0
  cijk(8,4,3)= -1.0d0 / 12.0d0
  cijk(8,4,5)=  1.0d0 / 12.0d0
  cijk(8,4,7)= -1.0d0 / 12.0d0
  cijk(8,5,1)= -1.0d0 / 12.0d0
  cijk(8,5,4)= -1.0d0 / 12.0d0
  cijk(8,5,6)=  1.0d0 / 12.0d0
  cijk(8,5,7)=  1.0d0 / 12.0d0
  cijk(8,6,5)= -1.0d0 / 12.0d0
  cijk(8,6,7)=  1.0d0 / 12.0d0
  cijk(8,7,3)=  1.0d0 / 12.0d0
  cijk(8,7,4)=  1.0d0 / 12.0d0
  cijk(8,7,5)= -1.0d0 / 12.0d0
  cijk(8,7,6)= -1.0d0 / 12.0d0



  return
end subroutine getcijk3d8nod

