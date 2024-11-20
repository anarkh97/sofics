C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C*                                                           *
C*  Purpose :  Construct ELement Stiffness Matrix for Beam   *
C*             element.  It stores only the upper triangular *
C*                                                           *
C*  Authors :  E. Pramono & Charbel Farhat                   *
C*                                                           *
C*            CSSC/ BOULDER 1990                             *
C*                                                           *
C*  Modified By : J. C. Chiou and P. Stern, 1/15/91          *
C*                                                           *
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C N O T E !! N O T E !! N O T E !! N O T E !! N O T E !! N O T E !!
C 
C   This routine is exactly the same as mstif6.f with only
C   one small modificatin on arguments.  Combines resmtx and eframe
C   to one argument, still call eframe.  This is to avoid the
C   messy argument passing from Domain class to Element class.
C   Po-Shu Chen Mar. 28 '94
C
C   Change the content of the argument eframe and remove elm.
C   eframe contains only a section, which is for the current element.
C   Po-Shu Chen Mar. 1 '95
C 
C N O T E !! N O T E !! N O T E !! N O T E !! N O T E !! N O T E !!
C
        subroutine modmstif6(area,e,elstif,eframe,ix,iy,iz,nu,
     &                       x,y,z,flg)
C*************************************************************C
C
C.... DECLARE GLOBAL VARIABLES
C
C.... REAL CONSTANTS
C
        real*8  area, e, nu, ix, iy, iz
C
C.... REAL ARRAYS
C
        real*8  eframe(9),elstif(12,12),x(2),y(2),z(2)
C
C.... LOCAL VARIABLE DECLARATIONS
C
        integer ii,jj,flg
C
        real*8 G,J,dx,dy,dz,length,length2,length3
        real*8 EIy, EIz
        real*8 vec(10),ke(12,12)
        real*8 u(3),v(3),w(3)
        real*8 tran(12,12),trant(12,12),ttke(12,12),stif(12,12)
C
C.... INITIALIZE LOCAL CONSTANTS
C
        G = E /(2.0d00*(1.0d00+nu))
        J = Ix
C
C.... COMPUTE THE LENGTH OF THE BEAM
C
        dx = x(2) - x(1)
        dy = y(2) - y(1)
        dz = z(2) - z(1)
        length  = dsqrt(dx*dx + dy*dy + dz*dz)

        if(length .eq. 0.0) then
          write(6,*) 'zero length in (Euler beam) modmstif6.f'
        endif

        length2 = length*length
        length3 = length2*length
        EIy = E*Iy
        EIz = E*Iz
C
        vec(1)   =  E * area / length
        vec(2)   = 12.0d00 * EIy / length3
        vec(3)   = 12.0d00 * EIz / length3
        vec(4)   =  6.0d00 * EIy / length2
        vec(5)   =  6.0d00 * EIz / length2
        vec(6)   =  G  * J  / length
        vec(7)   =  4.0d00 * EIy / length
        vec(8)   =  4.0d00 * EIz / length
        vec(9)   =  2.0d00 * EIy / length
        vec(10)  =  2.0d00 * EIz / length
C
C       Construct the Upper Triangle Stiffness Matrix
C
        call frame6(eframe,u,v,w)
C
        call transf(u,v,w,tran,trant)
C
        do 110 jj = 1, 12
          do 110 ii = 1, 12
            ke(ii,jj) = 0.0d00
 110    continue
C
C      write(6,*) "flg is", flg
C      write(6,*) "entries of stiffness matrix in local frame"
C      write(6,*) "vec(1) = ", vec(1)
C      write(6,*) "vec(2) = ", vec(2)
C      write(6,*) "vec(3) = ", vec(3)
C      write(6,*) "vec(4) = ", vec(4)
C      write(6,*) "vec(5) = ", vec(5)
C      write(6,*) "vec(6) = ", vec(6)
C      write(6,*) "vec(7) = ", vec(7)
C      write(6,*) "vec(8) = ", vec(8)
C      write(6,*) "vec(9) = ", vec(9)
C      write(6,*) "vec(10) = ", vec(10)
C
        ke(1,1)   = vec(1)
        ke(2,2)   = vec(3)
        ke(3,3)   = vec(2)
        ke(4,4)   = vec(6)
        ke(5,5)   = vec(7)
        ke(6,6)   = vec(8)
        ke(7,7)   = vec(1)
        ke(8,8)   = vec(3)
        ke(9,9)   = vec(2)
        ke(10,10) = vec(6)
        ke(11,11) = vec(7)
        ke(12,12) = vec(8)
        ke(5,3)   = -vec(4)
        ke(3,5)   = ke(5,3)
        ke(6,2)   = vec(5)
        ke(2,6)   = ke(6,2)
        ke(7,1)   = -vec(1)
        ke(1,7)   = ke(7,1)
        ke(8,2)   = -vec(3)
        ke(2,8)   = ke(8,2)
        ke(8,6)   = -vec(5)
        ke(6,8)   = ke(8,6)
        ke(9,3)   = -vec(2)
        ke(3,9)   = ke(9,3)
        ke(9,5)   = vec(4)
        ke(5,9)   = ke(9,5)
        ke(10,4)  = -vec(6)
        ke(4,10)  = ke(10,4)
        ke(11,3)  = -vec(4)
        ke(3,11)  = ke(11,3)
        ke(11,5)  = vec(9)
        ke(5,11)  = ke(11,5)
        ke(11,9)  = vec(4)
        ke(9,11)  = ke(11,9)
        ke(12,2)  = vec(5)
        ke(2,12)  = ke(12,2)
        ke(12,6)  = vec(10)
        ke(6,12)  = ke(12,6)
        ke(12,8)  = -vec(5)
        ke(8,12)  = ke(12,8)
C
C dgemm is a BLAS routine
C
        call dgemm('N','N',12,12,12,1.0d00,trant,12,ke,12,0.0d00,ttke,
     +  12)
        call dgemm('N','N',12,12,12,1.0d00,ttke,12,tran,12,0.0d00,stif,
     +  12)
C
C       KHP: MODIFICATION 8/21/98 for Nonlinear Corotational Method
C
C       if flg equals 1, return the transformed element stiffness matrix
C       else, return the local element stiffness matrix
C
        if(flg .eq. 1) then
          do 400 jj = 1, 12
            do 400 ii = 1, 12
               elstif(ii,jj) = stif(ii,jj)
 400      continue
        else
          do jj = 1, 12
            do ii = 1, 12
               elstif(ii,jj) = ke(ii,jj)
            enddo
          enddo
        endif
          
        return
        end
C=END FORTRAN
