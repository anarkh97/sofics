C $Id$

      subroutine cnodetriangle_cpproj(npairs,coords_node,
     *  coords_triangle,rmsnorm,ipushdr,ctrl,ctrcl,tola1,tola2,ctoln2 ) 
C
C Description: This routine computes the closest point projection for a 
C              three node triangle and a node.
C
C Memory:
C   npairs            i/-     number of pairs of faces/nodes to process
C   coords_node       i/-     node coordinates
C   coords_triangle   i/-     coordinates of face nodes
C   rmsnorm           i/-     face normals
C   ipushdr           i/o     pushback direction
C   ctrl              -/-     scratch array
C   ctrcl             -/o     contact information for each pair
C   tola1             i/-
C   tola2             i/-
C   ctoln2            i/-     search tolerance on normal gap

      implicit real*8 (a-h,o-z)
      include "search_parameters.par"

      dimension coords_node(n_v3d,npairs)
      dimension coords_triangle(n_v3d,3,npairs)
      dimension rmsnorm(n_v3d,npairs)
      dimension ipushdr(npairs)
   
      dimension ctrcl(ISIZCTRCL,npairs)
      dimension ctrl(3,8,npairs)

      do 100 i = 1,npairs
        if( ipushdr(i) .eq. inormal )then
          ctrcl(MSPARAM,i) = 0.0
c predicted slave node coordinates
          xsd = coords_node(k_v3dx,i)
          ysd = coords_node(k_v3dy,i)
          zsd = coords_node(k_v3dz,i)
c predicted triangular facet coordinates
          x1d = coords_triangle(k_v3dx,1,i)
          y1d = coords_triangle(K_v3dy,1,i)
          z1d = coords_triangle(k_v3dz,1,i)
          x2d = coords_triangle(k_v3dx,2,i)
          y2d = coords_triangle(k_v3dy,2,i)
          z2d = coords_triangle(k_v3dz,2,i)
          x3d = coords_triangle(k_v3dx,3,i)
          y3d = coords_triangle(k_v3dy,3,i)
          z3d = coords_triangle(k_v3dz,3,i)
C
C compute normal distance from the master surface to the slave node
C note that a negative distance implies the slave node is inside the
C face
          vx1s = xsd - x1d
          vy1s = ysd - y1d
          vz1s = zsd - z1d
C dot the vector from point 1 on master surface to slave node
C with the outward unit normal
          projn = vx1s*rmsnorm(k_v3dx,i) + 
     *            vy1s*rmsnorm(k_v3dy,i) +
     *            vz1s*rmsnorm(k_v3dz,i)
C find contact point
          xc0 = xsd - projn*rmsnorm(k_v3dx,i)
          yc0 = ysd - projn*rmsnorm(k_v3dy,i)
          zc0 = zsd - projn*rmsnorm(k_v3dz,i)
C determine if the contact point is inside the master surface triangle
C sub areas of slave node - to triangle corners
          vx1 = x1d - xc0
          vy1 = y1d - yc0
          vz1 = z1d - zc0
          vx2 = x2d - xc0
          vy2 = y2d - yc0
          vz2 = z2d - zc0
          vx3 = x3d - xc0
          vy3 = y3d - yc0
          vz3 = z3d - zc0
c area in triangle 1: 2-3-c, 2 2: 1-3-c, 3: 1-2-c
c   (Note: Actually 2.0*Area)
          a41 = ( vy2*vz3-vy3*vz2)*rmsnorm(k_v3dx,i) + 
     *          (-vx2*vz3+vx3*vz2)*rmsnorm(k_v3dy,i) + 
     *          ( vx2*vy3-vx3*vy2)*rmsnorm(k_v3dz,i)
          a42 = ( vy3*vz1-vy1*vz3)*rmsnorm(k_v3dx,i) +
     *          (-vx3*vz1+vx1*vz3)*rmsnorm(k_v3dy,i) + 
     *          ( vx3*vy1-vx1*vy3)*rmsnorm(k_v3dz,i)
          a43 = ( vy1*vz2-vy2*vz1)*rmsnorm(k_v3dx,i) + 
     *          (-vx1*vz2+vx2*vz1)*rmsnorm(k_v3dy,i) + 
     *          ( vx1*vy2-vx2*vy1)*rmsnorm(k_v3dz,i) 
          a44 = a41 + a42 + a43
c compute the minimum area and keep the "minimum" node
          if( a41 .lt. a42 )then
            if( a41 .lt. a43 )then
              a44r = a41
              imin_nd = 1
            else
              a44r = a43
              imin_nd = 3
            endif
          else
            if( a42 .lt. a43 )then
              a44r = a42
              imin_nd = 2
            else
              a44r = a43
              imin_nd = 3
            endif
          endif
          if( a44r.ge.-tola1*a44 .and. abs(projn).le.ctoln2 )then 
            ctrcl(MSPARAM,i)  = 1
            ctrcl(ICPOINTX,i) = xc0
            ctrcl(ICPOINTY,i) = yc0
            ctrcl(ICPOINTZ,i) = zc0
            ctrcl(IPENMAG,i)  = projn
            ctrcl(IPUSHX,i) = rmsnorm(k_v3dx,i)
            ctrcl(IPUSHY,i)  = rmsnorm(k_v3dy,i)
            ctrcl(IPUSHZ,i)  = rmsnorm(k_v3dz,i)
            ctrcl(INORMX,i)  = rmsnorm(k_v3dx,i)
            ctrcl(INORMY,i) = rmsnorm(k_v3dy,i)
            ctrcl(INORMZ,i) = rmsnorm(k_v3dz,i)
            ctrcl(ILOCATION,i) = iinside
          else if( abs(projn).le.ctoln2 )then 
c compute a "tangential distance"
c the area of the triangle is 0.5*base*height which can be written as
c     base = tang_dist = 2.0*area/e_len
c (Note: Because we computed 2.0*area above we will drop the 2.0 below)
c
c           compute the edge length
            if( imin_nd .eq. 1)then
              disx = x2d - x3d
              disy = y2d - y3d
              disz = z2d - z3d
              e_len = sqrt( disx*disx + disy*disy + disz*disz )
            else if( imin_nd .eq. 2 )then
              disx = x1d - x3d
              disy = y1d - y3d
              disz = z1d - z3d
              e_len = sqrt( disx*disx + disy*disy + disz*disz )
            else
              disx = x2d - x1d
              disy = y2d - y1d
              disz = z2d - z1d
              e_len = sqrt( disx*disx + disy*disy + disz*disz )
            endif
            tang_dist = abs(a44r)/e_len
            if( tang_dist .le. tola2 )then
              ipushdr(i) = iclose
C  FIX? WHY NOT DO THE PROJECTION TO THE CLOSEST EDGE HERE AND REMOVE THE
C  TWO SUBSEQUENT LOOPS? REJ 12/05
            endif
          endif
        endif
  100 continue


      do 200 i = 1,npairs
C     
C push back along along minimum distance to master surface
C     
        if( ipushdr(i) .eq. iclose )then
C predicted slave node coords
          xsd = coords_node(k_v3dx,i)
          ysd = coords_node(k_v3dy,i)
          zsd = coords_node(k_v3dz,i)
C predicted triangular facet coords
          x1d = coords_triangle(k_v3dx,1,i)
          y1d = coords_triangle(k_v3dy,1,i)
          z1d = coords_triangle(k_v3dz,1,i)
          x2d = coords_triangle(k_v3dx,2,i)
          y2d = coords_triangle(k_v3dy,2,i)
          z2d = coords_triangle(k_v3dz,2,i)
          x3d = coords_triangle(k_v3dx,3,i)
          y3d = coords_triangle(k_v3dy,3,i)
          z3d = coords_triangle(k_v3dz,3,i)
C edge vectors
          vx0   = x3d - x1d
          vy0   = y3d - y1d
          vz0   = z3d - z1d
          vx00  = x2d - x1d
          vy00  = y2d - y1d
          vz00  = z2d - z1d
          vx000 = x3d - x2d
          vy000 = y3d - y2d
          vz000 = z3d - z2d
C store facet normal
          ctrcl(INORMX,i) = rmsnorm(k_v3dx,i)
          ctrcl(INORMY,i) = rmsnorm(k_v3dy,i)
          ctrcl(INORMZ,i) = rmsnorm(k_v3dz,i)
C compute normal distance from the master surface to the slave node
C note that a negative distance implies the slave node is inside the
C face
          vx1s = xsd - x1d
          vy1s = ysd - y1d
          vz1s = zsd - z1d
C dot the vector from point 1 on master surface to slave node
C with the outward unit normal
          projn = vx1s*rmsnorm(k_v3dx,i) + 
     *            vy1s*rmsnorm(k_v3dy,i) + 
     *            vz1s*rmsnorm(k_v3dz,i)
C find contact point
          xc0 = xsd - projn*rmsnorm(k_v3dx,i)
          yc0 = ysd - projn*rmsnorm(k_v3dy,i)
          zc0 = zsd - projn*rmsnorm(k_v3dz,i)
C     
C projection of vector from node 1 to slave node onto the
C vector from node 1 to node 3, and the local coord, s
          proj13 = vx1s*vx0 + vy1s*vy0 + vz1s*vz0

C magnitude of vect. 1-3
          vmag13 = abs(vx0*vx0+vy0*vy0+vz0*vz0)
C local coord. (s) alonge 1-3
          s13 = proj13/vmag13
          s13 = max(zero,s13)
          s13 = min(one,s13)
C nearest point on line 1-3
          xc13 = x1d + s13*vx0
          yc13 = y1d + s13*vy0
          zc13 = z1d + s13*vz0
C vector from nearest point to slave node
          vcs13x = xsd-xc13
          vcs13y = ysd-yc13
          vcs13z = zsd-zc13
C distance from nearest point to slave node
          dis13 = sqrt(vcs13x*vcs13x+vcs13y*vcs13y+vcs13z*vcs13z)
C dot distance with surface normal to determine pos. or neg. dist.
          sdis = vcs13x*rmsnorm(k_v3dx,i) + 
     *           vcs13y*rmsnorm(k_v3dy,i) +
     *           vcs13z*rmsnorm(k_v3dz,i)
          sdis13 = sign(dis13,sdis)
          ctrl(1,1,i) = 1
          ctrl(1,2,i) = xc13
          ctrl(1,3,i) = yc13
          ctrl(1,4,i) = zc13
          ctrl(1,5,i) = sdis13
          ctrl(1,6,i) = vcs13x/sdis13
          ctrl(1,7,i) = vcs13y/sdis13
          ctrl(1,8,i) = vcs13z/sdis13
C     
C projection of vector from node 1 to slave node onto the
C vector from node 1 to node 2, and the local coord, s
          proj12 = vx1s*vx00 + vy1s*vy00 + vz1s*vz00
C magnitude of vect. 1-2
          vmag12 = abs(vx00*vx00+vy00*vy00+vz00*vz00)
C local coord. (s) alonge 1-2
          s12 = proj12/vmag12
          s12 = max(zero,s12)
          s12 = min(one,s12)
C nearest point on line 1-2
          xc12 = x1d + s12*vx00
          yc12 = y1d + s12*vy00
          zc12 = z1d + s12*vz00
C vector from nearest point to slave node
          vcs12x = xsd-xc12
          vcs12y = ysd-yc12
          vcs12z = zsd-zc12
C distance from nearest point to slave node
          dis12 = sqrt(vcs12x*vcs12x+vcs12y*vcs12y+vcs12z*vcs12z)
C dot distance with surface normal to determine pos. or neg. dist.
          sdis = vcs12x*rmsnorm(k_v3dx,i) + 
     *           vcs12y*rmsnorm(k_v3dy,i) + 
     *           vcs12z*rmsnorm(k_v3dz,i)
          sdis12 = sign(dis12,sdis)
          ctrl(2,1,i) = 1
          ctrl(2,2,i) = xc12
          ctrl(2,3,i) = yc12
          ctrl(2,4,i) = zc12
          ctrl(2,5,i) = sdis12
          ctrl(2,6,i) = vcs12x/sdis12
          ctrl(2,7,i) = vcs12y/sdis12
          ctrl(2,8,i) = vcs12z/sdis12
C
C projection of vector from node 2 to slave node onto the
C vector from node 2 to node 3, and the local coord, s
          vx1s = xsd - x2d
          vy1s = ysd - y2d
          vz1s = zsd - z2d
          proj23 = vx1s*vx000 + vy1s*vy000 + vz1s*vz000
          vmag23 = abs(vx000*vx000+vy000*vy000+vz000*vz000)
C local coord. (s) alonge 2-3
          s23 = proj23/vmag23
          s23 = max(zero,s23)
          s23 = min(one,s23)
C nearest point on line 2-3
          xc23 = x2d + s23*vx000
          yc23 = y2d + s23*vy000
          zc23 = z2d + s23*vz000
C vector from nearest point to slave node
          vcs23x = xsd-xc23
          vcs23y = ysd-yc23
          vcs23z = zsd-zc23
C distance from nearest point to slave node
          dis23 = sqrt(vcs23x*vcs23x+vcs23y*vcs23y+vcs23z*vcs23z)
C dot distance with surface normal to determine pos. or neg. dist.
          sdis = vcs23x*rmsnorm(k_v3dx,i) + 
     *           vcs23y*rmsnorm(k_v3dy,i) +
     *           vcs23z*rmsnorm(k_v3dz,i)
          sdis23 = sign(dis23,sdis)
          ctrl(3,1,i) = 1
          ctrl(3,2,i) = xc23
          ctrl(3,3,i) = yc23
          ctrl(3,4,i) = zc23
          ctrl(3,5,i) = sdis23
          ctrl(3,6,i) = vcs23x/sdis23
          ctrl(3,7,i) = vcs23y/sdis23
          ctrl(3,8,i) = vcs23z/sdis23
        endif
  200 continue

      do 300 i = 1,npairs
        
        if(ipushdr(i) .eq. iclose)then
          
          ismall = 1
          clmin = abs(ctrl(1,5,i))
          if(abs(ctrl(2,5,i)) .lt. clmin)then
            ismall = 2
            clmin = abs(ctrl(2,5,i))
          endif
          if(abs(ctrl(3,5,i)) .lt. clmin)then
            ismall = 3
            clmin = abs(ctrl(3,5,i))
          endif
          
          ctrcl(MSPARAM,i)   = ctrl(ismall,MSPARAM,i)
          ctrcl(ICPOINTX,i)  = ctrl(ismall,ICPOINTX,i)
          ctrcl(ICPOINTY,i)  = ctrl(ismall,ICPOINTY,i)
          ctrcl(ICPOINTZ,i)  = ctrl(ismall,ICPOINTZ,i)
          ctrcl(IPENMAG,i)   = ctrl(ismall,IPENMAG,i)
          ctrcl(IPUSHX,i)    = ctrl(ismall,IPUSHX,i)
          ctrcl(IPUSHY,i)    = ctrl(ismall,IPUSHY,i)
          ctrcl(IPUSHZ,i)    = ctrl(ismall,IPUSHZ,i)
          ctrcl(ILOCATION,i) = iout

        endif
  300 continue
C     
      return
      end

      subroutine cnodetriangle_cpproj_aug(npairs,
     *   coords_node_pre, coords_triangle_pre, rmsnorm_pre,
     *   coords_node_aug, coords_triangle_aug, rmsnorm_aug,
     *   ipushdr,ctrl,ctrcl,tola1,tola2,ctoln2 ) 
C
C Description: This routine computes the closest point projection for a 
C              three node triangle and a node.
C
C Memory:
C   npairs            i/-     number of pairs of faces/nodes to process
C   coords_node_pre       i/-     node coordinates
C   coords_triangle_pre   i/-     coordinates of face nodes
C   rmsnorm_pre           i/-     face normals
C   coords_node_aug      i/-     node coordinates in augmented configuration
C   coords_triangle_aug  i/-     coordinates of face nodes in augmented configuration
C   rmsnorm_aug          i/-     face normals in augmented configuration
C   ipushdr           i/o     pushback direction
C   ctrl              -/-     scratch array
C   ctrcl             -/o
C   tola1             i/-
C   tola2             i/-
C   ctoln2            i/-

      implicit real*8 (a-h,o-z)
      include "search_parameters.par"

      dimension coords_node_pre(n_v3d,npairs)
      dimension coords_triangle_pre(n_v3d,3,npairs)
      dimension rmsnorm_pre(n_v3d,npairs)
      dimension coords_node_aug(n_v3d,npairs)
      dimension coords_triangle_aug(n_v3d,3,npairs)
      dimension rmsnorm_aug(n_v3d,npairs)
      dimension ipushdr(npairs)
   
      dimension ctrcl(ISIZCTRCL,npairs)
      dimension ctrl(3,8,npairs)

      do 100 i = 1,npairs
        if( ipushdr(i) .eq. inormal )then
          ctrcl(MSPARAM,i) = 0.0
c Begin with a check for normal distance in the predicted configuration
c   predicted slave node coordinates
          xsd = coords_node_pre(k_v3dx,i)
          ysd = coords_node_pre(k_v3dy,i)
          zsd = coords_node_pre(k_v3dz,i)
c   predicted triangular facet coordinates
          x1d = coords_triangle_pre(k_v3dx,1,i)
          y1d = coords_triangle_pre(K_v3dy,1,i)
          z1d = coords_triangle_pre(k_v3dz,1,i)
C  compute normal distance from the master surface to the slave node
          vx1s = xsd - x1d
          vy1s = ysd - y1d
          vz1s = zsd - z1d
C dot the vector from point 1 on master surface to slave node
C with the outward unit normal
          projn = vx1s*rmsnorm_aug(k_v3dx,i) + 
     *            vy1s*rmsnorm_aug(k_v3dy,i) +
     *            vz1s*rmsnorm_aug(k_v3dz,i)
c if projn > ctoln2 the node is too far outside so we don't need to continue
          if( projn .le. ctoln2 )then
c augmented slave node coordinates
            xsd = coords_node_aug(k_v3dx,i)
            ysd = coords_node_aug(k_v3dy,i)
            zsd = coords_node_aug(k_v3dz,i)
c augmented triangular facet coordinates
            x1d = coords_triangle_aug(k_v3dx,1,i)
            y1d = coords_triangle_aug(K_v3dy,1,i)
            z1d = coords_triangle_aug(k_v3dz,1,i)
            x2d = coords_triangle_aug(k_v3dx,2,i)
            y2d = coords_triangle_aug(k_v3dy,2,i)
            z2d = coords_triangle_aug(k_v3dz,2,i)
            x3d = coords_triangle_aug(k_v3dx,3,i)
            y3d = coords_triangle_aug(k_v3dy,3,i)
            z3d = coords_triangle_aug(k_v3dz,3,i)
C
C  compute normal distance from the master surface to the slave node
C  note that a negative distance implies the slave node is inside the
C face
            vx1s = xsd - x1d
            vy1s = ysd - y1d
            vz1s = zsd - z1d
C dot the vector from point 1 on master surface to slave node
C with the outward unit normal
            projn = vx1s*rmsnorm_aug(k_v3dx,i) + 
     *              vy1s*rmsnorm_aug(k_v3dy,i) +
     *              vz1s*rmsnorm_aug(k_v3dz,i)
C find contact point
            xc0 = xsd - projn*rmsnorm_aug(k_v3dx,i)
            yc0 = ysd - projn*rmsnorm_aug(k_v3dy,i)
            zc0 = zsd - projn*rmsnorm_aug(k_v3dz,i)
C determine if the contact point is inside the master surface triangle
C sub areas of slave node - to triangle corners
            vx1 = x1d - xc0
            vy1 = y1d - yc0
            vz1 = z1d - zc0
            vx2 = x2d - xc0
            vy2 = y2d - yc0
            vz2 = z2d - zc0
            vx3 = x3d - xc0
            vy3 = y3d - yc0
            vz3 = z3d - zc0
c area in triangle 1: 2-3-c, 2 2: 1-3-c, 3: 1-2-c
c   (Note: Actually 2.0*Area)
            a41 = ( vy2*vz3-vy3*vz2)*rmsnorm_aug(k_v3dx,i) + 
     *            (-vx2*vz3+vx3*vz2)*rmsnorm_aug(k_v3dy,i) + 
     *            ( vx2*vy3-vx3*vy2)*rmsnorm_aug(k_v3dz,i)
            a42 = ( vy3*vz1-vy1*vz3)*rmsnorm_aug(k_v3dx,i) +
     *            (-vx3*vz1+vx1*vz3)*rmsnorm_aug(k_v3dy,i) + 
     *            ( vx3*vy1-vx1*vy3)*rmsnorm_aug(k_v3dz,i)
            a43 = ( vy1*vz2-vy2*vz1)*rmsnorm_aug(k_v3dx,i) + 
     *            (-vx1*vz2+vx2*vz1)*rmsnorm_aug(k_v3dy,i) + 
     *            ( vx1*vy2-vx2*vy1)*rmsnorm_aug(k_v3dz,i) 
            a44 = a41 + a42 + a43
            if( a41 .lt. a42 )then
              if( a41 .lt. a43 )then
                a44r = a41
                imin_nd = 1
              else
                a44r = a43
                imin_nd = 3
              endif
            else
              if( a42 .lt. a43 )then
                a44r = a42
                imin_nd = 2
              else
                a44r = a43
                imin_nd = 3
              endif
            endif

            if( a44r.ge.-tola1*a44 )then 
c re-compute contact gap and pushback based on predicted coords
c   predicted slave node coordinates
              xsd = coords_node_pre(k_v3dx,i)
              ysd = coords_node_pre(k_v3dy,i)
              zsd = coords_node_pre(k_v3dz,i)
c   predicted triangular facet coordinates
              x1d = coords_triangle_pre(k_v3dx,1,i)
              y1d = coords_triangle_pre(K_v3dy,1,i)
              z1d = coords_triangle_pre(k_v3dz,1,i)
              x2d = coords_triangle_pre(k_v3dx,2,i)
              y2d = coords_triangle_pre(k_v3dy,2,i)
              z2d = coords_triangle_pre(k_v3dz,2,i)
              x3d = coords_triangle_pre(k_v3dx,3,i)
              y3d = coords_triangle_pre(k_v3dy,3,i)
              z3d = coords_triangle_pre(k_v3dz,3,i)
c   compute global contact coordinates in the predicted configuration
c   assuming the local coords ( a41/a44 , a42/a44 , a43/a44 )
              xc0 = ( a41*x1d + a42*x2d + a43*x3d ) / a44
              yc0 = ( a41*y1d + a42*y2d + a43*y3d ) / a44
              zc0 = ( a41*z1d + a42*z2d + a43*z3d ) / a44
C   compute distance from the master surface to the slave node
              vx1s = xsd - xc0
              vy1s = ysd - yc0
              vz1s = zsd - zc0
c   re-compute gap & pushback dir
              projn = vx1s*vx1s + vy1s*vy1s + vz1s*vz1s
              if( projn .gt. 0 )then
                projn = sqrt( projn )           
                px =  vx1s/projn
                py =  vy1s/projn
                pz =  vz1s/projn
              else
                px = rmsnorm_aug(k_v3dx,i)
                py = rmsnorm_aug(k_v3dy,i)
                pz = rmsnorm_aug(k_v3dz,i)
              end if
c   follow correct sign convention
              pdn = px*rmsnorm_aug(k_v3dx,i) +
     *              py*rmsnorm_aug(k_v3dy,i) +
     *              pz*rmsnorm_aug(k_v3dz,i)
              if( pdn .lt. 0 )then
                projn = -projn
                px = -px 
                py = -py
                pz = -pz
              end if
C This check was removed to get scaling tolerances to work
              if( abs(projn).le.ctoln2 )then 
                ctrcl(MSPARAM,i)  = 1
                ctrcl(ICPOINTX,i) = xc0
                ctrcl(ICPOINTY,i) = yc0
                ctrcl(ICPOINTZ,i) = zc0
                ctrcl(IPENMAG,i)  = projn
                ctrcl(IPUSHX,i)   = px
                ctrcl(IPUSHY,i)   = py
                ctrcl(IPUSHZ,i)   = pz
                ctrcl(INORMX,i)   = rmsnorm_aug(k_v3dx,i)
                ctrcl(INORMY,i)   = rmsnorm_aug(k_v3dy,i)
                ctrcl(INORMZ,i)   = rmsnorm_aug(k_v3dz,i)
                ctrcl(ILOCATION,i) = iinside
              end if
            else if( abs(projn).le.ctoln2 )then 
c compute a "tangential distance"
c the area of the triangle is 0.5*base*height which can be written as
c     base = tang_dist = 2.0*area/e_len
c (Note: Because we computed 2.0*area above we will drop the 2.0 below)
c
c             compute the edge length
              if( imin_nd .eq. 1)then
                disx = x2d - x3d
                disy = y2d - y3d
                disz = z2d - z3d
                e_len = sqrt( disx*disx + disy*disy + disz*disz )
              else if( imin_nd .eq. 2 )then
                disx = x1d - x3d
                disy = y1d - y3d
                disz = z1d - z3d
                e_len = sqrt( disx*disx + disy*disy + disz*disz )
              else
                disx = x2d - x1d
                disy = y2d - y1d
                disz = z2d - z1d
                e_len = sqrt( disx*disx + disy*disy + disz*disz )
              endif
              tang_dist = abs(a44r)/e_len
              if( tang_dist .le. tola2 )then
                ipushdr(i) = iclose
c               write(*,*) "<> in cnodetriangle_cpproj: ",
c    &          "tang. dist. test triggered ", i
              endif
            endif
          endif
        endif
  100 continue


      do 200 i = 1,npairs
C     
C push back along along minimum distance to master surface
C     
        if( ipushdr(i) .eq. iclose )then
C augmented slave node coords
          xsd = coords_node_aug(k_v3dx,i)
          ysd = coords_node_aug(k_v3dy,i)
          zsd = coords_node_aug(k_v3dz,i)
C augmented triangular facet coords
          x1d = coords_triangle_aug(k_v3dx,1,i)
          y1d = coords_triangle_aug(k_v3dy,1,i)
          z1d = coords_triangle_aug(k_v3dz,1,i)
          x2d = coords_triangle_aug(k_v3dx,2,i)
          y2d = coords_triangle_aug(k_v3dy,2,i)
          z2d = coords_triangle_aug(k_v3dz,2,i)
          x3d = coords_triangle_aug(k_v3dx,3,i)
          y3d = coords_triangle_aug(k_v3dy,3,i)
          z3d = coords_triangle_aug(k_v3dz,3,i)
C edge vectors
          vx0   = x3d - x1d
          vy0   = y3d - y1d
          vz0   = z3d - z1d
          vx00  = x2d - x1d
          vy00  = y2d - y1d
          vz00  = z2d - z1d
          vx000 = x3d - x2d
          vy000 = y3d - y2d
          vz000 = z3d - z2d
C store facet normal in the predicted config
          ctrcl(INORMX,i) = rmsnorm_aug(k_v3dx,i)
          ctrcl(INORMY,i) = rmsnorm_aug(k_v3dy,i)
          ctrcl(INORMZ,i) = rmsnorm_aug(k_v3dz,i)
          px = rmsnorm_aug(k_v3dx,i)
          py = rmsnorm_aug(k_v3dy,i)
          pz = rmsnorm_aug(k_v3dz,i)
C compute normal distance from the master surface to the slave node
C note that a negative distance implies the slave node is inside the
C face
          vx1s = xsd - x1d
          vy1s = ysd - y1d
          vz1s = zsd - z1d
C dot the vector from point 1 on master surface to slave node
C with the outward unit normal
          projn = vx1s*rmsnorm_aug(k_v3dx,i) + 
     *            vy1s*rmsnorm_aug(k_v3dy,i) + 
     *            vz1s*rmsnorm_aug(k_v3dz,i)
C find contact point
          xc0 = xsd - projn*rmsnorm_aug(k_v3dx,i)
          yc0 = ysd - projn*rmsnorm_aug(k_v3dy,i)
          zc0 = zsd - projn*rmsnorm_aug(k_v3dz,i)
C     
C projection of vector from node 1 to slave node onto the
C vector from node 1 to node 3, and the local coord, s
          proj13 = vx1s*vx0 + vy1s*vy0 + vz1s*vz0

C magnitude of vect. 1-3
          vmag13 = abs(vx0*vx0+vy0*vy0+vz0*vz0)
C local coord. (s) along 1-3
          s13 = proj13/vmag13
          s13 = max(zero,s13)
          s13 = min(one,s13)
C nearest point on line 1-3
          xc13 = x1d + s13*vx0
          yc13 = y1d + s13*vy0
          zc13 = z1d + s13*vz0
C vector from nearest point to slave node
          vcs13x = xsd-xc13
          vcs13y = ysd-yc13
          vcs13z = zsd-zc13
C distance from nearest point to slave node
          dis13 = sqrt(vcs13x*vcs13x+vcs13y*vcs13y+vcs13z*vcs13z)
C dot distance with surface normal to determine pos. or neg. dist.
          sdis = vcs13x*rmsnorm_aug(k_v3dx,i) + 
     *           vcs13y*rmsnorm_aug(k_v3dy,i) +
     *           vcs13z*rmsnorm_aug(k_v3dz,i)
          sdis13 = sign(dis13,sdis)
          if (sdis13.gt.0.0) then
            px = vcs13x/sdis13
            py = vcs13y/sdis13
            pz = vcs13z/sdis13
          endif
          ctrl(1,1,i) = s13
          ctrl(1,2,i) = xc13
          ctrl(1,3,i) = yc13
          ctrl(1,4,i) = zc13
          ctrl(1,5,i) = sdis13
          ctrl(1,6,i) = px
          ctrl(1,7,i) = py
          ctrl(1,8,i) = pz
C     
C projection of vector from node 1 to slave node onto the
C vector from node 1 to node 2, and the local coord, s
          proj12 = vx1s*vx00 + vy1s*vy00 + vz1s*vz00
C magnitude of vect. 1-2
          vmag12 = abs(vx00*vx00+vy00*vy00+vz00*vz00)
C local coord. (s) alonge 1-2
          s12 = proj12/vmag12
          s12 = max(zero,s12)
          s12 = min(one,s12)
C nearest point on line 1-2
          xc12 = x1d + s12*vx00
          yc12 = y1d + s12*vy00
          zc12 = z1d + s12*vz00
C vector from nearest point to slave node
          vcs12x = xsd-xc12
          vcs12y = ysd-yc12
          vcs12z = zsd-zc12
C distance from nearest point to slave node
          dis12 = sqrt(vcs12x*vcs12x+vcs12y*vcs12y+vcs12z*vcs12z)
C dot distance with surface normal to determine pos. or neg. dist.
          sdis = vcs12x*rmsnorm_aug(k_v3dx,i) + 
     *           vcs12y*rmsnorm_aug(k_v3dy,i) + 
     *           vcs12z*rmsnorm_aug(k_v3dz,i)
          sdis12 = sign(dis12,sdis)
          if (sdis12.gt.0.0) then
            px = vcs12x/sdis12
            py = vcs12y/sdis12
            pz = vcs12z/sdis12
          endif
          ctrl(2,1,i) = s12
          ctrl(2,2,i) = xc12
          ctrl(2,3,i) = yc12
          ctrl(2,4,i) = zc12
          ctrl(2,5,i) = sdis12
          ctrl(2,6,i) = px
          ctrl(2,7,i) = py
          ctrl(2,8,i) = pz
C
C projection of vector from node 2 to slave node onto the
C vector from node 2 to node 3, and the local coord, s
          vx1s = xsd - x2d
          vy1s = ysd - y2d
          vz1s = zsd - z2d
          proj23 = vx1s*vx000 + vy1s*vy000 + vz1s*vz000
          vmag23 = abs(vx000*vx000+vy000*vy000+vz000*vz000)
C local coord. (s) alonge 2-3
          s23 = proj23/vmag23
          s23 = max(zero,s23)
          s23 = min(one,s23)
C nearest point on line 2-3
          xc23 = x2d + s23*vx000
          yc23 = y2d + s23*vy000
          zc23 = z2d + s23*vz000
C vector from nearest point to slave node
          vcs23x = xsd-xc23
          vcs23y = ysd-yc23
          vcs23z = zsd-zc23
C distance from nearest point to slave node
          dis23 = sqrt(vcs23x*vcs23x+vcs23y*vcs23y+vcs23z*vcs23z)
C dot distance with surface normal to determine pos. or neg. dist.
          sdis = vcs23x*rmsnorm_aug(k_v3dx,i) + 
     *           vcs23y*rmsnorm_aug(k_v3dy,i) +
     *           vcs23z*rmsnorm_aug(k_v3dz,i)
          sdis23 = sign(dis23,sdis)
          if (sdis23.gt.0.0) then
            px = vcs23x/sdis23
            py = vcs23y/sdis23
            pz = vcs23z/sdis23
          endif
          ctrl(3,1,i) = s23
          ctrl(3,2,i) = xc23
          ctrl(3,3,i) = yc23
          ctrl(3,4,i) = zc23
          ctrl(3,5,i) = sdis23
          ctrl(3,6,i) = px
          ctrl(3,7,i) = py
          ctrl(3,8,i) = pz
        endif
  200 continue

      do 300 i = 1,npairs
        
        if(ipushdr(i) .eq. iclose)then
          
          ismall = 1
          clmin = abs(ctrl(1,5,i))
          if(abs(ctrl(2,5,i)) .lt. clmin)then
            ismall = 2
            clmin = abs(ctrl(2,5,i))
          endif
          if(abs(ctrl(3,5,i)) .lt. clmin)then
            ismall = 3
            clmin = abs(ctrl(3,5,i))
          endif

C predicted slave node coords
          xsd = coords_node_pre(k_v3dx,i)
          ysd = coords_node_pre(k_v3dy,i)
          zsd = coords_node_pre(k_v3dz,i)
C predicted triangular facet coords
          x1d = coords_triangle_pre(k_v3dx,1,i)
          y1d = coords_triangle_pre(k_v3dy,1,i)
          z1d = coords_triangle_pre(k_v3dz,1,i)
          x2d = coords_triangle_pre(k_v3dx,2,i)
          y2d = coords_triangle_pre(k_v3dy,2,i)
          z2d = coords_triangle_pre(k_v3dz,2,i)
          x3d = coords_triangle_pre(k_v3dx,3,i)
          y3d = coords_triangle_pre(k_v3dy,3,i)
          z3d = coords_triangle_pre(k_v3dz,3,i)
C retrieve local coord along edge
          s = ctrl(ismall,1,i)

          px = rmsnorm_pre(k_v3dx,i)
          py = rmsnorm_pre(k_v3dy,i)
          pz = rmsnorm_pre(k_v3dz,i)
          if( ismall .eq. 1 )then
            xc13 = x1d - s*x1d + s*x3d
            yc13 = y1d - s*y1d + s*y3d
            zc13 = z1d - s*z1d + s*z3d
            vcs13x = xsd-xc13
            vcs13y = ysd-yc13
            vcs13z = zsd-zc13
            dis13 = sqrt(vcs13x*vcs13x+vcs13y*vcs13y+vcs13z*vcs13z)
            sdis = vcs13x*rmsnorm_pre(k_v3dx,i) + 
     *             vcs13y*rmsnorm_pre(k_v3dy,i) +
     *             vcs13z*rmsnorm_pre(k_v3dz,i)
            sdis13 = sign(dis13,sdis)
            if (sdis13.gt.0.0) then
              px = vcs13x/sdis13
              py = vcs13y/sdis13
              pz = vcs13z/sdis13
            endif
            ctrcl(MSPARAM,i)   = 1
            ctrcl(ICPOINTX,i)  = xc13
            ctrcl(ICPOINTY,i)  = yc13
            ctrcl(ICPOINTZ,i)  = zc13
            ctrcl(IPENMAG,i)   = sdis13
            ctrcl(IPUSHX,i)    = px
            ctrcl(IPUSHY,i)    = py
            ctrcl(IPUSHZ,i)    = pz
            ctrcl(ILOCATION,i) = iout
          else if( ismall .eq. 2 )then
            xc12 = x1d - s*x1d + s*x2d
            yc12 = y1d - s*y1d + s*y2d
            zc12 = z1d - s*z1d + s*z2d
            vcs12x = xsd-xc12
            vcs12y = ysd-yc12
            vcs12z = zsd-zc12
            dis12 = sqrt(vcs12x*vcs12x+vcs12y*vcs12y+vcs12z*vcs12z)
            sdis = vcs12x*rmsnorm_pre(k_v3dx,i) + 
     *             vcs12y*rmsnorm_pre(k_v3dy,i) + 
     *             vcs12z*rmsnorm_pre(k_v3dz,i)
            sdis12 = sign(dis12,sdis)
            if (sdis12.gt.0.0) then
              px = vcs12x/sdis12
              py = vcs12y/sdis12
              pz = vcs12z/sdis12
            endif
            ctrcl(MSPARAM,i)   = 1
            ctrcl(ICPOINTX,i)  = xc12
            ctrcl(ICPOINTY,i)  = yc12
            ctrcl(ICPOINTZ,i)  = zc12
            ctrcl(IPENMAG,i)   = sdis12
            ctrcl(IPUSHX,i)    = px
            ctrcl(IPUSHY,i)    = py
            ctrcl(IPUSHZ,i)    = pz
            ctrcl(ILOCATION,i) = iout
          else
            xc23 = x2d - s*x2d + s*x3d
            yc23 = y2d - s*y2d + s*y3d
            zc23 = z2d - s*z2d + s*z3d
            vcs23x = xsd-xc23
            vcs23y = ysd-yc23
            vcs23z = zsd-zc23
            dis23 = sqrt(vcs23x*vcs23x+vcs23y*vcs23y+vcs23z*vcs23z)
            sdis = vcs23x*rmsnorm_pre(k_v3dx,i) + 
     *             vcs23y*rmsnorm_pre(k_v3dy,i) +
     *             vcs23z*rmsnorm_pre(k_v3dz,i)
            sdis23 = sign(dis23,sdis)
            if (sdis23.gt.0.0) then
              px = vcs23x/sdis23
              py = vcs23y/sdis23
              pz = vcs23z/sdis23
            endif
            ctrcl(MSPARAM,i)   = 1
            ctrcl(ICPOINTX,i)  = xc23
            ctrcl(ICPOINTY,i)  = yc23
            ctrcl(ICPOINTZ,i)  = zc23
            ctrcl(IPENMAG,i)   = sdis23
            ctrcl(IPUSHX,i)    = px
            ctrcl(IPUSHY,i)    = py
            ctrcl(IPUSHZ,i)    = pz
            ctrcl(ILOCATION,i) = iout
          end if
        endif
  300 continue
C     
      return
      end


      subroutine cnodetriangle_movsrch_aug( npairs,
     *   coords_node_cur, coords_triangle_cur, rmsnorm_cur,
     *   coords_node_pre, coords_triangle_pre, rmsnorm_pre,
     *   coords_node_aug, coords_triangle_aug, rmsnorm_aug,
     *   ctrcl, ctoln2, ctolt1, ctolt2 )
C     
C***********************************************************************
C     
C Description:
C This subroutine calculates the dynamic intersection of a triangular
c face and a node.  The routine will return the status of the search in
c the MSPARAM location of ctrcl.  The possibilities are
c  
c   ctrcl(MSPARAM,i) = -1    Needs a static check
c   ctrcl(MSPARAM,i) = 0     No intersection
c   ctrcl(MSPARAM,i) = 1     Contact point defined 
C
C-----------------------------------------
C
C Arguments
C 
C  INPUT:
C   npairs                i/-    number of pairs in this workset
c   coords_node_cur       i/-    current configuration of the slave nodes
c   coords_node_pre       i/-    predicted configuration of the slave nodes
c   coords_triangle_cur   i/-    current configuration of the master nodes
c   coords_triangle_pre   i/-    predicted configuration of the master nodes
c   rmsnorm_cur           i/-    face normals for current configuration
c   rmsnorm_pre           i/-    face normals for predicted configuration
c   ctoln2                i/-    search tolerance on normal gap
c   ctolt2                i/-    serach toleracne on tangential gap
c
c  INPUT/OUTPUT:
c   ctrcl                 -/o    contact information for each pair
C
C***********************************************************************
C     
c
c

      implicit real*8 (a-h,o-z)
      include "search_parameters.par"

      dimension ctrcl(isizctrcl,npairs)

      dimension coords_node_cur(n_v3d,npairs)
      dimension coords_triangle_cur(n_v3d,3,npairs)
      dimension rmsnorm_cur(n_v3d,npairs)

      dimension coords_node_pre(n_v3d,npairs)
      dimension coords_triangle_pre(n_v3d,3,npairs)
      dimension rmsnorm_pre(n_v3d,npairs)

      dimension coords_node_aug(n_v3d,npairs)
      dimension coords_triangle_aug(n_v3d,3,npairs)
      dimension rmsnorm_aug(n_v3d,npairs)

      do 100 i = 1,npairs
c
c compute normal distance from the MS to SN in current configuration
c    NOTE: a negative distance implies the slave node is penetrating the face
        vfms0x = coords_node_cur(k_v3dx,i) -
     *           coords_triangle_cur(k_v3dx,3,i)
        vfms0y = coords_node_cur(k_v3dy,i) -
     *           coords_triangle_cur(k_v3dy,3,i)
        vfms0z = coords_node_cur(k_v3dz,i) -
     *           coords_triangle_cur(k_v3dz,3,i)
        dmag0  = vfms0x *rmsnorm_cur(k_v3dx,i) + 
     *           vfms0y *rmsnorm_cur(k_v3dy,i) + 
     *           vfms0z *rmsnorm_cur(k_v3dz,i)
c
c compute normal distance from the MS to SN in augmented configuration
c    NOTE: a negative distance implies the slave node is penetrating the face
        vfmspx = coords_node_aug(k_v3dx,i) -
     *           coords_triangle_aug(k_v3dx,3,i)
        vfmspy = coords_node_aug(k_v3dy,i) -
     *           coords_triangle_aug(k_v3dy,3,i)
        vfmspz = coords_node_aug(k_v3dz,i) -
     *           coords_triangle_aug(k_v3dz,3,i)
        dmagp  = vfmspx*rmsnorm_aug(k_v3dx,i) + 
     *           vfmspy*rmsnorm_aug(k_v3dy,i) + 
     *           vfmspz*rmsnorm_aug(k_v3dz,i)
        delmag = dmagp - dmag0
        delmaga = abs(delmag)
        if( delmaga .le. 1.e-10 )then
          rel_mot = 0.0
        else
          rel_mot = delmaga/max(abs(dmagp),abs(dmag0))
        endif
c if we are outside the search tolerance in both configurations, they don't
c interact.
        if( ((dmag0 .gt. ctoln2) .and. (dmagp .gt. ctoln2)) .or.
c if we are outside in both and inside in initial they are moving apart 
     *      ((dmag0 .ge. 0.0) .and. (dmagp .gt. ctoln2)) ) then
           msidf = 0
c if dmag0 & magp are the same (i.e., they aren't moving relative to
c one another, send to CPP
        else if( rel_mot .le. 1.e-8 )then
           msidf = -1
c dynamic intersection
        else
c compute normalized time to contact
        dtcn = -dmag0/delmag
c Compute the location of the slave node at contact time
        xc = coords_node_cur(k_v3dx,i) + 
     *       (coords_node_aug(k_v3dx,i)-coords_node_cur(k_v3dx,i))*dtcn
        yc = coords_node_cur(k_v3dy,i) +  
     *       (coords_node_aug(k_v3dy,i)-coords_node_cur(k_v3dy,i))*dtcn
        zc = coords_node_cur(k_v3dz,i) +  
     *       (coords_node_aug(k_v3dz,i)-coords_node_cur(k_v3dz,i))*dtcn
c position of master surface at contact
        x1c = coords_triangle_cur(k_v3dx,1,i) + dtcn*
     *        (coords_triangle_aug(k_v3dx,1,i)-
     *         coords_triangle_cur(k_v3dx,1,i))
        y1c = coords_triangle_cur(k_v3dy,1,i) + dtcn*
     *        (coords_triangle_aug(k_v3dy,1,i)-
     *         coords_triangle_cur(k_v3dy,1,i))
        z1c = coords_triangle_cur(k_v3dz,1,i) + dtcn*
     *        (coords_triangle_aug(k_v3dz,1,i)-
     *         coords_triangle_cur(k_v3dz,1,i))
        x2c = coords_triangle_cur(k_v3dx,2,i) + dtcn*
     *        (coords_triangle_aug(k_v3dx,2,i)-
     *         coords_triangle_cur(k_v3dx,2,i))
        y2c = coords_triangle_cur(k_v3dy,2,i) + dtcn*
     *        (coords_triangle_aug(k_v3dy,2,i)-
     *         coords_triangle_cur(k_v3dy,2,i))
        z2c = coords_triangle_cur(k_v3dz,2,i) + dtcn*
     *        (coords_triangle_aug(k_v3dz,2,i)-
     *         coords_triangle_cur(k_v3dz,2,i))
        x3c = coords_triangle_cur(k_v3dx,3,i) + dtcn*
     *        (coords_triangle_aug(k_v3dx,3,i)-
     *         coords_triangle_cur(k_v3dx,3,i))
        y3c = coords_triangle_cur(k_v3dy,3,i) + dtcn*
     *        (coords_triangle_aug(k_v3dy,3,i)-
     *         coords_triangle_cur(k_v3dy,3,i))
        z3c = coords_triangle_cur(k_v3dz,3,i) + dtcn*
     *        (coords_triangle_aug(k_v3dz,3,i)-
     *         coords_triangle_cur(k_v3dz,3,i))
c compute the unit normal at the contact time
        vxp0  = x3c - x1c
        vyp0  = y3c - y1c
        vzp0  = z3c - z1c
        vxp00 = x2c - x1c
        vyp00 = y2c - y1c
        vzp00 = z2c - z1c
        srfc1 = -vyp0*vzp00+vyp00*vzp0
        srfc2 =  vxp0*vzp00-vxp00*vzp0
        srfc3 = -vxp0*vyp00+vxp00*vyp0
        srfcn = sqrt(srfc1*srfc1+srfc2*srfc2+srfc3*srfc3)
        srfc1 = srfc1/srfcn
        srfc2 = srfc2/srfcn
        srfc3 = srfc3/srfcn
c compute local coords of contact point
c     determine if the contact point is inside the master surface triangle
c     sub areas of slave node - to triangle corners
        svx1 = x1c - xc
        svy1 = y1c - yc
        svz1 = z1c - zc
        svx2 = x2c - xc
        svy2 = y2c - yc
        svz2 = z2c - zc
        svx3 = x3c - xc
        svy3 = y3c - yc
        svz3 = z3c - zc
c area in triangle 1: 2-3-c, 2: 1-3-c, 3: 1-2-c
        a41 = ( svy2*svz3-svy3*svz2)*srfc1 + 
     *        (-svx2*svz3+svx3*svz2)*srfc2 + 
     *        ( svx2*svy3-svx3*svy2)*srfc3
        a42 = ( svy3*svz1-svy1*svz3)*srfc1 + 
     *        (-svx3*svz1+svx1*svz3)*srfc2 + 
     *        ( svx3*svy1-svx1*svy3)*srfc3
        a43 = ( svy1*svz2-svy2*svz1)*srfc1 + 
     *        (-svx1*svz2+svx2*svz1)*srfc2 + 
     *        ( svx1*svy2-svx2*svy1)*srfc3 
        a40 = a41 + a42 + a43        
c compute the minimum area and keep the "minimum" node
       if( a41 .lt. a42 )then
         if( a41 .lt. a43 )then
           a44r = a41
           imin_nd = 1
         else
           a44r = a43
           imin_nd = 3
         endif
       else
         if( a42 .lt. a43 )then
           a44r = a42
           imin_nd = 2
         else
           a44r = a43
           imin_nd = 3
         endif
       endif
c       compute the edge length
       if( imin_nd .eq. 1)then
         disx = x2c - x3c
         disy = y2c - y3c
         disz = z2c - z3c
         e_len = sqrt( disx*disx + disy*disy + disz*disz )
       else if( imin_nd .eq. 2 )then
         disx = x1c - x3c
         disy = y1c - y3c
         disz = z1c - z3c
         e_len = sqrt( disx*disx + disy*disy + disz*disz )
       else
         disx = x2c - x1c
         disy = y2c - y1c
         disz = z2c - z1c
         e_len = sqrt( disx*disx + disy*disy + disz*disz )
       endif
       tang_dist = abs(a44r)/e_len
C FIXED SHOULD THERE BE AN ABSOLUTE DISTANCE CHECK HERE I.E. COPY THE E_LEN 
C STUFF ABOVE (AND/OR THE EXISTING RELATIVE TEST) : REJ 12/05
        a40t = -ctolt1*a40
        msidf = 1
c if one of the triangle areas is negative -> no interaction
        if(( a41.lt.a40t .or. a42.lt.a40t .or. a43.lt.a40t )
     +  .and. (tang_dist .gt. ctolt2 )  )msidf = 0
c       if( a41.lt.a40t .or. a42.lt.a40t .or. a43.lt.a40t ) msidf = 0
c       if( a41.lt.a40t .or. a42.lt.a40t .or. a43.lt.a40t )
c    &   write(*,*) "<> in cnodetriangle_movsrch_aug: ",
c    &          "sub triangle test failed ", a40t, " > ", a41, a42, a43
c area coordinates
        A1 = a41/a40
        A2 = a42/a40
        A3 = a43/a40
c contact point in predicted configuration
        xc = A1*coords_triangle_pre(k_v3dx,1,i) 
     +     + A2*coords_triangle_pre(k_v3dx,2,i) 
     +     + A3*coords_triangle_pre(k_v3dx,3,i) 
        yc = A1*coords_triangle_pre(k_v3dy,1,i) 
     +     + A2*coords_triangle_pre(k_v3dy,2,i) 
     +     + A3*coords_triangle_pre(k_v3dy,3,i) 
        zc = A1*coords_triangle_pre(k_v3dz,1,i) 
     +     + A2*coords_triangle_pre(k_v3dz,2,i) 
     +     + A3*coords_triangle_pre(k_v3dz,3,i) 
c pushback direction
        pb_x = xc - coords_node_pre(k_v3dx,i)
        pb_y = yc - coords_node_pre(k_v3dy,i)
        pb_z = zc - coords_node_pre(k_v3dz,i)
        dmag = sqrt(pb_x*pb_x + pb_y*pb_y + pb_z*pb_z)

c
c  If no node motion, use the face normal for pushback direction
c
        if(dmag .eq. 0.0) then
          pb_x = srfc1
          pb_y = srfc2
          pb_z = srfc3
        else 
          pb_x = pb_x/dmag 
          pb_y = pb_y/dmag 
          pb_z = pb_z/dmag 
        endif
        pbdotn = pb_x*srfc1 + pb_y*srfc2 + pb_z*srfc3
        dmag = -sign(dmag,pbdotn) 
c
        ctrcl(ICPOINTX ,i) = xc
        ctrcl(ICPOINTY ,i) = yc
        ctrcl(ICPOINTZ ,i) = zc
        ctrcl(ICTIMC   ,i) = dtcn
        ctrcl(IPUSHX   ,i) = pb_x 
        ctrcl(IPUSHY   ,i) = pb_y 
        ctrcl(IPUSHZ   ,i) = pb_z 
        ctrcl(INORMX   ,i) = rmsnorm_aug(k_v3dx,i)
        ctrcl(INORMY   ,i) = rmsnorm_aug(k_v3dy,i)
        ctrcl(INORMZ   ,i) = rmsnorm_aug(k_v3dz,i)
        ctrcl(IPENMAG  ,i) = dmag
        ctrcl(ILOCATION,i) = iinside
        endif
        ctrcl(MSPARAM  ,i) = msidf

  100 continue

      return
      end













      subroutine cnodetriangle_cpproj_aug_noauto(npairs,
     *   coords_node_pre, coords_triangle_pre, rmsnorm_pre,
     *   coords_node_aug, coords_triangle_aug, rmsnorm_aug,
     *   ipushdr,ctrl,ctrcl,tola1,tola2,ctoln2 ) 
C
C Description: This routine computes the closest point projection for a 
C              three node triangle and a node.
C
C Memory:
C   npairs                i/-     number of pairs of faces/nodes to process
C   coords_node_pre       i/-     node coordinates
C   coords_triangle_pre   i/-     coordinates of face nodes
C   rmsnorm_pre           i/-     face normals
C   coords_node_aug       i/-     node coordinates in augmented configuration
C   coords_triangle_aug   i/-     coordinates of face nodes in augmented configuration
C   rmsnorm_aug           i/-     face normals in augmented configuration
C   ipushdr               i/o     pushback direction
C   ctrl                  -/-     scratch array
C   ctrcl                 -/o     contact information for each pair
C   tola1                 i/-
C   tola2                 i/-
C   ctoln2                i/-     search tolerance on normal gap

      implicit real*8 (a-h,o-z)
      include "search_parameters.par"
      
      dimension coords_node_pre(n_v3d,npairs)
      dimension coords_triangle_pre(n_v3d,3,npairs)
      dimension rmsnorm_pre(n_v3d,npairs)
      dimension coords_node_aug(n_v3d,npairs)
      dimension coords_triangle_aug(n_v3d,3,npairs)
      dimension rmsnorm_aug(n_v3d,npairs)
      dimension ipushdr(npairs)
   
      dimension ctrcl(ISIZCTRCL,npairs)
      dimension ctrl(3,8,npairs)
C     
      parameter ( zeron = -1.0e-5 )

      do 100 i = 1,npairs
        if( ipushdr(i) .eq. inormal )then
          ctrcl(MSPARAM,i) = 0.0
c Begin with a check for normal distance in the predicted configuration
c   predicted slave node coordinates
          xsd = coords_node_pre(k_v3dx,i)
          ysd = coords_node_pre(k_v3dy,i)
          zsd = coords_node_pre(k_v3dz,i)
c   predicted triangular facet coordinates
          x1d = coords_triangle_pre(k_v3dx,1,i)
          y1d = coords_triangle_pre(K_v3dy,1,i)
          z1d = coords_triangle_pre(k_v3dz,1,i)
C  compute normal distance from the master surface to the slave node
          vx1s = xsd - x1d
          vy1s = ysd - y1d
          vz1s = zsd - z1d
C dot the vector from point 1 on master surface to slave node
C with the outward unit normal
          projn = vx1s*rmsnorm_aug(k_v3dx,i) + 
     *            vy1s*rmsnorm_aug(k_v3dy,i) +
     *            vz1s*rmsnorm_aug(k_v3dz,i)
c if projn > ctoln2 the node is too far outside so we don't need to continue
          if( projn .le. ctoln2 )then
c augmented slave node coordinates
            xsd = coords_node_aug(k_v3dx,i)
            ysd = coords_node_aug(k_v3dy,i)
            zsd = coords_node_aug(k_v3dz,i)
c augmented triangular facet coordinates
            x1d = coords_triangle_aug(k_v3dx,1,i)
            y1d = coords_triangle_aug(K_v3dy,1,i)
            z1d = coords_triangle_aug(k_v3dz,1,i)
            x2d = coords_triangle_aug(k_v3dx,2,i)
            y2d = coords_triangle_aug(k_v3dy,2,i)
            z2d = coords_triangle_aug(k_v3dz,2,i)
            x3d = coords_triangle_aug(k_v3dx,3,i)
            y3d = coords_triangle_aug(k_v3dy,3,i)
            z3d = coords_triangle_aug(k_v3dz,3,i)
C
C  compute normal distance from the master surface to the slave node
C  note that a negative distance implies the slave node is inside the
C face
            vx1s = xsd - x1d
            vy1s = ysd - y1d
            vz1s = zsd - z1d
C dot the vector from point 1 on master surface to slave node
C with the outward unit normal
            projn = vx1s*rmsnorm_aug(k_v3dx,i) + 
     *              vy1s*rmsnorm_aug(k_v3dy,i) +
     *              vz1s*rmsnorm_aug(k_v3dz,i)
C find contact point
            xc0 = xsd - projn*rmsnorm_aug(k_v3dx,i)
            yc0 = ysd - projn*rmsnorm_aug(k_v3dy,i)
            zc0 = zsd - projn*rmsnorm_aug(k_v3dz,i)
C determine if the contact point is inside the master surface triangle
C sub areas of slave node - to triangle corners
            vx1 = x1d - xc0
            vy1 = y1d - yc0
            vz1 = z1d - zc0
            vx2 = x2d - xc0
            vy2 = y2d - yc0
            vz2 = z2d - zc0
            vx3 = x3d - xc0
            vy3 = y3d - yc0
            vz3 = z3d - zc0
c area in triangle 1: 2-3-c, 2 2: 1-3-c, 3: 1-2-c
c   (Note: Actually 2.0*Area)
            a41 = ( vy2*vz3-vy3*vz2)*rmsnorm_aug(k_v3dx,i) + 
     *            (-vx2*vz3+vx3*vz2)*rmsnorm_aug(k_v3dy,i) + 
     *            ( vx2*vy3-vx3*vy2)*rmsnorm_aug(k_v3dz,i)
            a42 = ( vy3*vz1-vy1*vz3)*rmsnorm_aug(k_v3dx,i) +
     *            (-vx3*vz1+vx1*vz3)*rmsnorm_aug(k_v3dy,i) + 
     *            ( vx3*vy1-vx1*vy3)*rmsnorm_aug(k_v3dz,i)
            a43 = ( vy1*vz2-vy2*vz1)*rmsnorm_aug(k_v3dx,i) + 
     *            (-vx1*vz2+vx2*vz1)*rmsnorm_aug(k_v3dy,i) + 
     *            ( vx1*vy2-vx2*vy1)*rmsnorm_aug(k_v3dz,i) 
            a44 = a41 + a42 + a43
            if( a41 .lt. a42 )then
              if( a41 .lt. a43 )then
                a44r = a41
                imin_nd = 1
              else
                a44r = a43
                imin_nd = 3
              endif
            else
              if( a42 .lt. a43 )then
                a44r = a42
                imin_nd = 2
              else
                a44r = a43
                imin_nd = 3
              endif
            endif

            if( a44r.ge.-tola1*a44 )then 
c re-compute contact gap and pushback based on predicted coords
c   predicted slave node coordinates
              xsd = coords_node_pre(k_v3dx,i)
              ysd = coords_node_pre(k_v3dy,i)
              zsd = coords_node_pre(k_v3dz,i)
c   predicted triangular facet coordinates
              x1d = coords_triangle_pre(k_v3dx,1,i)
              y1d = coords_triangle_pre(K_v3dy,1,i)
              z1d = coords_triangle_pre(k_v3dz,1,i)
              x2d = coords_triangle_pre(k_v3dx,2,i)
              y2d = coords_triangle_pre(k_v3dy,2,i)
              z2d = coords_triangle_pre(k_v3dz,2,i)
              x3d = coords_triangle_pre(k_v3dx,3,i)
              y3d = coords_triangle_pre(k_v3dy,3,i)
              z3d = coords_triangle_pre(k_v3dz,3,i)
c   compute global contact coordinates in the predicted configuration
c   assuming the local coords ( a41/a44 , a42/a44 , a43/a44 )
              xc0 = ( a41*x1d + a42*x2d + a43*x3d ) / a44
              yc0 = ( a41*y1d + a42*y2d + a43*y3d ) / a44
              zc0 = ( a41*z1d + a42*z2d + a43*z3d ) / a44
C   compute distance from the master surface to the slave node
              vx1s = xsd - xc0
              vy1s = ysd - yc0
              vz1s = zsd - zc0
c   re-compute gap & pushback dir
              projn = vx1s*vx1s + vy1s*vy1s + vz1s*vz1s
              if( projn .gt. 0 )then
                projn = sqrt( projn )           
                px =  vx1s/projn
                py =  vy1s/projn
                pz =  vz1s/projn
              else
                px = rmsnorm_aug(k_v3dx,i)
                py = rmsnorm_aug(k_v3dy,i)
                pz = rmsnorm_aug(k_v3dz,i)
              end if
c   follow correct sign convention
              pdn = px*rmsnorm_aug(k_v3dx,i) +
     *              py*rmsnorm_aug(k_v3dy,i) +
     *              pz*rmsnorm_aug(k_v3dz,i)
              if( pdn .lt. 0 )then
                projn = -projn
                px = -px 
                py = -py
                pz = -pz
              end if
              if( abs(projn).le.ctoln2 )then 
                ctrcl(MSPARAM,i)  = 1
                ctrcl(ICPOINTX,i) = xc0
                ctrcl(ICPOINTY,i) = yc0
                ctrcl(ICPOINTZ,i) = zc0
                ctrcl(IPENMAG,i)  = projn
                ctrcl(IPUSHX,i)   = px
                ctrcl(IPUSHY,i)   = py
                ctrcl(IPUSHZ,i)   = pz
                ctrcl(INORMX,i)   = rmsnorm_aug(k_v3dx,i)
                ctrcl(INORMY,i)   = rmsnorm_aug(k_v3dy,i)
                ctrcl(INORMZ,i)   = rmsnorm_aug(k_v3dz,i)
                ctrcl(ILOCATION,i) = iinside
              end if
            else if( abs(projn).le.ctoln2 )then 
c compute a "tangential distance"
c the area of the triangle is 0.5*base*height which can be written as
c     base = tang_dist = 2.0*area/e_len
c (Note: Because we computed 2.0*area above we will drop the 2.0 below)
c
c             compute the edge length
              if( imin_nd .eq. 1)then
                disx = x2d - x3d
                disy = y2d - y3d
                disz = z2d - z3d
                e_len = sqrt( disx*disx + disy*disy + disz*disz )
              else if( imin_nd .eq. 2 )then
                disx = x1d - x3d
                disy = y1d - y3d
                disz = z1d - z3d
                e_len = sqrt( disx*disx + disy*disy + disz*disz )
              else
                disx = x2d - x1d
                disy = y2d - y1d
                disz = z2d - z1d
                e_len = sqrt( disx*disx + disy*disy + disz*disz )
              endif
              tang_dist = abs(a44r)/e_len
              if( tang_dist .le. tola2 )then
                ipushdr(i) = iclose
              endif
            endif
          endif
        endif
  100 continue


      do 200 i = 1,npairs
C     
C push back along along minimum distance to master surface
C     
        if( ipushdr(i) .eq. iclose )then
C augmented slave node coords
          xsd = coords_node_aug(k_v3dx,i)
          ysd = coords_node_aug(k_v3dy,i)
          zsd = coords_node_aug(k_v3dz,i)
C augmented triangular facet coords
          x1d = coords_triangle_aug(k_v3dx,1,i)
          y1d = coords_triangle_aug(k_v3dy,1,i)
          z1d = coords_triangle_aug(k_v3dz,1,i)
          x2d = coords_triangle_aug(k_v3dx,2,i)
          y2d = coords_triangle_aug(k_v3dy,2,i)
          z2d = coords_triangle_aug(k_v3dz,2,i)
          x3d = coords_triangle_aug(k_v3dx,3,i)
          y3d = coords_triangle_aug(k_v3dy,3,i)
          z3d = coords_triangle_aug(k_v3dz,3,i)
C edge vectors
          vx0   = x3d - x1d
          vy0   = y3d - y1d
          vz0   = z3d - z1d
          vx00  = x2d - x1d
          vy00  = y2d - y1d
          vz00  = z2d - z1d
          vx000 = x3d - x2d
          vy000 = y3d - y2d
          vz000 = z3d - z2d
C store facet normal in the predicted config
          ctrcl(INORMX,i) = rmsnorm_aug(k_v3dx,i)
          ctrcl(INORMY,i) = rmsnorm_aug(k_v3dy,i)
          ctrcl(INORMZ,i) = rmsnorm_aug(k_v3dz,i)
          px = rmsnorm_aug(k_v3dx,i)
          py = rmsnorm_aug(k_v3dy,i)
          pz = rmsnorm_aug(k_v3dz,i)
C compute normal distance from the master surface to the slave node
C note that a negative distance implies the slave node is inside the
C face
          vx1s = xsd - x1d
          vy1s = ysd - y1d
          vz1s = zsd - z1d
C dot the vector from point 1 on master surface to slave node
C with the outward unit normal
          projn = vx1s*rmsnorm_aug(k_v3dx,i) + 
     *            vy1s*rmsnorm_aug(k_v3dy,i) + 
     *            vz1s*rmsnorm_aug(k_v3dz,i)
C find contact point
          xc0 = xsd - projn*rmsnorm_aug(k_v3dx,i)
          yc0 = ysd - projn*rmsnorm_aug(k_v3dy,i)
          zc0 = zsd - projn*rmsnorm_aug(k_v3dz,i)
C     
C projection of vector from node 1 to slave node onto the
C vector from node 1 to node 3, and the local coord, s
          proj13 = vx1s*vx0 + vy1s*vy0 + vz1s*vz0

C magnitude of vect. 1-3
          vmag13 = abs(vx0*vx0+vy0*vy0+vz0*vz0)
C local coord. (s) along 1-3
          s13 = proj13/vmag13
          s13 = max(zero,s13)
          s13 = min(one,s13)
C nearest point on line 1-3
          xc13 = x1d + s13*vx0
          yc13 = y1d + s13*vy0
          zc13 = z1d + s13*vz0
C vector from nearest point to slave node
          vcs13x = xsd-xc13
          vcs13y = ysd-yc13
          vcs13z = zsd-zc13
C distance from nearest point to slave node
          dis13 = sqrt(vcs13x*vcs13x+vcs13y*vcs13y+vcs13z*vcs13z)
C dot distance with surface normal to determine pos. or neg. dist.
          sdis = vcs13x*rmsnorm_aug(k_v3dx,i) + 
     *           vcs13y*rmsnorm_aug(k_v3dy,i) +
     *           vcs13z*rmsnorm_aug(k_v3dz,i)
          sdis13 = sign(dis13,sdis)
          if (sdis13.gt.0.0) then
            px = vcs13x/sdis13
            py = vcs13y/sdis13
            pz = vcs13z/sdis13
          endif
          ctrl(1,1,i) = s13
          ctrl(1,2,i) = xc13
          ctrl(1,3,i) = yc13
          ctrl(1,4,i) = zc13
          ctrl(1,5,i) = sdis13
          ctrl(1,6,i) = px
          ctrl(1,7,i) = py
          ctrl(1,8,i) = pz
C     
C projection of vector from node 1 to slave node onto the
C vector from node 1 to node 2, and the local coord, s
          proj12 = vx1s*vx00 + vy1s*vy00 + vz1s*vz00
C magnitude of vect. 1-2
          vmag12 = abs(vx00*vx00+vy00*vy00+vz00*vz00)
C local coord. (s) alonge 1-2
          s12 = proj12/vmag12
          s12 = max(zero,s12)
          s12 = min(one,s12)
C nearest point on line 1-2
          xc12 = x1d + s12*vx00
          yc12 = y1d + s12*vy00
          zc12 = z1d + s12*vz00
C vector from nearest point to slave node
          vcs12x = xsd-xc12
          vcs12y = ysd-yc12
          vcs12z = zsd-zc12
C distance from nearest point to slave node
          dis12 = sqrt(vcs12x*vcs12x+vcs12y*vcs12y+vcs12z*vcs12z)
C dot distance with surface normal to determine pos. or neg. dist.
          sdis = vcs12x*rmsnorm_aug(k_v3dx,i) + 
     *           vcs12y*rmsnorm_aug(k_v3dy,i) + 
     *           vcs12z*rmsnorm_aug(k_v3dz,i)
          sdis12 = sign(dis12,sdis)
          if (sdis12.gt.0.0) then
            px = vcs12x/sdis12
            py = vcs12y/sdis12
            pz = vcs12z/sdis12
          endif
          ctrl(2,1,i) = s12
          ctrl(2,2,i) = xc12
          ctrl(2,3,i) = yc12
          ctrl(2,4,i) = zc12
          ctrl(2,5,i) = sdis12
          ctrl(2,6,i) = px
          ctrl(2,7,i) = py
          ctrl(2,8,i) = pz
C
C projection of vector from node 2 to slave node onto the
C vector from node 2 to node 3, and the local coord, s
          vx1s = xsd - x2d
          vy1s = ysd - y2d
          vz1s = zsd - z2d
          proj23 = vx1s*vx000 + vy1s*vy000 + vz1s*vz000
          vmag23 = abs(vx000*vx000+vy000*vy000+vz000*vz000)
C local coord. (s) alonge 2-3
          s23 = proj23/vmag23
          s23 = max(zero,s23)
          s23 = min(one,s23)
C nearest point on line 2-3
          xc23 = x2d + s23*vx000
          yc23 = y2d + s23*vy000
          zc23 = z2d + s23*vz000
C vector from nearest point to slave node
          vcs23x = xsd-xc23
          vcs23y = ysd-yc23
          vcs23z = zsd-zc23
C distance from nearest point to slave node
          dis23 = sqrt(vcs23x*vcs23x+vcs23y*vcs23y+vcs23z*vcs23z)
C dot distance with surface normal to determine pos. or neg. dist.
          sdis = vcs23x*rmsnorm_aug(k_v3dx,i) + 
     *           vcs23y*rmsnorm_aug(k_v3dy,i) +
     *           vcs23z*rmsnorm_aug(k_v3dz,i)
          sdis23 = sign(dis23,sdis)
          if (sdis23.gt.0.0) then
            px = vcs23x/sdis23
            py = vcs23y/sdis23
            pz = vcs23z/sdis23
          endif
          ctrl(3,1,i) = s23
          ctrl(3,2,i) = xc23
          ctrl(3,3,i) = yc23
          ctrl(3,4,i) = zc23
          ctrl(3,5,i) = sdis23
          ctrl(3,6,i) = px
          ctrl(3,7,i) = py
          ctrl(3,8,i) = pz
        endif
  200 continue

      do 300 i = 1,npairs
        
        if(ipushdr(i) .eq. iclose)then
          
          ismall = 1
          clmin = abs(ctrl(1,5,i))
          if(abs(ctrl(2,5,i)) .lt. clmin)then
            ismall = 2
            clmin = abs(ctrl(2,5,i))
          endif
          if(abs(ctrl(3,5,i)) .lt. clmin)then
            ismall = 3
            clmin = abs(ctrl(3,5,i))
          endif

C predicted slave node coords
          xsd = coords_node_pre(k_v3dx,i)
          ysd = coords_node_pre(k_v3dy,i)
          zsd = coords_node_pre(k_v3dz,i)
C predicted triangular facet coords
          x1d = coords_triangle_pre(k_v3dx,1,i)
          y1d = coords_triangle_pre(k_v3dy,1,i)
          z1d = coords_triangle_pre(k_v3dz,1,i)
          x2d = coords_triangle_pre(k_v3dx,2,i)
          y2d = coords_triangle_pre(k_v3dy,2,i)
          z2d = coords_triangle_pre(k_v3dz,2,i)
          x3d = coords_triangle_pre(k_v3dx,3,i)
          y3d = coords_triangle_pre(k_v3dy,3,i)
          z3d = coords_triangle_pre(k_v3dz,3,i)
C retrieve local coord along edge
          s = ctrl(ismall,1,i)

          px = rmsnorm_pre(k_v3dx,i)
          py = rmsnorm_pre(k_v3dy,i)
          pz = rmsnorm_pre(k_v3dz,i)
          
          if( ismall .eq. 1 )then
            xc13 = x1d - s*x1d + s*x3d
            yc13 = y1d - s*y1d + s*y3d
            zc13 = z1d - s*z1d + s*z3d
            vcs13x = xsd-xc13
            vcs13y = ysd-yc13
            vcs13z = zsd-zc13
            dis13 = sqrt(vcs13x*vcs13x+vcs13y*vcs13y+vcs13z*vcs13z)
            sdis = vcs13x*rmsnorm_pre(k_v3dx,i) + 
     *             vcs13y*rmsnorm_pre(k_v3dy,i) +
     *             vcs13z*rmsnorm_pre(k_v3dz,i)
            sdis13 = sign(dis13,sdis)
            if (sdis13.gt.0.0) then
              px = vcs13x/sdis13
              py = vcs13y/sdis13
              pz = vcs13z/sdis13
            endif
            ctrcl(MSPARAM,i)   = 1
            ctrcl(ICPOINTX,i)  = xc13
            ctrcl(ICPOINTY,i)  = yc13
            ctrcl(ICPOINTZ,i)  = zc13
            ctrcl(IPENMAG,i)   = sdis13
            ctrcl(IPUSHX,i)    = px
            ctrcl(IPUSHY,i)    = py
            ctrcl(IPUSHZ,i)    = pz
            ctrcl(ILOCATION,i) = iout
          else if( ismall .eq. 2 )then
            xc12 = x1d - s*x1d + s*x2d
            yc12 = y1d - s*y1d + s*y2d
            zc12 = z1d - s*z1d + s*z2d
            vcs12x = xsd-xc12
            vcs12y = ysd-yc12
            vcs12z = zsd-zc12
            dis12 = sqrt(vcs12x*vcs12x+vcs12y*vcs12y+vcs12z*vcs12z)
            sdis = vcs12x*rmsnorm_pre(k_v3dx,i) + 
     *             vcs12y*rmsnorm_pre(k_v3dy,i) + 
     *             vcs12z*rmsnorm_pre(k_v3dz,i)
            sdis12 = sign(dis12,sdis)
            if (sdis12.gt.0.0) then
              px = vcs12x/sdis12
              py = vcs12y/sdis12
              pz = vcs12z/sdis12
            endif
            ctrcl(MSPARAM,i)   = 1
            ctrcl(ICPOINTX,i)  = xc12
            ctrcl(ICPOINTY,i)  = yc12
            ctrcl(ICPOINTZ,i)  = zc12
            ctrcl(IPENMAG,i)   = sdis12
            ctrcl(IPUSHX,i)    = px
            ctrcl(IPUSHY,i)    = py
            ctrcl(IPUSHZ,i)    = pz
            ctrcl(ILOCATION,i) = iout
          else
            xc23 = x2d - s*x2d + s*x3d
            yc23 = y2d - s*y2d + s*y3d
            zc23 = z2d - s*z2d + s*z3d
            vcs23x = xsd-xc23
            vcs23y = ysd-yc23
            vcs23z = zsd-zc23
            dis23 = sqrt(vcs23x*vcs23x+vcs23y*vcs23y+vcs23z*vcs23z)
            sdis = vcs23x*rmsnorm_pre(k_v3dx,i) + 
     *             vcs23y*rmsnorm_pre(k_v3dy,i) +
     *             vcs23z*rmsnorm_pre(k_v3dz,i)
            sdis23 = sign(dis23,sdis)
            if (sdis23.gt.0.0) then
              px = vcs23x/sdis23
              py = vcs23y/sdis23
              pz = vcs23z/sdis23
            endif
            ctrcl(MSPARAM,i)   = 1
            ctrcl(ICPOINTX,i)  = xc23
            ctrcl(ICPOINTY,i)  = yc23
            ctrcl(ICPOINTZ,i)  = zc23
            ctrcl(IPENMAG,i)   = sdis23
            ctrcl(IPUSHX,i)    = px
            ctrcl(IPUSHY,i)    = py
            ctrcl(IPUSHZ,i)    = pz
            ctrcl(ILOCATION,i) = iout
          end if
        endif
  300 continue
C     
      return
      end


      subroutine cnodetri_movsrch_aug_noauto( npairs,
     *   coords_node_cur, coords_triangle_cur, rmsnorm_cur,
     *   coords_node_pre, coords_triangle_pre, rmsnorm_pre,
     *   coords_node_aug, coords_triangle_aug, rmsnorm_aug,
     *   ctrcl, ctoln2, ctolt2 )
C     
C***********************************************************************
C     
C Description:
C This subroutine calculates the dynamic intersection of a triangular
c face and a node.  The routine will return the status of the search in
c the MSPARAM location of ctrcl.  The possibilities are
c  
c   ctrcl(MSPARAM,i) = -1    Needs a static check
c   ctrcl(MSPARAM,i) =  0    No intersection
c   ctrcl(MSPARAM,i) =  1    Contact point defined 
C
C-----------------------------------------
C
C Arguments
C 
C  INPUT:
C   npairs                i/-    number of pairs in this workset
c   coords_node_cur       i/-    current configuration of the slave nodes
c   coords_node_pre       i/-    predicted configuration of the slave nodes
c   coords_triangle_cur   i/-    current configuration of the master nodes
c   coords_triangle_pre   i/-    predicted configuration of the master nodes
c   rmsnorm_cur           i/-    face normals for current configuration
c   rmsnorm_pre           i/-    face normals for predicted configuration
c   ctoln2                i/-    search tolerance on normal gap
c   ctolt2                i/-    serach toleracne on tangential gap
c
c  INPUT/OUTPUT:
c   ctrcl                 -/o    contact information for each pair
C
C***********************************************************************
C     
c
c

      implicit real*8 (a-h,o-z)

      include "search_parameters.par"

      dimension ctrcl(isizctrcl,npairs)

      dimension coords_node_cur(n_v3d,npairs)
      dimension coords_triangle_cur(n_v3d,3,npairs)
      dimension rmsnorm_cur(n_v3d,npairs)

      dimension coords_node_pre(n_v3d,npairs)
      dimension coords_triangle_pre(n_v3d,3,npairs)
      dimension rmsnorm_pre(n_v3d,npairs)

      dimension coords_node_aug(n_v3d,npairs)
      dimension coords_triangle_aug(n_v3d,3,npairs)
      dimension rmsnorm_aug(n_v3d,npairs)

      parameter( tol1m3 = 1.D-3 )
c      parameter( ctolt1 = 1.D-6 )
      parameter( two = 2.D0 )

      ctolt1 = ctolt2

      do 100 i = 1,npairs
c
c compute normal distance from the MS to SN in current configuration
c    NOTE: a negative distance implies the slave node is penetrating the face
        vfms0x = coords_node_cur(k_v3dx,i) -
     *           coords_triangle_cur(k_v3dx,3,i)
        vfms0y = coords_node_cur(k_v3dy,i) -
     *           coords_triangle_cur(k_v3dy,3,i)
        vfms0z = coords_node_cur(k_v3dz,i) -
     *           coords_triangle_cur(k_v3dz,3,i)
        dmag0  = vfms0x *rmsnorm_cur(k_v3dx,i) + 
     *           vfms0y *rmsnorm_cur(k_v3dy,i) + 
     *           vfms0z *rmsnorm_cur(k_v3dz,i)
c
c compute normal distance from the MS to SN in augmented configuration
c    NOTE: a negative distance implies the slave node is penetrating the face
        vfmspx = coords_node_aug(k_v3dx,i) -
     *           coords_triangle_aug(k_v3dx,3,i)
        vfmspy = coords_node_aug(k_v3dy,i) -
     *           coords_triangle_aug(k_v3dy,3,i)
        vfmspz = coords_node_aug(k_v3dz,i) -
     *           coords_triangle_aug(k_v3dz,3,i)
        dmagp  = vfmspx*rmsnorm_aug(k_v3dx,i) + 
     *           vfmspy*rmsnorm_aug(k_v3dy,i) + 
     *           vfmspz*rmsnorm_aug(k_v3dz,i)
        delmag = dmagp - dmag0
        delmaga = abs(delmag)
        if( delmaga .le. 1.e-10 )then
          rel_mot = 0.0
        else
          rel_mot = delmaga/max(abs(dmagp),abs(dmag0))
        endif
c if we are outside the search tolerance in both configurations, they don't
c interact.
        if( ((dmag0 .gt. ctoln2) .and. (dmagp .gt. ctoln2)) .or.
c if we are outside in both and inside in initial they are moving apart 
     *      ((dmag0 .ge. 0.0) .and. (dmagp .gt. ctoln2)) ) then
           msidf = 0
c if dmag0 & magp are the same (i.e., they aren't moving relative to
c one another, send to CPP
        else if( rel_mot .le. 1.e-8 )then
           msidf = -1
c dynamic intersection
        else
c compute normalized time to contact
        dtcn = -dmag0/delmag
c Compute the location of the slave node at contact time
        xc = coords_node_cur(k_v3dx,i) + 
     *       (coords_node_aug(k_v3dx,i)-coords_node_cur(k_v3dx,i))*dtcn
        yc = coords_node_cur(k_v3dy,i) +  
     *       (coords_node_aug(k_v3dy,i)-coords_node_cur(k_v3dy,i))*dtcn
        zc = coords_node_cur(k_v3dz,i) +  
     *       (coords_node_aug(k_v3dz,i)-coords_node_cur(k_v3dz,i))*dtcn
c position of master surface at contact
        x1c = coords_triangle_cur(k_v3dx,1,i) + dtcn*
     *        (coords_triangle_aug(k_v3dx,1,i)-
     *         coords_triangle_cur(k_v3dx,1,i))
        y1c = coords_triangle_cur(k_v3dy,1,i) + dtcn*
     *        (coords_triangle_aug(k_v3dy,1,i)-
     *         coords_triangle_cur(k_v3dy,1,i))
        z1c = coords_triangle_cur(k_v3dz,1,i) + dtcn*
     *        (coords_triangle_aug(k_v3dz,1,i)-
     *         coords_triangle_cur(k_v3dz,1,i))
        x2c = coords_triangle_cur(k_v3dx,2,i) + dtcn*
     *        (coords_triangle_aug(k_v3dx,2,i)-
     *         coords_triangle_cur(k_v3dx,2,i))
        y2c = coords_triangle_cur(k_v3dy,2,i) + dtcn*
     *        (coords_triangle_aug(k_v3dy,2,i)-
     *         coords_triangle_cur(k_v3dy,2,i))
        z2c = coords_triangle_cur(k_v3dz,2,i) + dtcn*
     *        (coords_triangle_aug(k_v3dz,2,i)-
     *         coords_triangle_cur(k_v3dz,2,i))
        x3c = coords_triangle_cur(k_v3dx,3,i) + dtcn*
     *        (coords_triangle_aug(k_v3dx,3,i)-
     *         coords_triangle_cur(k_v3dx,3,i))
        y3c = coords_triangle_cur(k_v3dy,3,i) + dtcn*
     *        (coords_triangle_aug(k_v3dy,3,i)-
     *         coords_triangle_cur(k_v3dy,3,i))
        z3c = coords_triangle_cur(k_v3dz,3,i) + dtcn*
     *        (coords_triangle_aug(k_v3dz,3,i)-
     *         coords_triangle_cur(k_v3dz,3,i))
c compute the unit normal at the contact time
        vxp0  = x3c - x1c
        vyp0  = y3c - y1c
        vzp0  = z3c - z1c
        vxp00 = x2c - x1c
        vyp00 = y2c - y1c
        vzp00 = z2c - z1c
        srfc1 = -vyp0*vzp00+vyp00*vzp0
        srfc2 =  vxp0*vzp00-vxp00*vzp0
        srfc3 = -vxp0*vyp00+vxp00*vyp0
        srfcn = sqrt(srfc1*srfc1+srfc2*srfc2+srfc3*srfc3)
        srfc1 = srfc1/srfcn
        srfc2 = srfc2/srfcn
        srfc3 = srfc3/srfcn
c compute local coords of contact point
c     determine if the contact point is inside the master surface triangle
c     sub areas of slave node - to triangle corners
        svx1 = x1c - xc
        svy1 = y1c - yc
        svz1 = z1c - zc
        svx2 = x2c - xc
        svy2 = y2c - yc
        svz2 = z2c - zc
        svx3 = x3c - xc
        svy3 = y3c - yc
        svz3 = z3c - zc
c area in triangle 1: 2-3-c, 2: 1-3-c, 3: 1-2-c
        a41 = ( svy2*svz3-svy3*svz2)*srfc1 + 
     *        (-svx2*svz3+svx3*svz2)*srfc2 + 
     *        ( svx2*svy3-svx3*svy2)*srfc3
        a42 = ( svy3*svz1-svy1*svz3)*srfc1 + 
     *        (-svx3*svz1+svx1*svz3)*srfc2 + 
     *        ( svx3*svy1-svx1*svy3)*srfc3
        a43 = ( svy1*svz2-svy2*svz1)*srfc1 + 
     *        (-svx1*svz2+svx2*svz1)*srfc2 + 
     *        ( svx1*svy2-svx2*svy1)*srfc3 
        a40 = a41 + a42 + a43 

c++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c MWG: added 2/21 to make tangential tol absolute instead 
c of relative as per Reese's email
c compute the minimum area and keep the "minimum" node
ccc       if( a41 .lt. a42 )then
ccc         if( a41 .lt. a43 )then
ccc           a44r = a41
ccc           imin_nd = 1
ccc         else
ccc           a44r = a43
ccc           imin_nd = 3
ccc         endif
ccc       else
ccc         if( a42 .lt. a43 )then
ccc           a44r = a42
ccc           imin_nd = 2
ccc         else
ccc           a44r = a43
ccc           imin_nd = 3
ccc         endif
ccc       endif
cccc       compute the edge length
ccc       if( imin_nd .eq. 1)then
ccc         disx = x2c - x3c
ccc         disy = y2c - y3c
ccc         disz = z2c - z3c
ccc         e_len = sqrt( disx*disx + disy*disy + disz*disz )
ccc       else if( imin_nd .eq. 2 )then
ccc         disx = x1c - x3c
ccc         disy = y1c - y3c
ccc         disz = z1c - z3c
ccc         e_len = sqrt( disx*disx + disy*disy + disz*disz )
ccc       else
ccc         disx = x2c - x1c
ccc         disy = y2c - y1c
ccc         disz = z2c - z1c
ccc         e_len = sqrt( disx*disx + disy*disy + disz*disz )
ccc       endif
ccc       tang_dist = abs(a44r)/e_len
c--------------------------------------------------------
        a40t = -ctolt1*a40
        msidf = 1
c if one of the triangle areas is negative -> no interaction
c MWG: added 2/21 to make tangential tol absolute instead 
c of relative as per Reese's email
        if( a41.lt.a40t .or. a42.lt.a40t .or. a43.lt.a40t )msidf = 0
ccc        if(( a41.lt.a40t .or. a42.lt.a40t .or. a43.lt.a40t )
ccc     +     .and. (tang_dist .gt. ctolt2 )  )msidf = 0
c area coordinates
        A1 = a41/a40
        A2 = a42/a40
        A3 = a43/a40
c contact point in predicted configuration
        xc = A1*coords_triangle_pre(k_v3dx,1,i) 
     +     + A2*coords_triangle_pre(k_v3dx,2,i) 
     +     + A3*coords_triangle_pre(k_v3dx,3,i) 
        yc = A1*coords_triangle_pre(k_v3dy,1,i) 
     +     + A2*coords_triangle_pre(k_v3dy,2,i) 
     +     + A3*coords_triangle_pre(k_v3dy,3,i) 
        zc = A1*coords_triangle_pre(k_v3dz,1,i) 
     +     + A2*coords_triangle_pre(k_v3dz,2,i) 
     +     + A3*coords_triangle_pre(k_v3dz,3,i) 
c pushback direction
        pb_x = xc - coords_node_pre(k_v3dx,i)
        pb_y = yc - coords_node_pre(k_v3dy,i)
        pb_z = zc - coords_node_pre(k_v3dz,i)
        dmag = sqrt(pb_x*pb_x + pb_y*pb_y + pb_z*pb_z)

c
c  If no node motion, use the face normal for pushback direction
c
        if(dmag .eq. 0.0) then
          pb_x = srfc1
          pb_y = srfc2
          pb_z = srfc3
        else 
          pb_x = pb_x/dmag 
          pb_y = pb_y/dmag 
          pb_z = pb_z/dmag 
        endif
        pbdotn = pb_x*srfc1 + pb_y*srfc2 + pb_z*srfc3
        dmag = -sign(dmag,pbdotn) 
c
        ctrcl(ICPOINTX ,i) = xc
        ctrcl(ICPOINTY ,i) = yc
        ctrcl(ICPOINTZ ,i) = zc
        ctrcl(ICTIMC   ,i) = dtcn
        ctrcl(IPUSHX   ,i) = pb_x 
        ctrcl(IPUSHY   ,i) = pb_y 
        ctrcl(IPUSHZ   ,i) = pb_z 
        ctrcl(INORMX   ,i) = rmsnorm_aug(k_v3dx,i)
        ctrcl(INORMY   ,i) = rmsnorm_aug(k_v3dy,i)
        ctrcl(INORMZ   ,i) = rmsnorm_aug(k_v3dz,i)
        ctrcl(IPENMAG  ,i) = dmag
        ctrcl(ILOCATION,i) = iinside
        endif
        ctrcl(MSPARAM  ,i) = msidf

  100 continue

      return
      end


      subroutine cnodeline_cpproj(npairs,coords_node,coords_line,
     *                            rmsnorm,ipushdr,ctrl,ctrcl,
     *                            tola1,tola2,ctoln2 ) 
C
C Description: This routine computes the closest point projection 
C              for a two node line and a node.
C
C Memory:
C   npairs            i/-     number of pairs of faces/nodes to process
C   coords_node       i/-     node coordinates
C   coords_line       i/-     coordinates of face nodes
C   rmsnorm           i/-     face normals
C   ipushdr           i/o     pushback direction
C   ctrl              -/-     scratch array
C   ctrcl             -/o     contact information for each pair
C   tola1             i/-
C   tola2             i/-
C   ctoln2            i/-     search tolerance on normal gap

      implicit real*8 (a-h,o-z)
      include "search_parameters.par"

      dimension coords_node(n_v3d,npairs)
      dimension coords_line(n_v3d,3,npairs)
      dimension rmsnorm(n_v3d,npairs)
      dimension ipushdr(npairs)
   
      dimension ctrcl(ISIZCTRCL,npairs)
      dimension ctrl(3,8,npairs)

      do 100 i = 1,npairs
        if( ipushdr(i) .eq. inormal )then
          ctrcl(MSPARAM,i) = 0.0
c predicted slave node coordinates
          xsd = coords_node(k_v3dx,i)
          ysd = coords_node(k_v3dy,i)
          zsd = coords_node(k_v3dz,i)
c predicted line facet coordinates
          x1d = coords_line(k_v3dx,1,i)
          y1d = coords_line(K_v3dy,1,i)
          z1d = coords_line(k_v3dz,1,i)
          x2d = coords_line(k_v3dx,2,i)
          y2d = coords_line(k_v3dy,2,i)
          z2d = coords_line(k_v3dz,2,i)
C
C compute normal distance from the master surface to the slave node
C note that a negative distance implies the slave node is inside the
C face
          vx1s = xsd - x1d
          vy1s = ysd - y1d
          vz1s = zsd - z1d
C dot the vector from point 1 on master surface to slave node
C with the outward unit normal
          projn = vx1s*rmsnorm(k_v3dx,i) + 
     *            vy1s*rmsnorm(k_v3dy,i) +
     *            vz1s*rmsnorm(k_v3dz,i)
          if( abs(projn).le.ctoln2 )then 
C find contact point
            xc0 = xsd - projn*rmsnorm(k_v3dx,i)
            yc0 = ysd - projn*rmsnorm(k_v3dy,i)
            zc0 = zsd - projn*rmsnorm(k_v3dz,i)
C determine if the contact point is inside the master surface
            vx1  = xc0 - x1d
            vy1  = yc0 - y1d
            vz1  = zc0 - z1d
            vx2  = xc0 - x2d
            vy2  = yc0 - y2d
            vz2  = zc0 - z2d
            
            vx12 = x2d - x1d
            vy12 = y2d - y1d
            vz12 = z2d - z1d
            vx21 = x1d - x2d
            vy21 = y1d - y2d
            vz21 = z1d - z2d
            
            d1  = vx1*vx1 + vy1*vy1 + vz1*vz1
            d2  = vx2*vx2 + vy2*vy2 + vz2*vz2
            
            cos1c = vx1*vx12 + vy1*vy12 + vz1*vz12
            cos2c = vx2*vx21 + vy2*vy21 + vz2*vz21
          
            if( d1.le.d2 )then
              imin_nd = 1
            else
              imin_nd = 2
            endif
            
            if( cos1c.ge.0.0d0 .and. cos2c.ge.0.0d0 )then
              ctrcl(MSPARAM,i)   = 1
              ctrcl(ICPOINTX,i)  = xc0
              ctrcl(ICPOINTY,i)  = yc0
              ctrcl(ICPOINTZ,i)  = zc0
              ctrcl(IPENMAG,i)   = projn
              ctrcl(IPUSHX,i)    = rmsnorm(k_v3dx,i)
              ctrcl(IPUSHY,i)    = rmsnorm(k_v3dy,i)
              ctrcl(IPUSHZ,i)    = rmsnorm(k_v3dz,i)
              ctrcl(INORMX,i)    = rmsnorm(k_v3dx,i)
              ctrcl(INORMY,i)    = rmsnorm(k_v3dy,i)
              ctrcl(INORMZ,i)    = rmsnorm(k_v3dz,i)
              ctrcl(ILOCATION,i) = iinside
            else if( abs(projn).le.ctoln2 )then 
c compute a "tangential distance"
              if( imin_nd .eq. 1)then
                tang_dist = sqrt(d1)
              else
                tang_dist = sqrt(d2)
              endif
              if( tang_dist .le. tola2 )then
                ctrcl(MSPARAM,i)   = 1
                ctrcl(INORMX,i)    = rmsnorm(k_v3dx,i)
                ctrcl(INORMY,i)    = rmsnorm(k_v3dy,i)
                ctrcl(INORMZ,i)    = rmsnorm(k_v3dz,i)
                ctrcl(ILOCATION,i) = iout
                ipushdr(i)         = iclose
C     
C pushback along along minimum distance to master surface
C     
                if( imin_nd .eq. 1)then
                  projn = (xsd - x1d)*rmsnorm(k_v3dx,i) + 
     *                    (ysd - y1d)*rmsnorm(k_v3dy,i) +
     *                    (zsd - z1d)*rmsnorm(k_v3dz,i)
                  ctrcl(ICPOINTX,i)  = x1d
                  ctrcl(ICPOINTY,i)  = y1d
                  ctrcl(ICPOINTZ,i)  = z1d
                  ctrcl(IPENMAG,i)   = projn
                  ctrcl(IPUSHX,i)    = xsd - x1d
                  ctrcl(IPUSHY,i)    = ysd - y1d
                  ctrcl(IPUSHZ,i)    = zsd - z1d
                else
                  projn = (xsd - x2d)*rmsnorm(k_v3dx,i) + 
     *                    (ysd - y2d)*rmsnorm(k_v3dy,i) +
     *                    (zsd - z2d)*rmsnorm(k_v3dz,i)
                  ctrcl(ICPOINTX,i)  = x2d
                  ctrcl(ICPOINTY,i)  = y2d
                  ctrcl(ICPOINTZ,i)  = z2d
                  ctrcl(IPENMAG,i)   = projn
                  ctrcl(IPUSHX,i)    = xsd - x2d
                  ctrcl(IPUSHY,i)    = ysd - y2d
                  ctrcl(IPUSHZ,i)    = zsd - z2d
                endif
              endif
            endif
          endif
        endif
  100 continue
c
      return
      end
