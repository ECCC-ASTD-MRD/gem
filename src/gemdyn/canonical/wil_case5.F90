!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This library is free software; you can redist_8ribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! dist_8ributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END ---------------------------------

!**s/r wil_case5 - To setup Williamson Case 5: Zonal Flow over an isolated mountain (HEIGHT)

      subroutine wil_case5 (F_gz,F_mtn,F_minx,F_maxx,F_miny,F_maxy,F_nk)

      use gem_options
      use glb_ld
      use lun
      use ptopo
      use tdpack
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in)  :: F_minx,F_maxx,F_miny,F_maxy,F_nk
      real,    intent(out) :: F_gz(F_minx:F_maxx,F_miny:F_maxy,F_nk),F_mtn(F_minx:F_maxx,F_miny:F_maxy)

      !authors
      !     Abdessamad Qaddouri and Vivian Lee
      !
      !revision
      ! v5_00 - Tanguay M. - Clean Up
      !
      !object
      !==============================================================================
      !     To setup Williamson Case 5: Zonal Flow over an isolated mountain (HEIGHT)
      !     Williamson et al.,1992,JCP,102,211-224
      !==============================================================================

      integer :: i,j,k,g_i0,g_in,g_j0,g_jn,i0,in,j0,jn,zlist
      real(kind=REAL64) :: phi0_8,ubar_8,sina_8,cosa_8,phiamp_8, &
              rlon_8,rlat_8,sint_8,cost_8,                       &
              s_8(2,2),x_a_8,y_a_8,sinl_8,cosl_8,                &
              mounta_8,radius_8,dist_8,lambdc_8,thetc_8
      real :: mount(1-G_halox:G_ni+G_halox,1-G_haloy:G_nj+G_haloy), &
              picll(1-G_halox:G_ni+G_halox,1-G_haloy:G_nj+G_haloy), &
              mtloc(F_minx:F_maxx,F_miny:F_maxy), &
              gzloc(F_minx:F_maxx,F_miny:F_maxy)
!
!---------------------------------------------------------------------
!
      g_i0= 1-G_halox ; g_in= G_ni+G_halox
      g_j0= 1-G_haloy ; g_jn= G_nj+G_haloy

      i0= 1-G_halox ; in= l_ni+G_halox
      j0= 1-G_haloy ; jn= l_nj+G_haloy

      if (Lun_out>0) write(Lun_out,*) ''
      if (Lun_out>0) write(Lun_out,*) '--------------------------------------------'
      if (Lun_out>0) write(Lun_out,*) 'WILLIAMSON CASE5, Williamson et al. (1992)  '
      if (Lun_out>0) write(Lun_out,*) 'Zonal Flow over an isolated mountain        '
      if (Lun_out>0) write(Lun_out,*) '--------------------------------------------'

      phi0_8   = 5960.0d0
      ubar_8   = 20.0d0
      mounta_8 = 2000.0
      radius_8 = pi_8/9.0
      lambdc_8 = 0.5*pi_8
      thetc_8  = pi_8/6.0

      sina_8 = 0.0
      cosa_8 = 1.0

      phiamp_8 = rayt_8*omega_8*ubar_8 + ((ubar_8*ubar_8)/2.0)

      !Compute tracer for YIN
      !----------------------
      if (Ptopo_couleur==0) then

         do j=g_j0,g_jn

            rlat_8 = G_yg_8(j)

            cost_8 = cos(rlat_8)
            sint_8 = sin(rlat_8)

            do i=g_i0,g_in

               rlon_8 = G_xg_8(i)

               sinl_8 = sin(rlon_8)
               cosl_8 = cos(rlon_8)

               picll(i,j) = phi0_8-phiamp_8*sint_8**2/grav_8

               dist_8 = sqrt((rlon_8 -  lambdc_8)**2 + (rlat_8 - thetc_8)**2)

               if (dist_8 < radius_8) then
                   mount(i,j) = mounta_8*(1.0 - dist_8/radius_8)
               else
                   mount(i,j) = 0.0
               end if

            end do

         end do

      !Compute tracer for YAN
      !----------------------
      else

         do j=g_j0,g_jn

            do i=g_i0,g_in

               x_a_8 = G_xg_8(i)-acos(-1.D0)
               y_a_8 = G_yg_8(j)

               call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

               rlon_8 = rlon_8+acos(-1.D0)

               sint_8 = sin(rlat_8)
               cost_8 = cos(rlat_8)

               sinl_8 = sin(rlon_8)
               cosl_8 = cos(rlon_8)

               picll(i,j) = phi0_8-phiamp_8*sint_8**2/grav_8

               dist_8 = sqrt((rlon_8 -  lambdc_8)**2 + (rlat_8 - thetc_8)**2)

               if (dist_8 < radius_8) then
                   mount(i,j) = mounta_8*(1.0 - dist_8/radius_8)
               else
                   mount(i,j) = 0.0
               end if

            end do

         end do

      end if

      zlist = 1

      call glbdist_os (picll,gzloc,&
                       F_minx,F_maxx,F_miny,F_maxy,1,&
                       G_ni+G_halox,G_nj+G_haloy,zlist,1,1.0d0,0.d0)

      call glbdist_os (mount,mtloc,&
                       F_minx,F_maxx,F_miny,F_maxy,1,&
                       G_ni+G_halox,G_nj+G_haloy,zlist,1,1.0d0,0.d0)

      do k=1,F_nk
         F_gz(i0:in,j0:jn,k) = gzloc(i0:in,j0:jn)
      end do

      F_mtn(i0:in,j0:jn) = mtloc(i0:in,j0:jn)
!
!---------------------------------------------------------------------
!
      return
      end
