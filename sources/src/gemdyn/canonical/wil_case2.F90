!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END ---------------------------------

!**s/r wil_case2 - To setup Williamson Case 2: Steady state nonlinear geostrophic flow (GEOPOTENTIAL)

      subroutine wil_case2 (F_gz,F_minx,F_maxx,F_miny,F_maxy,F_nk)

      use wil_options
      use gem_options
      use tdpack

      use glb_ld
      use lun
      use ptopo
      implicit none

      integer F_minx,F_maxx,F_miny,F_maxy,F_nk
      real    F_gz(F_minx:F_maxx,F_miny:F_maxy,F_nk)

      !authors
      !     Abdessamad Qaddouri and Vivian Lee
      !
      !revision
      ! v5_00 - Tanguay M. - Clean Up
      !
      !object
      !=======================================================================================
      !     To setup Williamson Case 2: Steady state nonlinear geostrophic flow (GEOPOTENTIAL)
      !     Williamson et al.,1992,JCP,102,211-224
      !=======================================================================================


      !---------------------------------------------------------------

      integer i,j,k
      real*8  phi0_8,ubar_8,sina_8,cosa_8,phiamp_8, &
              rlon_8,rlat_8,sint_8,cost_8,          &
              s_8(2,2),x_a_8,y_a_8,sinl_8,cosl_8
      real    gzloc(F_minx:F_maxx,F_miny:F_maxy),picll(G_ni,G_nj)

      !---------------------------------------------------------------

      if (Lun_out>0) write(Lun_out,*) ''
      if (Lun_out>0) write(Lun_out,*) '--------------------------------------------'
      if (Lun_out>0) write(Lun_out,*) 'WILLIAMSON CASE2, Williamson et al. (1992)  '
      if (Lun_out>0) write(Lun_out,*) 'Steady state nonlinear geostrophic flow     '
      if (Lun_out>0) write(Lun_out,*) '--------------------------------------------'

      phi0_8   = 29400.0d0
      ubar_8   = (2.0*pi_8*rayt_8)/(12.0*24.0*3600.0)
      phiamp_8 = rayt_8*omega_8*ubar_8 + ((ubar_8*ubar_8)/2.0)

      sina_8 = sin(Williamson_alpha)
      cosa_8 = cos(Williamson_alpha)

      !Compute tracer for YIN
      !----------------------
      if (Ptopo_couleur==0) then

         do j=1,G_nj

            rlat_8 = G_yg_8(j)

            cost_8 = cos(rlat_8)
            sint_8 = sin(rlat_8)

            do i=1,G_ni

               rlon_8 = G_xg_8(i)

               sinl_8 = sin(rlon_8)
               cosl_8 = cos(rlon_8)

               picll(i,j) = (phi0_8-phiamp_8*(- cosl_8*cost_8*sina_8 +   &
                            sint_8*cosa_8)*(- cosl_8*cost_8*sina_8 +sint_8*cosa_8))/grav_8

            enddo

         enddo

      !Compute tracer for YAN
      !----------------------
      else

         do j=1,G_nj

            do i=1,G_ni

               x_a_8 = G_xg_8(i)-acos(-1.D0)
               y_a_8 = G_yg_8(j)

               call smat(s_8,rlon_8,rlat_8,x_a_8,y_a_8)

               rlon_8 = rlon_8+acos(-1.D0)

               sint_8 = sin(rlat_8)
               cost_8 = cos(rlat_8)

               sinl_8 = sin(rlon_8)
               cosl_8 = cos(rlon_8)

               picll(i,j) = (phi0_8-phiamp_8*(- cosl_8*cost_8*sina_8 + &
                            sint_8*cosa_8)*(- cosl_8*cost_8*sina_8 +sint_8*cosa_8)) /grav_8

            enddo

         enddo

      endif

      call glbdist (picll,G_ni,G_nj,gzloc,F_minx,l_maxx,l_miny,l_maxy,1,G_halox,G_haloy)

      do k=1,F_nk
         F_gz(1:l_ni,1:l_nj,k) = gzloc(1:l_ni,1:l_nj)
      enddo

      !---------------------------------------------------------------

      return
      end
