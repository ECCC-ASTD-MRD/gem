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

!**s/r get_density - Evaluate Fluid's density and mass

      subroutine get_density (F_density,F_mass,F_time,F_minx,F_maxx,F_miny,F_maxy,F_nk,F_k0)

      use dyn_fisl_options
      use dynkernel_options
      use geomh
      use glb_ld
      use gmm_itf_mod
      use gmm_vt1
      use gmm_vt0
      use metric
      use tdpack, only: rgasd_8, grav_8
      use ver

      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer,           intent(in) :: F_time                                     !Time 0 or Time 1
      integer,           intent(in) :: F_minx,F_maxx,F_miny,F_maxy                !Dimension H
      integer,           intent(in) :: F_k0                                       !Scope of operator
      integer,           intent(in) :: F_nk                                       !Number of vertical levels
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(out) :: F_density !Fluid's density
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(out) :: F_mass    !Fluid's mass

      !object
      !======================================
      !     Evaluate Fluid's density and mass
      !======================================

      integer :: i,j,k,istat
      logical :: GEM_P_L
      real, dimension(F_minx:F_maxx,F_miny:F_maxy)       :: pr_p0
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk+1):: pr_m,pr_t
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk)  :: sumq
      real, pointer, dimension(:,:,:) :: tr,w_qt
      real, pointer, dimension(:,:)   :: w2d
      character(len=1) :: timelevel_S
!
!---------------------------------------------------------------------
!
      GEM_P_L = trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_P'

      !Obtain momentum pressure levels at appropriate time
      !---------------------------------------------------
      if (GEM_P_L) then

         if (F_time == 0) istat = gmm_get(gmmk_st0_s,w2d)
         if (F_time == 1) istat = gmm_get(gmmk_st1_s,w2d)

         call calc_pressure ( pr_m, pr_t, pr_p0, w2d, F_minx, F_maxx, F_miny, F_maxy, F_nk )

         pr_m(1:l_ni,1:l_nj,F_nk+1) = pr_p0(1:l_ni,1:l_nj)

      else

         if (F_time == 0) istat = gmm_get(gmmk_qt0_s,w_qt)
         if (F_time == 1) istat = gmm_get(gmmk_qt1_s,w_qt)

         do k=1,F_nk+1
            pr_m(1:l_ni,1:l_nj,k) = exp( (w_qt(1:l_ni,1:l_nj,k)/(rgasd_8*Ver_Tstar_8%m(k))+lg_pstar(1:l_ni,1:l_nj,k)) )
         end do

      end if

      !Evaluate water tracers if dry mixing ratio
      !------------------------------------------
      sumq = 0.

      if (Schm_dry_mixing_ratio_L) then

         if (F_time == 1) timelevel_S = 'P'
         if (F_time == 0) timelevel_S = 'M'

         call sumhydro (sumq,F_minx,F_maxx,F_miny,F_maxy,F_nk,timelevel_S)

         istat = gmm_get('TR/HU:'//timelevel_S,tr)

         do k=1,l_nk
            sumq(1:l_ni,1:l_nj,k) = sumq(1:l_ni,1:l_nj,k) + tr(1:l_ni,1:l_nj,k)
         end do

      end if

      !Evaluate Fluid's density and mass
      !---------------------------------
      do k=F_k0,F_nk
         do j=1,l_nj
         do i=1,l_ni
            F_density(i,j,k) = + (pr_m(i,j,k+1) - pr_m(i,j,k)) * (1.-sumq(i,j,k)) * Ver_idz_8%t(k) / grav_8
         end do
         end do
      end do

      do k=F_k0,F_nk
         do j=1,l_nj
         do i=1,l_ni
            F_mass(i,j,k) = F_density(i,j,k) * geomh_area_8(i,j) * Ver_dz_8%t(k)
         end do
         end do
      end do
!
!---------------------------------------------------------------------
!
      return
      end subroutine get_density
