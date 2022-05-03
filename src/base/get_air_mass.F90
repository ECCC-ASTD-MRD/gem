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

!**s/r get_air_mass - Evaluate Air Mass at appropriate time

      subroutine get_air_mass (F_air_mass,F_time,F_minx,F_maxx,F_miny,F_maxy,F_nk,F_k0)

      use dyn_fisl_options
      use geomh
      use glb_ld
      use gmm_pw
      use tr3d
      use mem_tracers
      use tdpack, only: grav_8
      use ver

      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer,          intent(in) :: F_time                                       !Time 1/0
      integer,          intent(in) :: F_minx,F_maxx,F_miny,F_maxy                  !Dimension H
      integer,          intent(in) :: F_k0                                         !Scope of operator
      integer,          intent(in) :: F_nk                                         !Number of vertical levels
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(out) :: F_air_mass !Air_mass

      !object
      !==========================================
      !     Evaluate Air Mass at appropriate time
      !==========================================

      integer :: i,j,k
      real(kind=REAL64) :: dens,inv_grav
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk) :: sumq
      real, pointer, dimension(:,:,:) :: pr_m,tr
!
!---------------------------------------------------------------------
!
      !Obtain Momentum pressure levels at appropriate time
      !---------------------------------------------------
      if (F_time == 1) pr_m => pw_pm_plus
      if (F_time == 0) pr_m => pw_pm_moins

      !Evaluate water tracers at appropriate time if dry mixing ratio
      !--------------------------------------------------------------
      sumq = 0.

      if (Schm_dry_mixing_ratio_L) then

         if (F_time == 1) then
            call sumhydro (sumq,F_minx,F_maxx,F_miny,F_maxy,F_nk,Tr3d_ntr, trt1)
            tr=>tracers_P(Tr3d_hu)%pntr
         endif
         if (F_time == 0) then
            call sumhydro (sumq,F_minx,F_maxx,F_miny,F_maxy,F_nk,Tr3d_ntr, trt0)
            tr=>tracers_M(Tr3d_hu)%pntr
         endif

         do k=F_k0,l_nk
            sumq(1:l_ni,1:l_nj,k) = sumq(1:l_ni,1:l_nj,k) + tr(1:l_ni,1:l_nj,k)
         end do

      end if

      !Evaluate Air Mass
      !-----------------
      inv_grav= 1.d0 / grav_8
      F_air_mass(:,:,1:F_k0-1)= 0.
      do k=F_k0,F_nk
         F_air_mass(:,F_miny:0     ,k)= 0.
         F_air_mass(:,l_nj+1:F_maxy,k)= 0.
         do j=1,l_nj
         F_air_mass(F_minx:0     ,j,k)= 0.
         F_air_mass(l_ni+1:F_maxx,j,k)= 0.
         do i=1,l_ni
            dens= + (pr_m(i,j,k+1) - pr_m(i,j,k)) * (1.-sumq(i,j,k)) * Ver_idz_8%t(k) * inv_grav
            F_air_mass(i,j,k) = dens * geomh_area_8(i,j) * Ver_dz_8%t(k)
         end do
         end do
      end do
!
!---------------------------------------------------------------------
!
      return
      end subroutine get_air_mass
