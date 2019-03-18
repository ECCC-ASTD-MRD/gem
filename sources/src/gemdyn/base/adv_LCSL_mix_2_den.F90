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

!**s/r adv_LCSL_mix_2_den : F_action=1 => Mixing ratio to Density (on Advection grid)
!                           F_action=2 => Density to Mixing ratio

      subroutine adv_LCSL_mix_2_den (F_den_a,F_mix,F_aminx,F_amaxx,F_aminy,F_amaxy, &
                                     F_minx,F_maxx,F_miny,F_maxy,F_nk,F_k0,F_action)

      use glb_ld
      use adv_grid

      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer :: F_aminx,F_amaxx,F_aminy,F_amaxy !Dimension H Advection
      integer :: F_minx,F_maxx,F_miny,F_maxy     !Dimension H
      integer :: F_nk,F_k0                       !Scope of operator
      integer :: F_action                        !Action to be completed
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk) :: F_mix
      real, dimension(F_aminx:F_amaxx,F_aminy:F_amaxy,F_nk) :: F_den_a

      !object
      !==================================================================
      !     F_action=1 => Mixing ratio to Density  (on Advection grid)
      !                   I: F_mix   = Mixing ratio
      !                   O: F_den_a = Density (on Advection grid)
      !------------------------------------------------------------------
      !     F_action=2 => Density to Mixing ratio
      !                   I: F_mix = Density
      !                   O: F_mix = Mixing ratio
      !==================================================================

      !---------------------------------------------------------------------------

      integer i,j,k,time
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk) :: air_density,bidon,den_m

      !---------------------------------------------------------------------------

      if (F_action == 1) time = 1
      if (F_action == 2) time = 0

      call get_density (air_density,bidon,time,F_minx,F_maxx,F_miny,F_maxy,F_nk,F_k0)

      !Mixing ratio to Density (on Advection grid)
      !-------------------------------------------
      if (F_action == 1) then

!$omp parallel private(i,j,k) shared(den_m,air_density)
!$omp do
         do k=1,F_nk
         do j=1,l_nj
         do i=1,l_ni
            den_m(i,j,k) = F_mix(i,j,k) * air_density(i,j,k)
         end do
         end do
         end do
!$omp end do
!$omp end parallel

         !Transfer to advection grid
         !--------------------------
         call rpn_comm_xch_halox( den_m, F_minx, F_maxx,F_miny, F_maxy, &
               l_ni, l_nj, F_nk, adv_halox, adv_haloy, G_periodx, G_periody, &
               F_den_a, F_aminx, F_amaxx, F_aminy, F_amaxy, l_ni, 0)

      end if

      !Density to Mixing ratio
      !-----------------------
      if (F_action == 2) then

!$omp parallel private(i,j,k) shared(air_density)
!$omp do
         do k=1,F_nk
         do j=1,l_nj
         do i=1,l_ni
            F_mix(i,j,k) = F_mix(i,j,k) / air_density(i,j,k)
         end do
         end do
         end do
!$omp end do
!$omp end parallel

      end if

      return

      end subroutine adv_LCSL_mix_2_den
