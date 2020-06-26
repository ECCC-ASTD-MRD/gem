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

!**s/r set_post_tr - Resetting done at each timestep before calling adz_post_tr

      subroutine set_post_tr

      use adz_mem
      use adz_options
      use gem_timing
      use gmm_pw
      use gmm_tracers

      implicit none

#include <arch_specific.hf>

      !object
      !===============================================================
      !     Resetting done at each timestep before calling adz_post_tr
      !===============================================================

      integer :: i,j,k
      real, pointer, dimension (:,:,:) :: pt_0
      real, pointer, dimension (:,:)   :: p0_0
!
!     ---------------------------------------------------------------
!
      call gemtime_start (36, 'C_SETPOST', 33)

      !SET_POST_TR has just been launched
      !----------------------------------
      Adz_set_post_tr  = 1

      if (Adz_verbose>0) call stat_mass_tracers (1,"BEFORE ADVECTION")

      call gemtime_start (52, 'AIR_MASS', 36)

      !Reset Air Mass at TIME P and TIME M
      !-----------------------------------
      call get_air_mass (airm1,1,l_minx,l_maxx,l_miny,l_maxy,l_nk, Adz_k0)
      call get_air_mass (airm0,0,l_minx,l_maxx,l_miny,l_maxy,l_nk, Adz_k0)

      call gemtime_stop  (52)

      !Reset pr_t(k)/pr_s at TIME M (used in Bermejo-Conde)
      !----------------------------------------------------
      pt_0 => pw_pt_moins
      p0_0 => pw_p0_moins

      do k=1,l_nk
         do j=1,l_nj
            do i=1,l_ni
               pkps(i,j,k) = pt_0(i,j,k) / p0_0(i,j)
            end do
         end do
      end do

      call gemtime_stop (36)
!
!     ---------------------------------------------------------------
!
      return
      end
