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
      use dynkernel_options
      use gem_options
      use gmm_itf_mod
      use gmm_tracers
      use gmm_vt0
      use metric
      use rstr
      use tdpack, only : rgasd_8
      use tr3d
      use ver
      implicit none

#include <arch_specific.hf>

      !object
      !===============================================================
      !     Resetting done at each timestep before calling adz_post_tr
      !===============================================================

      integer :: err,i,j,k
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,l_nk)   :: bidon
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,l_nk+1) :: pr_m,pr_t,log_pr_m,log_pr_t
      real, dimension(l_minx:l_maxx,l_miny:l_maxy)        :: pr_p0
!
!     ---------------------------------------------------------------
!
      if (Adz_verbose>0) call stat_mass_tracers (1,"BEFORE ADVECTION")

      !Reset Air Mass at TIME 1 and TIME 0
      !-----------------------------------
      err = gmm_get(gmmk_airm1_s,airm1)
      err = gmm_get(gmmk_airm0_s,airm0)

      airm1 = 0. ; airm0 = 0. !J'en ai besoin... Comment m'y prendre autrement?

      call get_density ( bidon, airm1, 1,&
                         l_minx,l_maxx,l_miny,l_maxy,l_nk, Adz_k0 )
      call get_density ( bidon, airm0, 0,&
                         l_minx,l_maxx,l_miny,l_maxy,l_nk, Adz_k0 )

      !Fill Halo of Air Mass at TIME 0 (required in ILMC)
      !--------------------------------------------------
      call rpn_comm_xch_halo (airm0,l_minx,l_maxx,l_miny,l_maxy,&
            l_ni,l_nj,l_nk, G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

      !Reset pr_t(k)/pr_s at TIME 0 (used in Bermejo-Conde)
      !----------------------------------------------------
      if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_P') then

         err = gmm_get(gmmk_st0_s, st0)

         call calc_pressure (pr_m,pr_t,pr_p0,st0,l_minx,l_maxx,l_miny,l_maxy,l_nk)

      else if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H') then

         err = gmm_get(gmmk_qt0_s, qt0)

         do k=1,l_nk+1
            log_pr_m(1:l_ni,1:l_nj,k) = (qt0(1:l_ni,1:l_nj,k)/(rgasd_8*Ver_Tstar_8%m(k))+lg_pstar(1:l_ni,1:l_nj,k))
         end do

         pr_p0(1:l_ni,1:l_nj) = exp(log_pr_m(1:l_ni,1:l_nj,l_nk+1))

         do k=1,l_nk
            log_pr_t(1:l_ni,1:l_nj,k) = 0.5*(log_pr_m(1:l_ni,1:l_nj,k+1)+log_pr_m(1:l_ni,1:l_nj,k))
         end do

         do k=1,l_nk
            pr_t(1:l_ni,1:l_nj,k) = exp(log_pr_t(1:l_ni,1:l_nj,k))
         end do

      end if

      err = gmm_get(gmmk_pkps_s,pkps)

      do k=1,l_nk
         do j=1,l_nj
            do i=1,l_ni
               pkps(i,j,k) = pr_t(i,j,k) / pr_p0(i,j)
            end do
         end do
      end do
!
!     ---------------------------------------------------------------
!
      return
      end
