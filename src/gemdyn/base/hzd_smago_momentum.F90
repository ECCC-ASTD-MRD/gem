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

!**s/r hzd_smago_momentum   - applies horizontal diffusion based on the Smagorinsky approach

      subroutine hzd_smago_momentum()
      use gem_options
      use glb_ld
      use gmm_itf_mod
      use gmm_geof
      use gmm_vt0
      use hvdif_options
      use HORgrid_options
      use lun

      implicit none
#include <arch_specific.hf>

!
!Author:  Syed Husain
!
      integer :: istat
      logical :: smago_in_rhs_L, switch_on_wzd
!-------------------------------------------------------------------
!
      smago_in_rhs_L=.false.

      if( (hzd_smago_param <= 0.) .and. (hzd_smago_lnr(2) <=0.) .and. (.not. smago_in_rhs_L) ) return

      if (Lun_debug_L) write (Lun_out,1000)

      switch_on_wzd   = (Hzd_lnr <= 0.)

      istat = gmm_get(gmmk_ut0_s,ut0)
      istat = gmm_get(gmmk_vt0_s,vt0)
      istat = gmm_get(gmmk_zdt0_s,zdt0)
      istat = gmm_get(gmmk_tt0_s,tt0)
      istat = gmm_get(gmmk_wt0_s,wt0)

      call hzd_smago_in_split(ut0,vt0,wt0,tt0,zdt0, &
            l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.)


      if (Grd_yinyang_L) then
         call yyg_xchng_vec_uv2uv (ut0,vt0,l_minx,l_maxx,l_miny,l_maxy,G_nk)

         if (switch_on_wzd) then
            call yyg_xchng (zdt0, l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,G_nk,&
                            .false., 'CUBIC', .false.)
            call yyg_xchng (wt0 , l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,G_nk,&
                            .false., 'CUBIC', .false.)
         end if
      end if

 1000 format(3X,'SMAGO MOMENTUM DIFFUSION : (S/R HZD_SMAGO_MOMENTUM)')
!
!-------------------------------------------------------------------
!
      return
      end
