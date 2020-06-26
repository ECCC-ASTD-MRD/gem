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

!**s/r hzd_momentum - applies horizontal diffusion on zdt and possibly on u and v
!                     to smooth momentum components for trajectory calculations after
!                     a Crank-Nicholson step avoiding pole problems

      subroutine hzd_momentum
      use gmm_vt0
      use hzd_ctrl
      use HORgrid_options
      use hvdif_options
      use dyn_fisl_options
      use glb_ld
      use lun
      implicit none
#include <arch_specific.hf>

      logical, save :: switch_on_hzd= .true.
!
!-------------------------------------------------------------------
!
      if (Schm_hzdadw_L .and. switch_on_hzd) then
         if (Lun_debug_L) write (Lun_out,1000)

         call hzd_ctrl4 ( ut0, vt0, l_minx,l_maxx,l_miny,l_maxy,G_nk)
         call hzd_ctrl4 (zdt0, 'S', l_minx,l_maxx,l_miny,l_maxy,G_nk)

         if (Grd_yinyang_L) then
            call yyg_xchng (zdt0,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                            G_nk, .false., 'CUBIC', .false.)
            call yyg_xchng_vec_uv2uv (ut0,vt0,l_minx,l_maxx,&
                                      l_miny,l_maxy,G_nk)
         end if
      end if

      call hzd_smago_momentum()

1000  format(5X,'DIFFUSION ON U,V,ZD: (S/R HZD_MOMENTUM)')
!
!-------------------------------------------------------------------
!
      return
      end subroutine hzd_momentum
