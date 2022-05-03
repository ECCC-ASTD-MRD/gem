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

!**s/r hzd_smago_main   - applies horizontal diffusion based on the Smagorinsky approach

      subroutine hzd_smago_main()
      use gmm_vt1
      use dyn_fisl_options
      use glb_ld
      use lun
      use gmm_smag
      use hzd_mod
      use tr3d
      use step_options
      use hvdif_options
      use HORgrid_options
      use gmm_geof
      use mem_tracers
      implicit none
#include <arch_specific.hf>
!
!Author:  Syed Husain
!
      logical :: switch_on_THETA, switch_on_hu, switch_on_wzd
      logical yyblend, smago_in_rhs_L
!
!-------------------------------------------------------------------
!
      smago_in_rhs_L=.false.
      if( (hzd_smago_param <= 0.) .and. (hzd_smago_lnr(2) <=0.) .and. (.not. smago_in_rhs_L)) return

      switch_on_THETA = (Hzd_smago_prandtl > 0.) .and. (Hzd_lnr_theta <= 0.)
      switch_on_hu    = (Hzd_smago_prandtl_hu > 0.) .and. (Hzd_lnr_tr <= 0.)
      switch_on_wzd   = (Hzd_lnr <= 0.)

      if (Lun_debug_L) write (Lun_out,1000)

      call pw_update_GW()

      call hzd_smago_in_split(ut1,vt1,wt1,tt1,zdt1, &
               l_minx,l_maxx,l_miny,l_maxy,G_nk,.false.)

      if (Grd_yinyang_L) then
         call yyg_xchng_vec_uv2uv (ut1,vt1,l_minx,l_maxx,l_miny,l_maxy,G_nk)
         if (switch_on_THETA) then
            call yyg_xchng (tt1 , l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,G_nk,&
                            .false., 'CUBIC', .false.)
         end if

         if (switch_on_hu) then
            call yyg_xchng (tracers_P(Tr3d_hu)%pntr,l_minx,l_maxx,l_miny,l_maxy,&
                            l_ni,l_nj,G_nk,.false., 'CUBIC', .false.)
         end if

         if (switch_on_wzd) then
            call yyg_xchng (zdt1, l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,G_nk,&
                            .false., 'CUBIC', .false.)
            call yyg_xchng (wt1 , l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,G_nk,&
                            .false., 'CUBIC', .false.)
         end if

         yyblend= (Schm_nblendyy > 0)
         if (yyblend) then
             call yyg_blend (mod(Step_kount,Schm_nblendyy) == 0)
         end if
      end if

 1000 format(3X,'MAIN SMAGO DIFFUSION : (S/R HZD_SMAGO_MAIN)')
!
!-------------------------------------------------------------------
!
      return
      end
