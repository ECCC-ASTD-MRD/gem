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

      subroutine hzd_momentum_hlt ()
      use gmm_vt0
      use hzd_exp_hlt
      use HORgrid_options
      use hvdif_options
      use dyn_fisl_options
      use mem_tstp
      use glb_ld
      implicit none
!
!-------------------------------------------------------------------
!
      if (Schm_hzdadw_L) then

         call hzd_exp_deln ( ut0, Hzd_pwr, Hzd_lnR, WS1, l_minx,l_maxx,l_miny,l_maxy,G_nk)
         call hzd_exp_deln ( vt0, Hzd_pwr, Hzd_lnR, WS1, l_minx,l_maxx,l_miny,l_maxy,G_nk)
         call hzd_exp_deln (zdt0, Hzd_pwr, Hzd_lnR, WS1, l_minx,l_maxx,l_miny,l_maxy,G_nk)

         if (Grd_yinyang_L) then
!$omp single
            call yyg_xchng (zdt0,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                            G_nk, .false., 'CUBIC', .false.)
            call yyg_xchng_vec_uv2uv (ut0,vt0,l_minx,l_maxx,&
                                      l_miny,l_maxy,G_nk)
!$omp end single
         end if
      end if

!$omp single
      call hzd_smago_momentum()
!$omp end single
!
!-------------------------------------------------------------------
!
      return
      end subroutine hzd_momentum_hlt
