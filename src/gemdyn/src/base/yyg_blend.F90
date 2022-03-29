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
!**s/r yyg_blend - Blending ut1, vt1 and zdt1 over
!                  Yin-Yang total overlap region

      subroutine yyg_blend
      use dyn_fisl_options
      use gem_timing
      use step_options
      use gmm_vt1
      use glb_ld
      use yyg_param
      implicit none
#include <arch_specific.hf>
!
!----------------------------------------------------------------------
!
      if (Schm_nblendyy                 >  0) then
      if (mod(Step_kount,Schm_nblendyy) == 0) then

         call gemtime_start ( 7, 'YYG_BLEND', 0)

         call yyg_blend_sca ( zdt1, YYG_BLEN_q2q, &
                              l_minx,l_maxx,l_miny,l_maxy,G_nk )

         call yyg_blend_uv ( ut1, vt1, l_minx,l_maxx,l_miny,l_maxy,G_nk )

         call gemtime_stop (7)

      end if
      end if
!
!----------------------------------------------------------------------
!
      return
      end subroutine yyg_blend
