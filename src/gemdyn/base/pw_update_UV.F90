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

!**s/r pw_update_UV - Update physical unstaggered horizonal wind
!                     components pw_uu_plus and pw_vv_plus

      subroutine pw_update_UV()
      use glb_ld
      use gmm_itf_mod
      use gmm_pw
      use gmm_vt1
      use gem_timing
      implicit none
#include <arch_specific.hf>

      integer :: istat
!     ________________________________________________________________
!
      call gemtime_start ( 5, 'PW_UPDATE', 0)

      istat = gmm_get (gmmk_pw_uu_plus_s, pw_uu_plus)
      istat = gmm_get (gmmk_pw_vv_plus_s, pw_vv_plus)
      istat = gmm_get (gmmk_ut1_s       , ut1       )
      istat = gmm_get (gmmk_vt1_s       , vt1       )

      call hwnd_stag ( pw_uu_plus,pw_vv_plus, ut1,vt1, &
                       l_minx,l_maxx,l_miny,l_maxy,l_nk,.false.)

      call gemtime_stop (5)
!     ________________________________________________________________
!
      return
      end
