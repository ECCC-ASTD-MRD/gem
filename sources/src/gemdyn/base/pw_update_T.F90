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

!**s/r pw_update_T - Update physical quantities TT

      subroutine pw_update_T()
      use glb_ld
      use gmm_itf_mod
      use gmm_pw
      use gem_timing
      implicit none
#include <arch_specific.hf>

!author
!     Michel Desgagne - May 2010
!
!revision
! v4_14 - Desgagne, M.     - Initial revision

      integer :: istat
      real, pointer, dimension (:,:,:)  :: pw_tt  => null()
!     ________________________________________________________________
!
      call gemtime_start ( 5, 'PW_UPDATE', 0)
      istat = gmm_get (gmmk_pw_tt_plus_s, pw_tt    )
!
!     Compute temperature from virtual temperature
!     --------------------------------------------
!
      call tt2virt (pw_tt, .false., l_minx,l_maxx,l_miny,l_maxy,l_nk)
      call gemtime_stop (5)
!     ________________________________________________________________
!
      return
      end
