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
      use gmm_pw
      use outp
      use mem_tracers
      use tr3d
      use omp_timing
      implicit none
#include <arch_specific.hf>

!     ________________________________________________________________
!
      call gtmg_start ( 5, 'PW_UPDATE', 0)
!
!     Compute temperature from virtual temperature
!     --------------------------------------------
!
      call tt2virt (pw_tt_plus, .false.,l_minx,l_maxx,l_miny,l_maxy,l_nk)
      tdiag(:,:) = pw_tt_plus(:,:,G_nk)
      qdiag(:,:) = tracers_P(Tr3d_hu)%pntr(:,:,G_nk)

      call gtmg_stop (5)
!     ________________________________________________________________
!
      return
      end
