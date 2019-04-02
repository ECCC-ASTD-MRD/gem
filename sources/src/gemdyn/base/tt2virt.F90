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

!**s/r tt2virt - Update physical quantities
!
      subroutine tt2virt (F_t, F_tt2tv, Minx, Maxx, Miny, Maxy, Nk)
      use glb_ld
      use gmm_itf_mod
      use gmm_pw
      use gmm_vt1
      implicit none
#include <arch_specific.hf>
!
      logical, intent(in) :: F_tt2tv
      integer, intent(in) :: Minx, Maxx, Miny, Maxy, Nk
      real, dimension(Minx:Maxx, Miny:Maxy, Nk), intent(inout) :: F_t
!
!author
!     Michel Desgagne - Dec 2009
!
      integer :: istat
      real, pointer, dimension(:,:,:)     :: tvirt,tt,hu
      real, dimension(Minx:Maxx,Miny:Maxy,Nk) :: sumpqj
      character(1) :: timelevel_S
!     ________________________________________________________________
!
!     Compute temperature from virtual temperature when tt2tv = .false.
!     Compute virtual temperature from temperature when tt2tv = .true.
!     --------------------------------------------

      timelevel_S = 'P'

      nullify (tvirt, tt, hu)

      if (.not.F_tt2tv) then
         istat = gmm_get(gmmk_tt1_s, tvirt)
      else
         istat = gmm_get(gmmk_pw_tt_plus_s, tt)
      end if

      sumpqj = 0.

      call sumhydro(sumpqj, Minx, Maxx, Miny, Maxy, Nk, timelevel_S)

      istat= gmm_get('TR/HU:P', hu)

      if (.not.F_tt2tv) then
         call mfottvh2 (F_t, tvirt, hu,sumpqj, minx, maxx, miny, maxy, &
                        l_nk,1, l_ni,1, l_nj, F_tt2tv)
      else
         call mfottvh2 (tt, F_t, hu,sumpqj, minx, maxx, miny, maxy, &
                        l_nk,1, l_ni, 1, l_nj, F_tt2tv)
      end if
!
!     ________________________________________________________________
!
      return
      end
