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
      use gem_options
      use gmm_pw
      use gmm_vt1
      use tr3d
      use mem_tracers
      implicit none
#include <arch_specific.hf>

      logical, intent(in) :: F_tt2tv
      integer, intent(in) :: Minx, Maxx, Miny, Maxy, Nk
      real, dimension(Minx:Maxx, Miny:Maxy, Nk), intent(inout) :: F_t

      real, pointer, dimension(:,:,:)  :: tvirt,tt
      real, dimension(Minx:Maxx,Miny:Maxy,Nk) :: sumpqj
!     ________________________________________________________________
!
!     Compute temperature from virtual temperature when tt2tv = .false.
!     Compute virtual temperature from temperature when tt2tv = .true.
!     --------------------------------------------

      nullify (tvirt, tt)

      if (.not.F_tt2tv) then
         tvirt => tt1
      else
         tt => pw_tt_plus
      end if
     
      sumpqj = 0.

      call sumhydro (sumpqj, Minx, Maxx, Miny, Maxy, Nk,Tr3d_ntr, trt1)

      if (.not.F_tt2tv) then
         call mfottvh2 (F_t, tvirt, tracers_P(Tr3d_hu)%pntr, sumpqj,&
                        minx, maxx, miny, maxy, l_nk, &
                        1-G_halox, l_ni+G_halox, 1-G_haloy, l_nj+G_haloy, F_tt2tv)
      else
         call mfottvh2 (tt, F_t, tracers_P(Tr3d_hu)%pntr, sumpqj,&
                        minx, maxx, miny, maxy, l_nk, &
                        1-G_halox, l_ni+G_halox, 1-G_haloy, l_nj+G_haloy, F_tt2tv)
      end if
!
!     ________________________________________________________________
!
      return
      end
!**s/r tt2virt - Update physical quantities
!
      subroutine tt2virt2 (F_t, F_tt2tv, Minx, Maxx, Miny, Maxy, Nk,&
                           F_i0,F_in,F_j0,F_nj)
      use glb_ld
      use gmm_pw
      use gmm_vt1
      use tr3d
      use mem_tracers
      implicit none
#include <arch_specific.hf>

      logical, intent(in) :: F_tt2tv
      integer, intent(in) :: Minx,Maxx,Miny,Maxy, Nk, F_i0,F_in,F_j0,F_nj
      real, dimension(Minx:Maxx, Miny:Maxy, Nk), intent(inout) :: F_t

      real, pointer, dimension(:,:,:)  :: tvirt,tt
      real, dimension(Minx:Maxx,Miny:Maxy,Nk) :: sumpqj
!     ________________________________________________________________
!
!     Compute temperature from virtual temperature when tt2tv = .false.
!     Compute virtual temperature from temperature when tt2tv = .true.
!     --------------------------------------------

      nullify (tvirt, tt)

      if (.not.F_tt2tv) then
         tvirt => tt1
      else
         tt => pw_tt_plus
      end if

      sumpqj = 0.

      call sumhydro (sumpqj, Minx, Maxx, Miny, Maxy, Nk,Tr3d_ntr, trt1)

      if (.not.F_tt2tv) then
         call mfottvh2 (F_t, tvirt, tracers_P(Tr3d_hu)%pntr, sumpqj,&
               minx, maxx, miny, maxy, l_nk, F_i0,F_in,F_j0,F_nj,F_tt2tv)
      else
         call mfottvh2 (tt, F_t, tracers_P(Tr3d_hu)%pntr, sumpqj,&
           minx, maxx, miny, maxy, l_nk, F_i0,F_in,F_j0,F_nj, F_tt2tv)
      end if
!
!     ________________________________________________________________
!
      return
      end subroutine tt2virt2
