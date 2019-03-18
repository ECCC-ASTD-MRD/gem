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

      subroutine adz_cubic ( F_dest,  F_src, F_xyz, F_ni,F_nj,F_nk,&
                             F_minx,F_maxx,F_miny,F_maxy          ,&
                             F_i0, F_in, F_j0, F_jn               ,&
                             F_k0, F_lev_S, F_mono_L )
      use adz_mem
      implicit none

      character(len=1), intent(in) :: F_lev_S
      logical, intent(in) :: F_mono_L
      integer, intent(in) :: F_ni,F_nj,F_nk
      integer, intent(in) :: F_minx,F_maxx,F_miny, F_maxy
      integer, intent(in) :: F_i0, F_j0, F_in, F_jn, F_k0
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(in) :: F_src
      real, dimension(F_minx:F_maxx,F_miny:F_maxy,F_nk), intent(out) :: F_dest
      real, dimension(3,F_ni,F_nj,F_nk), intent(in) :: F_xyz

      integer NK
      real, dimension(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,F_nk) :: extended
      real, dimension(F_ni,F_nj,F_nk) :: wrkc
!
!---------------------------------------------------------------------
!
      call rpn_comm_xch_halox( F_src, l_minx, l_maxx, l_miny, l_maxy,&
            l_ni, l_nj, l_nk, Adz_halox, Adz_haloy, .false., .false.,&
            extended, Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy, l_ni, 0)

      NK = l_nk-F_k0+1
      if (F_mono_L) then
         call adz_tricub_lag3d_mono (wrkc, extended, F_xyz(1,1,1,F_k0),&
                                     Adz_2dnh*NK, F_lev_S)
      else
         call adz_tricub_lag3d      (wrkc, extended, F_xyz(1,1,1,F_k0),&
                                     Adz_2dnh*NK, F_lev_S)
      end if

      F_dest(F_i0:F_in,F_j0:F_jn,F_k0:l_nk) = &
        wrkc(F_i0:F_in,F_j0:F_jn,1:NK)
!
!---------------------------------------------------------------------
!
      return
      end subroutine adz_cubic
