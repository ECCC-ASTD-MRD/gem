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
!**s/r yyg_scalgeo - Interpolate and exchange non-halo geophysical data
!                    using local PEs

      subroutine yyg_scalgeo ( F_src, F_minx,F_maxx,F_miny,F_maxy, &
                               NK, Interp_S )
      use glb_ld
      implicit none
#include <arch_specific.hf>

      character(len=*) interp_S
      integer F_minx,F_maxx,F_miny,F_maxy,NK
      real F_src(F_minx:F_maxx,F_miny:F_maxy,NK)

!author
!           Abdessamad Qaddouri/V.Lee - October 2009


      real  tab_dstf(l_minx:l_maxx,l_miny:l_maxy,NK)
      character(len=32) :: UPinterp_S
!
!----------------------------------------------------------------------
!
      call low2up  (interp_S,UPinterp_S)

!     Copy original fields to a source and destination tables with
!     zeroed halo

      tab_dstf = 0.
      tab_dstf(1:l_ni,1:l_nj,1:Nk)= F_src(1:l_ni,1:l_nj,1:Nk)

      call yyg_int_xch_scal (tab_dstf, Nk, .false., trim(UPinterp_S), .false.)

      F_src(1:l_ni,1:l_nj,1:Nk)= tab_dstf(1:l_ni,1:l_nj,1:Nk)
!
!----------------------------------------------------------------------
!
      return
      end

