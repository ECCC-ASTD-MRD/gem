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

!**s/r yyg_xchng - Interpolate and exchange scalar fields (q2q)

      subroutine yyg_xchng ( F_src, Minx,Maxx,Miny,Maxy, F_ni, F_nj, Nk,&
                             mono_L, F_interpo_S , do_xch)
      use yyg_param
      implicit none
#include <arch_specific.hf>

      character(len=*), intent(in)  :: F_interpo_S
      logical, intent(in) :: mono_L, do_xch
      integer, intent(in) :: Minx,Maxx,Miny,Maxy, Nk, F_ni,F_nj
      real, intent(inout) :: F_src (Minx:Maxx,Miny:Maxy,Nk)
!
!----------------------------------------------------------------------
!
!     For the vast majority of scalar fields exchanges
      if (trim(F_interpo_S) == 'CUBIC') then
         call yyg_xchng_sca_q2q ( F_src, YYG_PILT_q2q              ,&
                                  Minx,Maxx,Miny,Maxy, F_ni,F_nj,Nk,&
                                  'CUBIC', mono_L, do_xch )

!     For horizontal filtering in the physics interface
      elseif (trim(F_interpo_S) == 'PHYSI') then
         call yyg_xchng_sca_q2q ( F_src, YYG_HALO_q2q              ,&
                                  Minx,Maxx,Miny,Maxy, F_ni,F_nj,Nk,&
                                  'CUBIC', mono_L, do_xch )

!     For geophysical fields in itf_phy_init
      else
         call yyg_xchng_sca_q2q ( F_src, YYG_NEAR_q2q              ,&
                                  Minx,Maxx,Miny,Maxy, F_ni,F_nj,Nk,&
                                  F_interpo_S, .false., do_xch )
      end if
!
!----------------------------------------------------------------------
!
      return
      end subroutine yyg_xchng

