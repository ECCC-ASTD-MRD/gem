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

!**s/r nest_bcs

      subroutine nest_bcs ( F_rhsu, F_rhsv, Minx, Maxx, Miny, Maxy, Nk )
      use gmm_nest
      use lam_options
      use glb_ld
      use cstv
      implicit none
#include <arch_specific.hf>
      integer, intent(in) :: Minx, Maxx, Miny, Maxy, Nk
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(out) :: &
                                               F_rhsu, F_rhsv
      integer i,j,k
!
!----------------------------------------------------------------------
!
      if (.not. Lam_ctebcs_L) call nest_intt

      call nest_bcs_t0 ()

!**************************************
! Apply HORIZONTAL BOUNDARY CONDITIONS
!**************************************

      if (l_west) then
         do k=1,l_nk
         do j= 1+pil_s, l_nj-pil_n
               F_rhsu (pil_w,j,k) = Cstv_invT_m_8 * nest_u(pil_w,j,k)
         end do
         end do
      end if

      if (l_east) then
         do k=1,l_nk
         do j= 1+pil_s, l_nj-pil_n
            F_rhsu (l_ni-pil_e,j,k) = Cstv_invT_m_8 * nest_u(l_ni-pil_e,j,k)
         end do
         end do
      end if

      if (l_south) then
         do k=1,l_nk
         do i= 1+pil_w, l_ni-pil_e
            F_rhsv (i,pil_s,k) = Cstv_invT_m_8 * nest_v(i,pil_s,k)
         end do
         end do
      end if

      if (l_north) then
         do k=1,l_nk
         do i= 1+pil_w, l_ni-pil_e
            F_rhsv (i,l_nj-pil_n,k) = Cstv_invT_m_8 * nest_v(i,l_nj-pil_n,k)
         end do
         end do
      end if
!
!----------------------------------------------------------------------
!
      return
      end
