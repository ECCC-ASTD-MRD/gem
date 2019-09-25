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

!**s/r sumhydrom - Sum over Hydrometeors (GMM)
!
      subroutine sumhydro (F_qh,minx,maxx,miny,maxy,nk,F_timelevel_S)
      use dyn_fisl_options
      use glb_ld
      use tr3d
      use gmm_itf_mod
      implicit none
#include <arch_specific.hf>

      integer minx,maxx,miny,maxy,nk
      character(1) F_timelevel_S
      real F_qh(minx:maxx,miny:maxy,nk)

      integer i, j, k, n,istat
      real, pointer, dimension(:,:,:)     :: tr
!     ________________________________________________________________
!
      F_qh = 0.

      if (.not.Schm_wload_L) return

!     Sum over Hydrometeors
      do n = 1, Tr3d_ntr
         if (Tr3d_wload (n)) then
            nullify (tr)
            istat = gmm_get('TR/'//trim(Tr3d_name_S(n))//':'//F_timelevel_S,tr)

            do k = 1, l_nk
               do j = 1, l_nj
                  do i = 1, l_ni
                     F_qh(i,j,k)=F_qh(i,j,k)+tr(i,j,k)
                  end do
               end do
            end do
         end if
      end do

!     ________________________________________________________________
!
      return
      end
