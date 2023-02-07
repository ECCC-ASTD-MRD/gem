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

      subroutine pw_update_UV_hlt()
      use glb_ld
      use gmm_pw
      use gmm_vt1
      use glb_ld
      use outp
      use omp_timing
      implicit none

      integer i,j
!     ________________________________________________________________
!
      call gtmg_start (5, 'PW_UPDATE', 0)

      call hwnd_stag_hlt ( pw_uu_plus,pw_vv_plus, ut1,vt1, &
                           l_minx,l_maxx,l_miny,l_maxy,l_nk,.false.)
!$omp do
      do j=l_miny, l_maxy
         do i=l_minx,l_maxx
            udiag(i,j) = pw_uu_plus(i,j,G_nk)
            vdiag(i,j) = pw_vv_plus(i,j,G_nk)
         end do
      end do
!$omp end do
      
      call gtmg_stop (5)
!     ________________________________________________________________
!
      return
      end subroutine pw_update_UV_hlt
