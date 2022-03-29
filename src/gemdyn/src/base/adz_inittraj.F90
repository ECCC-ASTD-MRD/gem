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

      subroutine adz_inittraj
      use glb_ld
      use adz_options
      use adz_mem

      integer i,j,k
!
!---------------------------------------------------------------------
!
      Adz_niter = max( 6, Adz_itraj )
      do k = Adz_k0m, l_nk
         do j = 1, l_nj
            do i = 1, l_ni
               Adz_pxyzm(1,i,j,k) = i + l_i0 - 1
               Adz_wpxyz(i,j,k,1) = Adz_pxyzm(1,i,j,k)
               Adz_pxyzm(2,i,j,k) = j + l_j0 - 1
               Adz_wpxyz(i,j,k,2) = Adz_pxyzm(2,i,j,k)
               Adz_pxyzm(3,i,j,k) = k
               Adz_wpxyz(i,j,k,3) = Adz_pxyzm(3,i,j,k)
            end do
         end do
      end do
!
!---------------------------------------------------------------------
!
      return
      end subroutine adz_inittraj
