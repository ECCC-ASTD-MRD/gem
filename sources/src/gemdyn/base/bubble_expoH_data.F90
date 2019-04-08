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
!**s/r bubble_expoH_data - generates initial condition for Robert's bubble
!                    experiment (Robert 1993 JAS) - EXPO height coord.
!
      subroutine bubble_expoH_data ( F_t,Mminx,Mmaxx,Mminy,Mmaxy,nk )
      use glb_ld
      use bubble_options
      implicit none
#include <arch_specific.hf>

      integer Mminx,Mmaxx,Mminy,Mmaxy,nk
      real F_t    (Mminx:Mmaxx,Mminy:Mmaxy,nk)

      integer :: i,j,k,ii
      real*8 :: theta, r,rad
!
!     ---------------------------------------------------------------
!
      if(.not.bubble_gaus_L) then

         do k=1,g_nk
            do j=1,l_nj
            do i=1,l_ni
               ii=i+l_i0-1
               theta = bubble_theta
               if ( (((ii)-bubble_ictr)**2 +((k)-bubble_kctr)**2) < bubble_rad**2 ) then
                   theta = theta + 0.5d0
               end if
               F_t(i,j,k) = theta
            end do
            end do
         end do

      else

      ! Gaussian

         do k=1,g_nk
            do j=1,l_nj
            do i=1,l_ni
               ii=i+l_i0-1
               r = sqrt( ( ((ii)-bubble_ictr) * dble(bubble_dx) )**2 + ( ((k)-bubble_kctr) * dble(bubble_dz) )**2)
               rad = bubble_rad*dble(bubble_dx)
               if ( r <= rad ) then
                   theta = bubble_theta + 0.5d0
               else
                   theta = bubble_theta + 0.5d0 * exp( -(r - rad)**2 / 100.d0**2 );
               end if
               F_t(i,j,k) = theta
            end do
            end do
         end do

      end if
!
!     -----------------------------------------------------------------
!
      return
      end subroutine bubble_expoH_data
